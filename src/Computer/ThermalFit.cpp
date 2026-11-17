// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright The XCSoar Project

#include "ThermalFit.hpp"
#include "Geo/Flat/FlatProjection.hpp"
#include "Geo/Math.hpp"
#include "Geo/SpeedVector.hpp"
#include "Math/Constants.hpp"
#include "Math/LowPassFilter.hpp"
#include "time/Cast.hxx"

#include <algorithm>
#include <cmath>

/* tuned on recorded XCVario flights, see the XCVarioPure AirMass
   analysis (RESULTS-2026-11-14-thermal.md) */

/** a sample is the mean vario of this long [s] */
static constexpr double SAMPLE_PERIOD = 0.5;
/** a sample's weight fades with its age [s], circle by circle while
    circling: fading with the age alone, the side flown last would
    count most and the centre would swing around with the heading */
static constexpr double MEMORY = 40;
/** the netto shows the air of this long ago [s], more while circling
    (the vario's damping, the transfer, and a GPS fix that is older
    than its arrival) */
static constexpr double LAG = 2;
static constexpr double CIRCLING_LAG = 0.5;
static constexpr double SEARCH_PERIOD = 1; // [s]
/** samples older than this are dropped [s] */
static constexpr double HISTORY = 120;

/** turning in a circle: faster than this [rad/s], on a radius below
    this [m].  Not the circling mode, which a vario's cruise/climb
    switch can force in straight flight and which stays on for a while
    after the last turn */
static constexpr double TURN_RATE_MIN = 0.0698; // 4 degrees/s
static constexpr double TURN_RADIUS_MAX = 200;
static constexpr unsigned MIN_SAMPLES = 8;

/** the widest thermal [m], each narrower one by 2^-1/4 */
static constexpr double RADIUS = 400;
static constexpr double NARROWER = 0.84089641525371454;
/** a sample's weight for the radius fades with its age [s] */
static constexpr double RADIUS_MEMORY = 120;
/** prior of the radius [m], log-normal: a few circles alone would fit
    nearly any width */
static constexpr double RADIUS_PRIOR = 160;
static constexpr double RADIUS_SIGMA = 0.5;
static constexpr double COARSE_STEP = 40; // [m]
static constexpr double FINE_STEP = 10; // [m]

/** vario error of a sample, as if each were independent [m/s] */
static constexpr double NOISE = 1;
/** priors of the air around and of the strength [m/s] */
static constexpr double AROUND_SIGMA = 1.5;
static constexpr double STRENGTH_SIGMA = 4;
/** the air around a thermal does not rise: steady lift all around the
    circle is a thermal centred on it [m/s] */
static constexpr double AROUND_MAX = 0.2;
static constexpr double STRENGTH_MAX = 8;

/** a thermal is found with this much lift above the air around [m/s]
    and log likelihood ratio against none; one found stays with less,
    and through searches without one */
static constexpr double STRENGTH_MIN = 0.5;
static constexpr double EVIDENCE_MIN = 6;
static constexpr double STRENGTH_KEEP = 0.25;
static constexpr double EVIDENCE_KEEP = 2;
static constexpr unsigned HOLD = 5;

/** share of the fitted strength given, the rest goes to the air
    around so the mean along the path stays: on real flights the lift
    differences fitted came half as large a circle later, the lift
    around a circle is mostly gusts.  The centre comes from the
    differences and is not held to this */
static constexpr double GIVEN = 0.4;

/** low pass of the thermal given [s]: centres that fit nearly alike
    take turns from search to search */
static constexpr double SMOOTH = 3;
/** the search stays around a thermal this close [m] */
static constexpr double FOLLOW = 250;
/** a centre this far off the column [m], or a radius twice as wide
    or half, has to fit better by one noise level */
static constexpr double CONTINUITY = 80;
static constexpr double SECTION_RADIUS = 0.7;

/** the strength is a state that drifts by this much in a minute and
    in 100 m of climb [m/s]; each search's fit updates it, and a
    strength off it has to fit better by its uncertainty */
static constexpr double STRENGTH_DRIFT = 4;
static constexpr double STRENGTH_DRIFT_HEIGHT = 2;
/** successive searches share most of their samples: each fit counts
    as this share of an independent one */
static constexpr double STRENGTH_OVERLAP = MEMORY / SEARCH_PERIOD;

/** the column's lean [m/m]: prior, the leans tried either side, and
    the share of the way to the best one taken per search */
static constexpr double LEAN_SIGMA = 0.5;
static constexpr double LEAN_STEP = 0.1;
static constexpr double LEAN_GAIN = 0.3;
static constexpr double LEAN_MAX = 2;

static constexpr unsigned RADII = 12;
static constexpr int COARSE = 13;
static constexpr int FINE = 7;

static constexpr double PRIOR_AROUND =
  NOISE * NOISE / (AROUND_SIGMA * AROUND_SIGMA);
static constexpr double PRIOR_STRENGTH =
  NOISE * NOISE / (STRENGTH_SIGMA * STRENGTH_SIGMA);

void
ThermalFit::Weights::Start() noexcept
{
  s1 = PRIOR_AROUND;
  sw = sww = 0;
}

void
ThermalFit::Weights::Add(unsigned k, double weight, double w) noexcept
{
  a[k] = weight;
  aw[k] = weight * w;
  s1 += weight;
  sw += weight * w;
  sww += weight * w * w;
}

void
ThermalFit::Reset() noexcept
{
  track_head = track_count = 0;
  head = count = 0;
  sum_w = sum_t = 0;
  sum_n = 0;
  valid = false;
  missed = 0;
  lean = FlatPoint(0, 0);
  profile = {};
}

GeoPoint
ThermalFit::TrackAt(TimeStamp time) const noexcept
{
  /* newest first */
  for (unsigned i = 0; i + 1 < track_count; ++i) {
    const auto &newer =
      track[(track_head + MAX_TRACK - 1 - i) % MAX_TRACK];
    const auto &older =
      track[(track_head + MAX_TRACK - 2 - i) % MAX_TRACK];
    if (older.first <= time) {
      const double span = ToFloatSeconds(newer.first - older.first);
      const double f = span > 0
        ? std::clamp(ToFloatSeconds(time - older.first) / span, 0., 1.)
        : 1.;
      return older.second.Interpolate(newer.second, f);
    }
  }

  return track[(track_head + MAX_TRACK - track_count) % MAX_TRACK].second;
}

void
ThermalFit::AddSample(TimeStamp time, double altitude, double lag) noexcept
{
  /* when and where the air of this reading was: it drifts with the
     wind from then */
  Sample &s = samples[head];
  s.time = time - FloatDuration{lag};
  s.location = TrackAt(s.time);
  s.altitude = altitude;
  s.w = sum_w / sum_n;
  head = (head + 1) % MAX_SAMPLES;
  count = std::min(count + 1, MAX_SAMPLES);
}

void
ThermalFit::Process(Angle turn_rate, bool circling_mode, double speed,
                    TimeStamp time,
                    const GeoPoint &location, double altitude,
                    double netto, const SpeedVector &wind,
                    ThermalLocatorInfo &therm) noexcept
{
  const double rate = std::fabs(turn_rate.Radians());
  const bool circling = rate > TURN_RATE_MIN &&
    speed < TURN_RADIUS_MAX * rate;

  if (track_count > 0 && (time < last_time ||
                          time - last_time > std::chrono::seconds{10}))
    Reset();

  /* the air drifts with the wind; drift per second */
  const GeoPoint drift = location -
    FindLatitudeLongitude(location, wind.bearing, wind.norm);

  if (track_count > 0) {
    const double dt = ToFloatSeconds(time - last_time);
    if (dt <= 0)
      return;

    if (valid) {
      found += drift * dt;
      centre += drift * dt;
    }
  } else {
    sum_start = time;
    last_search = time;
  }
  last_time = time;

  track[track_head] = {time, location};
  track_head = (track_head + 1) % MAX_TRACK;
  track_count = std::min(track_count + 1, MAX_TRACK);

  sum_w += netto;
  sum_t += ToFloatSeconds(time - sum_start);
  sum_n++;
  if (ToFloatSeconds(time - sum_start) >= SAMPLE_PERIOD) {
    /* at the mean time of the readings */
    AddSample(sum_start + FloatDuration{sum_t / sum_n}, altitude,
              circling ? LAG + CIRCLING_LAG : LAG);
    sum_w = sum_t = 0;
    sum_n = 0;
    sum_start = time;
  }

  while (count > 0 &&
         ToFloatSeconds(time - samples[(head + MAX_SAMPLES - count) %
                                       MAX_SAMPLES].time) > HISTORY)
    --count;

  if (count >= MIN_SAMPLES &&
      ToFloatSeconds(time - last_search) >= SEARCH_PERIOD) {
    last_search = time;
    Search(circling ? M_2PI / rate : 0, circling_mode, time, location,
           altitude, drift);
  }

  therm.estimate_valid = valid;
  therm.estimate_profile_valid = valid;
  therm.estimate_profile = profile;
  if (valid) {
    const FlatProjection projection(centre);
    const double scale = projection.GetApproximateScale();
    const FlatPoint shift = lean * ((altitude - found_altitude) / scale);
    therm.estimate_location = projection.Unproject(shift);
  }
}

/**
 * The lowest point of the parabola through three costs a step apart,
 * in steps from the middle one, within +/- limit; flat if it is not
 * convex.
 */
static double
ParabolaVertex(double below, double at, double above, double limit,
               double flat = 0.) noexcept
{
  const double bend = below + above - 2 * at;
  return bend > 1e-6
    ? std::clamp(0.5 * (below - above) / bend, -limit, limit)
    : flat;
}

ThermalFit::Fit
ThermalFit::Solve(double sg, double sgg, double swg, double radius,
                  const Weights &weights) noexcept
{
  /* weighted least squares of base + strength * g, both with a prior
     around zero */
  const double s1 = weights.s1, sw = weights.sw;
  const double held = sgg + PRIOR_STRENGTH;
  const double det = s1 * held - sg * sg;
  double strength = det > 1e-6
    ? std::clamp((s1 * swg - sg * sw) / det, 0., STRENGTH_MAX)
    : 0.;
  double around = (sw - sg * strength) / s1;
  double variance = NOISE * NOISE * (det > 1e-6 ? s1 / det : 1 / held);
  if (around > AROUND_MAX) {
    around = AROUND_MAX;
    strength = std::clamp((swg - around * sg) / held, 0., STRENGTH_MAX);
    variance = NOISE * NOISE / held;
  }

  const double cost = weights.sww - 2 * around * sw - 2 * strength * swg +
    around * around * s1 + 2 * around * strength * sg +
    strength * strength * held;
  return {cost, strength, radius, around, cost, sg / s1, variance};
}

void
ThermalFit::Lean(FlatPoint l,
                 std::array<FlatPoint, MAX_SAMPLES> &points) const noexcept
{
  for (unsigned k = 0; k < count; ++k)
    points[k] = xy[k] + l * (search_altitude - z[k]);
}

ThermalFit::Fit
ThermalFit::FitAt(FlatPoint c,
                  const std::array<FlatPoint, MAX_SAMPLES> &points,
                  const Weights &weights, bool hold_strength) const noexcept
{
  /* two chains a quarter octave apart, each narrower radius of a chain
     the square of the one before */
  double sg[RADII] = {}, sgg[RADII] = {}, swg[RADII] = {};
  for (unsigned k = 0; k < count; ++k) {
    const FlatPoint p = points[k] - c;
    const double d2 = p.DotProduct(p) / (RADIUS * RADIUS);
    for (unsigned chain = 0; chain < 2; ++chain) {
      double g = std::exp(chain == 0 ? -d2 : -d2 / (NARROWER * NARROWER));
      for (unsigned r = chain; r < RADII; r += 2, g *= g) {
        sg[r] += weights.a[k] * g;
        sgg[r] += weights.a[k] * g * g;
        swg[r] += weights.aw[k] * g;
      }
    }
  }

  Fit best{};
  best.total = INFINITY;
  double total[RADII];
  unsigned narrowest = 0;
  double radius = RADIUS;
  for (unsigned r = 0; r < RADII; ++r, radius *= NARROWER) {
    const Fit f = Solve(sg[r], sgg[r], swg[r], radius, weights);
    const double usual = std::log(radius / RADIUS_PRIOR) / RADIUS_SIGMA;
    total[r] = f.cost + usual * usual * NOISE * NOISE;
    if (valid) {
      const double wider = std::log(radius / section_radius) / SECTION_RADIUS;
      total[r] += wider * wider * NOISE * NOISE;
      /* not in the radius fit: a wider bell with less strength fits
         nearly as well, and held to the strength the radius would
         follow it */
      if (hold_strength) {
        const double stronger = f.strength - strength;
        total[r] += stronger * stronger / search_variance * NOISE * NOISE;
      }
    }

    if (total[r] < best.total) {
      best = f;
      best.total = total[r];
      narrowest = r;
    }
  }

  /* between the radii tried: the lowest point of the parabola through
     the best one and its neighbours */
  if (narrowest > 0 && narrowest < RADII - 1)
    best.radius *= std::pow(NARROWER,
                            ParabolaVertex(total[narrowest - 1],
                                           total[narrowest],
                                           total[narrowest + 1], 0.5));

  return best;
}

ThermalFit::Fit
ThermalFit::FitRadius(FlatPoint c, double radius) const noexcept
{
  double sg = 0, sgg = 0, swg = 0;
  for (unsigned k = 0; k < count; ++k) {
    const FlatPoint p = leaned[k] - c;
    const double g = std::exp(-p.DotProduct(p) / (radius * radius));
    sg += recent.a[k] * g;
    sgg += recent.a[k] * g * g;
    swg += recent.aw[k] * g;
  }
  return Solve(sg, sgg, swg, radius, recent);
}

double
ThermalFit::CostAt(FlatPoint c) const noexcept
{
  double cost = FitAt(c, leaned, recent, true).total;
  if (valid) {
    const FlatPoint off = c - column;
    cost += off.DotProduct(off) *
      (NOISE * NOISE / (CONTINUITY * CONTINUITY));
  }
  return cost;
}

FlatPoint
ThermalFit::FitLean(FlatPoint c) const noexcept
{
  /* a step towards the lean that fits the samples best, east and
     north; the fits do not use the new lean, so a wrong one cannot
     feed itself */
  FlatPoint result = lean;
  std::array<FlatPoint, MAX_SAMPLES> points;
  for (double FlatPoint::*axis : {&FlatPoint::x, &FlatPoint::y}) {
    double cost[3];
    for (int k = 0; k < 3; ++k) {
      FlatPoint l = lean;
      l.*axis += (k - 1) * LEAN_STEP;
      Lean(l, points);
      cost[k] = FitAt(c, points, recent, true).total +
        l.*axis * l.*axis * (NOISE * NOISE / (LEAN_SIGMA * LEAN_SIGMA));
    }

    const double to = ParabolaVertex(cost[0], cost[1], cost[2], 1.,
                                     cost[0] < cost[2] ? -1. : 1.);
    result.*axis = std::clamp(lean.*axis + LEAN_GAIN * LEAN_STEP * to,
                              -LEAN_MAX, LEAN_MAX);
  }
  return result;
}

void
ThermalFit::Search(double period, bool circling_mode, TimeStamp time,
                   const GeoPoint &location, double altitude,
                   const GeoPoint &drift) noexcept
{
  /* a new thermal is only found in a turn: lift crossed in straight
     flight is a place along the track with either side as probable,
     and most of it is a gust that no circle confirms */
  const bool circling = period > 0;
  if (!valid && !circling)
    return;

  const FlatProjection projection(location);
  const double scale = projection.GetApproximateScale();
  auto ToMetres = [&](const GeoPoint &p) {
    return projection.ProjectFloat(p) * scale;
  };
  auto FromMetres = [&](FlatPoint p) {
    return projection.Unproject(p * (1. / scale));
  };

  /* the samples drifted with the air to now, in metres from the
     aircraft; the centre of the circle flown is their mean over its
     last turn */
  recent.Start();
  wide.Start();
  FlatPoint circle(0, 0);
  unsigned in_circle = 0;
  for (unsigned k = 0; k < count; ++k) {
    const Sample &s = samples[(head + MAX_SAMPLES - 1 - k) % MAX_SAMPLES];
    const double age = ToFloatSeconds(time - s.time);
    xy[k] = ToMetres(s.location + drift * age);
    z[k] = s.altitude;
    w[k] = s.w;
    recent.Add(k, std::exp(-(period > 0 ? std::floor(age / period) * period
                                        : age) / MEMORY), w[k]);
    wide.Add(k, std::exp(-age / RADIUS_MEMORY), w[k]);
    if (age < period) {
      circle += xy[k];
      ++in_circle;
    }
  }

  search_altitude = altitude;
  Lean(lean, leaned);
  if (valid)
    search_variance = strength_variance +
      STRENGTH_DRIFT * STRENGTH_DRIFT *
      ToFloatSeconds(time - found_time) / 60 +
      STRENGTH_DRIFT_HEIGHT * STRENGTH_DRIFT_HEIGHT *
      std::fabs(altitude - found_altitude) / 100;

  FlatPoint mid(0, 0);
  if (in_circle > 0)
    mid = circle * (1. / in_circle);
  if (valid) {
    column = ToMetres(found) + lean * (altitude - found_altitude);
    if ((column - mid).Magnitude() < FOLLOW)
      mid = column;
  }

  /* coarse grid, then a fine one around its best cell, then the
     lowest point of the parabolas through the best fine cell and its
     neighbours: a centre held to the cells would stay while the
     thermal moves and then jump */
  auto Cell = [](FlatPoint middle, int size, double step, int i, int j) {
    return middle + FlatPoint(i - size / 2, j - size / 2) * step;
  };

  double coarse_best = INFINITY;
  int bi = 0, bj = 0;
  for (int i = 0; i < COARSE; ++i)
    for (int j = 0; j < COARSE; ++j) {
      const double cost = CostAt(Cell(mid, COARSE, COARSE_STEP, i, j));
      if (cost < coarse_best) {
        coarse_best = cost;
        bi = i;
        bj = j;
      }
    }

  const FlatPoint fine_mid = Cell(mid, COARSE, COARSE_STEP, bi, bj);
  double fine[FINE][FINE];
  int fi = 0, fj = 0;
  for (int i = 0; i < FINE; ++i)
    for (int j = 0; j < FINE; ++j) {
      fine[i][j] = CostAt(Cell(fine_mid, FINE, FINE_STEP, i, j));
      if (fine[i][j] < fine[fi][fj]) {
        fi = i;
        fj = j;
      }
    }

  FlatPoint best = Cell(fine_mid, FINE, FINE_STEP, fi, fj);
  if (fi > 0 && fi < FINE - 1)
    best.x += ParabolaVertex(fine[fi - 1][fj], fine[fi][fj],
                             fine[fi + 1][fj], 0.5) * FINE_STEP;
  if (fj > 0 && fj < FINE - 1)
    best.y += ParabolaVertex(fine[fi][fj - 1], fine[fi][fj],
                             fine[fi][fj + 1], 0.5) * FINE_STEP;

  const Fit fit = FitAt(best, leaned, recent, true);

  /* against no thermal: the air around, at most AROUND_MAX */
  const double calm = std::min(recent.sw / recent.s1, AROUND_MAX);
  const double evidence = (recent.sww - 2 * calm * recent.sw +
                           calm * calm * recent.s1 - fit.cost) /
    (2 * NOISE * NOISE);
  const bool inside = bi > 0 && bi < COARSE - 1 && bj > 0 && bj < COARSE - 1;

  const bool found_now = inside &&
    (valid
     ? fit.strength >= STRENGTH_KEEP && evidence >= EVIDENCE_KEEP
     : fit.strength >= STRENGTH_MIN && evidence >= EVIDENCE_MIN);
  if (!found_now) {
    /* while still turning the thermal stays where it was: the lift
       weakening or the circle drifting into sink does not end it; and
       until the circling mode ends, which records the climb with it */
    if (valid && (circling || circling_mode))
      return;

    valid = valid && ++missed <= HOLD;
    if (!valid)
      lean = FlatPoint(0, 0);
    return;
  }

  const Fit cross =
    FitRadius(best, FitAt(best, leaned, wide, false).radius);

  /* zero while not valid */
  if (valid)
    lean = FitLean(best);

  const double since = valid ? ToFloatSeconds(time - found_time) : 0.;
  const double a = valid ? since / (SMOOTH + since) : 1.;
  const FlatPoint before = valid
    ? ToMetres(centre) + lean * (altitude - found_altitude)
    : best;

  const double measured = cross.variance * STRENGTH_OVERLAP;
  if (valid) {
    const double gain = search_variance / (search_variance + measured);
    strength += (cross.strength - strength) * gain;
    strength_variance = (1 - gain) * search_variance;
  } else {
    strength = cross.strength;
    strength_variance = measured;
  }

  /* the strength given, the rest to the air around so the mean along
     the path stays */
  const double given = strength * GIVEN;
  const double around = cross.around + (cross.strength - given) * cross.mean;
  const LiftProfile now{around + given, around, cross.radius};
  profile.lift = LowPassFilter(profile.lift, now.lift, a);
  profile.around = LowPassFilter(profile.around, now.around, a);
  profile.radius = LowPassFilter(profile.radius, now.radius, a);

  found = FromMetres(best);
  centre = FromMetres(before + (best - before) * a);
  found_altitude = altitude;
  found_time = time;
  section_radius = cross.radius;
  missed = 0;
  valid = true;
}
