// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright The XCSoar Project

#include "Computer/ThermalFit.hpp"
#include "Geo/Flat/FlatProjection.hpp"
#include "Geo/SpeedVector.hpp"
#include "NMEA/ThermalLocator.hpp"
#include "TestUtil.hpp"

#include <cmath>

/* On noise free lift the first circle's centre is held by the
   column's continuity (a centre far off it has to fit better), so the
   centre found stays tens of metres off; on real flights the noise
   hides that. */

/**
 * Circle or fly straight at 1 Hz in a bell shaped thermal that drifts
 * with the wind, the netto showing the air 2.5 s late, and return how
 * far the estimate is from the thermal's centre [m]; -1 if there is
 * none.
 */
struct Flight {
  /** the thermal in the air [m east, m north], its radius and lift */
  double cx, cy, radius = 160, strength = 3;
  /** the wind [m/s], blowing towards the east */
  double wind = 0;
  double bank_rate = 0.33;
  double straight_rate = 0;

  const GeoPoint origin{Angle::Degrees(7), Angle::Degrees(51)};
  const FlatProjection projection{origin};
  ThermalFit fit;
  ThermalLocatorInfo info;
  double x = 0, y = 0, psi = 0, t = 0;
  double history[4] = {};

  Flight(double _cx, double _cy) noexcept : cx(_cx), cy(_cy) {
    info.Clear();
  }

  /** the air around the thermal [m/s] */
  double around = -0.5;

  double Lift(double px, double py) const noexcept {
    const double d2 = (px - cx) * (px - cx) + (py - cy) * (py - cy);
    return around + strength * std::exp(-d2 / (radius * radius));
  }

  /**
   * Turning at bank_rate while circling, else at straight_rate; the
   * computer's circling mode follows the turning unless given.
   */
  void Fly(double seconds, bool circling) noexcept {
    Fly(seconds, circling, circling);
  }

  void Fly(double seconds, bool circling, bool circling_mode) noexcept {
    const double rate = circling ? bank_rate : straight_rate;
    for (double end = t + seconds; t < end; t += 1) {
      psi += rate;
      x += 25 * std::sin(psi);
      y += 25 * std::cos(psi);
      for (unsigned i = 3; i > 0; --i)
        history[i] = history[i - 1];
      history[0] = Lift(x, y);
      const double netto = (history[2] + history[3]) / 2;

      const double scale = projection.GetApproximateScale();
      const GeoPoint ground =
        projection.Unproject(FlatPoint(x + wind * t, y) * (1. / scale));
      fit.Process(Angle::Radians(rate), circling_mode, 25,
                  TimeStamp{FloatDuration{t}},
                  ground, 1000, netto,
                  SpeedVector(Angle::Degrees(270), wind), info);
    }
  }

  double Error() const noexcept {
    if (!info.estimate_valid)
      return -1;
    const double scale = projection.GetApproximateScale();
    const FlatPoint p = projection.ProjectFloat(info.estimate_location) *
      scale;
    return std::hypot(p.x - (cx + wind * t), p.y - cy);
  }
};

static void
TestCircle()
{
  /* the circle flown, about 76 m wide, is centred at (77, -15) */
  Flight f(40, 80);
  f.Fly(70, true);
  ok(f.Error() >= 0 && f.Error() < 55, "circle off the centre: found, %.0f m",
     f.Error());
}

static void
TestWind()
{
  Flight f(130, 60);
  f.wind = 8;
  f.Fly(70, true);
  ok(f.Error() >= 0 && f.Error() < 70,
     "drifting with 8 m/s of wind: found, %.0f m", f.Error());
}

static void
TestStraight()
{
  Flight f(300, 0);
  f.psi = 1.5707963;
  f.Fly(20, false);
  ok1(f.Error() < 0);

  f.Fly(25, true);
  ok(f.Error() >= 0 && f.Error() < 60,
     "found once circling, %.0f m", f.Error());
}

static void
TestGentleTurn()
{
  /* gliding through lift in a wide turn, 4.6 degrees/s on a radius of
     more than 300 m: no thermal */
  Flight f(300, 0);
  f.psi = 1.5707963;
  f.straight_rate = 0.08;
  f.Fly(25, false);
  ok(f.Error() < 0, "gentle turn through lift: no thermal");
}

static void
TestDying()
{
  /* the thermal dies while the glider keeps turning: the centre stays
     while turning, and goes after a few seconds of straight flight */
  Flight f(40, 80);
  f.Fly(40, true);
  f.strength = 0;
  f.around = -1.5;
  f.Fly(40, true);
  ok1(f.Error() >= 0);
  f.Fly(30, false);
  ok1(f.Error() < 0);
}

static void
TestModeHold()
{
  /* leaving the thermal straight: the estimate stays until the
     circling mode ends, as the climb is recorded then, and goes a few
     seconds later */
  Flight f(40, 80);
  f.Fly(40, true);
  f.Fly(25, false, true);
  ok(f.Error() >= 0, "straight, still in circling mode: kept");
  f.Fly(30, false);
  ok(f.Error() < 0, "circling mode ended: gone");
}

static void
TestNoLift()
{
  Flight f(0, 0);
  f.strength = 0;
  unsigned found = 0;
  for (unsigned i = 0; i < 90; ++i) {
    f.Fly(1, i > 20);
    found += f.Error() >= 0;
  }
  ok(found == 0, "no lift: no thermal (%u s with one)", found);
}

int
main()
{
  plan_tests(10);
  TestCircle();
  TestWind();
  TestStraight();
  TestGentleTurn();
  TestDying();
  TestModeHold();
  TestNoLift();
  return exit_status();
}
