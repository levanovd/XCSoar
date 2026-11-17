// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright The XCSoar Project

#pragma once

#include "Geo/GeoPoint.hpp"
#include "Geo/Flat/FlatPoint.hpp"
#include "Math/Angle.hpp"
#include "time/Stamp.hpp"
#include "NMEA/ThermalLocator.hpp"

#include <array>

struct SpeedVector;

/**
 * Locates the centre of a thermal by fitting a bell shaped lift,
 * around + strength * exp(-d^2 / radius^2), to the netto vario of the
 * last two minutes at the places it was measured, in coordinates that
 * drift with the wind.  Every centre of a grid gets its best strength
 * and air around in closed form, so each has a likelihood; the most
 * probable one wins.  Unlike the lift weighted mean of the flown
 * positions, the centre can lie outside the circle flown.
 *
 * Over the height the thermal is one column, leaning by a slope that
 * is fitted to the samples too, and a centre or cross-section unlike
 * the column so far has to fit better.
 *
 * The radius comes from all the samples: those at many distances from
 * the centre show it, the last circles around it hardly do.  Strength
 * and air around are those of the recent samples at that radius.
 *
 * Needs the netto of a total energy vario: from the GPS altitude the
 * lift comes too late and too noisy, the fit is then no better than
 * the mean of the last circle.  A new thermal is only found in a
 * turn; one found stays while turning, also in sink, and through
 * short straight stretches.
 */
class ThermalFit {
public:
  ThermalFit() noexcept {
    Reset();
  }

  void Reset() noexcept;

  /**
   * @param turn_rate heading turn rate per second, low passed
   * @param circling_mode the computer's circling mode: a thermal
   * found stays while it lasts, so it is still there when the climb
   * is recorded at its end
   * @param speed airspeed, or ground speed without one [m/s]
   * @param altitude [m]
   * @param netto netto vario of a total energy vario [m/s]
   */
  void Process(Angle turn_rate, bool circling_mode, double speed,
               TimeStamp time,
               const GeoPoint &location, double altitude, double netto,
               const SpeedVector &wind,
               ThermalLocatorInfo &therm) noexcept;

private:
  static constexpr unsigned MAX_SAMPLES = 240;
  static constexpr unsigned MAX_TRACK = 32;

  struct Sample {
    TimeStamp time;
    /** where the air of this vario reading was */
    GeoPoint location;
    double altitude, w;
  };

  struct Fit {
    double cost, strength, radius, around;
    /** the cost with what an unlike cross-section has to fit better by */
    double total;
    /** the bell's mean over the samples, of its top */
    double mean;
    /** of the strength */
    double variance;
  };

  struct Weights {
    /** each sample's weight, and times its vario */
    std::array<double, MAX_SAMPLES> a, aw;
    double s1, sw, sww;

    void Start() noexcept;
    void Add(unsigned k, double weight, double w) noexcept;
  };

  /** recent positions, for where the late vario reading belongs */
  std::array<std::pair<TimeStamp, GeoPoint>, MAX_TRACK> track;
  unsigned track_head, track_count;

  std::array<Sample, MAX_SAMPLES> samples;
  unsigned head, count;

  /** the readings of the sample being collected: netto, and time
      since sum_start [s] */
  double sum_w, sum_t;
  unsigned sum_n;
  TimeStamp sum_start, last_time, last_search;

  /* the samples of one search, in metres east (x) and north (y) of
     the aircraft, drifted to now */
  std::array<FlatPoint, MAX_SAMPLES> xy;
  std::array<double, MAX_SAMPLES> z, w;
  /** xy at the search altitude, along the column's lean */
  std::array<FlatPoint, MAX_SAMPLES> leaned;
  /** for the centre, strength and air around; for the radius */
  Weights recent, wide;
  double search_altitude;
  FlatPoint column;

  /** the last search's centre, before the low pass */
  GeoPoint found;
  /** the centre given, after the low pass */
  GeoPoint centre;
  /** the altitude and time of the last search that found the thermal */
  double found_altitude;
  TimeStamp found_time;
  /** the column's lean: metres east and north per metre of altitude */
  FlatPoint lean;
  /** the radius of the last cross-section found [m] */
  double section_radius;
  /** the strength, a state that drifts slowly, and its variance, then
     and grown to the search */
  double strength, strength_variance, search_variance;
  /** the lift given, after the low pass */
  LiftProfile profile;
  bool valid;
  unsigned missed;

  GeoPoint TrackAt(TimeStamp time) const noexcept;
  void AddSample(TimeStamp time, double altitude, double lag) noexcept;
  /**
   * @param period the time of a turn [s], 0 when not circling
   */
  void Search(double period, bool circling_mode, TimeStamp time,
              const GeoPoint &location, double altitude,
              const GeoPoint &drift) noexcept;
  /** the samples at the search altitude, along the lean */
  void Lean(FlatPoint lean,
            std::array<FlatPoint, MAX_SAMPLES> &points) const noexcept;
  /**
   * @param points the samples, leaned
   * @param hold_strength whether a strength unlike the thermal's so
   * far has to fit better
   */
  Fit FitAt(FlatPoint c, const std::array<FlatPoint, MAX_SAMPLES> &points,
            const Weights &weights, bool hold_strength) const noexcept;
  Fit FitRadius(FlatPoint c, double radius) const noexcept;
  static Fit Solve(double sg, double sgg, double swg, double radius,
                   const Weights &weights) noexcept;
  double CostAt(FlatPoint c) const noexcept;
  FlatPoint FitLean(FlatPoint c) const noexcept;
};
