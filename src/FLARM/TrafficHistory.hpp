// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright The XCSoar Project

#pragma once

#include "FLARM/Id.hpp"
#include "Geo/GeoPoint.hpp"
#include "time/Stamp.hpp"

#include <map>
#include <deque>
#include <chrono>

/**
 * Stores historical positions for FLARM targets
 */
class FlarmTrafficHistory {
public:
  struct HistoryPoint {
    GeoPoint location;
    TimeStamp time;

    constexpr HistoryPoint(const GeoPoint &_location, TimeStamp _time) noexcept
      :location(_location), time(_time) {}
  };

  using TrailPoints = std::deque<HistoryPoint>;

private:
  std::map<FlarmId, TrailPoints> trails;

public:
  /**
   * Add a new position to the history
   * 
   * @param id FLARM ID
   * @param location Current location
   * @param time Current time
   */
  void AddPoint(FlarmId id, const GeoPoint &location, TimeStamp time) noexcept;

  /**
   * Get the trail for a specific FLARM target
   * 
   * @param id FLARM ID
   * @return Array of trail points or nullptr if no trail exists
   */
  const TrailPoints* GetTrail(FlarmId id) const noexcept;

  /**
   * Remove expired points (older than 15 minutes)
   * 
   * @param now Current time
   */
  void Expire(TimeStamp now) noexcept;
};
