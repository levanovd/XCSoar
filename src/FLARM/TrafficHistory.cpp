// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright The XCSoar Project

#include "TrafficHistory.hpp"

void
FlarmTrafficHistory::AddPoint(FlarmId id, const GeoPoint &location, 
                              TimeStamp time) noexcept
{
  // Get or create the trail for this ID
  auto &trail = trails[id];
  
  // Only add a point if it's been at least 15 seconds since the last one
  // or if this is the first point
  if (trail.empty() || 
      (time - trail.back().time) >= std::chrono::seconds(15)) {
    
    // If the trail has more than 60 points, remove the oldest point
    if (trail.size() >= 60)
      trail.pop_front();
    
    // Add the new point
    trail.push_back(HistoryPoint(location, time));
  }
}

const FlarmTrafficHistory::TrailPoints*
FlarmTrafficHistory::GetTrail(FlarmId id) const noexcept
{
  auto it = trails.find(id);
  if (it != trails.end())
    return &it->second;
  
  return nullptr;
}

void
FlarmTrafficHistory::Expire(TimeStamp now) noexcept
{
  // 15 minutes expiration time
  const auto max_age = std::chrono::minutes(15);
  
  for (auto &pair : trails) {
    auto &trail = pair.second;
    
    // Remove points older than 15 minutes
    while (!trail.empty() && (now - trail.front().time) > max_age)
      trail.pop_front();
  }
  
  // Remove empty trails
  for (auto it = trails.begin(); it != trails.end();) {
    if (it->second.empty())
      it = trails.erase(it);
    else
      ++it;
  }
}
