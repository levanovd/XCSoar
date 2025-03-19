// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright The XCSoar Project

#include "MapWindow.hpp"
#include "ui/canvas/Canvas.hpp"
#include "ui/canvas/Pen.hpp"
#include "Screen/Layout.hpp"
#include "Renderer/TrafficRenderer.hpp"
#include "FLARM/Friends.hpp"
#include "FLARM/Traffic.hpp"
#include "Look/TrafficLook.hpp"
#include "Tracking/SkyLines/Data.hpp"
#include "ui/canvas/Color.hpp"


/**
 * Draws the FLARM traffic icons onto the given canvas
 * @param canvas Canvas for drawing
 */
void
MapWindow::DrawFLARMTraffic(Canvas &canvas,
                            [[maybe_unused]] const PixelPoint aircraft_pos) const noexcept
{
  // Return if FLARM icons on moving map are disabled
  if (!GetMapSettings().show_flarm_on_map)
    return;

  // Return if FLARM data is not available
  const TrafficList &flarm = Basic().flarm.traffic;
  if (flarm.IsEmpty())
    return;

  const WindowProjection &projection = render_projection;

  // if zoomed in too far out, dont draw traffic since it will be too close to
  // the glider and so will be meaningless (serves only to clutter, cant help
  // the pilot)
  if (projection.GetMapScale() > 7300)
    return;

  // Update the FLARM traffic history
  const TimeStamp now = Basic().clock;
  flarm_traffic_history.Expire(now);
  
  // Process FLARM targets
  const FlarmTraffic *traffic = flarm.FirstTraffic();
  while (traffic != nullptr) {
    if (traffic->location_available) {
      // Add the current position to the history
      flarm_traffic_history.AddPoint(traffic->id, traffic->location, now);
      
      // Draw the FLARM target
      if (auto p = projection.GeoToScreenIfVisible(traffic->location)) {
        // Get the trail for this target
        const auto *trail = flarm_traffic_history.GetTrail(traffic->id);
        
        // Draw the traffic icon with trail
        auto color = FlarmFriends::GetFriendColor(traffic->id);
        TrafficRenderer::Draw(canvas, traffic_look, false, *traffic,
                              traffic->track - projection.GetScreenAngle(),
                              color, *p, trail, projection);
      }
    }
    
    // Move to next traffic
    traffic = flarm.NextTraffic(traffic);
  }
}

/**
 * Draws the GliderLink traffic icons onto the given canvas
 * @param canvas Canvas for drawing
 */
void
MapWindow::DrawGLinkTraffic([[maybe_unused]] Canvas &canvas) const noexcept
{
  // Original implementation
}

/**
 * Draws the teammate icon to the given canvas
 * @param canvas Canvas for drawing
 */
void
MapWindow::DrawTeammate(Canvas &canvas) const noexcept
{
  const TeamInfo &teamcode_info = Calculated();

  if (teamcode_info.teammate_available) {
    if (auto p = render_projection.GeoToScreenIfVisible(teamcode_info.teammate_location))
      traffic_look.teammate_icon.Draw(canvas, *p);
  }
}

#ifdef HAVE_SKYLINES_TRACKING
void
MapWindow::DrawSkyLinesTraffic([[maybe_unused]] Canvas &canvas) const noexcept
{
  // Original implementation
}
#endif
