#pragma once

#include <mapbox/geojsonvt/types.hpp>

namespace mapbox {
namespace geojsonvt {
namespace detail {

template <uint8_t I>
class clipper {
public:
    clipper(double k1_, double k2_, bool lineMetrics_ = false)
        : k1(k1_), k2(k2_), k1_metric(k1_), k2_metric(k2_), lineMetrics(lineMetrics_), z_(0), x_(0), y_(0) {}

    clipper(double k1_geom, double k2_geom, double k1_m, double k2_m, bool lineMetrics_ = false)
    : k1(k1_geom), k2(k2_geom), k1_metric(k1_m), k2_metric(k2_m), lineMetrics(lineMetrics_), z_(0), x_(0), y_(0) {}

    clipper(double k1_geom, double k2_geom, double k1_m, double k2_m, uint8_t z, uint32_t x, uint32_t y, bool lineMetrics_ = false)
    : k1(k1_geom), k2(k2_geom), k1_metric(k1_m), k2_metric(k2_m), lineMetrics(lineMetrics_), z_(z), x_(x), y_(y) {}

    const double k1;
    const double k2;
    const double k1_metric;
    const double k2_metric;
    const bool lineMetrics;
    const uint8_t z_;
    uint32_t x_, y_;

    vt_geometry operator()(const vt_empty& empty) const {
        return empty;
    }

    vt_geometry operator()(const vt_point& point) const {
        return point;
    }

    vt_geometry operator()(const vt_multi_point& points) const {
        vt_multi_point part;
        for (const auto& p : points) {
            const double ak = get<I>(p);
            if (ak >= k1 && ak <= k2)
                part.emplace_back(p);
        }
        return part;
    }

    vt_geometry operator()(const vt_line_string& line) const {
        vt_multi_line_string parts;
        clipLine(line, parts);
        if (parts.size() == 1)
            return parts[0];
        else
            return parts;
    }

    vt_geometry operator()(const vt_multi_line_string& lines) const {
        vt_multi_line_string parts;
        for (const auto& line : lines) {
            clipLine(line, parts);
        }
        if (parts.size() == 1)
            return parts[0];
        else
            return parts;
    }

    vt_geometry operator()(const vt_polygon& polygon) const {
        vt_polygon result;
        for (const auto& ring : polygon) {
            auto new_ring = clipRing(ring);
            if (!new_ring.empty())
                result.emplace_back(std::move(new_ring));
        }
        return result;
    }

    vt_geometry operator()(const vt_multi_polygon& polygons) const {
        vt_multi_polygon result;
        for (const auto& polygon : polygons) {
            vt_polygon p;
            for (const auto& ring : polygon) {
                auto new_ring = clipRing(ring);
                if (!new_ring.empty())
                    p.emplace_back(std::move(new_ring));
            }
            if (!p.empty())
                result.emplace_back(std::move(p));
        }
        return result;
    }

    vt_geometry operator()(const vt_geometry_collection& geometries) const {
        vt_geometry_collection result;
        for (const auto& geometry : geometries) {
            vt_geometry::visit(geometry,
                               [&](const auto& g) { result.emplace_back(this->operator()(g)); });
        }
        return result;
    }

private:
    vt_line_string newSlice(const vt_line_string& line) const {
        vt_line_string slice;
        slice.dist = line.dist;
        if (lineMetrics) {
            slice.segStart = line.segStart;
            slice.segEnd = line.segEnd;
        }
        return slice;
    }

    void clipLine(const vt_line_string& line, vt_multi_line_string& slices) const {
        const size_t len = line.size();
        double lineLen = line.segStart;
        double segLen = 0.0;
        double t_geom = 0.0;
        double t_metric = 0.0;


        if (len < 2)
            return;

        vt_line_string slice = newSlice(line);

        for (size_t i = 0; i < (len - 1); ++i) {
            const auto& a = line[i];
            const auto& b = line[i + 1];
            const double ak = get<I>(a);
            const double bk = get<I>(b);
            const bool isLastSeg = (i == (len - 2));

            if (lineMetrics) {
                segLen = ::hypot((b.x - a.x), (b.y - a.y));
            }

            if (ak < k1) {
                if (bk > k2) { // ---|-----|-->
                    t_geom = calc_progress<I>(a, b, k1);
                    slice.emplace_back(intersect<I>(a, b, k1, calc_progress<I>(a, b, k1)));
                    if (lineMetrics) {
                        t_metric =  calc_progress<I>(a, b, k1_metric);
                        slice.segStart = lineLen + segLen * t_metric;
                    }

                    t_geom = calc_progress<I>(a, b, k2);
                    slice.emplace_back(intersect<I>(a, b, k2, t_geom));
                    if (lineMetrics) {
                        t_metric = calc_progress<I>(a, b, k2_metric);
                        slice.segEnd = lineLen + segLen * t_metric;
                    }
                    slices.emplace_back(std::move(slice));

                    slice = newSlice(line);

                } else if (bk > k1) { // ---|-->  |
                    t_geom = calc_progress<I>(a, b, k1);
                    slice.emplace_back(intersect<I>(a, b, k1, t_geom));
                    if (lineMetrics) {
                        t_metric = calc_progress<I>(a, b, k1_metric);
                        slice.segStart = lineLen + segLen * t_metric;
                    }
                    if (isLastSeg) slice.emplace_back(b); // last point

                } else if (bk == k1 && !isLastSeg) { // --->|..  |
                    if (lineMetrics) {
                        slice.segStart = lineLen + segLen;
                    }
                    slice.emplace_back(b);
                }
            } else if (ak > k2) {
                if (bk < k1) { // <--|-----|---
                    t_geom = calc_progress<I>(a, b, k2);
                    slice.emplace_back(intersect<I>(a, b, k2, t_geom));
                    if (lineMetrics) {
                        t_metric = calc_progress<I>(a, b, k2_metric);
                        slice.segStart = lineLen + segLen * t_metric;
                    }

                    t_geom = calc_progress<I>(a, b, k1);
                    slice.emplace_back(intersect<I>(a, b, k1, t_geom));
                    if (lineMetrics) {
                        t_metric = calc_progress<I>(a, b, k1_metric);
                        slice.segEnd = lineLen + segLen * t_metric;
                    }

                    slices.emplace_back(std::move(slice));

                    slice = newSlice(line);

                } else if (bk < k2) { // |  <--|---
                    t_geom = calc_progress<I>(a, b, k2);
                    slice.emplace_back(intersect<I>(a, b, k2, t_geom));
                    if (lineMetrics) {
                        t_metric = calc_progress<I>(a, b, k2_metric);
                        slice.segStart = lineLen + segLen * t_metric;
                    }
                    if (isLastSeg) slice.emplace_back(b); // last point

                } else if (bk == k2 && !isLastSeg) { // |  ..|<---
                    if (lineMetrics) {
                        slice.segStart = lineLen + segLen;
                    }
                    slice.emplace_back(b);
                }
            } else {
                if (slice.empty() && lineMetrics) {
                    slice.segStart = lineLen;
                }
                slice.emplace_back(a);

                if (bk < k1) { // <--|---  |
                    t_geom = calc_progress<I>(a, b, k1);
                    slice.emplace_back(intersect<I>(a, b, k1, t_geom));
                    if (lineMetrics) {
                        t_metric = calc_progress<I>(a, b, k1_metric);
                        slice.segEnd = lineLen + segLen * t_metric;
                    }
                    slices.emplace_back(std::move(slice));
                    slice = newSlice(line);

                } else if (bk > k2) { // |  ---|-->
                    t_geom = calc_progress<I>(a, b, k2);
                    slice.emplace_back(intersect<I>(a, b, k2, t_geom));
                    if (lineMetrics) {
                        t_metric = calc_progress<I>(a, b, k2_metric);
                        slice.segEnd = lineLen + segLen * t_metric;
                    }
                    slices.emplace_back(std::move(slice));
                    slice = newSlice(line);

                } else if (isLastSeg) { // | --> |
                    slice.emplace_back(b);
                }
            }

            if (lineMetrics) lineLen += segLen;
        }

        if (!slice.empty()) { // add the final slice
            if (lineMetrics) slice.segEnd = lineLen;
            slices.emplace_back(std::move(slice));
        }
    }

    vt_linear_ring clipRing(const vt_linear_ring& ring) const {
        const size_t len = ring.size();
        vt_linear_ring slice;
        slice.area = ring.area;

        if (len < 2)
            return slice;

        for (size_t i = 0; i < (len - 1); ++i) {
            const auto& a = ring[i];
            const auto& b = ring[i + 1];
            const double ak = get<I>(a);
            const double bk = get<I>(b);

            if (ak < k1) {
                if (bk > k1) {
                    // ---|-->  |
                    slice.emplace_back(intersect<I>(a, b, k1, calc_progress<I>(a, b, k1)));
                    if (bk > k2)
                        // ---|-----|-->
                        slice.emplace_back(intersect<I>(a, b, k2, calc_progress<I>(a, b, k2)));
                    else if (i == len - 2)
                        slice.emplace_back(b); // last point
                }
            } else if (ak > k2) {
                if (bk < k2) { // |  <--|---
                    slice.emplace_back(intersect<I>(a, b, k2, calc_progress<I>(a, b, k2)));
                    if (bk < k1) // <--|-----|---
                        slice.emplace_back(intersect<I>(a, b, k1, calc_progress<I>(a, b, k1)));
                    else if (i == len - 2)
                        slice.emplace_back(b); // last point
                }
            } else {
                // | --> |
                slice.emplace_back(a);
                if (bk < k1)
                    // <--|---  |
                    slice.emplace_back(intersect<I>(a, b, k1, calc_progress<I>(a, b, k1)));
                else if (bk > k2)
                    // |  ---|-->
                    slice.emplace_back(intersect<I>(a, b, k2, calc_progress<I>(a, b, k2)));
            }
        }

        // close the polygon if its endpoints are not the same after clipping
        if (!slice.empty()) {
            const auto& first = slice.front();
            const auto& last = slice.back();
            if (first != last) {
                slice.emplace_back(first);
            }
        }

        return slice;
    }
};

/* clip features between two axis-parallel lines:
 *     |        |
 *  ___|___     |     /
 * /   |   \____|____/
 *     |        |
 */

    template <uint8_t I>
    inline vt_features clip(const vt_features& features,
                            const double k1_geom,
                            const double k2_geom,
                            const double k1_metric,
                            const double k2_metric,
                            const double minAll,
                            const double maxAll,
                            const bool lineMetrics,
                            const uint8_t z,
                            const uint32_t x,
                            const uint32_t y) {
        // Trivial accept for the entire feature set is only safe if lineMetrics are disabled.
        if (/*!lineMetrics &&*/ minAll >= k1_geom && maxAll < k2_geom)
            return features;

        if (maxAll < k1_geom || minAll >= k2_geom)
            return {};

        vt_features clipped;
        clipped.reserve(features.size());

        for (const auto& feature : features) {
            const auto& geom = feature.geometry;
            assert(feature.properties);
            const auto& props = feature.properties;
            const auto& id = feature.id;

            const double min = get<I>(feature.bbox.min);
            const double max = get<I>(feature.bbox.max);

            // Trivial accept for a single feature is also only safe if lineMetrics are disabled.
            if (/*!lineMetrics &&*/ min >= k1_geom && max < k2_geom) {
                clipped.emplace_back(feature);

            } else if (max < k1_geom || min >= k2_geom) { // Trivial reject is always safe.
                continue;

            } else { // Perform a detailed clip.
                const auto& clippedGeom = vt_geometry::visit(geom, clipper<I>{ k1_geom, k2_geom, k1_metric, k2_metric, z, x, y, lineMetrics });

                clippedGeom.match(
                    [&](const auto&) {
                        clipped.emplace_back(clippedGeom, props, id);
                    },
                    [&](const vt_multi_line_string& result) {
                        if (lineMetrics) {
                            for (const auto& segment : result) {
                                clipped.emplace_back(segment, props, id);
                            }
                        } else {
                            clipped.emplace_back(clippedGeom, props, id);
                        }
                    }
                );
            }
        }
        return clipped;
    }


    // Existing clip function (now calls the new one for backward compatibility)
    template <uint8_t I>
    inline vt_features clip(const vt_features& features,
                            const double k1,
                            const double k2,
                            const double minAll,
                            const double maxAll,
                            const bool lineMetrics) {
        // Pass 'k1' and 'k2' as both the geometry and metric boundaries
        return clip<I>(features, k1, k2, k1, k2, minAll, maxAll, lineMetrics, 0, 0, 0);
    }

} // namespace detail
} // namespace geojsonvt
} // namespace mapbox
