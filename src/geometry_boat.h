#ifndef GEOMETRY_BOAT_H
#define GEOMETRY_BOAT_H

#include "geometry.h"
#include <cmath>
#include <algorithm>

class Sailboat : public Geometry {
public:
    double L, B, T, mast_h, mast_r;
    double sail_chord, sail_low, sail_high, sail_camber, sail_thickness;
    double sail_angle;

    Sailboat(const Config& cfg) {
        L = cfg.L_ref;
        B = 0.25 * L;           // beam (width)
        T = 0.1 * L;            // draft (depth below waterline)
        mast_h = 0.8 * L;       // mast height
        mast_r = 0.015 * L;     // mast radius
        sail_chord = 0.3 * L;
        sail_low = 0.1 * L;
        sail_high = 0.75 * L;
        sail_camber = 0.08 * sail_chord;
        sail_thickness = 0.01 * L;
        sail_angle = cfg.sail_angle_deg * M_PI / 180.0;

        // Waterline at z = 0.3*NZ, body centered in x and y
        cx = cfg.NX * 0.35;
        cy = cfg.NY * 0.5;
        cz = cfg.NZ * 0.3;  // waterline position
    }

    double sdf(double x, double y, double z) const override {
        double d_hull = sdf_hull(x, y, z);
        double d_mast = sdf_mast(x, y, z);
        double d_sail = sdf_sail(x, y, z);

        return std::min({d_hull, d_mast, d_sail});
    }

    double reference_area(const Config& /*cfg*/) const override {
        // Sail area (projected)
        return sail_chord * (sail_high - sail_low);
    }

    std::string name() const override { return "Sailboat"; }

private:
    double sdf_hull(double x, double y, double z) const {
        // Half-ellipsoid below waterline (z <= 0 in body frame)
        double a = L / 2.0;
        double b = B / 2.0;
        double c = T;

        if (z > 0.0) {
            // Above waterline: thin deck
            double deck_height = 0.03 * L;
            if (z > deck_height) return 1e6;

            // Deck is same ellipse footprint as hull at waterline
            double r2 = (x / a) * (x / a) + (y / b) * (y / b);
            return (std::sqrt(r2) - 1.0) * std::min(a, b);
        }

        // Below waterline: half-ellipsoid
        double val = std::sqrt((x / a) * (x / a) + (y / b) * (y / b) + (z / c) * (z / c));
        return (val - 1.0) * std::min({a, b, c});
    }

    double sdf_mast(double x, double y, double z) const {
        double x_m = -0.1 * L;  // slightly forward of center

        if (z < 0.0 || z > mast_h) return 1e6;

        double r = std::sqrt((x - x_m) * (x - x_m) + y * y);
        return r - mast_r;
    }

    double sdf_sail(double x, double y, double z) const {
        if (z < sail_low || z > sail_high) return 1e6;

        double x_m = -0.1 * L;

        // Sail is attached to mast at (x_m, 0)
        // Rotate sail by sail_angle about z-axis from mast position
        double dx = x - x_m;
        double dy = y;

        // Transform into sail-local coords (along sail chord)
        double along = dx * std::cos(sail_angle) + dy * std::sin(sail_angle);
        double perp = -dx * std::sin(sail_angle) + dy * std::cos(sail_angle);

        if (along < 0.0 || along > sail_chord) return 1e6;

        // Parabolic camber
        double t_norm = along / sail_chord;
        double camber = 4.0 * sail_camber * t_norm * (1.0 - t_norm);

        return std::abs(perp - camber) - sail_thickness / 2.0;
    }
};

#endif
