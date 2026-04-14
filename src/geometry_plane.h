#ifndef GEOMETRY_PLANE_H
#define GEOMETRY_PLANE_H

#include "geometry.h"
#include <cmath>
#include <algorithm>

class Aircraft : public Geometry {
public:
    double L, D, span, c_root, c_tip;
    double h_span, h_chord, v_height, v_chord;

    Aircraft(const Config& cfg) {
        L = cfg.L_ref;
        D = 0.1 * L;           // fuselage diameter
        span = 2.5 * L;        // wing full span
        c_root = 0.3 * L;      // wing root chord
        c_tip = 0.15 * L;      // wing tip chord

        h_span = 0.6 * L;      // horizontal tail span
        h_chord = 0.12 * L;    // horizontal tail chord
        v_height = 0.25 * L;   // vertical tail height
        v_chord = 0.15 * L;    // vertical tail chord

        cx = cfg.NX * 0.35;
        cy = cfg.NY * 0.5;
        cz = cfg.NZ * 0.5;
    }

    double sdf(double x, double y, double z) const override {
        double d_fuse = sdf_fuselage(x, y, z);
        double d_wing = sdf_wing(x, y, z);
        double d_htail = sdf_htail(x, y, z);
        double d_vtail = sdf_vtail(x, y, z);

        return std::min({d_fuse, d_wing, d_htail, d_vtail});
    }

    double reference_area(const Config& /*cfg*/) const override {
        // Wing planform area (approximate trapezoidal)
        return span * (c_root + c_tip) * 0.5;
    }

    std::string name() const override { return "Aircraft (Plane)"; }

private:
    double sdf_fuselage(double x, double y, double z) const {
        double a = L / 2.0;   // semi-axis x
        double b = D / 2.0;   // semi-axis y, z
        double val = std::sqrt((x / a) * (x / a) + (y / b) * (y / b) + (z / b) * (z / b));
        return (val - 1.0) * b;  // scale to approximate Euclidean distance
    }

    double sdf_wing(double x, double y, double z) const {
        double half_span = span / 2.0;
        double abs_y = std::abs(y);

        if (abs_y < D / 2.0 || abs_y > half_span) return 1e6;

        // Taper: chord varies linearly from root to tip
        double t = (abs_y - D / 2.0) / (half_span - D / 2.0);
        double chord = c_root * (1.0 - t) + c_tip * t;

        // Wing leading edge: slight sweep (25 degrees)
        double sweep = 25.0 * M_PI / 180.0;
        double x_le = -0.05 * L + std::tan(sweep) * (abs_y - D / 2.0);

        double xn = (x - x_le) / chord;  // normalized chord position [0, 1]

        if (xn < 0.0 || xn > 1.0) return 1e6;

        // NACA 0012 half-thickness
        double thick = 0.12 * chord;
        double ht = 0.5 * thick * (2.98 * std::sqrt(xn) - 1.26 * xn
                     - 3.516 * xn * xn + 2.843 * xn * xn * xn - 1.015 * xn * xn * xn * xn);

        return std::abs(z) - ht;
    }

    double sdf_htail(double x, double y, double z) const {
        double x_tail = 0.42 * L;  // position behind center
        double abs_y = std::abs(y);

        if (abs_y > h_span / 2.0) return 1e6;

        double xn = (x - x_tail) / h_chord;
        if (xn < 0.0 || xn > 1.0) return 1e6;

        double thick = 0.08 * h_chord;
        double ht = 0.5 * thick * std::sqrt(4.0 * xn * (1.0 - xn));

        return std::abs(z) - ht;
    }

    double sdf_vtail(double x, double y, double z) const {
        double x_tail = 0.42 * L;

        if (z < 0.0 || z > v_height) return 1e6;

        double xn = (x - x_tail) / v_chord;
        if (xn < 0.0 || xn > 1.0) return 1e6;

        double thick = 0.08 * v_chord;
        double ht = 0.5 * thick * std::sqrt(4.0 * xn * (1.0 - xn));

        return std::abs(y) - ht;
    }
};

#endif
