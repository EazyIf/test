#ifndef GEOMETRY_CAR_H
#define GEOMETRY_CAR_H

#include "geometry.h"
#include <cmath>
#include <algorithm>

class AhmedBody : public Geometry {
public:
    double L, W, H, G, slant_angle, slant_len;

    AhmedBody(const Config& cfg) {
        L = cfg.L_ref;
        W = 0.389 * L;
        H = 0.288 * L;
        G = 0.05 * L;
        slant_angle = cfg.slant_angle_deg * M_PI / 180.0;
        slant_len = 0.222 * L;

        // Center body in domain
        cx = cfg.NX * 0.35;
        cy = cfg.NY * 0.5;
        cz = cfg.NZ * 0.35;  // center body at 35% of domain height
    }

    double sdf(double x, double y, double z) const override {
        // Body centered at origin in body frame
        // Body spans: x in [-L/2, L/2], y in [-W/2, W/2], z in [-H/2, H/2]
        // Ground clearance shifts z: body bottom at z_body = -H/2, actual z = z_body + G + H/2

        double dx = std::abs(x) - L / 2.0;
        double dy = std::abs(y) - W / 2.0;
        double dz_low = -(z + H / 2.0);     // below body
        double dz_high = z - H / 2.0;        // above body

        double d_box = std::max({dx, dy, dz_low, dz_high});

        // Rear slant: for x > L/2 - slant_len, top surface drops
        double x_slant_start = L / 2.0 - slant_len;
        if (x > x_slant_start) {
            double z_slant = H / 2.0 - std::tan(slant_angle) * (x - x_slant_start);
            double d_slant = z - z_slant;
            d_box = std::max(d_box, d_slant);
        }

        // Rounded front: smooth the leading edge with a cylindrical cap
        double front_radius = 0.1 * L;
        double x_round_start = -L / 2.0 + front_radius;
        if (x < x_round_start) {
            double r_xz = std::sqrt((x - x_round_start) * (x - x_round_start) +
                                     std::max(0.0, std::abs(z) - (H / 2.0 - front_radius)) *
                                     std::max(0.0, std::abs(z) - (H / 2.0 - front_radius)));
            double d_round = r_xz - front_radius;
            d_box = std::max(d_box, d_round);
        }

        // Ground plane: everything below z = -H/2 (after centering, ground is at z=-H/2 in body frame)
        // This is handled naturally by the box bottom

        return d_box;
    }

    double reference_area(const Config& /*cfg*/) const override {
        // Frontal area: W * H
        return W * H;
    }

    std::string name() const override { return "Ahmed Body (Car)"; }
};

#endif
