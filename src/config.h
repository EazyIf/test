#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <cmath>
#include <algorithm>

enum class VehicleType { CAR, PLANE, SAILBOAT };

struct Config {
    VehicleType vehicle = VehicleType::CAR;
    int NX = 240, NY = 100, NZ = 80;
    int num_timesteps = 20000;
    double Re = 500.0;
    double angle_of_attack_deg = 0.0;
    double yaw_deg = 0.0;
    double U_LB = 0.05;
    int vtk_interval = 2000;
    int log_interval = 100;
    std::string output_dir = "./output";

    // Vehicle-specific
    double slant_angle_deg = 25.0;  // Ahmed body
    double sail_angle_deg = 15.0;   // Sailboat

    // Derived (computed by compute_derived())
    double nu_LB = 0.0;
    double tau = 0.0;
    int L_ref = 0;
    double A_ref = 0.0;

    void compute_derived() {
        // L_ref scales proportionally with grid size
        switch (vehicle) {
            case VehicleType::CAR:      L_ref = NX / 3; break;
            case VehicleType::PLANE:    L_ref = NX * 3 / 10; break;
            case VehicleType::SAILBOAT: L_ref = NX / 3; break;
        }

        // Auto-set U_LB: balance stability (tau >= 0.52) and compressibility (Ma < 0.1)
        // tau = 3*U*L/Re + 0.5
        // Target tau = 0.52 → U = 0.02*Re/(3*L)
        // Ma limit: U < 0.057
        double U_for_tau = 0.02 * Re / (3.0 * L_ref);
        double max_U_mach = 0.057;
        U_LB = std::min(U_for_tau, max_U_mach);
        if (U_LB < 0.005) U_LB = 0.005;

        nu_LB = U_LB * L_ref / Re;
        tau = 3.0 * nu_LB + 0.5;
    }

    void set_grid_for_vehicle(const std::string& resolution) {
        double scale = 1.0;
        if (resolution == "coarse") scale = 0.5;
        else if (resolution == "medium") scale = 0.75;
        else scale = 1.0; // fine

        switch (vehicle) {
            case VehicleType::CAR:
                NX = (int)(240 * scale); NY = (int)(100 * scale); NZ = (int)(80 * scale);
                break;
            case VehicleType::PLANE:
                NX = (int)(200 * scale); NY = (int)(120 * scale); NZ = (int)(100 * scale);
                break;
            case VehicleType::SAILBOAT:
                NX = (int)(180 * scale); NY = (int)(100 * scale); NZ = (int)(80 * scale);
                break;
        }
    }

    double memory_estimate_mb() const {
        size_t N = (size_t)NX * NY * NZ;
        // 2 distribution arrays (19 doubles each) + rho + ux,uy,uz + cell_type
        size_t bytes = 2 * 19 * N * sizeof(double) + 4 * N * sizeof(double) + N;
        return bytes / (1024.0 * 1024.0);
    }
};

#endif
