#include "config.h"
#include "lattice.h"
#include "solver.h"
#include "boundary.h"
#include "geometry.h"
#include "forces.h"
#include "vtk_writer.h"
#include "csv_logger.h"
#include "tui.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include <sys/stat.h>

static double read_double(const std::string& prompt, double default_val) {
    printf("%s [%.4g]: ", prompt.c_str(), default_val);
    std::string line;
    std::getline(std::cin, line);
    if (line.empty()) return default_val;
    return std::atof(line.c_str());
}

static int read_int(const std::string& prompt, int default_val) {
    printf("%s [%d]: ", prompt.c_str(), default_val);
    std::string line;
    std::getline(std::cin, line);
    if (line.empty()) return default_val;
    return std::atoi(line.c_str());
}

static std::string read_string(const std::string& prompt, const std::string& default_val) {
    printf("%s [%s]: ", prompt.c_str(), default_val.c_str());
    std::string line;
    std::getline(std::cin, line);
    if (line.empty()) return default_val;
    return line;
}

static Config interactive_menu() {
    Config cfg;

    printf("\n");
    printf("========================================\n");
    printf("   3D Aerodynamics Simulator (LBM)\n");
    printf("   Lattice Boltzmann D3Q19 + TRT\n");
    printf("========================================\n\n");

    printf("Select vehicle type:\n");
    printf("  [1] Car (Ahmed Body)\n");
    printf("  [2] Plane (Aircraft)\n");
    printf("  [3] Sailboat\n\n");

    int choice = read_int("Enter choice", 1);
    switch (choice) {
        case 1: cfg.vehicle = VehicleType::CAR; break;
        case 2: cfg.vehicle = VehicleType::PLANE; break;
        case 3: cfg.vehicle = VehicleType::SAILBOAT; break;
        default:
            printf("Invalid choice, defaulting to Car.\n");
            cfg.vehicle = VehicleType::CAR;
    }

    printf("\n--- Initial Configuration ---\n");
    cfg.Re = read_double("Reynolds number", 150.0);

    switch (cfg.vehicle) {
        case VehicleType::CAR:
            cfg.yaw_deg = read_double("Yaw angle (degrees)", 0.0);
            cfg.slant_angle_deg = read_double("Rear slant angle (degrees)", 25.0);
            break;
        case VehicleType::PLANE:
            cfg.angle_of_attack_deg = read_double("Angle of attack (degrees)", 5.0);
            cfg.yaw_deg = read_double("Yaw angle (degrees)", 0.0);
            break;
        case VehicleType::SAILBOAT:
            cfg.yaw_deg = read_double("Wind yaw angle (degrees)", 0.0);
            cfg.sail_angle_deg = read_double("Sail angle (degrees)", 15.0);
            break;
    }

    printf("\n--- Grid Resolution ---\n");
    printf("  [1] Coarse  (fast, best for interactive)\n");
    printf("  [2] Medium  (balanced)\n");
    printf("  [3] Fine    (high fidelity, slower)\n");
    int res_choice = read_int("Resolution", 1);
    std::string resolution = "coarse";
    if (res_choice == 2) resolution = "medium";
    else if (res_choice == 3) resolution = "fine";
    cfg.set_grid_for_vehicle(resolution);
    cfg.output_dir = read_string("Output directory", "./output");

    cfg.compute_derived();

    printf("\n  Grid: %d x %d x %d (%d nodes, ~%.0f MB)\n",
           cfg.NX, cfg.NY, cfg.NZ, cfg.NX * cfg.NY * cfg.NZ, cfg.memory_estimate_mb());
    printf("  tau = %.4f, Ma = %.4f\n\n", cfg.tau, cfg.U_LB / 0.5774);

    std::string proceed = read_string("Launch interactive simulation? (Y/n)", "Y");
    if (proceed == "n" || proceed == "N") {
        printf("Aborted.\n");
        std::exit(0);
    }

    return cfg;
}

static void do_revoxelize(Config& cfg, Lattice& lat, std::unique_ptr<Geometry>& geom,
                          ForceComputer& fc) {
    std::fill(lat.cell.begin(), lat.cell.end(), CellType::FLUID);
    geom = create_geometry(cfg);
    geom->voxelize(lat, cfg);
    cfg.A_ref = geom->reference_area(cfg);
    lat.initialize(1.0, cfg.U_LB, 0.0, 0.0);
    fc.find_surface_cells(lat);
}

int main() {
    Config cfg = interactive_menu();

    mkdir(cfg.output_dir.c_str(), 0755);

    printf("Allocating lattice...\n");
    Lattice lat;
    lat.allocate(cfg.NX, cfg.NY, cfg.NZ);

    printf("Creating geometry...\n");
    auto geom = create_geometry(cfg);
    geom->voxelize(lat, cfg);
    cfg.A_ref = geom->reference_area(cfg);

    printf("Initializing flow field...\n");
    lat.initialize(1.0, cfg.U_LB, 0.0, 0.0);

    Solver solver;
    BoundaryConditions bc;
    ForceComputer fc;
    fc.find_surface_cells(lat);
    VTKWriter vtk;
    CSVLogger logger;
    logger.open(cfg.output_dir + "/forces.csv");

    // Residual tracking
    std::vector<double> ux_old(lat.ux.begin(), lat.ux.end());
    std::vector<double> uy_old(lat.uy.begin(), lat.uy.end());
    std::vector<double> uz_old(lat.uz.begin(), lat.uz.end());

    // ASCII flow viz buffer
    const int VIZ_COLS = 60, VIZ_ROWS = 20;
    std::vector<char> flow_buf(VIZ_COLS * VIZ_ROWS, ' ');

    TUIState tui;
    auto params = build_param_list(cfg);

    // Enter raw terminal mode
    printf("\033[2J");  // clear screen
    RawMode raw_mode;

    // Initial render
    render_flow_slice(lat, cfg.NY / 2, VIZ_COLS, VIZ_ROWS, flow_buf.data(), cfg.U_LB * 2.0);
    draw_screen(tui, cfg, params, flow_buf.data(), VIZ_COLS, VIZ_ROWS);

    int force_interval = 20;  // compute forces every N steps
    int render_interval = 5;  // redraw screen every N steps

    while (true) {
        // --- 1. Poll keyboard input ---
        KeyAction key = poll_key();

        switch (key) {
            case KeyAction::QUIT:
                goto cleanup;

            case KeyAction::PAUSE:
                tui.paused = !tui.paused;
                tui.status_msg = tui.paused ? "Paused. Press SPACE to resume." : "";
                break;

            case KeyAction::UP:
                tui.selected_param = (tui.selected_param - 1 + (int)params.size()) % (int)params.size();
                break;

            case KeyAction::DOWN:
                tui.selected_param = (tui.selected_param + 1) % (int)params.size();
                break;

            case KeyAction::LEFT:
            case KeyAction::RIGHT: {
                auto& p = params[tui.selected_param];
                double delta = (key == KeyAction::RIGHT) ? p.step_size : -p.step_size;
                *p.value = std::clamp(*p.value + delta, p.min_val, p.max_val);

                if (p.requires_revoxelize) {
                    tui.status_msg = "Re-voxelizing geometry...";
                    draw_screen(tui, cfg, params, flow_buf.data(), VIZ_COLS, VIZ_ROWS);
                    cfg.compute_derived();
                    do_revoxelize(cfg, lat, geom, fc);
                    ux_old.assign(lat.ux.begin(), lat.ux.end());
                    uy_old.assign(lat.uy.begin(), lat.uy.end());
                    uz_old.assign(lat.uz.begin(), lat.uz.end());
                    tui.step = 0;
                    tui.Cd = tui.Cl = tui.Cs = 0;
                    tui.residual = 0;
                    tui.unstable = false;
                    tui.status_msg = "";
                } else {
                    // Re-only: just recompute tau
                    cfg.compute_derived();
                    tui.status_msg = "";
                }
                break;
            }

            case KeyAction::RESET:
                lat.initialize(1.0, cfg.U_LB, 0.0, 0.0);
                ux_old.assign(lat.ux.begin(), lat.ux.end());
                uy_old.assign(lat.uy.begin(), lat.uy.end());
                uz_old.assign(lat.uz.begin(), lat.uz.end());
                tui.step = 0;
                tui.Cd = tui.Cl = tui.Cs = 0;
                tui.residual = 0;
                tui.unstable = false;
                tui.status_msg = "Flow field reset.";
                break;

            case KeyAction::VEHICLE: {
                int v = ((int)cfg.vehicle + 1) % 3;
                cfg.vehicle = (VehicleType)v;
                cfg.compute_derived();
                // Reset vehicle-specific defaults
                switch (cfg.vehicle) {
                    case VehicleType::CAR:
                        cfg.angle_of_attack_deg = 0; cfg.slant_angle_deg = 25; break;
                    case VehicleType::PLANE:
                        cfg.angle_of_attack_deg = 5; cfg.yaw_deg = 0; break;
                    case VehicleType::SAILBOAT:
                        cfg.sail_angle_deg = 15; cfg.yaw_deg = 0; break;
                }
                do_revoxelize(cfg, lat, geom, fc);
                ux_old.assign(lat.ux.begin(), lat.ux.end());
                uy_old.assign(lat.uy.begin(), lat.uy.end());
                uz_old.assign(lat.uz.begin(), lat.uz.end());
                params = build_param_list(cfg);
                tui.selected_param = 0;
                tui.step = 0;
                tui.Cd = tui.Cl = tui.Cs = 0;
                tui.residual = 0;
                tui.unstable = false;
                tui.status_msg = "Switched vehicle.";
                break;
            }

            case KeyAction::SPEED_UP:
                tui.render_skip = std::min(tui.render_skip * 2, 64);
                break;

            case KeyAction::SPEED_DOWN:
                tui.render_skip = std::max(tui.render_skip / 2, 1);
                break;

            default: break;
        }

        // --- 2. Run simulation step (unless paused or unstable) ---
        if (!tui.paused && !tui.unstable) {
            solver.collide(lat, cfg.tau);
            solver.stream(lat);
            bc.apply_all(lat, cfg);
            solver.compute_macroscopic(lat);
            tui.step++;

            // Compute forces periodically
            if (tui.step % force_interval == 0) {
                ForceResult fr = fc.compute(lat, cfg.tau, cfg.U_LB, cfg.A_ref);
                tui.Cd = fr.Cd;
                tui.Cl = fr.Cl;
                tui.Cs = fr.Cs;
                logger.log(tui.step, fr);

                // Compute residual + max velocity
                double sum_diff = 0.0, sum_norm = 0.0;
                double max_v2 = 0.0;
                size_t N = (size_t)cfg.NX * cfg.NY * cfg.NZ;
                for (size_t n = 0; n < N; ++n) {
                    if (lat.cell[n] == CellType::SOLID) continue;
                    double dux = lat.ux[n] - ux_old[n];
                    double duy = lat.uy[n] - uy_old[n];
                    double duz = lat.uz[n] - uz_old[n];
                    sum_diff += dux * dux + duy * duy + duz * duz;
                    sum_norm += lat.ux[n] * lat.ux[n] + lat.uy[n] * lat.uy[n] + lat.uz[n] * lat.uz[n];
                    double v2 = lat.ux[n]*lat.ux[n] + lat.uy[n]*lat.uy[n] + lat.uz[n]*lat.uz[n];
                    if (v2 > max_v2) max_v2 = v2;
                }
                tui.max_vel = std::sqrt(max_v2);
                tui.residual = std::sqrt(sum_diff / (sum_norm + 1e-30));

                ux_old.assign(lat.ux.begin(), lat.ux.end());
                uy_old.assign(lat.uy.begin(), lat.uy.end());
                uz_old.assign(lat.uz.begin(), lat.uz.end());

                if (tui.max_vel > 0.3 || std::isnan(tui.max_vel)) {
                    tui.unstable = true;
                    tui.status_msg = "UNSTABLE! Lower Re or press [r] to reset.";
                }
            }
        } else if (tui.paused) {
            usleep(16000);  // save CPU when paused
        }

        // --- 3. Render screen periodically ---
        if (tui.step % (render_interval * tui.render_skip) == 0 || tui.paused || key != KeyAction::NONE) {
            render_flow_slice(lat, cfg.NY / 2, VIZ_COLS, VIZ_ROWS, flow_buf.data(), cfg.U_LB * 2.0);
            draw_screen(tui, cfg, params, flow_buf.data(), VIZ_COLS, VIZ_ROWS);
        }

        // --- 4. Write VTK periodically ---
        if (!tui.paused && tui.step > 0 && tui.step % cfg.vtk_interval == 0) {
            std::string suffix = std::to_string(tui.step);
            vtk.write_slice_xz(lat, cfg.NY / 2, cfg.output_dir + "/slice_xz_" + suffix + ".vtk");
        }
    }

cleanup:
    // Write final VTK
    vtk.write_slice_xz(lat, cfg.NY / 2, cfg.output_dir + "/slice_xz_final.vtk");
    vtk.write_surface(lat, cfg.output_dir + "/surface_final.vtk");
    logger.close();

    // Terminal is restored by RawMode destructor
    printf("\033[?25h\033[2J\033[H");
    printf("Simulation ended at step %d.\n", tui.step);
    printf("Final Cd = %.4f, Cl = %.4f, Cs = %.4f\n", tui.Cd, tui.Cl, tui.Cs);
    printf("Output written to: %s/\n", cfg.output_dir.c_str());

    return 0;
}
