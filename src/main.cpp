#include "config.h"
#include "lattice.h"
#include "solver.h"
#include "boundary.h"
#include "geometry.h"
#include "forces.h"
#include "vtk_writer.h"
#include "csv_logger.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <string>
#include <iostream>
#include <vector>
#include <numeric>
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

    printf("\n--- Vehicle Configuration ---\n");
    cfg.Re = read_double("Reynolds number", 500.0);

    switch (cfg.vehicle) {
        case VehicleType::CAR:
            cfg.angle_of_attack_deg = 0.0;
            cfg.yaw_deg = read_double("Yaw angle (degrees)", 0.0);
            cfg.slant_angle_deg = read_double("Rear slant angle (degrees)", 25.0);
            break;
        case VehicleType::PLANE:
            cfg.angle_of_attack_deg = read_double("Angle of attack (degrees)", 5.0);
            cfg.yaw_deg = read_double("Yaw angle (degrees)", 0.0);
            break;
        case VehicleType::SAILBOAT:
            cfg.angle_of_attack_deg = 0.0;
            cfg.yaw_deg = read_double("Wind yaw angle (degrees)", 0.0);
            cfg.sail_angle_deg = read_double("Sail angle (degrees)", 15.0);
            break;
    }

    printf("\n--- Grid Resolution ---\n");
    printf("  [1] Coarse  (fast, lower accuracy)\n");
    printf("  [2] Medium  (balanced)\n");
    printf("  [3] Fine    (high fidelity, slower)\n");
    int res_choice = read_int("Resolution", 3);
    std::string resolution = "fine";
    if (res_choice == 1) resolution = "coarse";
    else if (res_choice == 2) resolution = "medium";
    cfg.set_grid_for_vehicle(resolution);

    printf("\n--- Simulation Settings ---\n");
    cfg.num_timesteps = read_int("Number of timesteps", 20000);
    cfg.vtk_interval = read_int("VTK output interval", 2000);
    cfg.log_interval = read_int("Force log interval", 100);
    cfg.output_dir = read_string("Output directory", "./output");

    // Compute derived quantities
    cfg.compute_derived();

    // Stability check
    if (cfg.tau < 0.505) {
        printf("\n  WARNING: tau = %.4f is very close to 0.5 (stability limit).\n", cfg.tau);
        printf("  Consider reducing Reynolds number or increasing resolution.\n");
        double max_re = cfg.U_LB * cfg.L_ref / ((0.505 - 0.5) / 3.0);
        printf("  Maximum stable Re for this resolution: ~%.0f\n", max_re);
    }

    double mach = cfg.U_LB / std::sqrt(1.0 / 3.0);

    printf("\n========================================\n");
    printf("   Configuration Summary\n");
    printf("========================================\n");

    const char* vname = "";
    switch (cfg.vehicle) {
        case VehicleType::CAR: vname = "Ahmed Body (Car)"; break;
        case VehicleType::PLANE: vname = "Aircraft (Plane)"; break;
        case VehicleType::SAILBOAT: vname = "Sailboat"; break;
    }

    printf("  Vehicle:       %s\n", vname);
    printf("  Grid:          %d x %d x %d (%d nodes)\n", cfg.NX, cfg.NY, cfg.NZ, cfg.NX * cfg.NY * cfg.NZ);
    printf("  Memory:        ~%.0f MB\n", cfg.memory_estimate_mb());
    printf("  Reynolds:      %.0f\n", cfg.Re);
    printf("  tau:           %.6f\n", cfg.tau);
    printf("  nu (lattice):  %.6e\n", cfg.nu_LB);
    printf("  Mach number:   %.4f\n", mach);
    printf("  L_ref:         %d lattice units\n", cfg.L_ref);
    printf("  Timesteps:     %d\n", cfg.num_timesteps);
    printf("  VTK interval:  %d\n", cfg.vtk_interval);
    printf("  Log interval:  %d\n", cfg.log_interval);
    printf("  Output dir:    %s\n", cfg.output_dir.c_str());
    printf("========================================\n\n");

    std::string proceed = read_string("Proceed? (Y/n)", "Y");
    if (proceed == "n" || proceed == "N") {
        printf("Aborted.\n");
        std::exit(0);
    }

    return cfg;
}

int main() {
    Config cfg = interactive_menu();

    // Create output directory
    mkdir(cfg.output_dir.c_str(), 0755);

    printf("\n[1/5] Allocating lattice (%d x %d x %d)...\n", cfg.NX, cfg.NY, cfg.NZ);
    Lattice lat;
    lat.allocate(cfg.NX, cfg.NY, cfg.NZ);

    printf("[2/5] Creating geometry and voxelizing...\n");
    auto geom = create_geometry(cfg);
    geom->voxelize(lat, cfg);
    cfg.A_ref = geom->reference_area(cfg);
    printf("  Reference area: %.2f lattice units^2\n", cfg.A_ref);

    printf("[3/5] Initializing flow field...\n");
    lat.initialize(1.0, cfg.U_LB, 0.0, 0.0);

    printf("[4/5] Setting up solver...\n");
    Solver solver;
    BoundaryConditions bc;
    ForceComputer fc;
    fc.find_surface_cells(lat);
    VTKWriter vtk;
    CSVLogger logger;
    logger.open(cfg.output_dir + "/forces.csv");

    // Storage for steady-state averaging
    std::vector<double> cd_history, cl_history, cs_history;

    printf("[5/5] Starting simulation...\n\n");
    printf("  %-8s  %-12s  %-12s  %-12s  %-12s\n", "Step", "Cd", "Cl", "Cs", "Residual");
    printf("  %-8s  %-12s  %-12s  %-12s  %-12s\n", "----", "--", "--", "--", "--------");

    // For residual computation
    std::vector<double> ux_old(lat.ux.begin(), lat.ux.end());
    std::vector<double> uy_old(lat.uy.begin(), lat.uy.end());
    std::vector<double> uz_old(lat.uz.begin(), lat.uz.end());

    clock_t t_start = clock();

    bool diverged = false;

    for (int t = 0; t < cfg.num_timesteps; ++t) {
        solver.collide(lat, cfg.tau);
        solver.stream(lat);
        bc.apply_all(lat, cfg);
        solver.compute_macroscopic(lat);

        if (t % cfg.log_interval == 0) {
            // Stress integration force computation (uses macroscopic quantities)
            ForceResult fr = fc.compute(lat, cfg.tau, cfg.U_LB, cfg.A_ref);
            logger.log(t, fr);

            cd_history.push_back(fr.Cd);
            cl_history.push_back(fr.Cl);
            cs_history.push_back(fr.Cs);

            // Compute residual and check stability
            double sum_diff = 0.0, sum_norm = 0.0;
            double max_vel = 0.0;
            size_t N = (size_t)cfg.NX * cfg.NY * cfg.NZ;
            for (size_t n = 0; n < N; ++n) {
                if (lat.cell[n] == CellType::SOLID) continue;
                double dux = lat.ux[n] - ux_old[n];
                double duy = lat.uy[n] - uy_old[n];
                double duz = lat.uz[n] - uz_old[n];
                sum_diff += dux * dux + duy * duy + duz * duz;
                sum_norm += lat.ux[n] * lat.ux[n] + lat.uy[n] * lat.uy[n] + lat.uz[n] * lat.uz[n];
                double vm = lat.ux[n]*lat.ux[n] + lat.uy[n]*lat.uy[n] + lat.uz[n]*lat.uz[n];
                if (vm > max_vel) max_vel = vm;
            }
            max_vel = std::sqrt(max_vel);
            double residual = std::sqrt(sum_diff / (sum_norm + 1e-30));

            printf("  %-8d  %-12.6f  %-12.6f  %-12.6f  %-12.4e  v_max=%.4f\n",
                   t, fr.Cd, fr.Cl, fr.Cs, residual, max_vel);

            if (max_vel > 0.3 || std::isnan(max_vel)) {
                printf("\n  !! SIMULATION UNSTABLE (max velocity = %.4f) !!\n", max_vel);
                printf("  !! Try: lower Re, coarser grid, or smaller U_LB !!\n\n");
                diverged = true;
                break;
            }

            ux_old.assign(lat.ux.begin(), lat.ux.end());
            uy_old.assign(lat.uy.begin(), lat.uy.end());
            uz_old.assign(lat.uz.begin(), lat.uz.end());
        }

        if (t > 0 && t % cfg.vtk_interval == 0) {
            std::string suffix = std::to_string(t);
            vtk.write_slice_xz(lat, cfg.NY / 2, cfg.output_dir + "/slice_xz_" + suffix + ".vtk");
            vtk.write_slice_xy(lat, cfg.NZ / 2, cfg.output_dir + "/slice_xy_" + suffix + ".vtk");
            printf("  >> Wrote VTK slices at step %d\n", t);
        }
    }

    clock_t t_end = clock();
    double elapsed = (double)(t_end - t_start) / CLOCKS_PER_SEC;

    // Final VTK outputs
    if (diverged) {
        printf("Simulation stopped early due to instability.\n");
        printf("Output files may contain invalid data.\n\n");
    }
    printf("Writing final output files...\n");
    vtk.write_3d_field(lat, cfg.output_dir + "/field_final.vtk");
    vtk.write_surface(lat, cfg.output_dir + "/surface_final.vtk");
    vtk.write_slice_xz(lat, cfg.NY / 2, cfg.output_dir + "/slice_xz_final.vtk");
    vtk.write_slice_xy(lat, cfg.NZ / 2, cfg.output_dir + "/slice_xy_final.vtk");
    logger.close();

    // Compute steady-state averages (last 25% of data)
    size_t avg_start = cd_history.size() * 3 / 4;
    size_t avg_count = cd_history.size() - avg_start;

    double cd_avg = 0, cl_avg = 0, cs_avg = 0;
    double cd_var = 0, cl_var = 0, cs_var = 0;

    if (avg_count > 0) {
        for (size_t i = avg_start; i < cd_history.size(); ++i) {
            cd_avg += cd_history[i];
            cl_avg += cl_history[i];
            cs_avg += cs_history[i];
        }
        cd_avg /= avg_count;
        cl_avg /= avg_count;
        cs_avg /= avg_count;

        for (size_t i = avg_start; i < cd_history.size(); ++i) {
            cd_var += (cd_history[i] - cd_avg) * (cd_history[i] - cd_avg);
            cl_var += (cl_history[i] - cl_avg) * (cl_history[i] - cl_avg);
            cs_var += (cs_history[i] - cs_avg) * (cs_history[i] - cs_avg);
        }
        cd_var = std::sqrt(cd_var / avg_count);
        cl_var = std::sqrt(cl_var / avg_count);
        cs_var = std::sqrt(cs_var / avg_count);
    }

    // Find max velocity
    double max_vel = 0;
    double total_mass = 0;
    size_t N = (size_t)cfg.NX * cfg.NY * cfg.NZ;
    for (size_t n = 0; n < N; ++n) {
        double vm = std::sqrt(lat.ux[n]*lat.ux[n] + lat.uy[n]*lat.uy[n] + lat.uz[n]*lat.uz[n]);
        if (vm > max_vel) max_vel = vm;
        total_mass += lat.rho[n];
    }

    int hours = (int)(elapsed / 3600);
    int mins = (int)((elapsed - hours * 3600) / 60);
    int secs = (int)(elapsed - hours * 3600 - mins * 60);

    printf("\n========================================\n");
    printf("   SIMULATION COMPLETE\n");
    printf("========================================\n");
    printf("  Vehicle:       %s\n", geom->name().c_str());
    printf("  Reynolds:      %.0f\n", cfg.Re);
    printf("  Grid:          %d x %d x %d (%d nodes)\n", cfg.NX, cfg.NY, cfg.NZ, cfg.NX * cfg.NY * cfg.NZ);
    printf("  Timesteps:     %d\n", cfg.num_timesteps);
    printf("  Runtime:       %dh %dm %ds\n", hours, mins, secs);

    printf("\n--- Force Coefficients (averaged last 25%%) ---\n");
    printf("  Drag  (Cd) = %.6f +/- %.6f\n", cd_avg, cd_var);
    printf("  Lift  (Cl) = %.6f +/- %.6f\n", cl_avg, cl_var);
    printf("  Side  (Cs) = %.6f +/- %.6f\n", cs_avg, cs_var);

    printf("\n--- Flow Statistics ---\n");
    printf("  Max velocity:   %.4f (Mach = %.4f)\n", max_vel, max_vel / std::sqrt(1.0/3.0));
    printf("  Total mass:     %.2f (expected: %.2f)\n", total_mass, (double)N);

    printf("\n--- Output Files ---\n");
    printf("  Force log:     %s/forces.csv\n", cfg.output_dir.c_str());
    printf("  VTK 3D field:  %s/field_final.vtk\n", cfg.output_dir.c_str());
    printf("  VTK surface:   %s/surface_final.vtk\n", cfg.output_dir.c_str());
    printf("  VTK slices:    %s/slice_*.vtk\n", cfg.output_dir.c_str());
    printf("========================================\n\n");

    return 0;
}
