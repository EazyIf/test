#include "config.h"
#include "lattice.h"
#include "solver.h"
#include "boundary.h"
#include "geometry.h"
#include "forces.h"
#include "vtk_writer.h"
#include "csv_logger.h"
#include "webserver.h"
#include "web_ui.h"


#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <memory>
#include <sys/stat.h>
#include <unistd.h>
#include <mutex>
#include <atomic>

// --- Shared simulation state for web endpoints ---
struct SimState {
    std::atomic<bool> paused{false};
    std::atomic<bool> unstable{false};
    int step = 0;
    double Cd = 0, Cl = 0, Cs = 0;
    double residual = 0, max_vel = 0;

    // Pending commands from web UI
    std::atomic<bool> cmd_reset{false};
    std::atomic<bool> cmd_vehicle{false};
    std::atomic<bool> cmd_pause{false};
    std::atomic<bool> param_changed{false};
    std::atomic<bool> geom_changed{false};

    double new_re = 0;
    double new_yaw = 0;
    double new_aoa = 0;
    double new_slant = 0;
    double new_sail = 0;
};

static SimState g_sim;

// --- JSON helpers ---
static std::string json_state(const Config& cfg, const Lattice& lat) {
    std::ostringstream ss;
    ss << std::fixed;
    ss << "{";
    ss << "\"step\":" << g_sim.step;
    ss << ",\"cd\":" << g_sim.Cd;
    ss << ",\"cl\":" << g_sim.Cl;
    ss << ",\"cs\":" << g_sim.Cs;
    ss << ",\"residual\":" << g_sim.residual;
    ss << ",\"max_vel\":" << g_sim.max_vel;
    ss << ",\"tau\":" << cfg.tau;
    ss << ",\"re\":" << cfg.Re;
    ss << ",\"nx\":" << cfg.NX << ",\"ny\":" << cfg.NY << ",\"nz\":" << cfg.NZ;
    ss << ",\"paused\":" << (g_sim.paused ? "true" : "false");
    ss << ",\"unstable\":" << (g_sim.unstable ? "true" : "false");

    const char* vname = "car";
    switch (cfg.vehicle) {
        case VehicleType::CAR: vname = "car"; break;
        case VehicleType::PLANE: vname = "plane"; break;
        case VehicleType::SAILBOAT: vname = "sailboat"; break;
    }
    ss << ",\"vehicle\":\"" << vname << "\"";

    // Center of geometry
    ss << ",\"cx\":" << cfg.NX * 0.35;
    ss << ",\"cy\":" << cfg.NY * 0.5;
    ss << ",\"cz\":" << cfg.NZ * 0.35;

    // Sampled flow field (stride for performance)
    int stride = std::max(3, cfg.NX / 30);
    ss << ",\"flow\":[";
    bool first = true;
    for (int x = stride/2; x < cfg.NX; x += stride) {
        for (int y = cfg.NY/2; y <= cfg.NY/2; y += stride) { // single y-slice
            for (int z = stride/2; z < cfg.NZ; z += stride) {
                size_t n = lat.idx(x, y, z);
                if (lat.cell[n] == CellType::SOLID) continue;
                double vmag = std::sqrt(lat.ux[n]*lat.ux[n]+lat.uy[n]*lat.uy[n]+lat.uz[n]*lat.uz[n]);
                if (vmag < cfg.U_LB * 0.01) continue;
                if (!first) ss << ",";
                ss << "[" << x << "," << y << "," << z << ","
                   << lat.ux[n] << "," << lat.uy[n] << "," << lat.uz[n] << "]";
                first = false;
            }
        }
    }
    ss << "]";

    ss << "}";
    return ss.str();
}

static std::string json_mesh(const Config& cfg, const Lattice& lat) {
    std::ostringstream ss;
    ss << std::fixed;
    ss << "{\"cx\":" << cfg.NX * 0.35
       << ",\"cy\":" << cfg.NY * 0.5
       << ",\"cz\":" << cfg.NZ * 0.35;
    ss << ",\"body\":[";

    bool first = true;
    for (int x = 0; x < lat.NX; ++x) {
        for (int y = 0; y < lat.NY; ++y) {
            for (int z = 0; z < lat.NZ; ++z) {
                size_t n = lat.idx(x, y, z);
                if (lat.cell[n] != CellType::FLUID) continue;

                // Check if adjacent to solid (surface cell)
                bool surface = false;
                double rnx = 0, rny = 0, rnz = 0;
                for (int i = 1; i <= 6; ++i) {
                    int xn = x+ex[i], yn = y+ey[i], zn = z+ez[i];
                    if (xn < 0 || xn >= lat.NX || yn < 0 || yn >= lat.NY || zn < 0 || zn >= lat.NZ) continue;
                    if (lat.cell[lat.idx(xn,yn,zn)] == CellType::SOLID) {
                        surface = true;
                        rnx += ex[i]; rny += ey[i]; rnz += ez[i];
                    }
                }
                if (!surface) continue;

                double mag = std::sqrt(rnx*rnx + rny*rny + rnz*rnz);
                if (mag > 0) { rnx /= mag; rny /= mag; rnz /= mag; }
                double p = lat.rho[n] / 3.0;

                if (!first) ss << ",";
                ss << "[" << x << "," << y << "," << z << ","
                   << rnx << "," << rny << "," << rnz << "," << p << "]";
                first = false;
            }
        }
    }

    ss << "]}";
    return ss.str();
}

static void handle_control(const std::string& body, Config& cfg) {
    // Minimal JSON parsing for control commands
    if (body.find("\"reset\"") != std::string::npos) {
        g_sim.cmd_reset = true;
    }
    if (body.find("\"pause\"") != std::string::npos) {
        g_sim.cmd_pause = true;
    }
    if (body.find("\"vehicle\"") != std::string::npos) {
        g_sim.cmd_vehicle = true;
    }
    if (body.find("\"param\"") != std::string::npos) {
        // Extract name and value
        // Find "value": number
        size_t vpos = body.find("\"value\":");
        size_t npos = body.find("\"name\":\"");
        if (vpos != std::string::npos && npos != std::string::npos) {
            double val = std::atof(body.c_str() + vpos + 8);
            std::string name = body.substr(npos + 8);
            name = name.substr(0, name.find('"'));

            if (name == "re") {
                cfg.Re = val;
                cfg.compute_derived();
                g_sim.param_changed = true;
            } else if (name == "yaw") {
                cfg.yaw_deg = val;
                g_sim.geom_changed = true;
            } else if (name == "aoa") {
                cfg.angle_of_attack_deg = val;
                g_sim.geom_changed = true;
            } else if (name == "slant") {
                cfg.slant_angle_deg = val;
                g_sim.geom_changed = true;
            } else if (name == "sail") {
                cfg.sail_angle_deg = val;
                g_sim.geom_changed = true;
            }
        }
    }
}

// --- Setup menu ---
static double read_double(const std::string& prompt, double def) {
    printf("%s [%.4g]: ", prompt.c_str(), def);
    std::string line; std::getline(std::cin, line);
    return line.empty() ? def : std::atof(line.c_str());
}
static int read_int(const std::string& prompt, int def) {
    printf("%s [%d]: ", prompt.c_str(), def);
    std::string line; std::getline(std::cin, line);
    return line.empty() ? def : std::atoi(line.c_str());
}

static Config setup_menu() {
    Config cfg;
    printf("\n========================================\n");
    printf("   3D Aerodynamics Simulator\n");
    printf("   Web-based 3D viewer + LBM engine\n");
    printf("========================================\n\n");

    printf("Select vehicle:\n  [1] Car  [2] Plane  [3] Sailboat\n");
    int c = read_int("Choice", 1);
    cfg.vehicle = (c==2) ? VehicleType::PLANE : (c==3) ? VehicleType::SAILBOAT : VehicleType::CAR;

    cfg.Re = read_double("Reynolds number", 150.0);

    printf("\nGrid: [1] Coarse (fast)  [2] Medium  [3] Fine\n");
    int r = read_int("Resolution", 1);
    cfg.set_grid_for_vehicle(r==2 ? "medium" : r==3 ? "fine" : "coarse");
    cfg.compute_derived();

    printf("\n  Grid: %dx%dx%d  tau=%.4f  Ma=%.3f\n\n",
           cfg.NX, cfg.NY, cfg.NZ, cfg.tau, cfg.U_LB/0.5774);
    return cfg;
}

static void do_revoxelize(Config& cfg, Lattice& lat, std::unique_ptr<Geometry>& geom,
                          ForceComputer& fc) {
    std::fill(lat.cell.begin(), lat.cell.end(), CellType::FLUID);
    cfg.compute_derived();
    geom = create_geometry(cfg);
    geom->voxelize(lat, cfg);
    cfg.A_ref = geom->reference_area(cfg);
    lat.initialize(1.0, cfg.U_LB, 0.0, 0.0);
    fc.find_surface_cells(lat);
}

int main() {
    Config cfg = setup_menu();
    mkdir("./output", 0755);

    printf("Allocating %dx%dx%d lattice (%.0f MB)...\n",
           cfg.NX, cfg.NY, cfg.NZ, cfg.memory_estimate_mb());
    Lattice lat;
    lat.allocate(cfg.NX, cfg.NY, cfg.NZ);

    auto geom = create_geometry(cfg);
    geom->voxelize(lat, cfg);
    cfg.A_ref = geom->reference_area(cfg);
    lat.initialize(1.0, cfg.U_LB, 0.0, 0.0);

    Solver solver;
    BoundaryConditions bc;
    ForceComputer fc;
    fc.find_surface_cells(lat);

    // Residual tracking
    std::vector<double> ux_old(lat.ux.begin(), lat.ux.end());
    std::vector<double> uy_old(lat.uy.begin(), lat.uy.end());
    std::vector<double> uz_old(lat.uz.begin(), lat.uz.end());

    g_sim.new_re = cfg.Re;

    // --- Start web server ---
    WebServer web;
    int port = 8080;
    while (!web.start(port) && port < 8090) port++;
    if (port >= 8090) { printf("Could not bind to any port!\n"); return 1; }

    web.route("/", [](const HttpRequest&) {
        return http_response("text/html", HTML_PAGE);
    });

    web.route("/state", [&](const HttpRequest&) {
        return http_response("application/json", json_state(cfg, lat));
    });

    web.route("/mesh", [&](const HttpRequest&) {
        return http_response("application/json", json_mesh(cfg, lat));
    });

    web.route("/control", [&](const HttpRequest& req) {
        handle_control(req.body, cfg);
        return http_response("application/json", "{\"ok\":true}");
    });

    printf("\n");
    printf("=========================================\n");
    printf("  3D viewer running at:\n");
    printf("  http://localhost:%d\n", port);
    printf("=========================================\n");
    printf("  Open in your browser to see the 3D view.\n");
    printf("  Press Ctrl+C to stop.\n\n");

    // --- Main simulation loop ---
    int force_interval = 20;

    while (true) {
        // Handle web requests (non-blocking)
        web.poll();

        // Handle commands from web UI
        if (g_sim.cmd_pause.exchange(false)) {
            g_sim.paused = !g_sim.paused;
        }
        if (g_sim.cmd_reset.exchange(false)) {
            lat.initialize(1.0, cfg.U_LB, 0.0, 0.0);
            ux_old.assign(lat.ux.begin(), lat.ux.end());
            uy_old.assign(lat.uy.begin(), lat.uy.end());
            uz_old.assign(lat.uz.begin(), lat.uz.end());
            g_sim.step = 0;
            g_sim.Cd = g_sim.Cl = g_sim.Cs = 0;
            g_sim.unstable = false;
        }
        if (g_sim.cmd_vehicle.exchange(false)) {
            int v = ((int)cfg.vehicle + 1) % 3;
            cfg.vehicle = (VehicleType)v;
            switch (cfg.vehicle) {
                case VehicleType::CAR: cfg.angle_of_attack_deg=0; cfg.slant_angle_deg=25; break;
                case VehicleType::PLANE: cfg.angle_of_attack_deg=5; cfg.yaw_deg=0; break;
                case VehicleType::SAILBOAT: cfg.sail_angle_deg=15; cfg.yaw_deg=0; break;
            }
            do_revoxelize(cfg, lat, geom, fc);
            ux_old.assign(lat.ux.begin(), lat.ux.end());
            uy_old.assign(lat.uy.begin(), lat.uy.end());
            uz_old.assign(lat.uz.begin(), lat.uz.end());
            g_sim.step = 0;
            g_sim.unstable = false;
        }
        if (g_sim.geom_changed.exchange(false)) {
            do_revoxelize(cfg, lat, geom, fc);
            ux_old.assign(lat.ux.begin(), lat.ux.end());
            uy_old.assign(lat.uy.begin(), lat.uy.end());
            uz_old.assign(lat.uz.begin(), lat.uz.end());
            g_sim.step = 0;
            g_sim.unstable = false;
        }

        // Run simulation step
        if (!g_sim.paused && !g_sim.unstable) {
            solver.collide(lat, cfg.tau);
            solver.stream(lat);
            bc.apply_all(lat, cfg);
            solver.compute_macroscopic(lat);
            g_sim.step++;

            if (g_sim.step % force_interval == 0) {
                ForceResult fr = fc.compute(lat, cfg.tau, cfg.U_LB, cfg.A_ref);
                g_sim.Cd = fr.Cd;
                g_sim.Cl = fr.Cl;
                g_sim.Cs = fr.Cs;

                double sum_diff = 0, sum_norm = 0, max_v2 = 0;
                size_t N = (size_t)cfg.NX * cfg.NY * cfg.NZ;
                for (size_t n = 0; n < N; ++n) {
                    if (lat.cell[n] == CellType::SOLID) continue;
                    double dux = lat.ux[n]-ux_old[n], duy = lat.uy[n]-uy_old[n], duz = lat.uz[n]-uz_old[n];
                    sum_diff += dux*dux + duy*duy + duz*duz;
                    sum_norm += lat.ux[n]*lat.ux[n] + lat.uy[n]*lat.uy[n] + lat.uz[n]*lat.uz[n];
                    double v2 = lat.ux[n]*lat.ux[n] + lat.uy[n]*lat.uy[n] + lat.uz[n]*lat.uz[n];
                    if (v2 > max_v2) max_v2 = v2;
                }
                g_sim.max_vel = std::sqrt(max_v2);
                g_sim.residual = std::sqrt(sum_diff / (sum_norm + 1e-30));

                ux_old.assign(lat.ux.begin(), lat.ux.end());
                uy_old.assign(lat.uy.begin(), lat.uy.end());
                uz_old.assign(lat.uz.begin(), lat.uz.end());

                if (g_sim.max_vel > 0.3 || std::isnan(g_sim.max_vel)) {
                    g_sim.unstable = true;
                }
            }
        } else {
            usleep(10000); // 10ms idle when paused
        }
    }

    web.stop();
    return 0;
}
