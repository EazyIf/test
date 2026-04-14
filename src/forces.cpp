#include "forces.h"
#include <cstdio>
#include <cmath>

void ForceComputer::find_surface_cells(const Lattice& lat) {
    surface_.clear();

    for (int x = 0; x < lat.NX; ++x) {
        for (int y = 0; y < lat.NY; ++y) {
            for (int z = 0; z < lat.NZ; ++z) {
                size_t n = lat.idx(x, y, z);
                if (lat.cell[n] != CellType::FLUID) continue;

                // Compute surface normal from axis-aligned solid neighbors only
                double rnx = 0, rny = 0, rnz = 0;
                for (int i = 1; i <= 6; ++i) {  // axis-aligned directions only
                    int xn = x + ex[i], yn = y + ey[i], zn = z + ez[i];
                    if (xn < 0 || xn >= lat.NX || yn < 0 || yn >= lat.NY || zn < 0 || zn >= lat.NZ)
                        continue;
                    if (lat.cell[lat.idx(xn, yn, zn)] == CellType::SOLID) {
                        rnx += ex[i];
                        rny += ey[i];
                        rnz += ez[i];
                    }
                }

                double mag = std::sqrt(rnx * rnx + rny * rny + rnz * rnz);
                if (mag < 1e-10) continue;

                // Unit area per surface cell, normal from axis-aligned neighbors
                surface_.push_back({x, y, z, rnx / mag, rny / mag, rnz / mag, 1.0});
            }
        }
    }

    printf("  Found %zu surface cells for force computation\n", surface_.size());
}

ForceResult ForceComputer::compute(const Lattice& lat, double tau, double U_in, double A_ref) {
    ForceResult fr;
    const double cs2 = 1.0 / 3.0;
    const double visc_factor = 1.0 - 0.5 / tau;
    const size_t N = (size_t)lat.NX * lat.NY * lat.NZ;

    for (const auto& sc : surface_) {
        size_t n = lat.idx(sc.x, sc.y, sc.z);
        double rho = lat.rho[n];
        double ux = lat.ux[n];
        double uy = lat.uy[n];
        double uz = lat.uz[n];

        // Pressure force: -(p - p_ref) * n, where p = rho*cs^2, p_ref = 1*cs^2
        double dp = (rho - 1.0) * cs2;

        // Non-equilibrium stress tensor: Pi^neq_ab = sum_i (f_i - f_i^eq) * e_ia * e_ib
        double Pxx = 0, Pxy = 0, Pxz = 0;
        double Pyy = 0, Pyz = 0, Pzz = 0;

        for (int i = 0; i < Q; ++i) {
            double fi = lat.f[(size_t)i * N + n];
            double eu = ex[i] * ux + ey[i] * uy + ez[i] * uz;
            double usq = ux * ux + uy * uy + uz * uz;
            double feq = w[i] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * usq);
            double fneq = fi - feq;

            Pxx += fneq * ex[i] * ex[i];
            Pxy += fneq * ex[i] * ey[i];
            Pxz += fneq * ex[i] * ez[i];
            Pyy += fneq * ey[i] * ey[i];
            Pyz += fneq * ey[i] * ez[i];
            Pzz += fneq * ez[i] * ez[i];
        }

        // Total stress on body: F = (dp * n + visc_factor * Pi^neq * n) * area
        double fx = (dp * sc.nx + visc_factor * (Pxx * sc.nx + Pxy * sc.ny + Pxz * sc.nz)) * sc.area;
        double fy = (dp * sc.ny + visc_factor * (Pxy * sc.nx + Pyy * sc.ny + Pyz * sc.nz)) * sc.area;
        double fz = (dp * sc.nz + visc_factor * (Pxz * sc.nx + Pyz * sc.ny + Pzz * sc.nz)) * sc.area;

        fr.Fx += fx;
        fr.Fy += fy;
        fr.Fz += fz;
    }

    double denom = 0.5 * 1.0 * U_in * U_in * A_ref;
    if (denom > 1e-15) {
        fr.Cd = fr.Fx / denom;
        fr.Cl = fr.Fz / denom;
        fr.Cs = fr.Fy / denom;
    }

    return fr;
}
