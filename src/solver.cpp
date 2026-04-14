#include "solver.h"
#include <cmath>
#include <algorithm>

void Solver::collide(Lattice& lat, double tau) {
    // TRT (Two-Relaxation-Time) collision operator
    // Magic parameter Lambda = (tau+ - 0.5)(tau- - 0.5) = 1/4
    const double tau_plus = tau;
    const double tau_minus = 0.5 + 1.0 / (4.0 * (tau_plus - 0.5));

    const double inv_tau_p = 1.0 / tau_plus;
    const double inv_tau_m = 1.0 / tau_minus;

    const size_t N = (size_t)lat.NX * lat.NY * lat.NZ;

    for (size_t n = 0; n < N; ++n) {
        if (lat.cell[n] == CellType::SOLID) continue;

        double r = lat.rho[n];
        double vx = lat.ux[n];
        double vy = lat.uy[n];
        double vz = lat.uz[n];

        // Compute equilibrium for all directions
        double feq[Q];
        for (int i = 0; i < Q; ++i)
            feq[i] = equilibrium(i, r, vx, vy, vz);

        // Process rest direction (i=0, self-paired) with BGK
        {
            size_t fi = n;  // i=0, fidx = 0*N + n
            lat.f[fi] -= inv_tau_p * (lat.f[fi] - feq[0]);
        }

        // Process paired directions: (1,2), (3,4), (5,6), (7,10), (8,9),
        // (11,14), (12,13), (15,18), (16,17)
        // For each pair (i, j=opp[i]) with i < j, compute both simultaneously
        static const int pairs[][2] = {
            {1,2}, {3,4}, {5,6}, {7,10}, {8,9},
            {11,14}, {12,13}, {15,18}, {16,17}
        };

        for (const auto& p : pairs) {
            int i = p[0], j = p[1];
            size_t fi_idx = (size_t)i * N + n;
            size_t fj_idx = (size_t)j * N + n;

            double fi = lat.f[fi_idx];
            double fj = lat.f[fj_idx];

            // Symmetric (even) and antisymmetric (odd) non-equilibrium parts
            double neq_plus  = 0.5 * (fi + fj) - 0.5 * (feq[i] + feq[j]);
            double neq_minus = 0.5 * (fi - fj) - 0.5 * (feq[i] - feq[j]);

            // Update both directions simultaneously
            lat.f[fi_idx] = fi - inv_tau_p * neq_plus - inv_tau_m * neq_minus;
            lat.f[fj_idx] = fj - inv_tau_p * neq_plus + inv_tau_m * neq_minus;
        }
    }
}

void Solver::stream(Lattice& lat) {
    const int NX = lat.NX, NY = lat.NY, NZ = lat.NZ;

    for (int x = 0; x < NX; ++x) {
        for (int y = 0; y < NY; ++y) {
            for (int z = 0; z < NZ; ++z) {
                size_t n = lat.idx(x, y, z);
                if (lat.cell[n] == CellType::SOLID) {
                    for (int i = 0; i < Q; ++i)
                        lat.f_new[lat.fidx(i, x, y, z)] = lat.f[lat.fidx(i, x, y, z)];
                    continue;
                }

                for (int i = 0; i < Q; ++i) {
                    int xs = x - ex[i];
                    int ys = y - ey[i];
                    int zs = z - ez[i];

                    if (xs < 0 || xs >= NX || ys < 0 || ys >= NY || zs < 0 || zs >= NZ) {
                        lat.f_new[lat.fidx(i, x, y, z)] = lat.f[lat.fidx(i, x, y, z)];
                    } else if (lat.cell[lat.idx(xs, ys, zs)] == CellType::SOLID) {
                        lat.f_new[lat.fidx(i, x, y, z)] = lat.f[lat.fidx(opp[i], x, y, z)];
                    } else {
                        lat.f_new[lat.fidx(i, x, y, z)] = lat.f[lat.fidx(i, xs, ys, zs)];
                    }
                }
            }
        }
    }

    std::swap(lat.f, lat.f_new);
}

void Solver::compute_macroscopic(Lattice& lat) {
    const size_t N = (size_t)lat.NX * lat.NY * lat.NZ;

    for (size_t n = 0; n < N; ++n) {
        if (lat.cell[n] == CellType::SOLID) {
            lat.rho[n] = 1.0;
            lat.ux[n] = lat.uy[n] = lat.uz[n] = 0.0;
            continue;
        }

        double r = 0.0, vx = 0.0, vy = 0.0, vz = 0.0;
        for (int i = 0; i < Q; ++i) {
            double fi = lat.f[(size_t)i * N + n];
            r += fi;
            vx += fi * ex[i];
            vy += fi * ey[i];
            vz += fi * ez[i];
        }

        lat.rho[n] = r;
        double inv_r = 1.0 / r;
        lat.ux[n] = vx * inv_r;
        lat.uy[n] = vy * inv_r;
        lat.uz[n] = vz * inv_r;
    }
}
