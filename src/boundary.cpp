#include "boundary.h"

void BoundaryConditions::apply_inlet_zou_he(Lattice& lat, double U_in) {
    const int NY = lat.NY, NZ = lat.NZ;
    const int x = 0;

    for (int y = 0; y < NY; ++y) {
        for (int z = 0; z < NZ; ++z) {
            size_t n = lat.idx(x, y, z);
            if (lat.cell[n] != CellType::INLET) continue;

            // Known distributions at inlet (x=0)
            // Unknown: f1(+x), f7(+x+y), f9(+x-y), f11(+x+z), f13(+x-z)
            // Known inward from -x: f2, f8, f10, f12, f14
            double f0  = lat.f[lat.fidx(0, x, y, z)];
            double f2  = lat.f[lat.fidx(2, x, y, z)];
            double f3  = lat.f[lat.fidx(3, x, y, z)];
            double f4  = lat.f[lat.fidx(4, x, y, z)];
            double f5  = lat.f[lat.fidx(5, x, y, z)];
            double f6  = lat.f[lat.fidx(6, x, y, z)];
            double f8  = lat.f[lat.fidx(8, x, y, z)];
            double f10 = lat.f[lat.fidx(10, x, y, z)];
            double f12 = lat.f[lat.fidx(12, x, y, z)];
            double f14 = lat.f[lat.fidx(14, x, y, z)];
            double f15 = lat.f[lat.fidx(15, x, y, z)];
            double f16 = lat.f[lat.fidx(16, x, y, z)];
            double f17 = lat.f[lat.fidx(17, x, y, z)];
            double f18 = lat.f[lat.fidx(18, x, y, z)];

            double rho_minus = f2 + f8 + f10 + f12 + f14;
            double rho_zero = f0 + f3 + f4 + f5 + f6 + f15 + f16 + f17 + f18;

            double rho = (rho_zero + 2.0 * rho_minus) / (1.0 - U_in);

            double ru = rho * U_in;

            // D3Q19 Zou-He (Hecht & Harting 2010)
            // Correction terms for transverse momentum balance
            double Ny = 0.5 * (f3 + f15 + f17 - f4 - f16 - f18);
            double Nz = 0.5 * (f5 + f15 + f16 - f6 - f17 - f18);

            lat.f[lat.fidx(1, x, y, z)]  = f2  + (1.0/3.0) * ru;
            lat.f[lat.fidx(7, x, y, z)]  = f10 + (1.0/6.0) * ru - Ny;
            lat.f[lat.fidx(9, x, y, z)]  = f8  + (1.0/6.0) * ru + Ny;
            lat.f[lat.fidx(11, x, y, z)] = f14 + (1.0/6.0) * ru - Nz;
            lat.f[lat.fidx(13, x, y, z)] = f12 + (1.0/6.0) * ru + Nz;

            lat.rho[n] = rho;
            lat.ux[n] = U_in;
            lat.uy[n] = 0.0;
            lat.uz[n] = 0.0;
        }
    }
}

void BoundaryConditions::apply_outlet_extrapolation(Lattice& lat) {
    const int NX = lat.NX, NY = lat.NY, NZ = lat.NZ;
    const int x = NX - 1;
    const int xs = NX - 2;

    for (int y = 0; y < NY; ++y) {
        for (int z = 0; z < NZ; ++z) {
            size_t n = lat.idx(x, y, z);
            if (lat.cell[n] != CellType::OUTLET) continue;

            // Zero-gradient extrapolation for all distributions
            for (int i = 0; i < Q; ++i) {
                lat.f[lat.fidx(i, x, y, z)] = lat.f[lat.fidx(i, xs, y, z)];
            }
        }
    }
}

void BoundaryConditions::apply_free_slip_y(Lattice& lat) {
    const int NX = lat.NX, NY = lat.NY, NZ = lat.NZ;

    // y = 0 face — skip inlet/outlet/solid cells to avoid overwriting their BCs
    for (int x = 1; x < NX - 1; ++x) {
        for (int z = 0; z < NZ; ++z) {
            size_t n = lat.idx(x, 0, z);
            if (lat.cell[n] != CellType::FLUID) continue;
            for (int i = 0; i < Q; ++i) {
                if (ey[i] < 0) {
                    lat.f[lat.fidx(mirror_y[i], x, 0, z)] = lat.f[lat.fidx(i, x, 0, z)];
                }
            }
        }
    }

    // y = NY-1 face
    for (int x = 1; x < NX - 1; ++x) {
        for (int z = 0; z < NZ; ++z) {
            size_t n = lat.idx(x, NY - 1, z);
            if (lat.cell[n] != CellType::FLUID) continue;
            for (int i = 0; i < Q; ++i) {
                if (ey[i] > 0) {
                    lat.f[lat.fidx(mirror_y[i], x, NY - 1, z)] = lat.f[lat.fidx(i, x, NY - 1, z)];
                }
            }
        }
    }
}

void BoundaryConditions::apply_free_slip_z(Lattice& lat) {
    const int NX = lat.NX, NY = lat.NY, NZ = lat.NZ;

    // z = 0 face — skip inlet/outlet/solid cells
    for (int x = 1; x < NX - 1; ++x) {
        for (int y = 0; y < NY; ++y) {
            size_t n = lat.idx(x, y, 0);
            if (lat.cell[n] != CellType::FLUID) continue;
            for (int i = 0; i < Q; ++i) {
                if (ez[i] < 0) {
                    lat.f[lat.fidx(mirror_z[i], x, y, 0)] = lat.f[lat.fidx(i, x, y, 0)];
                }
            }
        }
    }

    // z = NZ-1 face
    for (int x = 1; x < NX - 1; ++x) {
        for (int y = 0; y < NY; ++y) {
            size_t n = lat.idx(x, y, NZ - 1);
            if (lat.cell[n] != CellType::FLUID) continue;
            for (int i = 0; i < Q; ++i) {
                if (ez[i] > 0) {
                    lat.f[lat.fidx(mirror_z[i], x, y, NZ - 1)] = lat.f[lat.fidx(i, x, y, NZ - 1)];
                }
            }
        }
    }
}

void BoundaryConditions::apply_all(Lattice& lat, const Config& cfg) {
    apply_inlet_zou_he(lat, cfg.U_LB);
    apply_outlet_extrapolation(lat);
    apply_free_slip_y(lat);
    apply_free_slip_z(lat);
}
