#include "lattice.h"

void Lattice::allocate(int nx, int ny, int nz) {
    NX = nx; NY = ny; NZ = nz;
    size_t N = (size_t)NX * NY * NZ;

    f.resize((size_t)Q * N, 0.0);
    f_new.resize((size_t)Q * N, 0.0);
    cell.resize(N, CellType::FLUID);
    rho.resize(N, 1.0);
    ux.resize(N, 0.0);
    uy.resize(N, 0.0);
    uz.resize(N, 0.0);
}

void Lattice::initialize(double rho0, double ux0, double uy0, double uz0) {
    for (int x = 0; x < NX; ++x) {
        for (int y = 0; y < NY; ++y) {
            for (int z = 0; z < NZ; ++z) {
                size_t n = idx(x, y, z);
                rho[n] = rho0;

                double vx = (cell[n] == CellType::SOLID) ? 0.0 : ux0;
                double vy = (cell[n] == CellType::SOLID) ? 0.0 : uy0;
                double vz = (cell[n] == CellType::SOLID) ? 0.0 : uz0;

                ux[n] = vx;
                uy[n] = vy;
                uz[n] = vz;

                for (int i = 0; i < Q; ++i) {
                    double eu = ex[i] * vx + ey[i] * vy + ez[i] * vz;
                    double usq = vx * vx + vy * vy + vz * vz;
                    f[fidx(i, x, y, z)] = w[i] * rho0 * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * usq);
                }
            }
        }
    }
}
