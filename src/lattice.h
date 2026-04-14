#ifndef LATTICE_H
#define LATTICE_H

#include <vector>
#include <cstdint>
#include <cstddef>

constexpr int Q = 19;

constexpr int ex[Q] = {0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0};
constexpr int ey[Q] = {0, 0, 0, 1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 1,-1, 1,-1};
constexpr int ez[Q] = {0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1};

constexpr double w[Q] = {
    1.0/3.0,
    1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0,
    1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,
    1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,
    1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0
};

constexpr int opp[Q] = {0, 2,1, 4,3, 6,5, 10,9,8,7, 14,13,12,11, 18,17,16,15};

constexpr int mirror_y[Q] = {0,1,2, 4,3, 5,6, 9,10,7,8, 11,12,13,14, 16,15,18,17};
constexpr int mirror_z[Q] = {0,1,2, 3,4, 6,5, 7,8,9,10, 13,14,11,12, 17,18,15,16};

enum class CellType : uint8_t { FLUID = 0, SOLID = 1, INLET = 2, OUTLET = 3 };

struct Lattice {
    int NX, NY, NZ;
    std::vector<double> f;
    std::vector<double> f_new;
    std::vector<CellType> cell;
    std::vector<double> rho;
    std::vector<double> ux, uy, uz;

    inline size_t idx(int x, int y, int z) const {
        return (size_t)x * NY * NZ + (size_t)y * NZ + z;
    }

    inline size_t fidx(int i, int x, int y, int z) const {
        return (size_t)i * NX * NY * NZ + idx(x, y, z);
    }

    void allocate(int nx, int ny, int nz);
    void initialize(double rho0, double ux0, double uy0, double uz0);
};

#endif
