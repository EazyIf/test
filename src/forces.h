#ifndef FORCES_H
#define FORCES_H

#include "lattice.h"
#include <vector>

struct ForceResult {
    double Fx = 0, Fy = 0, Fz = 0;
    double Cd = 0, Cl = 0, Cs = 0;
};

class ForceComputer {
public:
    void find_surface_cells(const Lattice& lat);
    ForceResult compute(const Lattice& lat, double tau, double U_in, double A_ref);

private:
    struct SurfaceCell {
        int x, y, z;
        double nx, ny, nz;  // surface normal (pointing into body)
        double area;         // effective surface area
    };
    std::vector<SurfaceCell> surface_;
};

#endif
