#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "lattice.h"
#include "config.h"

class BoundaryConditions {
public:
    void apply_inlet_zou_he(Lattice& lat, double U_in);
    void apply_outlet_extrapolation(Lattice& lat);
    void apply_free_slip_y(Lattice& lat);
    void apply_free_slip_z(Lattice& lat);
    void apply_all(Lattice& lat, const Config& cfg);
};

#endif
