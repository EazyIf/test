#ifndef SOLVER_H
#define SOLVER_H

#include "lattice.h"

class Solver {
public:
    void collide(Lattice& lat, double tau);
    void stream(Lattice& lat);
    void compute_macroscopic(Lattice& lat);

private:
    inline double equilibrium(int i, double rho, double ux, double uy, double uz) const {
        double eu = ex[i] * ux + ey[i] * uy + ez[i] * uz;
        double usq = ux * ux + uy * uy + uz * uz;
        return w[i] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * usq);
    }
};

#endif
