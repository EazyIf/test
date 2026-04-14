#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "lattice.h"
#include "config.h"
#include <memory>
#include <string>

class Geometry {
public:
    virtual ~Geometry() = default;

    // Signed distance function: negative = inside solid
    virtual double sdf(double x, double y, double z) const = 0;

    // Populate cell types in the lattice based on SDF
    void voxelize(Lattice& lat, const Config& cfg);

    // Compute reference frontal area (for force coefficients)
    virtual double reference_area(const Config& cfg) const = 0;

    // Name for display
    virtual std::string name() const = 0;

    // Center of body in lattice coordinates
    double cx, cy, cz;
};

std::unique_ptr<Geometry> create_geometry(const Config& cfg);

#endif
