#ifndef VTK_WRITER_H
#define VTK_WRITER_H

#include "lattice.h"
#include <string>

class VTKWriter {
public:
    void write_3d_field(const Lattice& lat, const std::string& filename);
    void write_slice_xz(const Lattice& lat, int y_slice, const std::string& filename);
    void write_slice_xy(const Lattice& lat, int z_slice, const std::string& filename);
    void write_surface(const Lattice& lat, const std::string& filename);
};

#endif
