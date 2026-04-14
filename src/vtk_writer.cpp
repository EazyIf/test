#include "vtk_writer.h"
#include <fstream>
#include <cmath>
#include <cstdio>

void VTKWriter::write_3d_field(const Lattice& lat, const std::string& filename) {
    std::ofstream out(filename);
    if (!out) { printf("  Warning: could not open %s\n", filename.c_str()); return; }

    int NX = lat.NX, NY = lat.NY, NZ = lat.NZ;
    size_t N = (size_t)NX * NY * NZ;

    out << "# vtk DataFile Version 3.0\n";
    out << "LBM 3D Aerodynamics Simulation\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS " << NX << " " << NY << " " << NZ << "\n";
    out << "ORIGIN 0 0 0\n";
    out << "SPACING 1 1 1\n";
    out << "POINT_DATA " << N << "\n";

    // Velocity magnitude
    out << "SCALARS velocity_magnitude double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int z = 0; z < NZ; ++z)
        for (int y = 0; y < NY; ++y)
            for (int x = 0; x < NX; ++x) {
                size_t n = lat.idx(x, y, z);
                double vm = std::sqrt(lat.ux[n]*lat.ux[n] + lat.uy[n]*lat.uy[n] + lat.uz[n]*lat.uz[n]);
                out << vm << "\n";
            }

    // Pressure
    out << "SCALARS pressure double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int z = 0; z < NZ; ++z)
        for (int y = 0; y < NY; ++y)
            for (int x = 0; x < NX; ++x) {
                size_t n = lat.idx(x, y, z);
                out << lat.rho[n] / 3.0 << "\n";
            }

    // Velocity vectors
    out << "VECTORS velocity double\n";
    for (int z = 0; z < NZ; ++z)
        for (int y = 0; y < NY; ++y)
            for (int x = 0; x < NX; ++x) {
                size_t n = lat.idx(x, y, z);
                out << lat.ux[n] << " " << lat.uy[n] << " " << lat.uz[n] << "\n";
            }

    // Cell type
    out << "SCALARS cell_type int 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int z = 0; z < NZ; ++z)
        for (int y = 0; y < NY; ++y)
            for (int x = 0; x < NX; ++x) {
                out << static_cast<int>(lat.cell[lat.idx(x, y, z)]) << "\n";
            }

    printf("  Wrote 3D field: %s\n", filename.c_str());
}

void VTKWriter::write_slice_xz(const Lattice& lat, int y_slice, const std::string& filename) {
    std::ofstream out(filename);
    if (!out) return;

    int NX = lat.NX, NZ = lat.NZ;

    out << "# vtk DataFile Version 3.0\n";
    out << "LBM XZ slice at y=" << y_slice << "\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS " << NX << " 1 " << NZ << "\n";
    out << "ORIGIN 0 " << y_slice << " 0\n";
    out << "SPACING 1 1 1\n";
    out << "POINT_DATA " << NX * NZ << "\n";

    out << "SCALARS velocity_magnitude double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int z = 0; z < NZ; ++z)
        for (int x = 0; x < NX; ++x) {
            size_t n = lat.idx(x, y_slice, z);
            double vm = std::sqrt(lat.ux[n]*lat.ux[n] + lat.uy[n]*lat.uy[n] + lat.uz[n]*lat.uz[n]);
            out << vm << "\n";
        }

    out << "SCALARS pressure double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int z = 0; z < NZ; ++z)
        for (int x = 0; x < NX; ++x) {
            size_t n = lat.idx(x, y_slice, z);
            out << lat.rho[n] / 3.0 << "\n";
        }

    out << "VECTORS velocity double\n";
    for (int z = 0; z < NZ; ++z)
        for (int x = 0; x < NX; ++x) {
            size_t n = lat.idx(x, y_slice, z);
            out << lat.ux[n] << " " << lat.uy[n] << " " << lat.uz[n] << "\n";
        }

    out << "SCALARS cell_type int 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int z = 0; z < NZ; ++z)
        for (int x = 0; x < NX; ++x)
            out << static_cast<int>(lat.cell[lat.idx(x, y_slice, z)]) << "\n";
}

void VTKWriter::write_slice_xy(const Lattice& lat, int z_slice, const std::string& filename) {
    std::ofstream out(filename);
    if (!out) return;

    int NX = lat.NX, NY = lat.NY;

    out << "# vtk DataFile Version 3.0\n";
    out << "LBM XY slice at z=" << z_slice << "\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS " << NX << " " << NY << " 1\n";
    out << "ORIGIN 0 0 " << z_slice << "\n";
    out << "SPACING 1 1 1\n";
    out << "POINT_DATA " << NX * NY << "\n";

    out << "SCALARS velocity_magnitude double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int y = 0; y < NY; ++y)
        for (int x = 0; x < NX; ++x) {
            size_t n = lat.idx(x, y, z_slice);
            double vm = std::sqrt(lat.ux[n]*lat.ux[n] + lat.uy[n]*lat.uy[n] + lat.uz[n]*lat.uz[n]);
            out << vm << "\n";
        }

    out << "SCALARS pressure double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int y = 0; y < NY; ++y)
        for (int x = 0; x < NX; ++x) {
            size_t n = lat.idx(x, y, z_slice);
            out << lat.rho[n] / 3.0 << "\n";
        }

    out << "VECTORS velocity double\n";
    for (int y = 0; y < NY; ++y)
        for (int x = 0; x < NX; ++x) {
            size_t n = lat.idx(x, y, z_slice);
            out << lat.ux[n] << " " << lat.uy[n] << " " << lat.uz[n] << "\n";
        }

    out << "SCALARS cell_type int 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int y = 0; y < NY; ++y)
        for (int x = 0; x < NX; ++x)
            out << static_cast<int>(lat.cell[lat.idx(x, y, z_slice)]) << "\n";
}

void VTKWriter::write_surface(const Lattice& lat, const std::string& filename) {
    // Collect surface cells (fluid cells adjacent to at least one solid cell)
    struct SurfPoint { int x, y, z; double pressure; };
    std::vector<SurfPoint> points;

    for (int x = 0; x < lat.NX; ++x)
        for (int y = 0; y < lat.NY; ++y)
            for (int z = 0; z < lat.NZ; ++z) {
                size_t n = lat.idx(x, y, z);
                if (lat.cell[n] != CellType::FLUID) continue;

                bool near_solid = false;
                for (int i = 1; i < Q; ++i) {
                    int xn = x + ex[i], yn = y + ey[i], zn = z + ez[i];
                    if (xn >= 0 && xn < lat.NX && yn >= 0 && yn < lat.NY &&
                        zn >= 0 && zn < lat.NZ &&
                        lat.cell[lat.idx(xn, yn, zn)] == CellType::SOLID) {
                        near_solid = true;
                        break;
                    }
                }

                if (near_solid) {
                    points.push_back({x, y, z, lat.rho[n] / 3.0});
                }
            }

    std::ofstream out(filename);
    if (!out) return;

    out << "# vtk DataFile Version 3.0\n";
    out << "Surface pressure data\n";
    out << "ASCII\n";
    out << "DATASET POLYDATA\n";
    out << "POINTS " << points.size() << " double\n";
    for (const auto& p : points)
        out << p.x << " " << p.y << " " << p.z << "\n";

    out << "POINT_DATA " << points.size() << "\n";
    out << "SCALARS surface_pressure double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (const auto& p : points)
        out << p.pressure << "\n";

    printf("  Wrote surface data: %s (%zu points)\n", filename.c_str(), points.size());
}
