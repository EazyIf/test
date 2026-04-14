#include "geometry.h"
#include "geometry_car.h"
#include "geometry_plane.h"
#include "geometry_boat.h"

void Geometry::voxelize(Lattice& lat, const Config& cfg) {
    int solid_count = 0;

    for (int x = 0; x < lat.NX; ++x) {
        for (int y = 0; y < lat.NY; ++y) {
            for (int z = 0; z < lat.NZ; ++z) {
                double px = x - cx;
                double py = y - cy;
                double pz = z - cz;

                // Apply yaw rotation (about z-axis) — inverse rotation
                double yaw = cfg.yaw_deg * M_PI / 180.0;
                double rx = px * cos(yaw) + py * sin(yaw);
                double ry = -px * sin(yaw) + py * cos(yaw);

                // Apply angle of attack (about y-axis) — inverse rotation
                double aoa = cfg.angle_of_attack_deg * M_PI / 180.0;
                double rx2 = rx * cos(aoa) + pz * sin(aoa);
                double rz = -rx * sin(aoa) + pz * cos(aoa);

                double d = sdf(rx2, ry, rz);
                size_t n = lat.idx(x, y, z);

                if (d < 0.0) {
                    lat.cell[n] = CellType::SOLID;
                    solid_count++;
                }
            }
        }
    }

    // Mark inlet and outlet
    for (int y = 0; y < lat.NY; ++y) {
        for (int z = 0; z < lat.NZ; ++z) {
            size_t n0 = lat.idx(0, y, z);
            if (lat.cell[n0] != CellType::SOLID)
                lat.cell[n0] = CellType::INLET;

            size_t n1 = lat.idx(lat.NX - 1, y, z);
            if (lat.cell[n1] != CellType::SOLID)
                lat.cell[n1] = CellType::OUTLET;
        }
    }

    printf("  Voxelization: %d solid cells out of %d total (%.1f%%)\n",
           solid_count, lat.NX * lat.NY * lat.NZ,
           100.0 * solid_count / ((double)lat.NX * lat.NY * lat.NZ));
}

std::unique_ptr<Geometry> create_geometry(const Config& cfg) {
    switch (cfg.vehicle) {
        case VehicleType::CAR:
            return std::make_unique<AhmedBody>(cfg);
        case VehicleType::PLANE:
            return std::make_unique<Aircraft>(cfg);
        case VehicleType::SAILBOAT:
            return std::make_unique<Sailboat>(cfg);
    }
    return nullptr;
}
