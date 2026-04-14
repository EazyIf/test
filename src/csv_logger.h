#ifndef CSV_LOGGER_H
#define CSV_LOGGER_H

#include "forces.h"
#include <fstream>
#include <string>

class CSVLogger {
public:
    CSVLogger() = default;

    void open(const std::string& filename) {
        file_.open(filename);
        file_ << "timestep,Fx,Fy,Fz,Cd,Cl,Cs\n";
    }

    void log(int timestep, const ForceResult& fr) {
        if (!file_.is_open()) return;
        file_ << timestep << ","
              << fr.Fx << "," << fr.Fy << "," << fr.Fz << ","
              << fr.Cd << "," << fr.Cl << "," << fr.Cs << "\n";
        file_.flush();
    }

    void close() { if (file_.is_open()) file_.close(); }

private:
    std::ofstream file_;
};

#endif
