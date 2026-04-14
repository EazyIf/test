#ifndef TUI_H
#define TUI_H

#include "config.h"
#include "lattice.h"
#include <termios.h>
#include <unistd.h>
#include <sys/ioctl.h>
#include <csignal>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

// --- Terminal raw mode (RAII) ---

static termios g_original_termios;
static bool g_raw_mode_active = false;

static void restore_terminal() {
    if (g_raw_mode_active) {
        tcsetattr(STDIN_FILENO, TCSANOW, &g_original_termios);
        printf("\033[?25h\033[0m\n");
        g_raw_mode_active = false;
    }
}

static void sigint_handler(int) {
    restore_terminal();
    _exit(0);
}

struct RawMode {
    RawMode() {
        tcgetattr(STDIN_FILENO, &g_original_termios);
        termios raw = g_original_termios;
        raw.c_lflag &= ~(ICANON | ECHO);
        raw.c_cc[VMIN] = 0;
        raw.c_cc[VTIME] = 0;
        tcsetattr(STDIN_FILENO, TCSANOW, &raw);
        g_raw_mode_active = true;
        signal(SIGINT, sigint_handler);
        atexit(restore_terminal);
    }
    ~RawMode() { restore_terminal(); }
};

// --- Key input ---

enum class KeyAction {
    NONE, UP, DOWN, LEFT, RIGHT,
    PAUSE, RESET, VEHICLE, QUIT,
    SPEED_UP, SPEED_DOWN
};

static KeyAction poll_key() {
    char c;
    if (read(STDIN_FILENO, &c, 1) != 1) return KeyAction::NONE;

    if (c == '\033') {
        char seq[2];
        if (read(STDIN_FILENO, &seq[0], 1) != 1) return KeyAction::NONE;
        if (read(STDIN_FILENO, &seq[1], 1) != 1) return KeyAction::NONE;
        if (seq[0] == '[') {
            switch (seq[1]) {
                case 'A': return KeyAction::UP;
                case 'B': return KeyAction::DOWN;
                case 'C': return KeyAction::RIGHT;
                case 'D': return KeyAction::LEFT;
            }
        }
        return KeyAction::NONE;
    }

    switch (c) {
        case ' ':           return KeyAction::PAUSE;
        case 'r': case 'R': return KeyAction::RESET;
        case 'v': case 'V': return KeyAction::VEHICLE;
        case 'q': case 'Q': return KeyAction::QUIT;
        case '+': case '=': return KeyAction::SPEED_UP;
        case '-':           return KeyAction::SPEED_DOWN;
    }
    return KeyAction::NONE;
}

// --- Adjustable parameter ---

struct AdjustableParam {
    const char* label;
    double* value;
    double step_size;
    double min_val, max_val;
    bool requires_revoxelize;
};

static std::vector<AdjustableParam> build_param_list(Config& cfg) {
    std::vector<AdjustableParam> params;
    params.push_back({"Reynolds", &cfg.Re, 25.0, 25.0, 2000.0, false});

    switch (cfg.vehicle) {
        case VehicleType::CAR:
            params.push_back({"Yaw (deg)", &cfg.yaw_deg, 1.0, -30.0, 30.0, true});
            params.push_back({"Slant (deg)", &cfg.slant_angle_deg, 1.0, 5.0, 40.0, true});
            break;
        case VehicleType::PLANE:
            params.push_back({"AoA (deg)", &cfg.angle_of_attack_deg, 1.0, -10.0, 20.0, true});
            params.push_back({"Yaw (deg)", &cfg.yaw_deg, 1.0, -20.0, 20.0, true});
            break;
        case VehicleType::SAILBOAT:
            params.push_back({"Yaw (deg)", &cfg.yaw_deg, 1.0, -30.0, 30.0, true});
            params.push_back({"Sail (deg)", &cfg.sail_angle_deg, 1.0, -45.0, 45.0, true});
            break;
    }
    return params;
}

// --- TUI state ---

struct TUIState {
    int selected_param = 0;
    bool paused = false;
    int render_skip = 1;
    int step = 0;

    double Cd = 0, Cl = 0, Cs = 0;
    double residual = 0, max_vel = 0;
    bool unstable = false;

    std::string status_msg;
};

// --- ASCII flow visualization ---

static const char palette[] = " .:-=+*#%@";
static const int palette_len = 10;

static void render_flow_slice(const Lattice& lat, int y_slice,
                              int viz_cols, int viz_rows, char* buf, double v_scale) {
    double x_scale = (double)lat.NX / viz_cols;
    double z_scale = (double)lat.NZ / viz_rows;

    for (int r = 0; r < viz_rows; ++r) {
        int z = lat.NZ - 1 - (int)(r * z_scale);
        if (z < 0) z = 0;
        if (z >= lat.NZ) z = lat.NZ - 1;

        for (int c = 0; c < viz_cols; ++c) {
            int x = (int)(c * x_scale);
            if (x >= lat.NX) x = lat.NX - 1;

            size_t n = lat.idx(x, y_slice, z);

            if (lat.cell[n] == CellType::SOLID) {
                buf[r * viz_cols + c] = '#';
            } else {
                double vmag = std::sqrt(lat.ux[n]*lat.ux[n] + lat.uy[n]*lat.uy[n] + lat.uz[n]*lat.uz[n]);
                int idx = (int)(vmag / v_scale * (palette_len - 1));
                if (idx >= palette_len) idx = palette_len - 1;
                if (idx < 0) idx = 0;
                buf[r * viz_cols + c] = palette[idx];
            }
        }
    }
}

// --- Screen drawing ---

static void draw_screen(const TUIState& st, const Config& cfg,
                        const std::vector<AdjustableParam>& params,
                        const char* flow_buf, int viz_cols, int viz_rows) {
    // Get terminal size
    struct winsize ws;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &ws);
    int term_w = ws.ws_col > 0 ? ws.ws_col : 80;

    std::string out;
    out.reserve(16384);

    out += "\033[?25l\033[H";  // hide cursor, home

    char line[256];

    // --- Title ---
    const char* status = st.paused ? "\033[1;33m PAUSED \033[0m" :
                         st.unstable ? "\033[1;31m UNSTABLE \033[0m" :
                         "\033[1;32m RUNNING \033[0m";
    snprintf(line, sizeof(line), "\033[1;37m  3D Aerodynamics Simulator \033[0m [%s]", status);
    out += line;
    out += "\033[K\n";

    // --- Step / Speed ---
    snprintf(line, sizeof(line), "  Step: %-8d  Speed: x%-3d",
             st.step, st.render_skip);
    out += line;
    out += "\033[K\n\033[K\n";

    // --- Vehicle ---
    const char* vname = "";
    switch (cfg.vehicle) {
        case VehicleType::CAR:      vname = "Car (Ahmed Body)"; break;
        case VehicleType::PLANE:    vname = "Aircraft"; break;
        case VehicleType::SAILBOAT: vname = "Sailboat"; break;
    }
    snprintf(line, sizeof(line), "  \033[1mVehicle:\033[0m %s    Grid: %dx%dx%d   L=%d",
             vname, cfg.NX, cfg.NY, cfg.NZ, cfg.L_ref);
    out += line;
    out += "\033[K\n\033[K\n";

    // --- Forces ---
    out += "  \033[1;4mAerodynamic Forces\033[0m\033[K\n";
    snprintf(line, sizeof(line), "  Cd = \033[1;36m%10.4f\033[0m   Cl = \033[1;36m%10.4f\033[0m   Cs = \033[1;36m%10.4f\033[0m",
             st.Cd, st.Cl, st.Cs);
    out += line;
    out += "\033[K\n";
    snprintf(line, sizeof(line), "  Residual: %.2e    Max vel: %.4f (Ma=%.3f)    tau: %.4f",
             st.residual, st.max_vel, st.max_vel / 0.5774, cfg.tau);
    out += line;
    out += "\033[K\n\033[K\n";

    // --- Parameters ---
    out += "  \033[1;4mParameters\033[0m  (\033[33m\xe2\x86\x91\xe2\x86\x93\033[0m select  \033[33m\xe2\x86\x90\xe2\x86\x92\033[0m adjust)\033[K\n";
    for (int i = 0; i < (int)params.size(); ++i) {
        if (i == st.selected_param) {
            snprintf(line, sizeof(line), "  \033[1;33m> %-14s  [ %8.1f ]\033[0m %s",
                     params[i].label, *params[i].value,
                     params[i].requires_revoxelize ? "\033[2m(restart)\033[0m" : "");
        } else {
            snprintf(line, sizeof(line), "    %-14s  [ %8.1f ] %s",
                     params[i].label, *params[i].value,
                     params[i].requires_revoxelize ? "\033[2m(restart)\033[0m" : "");
        }
        out += line;
        out += "\033[K\n";
    }
    out += "\033[K\n";

    // --- Flow visualization ---
    snprintf(line, sizeof(line), "  \033[1;4mFlow Field\033[0m  (XZ slice at y=%d)  \033[2m[ .:-=+*#%%@ ]\033[0m",
             cfg.NY / 2);
    out += line;
    out += "\033[K\n";

    int max_viz_w = std::min(viz_cols, term_w - 4);
    for (int r = 0; r < viz_rows; ++r) {
        out += "  ";
        out.append(flow_buf + r * viz_cols, max_viz_w);
        out += "\033[K\n";
    }
    out += "\033[K\n";

    // --- Help bar ---
    out += "  \033[7m SPACE \033[0m pause  \033[7m r \033[0m reset  \033[7m v \033[0m vehicle  \033[7m q \033[0m quit  \033[7m +/- \033[0m speed\033[K\n";

    // Status message
    if (!st.status_msg.empty()) {
        out += "  \033[1;33m" + st.status_msg + "\033[0m\033[K\n";
    } else {
        out += "\033[K\n";
    }

    // Clear remaining lines
    out += "\033[J";

    (void)!write(STDOUT_FILENO, out.data(), out.size());
}

#endif
