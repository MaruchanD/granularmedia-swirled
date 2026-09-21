#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>
#include <algorithm>
#include <string>
#include <utility>
#include <cstdlib>

void save_graph_svg(const std::vector<double>& t,
                    const std::vector<double>& values,
                    double tau,
                    const std::string& filename,
                    const std::string& title,
                    const std::string& ylabel) {
    if (t.empty() || values.empty() || t.size() != values.size()) return;

    const int W = 1100, H = 650;
    const int left = 90, right = 40, top = 60, bottom = 80;

    double xmin = t.front();
    double xmax = t.back();
    double ymin = *std::min_element(values.begin(), values.end());
    double ymax = *std::max_element(values.begin(), values.end());

    if (xmax <= xmin) xmax = xmin + 1.0;
    if (ymax <= ymin) ymax = ymin + 1.0;

    double margen = 0.05 * (ymax - ymin);
    ymin -= margen;
    ymax += margen;

    auto xmap = [&](double x) {
        return left + (x - xmin) / (xmax - xmin) * (W - left - right);
    };

    auto ymap = [&](double y) {
        return H - bottom - (y - ymin) / (ymax - ymin) *
               (H - top - bottom);
    };

    std::ofstream svg(filename);
    if (!svg) return;

    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
        << "width=\"" << W << "\" height=\"" << H << "\">\n";
    svg << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

    // ejes principales
    svg << "<line x1=\"" << left << "\" y1=\"" << top
        << "\" x2=\"" << left << "\" y2=\"" << H - bottom
        << "\" stroke=\"black\" stroke-width=\"1.5\"/>\n";

    svg << "<line x1=\"" << left << "\" y1=\"" << H - bottom
        << "\" x2=\"" << W - right << "\" y2=\"" << H - bottom
        << "\" stroke=\"black\" stroke-width=\"1.5\"/>\n";

    // escala horizontal: ticks + etiquetas
    int nx = 6;
    for (int i = 0; i <= nx; ++i) {
        double xv = xmin + (xmax - xmin) * i / nx;
        double px = xmap(xv);

        svg << "<line x1=\"" << px << "\" y1=\"" << H - bottom
            << "\" x2=\"" << px << "\" y2=\"" << H - bottom + 6
            << "\" stroke=\"black\" stroke-width=\"1\"/>\n";

        svg << "<text x=\"" << px << "\" y=\"" << H - bottom + 18
            << "\" font-size=\"10\" text-anchor=\"middle\">"
            << std::fixed << std::setprecision(2) << xv
            << "</text>\n";
    }

    // escala vertical: ticks + etiquetas
    int ny = 6;
    for (int i = 0; i <= ny; ++i) {
        double yv = ymin + (ymax - ymin) * (ny - i) / ny;
        double py = ymap(yv);

        svg << "<line x1=\"" << left - 6 << "\" y1=\"" << py
            << "\" x2=\"" << left << "\" y2=\"" << py
            << "\" stroke=\"black\" stroke-width=\"1\"/>\n";

        svg << "<text x=\"" << left - 10 << "\" y=\"" << py + 3
            << "\" font-size=\"10\" text-anchor=\"end\">"
            << std::fixed << std::setprecision(2) << yv
            << "</text>\n";
    }

    // curva
    svg << "<polyline fill=\"none\" stroke=\"blue\" stroke-width=\"2\" points=\"";
    for (size_t i = 0; i < t.size(); ++i)
        svg << xmap(t[i]) << "," << ymap(values[i]) << " ";
    svg << "\"/>\n";

    // tiempo de relajación
    if (tau >= xmin && tau <= xmax) {
        double xtau = xmap(tau);

        svg << "<line x1=\"" << xtau << "\" y1=\"" << top
            << "\" x2=\"" << xtau << "\" y2=\"" << H - bottom
            << "\" stroke=\"red\" stroke-width=\"2\" "
            << "stroke-dasharray=\"7,5\"/>\n";

        svg << "<text x=\"" << xtau + 8 << "\" y=\"" << top + 22
            << "\" fill=\"red\" font-size=\"12\">tau = " << tau << " s</text>\n";
    }

    // título y etiquetas
    svg << "<text x=\"" << W / 2 << "\" y=\"30\" "
        << "text-anchor=\"middle\" font-size=\"20\">"
        << title << "</text>\n";

    svg << "<text x=\"" << W / 2 << "\" y=\"" << H - 25
        << "\" text-anchor=\"middle\" font-size=\"12\">Tiempo t [s]</text>\n";

    svg << "<text x=\"22\" y=\"" << H / 2
        << "\" transform=\"rotate(-90 22," << H / 2
        << ")\" text-anchor=\"middle\" font-size=\"12\">" << ylabel << "</text>\n";

    svg << "</svg>\n";
}

struct Params {
    int    n          = 40;     // número de partículas
    double k_n        = 5000.0;  // rigidez de repulsión elástica interpartícula [N/m]
    double gamma_n    = 2.5;     // amortiguamiento normal viscoso interpartícula [kg/s]
    double r_part     = 0.15;    // radio físico de la partícula [m]
    double k          = 1000.0;  // rigidez de la pared blanda del plato [N/m]
    double R_plato    = 6.0;     // radio del plato [m]
    double R          = 1.5;     // distancia del origen al centro del plato [m]
    double Omega      = 1.2;     // velocidad angular orbital del plato [rad/s]
    double t_max      = 100.0;   // tiempo total de simulación [s]
    double dt         = 0.0005;  // paso de integración [s]
    int    save_every = 40;     // guardar 1 de cada 'save_every' pasos
};

// Potencial de pared blanda del plato
inline double Vborde(double rho, const Params& p) {
    double d = rho - p.R_plato;
    return (d > 0.0) ? 0.5 * p.k * d * d : 0.0;
}

inline double dVborde_drho(double rho, const Params& p) {
    double d = rho - p.R_plato;
    return (d > 0.0) ? p.k * d : 0.0;
}

// Términos de arrastre cinemático del plato orbital
inline void driftAB(double phi, double t, const Params& p, double& A, double& B) {
    double psi = phi - p.Omega * t;
    double ROmega = p.R * p.Omega;
    A = ROmega * std::sin(psi);
    B = ROmega * std::cos(psi);
}

// Estado local: velocidades reales y momentos efectivos
void local_state(double t, const std::vector<double>& y, const std::vector<double>& m,
                 const Params& p, std::vector<double>& rho_dot, std::vector<double>& phi_dot,
                 std::vector<double>& vx, std::vector<double>& vy,
                 std::vector<double>& Prho_eff, std::vector<double>& Pphi_eff) {
    int n = p.n;
    for (int i = 0; i < n; ++i) {
        double rho = y[4*i+0], phi = y[4*i+1], prho = y[4*i+2], pphi = y[4*i+3], mi = m[i];
        double A, B;
        driftAB(phi, t, p, A, B);
        Prho_eff[i] = prho - mi * A;
        Pphi_eff[i] = pphi - mi * rho * B;
        rho_dot[i] = Prho_eff[i] / mi;
        phi_dot[i] = Pphi_eff[i] / (mi * rho * rho);
        vx[i] = rho_dot[i] * std::cos(phi) - rho * phi_dot[i] * std::sin(phi);
        vy[i] = rho_dot[i] * std::sin(phi) + rho * phi_dot[i] * std::cos(phi);
    }
}

struct Vec2D {
    double x = 0.0;
    double y = 0.0;
};

inline Vec2D operator+(const Vec2D& a, const Vec2D& b) {
    return {a.x + b.x, a.y + b.y};
}

inline Vec2D operator-(const Vec2D& a, const Vec2D& b) {
    return {a.x - b.x, a.y - b.y};
}

inline Vec2D operator*(double s, const Vec2D& v) {
    return {s * v.x, s * v.y};
}

inline Vec2D operator*(const Vec2D& v, double s) {
    return {v.x * s, v.y * s};
}

inline Vec2D operator/(const Vec2D& v, double s) {
    return {v.x / s, v.y / s};
}

inline double dot(const Vec2D& a, const Vec2D& b) {
    return a.x * b.x + a.y * b.y;
}

inline double norm_sq(const Vec2D& v) {
    return dot(v, v);
}

inline double norm(const Vec2D& v) {
    return std::hypot(v.x, v.y);
}

inline Vec2D calcular_fuerzas_contacto(int i,
                                const std::vector<double>& y,
                                const std::vector<double>& vx,
                                const std::vector<double>& vy,
                                const Params& p) {
    Vec2D F = {0.0, 0.0};
    double rho_i = y[4 * i + 0];
    double phi_i = y[4 * i + 1];

    double qxi = rho_i * std::cos(phi_i);
    double qyi = rho_i * std::sin(phi_i);
    double suma_radios = 2.0 * p.r_part;

    for (int j = 0; j < p.n; ++j) {
        if (i == j) continue;

        double rho_j = y[4 * j + 0];
        double phi_j = y[4 * j + 1];
        double qxj = rho_j * std::cos(phi_j);
        double qyj = rho_j * std::sin(phi_j);

        double dx = qxi - qxj;
        double dy = qyi - qyj;
        double dist = std::hypot(dx, dy);

        if (dist < suma_radios && dist > 0.0) {
            double delta = suma_radios - dist;
            double nx = dx / dist;
            double ny = dy / dist;

            double dvx = vx[i] - vx[j];
            double dvy = vy[i] - vy[j];
            double v_rel_n = dvx * nx + dvy * ny;

            double F_mag = p.k_n * delta - p.gamma_n * v_rel_n;

            F.x += F_mag * nx;
            F.y += F_mag * ny;
        }
    }
    return F;
}

std::vector<double> derivatives(double t, const std::vector<double>& y,
                                const std::vector<double>& m, const Params& p) {
    int n = p.n;
    std::vector<double> dydt(4 * n);
    std::vector<double> rho_dot(n), phi_dot(n), vx(n), vy(n), Prho_eff(n), Pphi_eff(n);

    local_state(t, y, m, p, rho_dot, phi_dot, vx, vy, Prho_eff, Pphi_eff);

    for (int i = 0; i < n; ++i) {
        double rho = y[4*i+0], phi = y[4*i+1], mi = m[i];
        double A, B;
        driftAB(phi, t, p, A, B);
        double cphi = std::cos(phi), sphi = std::sin(phi);

        Vec2D F_dem = calcular_fuerzas_contacto(i, y, vx, vy, p);

        double Fx_pol = F_dem.x * cphi + F_dem.y * sphi;
        double Fy_pol = -F_dem.x * sphi + F_dem.y * cphi;

        double prho_dot = B * Pphi_eff[i] / (rho * rho)
                        + (Pphi_eff[i] * Pphi_eff[i]) / (mi * rho * rho * rho)
                        - dVborde_drho(rho, p)
                        + Fx_pol;

        double pphi_dot = B * Prho_eff[i] - A * Pphi_eff[i] / rho
                        + rho * Fy_pol;

        dydt[4*i+0] = rho_dot[i];
        dydt[4*i+1] = phi_dot[i];
        dydt[4*i+2] = prho_dot;
        dydt[4*i+3] = pphi_dot;
    }
    return dydt;
}

void euler_richardson_step(double t, std::vector<double>& y, const std::vector<double>& m,
                          const Params& p) {
    double dt = p.dt;
    auto f1 = derivatives(t, y, m, p);
    std::vector<double> y_mid(y.size());
    for (size_t j = 0; j < y.size(); ++j) y_mid[j] = y[j] + 0.5 * dt * f1[j];
    auto f_mid = derivatives(t + 0.5 * dt, y_mid, m, p);
    for (size_t j = 0; j < y.size(); ++j) y[j] += dt * f_mid[j];
}

// Energía Potencial de deformación por contacto elástico interpartícula
double V_contacto(const std::vector<double>& y, const Params& p) {
    double V_tot = 0.0;
    int n = p.n;
    double suma_radios = 2.0 * p.r_part;

    for (int i = 0; i < n; ++i) {
        double xi = y[4 * i + 0] * std::cos(y[4 * i + 1]);
        double yi = y[4 * i + 0] * std::sin(y[4 * i + 1]);

        for (int j = i + 1; j < n; ++j) {
            double xj = y[4 * j + 0] * std::cos(y[4 * j + 1]);
            double yj = y[4 * j + 0] * std::sin(y[4 * j + 1]);

            double dx = xi - xj;
            double dy = yi - yj;
            double dist = std::hypot(dx, dy);

            if (dist < suma_radios && dist > 0.0) {
                double delta = suma_radios - dist;
                V_tot += 0.5 * p.k_n * delta * delta;
            }
        }
    }
    return V_tot;
}

// Hamiltoniano Canónico Total
double hamiltonian_canonical(double t, const std::vector<double>& y,
                            const std::vector<double>& m, const Params& p) {
    int n = p.n;
    double H = 0.0;
    for (int i = 0; i < n; ++i) {
        double rho = y[4*i+0], phi = y[4*i+1], prho = y[4*i+2], pphi = y[4*i+3], mi = m[i];
        double A, B;
        driftAB(phi, t, p, A, B);
        double Prho_eff = prho - mi * A;
        double Pphi_eff = pphi - mi * rho * B;
        H += Prho_eff * Prho_eff / (2.0 * mi)
           + Pphi_eff * Pphi_eff / (2.0 * mi * rho * rho)
           - 0.5 * mi * p.R * p.R * p.Omega * p.Omega
           + Vborde(rho, p);
    }
    return H + V_contacto(y, p);
}

// Energía Física Real Total
double physical_energy(double t, const std::vector<double>& y,
                      const std::vector<double>& m, const Params& p) {
    int n = p.n;
    double ROmega = p.R * p.Omega;
    double E = 0.0;
    std::vector<double> rho_dot(n), phi_dot(n), vx(n), vy(n), Prho_eff(n), Pphi_eff(n);
    local_state(t, y, m, p, rho_dot, phi_dot, vx, vy, Prho_eff, Pphi_eff);
    for (int i = 0; i < n; ++i) {
        double rho = y[4*i+0], phi = y[4*i+1], mi = m[i];
        double psi = phi - p.Omega * t;

        double T = 0.5 * mi * rho_dot[i] * rho_dot[i]
                 + 0.5 * mi * rho * rho * phi_dot[i] * phi_dot[i]
                 + mi * ROmega * (rho_dot[i] * std::sin(psi) + rho * phi_dot[i] * std::cos(psi))
                 + 0.5 * mi * ROmega * ROmega;

        E += T + Vborde(rho, p);
    }
    return E + V_contacto(y, p);
}

std::vector<double> moving_average(const std::vector<double>& x, int window_n) {
    int N = (int)x.size();
    std::vector<double> out(N);
    int half = window_n / 2;
    for (int i = 0; i < N; ++i) {
        int lo = std::max(0, i - half);
        int hi = std::min(N - 1, i + half);
        double s = 0.0;
        for (int j = lo; j <= hi; ++j) s += x[j];
        out[i] = s / (hi - lo + 1);
    }
    return out;
}

double period_envelope_relaxation_time(const std::vector<double>& tvec,
                                       const std::vector<double>& x, double T_orb,
                                       double tol_rel = 0.02) {
    double t0 = tvec.front(), t1 = tvec.back();
    int n_cycles = (int)std::floor((t1 - t0) / T_orb);
    if (n_cycles < 3) return -1.0;

    std::vector<double> cycle_start(n_cycles), peak(n_cycles), trough(n_cycles);
    size_t idx = 0;
    for (int c = 0; c < n_cycles; ++c) {
        double lo = t0 + c * T_orb, hi = t0 + (c + 1) * T_orb;
        double mx = -1e300, mn = 1e300;
        while (idx < tvec.size() && tvec[idx] < hi) {
            if (tvec[idx] >= lo) {
                mx = std::max(mx, x[idx]);
                mn = std::min(mn, x[idx]);
            }
            ++idx;
        }
        while (idx > 0 && tvec[idx - 1] >= hi) --idx;
        cycle_start[c] = lo;
        peak[c] = mx;
        trough[c] = mn;
    }

    auto settling_cycle = [&](const std::vector<double>& seq) -> int {
        int N = (int)seq.size();
        int tail = std::max(2, N / 4);
        double ref = 0.0;
        for (int i = N - tail; i < N; ++i) ref += seq[i];
        ref /= tail;

        double amp = *std::max_element(seq.begin(), seq.end())
                   - *std::min_element(seq.begin(), seq.end());
        double band = tol_rel * std::max(amp, 1e-9);

        int last_outside = -1;
        for (int i = 0; i < N; ++i) {
            if (std::fabs(seq[i] - ref) > band) last_outside = i;
        }
        return (last_outside == -1) ? 0 : std::min(last_outside + 1, N - 1);
    };

    int c_peak = settling_cycle(peak);
    int c_trough = settling_cycle(trough);
    int c_final = std::max(c_peak, c_trough);
    return cycle_start[c_final];
}

double wall_potential_from_state(const std::vector<double>& y, const Params& p) {
    double V = 0.0;
    for (int i = 0; i < p.n; ++i) {
        double rho = y[4 * i + 0];
        double d = rho - p.R_plato;
        if (d > 0.0) {
            V += 0.5 * p.k * d * d;
        }
    }
    return V;
}

double kinetic_energy(double t,
                      const std::vector<double>& y,
                      const std::vector<double>& m,
                      const Params& p) {
    int n = p.n;
    std::vector<double> rho_dot(n), phi_dot(n), vx(n), vy(n), Prho_eff(n), Pphi_eff(n);
    local_state(t, y, m, p, rho_dot, phi_dot, vx, vy, Prho_eff, Pphi_eff);

    double K = 0.0;
    for (int i = 0; i < n; ++i) {
        double mi = m[i];
        K += 0.5 * mi * (vx[i] * vx[i] + vy[i] * vy[i]);
    }
    return K;
}

double cluster_angular_velocity(double t,
                               const std::vector<double>& y,
                               const std::vector<double>& m,
                               const Params& p) {
    int n = p.n;
    std::vector<double> rho_dot(n), phi_dot(n), vx(n), vy(n), Prho_eff(n), Pphi_eff(n);
    local_state(t, y, m, p, rho_dot, phi_dot, vx, vy, Prho_eff, Pphi_eff);

    double M = 0.0;
    Vec2D cm = {0.0, 0.0};

    for (int i = 0; i < n; ++i) {
        double rho = y[4*i+0], phi = y[4*i+1], mi = m[i];
        double x = rho * std::cos(phi);
        double yv = rho * std::sin(phi);
        cm.x += mi * x;
        cm.y += mi * yv;
        M += mi;
    }

    if (M < 1e-12) return 0.0;
    cm.x /= M; cm.y /= M;

    double Lz = 0.0;
    double I = 0.0;

    for (int i = 0; i < n; ++i) {
        double rho = y[4*i+0], phi = y[4*i+1], mi = m[i];
        double x = rho * std::cos(phi);
        double yv = rho * std::sin(phi);

        Vec2D r_rel = {x - cm.x, yv - cm.y};
        Vec2D v_rel = {vx[i], vy[i]};

        double r2 = norm_sq(r_rel);
        I += mi * r2;

        Lz += mi * (r_rel.x * v_rel.y - r_rel.y * v_rel.x);
    }

    if (std::abs(I) < 1e-12) return 0.0;
    return Lz / I;
}

void save_energy_svg(const std::vector<double>& t,
                     const std::vector<double>& Kvec,
                     const std::vector<double>& Wvec,
                     const std::vector<double>& Etotal,
                     double tau,
                     const std::string& filename) {
    const int W = 1100, H = 650;
    const int left = 90, right = 40, top = 60, bottom = 80;

    if (t.empty() || Kvec.empty() || Wvec.empty() || Etotal.empty()) return;
    if (t.size() != Kvec.size() || t.size() != Wvec.size() || t.size() != Etotal.size()) return;

    double xmin = t.front();
    double xmax = t.back();

    double ymin = std::min({*std::min_element(Kvec.begin(), Kvec.end()),
                            *std::min_element(Wvec.begin(), Wvec.end()),
                            *std::min_element(Etotal.begin(), Etotal.end())});
    double ymax = std::max({*std::max_element(Kvec.begin(), Kvec.end()),
                            *std::max_element(Wvec.begin(), Wvec.end()),
                            *std::max_element(Etotal.begin(), Etotal.end())});

    if (ymax <= ymin) ymax = ymin + 1.0;
    double margen = 0.05 * (ymax - ymin);
    ymin -= margen;
    ymax += margen;

    auto xmap = [&](double x) {
        return left + (x - xmin) / (xmax - xmin) * (W - left - right);
    };

    auto ymap = [&](double y) {
        return H - bottom - (y - ymin) / (ymax - ymin) * (H - top - bottom);
    };

    std::ofstream svg(filename);
    if (!svg) return;

    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << W
        << "\" height=\"" << H << "\">\n";
    svg << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

    svg << "<line x1=\"" << left << "\" y1=\"" << top
        << "\" x2=\"" << left << "\" y2=\"" << H - bottom
        << "\" stroke=\"black\" stroke-width=\"1.5\"/>\n";
    svg << "<line x1=\"" << left << "\" y1=\"" << H - bottom
        << "\" x2=\"" << W - right << "\" y2=\"" << H - bottom
        << "\" stroke=\"black\" stroke-width=\"1.5\"/>\n";

    int nx = 6;
    for (int i = 0; i <= nx; ++i) {
        double xv = xmin + (xmax - xmin) * i / nx;
        double px = xmap(xv);

        svg << "<line x1=\"" << px << "\" y1=\"" << H - bottom
            << "\" x2=\"" << px << "\" y2=\"" << H - bottom + 6
            << "\" stroke=\"black\" stroke-width=\"1\"/>\n";
        svg << "<text x=\"" << px << "\" y=\"" << H - bottom + 18
            << "\" font-size=\"10\" text-anchor=\"middle\">"
            << std::fixed << std::setprecision(2) << xv
            << "</text>\n";
    }

    int ny = 6;
    for (int i = 0; i <= ny; ++i) {
        double yv = ymin + (ymax - ymin) * (ny - i) / ny;
        double py = ymap(yv);

        svg << "<line x1=\"" << left - 6 << "\" y1=\"" << py
            << "\" x2=\"" << left << "\" y2=\"" << py
            << "\" stroke=\"black\" stroke-width=\"1\"/>\n";
        svg << "<text x=\"" << left - 10 << "\" y=\"" << py + 3
            << "\" font-size=\"10\" text-anchor=\"end\">"
            << std::fixed << std::setprecision(2) << yv
            << "</text>\n";
    }

    svg << "<polyline fill=\"none\" stroke=\"blue\" stroke-width=\"2\" points=\"";
    for (size_t i = 0; i < t.size(); ++i) {
        svg << xmap(t[i]) << "," << ymap(Kvec[i]) << " ";
    }
    svg << "\"/>\n";

    svg << "<polyline fill=\"none\" stroke=\"green\" stroke-width=\"2\" points=\"";
    for (size_t i = 0; i < t.size(); ++i) {
        svg << xmap(t[i]) << "," << ymap(Wvec[i]) << " ";
    }
    svg << "\"/>\n";

    svg << "<polyline fill=\"none\" stroke=\"red\" stroke-width=\"2\" points=\"";
    for (size_t i = 0; i < t.size(); ++i) {
        svg << xmap(t[i]) << "," << ymap(Etotal[i]) << " ";
    }
    svg << "\"/>\n";

    if (tau >= xmin && tau <= xmax) {
        double xtau = xmap(tau);
        svg << "<line x1=\"" << xtau << "\" y1=\"" << top
            << "\" x2=\"" << xtau << "\" y2=\"" << H - bottom
            << "\" stroke=\"red\" stroke-width=\"2\" stroke-dasharray=\"7,5\"/>\n";
        svg << "<text x=\"" << xtau + 8 << "\" y=\"" << top + 22
            << "\" fill=\"red\" font-size=\"12\">tau = "
            << std::fixed << std::setprecision(2) << tau << " s</text>\n";
    }

    svg << "<line x1=\"120\" y1=\"55\" x2=\"150\" y2=\"55\" stroke=\"blue\" stroke-width=\"3\"/>\n";
    svg << "<text x=\"160\" y=\"60\" font-size=\"14\" fill=\"blue\">Energía cinética</text>\n";

    svg << "<line x1=\"330\" y1=\"55\" x2=\"360\" y2=\"55\" stroke=\"green\" stroke-width=\"3\"/>\n";
    svg << "<text x=\"370\" y=\"60\" font-size=\"14\" fill=\"green\">Energía del contenedor</text>\n";

    svg << "<line x1=\"660\" y1=\"55\" x2=\"690\" y2=\"55\" stroke=\"red\" stroke-width=\"3\"/>\n";
    svg << "<text x=\"700\" y=\"60\" font-size=\"14\" fill=\"red\">Energía total</text>\n";

    svg << "<text x=\"" << W / 2 << "\" y=\"25\" text-anchor=\"middle\" font-size=\"18\">"
        << "Energía del sistema</text>\n";
    svg << "<text x=\"" << W / 2 << "\" y=\"" << H - 20
        << "\" text-anchor=\"middle\" font-size=\"12\">Tiempo t [s]</text>\n";
    svg << "<text x=\"20\" y=\"" << H / 2
        << "\" transform=\"rotate(-90 20," << H / 2
        << ")\" text-anchor=\"middle\" font-size=\"12\">Energía [J]</text>\n";
    svg << "</svg>\n";
}

void save_frequency_svg(const std::vector<double>& t,
                        const std::vector<double>& omega_plato,
                        const std::vector<double>& omega_cluster,
                        const std::string& filename) {
    const int W = 1100, H = 650;
    const int left = 90, right = 40, top = 60, bottom = 80;

    if (t.empty() || omega_plato.empty() || omega_cluster.empty()) return;
    if (t.size() != omega_plato.size() || t.size() != omega_cluster.size()) return;

    double xmin = t.front();
    double xmax = t.back();
    double ymin = std::min(*std::min_element(omega_plato.begin(), omega_plato.end()),
                          *std::min_element(omega_cluster.begin(), omega_cluster.end()));
    double ymax = std::max(*std::max_element(omega_plato.begin(), omega_plato.end()),
                          *std::max_element(omega_cluster.begin(), omega_cluster.end()));

    if (ymax <= ymin) ymax = ymin + 1.0;
    double margen = 0.05 * (ymax - ymin);
    ymin -= margen;
    ymax += margen;

    auto xmap = [&](double x) {
        return left + (x - xmin) / (xmax - xmin) * (W - left - right);
    };

    auto ymap = [&](double y) {
        return H - bottom - (y - ymin) / (ymax - ymin) * (H - top - bottom);
    };

    std::ofstream svg(filename);
    if (!svg) return;

    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << W << "\" height=\"" << H << "\">\n";
    svg << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

    svg << "<line x1=\"" << left << "\" y1=\"" << top
        << "\" x2=\"" << left << "\" y2=\"" << H - bottom
        << "\" stroke=\"black\" stroke-width=\"1.5\"/>\n";
    svg << "<line x1=\"" << left << "\" y1=\"" << H - bottom
        << "\" x2=\"" << W - right << "\" y2=\"" << H - bottom
        << "\" stroke=\"black\" stroke-width=\"1.5\"/>\n";

    int nx = 6;
    for (int i = 0; i <= nx; ++i) {
        double xv = xmin + (xmax - xmin) * i / nx;
        double px = xmap(xv);

        svg << "<line x1=\"" << px << "\" y1=\"" << H - bottom
            << "\" x2=\"" << px << "\" y2=\"" << H - bottom + 6
            << "\" stroke=\"black\" stroke-width=\"1\"/>\n";
        svg << "<text x=\"" << px << "\" y=\"" << H - bottom + 18
            << "\" font-size=\"10\" text-anchor=\"middle\">"
            << std::fixed << std::setprecision(2) << xv
            << "</text>\n";
    }

    int ny = 6;
    for (int i = 0; i <= ny; ++i) {
        double yv = ymin + (ymax - ymin) * (ny - i) / ny;
        double py = ymap(yv);

        svg << "<line x1=\"" << left - 6 << "\" y1=\"" << py
            << "\" x2=\"" << left << "\" y2=\"" << py
            << "\" stroke=\"black\" stroke-width=\"1\"/>\n";
        svg << "<text x=\"" << left - 10 << "\" y=\"" << py + 3
            << "\" font-size=\"10\" text-anchor=\"end\">"
            << std::fixed << std::setprecision(2) << yv
            << "</text>\n";
    }

    svg << "<polyline fill=\"none\" stroke=\"blue\" stroke-width=\"2\" points=\"";
    for (size_t i = 0; i < t.size(); ++i) {
        svg << xmap(t[i]) << "," << ymap(omega_plato[i]) << " ";
    }
    svg << "\"/>\n";

    svg << "<polyline fill=\"none\" stroke=\"red\" stroke-width=\"2\" points=\"";
    for (size_t i = 0; i < t.size(); ++i) {
        svg << xmap(t[i]) << "," << ymap(omega_cluster[i]) << " ";
    }
    svg << "\"/>\n";

    svg << "<line x1=\"120\" y1=\"55\" x2=\"150\" y2=\"55\" stroke=\"blue\" stroke-width=\"3\"/>\n";
    svg << "<text x=\"160\" y=\"60\" font-size=\"14\" fill=\"blue\">Omega plato</text>\n";

    svg << "<line x1=\"360\" y1=\"55\" x2=\"390\" y2=\"55\" stroke=\"red\" stroke-width=\"3\"/>\n";
    svg << "<text x=\"400\" y=\"60\" font-size=\"14\" fill=\"red\">Omega cluster</text>\n";

    svg << "<text x=\"" << W / 2 << "\" y=\"25\" text-anchor=\"middle\" font-size=\"18\">"
        << "Frecuencia angular</text>\n";
    svg << "<text x=\"" << W / 2 << "\" y=\"" << H - 20
        << "\" text-anchor=\"middle\" font-size=\"12\">Tiempo t [s]</text>\n";
    svg << "<text x=\"20\" y=\"" << H / 2
        << "\" transform=\"rotate(-90 20," << H / 2
        << ")\" text-anchor=\"middle\" font-size=\"12\">Omega [rad/s]</text>\n";
    svg << "</svg>\n";
}

int main() {
    Params p;

    std::mt19937 rng(42);
    std::uniform_real_distribution<double> U_rho(0.25 * p.R_plato, 0.85 * p.R_plato);
    std::uniform_real_distribution<double> U_phi(0.0, 2.0 * std::acos(-1.0));
    std::uniform_real_distribution<double> U_vrad(-0.8, 0.8);
    std::uniform_real_distribution<double> U_omega(1.0, 3.0);
    std::uniform_real_distribution<double> U_mass(0.8, 1.2);

    std::vector<double> m(p.n), y(4 * p.n);

    for (int i = 0; i < p.n; ++i) {
        double mi = U_mass(rng);
        double rho0 = U_rho(rng);
        double phi0 = U_phi(rng);
        double vrad0 = U_vrad(rng);
        double om0 = U_omega(rng);

        double A0, B0;
        driftAB(phi0, 0.0, p, A0, B0);

        m[i] = mi;
        y[4*i + 0] = rho0;
        y[4*i + 1] = phi0;
        y[4*i + 2] = mi * vrad0 + mi * A0;
        y[4*i + 3] = mi * rho0 * rho0 * om0 + mi * rho0 * B0;
    }

    int nsteps = (int)std::round(p.t_max / p.dt);

    std::vector<double> tvec;
    std::vector<double> Ekin_vec, Ewall_vec, Etotal_vec;
    std::vector<double> Oplato_vec, Ocluster_vec;

    tvec.reserve(nsteps / p.save_every + 1);
    Ekin_vec.reserve(nsteps / p.save_every + 1);
    Ewall_vec.reserve(nsteps / p.save_every + 1);
    Etotal_vec.reserve(nsteps / p.save_every + 1);
    Oplato_vec.reserve(nsteps / p.save_every + 1);
    Ocluster_vec.reserve(nsteps / p.save_every + 1);

    double t = 0.0;
    for (int step = 0; step <= nsteps; ++step) {
        if (step % p.save_every == 0) {
            tvec.push_back(t);

            double K = kinetic_energy(t, y, m, p);
            double Vwall = wall_potential_from_state(y, p);
            double Etot = K + Vwall + V_contacto(y, p);

            Ekin_vec.push_back(K);
            Ewall_vec.push_back(Vwall);
            Etotal_vec.push_back(Etot);

            Oplato_vec.push_back(p.Omega);
            Ocluster_vec.push_back(cluster_angular_velocity(t, y, m, p));
        }
        if (step == nsteps) break;
        euler_richardson_step(t, y, m, p);
        t += p.dt;
    }

    double T_orb = 2.0 * std::acos(-1.0) / std::fabs(p.Omega);
    double sample_dt = p.dt * p.save_every;

    // Cálculo del tiempo de relajación usando la energía del sistema
    double tau = period_envelope_relaxation_time(tvec, Etotal_vec, T_orb, 0.02);

    std::ofstream fout("dem_spring_dashpot.csv");
    fout << "t,E_kinetica,E_contenedor,E_total,Omega_plato,Omega_cluster,tau\n";
    for (size_t i = 0; i < tvec.size(); ++i) {
        fout << std::setprecision(10)
             << tvec[i] << ","
             << Ekin_vec[i] << ","
             << Ewall_vec[i] << ","
             << Etotal_vec[i] << ","
             << Oplato_vec[i] << ","
             << Ocluster_vec[i] << ","
             << tau << "\n";
    }
    fout.close();

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Simulación completada.\n";
    std::cout << "Tiempo de relajación tau = " << tau << " s\n";
    std::cout << "Omega_plato final = " << Oplato_vec.back() << " rad/s\n";
    std::cout << "Omega_cluster final = " << Ocluster_vec.back() << " rad/s\n";
    std::cout << "Error de acople = " << std::abs(Ocluster_vec.back() - Oplato_vec.back()) << " rad/s\n";

    save_energy_svg(tvec, Ekin_vec, Ewall_vec, Etotal_vec, tau, "energia_sistema.svg");
    save_frequency_svg(tvec, Oplato_vec, Ocluster_vec, "frecuencia_cluster_vs_plato.svg");

#ifdef _WIN32
    std::system("start energia_sistema.svg");
    std::system("start frecuencia_cluster_vs_plato.svg");
#else
    std::system("xdg-open energia_sistema.svg");
    std::system("xdg-open frecuencia_cluster_vs_plato.svg");
#endif

    return 0;
}