#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>
#include <algorithm>
#include <string>
#include <cstdlib>

// Estructura vectorial 2D para sustituir StaticArrays de Julia
struct Vec2D {
    double x = 0.0, y = 0.0;

    Vec2D operator+(const Vec2D& o) const { return {x + o.x, y + o.y}; }
    Vec2D operator-(const Vec2D& o) const { return {x - o.x, y - o.y}; }
    Vec2D operator*(double s) const { return {x * s, y * s}; }
    Vec2D operator/(double s) const { return {x / s, y / s}; }
    Vec2D& operator+=(const Vec2D& o) { x += o.x; y += o.y; return *this; }
    Vec2D& operator-=(const Vec2D& o) { x -= o.x; y -= o.y; return *this; }
};

inline Vec2D operator*(double s, const Vec2D& v) { return {s * v.x, s * v.y}; }
inline double dot(const Vec2D& a, const Vec2D& b) { return a.x * b.x + a.y * b.y; }
inline double norm_sq(const Vec2D& v) { return v.x * v.x + v.y * v.y; }
inline double norm(const Vec2D& v) { return std::sqrt(norm_sq(v)); }

struct Particle {
    Vec2D r;       // Posición
    Vec2D v;       // Velocidad
    Vec2D a;       // Aceleración
    double mass;   // Masa
    double radius; // Radio
};

struct Params {
    int n = 40;                  // Número de partículas (igualado al código en C++)
    double R_plato = 6.0;         // Radio del contenedor [m]
    double R = 1.5;               // Amplitud de excitación / Radio orbital [m]
    double Omega = 1.2;           // Velocidad angular orbital [rad/s]
    double kn_wall = 1000.0;      // Rigidez de la pared [N/m] (igualado a p.k)
    double gamma_wall = 0.5;      // Coeficiente de amortiguamiento pared
    double kn = 1.0e4;            // Rigidez normal entre partículas
    double gamma_n = 0.5;         // Amortiguamiento/Fricción entre partículas (igual a p.gamma_)
    double tau = 1.0e-6;          // Rampa instantánea para encendido desde t=0
    double radio_particula = 0.15;// Radio de la partícula [m]
    double dt = 1.0e-4;           // Paso de tiempo [s]
    double t_max = 100.0;          // Tiempo total de simulación [s]
    int save_every = 100;         // Guardado de frames cada N pasos
};

void generar_configuracion(std::vector<Particle>& particles, const Params& p) {
    std::mt19937 rng(42);
    
    // Distribuciones idénticas a la versión hamiltoniana
    double r_min = 0.25 * p.R_plato;
    double r_max = 0.85 * p.R_plato;
    std::uniform_real_distribution<double> U_r2(r_min * r_min, r_max * r_max);
    std::uniform_real_distribution<double> U_phi(0.0, 2.0 * std::acos(-1.0));
    std::uniform_real_distribution<double> U_vrad(-0.8, 0.8);
    std::uniform_real_distribution<double> U_omega(1.0, 3.0);
    std::uniform_real_distribution<double> U_mass(0.8, 1.2);

    // 1. Asignación espacial en la corona circular
    for (int i = 0; i < p.n; ++i) {
        particles[i].mass = U_mass(rng);
        particles[i].radius = p.radio_particula;

        double rho0 = std::sqrt(U_r2(rng));
        double phi0 = U_phi(rng);
        particles[i].r = {rho0 * std::cos(phi0), rho0 * std::sin(phi0)};

        // Asignación de velocidades iniciales radiales y angulares aleatorias
        double vrad0 = U_vrad(rng);
        double om0 = U_omega(rng);
        double vx0 = vrad0 * std::cos(phi0) - rho0 * om0 * std::sin(phi0);
        double vy0 = vrad0 * std::sin(phi0) + rho0 * om0 * std::cos(phi0);
        particles[i].v = {vx0, vy0};
    }

    // 2. Relajación de superposiciones iniciales
    double diametro = 2.0 * p.radio_particula;
    double diam_cuadrado = diametro * diametro;
    std::vector<Vec2D> desplazamientos(p.n);

    for (int iter = 0; iter < 20000; ++iter) {
        double max_superposicion = 0.0;
        std::fill(desplazamientos.begin(), desplazamientos.end(), Vec2D{0.0, 0.0});

        for (int i = 0; i < p.n - 1; ++i) {
            for (int j = i + 1; j < p.n; ++j) {
                Vec2D dr = particles[i].r - particles[j].r;
                double dist_sq = norm_sq(dr);
                if (dist_sq < diam_cuadrado && dist_sq > 0.0) {
                    double dist = std::sqrt(dist_sq);
                    double superposicion = diametro - dist;
                    max_superposicion = std::max(max_superposicion, superposicion);
                    Vec2D f = (dr / dist) * (superposicion * 0.5);
                    desplazamientos[i] += f;
                    desplazamientos[j] -= f;
                }
            }
        }

        for (int i = 0; i < p.n; ++i) {
            particles[i].r += desplazamientos[i] * 0.1;
        }

        if (max_superposicion < 1.0e-6) break;
    }
}

void excitacion_orbital_rampa(std::vector<Particle>& particles, double tiempo, const Params& p) {
    double omega = p.Omega;
    double factor_exp = std::exp(-tiempo / p.tau);

    double A_t = p.R * (1.0 - factor_exp * (1.0 + tiempo / p.tau));
    double v_A = p.R * (tiempo / (p.tau * p.tau)) * factor_exp;
    double a_A = (p.R / (p.tau * p.tau)) * factor_exp * (1.0 - tiempo / p.tau);

    double aceleracion_x = -(a_A * std::cos(omega * tiempo) - 2.0 * v_A * omega * std::sin(omega * tiempo) - A_t * omega * omega * std::cos(omega * tiempo));
    double aceleracion_y = -(a_A * std::sin(omega * tiempo) + 2.0 * v_A * omega * std::cos(omega * tiempo) - A_t * omega * omega * std::sin(omega * tiempo));
    Vec2D aceleracion_inercial = {aceleracion_x, aceleracion_y};

    for (auto& part : particles) {
        part.a += aceleracion_inercial;
    }
}

void contenedor_circular(std::vector<Particle>& particles, const Params& p) {
    for (auto& part : particles) {
        double d = norm(part.r);
        double delta = d + part.radius - p.R_plato;

        if (delta > 0.0 && d > 0.0) {
            Vec2D n_wall = part.r * (-1.0 / d);
            double vn = dot(part.v, n_wall);
            double Fn_mag = std::max(p.kn_wall * delta - p.gamma_wall * vn, 0.0);
            part.a += (n_wall * Fn_mag) / part.mass;
        }
    }
}

void contacto_particulas(std::vector<Particle>& particles, const Params& p) {
    int num_p = particles.size();
    for (int i = 0; i < num_p; ++i) {
        for (int j = i + 1; j < num_p; ++j) {
            Vec2D r_ij = particles[i].r - particles[j].r;
            double d = norm(r_ij);
            double suma_radios = particles[i].radius + particles[j].radius;

            if (d < suma_radios && d > 0.0) {
                double delta = suma_radios - d;
                Vec2D n_ij = r_ij / d;
                Vec2D v_ij = particles[i].v - particles[j].v;
                double vn = dot(v_ij, n_ij);

                double Fn_mag = std::max(p.kn * delta - p.gamma_n * vn, 0.0);
                Vec2D Fn_vec = n_ij * Fn_mag;

                particles[i].a += Fn_vec / particles[i].mass;
                particles[j].a -= Fn_vec / particles[j].mass;
            }
        }
    }
}

void fuerza_total(std::vector<Particle>& particles, double tiempo, const Params& p) {
    for (auto& part : particles) {
        part.a = {0.0, 0.0};
    }
    excitacion_orbital_rampa(particles, tiempo, p);
    contenedor_circular(particles, p);
    contacto_particulas(particles, p);
}

void velocity_verlet_step(std::vector<Particle>& particles, double dt, double tiempo, const Params& p) {
    for (auto& part : particles) {
        part.r += part.v * dt + part.a * (0.5 * dt * dt);
        part.v += part.a * (0.5 * dt);
    }

    fuerza_total(particles, tiempo + dt, p);

    for (auto& part : particles) {
        part.v += part.a * (0.5 * dt);
    }
}

void guardar_frame_xyz(std::ofstream& archivo, const std::vector<Particle>& particles, double tiempo, double radio_contenedor) {
    archivo << particles.size() + 1 << "\n";
    archivo << "Properties=species:S:1:pos:R:3:radius:R:1 Time=" << tiempo << "\n";
    for (const auto& part : particles) {
        archivo << "Granulo " << part.r.x << " " << part.r.y << " 0.0 " << part.radius << "\n";
    }
    archivo << "Contenedor 0.0 0.0 0.0 " << radio_contenedor << "\n";
}

// --------- Energía del sistema y gráfica ---------
double wall_potential_energy(const std::vector<Particle>& particles, const Params& p) {
    double V = 0.0;
    for (const auto& part : particles) {
        double d = norm(part.r);
        double delta = d + part.radius - p.R_plato;
        if (delta > 0.0) {
            V += 0.5 * p.kn_wall * delta * delta;
        }
    }
    return V;
}

double pair_potential_energy(const std::vector<Particle>& particles, const Params& p) {
    double V = 0.0;
    int n = static_cast<int>(particles.size());
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            Vec2D dr = particles[i].r - particles[j].r;
            double d = norm(dr);
            double suma_radios = particles[i].radius + particles[j].radius;

            if (d < suma_radios && d > 0.0) {
                double delta = suma_radios - d;
                V += 0.5 * p.kn * delta * delta;
            }
        }
    }
    return V;
}

double kinetic_energy(const std::vector<Particle>& particles) {
    double K = 0.0;
    for (const auto& part : particles) {
        K += 0.5 * part.mass * norm_sq(part.v);
    }
    return K;
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

// Tiempo de relajación usando envolvente de la energía total
double period_envelope_relaxation_time(const std::vector<double>& tvec,
                                       const std::vector<double>& x,
                                       double T_orb,
                                       double tol_rel = 0.02) {
    if (tvec.size() < 3 || x.size() != tvec.size()) return -1.0;

    double t0 = tvec.front();
    double t1 = tvec.back();
    int n_cycles = (int)std::floor((t1 - t0) / T_orb);

    if (n_cycles < 3) return -1.0;

    std::vector<double> cycle_max(n_cycles);
    std::vector<double> cycle_min(n_cycles);

    size_t idx = 0;
    for (int c = 0; c < n_cycles; ++c) {
        double lo = t0 + c * T_orb;
        double hi = t0 + (c + 1) * T_orb;

        double mx = -1e300;
        double mn = 1e300;

        while (idx < tvec.size() && tvec[idx] < hi) {
            if (tvec[idx] >= lo) {
                mx = std::max(mx, x[idx]);
                mn = std::min(mn, x[idx]);
            }
            ++idx;
        }

        cycle_max[c] = mx;
        cycle_min[c] = mn;
    }

    double amp0 = cycle_max[0] - cycle_min[0];
    double ref = cycle_max.back();

    for (int c = 0; c < n_cycles; ++c) {
        double amp = cycle_max[c] - cycle_min[c];
        if (std::fabs((amp / std::max(amp0, 1e-9)) - 1.0) < tol_rel) {
            return tvec[std::min((size_t)c, tvec.size() - 1)];
        }
    }

    if (!tvec.empty()) return tvec.back();
    return -1.0;
}

// tau = primer instante en que la oscilación entra en régimen estable
// Es decir: la amplitud por ciclo deja de decay/variar y se mantiene casi constante
double stable_oscillation_start_time(const std::vector<double>& tvec,
                                    const std::vector<double>& x,
                                    double T_orb,
                                    double tol_rel = 0.10,
                                    int required_cycles = 3) {
    if (tvec.size() < 6 || x.size() != tvec.size()) return -1.0;

    std::vector<double> cycle_amp;
    std::vector<double> cycle_start;

    double t0 = tvec.front();
    size_t idx = 0;

    while (idx < tvec.size()) {
        double lo = t0 + cycle_amp.size() * T_orb;
        double hi = lo + T_orb;

        double mx = -1e300;
        double mn = 1e300;

        while (idx < tvec.size() && tvec[idx] < hi) {
            if (tvec[idx] >= lo) {
                mx = std::max(mx, x[idx]);
                mn = std::min(mn, x[idx]);
            }
            ++idx;
        }

        if (mx > -1e250 && mn < 1e250) {
            cycle_amp.push_back(mx - mn);
            cycle_start.push_back(lo);
        }
    }

    if (cycle_amp.size() < (size_t)(required_cycles + 2)) return -1.0;

    std::vector<double> rel_change(cycle_amp.size() - 1);
    for (size_t i = 1; i < cycle_amp.size(); ++i) {
        double denom = std::max(std::max(cycle_amp[i], cycle_amp[i - 1]), 1e-12);
        rel_change[i - 1] = std::abs(cycle_amp[i] - cycle_amp[i - 1]) / denom;
    }

    for (size_t i = required_cycles; i < rel_change.size(); ++i) {
        bool stable = true;
        for (size_t j = i - required_cycles + 1; j <= i; ++j) {
            if (rel_change[j] > tol_rel) {
                stable = false;
                break;
            }
        }

        if (stable) {
            return cycle_start[i - required_cycles + 1];
        }
    }

    return tvec.back();
}

double cluster_angular_velocity(const std::vector<Particle>& particles) {
    double M = 0.0;
    Vec2D cm = {0.0, 0.0};

    for (const auto& part : particles) {
        M += part.mass;
        cm = cm + part.r * part.mass;
    }

    if (M <= 1e-12) return 0.0;
    cm = cm / M;

    double Lz = 0.0;
    double I = 0.0;

    for (const auto& part : particles) {
        Vec2D r_rel = part.r - cm;
        Vec2D v_rel = part.v;
        double r2 = norm_sq(r_rel);
        I += part.mass * r2;
        Lz += part.mass * (r_rel.x * v_rel.y - r_rel.y * v_rel.x);
    }

    if (std::abs(I) < 1e-12) return 0.0;
    return Lz / I;
}

void save_energy_svg(const std::vector<double>& t,
                     const std::vector<double>& Ke,
                     const std::vector<double>& Wp,
                     const std::vector<double>& Etot,
                     double tau,
                     const std::string& filename) {
    const int W = 1100, H = 650;
    const int left = 90, right = 40, top = 60, bottom = 80;

    if (t.empty() || Ke.empty() || Wp.empty() || Etot.empty()) return;
    if (t.size() != Ke.size() || t.size() != Wp.size() || t.size() != Etot.size()) return;

    double xmin = t.front();
    double xmax = t.back();

    double ymin = std::min({*std::min_element(Ke.begin(), Ke.end()),
                            *std::min_element(Wp.begin(), Wp.end()),
                            *std::min_element(Etot.begin(), Etot.end())});
    double ymax = std::max({*std::max_element(Ke.begin(), Ke.end()),
                            *std::max_element(Wp.begin(), Wp.end()),
                            *std::max_element(Etot.begin(), Etot.end())});

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
        svg << xmap(t[i]) << "," << ymap(Ke[i]) << " ";
    }
    svg << "\"/>\n";

    svg << "<polyline fill=\"none\" stroke=\"green\" stroke-width=\"2\" points=\"";
    for (size_t i = 0; i < t.size(); ++i) {
        svg << xmap(t[i]) << "," << ymap(Wp[i]) << " ";
    }
    svg << "\"/>\n";

    svg << "<polyline fill=\"none\" stroke=\"red\" stroke-width=\"2\" points=\"";
    for (size_t i = 0; i < t.size(); ++i) {
        svg << xmap(t[i]) << "," << ymap(Etot[i]) << " ";
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
    svg << "<text x=\"370\" y=\"60\" font-size=\"14\" fill=\"green\">Energía del plato</text>\n";

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
    std::vector<Particle> particles(p.n);
    generar_configuracion(particles, p);

    std::string archivo_salida = "giro_swirling_dem.xyz";
    std::ofstream fout(archivo_salida);

    int pasos = static_cast<int>(std::round(p.t_max / p.dt));
    std::cout << "Iniciando simulacion DEM en C++...\n";

    std::vector<double> tiempos, energy_kinetic, energy_plato, energy_total;
    std::vector<double> omega_plato, omega_cluster;

    std::ofstream energy_file("energia_sistema.csv");
    energy_file << "t,E_kinetica,E_plato,E_total,Omega_plato,Omega_cluster\n";

    for (int paso = 1; paso <= pasos; ++paso) {
        double t_actual = paso * p.dt;
        velocity_verlet_step(particles, p.dt, t_actual, p);

        if (paso % p.save_every == 0) {
            guardar_frame_xyz(fout, particles, t_actual, p.R_plato);

            double K = kinetic_energy(particles);
            double Vw = wall_potential_energy(particles, p);
            double Vp = pair_potential_energy(particles, p);
            double E_total = K + Vw + Vp;

            double omega_cluster_inst = cluster_angular_velocity(particles);

            tiempos.push_back(t_actual);
            energy_kinetic.push_back(K);
            energy_plato.push_back(Vw);
            energy_total.push_back(E_total);

            omega_plato.push_back(p.Omega);
            omega_cluster.push_back(omega_cluster_inst);

            energy_file << std::setprecision(12)
                       << t_actual << ","
                       << K << ","
                       << Vw << ","
                       << E_total << ","
                       << p.Omega << ","
                       << omega_cluster_inst << "\n";
        }
    }

    energy_file.close();

    double T_orb = 2.0 * std::acos(-1.0) / std::fabs(p.Omega);
    double tau = stable_oscillation_start_time(tiempos, energy_total, T_orb, 0.10, 3);

    save_energy_svg(tiempos, energy_kinetic, energy_plato, energy_total, tau, "energia_sistema.svg");
    save_frequency_svg(tiempos, omega_plato, omega_cluster, "frecuencia_cluster_vs_plato.svg");

    fout.close();

    std::cout << "Simulación terminada. Frames en " << archivo_salida << "\n";
    std::cout << "Energías guardadas en energia_sistema.csv y energia_sistema.svg\n";
    std::cout << "Frecuencias guardadas en frecuencia_cluster_vs_plato.svg\n";
    std::cout << "Tiempo de relajación tau = " << tau << " s\n";
    std::cout << "Omega_plato final = " << omega_plato.back() << " rad/s\n";
    std::cout << "Omega_cluster final = " << omega_cluster.back() << " rad/s\n";
    std::cout << "Error de acople = " << std::abs(omega_cluster.back() - omega_plato.back()) << " rad/s\n";

#ifdef _WIN32
    std::system("start energia_sistema.svg");
    std::system("start frecuencia_cluster_vs_plato.svg");
#else
    std::system("xdg-open energia_sistema.svg");
    std::system("xdg-open frecuencia_cluster_vs_plato.svg");
#endif

    return 0;
}