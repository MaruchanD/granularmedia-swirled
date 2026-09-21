// =====================================================================================
//  swirling_rayleigh_omega.cpp
//
//  Version MAS GENERAL de la simulacion: Omega != 0 (el plato orbita el laboratorio)
//  Y disipacion introducida de forma EXPLICITA mediante la funcion de Rayleigh F
//  (Ec. 18), 
//  Con psi_i = phi_i - Omega t,  A_i = R Omega sin(psi_i) = R Omega (Vhat.rhat_i),
//  B_i = R Omega cos(psi_i) = R Omega (Vhat.phihat_i):
//
//      H_i = (prho_i - m_i A_i)^2 / (2 m_i)
//          + (pphi_i - m_i rho_i B_i)^2 / (2 m_i rho_i^2)
//          - (1/2) m_i R^2 Omega^2
//          + Vborde(rho_i)
//
//  A diferencia de la version con Caldirola-Kanai, AQUI NO HAY factores e^{+-gamma t}:
//  el termino "-(1/2) m_i R^2 Omega^2" es una simple CONSTANTE (no diverge con t).
//  Sin embargo, H SIGUE siendo explicitamente dependiente del tiempo, porque A_i y B_i
//  dependen de psi_i = phi_i - Omega t: esto es una dependencia temporal FISICA real
//  (el plato orbita, su orientacion relativa a cada particula cambia con t), no un
//  artefacto matematico como en Caldirola-Kanai.
//
// -------------------------------------------------------------------------------------
//  2) FUNCION DE DISIPACION DE RAYLEIGH (Ec. 18) -- SIN CAMBIOS respecto al caso Omega=0
// -------------------------------------------------------------------------------------
//  La velocidad de arrastre ROmega*Vhat es COMUN a todas las
//  particulas (mismo plato), por lo que se CANCELA idénticamente al calcular
//  velocidades relativas entre pares. Por eso F depende solo de la velocidad LOCAL
//  vrel_i = rho_dot_i rhat_i + rho_i phi_dot_i phihat_i, exactamente como en el caso
//  Omega=0:
//
//      F = (gamma/4) Sum_i Sum_j |vrel_i - vrel_j|^2
//
//      dF/d(rho_dot_i) = gamma * Sum_{j!=i} (vrel_i-vrel_j).rhat_i
//      dF/d(phi_dot_i) = gamma * rho_i * Sum_{j!=i} (vrel_i-vrel_j).phihat_i
//
//  (formulas validadas numericamente por diferencias finitas;
//  no dependen de R ni de Omega, solo de rho_i, phi_i, rho_dot_i, phi_dot_i).
//
// -------------------------------------------------------------------------------------
//  3) ECUACIONES DE MOVIMIENTO COMPLETAS (Omega != 0 + Rayleigh, sin CK)
// -------------------------------------------------------------------------------------
//
//      rho_dot_i = Prho_eff_i / m_i
//      phi_dot_i = Pphi_eff_i / (m_i rho_i^2)
//
//      prho_dot_i =  B_i Pphi_eff_i / rho_i^2  +  Pphi_eff_i^2/(m_i rho_i^3)
//                    - k(rho_i-R_PLATO)*[rho_i>R_PLATO]
//                    - dF/d(rho_dot_i)
//
//      pphi_dot_i =  B_i Prho_eff_i  -  A_i Pphi_eff_i / rho_i
//                    - dF/d(phi_dot_i)
//
//  donde Prho_eff_i = prho_i - m_i A_i ,  Pphi_eff_i = pphi_i - m_i rho_i B_i.
//
// -------------------------------------------------------------------------------------
//  4) QUE SE GRAFICA Y COMO SE MIDE LA RELAJACION
// -------------------------------------------------------------------------------------
//  Igual que con Caldirola-Kanai + Omega!=0, el Hamiltoniano CANONICO H(t) no es la
//  energia fisica real (por la dependencia explicita del tiempo). Por eso se calculan
//  y grafican DOS cantidades:
//    - H_CK-like(t): el Hamiltoniano canonico H = Sum H_i de arriba (diagnostico).
//    - E_fis(t): la energia mecanica real T+V (Ec. 13, con el termino cruzado R*Omega
//      incluido), calculada con las velocidades fisicas verdaderas.
//  Como el plato sigue orbitando para siempre, el sistema no relaja al reposo sino a
//  un regimen dinamico periodico/cuasi-periodico. El "tiempo de relajacion" se estima
//  comparando, ciclo orbital a ciclo orbital, el PICO y el VALLE de E_fis(t) dentro
//  de cada periodo (envolvente gruesa) contra el promedio de los ultimos ciclos: es
//  mas robusto que comparar E(t) fase a fase, porque no se ve afectado por el "beat"
//  entre la frecuencia orbital y las frecuencias propias de cada particula.
//  Integracion: Euler-Richardson (punto medio, 2do orden).
//  Salida: swirling_rayleigh_omega.csv  (columnas: t,H_canonico,E_fis,E_fis_envolvente)
// =====================================================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>
#include <algorithm>
#include <string>

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

    // Eje vertical izquierdo
    svg << "<line x1=\"" << left << "\" y1=\"" << top
        << "\" x2=\"" << left << "\" y2=\"" << H - bottom
        << "\" stroke=\"black\" stroke-width=\"1.5\"/>\n";

    // Eje horizontal inferior
    svg << "<line x1=\"" << left << "\" y1=\"" << H - bottom
        << "\" x2=\"" << W - right << "\" y2=\"" << H - bottom
        << "\" stroke=\"black\" stroke-width=\"1.5\"/>\n";

    // Escala horizontal
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

    // Escala vertical
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

    // Curva principal
    svg << "<polyline fill=\"none\" stroke=\"blue\" stroke-width=\"2\" points=\"";
    for (size_t i = 0; i < t.size(); ++i)
        svg << xmap(t[i]) << "," << ymap(values[i]) << " ";
    svg << "\"/>\n";

    // Línea del tiempo de relajación
    if (tau >= xmin && tau <= xmax) {
        double xtau = xmap(tau);

        svg << "<line x1=\"" << xtau << "\" y1=\"" << top
            << "\" x2=\"" << xtau << "\" y2=\"" << H - bottom
            << "\" stroke=\"red\" stroke-width=\"2\" "
            << "stroke-dasharray=\"7,5\"/>\n";

        svg << "<text x=\"" << xtau + 8 << "\" y=\"" << top + 22
            << "\" fill=\"red\" font-size=\"12\">tau = " << tau << " s</text>\n";
    }

    // Título y etiquetas
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
    int    n         = 5;      // numero de particulas
    double gamma_    = 0.5;    // constante de acoplamiento de la friccion de Rayleigh
    double k         = 1000.0;   // rigidez de la pared blanda del plato [N/m]
    double R_plato   = 6.0;    // radio del plato [m]
    double R         = 1.5;    // distancia del origen del laboratorio al centro del plato [m]
    double Omega     = 1.2;    // velocidad angular orbital del plato [rad/s] (!= 0)
    double t_max     = 100.0;   // tiempo total de simulacion [s] (~7-8 periodos orbitales)
    double dt        = 0.0005; // paso de integracion [s]
    int    save_every = 40;    // guardar 1 de cada 'save_every' pasos
};

inline double Vborde(double rho, const Params& p) {
    double d = rho - p.R_plato;
    return (d > 0.0) ? 0.5*p.k*d*d : 0.0;
}
inline double dVborde_drho(double rho, const Params& p) {
    double d = rho - p.R_plato;
    return (d > 0.0) ? p.k*d : 0.0;
}

// A_i(t) = R*Omega*sin(psi_i), B_i(t) = R*Omega*cos(psi_i), psi_i = phi_i - Omega*t
inline void driftAB(double phi, double t, const Params& p, double& A, double& B) {
    double psi = phi - p.Omega * t;
    double ROmega = p.R * p.Omega;
    A = ROmega * std::sin(psi);
    B = ROmega * std::cos(psi);
}

// Calcula, para todas las particulas, rho_dot_i, phi_dot_i (velocidades fisicas
// locales) y sus componentes cartesianas locales vx_i,vy_i (usadas por F).
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
        vx[i] = rho_dot[i]*std::cos(phi) - rho*phi_dot[i]*std::sin(phi);
        vy[i] = rho_dot[i]*std::sin(phi) + rho*phi_dot[i]*std::cos(phi);
    }
}

std::vector<double> derivatives(double t, const std::vector<double>& y,
                                 const std::vector<double>& m, const Params& p) {
    int n = p.n;
    std::vector<double> dydt(4*n);
    std::vector<double> rho_dot(n), phi_dot(n), vx(n), vy(n), Prho_eff(n), Pphi_eff(n);
    local_state(t, y, m, p, rho_dot, phi_dot, vx, vy, Prho_eff, Pphi_eff);

    for (int i = 0; i < n; ++i) {
        double rho = y[4*i+0], phi = y[4*i+1], mi = m[i];
        double A, B;
        driftAB(phi, t, p, A, B);
        double cphi = std::cos(phi), sphi = std::sin(phi);

        // Fuerza generalizada disipativa (Rayleigh), acoplando con las demas particulas
        double sum_r = 0.0, sum_phi = 0.0;
        for (int j = 0; j < n; ++j) {
            if (j == i) continue;
            double dvx = vx[i]-vx[j], dvy = vy[i]-vy[j];
            sum_r   += dvx*cphi + dvy*sphi;
            sum_phi += -dvx*sphi + dvy*cphi;
        }
        double dF_drhodot = p.gamma_ * sum_r;
        double dF_dphidot = p.gamma_ * rho * sum_phi;

        double prho_dot = B*Pphi_eff[i]/(rho*rho) + (Pphi_eff[i]*Pphi_eff[i])/(mi*rho*rho*rho)
                           - dVborde_drho(rho, p) - dF_drhodot;
        double pphi_dot = B*Prho_eff[i] - A*Pphi_eff[i]/rho - dF_dphidot;

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
    for (size_t j = 0; j < y.size(); ++j) y_mid[j] = y[j] + 0.5*dt*f1[j];
    auto f_mid = derivatives(t + 0.5*dt, y_mid, m, p);
    for (size_t j = 0; j < y.size(); ++j) y[j] += dt*f_mid[j];
}

// Hamiltoniano canonico (Ec. 17, sin CK). Diagnostico: NO es la energia fisica real.
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
        H += Prho_eff*Prho_eff/(2*mi) + Pphi_eff*Pphi_eff/(2*mi*rho*rho)
             - 0.5*mi*p.R*p.R*p.Omega*p.Omega + Vborde(rho, p);
    }
    return H;
}

// Energia mecanica fisica real T+V (Ec. 13 con el termino cruzado R*Omega incluido)
double physical_energy(double t, const std::vector<double>& y,
                        const std::vector<double>& m, const Params& p) {
    int n = p.n;
    double ROmega = p.R * p.Omega;
    double E = 0.0;
    std::vector<double> rho_dot(n), phi_dot(n), vx(n), vy(n), Prho_eff(n), Pphi_eff(n);
    local_state(t, y, m, p, rho_dot, phi_dot, vx, vy, Prho_eff, Pphi_eff);
    for (int i = 0; i < n; ++i) {
        double rho = y[4*i+0], phi = y[4*i+1], mi = m[i];
        double psi = phi - p.Omega*t;
        double T = 0.5*mi*rho_dot[i]*rho_dot[i]
                 + 0.5*mi*rho*rho*phi_dot[i]*phi_dot[i]
                 + mi*ROmega*(rho_dot[i]*std::sin(psi) + rho*phi_dot[i]*std::cos(psi))
                 + 0.5*mi*ROmega*ROmega;
        E += T + Vborde(rho, p);
    }
    return E;
}

double kinetic_energy(double t, const std::vector<double>& y, const std::vector<double>& m, const Params& p) {
    int n = p.n;
    std::vector<double> rho_dot(n), phi_dot(n), vx(n), vy(n), Prho_eff(n), Pphi_eff(n);
    local_state(t, y, m, p, rho_dot, phi_dot, vx, vy, Prho_eff, Pphi_eff);

    double K = 0.0;
    for (int i = 0; i < n; ++i) {
        K += 0.5 * m[i] * (vx[i]*vx[i] + vy[i]*vy[i]);
    }
    return K;
}

double wall_energy(const std::vector<double>& y, const std::vector<double>& m, const Params& p) {
    (void)m;
    double Vw = 0.0;
    for (int i = 0; i < p.n; ++i) {
        double rho = y[4*i+0];
        double d = rho - p.R_plato;
        if (d > 0.0) Vw += 0.5 * p.k * d * d;
    }
    return Vw;
}

double cluster_angular_velocity(const std::vector<double>& y, const std::vector<double>& m, double t, const Params& p) {
    int n = p.n;
    std::vector<double> rho_dot(n), phi_dot(n), vx(n), vy(n), Prho_eff(n), Pphi_eff(n);
    local_state(t, y, m, p, rho_dot, phi_dot, vx, vy, Prho_eff, Pphi_eff);

    double M = 0.0;
    double cmx = 0.0, cmy = 0.0;
    for (int i = 0; i < n; ++i) {
        double rho = y[4*i+0], phi = y[4*i+1];
        cmx += m[i] * rho * std::cos(phi);
        cmy += m[i] * rho * std::sin(phi);
        M += m[i];
    }
    if (M <= 1e-12) return 0.0;
    cmx /= M; cmy /= M;

    double Lz = 0.0, I = 0.0;
    for (int i = 0; i < n; ++i) {
        double rho = y[4*i+0], phi = y[4*i+1];
        double x = rho * std::cos(phi);
        double yv = rho * std::sin(phi);
        double rx = x - cmx;
        double ry = yv - cmy;

        Lz += m[i] * (rx * vy[i] - ry * vx[i]);
        I += m[i] * (rx * rx + ry * ry);
    }
    if (std::abs(I) < 1e-12) return 0.0;
    return Lz / I;
}

std::vector<double> moving_average(const std::vector<double>& x, int window_n) {
    int N = (int)x.size();
    std::vector<double> out(N);
    int half = window_n / 2;
    for (int i = 0; i < N; ++i) {
        int lo = std::max(0, i-half), hi = std::min(N-1, i+half);
        double s = 0.0;
        for (int j = lo; j <= hi; ++j) s += x[j];
        out[i] = s/(hi-lo+1);
    }
    return out;
}

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
        y[4*i+0] = rho0;
        y[4*i+1] = phi0;
        y[4*i+2] = mi * vrad0 + mi * A0;
        y[4*i+3] = mi * rho0 * rho0 * om0 + mi * rho0 * B0;
    }

    int nsteps = (int)std::round(p.t_max / p.dt);

    std::vector<double> tvec;
    std::vector<double> Hvec, Ephys_vec;
    std::vector<double> Ekin_vec, Eplato_vec, Etotal_vec;
    std::vector<double> Omega_plato_vec, Omega_cluster_vec;

    tvec.reserve(nsteps / p.save_every + 1);
    Hvec.reserve(nsteps / p.save_every + 1);
    Ephys_vec.reserve(nsteps / p.save_every + 1);
    Ekin_vec.reserve(nsteps / p.save_every + 1);
    Eplato_vec.reserve(nsteps / p.save_every + 1);
    Etotal_vec.reserve(nsteps / p.save_every + 1);
    Omega_plato_vec.reserve(nsteps / p.save_every + 1);
    Omega_cluster_vec.reserve(nsteps / p.save_every + 1);

    double t = 0.0;
    for (int step = 0; step <= nsteps; ++step) {
        if (step % p.save_every == 0) {
            tvec.push_back(t);

            double K = kinetic_energy(t, y, m, p);
            double Vwall = wall_energy(y, m, p);
            double E_total = K + Vwall;

            Hvec.push_back(hamiltonian_canonical(t, y, m, p));
            Ephys_vec.push_back(physical_energy(t, y, m, p));
            Ekin_vec.push_back(K);
            Eplato_vec.push_back(Vwall);
            Etotal_vec.push_back(E_total);

            Omega_plato_vec.push_back(p.Omega);
            Omega_cluster_vec.push_back(cluster_angular_velocity(y, m, t, p));
        }

        if (step == nsteps) break;
        euler_richardson_step(t, y, m, p);
        t += p.dt;
    }

    double T_orb = 2.0 * std::acos(-1.0) / std::fabs(p.Omega);
    double tau = stable_oscillation_start_time(tvec, Etotal_vec, T_orb, 0.10, 3);

    std::ofstream fout("swirling_rayleigh_omega.csv");
    fout << "t,H_canonico,E_fis,E_kinetica,E_plato,E_total,Omega_plato,Omega_cluster,tau\n";
    for (size_t i = 0; i < tvec.size(); ++i) {
        fout << std::setprecision(10)
             << tvec[i] << ","
             << Hvec[i] << ","
             << Ephys_vec[i] << ","
             << Ekin_vec[i] << ","
             << Eplato_vec[i] << ","
             << Etotal_vec[i] << ","
             << Omega_plato_vec[i] << ","
             << Omega_cluster_vec[i] << ","
             << tau << "\n";
    }
    fout.close();

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Simulacion completada.\n";
    std::cout << "tau = " << tau << " s\n";
    std::cout << "Omega_plato final = " << Omega_plato_vec.back() << " rad/s\n";
    std::cout << "Omega_cluster final = " << Omega_cluster_vec.back() << " rad/s\n";
    std::cout << "Error de acople = " << std::abs(Omega_cluster_vec.back() - Omega_plato_vec.back()) << " rad/s\n";

    save_energy_svg(tvec, Ekin_vec, Eplato_vec, Etotal_vec, tau, "energia_sistema.svg");
    save_frequency_svg(tvec, Omega_plato_vec, Omega_cluster_vec, "frecuencia_cluster_vs_plato.svg");

#ifdef _WIN32
    std::system("start energia_sistema.svg");
    std::system("start frecuencia_cluster_vs_plato.svg");
#else
    std::system("xdg-open energia_sistema.svg");
    std::system("xdg-open frecuencia_cluster_vs_plato.svg");
#endif

    return 0;
}
