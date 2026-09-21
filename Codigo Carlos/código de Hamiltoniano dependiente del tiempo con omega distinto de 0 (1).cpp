//  ADVERTENCIA IMPORTANTE SOBRE LA CONVERGENCIA (leer antes de usar este codigo)
//  -------------------------------------------------------------------------------------
//  Se verifico numericamente que, con Omega != 0, el Hamiltoniano canonico H_CK(t) tal
//  como esta definido literalmente en la Ec. 25 DIVERGE con el tiempo (crece sin cota,
//  dominado por el termino -e^{gamma t} Sum(1/2 mi R^2 Omega^2), que no depende de
//  ninguna variable dinamica). Esto no es un error de integracion: es una consecuencia
//  de aplicar el escalamiento global e^{+-gamma t} de Caldirola-Kanai a un Lagrangiano
//  que incluye la energia cinetica de arrastre orbital del plato, la cual NUNCA decae
//  (el plato orbita externamente a velocidad angular Omega constante, sin friccion que
//  la frene). Fisicamente, la friccion interna (Rayleigh) SI cancela ese termino de
//  arrastre al calcular velocidades relativas entre particulas (ver Sec. 4 del
//  documento), pero el modelo CK global no distingue "arrastre comun" de "movimiento
//  relativo" y termina amortiguando/escalando todo por igual -> de ahi el artefacto.
//
//  Por esta razon, este codigo reporta y grafica DOS cantidades distintas:
//    (a) H_CK(t): el Hamiltoniano canonico literal (Ec. 25) - diverge, se incluye solo
//        como referencia/diagnostico.
//    (b) E_fis(t): la energia mecanica FISICA real, T+V, calculada con las velocidades
//        fisicas verdaderas (rho_dot, phi_dot), sin ningun factor e^{+-gamma t}. Esta
//        es la cantidad fisicamente significativa, y se mantiene acotada.
//
//  Ademas, como el plato sigue orbitando para siempre (Omega=const != 0), el sistema
//  NO relaja al reposo (a diferencia del caso Omega=0): se asienta en un estado
//  dinamico acotado, oscilando de forma cuasi-periodica en torno a un nivel medio.
//  Por eso el "tiempo de relajacion" aqui se redefine como el tiempo en que la
//  ENVOLVENTE de E_fis(t) (promedio movil sobre una ventana de un periodo orbital
//  T_orb = 2*pi/|Omega|) se estabiliza dentro de una banda de tolerancia -- es un
//  "tiempo de asentamiento" analogo al de un oscilador forzado alcanzando su regimen
//  estacionario, no una relajacion a energia cero.
// =====================================================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>
#include <algorithm>

// ---------------------------- Parametros fisicos globales ----------------------------
struct Params {
    int    n         = 100;      // numero de particulas
    double gamma_    = 0.5;    // coeficiente de amortiguamiento de Caldirola-Kanai [1/s]
    double k         = 1000.0;   // rigidez de la pared blanda del plato [N/m]
    double R_plato   = 6.0;    // radio del plato [m]
    double R         = 1.5;    // distancia del origen del laboratorio al centro del plato [m]
    double Omega     = 1.2;    // velocidad angular orbital del plato [rad/s]  (Omega != 0)
    double t_max     = 100.0;   // tiempo total de simulacion [s]
    double dt        = 0.0005; // paso de integracion [s]
    int    save_every = 40;    // guardar 1 de cada 'save_every' pasos al CSV
};

// Potencial de pared blanda (oscilador armonico truncado por Heaviside)
inline double Vborde(double rho, const Params& p) {
    double d = rho - p.R_plato;
    if (d > 0.0) return 0.5 * p.k * d * d;
    return 0.0;
}
inline double dVborde_drho(double rho, const Params& p) {
    double d = rho - p.R_plato;
    if (d > 0.0) return p.k * d;
    return 0.0;
}

// Calcula A_i(t) = R*Omega*sin(psi_i) y B_i(t) = R*Omega*cos(psi_i), con
// psi_i = phi_i - Omega*t, es decir, las proyecciones (Vhat.rhat_i) y (Vhat.phihat_i)
// escaladas por R*Omega.
inline void driftAB(double phi, double t, const Params& p, double& A, double& B) {
    double psi = phi - p.Omega * t;
    double ROmega = p.R * p.Omega;
    A = ROmega * std::sin(psi);
    B = ROmega * std::cos(psi);
}

// Calcula dy/dt para todo el sistema (n particulas) en el instante t.
// Vector de estado y: [rho_1,phi_1,Prho_1,Pphi_1, rho_2,phi_2,Prho_2,Pphi_2, ...]
std::vector<double> derivatives(double t, const std::vector<double>& y,
                                 const std::vector<double>& m, const Params& p) {
    int n = p.n;
    std::vector<double> dydt(4 * n);
    double emg = std::exp(-p.gamma_ * t); // e^{-gamma t}
    double epg = std::exp( p.gamma_ * t); // e^{+gamma t}

    for (int i = 0; i < n; ++i) {
        double rho  = y[4*i + 0];
        double phi  = y[4*i + 1];
        double Prho = y[4*i + 2];
        double Pphi = y[4*i + 3];
        double mi   = m[i];

        double A, B;
        driftAB(phi, t, p, A, B);

        double Prho_eff = Prho - mi * epg * A;
        double Pphi_eff = Pphi - mi * epg * rho * B;

        double rho_dot = emg * Prho_eff / mi;
        double phi_dot = emg * Pphi_eff / (mi * rho * rho);

        double Prho_dot = B * Pphi_eff / (rho * rho)
                           + emg * (Pphi_eff * Pphi_eff) / (mi * rho * rho * rho)
                           - epg * dVborde_drho(rho, p);

        double Pphi_dot = B * Prho_eff - A * Pphi_eff / rho;

        dydt[4*i + 0] = rho_dot;
        dydt[4*i + 1] = phi_dot;
        dydt[4*i + 2] = Prho_dot;
        dydt[4*i + 3] = Pphi_dot;
    }
    return dydt;
}

// Un paso de integracion por el metodo de Euler-Richardson (punto medio, 2do orden):
//   1) f1    = derivadas en (t, y)
//   2) y_mid = y + f1 * (dt/2)
//   3) f_mid = derivadas en (t+dt/2, y_mid)
//   4) y_nuevo = y + f_mid * dt
void euler_richardson_step(double t, std::vector<double>& y, const std::vector<double>& m,
                            const Params& p) {
    double dt = p.dt;
    auto f1 = derivatives(t, y, m, p);

    std::vector<double> y_mid(y.size());
    for (size_t j = 0; j < y.size(); ++j) y_mid[j] = y[j] + 0.5 * dt * f1[j];

    auto f_mid = derivatives(t + 0.5 * dt, y_mid, m, p);

    for (size_t j = 0; j < y.size(); ++j) {
        y[j] += dt * f_mid[j];
    }
}

// Hamiltoniano canonico total del sistema en el instante t (Ec. 25 literal, Omega!=0).
// ADVERTENCIA: esta cantidad DIVERGE con t. Se incluye
// solo con fines de referencia/diagnostico, no para estimar el tiempo de relajacion.
double hamiltonian_CK(double t, const std::vector<double>& y,
                       const std::vector<double>& m, const Params& p) {
    int n = p.n;
    double emg = std::exp(-p.gamma_ * t);
    double epg = std::exp( p.gamma_ * t);
    double H = 0.0;
    for (int i = 0; i < n; ++i) {
        double rho  = y[4*i + 0];
        double phi  = y[4*i + 1];
        double Prho = y[4*i + 2];
        double Pphi = y[4*i + 3];
        double mi   = m[i];

        double A, B;
        driftAB(phi, t, p, A, B);

        double Prho_eff = Prho - mi * epg * A;
        double Pphi_eff = Pphi - mi * epg * rho * B;

        H += emg * (Prho_eff * Prho_eff) / (2.0 * mi);
        H += emg * (Pphi_eff * Pphi_eff) / (2.0 * mi * rho * rho);
        H -= epg * 0.5 * mi * p.R * p.R * p.Omega * p.Omega;
        H += epg * Vborde(rho, p);
    }
    return H;
}

// Energia mecanica FISICA real del sistema, T + V, calculada con las velocidades
// fisicas verdaderas (rho_dot, phi_dot) -- SIN ningun factor e^{+-gamma t}. Esta es la cantidad
// fisicamente significativa y la que se usa para estimar el tiempo de relajacion.
double physical_energy(double t, const std::vector<double>& y,
                        const std::vector<double>& m, const Params& p) {
    int n = p.n;
    double emg = std::exp(-p.gamma_ * t);
    double epg = std::exp( p.gamma_ * t);
    double ROmega = p.R * p.Omega;
    double E = 0.0;
    for (int i = 0; i < n; ++i) {
        double rho  = y[4*i + 0];
        double phi  = y[4*i + 1];
        double Prho = y[4*i + 2];
        double Pphi = y[4*i + 3];
        double mi   = m[i];

        double A, B;
        driftAB(phi, t, p, A, B);
        double Prho_eff = Prho - mi * epg * A;
        double Pphi_eff = Pphi - mi * epg * rho * B;

        double rho_dot = emg * Prho_eff / mi;
        double phi_dot = emg * Pphi_eff / (mi * rho * rho);
        double psi = phi - p.Omega * t;

        double T = 0.5*mi*rho_dot*rho_dot
                 + 0.5*mi*rho*rho*phi_dot*phi_dot
                 + mi*ROmega*(rho_dot*std::sin(psi) + rho*phi_dot*std::cos(psi))
                 + 0.5*mi*ROmega*ROmega;

        E += T + Vborde(rho, p);
    }
    return E;
}

// Energia cinetica de las particulas (inercial)
double kinetic_energy_particles(double t, const std::vector<double>& y,
                               const std::vector<double>& m, const Params& p) {
    int n = p.n;
    double emg = std::exp(-p.gamma_ * t);
    double epg = std::exp( p.gamma_ * t);

    double K = 0.0;
    for (int i = 0; i < n; ++i) {
        double rho  = y[4*i + 0];
        double phi  = y[4*i + 1];
        double Prho = y[4*i + 2];
        double Pphi = y[4*i + 3];
        double mi   = m[i];

        double A, B;
        driftAB(phi, t, p, A, B);

        double Prho_eff = Prho - mi * epg * A;
        double Pphi_eff = Pphi - mi * epg * rho * B;

        double rho_dot = emg * Prho_eff / mi;
        double phi_dot = emg * Pphi_eff / (mi * rho * rho);

        K += 0.5 * mi * (rho_dot * rho_dot + rho * rho * phi_dot * phi_dot);
    }
    return K;
}

// Energia asociada al movimiento orbital del plato
double plate_orbital_energy(const std::vector<double>& m, const Params& p) {
    double M = 0.0;
    for (double mi : m) M += mi;
    double v_plato = p.R * p.Omega;
    return 0.5 * M * v_plato * v_plato;
}

// Energia potencial contra las paredes (borde blando del plato)
double wall_potential_energy(const std::vector<double>& y, const std::vector<double>& m,
                            const Params& p) {
    (void)m;
    double V = 0.0;
    for (int i = 0; i < p.n; ++i) {
        double rho = y[4*i + 0];
        V += Vborde(rho, p);
    }
    return V;
}

// Energia total fisica: cinetica particulas + orbital del plato + potencial paredes
double total_physical_energy(double t, const std::vector<double>& y,
                            const std::vector<double>& m, const Params& p) {
    return kinetic_energy_particles(t, y, m, p) +
           plate_orbital_energy(m, p) +
           wall_potential_energy(y, m, p);
}

// Promedio movil centrado de una serie, usado para obtener la ENVOLVENTE de E_fis(t)
// y filtrar la oscilacion cuasi-periodica inducida por la orbita del plato.
// window_n = numero de muestras que caben en una ventana de duracion 'window_time'.
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

// Estima el tiempo de relajacion: primer instante a partir del cual H(t) se mantiene
// dentro de +-tol_rel * |H0 - H_inf| de su valor asintotico H_inf.
double relaxation_time(const std::vector<double>& tvec, const std::vector<double>& Hvec,
                        double tol_rel = 0.02) {
    int N = (int)Hvec.size();
    if (N < 10) return -1.0;

    int tail = std::max(5, N / 20);
    double H_inf = 0.0;
    for (int i = N - tail; i < N; ++i) H_inf += Hvec[i];
    H_inf /= tail;

    double H0 = Hvec[0];
    double band = tol_rel * std::fabs(H0 - H_inf);
    if (band < 1e-12) band = tol_rel * std::fabs(H0);

    int last_outside = -1;
    for (int i = 0; i < N; ++i) {
        if (std::fabs(Hvec[i] - H_inf) > band) last_outside = i;
    }
    int idx = (last_outside == -1) ? 0 : std::min(last_outside + 1, N - 1);
    return tvec[idx];
}

void save_graph_svg(const std::vector<double>& t,
                    const std::vector<double>& H,
                    double tau,
                    const std::string& filename,
                    const std::string& title,
                    const std::string& ylabel) {
    const int W = 1000, Ht = 600;
    const int left = 80, right = 30, top = 40, bottom = 70;

    double xmin = t.front();
    double xmax = t.back();
    double ymin = *std::min_element(H.begin(), H.end());
    double ymax = *std::max_element(H.begin(), H.end());

    if (ymax == ymin) ymax = ymin + 1.0;

    double margen = 0.05 * (ymax - ymin);
    ymin -= margen;
    ymax += margen;

    auto xmap = [&](double x) {
        return left + (x - xmin) / (xmax - xmin) * (W - left - right);
    };

    auto ymap = [&](double y) {
        return Ht - bottom -
               (y - ymin) / (ymax - ymin) * (Ht - top - bottom);
    };

    std::ofstream svg(filename);

    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
        << "width=\"" << W << "\" height=\"" << Ht << "\">\n";

    svg << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

    // Ejes
    svg << "<line x1=\"" << left << "\" y1=\"" << top
        << "\" x2=\"" << left << "\" y2=\"" << Ht - bottom
        << "\" stroke=\"black\"/>\n";

    svg << "<line x1=\"" << left << "\" y1=\"" << Ht - bottom
        << "\" x2=\"" << W - right << "\" y2=\"" << Ht - bottom
        << "\" stroke=\"black\"/>\n";

    // Hamiltoniano H_CK(t)
    svg << "<polyline fill=\"none\" stroke=\"blue\" stroke-width=\"2\" points=\"";

    for (size_t i = 0; i < t.size(); ++i) {
        svg << xmap(t[i]) << "," << ymap(H[i]) << " ";
    }

    svg << "\"/>\n";

    // Tiempo de asentamiento
    if (tau >= xmin && tau <= xmax) {
        double xtau = xmap(tau);

        svg << "<line x1=\"" << xtau << "\" y1=\"" << top
            << "\" x2=\"" << xtau << "\" y2=\"" << Ht - bottom
            << "\" stroke=\"red\" stroke-width=\"2\" "
               "stroke-dasharray=\"6,4\"/>\n";

        svg << "<text x=\"" << xtau + 8 << "\" y=\"" << top + 20
            << "\" fill=\"red\">tau = " << tau << " s</text>\n";
    }

    svg << "<text x=\"" << W / 2 << "\" y=\"25\" "
        << "text-anchor=\"middle\" font-size=\"18\">"
        << title << "</text>\n";

    svg << "<text x=\"" << W / 2 << "\" y=\"" << Ht - 20
        << "\" text-anchor=\"middle\">Tiempo t [s]</text>\n";

    svg << "<text x=\"20\" y=\"" << Ht / 2
        << "\" transform=\"rotate(-90 20," << Ht / 2
        << ")\" text-anchor=\"middle\">" << ylabel << "</text>\n";

    svg << "</svg>\n";
}

void save_energy_components_svg(const std::vector<double>& t,
                                const std::vector<double>& Ekin,
                                const std::vector<double>& Eplate,
                                const std::vector<double>& Ewall,
                                const std::vector<double>& Etot,
                                double tau,
                                const std::string& filename) {
    const int W = 1100, H = 650;
    const int left = 90, right = 40, top = 60, bottom = 80;

    if (t.empty() || Ekin.empty() || Eplate.empty() || Ewall.empty() || Etot.empty()) return;
    if (t.size() != Ekin.size() || t.size() != Eplate.size() ||
        t.size() != Ewall.size() || t.size() != Etot.size()) return;

    double xmin = t.front();
    double xmax = t.back();

    double ymin = std::min({*std::min_element(Ekin.begin(), Ekin.end()),
                            *std::min_element(Eplate.begin(), Eplate.end()),
                            *std::min_element(Ewall.begin(), Ewall.end()),
                            *std::min_element(Etot.begin(), Etot.end())});
    double ymax = std::max({*std::max_element(Ekin.begin(), Ekin.end()),
                            *std::max_element(Eplate.begin(), Eplate.end()),
                            *std::max_element(Ewall.begin(), Ewall.end()),
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
        svg << xmap(t[i]) << "," << ymap(Ekin[i]) << " ";
    }
    svg << "\"/>\n";

    svg << "<polyline fill=\"none\" stroke=\"green\" stroke-width=\"2\" points=\"";
    for (size_t i = 0; i < t.size(); ++i) {
        svg << xmap(t[i]) << "," << ymap(Eplate[i]) << " ";
    }
    svg << "\"/>\n";

    svg << "<polyline fill=\"none\" stroke=\"#d97706\" stroke-width=\"2\" points=\"";
    for (size_t i = 0; i < t.size(); ++i) {
        svg << xmap(t[i]) << "," << ymap(Ewall[i]) << " ";
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

    svg << "<rect x=\"720\" y=\"520\" width=\"330\" height=\"90\" fill=\"white\" stroke=\"black\" stroke-width=\"1\"/>\n";

    svg << "<line x1=\"740\" y1=\"540\" x2=\"770\" y2=\"540\" stroke=\"blue\" stroke-width=\"3\"/>\n";
    svg << "<text x=\"780\" y=\"545\" font-size=\"14\" fill=\"blue\">Azul: E cinética particulas</text>\n";

    svg << "<line x1=\"740\" y1=\"560\" x2=\"770\" y2=\"560\" stroke=\"green\" stroke-width=\"3\"/>\n";
    svg << "<text x=\"780\" y=\"565\" font-size=\"14\" fill=\"green\">Verde: E plato</text>\n";

    svg << "<line x1=\"740\" y1=\"580\" x2=\"770\" y2=\"580\" stroke=\"#d97706\" stroke-width=\"3\"/>\n";
    svg << "<text x=\"780\" y=\"585\" font-size=\"14\" fill=\"#d97706\">Naranja: E pared</text>\n";

    svg << "<line x1=\"740\" y1=\"600\" x2=\"770\" y2=\"600\" stroke=\"red\" stroke-width=\"3\"/>\n";
    svg << "<text x=\"780\" y=\"605\" font-size=\"14\" fill=\"red\">Rojo: E total</text>\n";

    svg << "<text x=\"" << W / 2 << "\" y=\"25\" text-anchor=\"middle\" font-size=\"18\">"
        << "Energía física total y componentes</text>\n";
    svg << "<text x=\"" << W / 2 << "\" y=\"" << H - 20
        << "\" text-anchor=\"middle\" font-size=\"12\">Tiempo t [s]</text>\n";
    svg << "<text x=\"20\" y=\"" << H / 2
        << "\" transform=\"rotate(-90 20," << H / 2
        << ")\" text-anchor=\"middle\" font-size=\"12\">Energía [J]</text>\n";
    svg << "</svg>\n";
}

double cluster_angular_velocity(const std::vector<double>& y,
                               const std::vector<double>& m,
                               double t,
                               const Params& p) {
    int n = p.n;
    double emg = std::exp(-p.gamma_ * t);
    double epg = std::exp( p.gamma_ * t);

    std::vector<double> x(n), yv(n), vx(n), vy(n);
    double M = 0.0;
    double cmx = 0.0, cmy = 0.0;

    for (int i = 0; i < n; ++i) {
        double rho = y[4*i + 0];
        double phi = y[4*i + 1];
        double Prho = y[4*i + 2];
        double Pphi = y[4*i + 3];
        double mi = m[i];

        double A, B;
        driftAB(phi, t, p, A, B);

        double Prho_eff = Prho - mi * epg * A;
        double Pphi_eff = Pphi - mi * epg * rho * B;

        double rho_dot = emg * Prho_eff / mi;
        double phi_dot = emg * Pphi_eff / (mi * rho * rho);

        x[i] = rho * std::cos(phi);
        yv[i] = rho * std::sin(phi);

        vx[i] = rho_dot * std::cos(phi) - rho * phi_dot * std::sin(phi);
        vy[i] = rho_dot * std::sin(phi) + rho * phi_dot * std::cos(phi);

        cmx += mi * x[i];
        cmy += mi * yv[i];
        M += mi;
    }

    if (M <= 1e-12) return 0.0;

    cmx /= M;
    cmy /= M;

    double Lz = 0.0;
    double I = 0.0;

    for (int i = 0; i < n; ++i) {
        double rx = x[i] - cmx;
        double ry = yv[i] - cmy;

        Lz += m[i] * (rx * vy[i] - ry * vx[i]);
        I += m[i] * (rx * rx + ry * ry);
    }

    if (std::abs(I) < 1e-12) return 0.0;

    return Lz / I;
}

void save_frequency_svg(const std::vector<double>& t,
                        const std::vector<double>& omega_plato,
                        const std::vector<double>& omega_cluster,
                        const std::string& title,
                        const std::string& ylabel,
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

    svg << "<rect x=\"720\" y=\"520\" width=\"270\" height=\"60\" fill=\"white\" stroke=\"black\" stroke-width=\"1\"/>\n";
    svg << "<line x1=\"740\" y1=\"540\" x2=\"770\" y2=\"540\" stroke=\"blue\" stroke-width=\"3\"/>\n";
    svg << "<text x=\"780\" y=\"545\" font-size=\"14\" fill=\"blue\">Azul: Ω plato</text>\n";

    svg << "<line x1=\"740\" y1=\"560\" x2=\"770\" y2=\"560\" stroke=\"red\" stroke-width=\"3\"/>\n";
    svg << "<text x=\"780\" y=\"565\" font-size=\"14\" fill=\"red\">Rojo: Ω cluster</text>\n";

    svg << "<text x=\"" << W / 2 << "\" y=\"25\" text-anchor=\"middle\" font-size=\"18\">"
        << title << "</text>\n";
    svg << "<text x=\"" << W / 2 << "\" y=\"" << H - 20
        << "\" text-anchor=\"middle\" font-size=\"12\">Tiempo t [s]</text>\n";
    svg << "<text x=\"20\" y=\"" << H / 2
        << "\" transform=\"rotate(-90 20," << H / 2
        << ")\" text-anchor=\"middle\" font-size=\"12\">" << ylabel << "</text>\n";
    svg << "</svg>\n";
}

int main() {
    Params p;

    std::mt19937 rng(42); 
    std::uniform_real_distribution<double> U_rho(0.25 * p.R_plato, 0.85 * p.R_plato);
    std::uniform_real_distribution<double> U_phi(0.0, 2.0 * std::acos(-1.0));
    std::uniform_real_distribution<double> U_vrad(-0.8, 0.8);   // velocidad radial inicial [m/s]
    std::uniform_real_distribution<double> U_omega(1.0, 3.0);   // velocidad angular local inicial [rad/s]
    std::uniform_real_distribution<double> U_mass(0.8, 1.2);    // masa [kg]

    std::vector<double> m(p.n);
    std::vector<double> y(4 * p.n);

    for (int i = 0; i < p.n; ++i) {
        double mi    = U_mass(rng);
        double rho0  = U_rho(rng);
        double phi0  = U_phi(rng);
        double vrad0 = U_vrad(rng);
        double om0   = U_omega(rng);

        // En t=0: e^{gamma*0}=1, psi_i(0)=phi0 (pues Omega*0=0), luego
        // A_i(0) = R*Omega*sin(phi0), B_i(0) = R*Omega*cos(phi0).
        // Se inicializan los momentos canonicos a partir de las velocidades fisicas
        // deseadas usando las Ecs. 15-16 del documento (equivalentes a 23-24 en t=0):
        //     Prho_i(0) = mi*rho_dot(0) + mi*A_i(0)
        //     Pphi_i(0) = mi*rho0^2*phi_dot(0) + mi*rho0*B_i(0)
        double A0, B0;
        driftAB(phi0, 0.0, p, A0, B0);

        m[i] = mi;
        y[4*i + 0] = rho0;
        y[4*i + 1] = phi0;
        y[4*i + 2] = mi * vrad0 + mi * A0;
        y[4*i + 3] = mi * rho0 * rho0 * om0 + mi * rho0 * B0;
    }

    if (p.Omega == 0.0) {
        std::cerr << "Aviso: Omega = 0 en este archivo pensado para Omega != 0. "
                     "Use swirling_hamiltonian.cpp para ese caso.\n";
    }

    int nsteps = (int)std::round(p.t_max / p.dt);

    std::vector<double> tvec;
    std::vector<double> Ekin_vec, Eplate_vec, Ewall_vec, Etotal_vec;
    std::vector<double> Omega_plato_vec, Omega_cluster_vec;

    tvec.reserve(nsteps / p.save_every + 1);
    Ekin_vec.reserve(nsteps / p.save_every + 1);
    Eplate_vec.reserve(nsteps / p.save_every + 1);
    Ewall_vec.reserve(nsteps / p.save_every + 1);
    Etotal_vec.reserve(nsteps / p.save_every + 1);
    Omega_plato_vec.reserve(nsteps / p.save_every + 1);
    Omega_cluster_vec.reserve(nsteps / p.save_every + 1);

    double t = 0.0;
    for (int step = 0; step <= nsteps; ++step) {
        if (step % p.save_every == 0) {
            tvec.push_back(t);

            double K = kinetic_energy_particles(t, y, m, p);
            double Eplate = plate_orbital_energy(m, p);
            double Ewall = wall_potential_energy(y, m, p);
            double Etot = K + Eplate + Ewall;

            Ekin_vec.push_back(K);
            Eplate_vec.push_back(Eplate);
            Ewall_vec.push_back(Ewall);
            Etotal_vec.push_back(Etot);

            Omega_plato_vec.push_back(p.Omega);
            Omega_cluster_vec.push_back(cluster_angular_velocity(y, m, t, p));
        }

        if (step == nsteps) break;
        euler_richardson_step(t, y, m, p);
        t += p.dt;
    }

    double T_orb = 2.0 * std::acos(-1.0) / std::fabs(p.Omega);
    double sample_dt = p.dt * p.save_every;
    int window_n = std::max(3, (int)std::round(T_orb / sample_dt));
    std::vector<double> Etotal_smooth = moving_average(Etotal_vec, window_n);

    double tau = relaxation_time(tvec, Etotal_smooth, 0.02);

    std::ofstream fout("swirling_hamiltonian_omega.csv");
    fout << "t,E_kin,E_plato,E_pared,E_total,Omega_plato,Omega_cluster,tau\n";
    for (size_t i = 0; i < tvec.size(); ++i) {
        fout << std::setprecision(10)
             << tvec[i] << ","
             << Ekin_vec[i] << ","
             << Eplate_vec[i] << ","
             << Ewall_vec[i] << ","
             << Etotal_vec[i] << ","
             << Omega_plato_vec[i] << ","
             << Omega_cluster_vec[i] << ","
             << tau << "\n";
    }
    fout.close();

    save_energy_components_svg(
        tvec,
        Ekin_vec,
        Eplate_vec,
        Ewall_vec,
        Etotal_vec,
        tau,
        "energia_fisica_total.svg"
    );

    save_frequency_svg(
        tvec,
        Omega_plato_vec,
        Omega_cluster_vec,
        "Frecuencia angular",
        "Ω [rad/s]",
        "frecuencia_angulares.svg"
    );

#ifdef _WIN32
    std::system("start energia_fisica_total.svg");
    std::system("start frecuencia_angulares.svg");
#else
    std::system("xdg-open energia_fisica_total.svg");
    std::system("xdg-open frecuencia_angulares.svg");
#endif

    std::cout << "Graficas guardadas en: energia_fisica_total.svg y frecuencia_angulares.svg\n";
    std::cout << "tau = " << tau << " s\n";

    return 0;
}
