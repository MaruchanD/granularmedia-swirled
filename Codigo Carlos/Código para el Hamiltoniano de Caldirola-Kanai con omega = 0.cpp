#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>
#include <algorithm>
// ---------------------------- Parametros fisicos globales ----------------------------
struct Params {
    int    n        = 6;       // numero de particulas
    double gamma_    = 0.5;    // coeficiente de amortiguamiento de Caldirola-Kanai [1/s]
    double k         = 60.0;   // rigidez de la pared blanda del plato [N/m]
    double R_plato   = 1.0;    // radio del plato [m]
    double t_max     = 20.0;   // tiempo total de simulacion [s]
    double dt        = 0.0005; // paso de integracion RK4 [s]
    int    save_every = 40;    // guardar 1 de cada 'save_every' pasos al CSV
};

// Estado de una particula: (rho, phi, Prho, Pphi). Las masas se guardan aparte.
// El vector de estado global y tiene tamano 4*n, ordenado por particula.

// Potencial de pared blanda (oscilador armonico truncado por Heaviside)
inline double Vborde(double rho, const Params& p) {
    double d = rho - p.R_plato;
    if (d > 0.0) return 0.5 * p.k * d * d;
    return 0.0;
}

// Derivada del potencial de pared respecto a rho (fuerza generalizada, sin signo)
inline double dVborde_drho(double rho, const Params& p) {
    double d = rho - p.R_plato;
    if (d > 0.0) return p.k * d;
    return 0.0;
}

// Calcula dy/dt para todo el sistema (n particulas) en el instante t
std::vector<double> derivatives(double t, const std::vector<double>& y,
                                 const std::vector<double>& m, const Params& p) {
    int n = p.n;
    std::vector<double> dydt(4 * n);
    double emg = std::exp(-p.gamma_ * t); // e^{-gamma t}
    double epg = std::exp( p.gamma_ * t); // e^{+gamma t}

    for (int i = 0; i < n; ++i) {
        double rho  = y[4*i + 0];
        // double phi = y[4*i + 1]; // no se usa explicitamente (coordenada ciclica)
        double Prho = y[4*i + 2];
        double Pphi = y[4*i + 3];
        double mi   = m[i];

        double rho_dot  = emg * Prho / mi;
        double phi_dot   = emg * Pphi / (mi * rho * rho);
        double Prho_dot  = emg * (Pphi * Pphi) / (mi * rho * rho * rho)
                            - epg * dVborde_drho(rho, p);
        double Pphi_dot  = 0.0; // phi ciclica: momento angular local conservado

        dydt[4*i + 0] = rho_dot;
        dydt[4*i + 1] = phi_dot;
        dydt[4*i + 2] = Prho_dot;
        dydt[4*i + 3] = Pphi_dot;
    }
    return dydt;
}


void euler_richardson_step(double t, std::vector<double>& y, const std::vector<double>& m,
                            const Params& p) {
    double dt = p.dt;

    // Paso 1: pendiente en el instante inicial del intervalo
    auto f1 = derivatives(t, y, m, p);

    // Paso 2: estado estimado en el punto medio del intervalo temporal
    std::vector<double> y_mid(y.size());
    for (size_t j = 0; j < y.size(); ++j) y_mid[j] = y[j] + 0.5 * dt * f1[j];

    // Paso 3: pendiente evaluada en el punto medio (t+dt/2, y_mid)
    auto f_mid = derivatives(t + 0.5 * dt, y_mid, m, p);

    // Paso 4: se avanza el paso completo usando la pendiente del punto medio
    for (size_t j = 0; j < y.size(); ++j) {
        y[j] += dt * f_mid[j];
    }
}

// Hamiltoniano total del sistema en el instante t (Ec. 25 con Omega = 0)
double hamiltonian(double t, const std::vector<double>& y,
                    const std::vector<double>& m, const Params& p) {
    int n = p.n;
    double emg = std::exp(-p.gamma_ * t);
    double epg = std::exp( p.gamma_ * t);
    double H = 0.0;
    for (int i = 0; i < n; ++i) {
        double rho  = y[4*i + 0];
        double Prho = y[4*i + 2];
        double Pphi = y[4*i + 3];
        double mi   = m[i];

        H += emg * (Prho * Prho) / (2.0 * mi);
        H += emg * (Pphi * Pphi) / (2.0 * mi * rho * rho);
        H += epg * Vborde(rho, p);
    }
    return H;
}

// Estima el tiempo de relajacion: primer instante a partir del cual H(t) se mantiene
// dentro de +-tol_rel * H0 de su valor asintotico H_inf (criterio de "settling time").
double relaxation_time(const std::vector<double>& tvec, const std::vector<double>& Hvec,
                        double tol_rel = 0.02) {
    int N = (int)Hvec.size();
    if (N < 10) return -1.0;

    // Valor asintotico: promedio del ultimo 5% de la serie
    int tail = std::max(5, N / 20);
    double H_inf = 0.0;
    for (int i = N - tail; i < N; ++i) H_inf += Hvec[i];
    H_inf /= tail;

    double H0 = Hvec[0];
    double band = tol_rel * std::fabs(H0 - H_inf);
    if (band < 1e-12) band = tol_rel * std::fabs(H0); // fallback si H0 ~ H_inf

    // Recorremos desde el final hacia atras: buscamos el ultimo indice que se sale
    // de la banda. El tiempo de relajacion es el instante inmediatamente posterior.
    int last_outside = -1;
    for (int i = 0; i < N; ++i) {
        if (std::fabs(Hvec[i] - H_inf) > band) last_outside = i;
    }
    int idx = (last_outside == -1) ? 0 : std::min(last_outside + 1, N - 1);
    return tvec[idx];
}


void save_graph_svg(const std::vector<double>& t,
                    const std::vector<double>& H,
                    double tau) {
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
        return Ht - bottom - (y - ymin) / (ymax - ymin) * (Ht - top - bottom);
    };

    std::ofstream svg("relajacion_hamiltoniano.svg");

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

    // Curva H(t)
    svg << "<polyline fill=\"none\" stroke=\"blue\" stroke-width=\"2\" points=\"";

    for (size_t i = 0; i < t.size(); ++i) {
        svg << xmap(t[i]) << "," << ymap(H[i]) << " ";
    }

    svg << "\"/>\n";

    // Línea vertical del tiempo de relajación
    if (tau >= xmin && tau <= xmax) {
        double xtau = xmap(tau);

        svg << "<line x1=\"" << xtau << "\" y1=\"" << top
            << "\" x2=\"" << xtau << "\" y2=\"" << Ht - bottom
            << "\" stroke=\"red\" stroke-width=\"2\" stroke-dasharray=\"6,4\"/>\n";

        svg << "<text x=\"" << xtau + 8 << "\" y=\"" << top + 20
            << "\" fill=\"red\">tau = " << tau << " s</text>\n";
    }

    // Etiquetas
    svg << "<text x=\"" << W / 2 << "\" y=\"25\" text-anchor=\"middle\" "
        << "font-size=\"18\">Relajacion del Hamiltoniano</text>\n";

    svg << "<text x=\"" << W / 2 << "\" y=\"" << Ht - 20
        << "\" text-anchor=\"middle\">Tiempo t [s]</text>\n";

    svg << "<text x=\"20\" y=\"" << Ht / 2
        << "\" transform=\"rotate(-90 20," << Ht / 2
        << ")\" text-anchor=\"middle\">H(t)</text>\n";

    svg << "</svg>\n";
}



int main() {
    Params p;

    std::mt19937 rng(42); // semilla fija -> resultados reproducibles
    std::uniform_real_distribution<double> U_rho(0.25 * p.R_plato, 0.85 * p.R_plato);
    std::uniform_real_distribution<double> U_phi(
        0.0, 2.0 * std::acos(-1.0)
    );
    std::uniform_real_distribution<double> U_vrad(-0.8, 0.8);   // velocidad radial inicial [m/s]
    std::uniform_real_distribution<double> U_omega(1.0, 3.0);   // velocidad angular inicial [rad/s]
    std::uniform_real_distribution<double> U_mass(0.8, 1.2);    // masa [kg]

    std::vector<double> m(p.n);
    std::vector<double> y(4 * p.n);

    for (int i = 0; i < p.n; ++i) {
        double mi   = U_mass(rng);
        double rho0 = U_rho(rng);
        double phi0 = U_phi(rng);
        double vrad0 = U_vrad(rng);
        double om0   = U_omega(rng);

        m[i] = mi;
        y[4*i + 0] = rho0;
        y[4*i + 1] = phi0;
        // En t=0, e^{gamma*0}=1, por lo que Prho = m*rho_dot y Pphi = m*rho^2*phi_dot
        y[4*i + 2] = mi * vrad0;
        y[4*i + 3] = mi * rho0 * rho0 * om0;
    }

    int nsteps = (int)std::round(p.t_max / p.dt);

    std::vector<double> tvec, Hvec;
    tvec.reserve(nsteps / p.save_every + 1);
    Hvec.reserve(nsteps / p.save_every + 1);

    std::ofstream fout("swirling_hamiltonian.csv");
    fout << "t,H\n";

    double t = 0.0;
    for (int step = 0; step <= nsteps; ++step) {
        if (step % p.save_every == 0) {
            double H = hamiltonian(t, y, m, p);
            tvec.push_back(t);
            Hvec.push_back(H);
            fout << std::setprecision(10) << t << "," << H << "\n";
        }
        if (step == nsteps) break;
        euler_richardson_step(t, y, m, p);
        t += p.dt;
    }
    fout.close();

double tau = relaxation_time(tvec, Hvec, 0.02);

std::cout << "Simulacion completada.\n";
std::cout << "  Particulas (n)          : " << p.n << "\n";
std::cout << "  gamma (Caldirola-Kanai): " << p.gamma_ << " 1/s\n";
std::cout << "  H(0)                    : " << Hvec.front() << "\n";
std::cout << "  H_inf (asintotico)      : " << Hvec.back() << "\n";
std::cout << "  Tiempo de relajacion    : " << tau << " s\n";
std::cout << "  Datos guardados en      : swirling_hamiltonian.csv\n";

save_graph_svg(tvec, Hvec, tau);
std::cout << "Grafica guardada en      : relajacion_hamiltoniano.svg\n";

    return 0;
}