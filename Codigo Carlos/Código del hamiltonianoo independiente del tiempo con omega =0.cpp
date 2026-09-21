//  Con Omega=0, el Hamiltoniano conservativo (sin disipar) es:
//
//      H = Sum_i [ prho_i^2/(2 m_i) + pphi_i^2/(2 m_i rho_i^2) + Vborde(rho_i) ]
//
//  FUNCION DE DISIPACION DE RAYLEIGH (Ec. 18, con Omega=0)
//  Con Omega=0 la velocidad de arrastre del plato es cero, asi que la velocidad local
//  de la particula i coincide con su velocidad total. En coordenadas cartesianas del
//  sistema local (no rotante) x',y':
//
//      vx_i = rho_dot_i cos(phi_i) - rho_i phi_dot_i sin(phi_i)
//      vy_i = rho_dot_i sin(phi_i) + rho_i phi_dot_i cos(phi_i)
//
//      F = (gamma/4) Sum_i Sum_j [ (vx_i-vx_j)^2 + (vy_i-vy_j)^2 ]
//
// -------------------------------------------------------------------------------------
//  ECUACIONES DE MOVIMIENTO FINALES (Omega = 0)
// -------------------------------------------------------------------------------------
//      rho_dot_i =  prho_i / m_i
//      phi_dot_i =  pphi_i / (m_i rho_i^2)
//
//      prho_dot_i =  pphi_i^2/(m_i rho_i^3) - k(rho_i-R_PLATO)*[rho_i>R_PLATO]
//                    - dF/d(rho_dot_i)
//
//      pphi_dot_i =  - dF/d(phi_dot_i)
//
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>
#include <algorithm>

// ---------------------------- Parametros fisicos globales ----------------------------
struct Params {
    int    n         = 6;      // numero de particulas (igual que en las versiones anteriores)
    double gamma_    = 0.5;    // coeficiente de friccion de Rayleigh [misma magnitud que
                                // el "gamma" usado en Caldirola-Kanai, para comparar]
    double k         = 60.0;   // rigidez de la pared blanda del plato [N/m]
    double R_plato   = 1.0;    // radio del plato [m]
    double t_max     = 20.0;   // tiempo total de simulacion [s]
    double dt        = 0.0005; // paso de integracion [s]
    int    save_every = 40;    // guardar 1 de cada 'save_every' pasos al CSV
};

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

// Calcula, para el estado actual, las velocidades cartesianas locales v=(vx,vy) de
// cada particula (usadas por la funcion de disipacion de Rayleigh).
void local_velocities(const std::vector<double>& y, const std::vector<double>& m,
                       int n, std::vector<double>& vx, std::vector<double>& vy,
                       std::vector<double>& rho_dot, std::vector<double>& phi_dot) {
    for (int i = 0; i < n; ++i) {
        double rho = y[4*i+0], phi = y[4*i+1], prho = y[4*i+2], pphi = y[4*i+3];
        double mi = m[i];
        rho_dot[i] = prho / mi;
        phi_dot[i] = pphi / (mi * rho * rho);
        vx[i] = rho_dot[i]*std::cos(phi) - rho*phi_dot[i]*std::sin(phi);
        vy[i] = rho_dot[i]*std::sin(phi) + rho*phi_dot[i]*std::cos(phi);
    }
}

// Funcion de disipacion de Rayleigh total F (Ec. 18, Omega=0)
double rayleigh_F(const std::vector<double>& vx, const std::vector<double>& vy,
                   int n, double gamma) {
    double F = 0.0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            double dvx = vx[i]-vx[j], dvy = vy[i]-vy[j];
            F += dvx*dvx + dvy*dvy;
        }
    return 0.25 * gamma * F;
}

// dy/dt para todo el sistema. Vector de estado: [rho_1,phi_1,prho_1,pphi_1, ...]
std::vector<double> derivatives(const std::vector<double>& y, const std::vector<double>& m,
                                 const Params& p) {
    int n = p.n;
    std::vector<double> dydt(4*n);
    std::vector<double> vx(n), vy(n), rho_dot(n), phi_dot(n);
    local_velocities(y, m, n, vx, vy, rho_dot, phi_dot);

    for (int i = 0; i < n; ++i) {
        double rho = y[4*i+0], phi = y[4*i+1], pphi = y[4*i+3], mi = m[i];
        double cphi = std::cos(phi), sphi = std::sin(phi);

        double sum_r = 0.0, sum_phi = 0.0;
        for (int j = 0; j < n; ++j) {
            if (j == i) continue;
            double dvx = vx[i]-vx[j], dvy = vy[i]-vy[j];
            sum_r   += dvx*cphi + dvy*sphi;          // (v_i-v_j).rhat_i
            sum_phi += -dvx*sphi + dvy*cphi;         // (v_i-v_j).phihat_i
        }
        double dF_drhodot = p.gamma_ * sum_r;
        double dF_dphidot = p.gamma_ * rho * sum_phi;

        dydt[4*i+0] = rho_dot[i];
        dydt[4*i+1] = phi_dot[i];
        dydt[4*i+2] = (pphi*pphi)/(mi*rho*rho*rho) - dVborde_drho(rho, p) - dF_drhodot;
        dydt[4*i+3] = -dF_dphidot;
    }
    return dydt;
}

// Paso de Euler-Richardson (punto medio, 2do orden)
void euler_richardson_step(std::vector<double>& y, const std::vector<double>& m,
                            const Params& p) {
    double dt = p.dt;
    auto f1 = derivatives(y, m, p);
    std::vector<double> y_mid(y.size());
    for (size_t j = 0; j < y.size(); ++j) y_mid[j] = y[j] + 0.5*dt*f1[j];
    auto f_mid = derivatives(y_mid, m, p);
    for (size_t j = 0; j < y.size(); ++j) y[j] += dt*f_mid[j];
}

// Hamiltoniano conservativo (T+V), independiente del tiempo (Omega=0)
double hamiltonian(const std::vector<double>& y, const std::vector<double>& m,
                    const Params& p) {
    int n = p.n;
    double H = 0.0;
    for (int i = 0; i < n; ++i) {
        double rho = y[4*i+0], prho = y[4*i+2], pphi = y[4*i+3], mi = m[i];
        H += (prho*prho)/(2*mi) + (pphi*pphi)/(2*mi*rho*rho) + Vborde(rho, p);
    }
    return H;
}

// Momento angular total del sistema (respecto al centro del plato): Lz = Sum m_i rho_i^2 phi_dot_i = Sum pphi_i
double angular_momentum_total(const std::vector<double>& y, int n) {
    double Lz = 0.0;
    for (int i = 0; i < n; ++i) Lz += y[4*i+3];
    return Lz;
}

double relaxation_time(const std::vector<double>& tvec, const std::vector<double>& Hvec,
                        double tol_rel = 0.02) {
    int N = (int)Hvec.size();
    if (N < 10) return -1.0;
    int tail = std::max(5, N/20);
    double H_inf = 0.0;
    for (int i = N-tail; i < N; ++i) H_inf += Hvec[i];
    H_inf /= tail;
    double H0 = Hvec[0];
    double band = tol_rel * std::fabs(H0 - H_inf);
    if (band < 1e-12) band = tol_rel * std::fabs(H0);
    int last_outside = -1;
    for (int i = 0; i < N; ++i) if (std::fabs(Hvec[i]-H_inf) > band) last_outside = i;
    int idx = (last_outside == -1) ? 0 : std::min(last_outside+1, N-1);
    return tvec[idx];
}

void save_graph_svg(const std::vector<double>& t,
                    const std::vector<double>& H,
                    const std::vector<double>& F) {
    const int W = 1100, Ht = 650;
    const int left = 90, right = 90, top = 50, bottom = 80;

    double xmin = t.front();
    double xmax = t.back();

    double Hmin = *std::min_element(H.begin(), H.end());
    double Hmax = *std::max_element(H.begin(), H.end());
    double Fmin = *std::min_element(F.begin(), F.end());
    double Fmax = *std::max_element(F.begin(), F.end());

    if (Hmax == Hmin) Hmax = Hmin + 1.0;
    if (Fmax == Fmin) Fmax = Fmin + 1.0;

    double margenH = 0.05 * (Hmax - Hmin);
    double margenF = 0.05 * (Fmax - Fmin);

    Hmin -= margenH;
    Hmax += margenH;
    Fmin = std::max(0.0, Fmin - margenF);
    Fmax += margenF;

    auto xmap = [&](double x) {
        return left + (x - xmin) / (xmax - xmin) *
               (W - left - right);
    };

    auto yH = [&](double y) {
        return Ht - bottom -
               (y - Hmin) / (Hmax - Hmin) *
               (Ht - top - bottom);
    };

    auto yF = [&](double y) {
        return Ht - bottom -
               (y - Fmin) / (Fmax - Fmin) *
               (Ht - top - bottom);
    };

    std::ofstream svg("hamiltoniano_y_disipacion_vs_tiempo.svg");

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

    // Curva del Hamiltoniano
    svg << "<polyline fill=\"none\" stroke=\"blue\" stroke-width=\"2\" points=\"";
    for (size_t i = 0; i < t.size(); ++i)
        svg << xmap(t[i]) << "," << yH(H[i]) << " ";
    svg << "\"/>\n";

    // Curva disipativa F(t)
    svg << "<polyline fill=\"none\" stroke=\"red\" stroke-width=\"2\" points=\"";
    for (size_t i = 0; i < t.size(); ++i)
        svg << xmap(t[i]) << "," << yF(F[i]) << " ";
    svg << "\"/>\n";

    // Título
    svg << "<text x=\"" << W / 2 << "\" y=\"25\" "
        << "text-anchor=\"middle\" font-size=\"18\">"
        << "Hamiltoniano y funcion disipativa</text>\n";

    // Etiquetas
    svg << "<text x=\"" << W / 2 << "\" y=\"" << Ht - 20
        << "\" text-anchor=\"middle\">Tiempo t [s]</text>\n";

    svg << "<text x=\"20\" y=\"" << Ht / 2
        << "\" transform=\"rotate(-90 20," << Ht / 2
        << ")\" text-anchor=\"middle\" fill=\"blue\">H(t)</text>\n";

    svg << "<text x=\"" << W - 20 << "\" y=\"" << Ht / 2
        << "\" transform=\"rotate(90 " << W - 20 << "," << Ht / 2
        << ")\" text-anchor=\"middle\" fill=\"red\">F(t)</text>\n";

    // Leyenda
    svg << "<line x1=\"120\" y1=\"55\" x2=\"150\" y2=\"55\" "
        << "stroke=\"blue\" stroke-width=\"3\"/>\n";
    svg << "<text x=\"160\" y=\"60\" fill=\"blue\">Hamiltoniano H(t)</text>\n";

    svg << "<line x1=\"350\" y1=\"55\" x2=\"380\" y2=\"55\" "
        << "stroke=\"red\" stroke-width=\"3\"/>\n";
    svg << "<text x=\"390\" y=\"60\" fill=\"red\">Disipacion F(t)</text>\n";

    svg << "</svg>\n";
}

int main() {
    Params p;

    // Mismas distribuciones, mismo orden y misma semilla que en swirling_hamiltonian.cpp
    // (version Caldirola-Kanai, Omega=0) para que H(0) coincida exactamente y la
    // comparacion entre ambos metodos sea directa.
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> U_rho(0.25*p.R_plato, 0.85*p.R_plato);
    std::uniform_real_distribution<double> U_phi(
        0.0, 2.0 * std::acos(-1.0)
    );
    std::uniform_real_distribution<double> U_vrad(-0.8, 0.8);
    std::uniform_real_distribution<double> U_omega(1.0, 3.0);
    std::uniform_real_distribution<double> U_mass(0.8, 1.2);

    std::vector<double> m(p.n), y(4*p.n);
    for (int i = 0; i < p.n; ++i) {
        double mi = U_mass(rng), rho0 = U_rho(rng), phi0 = U_phi(rng);
        double vrad0 = U_vrad(rng), om0 = U_omega(rng);
        m[i] = mi;
        y[4*i+0] = rho0;
        y[4*i+1] = phi0;
        y[4*i+2] = mi*vrad0;              // prho_i(0) = m rho_dot(0)  (sin factor CK)
        y[4*i+3] = mi*rho0*rho0*om0;      // pphi_i(0) = m rho^2 phi_dot(0)
    }

    int nsteps = (int)std::round(p.t_max/p.dt);
    std::vector<double> tvec, Hvec, Fvec, Lzvec;
    tvec.reserve(nsteps/p.save_every+1);
    Hvec.reserve(nsteps/p.save_every+1);
    Fvec.reserve(nsteps/p.save_every+1);
    Lzvec.reserve(nsteps/p.save_every+1);

    std::ofstream fout("swirling_rayleigh.csv");
    fout << "t,H,F,Lz_total\n";

    double t = 0.0;
    for (int step = 0; step <= nsteps; ++step) {
        if (step % p.save_every == 0) {
            std::vector<double> vx(p.n), vy(p.n), rho_dot(p.n), phi_dot(p.n);
            local_velocities(y, m, p.n, vx, vy, rho_dot, phi_dot);
            double H = hamiltonian(y, m, p);
            double F = rayleigh_F(vx, vy, p.n, p.gamma_);
            double Lz = angular_momentum_total(y, p.n);
            tvec.push_back(t); Hvec.push_back(H); Fvec.push_back(F); Lzvec.push_back(Lz);
            fout << std::setprecision(10) << t << "," << H << "," << F << "," << Lz << "\n";
        }
        if (step == nsteps) break;
        euler_richardson_step(y, m, p);
        t += p.dt;
    }
    fout.close();

    double tau = relaxation_time(tvec, Hvec, 0.02);

    // Verificacion de la identidad dH/dt = -2F, usando diferencias finitas sobre la
    // serie guardada de H(t) (chequeo grueso, con el paso de guardado, no con dt).
    double max_check_err = 0.0;
    for (size_t i = 1; i+1 < tvec.size(); ++i) {
        double dHdt_num = (Hvec[i+1]-Hvec[i-1]) / (tvec[i+1]-tvec[i-1]);
        double pred = -2.0*Fvec[i];
        max_check_err = std::max(max_check_err, std::fabs(dHdt_num - pred));
    }

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Simulacion completada (Omega=0, disipacion de Rayleigh explicita).\n";
    std::cout << "  Particulas (n)         : " << p.n << "\n";
    std::cout << "  gamma (Rayleigh)       : " << p.gamma_ << "\n";
    std::cout << "  H(0)                   : " << Hvec.front() << "\n";
    std::cout << "  H_inf (asintotico)     : " << Hvec.back() << "\n";
    std::cout << "  Tiempo de relajacion   : " << tau << " s (banda del 2%)\n";
    std::cout << "  Lz_total(0)            : " << Lzvec.front() << "\n";
    std::cout << "  Lz_total(t_max)        : " << Lzvec.back()
               << "   (NO tiene por que conservarse: la friccion es antisimetrica mas\n"
               << "                            no es una fuerza central, asi que puede\n"
               << "                            ejercer torque neto; ver nota en el codigo)\n";
    std::cout << "  max |dH/dt - (-2F)|    : " << max_check_err
               << "  (verificacion de la identidad de disipacion de Rayleigh)\n";
    std::cout << "Datos guardados en swirling_rayleigh.csv\n";

    save_graph_svg(tvec, Hvec, Fvec);

    std::cout << "Grafica guardada en "
              << "hamiltoniano_y_disipacion_vs_tiempo.svg\n";

    return 0;
}