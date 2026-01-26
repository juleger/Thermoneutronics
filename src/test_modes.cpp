#include "test_modes.h"
#include "config.h"
#include "mesh.h"
#include "scheme.h"
#include "material.h"
#include <iostream>
#include <vector>
#include <cmath>

// Le but de ce test est de vérifier que les modes propres analytiques sont bien reproduits par le solveur neutronique.
// On considère un domaine carré unique avec conditions de Neumann sur toutes les faces, et un matériau uniforme (UO2).
// On initialise le champ de flux neutronique avec une combinaison de modes propres analytiques, et on laisse évoluer la simulation.
// Au bout de tf, on compare le résultat numérique avec la solution analytique obtenue par la décroissance exponentielle des modes propres.
// La thermique n'est pas prise en compte dans ce test, on se concentre uniquement sur le solveur neutronique qui est le plus pertinent (c'est lui qui influence la thermique et pas l'inverse)

void runEigenmodeTest() {
    std::cout << "=== Test Eigenmodes (Neumann, UO2 uniforme) ===" << std::endl;
    // Config de test: carré unique, UO2 partout, Neumann
    SimulationConfig cfg;
    cfg.N = 1;
    cfg.cells_per_pitch = 200;
    cfg.pitch = 1.0; // longueur L
    cfg.bc_type = BCType::NEUMANN;
    cfg.f_ins = 0.0;
    cfg.guide_positions.clear();
    cfg.r_guide = 0.0;
    cfg.r_fuel = 10.0 * cfg.pitch; // Afin d'avoir seulement du combustible dans le domaine
    cfg.tf = 1;
    cfg.dt = 5e-4;
    cfg.initial_phi = 1.0;

    Mesh mesh(cfg);
    ImplicitEulerScheme scheme(&mesh, cfg);
    scheme.Initialize();

    const auto& cells = mesh.getCells();
    int Ncells = mesh.getNumCells();
    double L = cfg.N * cfg.pitch;

    const Material& mat = cells[0].mat;
    double D = mat.D;
    double Sigma_a = mat.Sigma_a;
    double Sigma_f = mat.Sigma_f;
    double nu = mat.nu;

    std::vector<std::pair<int,int>> modes;
    modes.push_back({0,1});
    modes.push_back({1,0});
    modes.push_back({1,1});

    for (size_t idx_mode = 0; idx_mode < modes.size(); ++idx_mode) {
        int n = modes[idx_mode].first;
        int m = modes[idx_mode].second;
        // On initialise le mode propre analytique
        // phi(,x,y) = cos(n*pi*x/L) * cos(m*pi*y/L) 
        Eigen::VectorXd phi0(Ncells);
        for (int c = 0; c < Ncells; ++c) {
            double x = cells[c].x; // [0, L]
            double y = cells[c].y; // [0, L]
            double valx = (n == 0) ? 1.0 : std::cos(M_PI * n * x / L);
            double valy = (m == 0) ? 1.0 : std::cos(M_PI * m * y / L);
            phi0(c) = valx * valy;
        }
        scheme.setPhi(phi0);

        int steps = static_cast<int>(cfg.tf / cfg.dt);
        for (int n = 0; n < steps; ++n) {
            scheme.Advance();
        }

        double k2 = (n == 0 ? 0.0 : (M_PI * n / L) * (M_PI * n / L)) +
                    (m == 0 ? 0.0 : (M_PI * m / L) * (M_PI * m / L));
        double lambda = -D * k2 - Sigma_a + nu * Sigma_f;
        double amp = std::exp(lambda * cfg.tf);
        Eigen::VectorXd phi_exact = amp * phi0;
        const Eigen::VectorXd& phi_num = scheme.getPhi();

        if (phi_exact.norm() == 0.0) {
            continue;
        }

        double errL2 = (phi_num - phi_exact).norm() / phi_exact.norm();

        std::cout << "Mode (n=" << n << ", m=" << m << ") : Erreur L2 relative = " << errL2 << std::endl;
    }

    std::cout << "=== Fin du test eigenmodes ===" << std::endl;
}
