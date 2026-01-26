#include "material.h"
#include "mesh.h"
#include "scheme.h"
#include "config.h"
#include <iostream>
#include <fstream>
#include <memory>
#include <chrono>
#include "test_modes.h"
#include <limits>

using namespace std;

// Ce projet simule la résolution d'équations de diffusion neutronique et thermique dans un assemblage nucléaire 2D carré.
// Le maillage est construit en fonction de la configuration spécifiée, avec des matériaux assignés selon la géométrie de l'assemblage.
// Le projet se veut être physiquement réaliste (au niveau des grandeurs, des matériaux, etc...) tout en restant simple et accessible pour notre niveau.
// Les équations omettent un certain nombre d'effets très complexes dans les réacteurs réels :
// neutrons multi-groupes, effet Doppler, convection thermique, géométrie 3D, dilatation thermique, etc...

int main(int argc, char* argv[]) {
    cout << "=== Simulation Thermique et Neutronique d'un assemblage nucléaire ===" << endl;
    string config_file = "../config/config.txt";
    if (argc > 1) {
        config_file = argv[1];
    }
    if (config_file == string("test_eigenmode")) {
        runEigenmodeTest();
        return 0;
    }
    SimulationConfig config;
    config.readFromFile(config_file);
    config.printSummary(config_file);
    
    Mesh mesh(config);
    writeMeshMaterialsVTK(mesh, config, "../results/mesh_materials.vtk"); // Export du maillage avec matériaux pour visualisation
    
    // On initialise le schéma numérique
    ImplicitEulerScheme scheme(&mesh, config);
    scheme.Initialize();
    cout << "Début de la simulation..." << endl;
    // k_inf théorique (Neumann) via eigenvalue (généralisé):
    cout << "k_inf théorique (Neumann, eigen): " << scheme.computeKinfEigen() << endl;
    // k_inf théorique (Neumann, uniforme): sum(nu*Sigma_f)/sum(Sigma_a)
    {
        const auto& cells = mesh.getCells();
        double sum_prod_uniform = 0.0, sum_abs_uniform = 0.0;
        for (const auto& c : cells) {
            sum_prod_uniform += c.mat.nu * c.mat.Sigma_f;
            sum_abs_uniform  += c.mat.Sigma_a;
        }
        double kinf_uniform = (sum_abs_uniform > 0.0) ? (sum_prod_uniform / sum_abs_uniform) : std::numeric_limits<double>::infinity();
        cout << "k_inf théorique (Neumann, uniforme): " << kinf_uniform << endl;
    }
    // k_inf théorique (Neumann) – deux variantes:
    // 1) Uniforme (flux supposé uniforme): sum(nu*Sigma_f)/sum(Sigma_a)
    // 2) Pondéré par le flux initial:      sum(nu*Sigma_f*phi0)/sum(Sigma_a*phi0)
    {
        const auto& cells = mesh.getCells();
        // Uniforme
        double sum_prod_uniform = 0.0, sum_abs_uniform = 0.0;
        for (const auto& c : cells) {
            sum_prod_uniform += c.mat.nu * c.mat.Sigma_f;
            sum_abs_uniform  += c.mat.Sigma_a;
        }
        double kinf_uniform = (sum_abs_uniform > 0.0) ? (sum_prod_uniform / sum_abs_uniform) : std::numeric_limits<double>::infinity();

        // Pondéré par le flux initial
        const auto& phi0 = scheme.getPhi();
        double sum_prod_phi = 0.0, sum_abs_phi = 0.0;
        for (size_t i = 0; i < cells.size(); ++i) {
            const auto& c = cells[i];
            sum_prod_phi += c.mat.nu * c.mat.Sigma_f * phi0(i);
            sum_abs_phi  += c.mat.Sigma_a * phi0(i);
        }
        double kinf_weighted_init = (sum_abs_phi > 0.0) ? (sum_prod_phi / sum_abs_phi) : std::numeric_limits<double>::infinity();
        cout << "k_inf théorique (Neumann, uniforme): " << kinf_uniform << endl;
        cout << "k_inf théorique (Neumann, pondéré par flux initial): " << kinf_weighted_init << endl;
    }
    
    int n_steps = static_cast<int>(config.tf / config.dt);
    int output_interval = (config.n_out > 0) ? std::max(1, n_steps / config.n_out) : std::max(1, n_steps);

    auto start = std::chrono::high_resolution_clock::now();

    // Prepare flux time series export (mean flux and k_inf proxy)
    std::string seriesName = "../results/flux_N" + std::to_string(config.N) + "_c" + std::to_string(config.N * config.cells_per_pitch) + "_fins" + std::to_string((int)config.f_ins) + ".dat";
    std::ofstream fluxSeries(seriesName);
    fluxSeries << "# t\tmean_phi\tmean_T\tk_inf_proxy\n";

    double kinf_sum = 0.0; int kinf_count = 0;
    for(int n = 0; n < n_steps; ++n) {
        scheme.Advance();
        // Record mean flux and k_inf at each step for quick validation
        double tcur = (n+1) * config.dt;
        double kinf_cur = scheme.computeKinf();
        double mean_T = scheme.getT().mean();
        fluxSeries << tcur << "\t" << scheme.getPhi().mean() << "\t" << mean_T << "\t" << kinf_cur << "\n";
        kinf_sum += kinf_cur; kinf_count++;
        if(n % output_interval == 0) {
            std::string filename = "../results/solution_N" + std::to_string(config.N) + "_c" + std::to_string(config.N * config.cells_per_pitch) + "_fins" + std::to_string((int)config.f_ins) + "_t" + std::to_string(n) + ".vtk";
            scheme.exportResults(filename);   
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    cout << "Simulation completée en " << duration.count() << " ms" << endl;
    fluxSeries.close();
    if (kinf_count > 0) {
        cout << "k_inf moyen numérique (sur pas de temps): " << (kinf_sum/kinf_count) << endl;
    }
    
    // Métriques simples pour avoir une idée immédiate des résultats
    cout << "Moyenne du flux neutronique phi dans le domaine : " << scheme.getPhi().mean() << endl;
    cout << "Moyenne de la température T dans le domaine : " << scheme.getT().mean() << endl;
    cout << "Fichier série flux/k_inf : " << seriesName << endl;
    // k_inf théorique (Neumann, eigen) à la fin de la simulation
    cout << "k_inf théorique (Neumann, eigen, fin): " << scheme.computeKinfEigen() << endl;
    // Diagnostic final: contributions pondérées par flux par matériau et k_inf final
    {
        const auto& cells = mesh.getCells();
        const auto& phi = scheme.getPhi();
        double prod_all = 0.0, abs_all = 0.0;
        double prod_fuel = 0.0, abs_fuel = 0.0, phi_fuel = 0.0;
        double abs_water = 0.0, phi_water = 0.0;
        double abs_rod = 0.0, phi_rod = 0.0;
        for (size_t i = 0; i < cells.size(); ++i) {
            const auto& c = cells[i];
            double ph = phi(i);
            double p = c.mat.nu * c.mat.Sigma_f * ph;
            double a = c.mat.Sigma_a * ph;
            prod_all += p; abs_all += a;
            if (c.mat.name == std::string("UO2")) {
                prod_fuel += p; abs_fuel += a; phi_fuel += ph;
            } else if (c.mat.name == std::string("H2O")) {
                abs_water += a; phi_water += ph;
            } else if (c.mat.name == std::string("ControlRod") || c.mat.name == std::string("B4C")) {
                abs_rod += a; phi_rod += ph;
            }
        }
        double kinf_weighted_final = (abs_all > 0.0) ? (prod_all / abs_all) : std::numeric_limits<double>::infinity();
        cout << "k_inf numérique (état final, pondéré flux): " << kinf_weighted_final << endl;
        cout << "Flux par matériau (∑phi): fuel=" << phi_fuel << ", eau=" << phi_water << ", barre=" << phi_rod << endl;
        cout << "Absorption pondérée (∑Σ_a·phi): fuel=" << abs_fuel << ", eau=" << abs_water << ", barre=" << abs_rod << endl;
    }

    

    
    return 0;
}
 