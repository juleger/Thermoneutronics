#ifndef MATERIAL_H
#define MATERIAL_H

#include <string>
#include "config.h"

class Mesh;

struct Material {
    std::string name;
    double rho;      // Densité
    double Cp;       // Capacité calorifique massique
    double lambda;   // Conductivité thermique
    double D;        // Coefficient de diffusion neutronique
    double Sigma_a;  // Section efficace d'absorption
    double Sigma_f;  // Section efficace de fission
    double E_f;      // Énergie libérée par fission
    double nu;
    double chi;

    Material()
        : name(""), rho(0), Cp(0), lambda(0), D(0), Sigma_a(0), Sigma_f(0), E_f(0), nu(0), chi(0) {}

    Material(const std::string& name_, double rho_, double Cp_, double lambda_, double D_, double Sigma_a_, double Sigma_f_, double E_f_, double nu_, double chi_);

    static Material interpolate(const Material& mat1, const Material& mat2, double factor, const std::string& name_);
};

void createStandardMaterials(Material& fuel, Material& water, Material& boron);

void writeMeshMaterialsVTK(const Mesh& mesh, const SimulationConfig& config, const std::string& filename);

#endif