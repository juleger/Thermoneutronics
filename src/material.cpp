#include "material.h"
#include "mesh.h"
#include <algorithm>
#include <fstream>
#include <map>
#include <iostream>

// Cette classe représente un matériau avec ses propriétés thermiques et neutroniques, et nucléaire.
// Elle inclut une méthode statique pour interpoler entre deux matériaux (notamment pour les barres de contrôle, selon leur facteur d'insertion).
using namespace std;

Material::Material(const std::string& name, double rho, double Cp, double lambda, double D, double Sigma_a, double Sigma_f, double E_f, double nu, double chi)
    : name(name), rho(rho), Cp(Cp), lambda(lambda), D(D), Sigma_a(Sigma_a), Sigma_f(Sigma_f), E_f(E_f), nu(nu), chi(chi) {}

Material Material::interpolate(const Material& mat1, const Material& mat2, double factor, const string& name) {
    double f = max(0.0, min(1.0, factor));
    double w1 = 1.0 - f;
    double w2 = f;

    return Material(
        name,
        mat1.rho * w1 + mat2.rho * w2,
        mat1.Cp * w1 + mat2.Cp * w2,
        mat1.lambda * w1 + mat2.lambda * w2,
        mat1.D * w1 + mat2.D * w2,
        mat1.Sigma_a * w1 + mat2.Sigma_a * w2,
        mat1.Sigma_f * w1 + mat2.Sigma_f * w2,
        mat1.E_f * w1 + mat2.E_f * w2,
        mat1.nu * w1 + mat2.nu * w2,
        mat1.chi * w1 + mat2.chi * w2
    );
}

void createStandardMaterials(Material& fuel, Material& water, Material& boron) {
    fuel = Material(
        "UO2",     // formule chimique
        10500.0,   // masse volumique
        1040.0,    // capacité calorifique
        5.0,       // conductivité thermique
        1.0,       // coefficient de diffusion neutronique
        0.1,       // section efficace d'absorption
        0.08,      // section efficace de fission
        3.204e-11, // J, énergie libérée par fission
        2.43,      // nombre de neutrons moyen par fission
        1.0        // spectre de fission (chi)
    );

    water = Material(
        "H2O",
        700.0,
        4200.0,
        0.6,
        0.16,  
        0.019,   
        0.0,
        0.0,
        0.0,
        0.0
    );

    boron = Material(
        "B4C",
        2000.0,
        1000.0,
        10.0,
        0.1,
        3.5,
        0.0,
        0.0,
        0.0,
        0.0
    );
}

void writeMeshMaterialsVTK(const Mesh& mesh, const SimulationConfig& config, const std::string& filename) {
    const auto& cells = mesh.getCells();
    size_t N_cells = cells.size();
    int nx = config.N * config.cells_per_pitch;
    int ny = nx;
    double pitch = config.pitch;
    double dx = pitch / config.cells_per_pitch;

    std::ofstream ofs(filename);
    if (!ofs) return;

    ofs << "# vtk DataFile Version 3.0\n";
    ofs << "Mesh materials\n";
    ofs << "ASCII\n";
    ofs << "DATASET POLYDATA\n";

    // Points at grid intersections
    int n_points = (nx + 1) * (ny + 1);
    ofs << "POINTS " << n_points << " float\n";
    for (int i = 0; i <= nx; ++i) {
        for (int j = 0; j <= ny; ++j) {
            double x = i * dx;
            double y = j * dx;
            ofs << x << " " << y << " 0\n";
        }
    }

    ofs << "POLYGONS " << N_cells << " " << 5 * N_cells << "\n";
    for (int ii = 0; ii < nx; ++ii) {
        for (int jj = 0; jj < ny; ++jj) {
            int p0 = ii * (ny + 1) + jj;
            int p1 = (ii + 1) * (ny + 1) + jj;
            int p2 = (ii + 1) * (ny + 1) + (jj + 1);
            int p3 = ii * (ny + 1) + (jj + 1);
            ofs << "4 " << p0 << " " << p1 << " " << p2 << " " << p3 << "\n";
        }
    }

    std::map<std::string,int> idmap;
    int nextId = 0;
    for (const auto& c : cells) {
        if (idmap.find(c.mat.name) == idmap.end()) idmap[c.mat.name] = nextId++;
    }

    ofs << "CELL_DATA " << N_cells << "\n";
    ofs << "SCALARS material_id int 1\n";
    ofs << "LOOKUP_TABLE default\n";
    for (const auto& c : cells) ofs << idmap[c.mat.name] << "\n";

    ofs.close();
    std::cout << "Fichier VTK généré : " << filename << std::endl;
}