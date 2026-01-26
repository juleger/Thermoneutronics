#include "mesh.h"
#include "material.h"
#include <iostream>
#include <cmath>
#include <algorithm>

using namespace std;

// Ce fichier implémente la construction du maillage 2D carré avec les matériaux assignés selon la configuration spécifiée.

Mesh::Mesh(const SimulationConfig& Config) {
    dim = Config.dimension;
    int N = Config.N;
    double pitch = Config.pitch;
    double r_fuel = Config.r_fuel;
    double r_guide = Config.r_guide;
    int cells_per_pitch = Config.cells_per_pitch;
    double f_ins = Config.f_ins; // Facteur d'insertion des barres, influe seulement sur les guides 

    N_cells = N * N * cells_per_pitch * cells_per_pitch;
    this->cells_per_pitch = cells_per_pitch;
    cells.reserve(N_cells);

    Material fuel, water, boron;
    createStandardMaterials(fuel, water, boron);
    Material control_rod = Material::interpolate(water, boron, f_ins, "ControlRod");

    double dx = pitch / cells_per_pitch;
    
    for (int ii = 0; ii < N * cells_per_pitch; ++ii) {
        for (int jj = 0; jj < N * cells_per_pitch; ++jj) {
            double x = (ii + 0.5) * dx;
            double y = (jj + 0.5) * dx;

            int assembly_i = ii / cells_per_pitch;
            int assembly_j = jj / cells_per_pitch;
            double r = sqrt(pow(x - (assembly_i + 0.5) * pitch, 2) + pow(y - (assembly_j + 0.5) * pitch, 2));
            bool is_guide = false;
            for (const auto& pos : Config.guide_positions) {
                if (assembly_i == pos.first && assembly_j == pos.second) {
                    is_guide = true;
                    break;
                }
            }

            Material cell_material = water;
            if (is_guide) {
                if (r <= r_guide) {
                    cell_material = control_rod;
                }
            } else {
                if (r <= r_fuel) {
                    cell_material = fuel;
                }
            }

            bool isBoundary = (ii == 0 || jj == 0 || ii == N * cells_per_pitch - 1 || jj == N * cells_per_pitch - 1);
            int cell_id = ii * N * cells_per_pitch + jj;
            vector<int> neighbors;
            
            if (ii > 0) neighbors.push_back((ii - 1) * N * cells_per_pitch + jj);
            if (ii < N * cells_per_pitch - 1) neighbors.push_back((ii + 1) * N * cells_per_pitch + jj);
            if (jj > 0) neighbors.push_back(ii * N * cells_per_pitch + (jj - 1));
            if (jj < N * cells_per_pitch - 1) neighbors.push_back(ii * N * cells_per_pitch + (jj + 1));
            Cell cell(cell_id, x, y, ii, jj, cell_material, isBoundary, neighbors);
            cells.push_back(cell);
        }
    }
}