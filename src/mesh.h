#ifndef MESH_H
#define MESH_H

#include <vector>
#include "material.h"
#include "config.h"
 
struct Cell {
    // Maille 2D en cartésien
    int id;           // ID unique
    double x, y;      // Coordonnées centrales
    int i, j;         // Indices dans la grille

    Material mat;      // Matériau de la cellule
    bool isBoundary = false; 
    std::vector<int> neighbors; // IDs des cellules voisines
    Cell(int id_, double x_, double y_, int i_, int j_, const Material& mat_, bool isBoundary_, const std::vector<int>& neighbors_)
        : id(id_), x(x_), y(y_), i(i_), j(j_), mat(mat_), isBoundary(isBoundary_), neighbors(neighbors_) {}

}; 

class Mesh {
    private:
        std::vector<Cell> cells;
        int N_cells;
        int dim;
        int cells_per_pitch;
    public:
        // Constructeur
        Mesh(const SimulationConfig& config);
        ~Mesh() {}

        // Getters
        const std::vector<Cell>& getCells() const { return cells; }
        int getNumCells() const { return N_cells; }
        int getCellsPerPitch() const { return cells_per_pitch; }

        //Setters
        void addCell(const Cell& cell) { cells.push_back(cell); }

};
#endif