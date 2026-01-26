#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <vector>

enum SolverType {NEUTRONIC, COUPLED};
enum SchemeType {IMPLICIT_EULER, CRANK_NICOLSON}; // CK pas implementé pour l'instant
enum BCType {DIRICHLET, NEUMANN, ROBIN}; // Pareil pour Robin

class SimulationConfig {
    public:
        SolverType solver_type = SolverType::NEUTRONIC;
        SchemeType scheme_type = SchemeType::IMPLICIT_EULER;
        BCType bc_type = BCType::DIRICHLET;
        int dimension = 2;
        int N = 3;
        double pitch = 1.26;
        double r_fuel = 0.41;
        double r_guide = 0.56;
        std::vector<std::pair<int, int>> guide_positions = {{1,1}};
        int cells_per_pitch = 10;
        double f_ins = 0.0;
        double tf = 100.0;
        double dt = 1.0;
        int n_out = 100;
        double initial_phi = 1e8;
        double initial_T = 300.0;
        std::string length_unit = "cm";
        
        void readFromFile(const std::string& filename);
        void printSummary(const std::string& filename) const;
};

#endif
