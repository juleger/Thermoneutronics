#include "scheme.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <limits>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

using namespace std;
using namespace Eigen;

// Ce fichier implémente les schémas numériques pour la résolution couplée des équations neutroniques et thermiques,
// il est complémentaire avec solver.cpp qui implémente les solveurs individuels. Le fichier permet aussi d'exporter les résultats.

Scheme::Scheme(Mesh* mesh, const SimulationConfig& config)
    : mesh(mesh), config(config), t(0.0), dt(config.dt), neutronic_solver(nullptr), thermal_solver(nullptr) {
    int N = mesh->getNumCells();
    phi.resize(N);
    T.resize(N);
}

void Scheme::Initialize() {
    int N = mesh->getNumCells();
    phi.setConstant(config.initial_phi);
    T.setConstant(config.initial_T);
    const auto& cells = mesh->getCells();
    for (int i = 0; i < N; ++i) {
        if (cells[i].isBoundary && config.bc_type == BCType::DIRICHLET) phi(i) = 0.0;
    }
    double length_scale = (config.length_unit == "cm") ? 0.01 : 1.0;
    neutronic_solver = new NeutronicSolver(mesh, dt, phi);
    thermal_solver = new ThermalSolver(mesh, dt, T, phi, length_scale);
}

void Scheme::setPhi(const Eigen::VectorXd& newPhi) {
    phi = newPhi;
}

void Scheme::exportResults(const std::string& filename) {
    std::ofstream file(filename);
    file << "# vtk DataFile Version 3.0\n";
    file << "Neutronic and Thermal Results\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";
    const auto& cells = mesh->getCells();
    int N = cells.size();

    // Construire un maillage surfacique en QUADs sur une grille régulière
    int nx = config.N * config.cells_per_pitch;
    int ny = nx;
    double pitch = config.pitch;
    double dx = pitch / config.cells_per_pitch;

    // Points aux intersections de la grille
    int n_points = (nx + 1) * (ny + 1);
    file << "POINTS " << n_points << " float\n";
    for (int i = 0; i <= nx; ++i) {
        for (int j = 0; j <= ny; ++j) {
            double x = i * dx;
            double y = j * dx;
            file << x << " " << y << " 0.0\n";
        }
    }

    // Cellules QUAD
    file << "CELLS " << N << " " << 5 * N << "\n";
    for (int ii = 0; ii < nx; ++ii) {
        for (int jj = 0; jj < ny; ++jj) {
            int p0 = ii * (ny + 1) + jj;
            int p1 = (ii + 1) * (ny + 1) + jj;
            int p2 = (ii + 1) * (ny + 1) + (jj + 1);
            int p3 = ii * (ny + 1) + (jj + 1);
            file << "4 " << p0 << " " << p1 << " " << p2 << " " << p3 << "\n";
        }
    }

    // Types de cellules : 9 pour QUAD
    file << "CELL_TYPES " << N << "\n";
    for (int i = 0; i < N; ++i) {
        file << "9\n";
    }

    // Données par cellule
    file << "CELL_DATA " << N << "\n";
    file << "SCALARS Neutron_Flux float 1\n";
    file << "LOOKUP_TABLE default\n";   
    for (int i = 0; i < N; ++i) {
        file << phi(i) << "\n";
    }
    file << "SCALARS Temperature float 1\n";
    file << "LOOKUP_TABLE default\n";   
    for (int i = 0; i < N; ++i) {
        file << T(i) << "\n";       
    }
    file.close();
}

double Scheme::computeKinf() const {
    const auto& cells = mesh->getCells();
    const int N = cells.size();
    double production = 0.0;
    double absorption = 0.0;
    for (int i = 0; i < N; ++i) {
        const Cell& c = cells[i];
        const double phi_i = phi(i);
        production += c.mat.nu * c.mat.Sigma_f * phi_i;
        absorption += c.mat.Sigma_a * phi_i;
    }
    if (absorption <= 0.0) return std::numeric_limits<double>::infinity();
    return production / absorption;
}

double Scheme::computeKinfEigen(int maxIters, double tol) const {
    // Build M = (-D ∇² + Σ_a) and F = (ν Σ_f) diag, then power iteration on M^{-1} F
    const auto& cells = mesh->getCells();
    const int N = static_cast<int>(cells.size());
    const double dx = config.pitch / config.cells_per_pitch; // grid spacing

    std::vector<Eigen::Triplet<double>> tripM;
    std::vector<Eigen::Triplet<double>> tripF;
    tripM.reserve(N * 5);
    tripF.reserve(N);

    for (int i = 0; i < N; ++i) {
        const Cell& ci = cells[i];
        double Di = ci.mat.D;
        double Sigma_a = ci.mat.Sigma_a;
        double sum_lap = 0.0;
        for (int j : ci.neighbors) {
            const Cell& cj = cells[j];
            double Dj = cj.mat.D;
            double lap_ij = (2.0 * Di * Dj / (Di + Dj)) / (dx * dx);
            // For M = -D ∇² + Σ_a, off-diagonal is -lap_ij
            tripM.emplace_back(i, j, -lap_ij);
            sum_lap += lap_ij;
        }
        // Diagonal of M: +sum_lap + Σ_a
        tripM.emplace_back(i, i, sum_lap + Sigma_a);
        // Diagonal of F: ν Σ_f
        tripF.emplace_back(i, i, ci.mat.nu * ci.mat.Sigma_f);
    }

    Eigen::SparseMatrix<double> M(N, N);
    Eigen::SparseMatrix<double> F(N, N);
    M.setFromTriplets(tripM.begin(), tripM.end());
    F.setFromTriplets(tripF.begin(), tripF.end());

    Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
    lu.analyzePattern(M);
    lu.factorize(M);

    Eigen::VectorXd v = Eigen::VectorXd::Ones(N);
    double k_prev = 0.0;
    for (int it = 0; it < maxIters; ++it) {
        // y = F v
        Eigen::VectorXd y = F * v;
        // solve M x = y
        Eigen::VectorXd x = lu.solve(y);
        // normalize
        double nrm = x.norm();
        if (nrm <= 0.0) break;
        v = x / nrm;
        // Rayleigh quotient for generalized eigenvalue: k = (v^T F v) / (v^T M v)
        double num = v.dot(F * v);
        double den = v.dot(M * v);
        double k = (den > 0.0) ? (num / den) : std::numeric_limits<double>::infinity();
        if (std::abs(k - k_prev) / (std::abs(k_prev) + 1e-12) < tol) {
            return k;
        }
        k_prev = k;
    }
    return k_prev;
}

void ImplicitEulerScheme::Advance() {
    neutronic_solver->solveStep();
    thermal_solver->solveStep();
    t += dt;
}