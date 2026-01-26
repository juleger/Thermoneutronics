#include "solver.h"
#include <vector>
#include <iostream>

// L'équation de la diffusion neutronique s'écrit sous la forme :
//      ∂φ/∂t = D ∇²φ - Σ_a φ + ν Σ_f φ
// L'équation de la thermique qui en découle est :
//      ρ Cp ∂T/∂t = ∇·(λ ∇T) + Q_thermal avec Q_thermal = E_f ν Σ_f φ
// Le schéma d'Euler implicite sous forme matricielle est utilisée, ce fichier permet donc de construire les matrices associées et de les résoudre.

NeutronicSolver::NeutronicSolver(Mesh* mesh, double dt, Eigen::VectorXd& phi)
    : Solver(mesh, dt, phi) {
    dx = mesh->getCells()[1].y - mesh->getCells()[0].y;
    A.resize(mesh->getNumCells(), mesh->getNumCells());
    buildMatrix();
    // Precompute and factorize the constant system matrix LHS = I - dt*A
    LHS.resize(A.rows(), A.cols());
    LHS.setIdentity();
    LHS -= dt * A;
    lu_solver.analyzePattern(LHS);
    lu_solver.factorize(LHS);
}

void NeutronicSolver::buildMatrix() {
    // Construction de la matrice A pour le schéma implicite
    const auto& cells = mesh->getCells();
    int N = cells.size();
    std::vector<Eigen::Triplet<double>> triplets;
    for (int i = 0; i < N; ++i) {
        const Cell& c = cells[i];
        double D = c.mat.D;
        double Sigma_a = c.mat.Sigma_a;
        double Sigma_f = c.mat.nu * c.mat.Sigma_f;
        double sum_lap = 0.0;
        // Contribution des voisins pour le terme laplacien
        for (int j : c.neighbors) {
            double D_j = cells[j].mat.D;
            double lap_coeff_ij = (2 * D * D_j / (D + D_j)) / (dx * dx);
            triplets.push_back(Eigen::Triplet<double>(i, j, lap_coeff_ij));
            sum_lap += lap_coeff_ij;
        }
        double diag = -sum_lap - Sigma_a + Sigma_f;
        triplets.push_back(Eigen::Triplet<double>(i, i, diag));
    }
    // Assemblage final de la matrice A
    A.setFromTriplets(triplets.begin(), triplets.end());
}

void NeutronicSolver::solveStep() {
    // Solve (I - dt A) φ^{n+1} = φ^{n} using pre-factorized LHS
    Eigen::VectorXd rhs = u;
    u = lu_solver.solve(rhs);
}

ThermalSolver::ThermalSolver(Mesh* mesh, double dt, Eigen::VectorXd& T, const Eigen::VectorXd& phi, double length_scale)
    : Solver(mesh, dt, T), length_scale(length_scale), phi(phi) {
    dx = mesh->getCells()[1].y - mesh->getCells()[0].y;
    A.resize(mesh->getNumCells(), mesh->getNumCells());
    Q_thermal.resize(mesh->getNumCells());
    buildMatrix();
    lu_solver.compute(A);
}

void ThermalSolver::buildMatrix() {
    // Construction de la matrice A pour le schéma implicite
    const auto& cells = mesh->getCells();
    int N = cells.size();
    std::vector<Eigen::Triplet<double>> triplets;
    for (int i = 0; i < N; ++i) {
        const Cell& c = cells[i];
        double lambda = c.mat.lambda;
        double rho = c.mat.rho;
        double Cp = c.mat.Cp;
        double dx_m = dx * length_scale;
        double V = dx_m * dx_m;
        double sum_lap = 0.0;
        for (int j : c.neighbors) {
            double lambda_j = cells[j].mat.lambda;
            double lap_coeff_ij = dt * (2 * lambda * lambda_j / (lambda + lambda_j)) / (dx_m * dx_m);
            triplets.push_back(Eigen::Triplet<double>(i, j, -lap_coeff_ij));
            sum_lap += lap_coeff_ij;
        }
        double diag = rho * Cp * V + sum_lap;
        triplets.push_back(Eigen::Triplet<double>(i, i, diag));
    }
    A.setFromTriplets(triplets.begin(), triplets.end());
}

void ThermalSolver::solveStep() {
    // Calcul de la source thermique Q_thermal à partir du flux neutronique phi
    const auto& cells = mesh->getCells();
    int N = cells.size();
    for (int i = 0; i < N; ++i) {
        const Cell& c = cells[i];
        double Sigma_f_SI = c.mat.Sigma_f / length_scale;
        double phi_SI = phi(i) / (length_scale * length_scale);
        Q_thermal(i) = c.mat.E_f * c.mat.nu * Sigma_f_SI * phi_SI;
    }
    Eigen::VectorXd b(N);
    // Construction du second membre
    for (int i = 0; i < N; ++i) {
        const Cell& c = cells[i];
        double dx_m = dx * length_scale;
        double V = dx_m * dx_m;
        b(i) = c.mat.rho * c.mat.Cp * V * u(i) + dt * Q_thermal(i) * V;
    }
    // Résolution du système linéaire A T^{n+1} = b
    u = lu_solver.solve(b);
}