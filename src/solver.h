#ifndef SOLVER_H
#define SOLVER_H

#include "mesh.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

class Solver {
protected:
    Mesh* mesh;
    double dt;
    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd& u;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> lu_solver;

public:
    Solver(Mesh* mesh, double dt, Eigen::VectorXd& u) : mesh(mesh), dt(dt), u(u) {}
    virtual ~Solver() {}

    virtual void buildMatrix() = 0;
    virtual void solveStep() = 0;
};

class NeutronicSolver : public Solver {
private:
    double dx;
    Eigen::SparseMatrix<double> LHS;

public:
    NeutronicSolver(Mesh* mesh, double dt, Eigen::VectorXd& phi);
    void buildMatrix() override;
    void solveStep() override;
};

class ThermalSolver : public Solver {
private:
    double dx;
    double length_scale;
    Eigen::VectorXd Q_thermal;
    const Eigen::VectorXd& phi;

public:
    ThermalSolver(Mesh* mesh, double dt, Eigen::VectorXd& T, const Eigen::VectorXd& phi, double length_scale);
    void buildMatrix() override;
    void solveStep() override;
};

#endif