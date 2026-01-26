#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "mesh.h"
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

class BoundaryCondition {
        //On identifie sur quel(s) bord(s) la condition est appliquée
        std::vector<int> boundary_ids; 
    
    public:
        BoundaryCondition() {}
        virtual ~BoundaryCondition() {}
        virtual void Apply(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b, const Mesh& mesh) = 0;
};

class DirichletBoundaryCondition : public BoundaryCondition {
    private:
        double value;
        
    public:
        DirichletBoundaryCondition(double value_) : value(value_) {}
        void Apply(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b, const Mesh& mesh) override;
};

class NeumannBoundaryCondition : public BoundaryCondition {
    private:
        double flux;
    public:
        NeumannBoundaryCondition(double flux_) : flux(flux_) {}
        void Apply(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b, const Mesh& mesh) override;
};
#endif