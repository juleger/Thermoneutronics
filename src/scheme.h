#ifndef SCHEME_H
#define SCHEME_H

#include "mesh.h"
#include "material.h"
#include "config.h"
#include "solver.h"
#include <Eigen/Dense>
#include <vector>

class Scheme {
    protected:
        Mesh* mesh;
        const SimulationConfig& config;
        double t;
        double dt;
        Eigen::VectorXd phi;
        Eigen::VectorXd T;

        NeutronicSolver* neutronic_solver;
        ThermalSolver* thermal_solver;
        
    public:

        Scheme(Mesh* mesh, const SimulationConfig& config);
        virtual ~Scheme() {
            delete neutronic_solver;
            delete thermal_solver;
        }

        void Initialize();

        virtual void Advance() = 0;

        void exportResults(const std::string& filename);
        
        
        const Eigen::VectorXd& getPhi() const { return phi; }
        const Eigen::VectorXd& getT() const { return T; }
        void setPhi(const Eigen::VectorXd& newPhi);

        // Calcule un proxy k_inf instantané: sum(nu*Sigma_f*phi)/sum(Sigma_a*phi)
        double computeKinf() const;

        // Calcule k_inf par résolution de l'équation aux valeurs propres stationnaire (Neumann):
        // (-D ∇² + Σ_a) φ = (1/k) ν Σ_f φ  → k via itération de puissance.
        double computeKinfEigen(int maxIters = 100, double tol = 1e-8) const;
};

class ImplicitEulerScheme : public Scheme {
    public:
    ImplicitEulerScheme(Mesh* mesh, const SimulationConfig& config) : Scheme(mesh, config) {}
        void Advance() override;
};

#endif