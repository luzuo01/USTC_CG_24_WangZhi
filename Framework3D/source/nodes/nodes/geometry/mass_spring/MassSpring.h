#pragma once 
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <set>
#include "utils.h"
#include <chrono>
#include <cassert>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <vector>
#include <iostream>

#define TIC(name) auto start_##name = std::chrono::high_resolution_clock::now(); 
#define TOC(name) \
    auto end_##name = std::chrono::high_resolution_clock::now(); \
    if(enable_time_profiling) std::cout << "Time taken by " << #name << ": " << std::chrono::duration_cast<std::chrono::microseconds>(end_##name- start_##name).count() << " microseconds\n"; \


namespace USTC_CG::node_mass_spring {

typedef OpenMesh::PolyMesh_ArrayKernelT<> SpringMesh;
using namespace Eigen;
using Edge = std::pair<int, int>;
using EdgeList = std::vector<Edge>;
using MatrixXd = Eigen::MatrixXd;
using SparseMatrix_d = Eigen::SparseMatrix<double>;
using Trip_d = Eigen::Triplet<double>;

class MassSpring {
   public:
    MassSpring() = default;

    virtual ~MassSpring() = default;

    enum TimeIntegrator { IMPLICIT_EULER = 0, SEMI_IMPLICIT_EULER = 1 };

    MassSpring(const Eigen::MatrixXd &X, const EdgeSet &E);
    MassSpring(const Eigen::MatrixXd &X, const EdgeSet &E, Eigen::MatrixXi F);

    virtual void step();
    void reset();

    void dihedralConstraint(double mass_per_vertex);

    // energy related function
    virtual double computeEnergy(double stiffness);
    virtual Eigen::MatrixXd computeGrad(double stiffness);
    virtual Eigen::SparseMatrix<double> computeHessianSparse(double stiffness);

    // make matrix positive definite
    // Eigen::SparseMatrix<double> makeSPD(const Eigen::SparseMatrix<double> &A);
    bool checkSPD(const Eigen::SparseMatrix<double> &A);

    Eigen::MatrixXd getVelocity() const
    {
        return vel;
    }
    Eigen::MatrixXd getX() const
    {
        return X;
    }

    // Detect collision and compute the penalty-based collision force with given sphere
    Eigen::MatrixXd getSphereCollisionForce(Eigen::Vector3d center, double radius);

    bool set_dirichlet_bc_mask(const std::vector<bool>& mask);
    bool update_dirichlet_bc_vertices(const MatrixXd &control_vertices); 
    bool init_dirichlet_bc_vertices_control_pair(const MatrixXd &control_vertices,
                                      const std::vector<bool>& control_mask);

    // Simulation parameters
    double stiffness = 100.0;
    double damping = 0.995;
    enum TimeIntegrator time_integrator = IMPLICIT_EULER;
    double mass = 1.0;  // total mass of the mesh
    double h = 1e-2;    // time step
    Eigen::Vector3d gravity = { 0, 0, -9.8 };
    Eigen::Vector3d wind_ext_acc = { 0, 0, 0}; // (HW TODO) feel free to change the wind acceleration

    // (HW Optional) sphere collision parameters
    double collision_penalty_k = 10000.0;
    double collision_scale_factor = 1.1; 
    Eigen::Vector3f sphere_center = Eigen::Vector3f(0, -0.5, 0.2);
    double sphere_radius = 0.4;

    // Useful switches
    bool enable_sphere_collision = false;
    bool enable_time_profiling = false;
    bool enable_make_SPD = false;
    bool enable_check_SPD = false;
    bool enable_damping = true;
    bool enable_debug_output = false;

   protected:
    Eigen::MatrixXd init_X;  // For reset
    Eigen::MatrixXd X;
    Eigen::MatrixXd vel;
    EdgeSet E;
    Eigen::MatrixXi F;
    SpringMesh MyMesh;
    std::vector<double> E_rest_length;
    std::vector<bool>
        dirichlet_bc_mask;  // mask for marking fixed points (Dirichlet boundary condition)
    std::vector<std::pair<int, int>> dirichlet_bc_control_pair;
};
}  // namespace USTC_CG::node_mass_spring