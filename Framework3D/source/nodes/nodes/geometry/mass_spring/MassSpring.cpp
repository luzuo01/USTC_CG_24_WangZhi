#include "MassSpring.h"
#include <iostream>

namespace USTC_CG::node_mass_spring {
MassSpring::MassSpring(const Eigen::MatrixXd& X, const EdgeSet& E)
{
    this->X = this->init_X = X;
    this->vel = Eigen::MatrixXd::Zero(X.rows(), X.cols());
    this->E = E;

    std::cout << "number of edges: " << E.size() << std::endl;
    std::cout << "init mass spring" << std::endl;

    // Compute the rest pose edge length
    for (const auto& e : E) {
        Eigen::Vector3d x0 = X.row(e.first);
        Eigen::Vector3d x1 = X.row(e.second);
        this->E_rest_length.push_back((x0 - x1).norm());
    }

    // Initialize the mask for Dirichlet boundary condition
    dirichlet_bc_mask.resize(X.rows(), false);

    // (HW_TODO) Fix two vertices, feel free to modify this 
    unsigned n_fix = sqrt(X.rows());  // Here we assume the cloth is square
    dirichlet_bc_mask[0] = true;
    dirichlet_bc_mask[n_fix - 1] = true;
}

void MassSpring::step()
{
    Eigen::Vector3d acceleration_ext = gravity + wind_ext_acc;

    unsigned n_vertices = X.rows();

    const auto I = Eigen::MatrixXd::Identity(3 * n_vertices, 3 * n_vertices);

    // The reason to not use 1.0 as mass per vertex: the cloth gets heavier as we increase the resolution
    double mass_per_vertex =
        mass / n_vertices; 

    //----------------------------------------------------
    // (HW Optional) Bonus part: Sphere collision
    Eigen::MatrixXd acceleration_collision =
        getSphereCollisionForce(sphere_center.cast<double>(), sphere_radius);
    //----------------------------------------------------

    if (time_integrator == IMPLICIT_EULER) {
        // Implicit Euler
        TIC(step)

        // (HW TODO) 
        auto H_elastic = computeHessianSparse(stiffness);  // size = [nx3, nx3]
        // then we get H_g from H_elastic
        Eigen::SparseMatrix<double> H_g = mass_per_vertex * I / h / h + H_elastic;
        // translate the diagonal to ensure positive definite-ness
        double epsilon = 0.02;
        for (size_t i = 0; i < 3 * n_vertices; i++)
        {
            if (!dirichlet_bc_mask[i])
            {
                H_g.coeffRef(i, i) += epsilon;
            }
        }
        
        if (checkSPD(H_g))
        {
        // compute Y 
        Eigen::MatrixXd Y = X + h * vel;
        Y.rowwise() += h * h * acceleration_ext.transpose();

        // compute grad_g
        Eigen::MatrixXd grad_g = mass_per_vertex * (X - Y) / h / h + computeGrad(stiffness);

        // then we fix certain points 
        for (size_t i = 0; i < n_vertices; i++)
        {
            if (dirichlet_bc_mask[i])
            {
                Eigen::Vector3d vec = Eigen::Vector3d::Zero();
                grad_g.row(i) =  vec;
            }
        }

        // flatten grad_g
        VectorXd grad_g_flatten = flatten(grad_g);
        VectorXd X_flatten = flatten(X);


        // Solve Newton's search direction with linear solver 
        SparseLU<SparseMatrix<double>> solver;
        solver.compute(H_g);
        Eigen::VectorXd delta_X_flatten = solver.solve((-1.0) * grad_g_flatten);

        // update X and vel 
        X_flatten += delta_X_flatten;
        Eigen::MatrixXd X_prev = X;
        X = unflatten(X_flatten);
        vel = (X - X_prev) / h;
        } else {
            std::cerr << "The Hessian is not positive definite!" << std::endl;
            return;
        }
        
        TOC(step)
    }
    else if (time_integrator == SEMI_IMPLICIT_EULER) {

        // Semi-implicit Euler
        Eigen::MatrixXd acceleration = -computeGrad(stiffness) / mass_per_vertex;
        acceleration.rowwise() += acceleration_ext.transpose();
        // modify acceleration first to fix certain points 
        for (size_t i = 0; i < X.rows(); i++)
        {
            if (dirichlet_bc_mask[i])
            {
                Eigen::Vector3d vec = Eigen::Vector3d::Zero();
                acceleration.row(i) = vec;
            }
        }
        // -----------------------------------------------
        // (HW Optional)
        if (enable_sphere_collision) {
            acceleration += acceleration_collision;
        }
        // -----------------------------------------------
        vel += h * acceleration;
        X += h * vel;
        vel *= damping;
    }
    else {
        std::cerr << "Unknown time integrator!" << std::endl;
        return;
    }
}

// There are different types of mass spring energy:
// For this homework we will adopt Prof. Huamin Wang's energy definition introduced in GAMES103
// course Lecture 2 E = 0.5 * stiffness * sum_{i=1}^{n} (||x_i - x_j|| - l)^2 There exist other
// types of energy definition, e.g., Prof. Minchen Li's energy definition
// https://www.cs.cmu.edu/~15769-f23/lec/3_Mass_Spring_Systems.pdf
double MassSpring::computeEnergy(double stiffness)
{
    double sum = 0.;
    unsigned i = 0;
    for (const auto& e : E) {
        auto diff = X.row(e.first) - X.row(e.second);
        auto l = E_rest_length[i];
        sum += 0.5 * stiffness * std::pow((diff.norm() - l), 2);
        i++;
    }
    return sum;
}

Eigen::MatrixXd MassSpring::computeGrad(double stiffness)
{
    Eigen::MatrixXd g = Eigen::MatrixXd::Zero(X.rows(), X.cols());
    unsigned i = 0;
    for (const auto& e : E) {
        auto diff = X.row(e.first) - X.row(e.second);
        auto l = E_rest_length[i];
        double length = diff.norm();
        double scale = stiffness * (length - l) / length;
        auto vec = scale * diff;
        g.row(e.first) += vec;
        g.row(e.second) -= vec;
        i++;
    }
    return g;
}

Eigen::SparseMatrix<double> MassSpring::computeHessianSparse(double stiffness)
{
    unsigned n_vertices = X.rows();
    Eigen::SparseMatrix<double> H(n_vertices * 3, n_vertices * 3);

    unsigned i = 0;
    auto k = stiffness;
    const auto I = Eigen::MatrixXd::Identity(3, 3);
    for (const auto& e : E) {
        // the hessians of e.first and e.second are both obtained from
        // the hessian of e
        MatrixXd Hessian_e = Eigen::MatrixXd::Zero(3, 3);
        Eigen::Vector3d diff = X.row(e.first) - X.row(e.second);
        auto l = E_rest_length[i];
        double length = diff.norm();
        Hessian_e += stiffness * diff * diff.transpose() / length / length;
        if (length >= l)
        {
            Hessian_e += stiffness * (1 - l / length) * (I - diff * diff.transpose() / length / length);
        }
        // then we update H
        for (size_t i = 3 * e.first; i < 3 * e.first + 3; i++)
        {
            for (size_t j = 3 * e.first; j < 3 * e.first + 3; j++)
            {
                double val = Hessian_e.coeff(i - 3 * e.first, j - 3 * e.first);
                H.coeffRef(i, j) += val;
            } 
        }
        for (size_t i = 3 * e.second; i < 3 * e.second + 3; i++)
        {
            for (size_t j = 3 * e.second; j < 3 * e.second + 3; j++)
            {
                double val = Hessian_e.coeff(i - 3 * e.second, j - 3 * e.second);
                H.coeffRef(i, j) += val;
            } 
        }
        for (size_t i = 3 * e.first; i < 3 * e.first + 3; i++)
        {
            for (size_t j = 3 * e.second; j < 3 * e.second + 3; j++)
            {
                double val = Hessian_e.coeff(i - 3 * e.first, j - 3 * e.second);
                H.coeffRef(i, j) -= val;
            } 
        }
        for (size_t i = 3 * e.second; i < 3 * e.second + 3; i++)
        {
            for (size_t j = 3 * e.first; j < 3 * e.first + 3; j++)
            {
                double val = Hessian_e.coeff(i - 3 * e.second, j - 3 * e.first);
                H.coeffRef(i, j) -= val;
            } 
        }
        i++;
    }
    // fix points
    for (int i = 0; i < n_vertices; i++)
    {
        if (dirichlet_bc_mask[i])
        {
            for (size_t a = 0; a < 3; a++)
            {
                for (size_t b = 0; b < 3; b++)
                {
                    H.coeffRef(3 * i + a, 3 * i + b) = I.coeff(a, b);
                }
                
            }
        }
    }
    

    // try the simplest way to make it positive definite
    H.makeCompressed();
    return H;
}


bool MassSpring::checkSPD(const Eigen::SparseMatrix<double>& A)
{
    // Eigen::SimplicialLDLT<SparseMatrix_d> ldlt(A);
    // return ldlt.info() == Eigen::Success;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
    auto eigen_values = es.eigenvalues();
    return eigen_values.minCoeff() >= 1e-10;
}

void MassSpring::reset()
{
    std::cout << "reset" << std::endl;
    this->X = this->init_X;
    this->vel.setZero();
}

// ----------------------------------------------------------------------------------
// (HW Optional) Bonus part
Eigen::MatrixXd MassSpring::getSphereCollisionForce(Eigen::Vector3d center, double radius)
{
    Eigen::MatrixXd force = Eigen::MatrixXd::Zero(X.rows(), X.cols());
    for (int i = 0; i < X.rows(); i++) {
       // (HW Optional) Implement penalty-based force here 
    }
    return force;
}
// ----------------------------------------------------------------------------------
 
bool MassSpring::set_dirichlet_bc_mask(const std::vector<bool>& mask)
{
	if (mask.size() == X.rows())
	{
		dirichlet_bc_mask = mask;
		return true;
	}
	else
		return false;
}

bool MassSpring::update_dirichlet_bc_vertices(const MatrixXd &control_vertices)
{
   for (int i = 0; i < dirichlet_bc_control_pair.size(); i++)
   {
       int idx = dirichlet_bc_control_pair[i].first;
	   int control_idx = dirichlet_bc_control_pair[i].second;
	   X.row(idx) = control_vertices.row(control_idx);
   }

   return true; 
}

bool MassSpring::init_dirichlet_bc_vertices_control_pair(const MatrixXd &control_vertices,
    const std::vector<bool>& control_mask)
{
    
	if (control_mask.size() != control_vertices.rows())
			return false; 

   // TODO: optimize this part from O(n) to O(1)
   // First, get selected_control_vertices
   std::vector<VectorXd> selected_control_vertices; 
   std::vector<int> selected_control_idx; 
   for (int i = 0; i < control_mask.size(); i++)
   {
       if (control_mask[i])
       {
			selected_control_vertices.push_back(control_vertices.row(i));
            selected_control_idx.push_back(i);
		}
   }

   // Then update mass spring fixed vertices 
   for (int i = 0; i < dirichlet_bc_mask.size(); i++)
   {
       if (dirichlet_bc_mask[i])
       {
           // O(n^2) nearest point search, can be optimized
           // -----------------------------------------
           int nearest_idx = 0;
           double nearst_dist = 1e6; 
           VectorXd X_i = X.row(i);
           for (int j = 0; j < selected_control_vertices.size(); j++)
           {
               double dist = (X_i - selected_control_vertices[j]).norm();
               if (dist < nearst_dist)
               {
				   nearst_dist = dist;
				   nearest_idx = j;
			   }
           }
           //-----------------------------------------
           
		   X.row(i) = selected_control_vertices[nearest_idx];
           dirichlet_bc_control_pair.push_back(std::make_pair(i, selected_control_idx[nearest_idx]));
	   }
   }

   return true; 
}

}  // namespace USTC_CG::node_mass_spring

