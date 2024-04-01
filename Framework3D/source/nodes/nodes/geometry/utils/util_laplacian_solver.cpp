#include <iostream>
#include <cmath>
#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "util_laplacian_solver.h"

namespace USTC_CG {
void Laplace_solver::initialize(shared_ptr<PolyMesh> our_mesh, int weight_type) {
    this->our_mesh = our_mesh;
    boundary_detect(our_mesh);
    if (weight_type == 0) {
        set_uniform_weights();
    } else if (weight_type == 1) {
        set_uniform_weights();
        //set_cotangent_weights();
    } else {
        set_uniform_weights();
    }
    assemble_A();
    std::cout << "the matrix A is assembled" << std::endl;
    assemble_targets();
    std::cout << "target vectors are assembled" << std::endl;
    std::cout<< "init success" << std::endl;
}
void Laplace_solver::detect(shared_ptr<PolyMesh> our_mesh) {
    this->our_mesh = our_mesh;
    boundary_detect(our_mesh);
    std::cout << "boundary detected" << std::endl;
}
std::map<int, PolyMesh::Point> Laplace_solver::solve() {
    std::map<int, PolyMesh::Point> rst;
    SparseLU<SparseMatrix<double>> solver;
    solver.compute(A);
    if (solver.info() != Success) {
        throw std::runtime_error ("solver error");
    }
    vector<VectorXd> solutions(dimension);
    solutions[0] = solver.solve(targets[0]);
    solutions[1] = solver.solve(targets[1]);
    solutions[2] = solver.solve(targets[2]);
    for (int i = 0; i < n_internal_vertex; i++) {
        PolyMesh::Point point(solutions[0][i], solutions[1][i], solutions[2][i]);
        rst[relabel[i]] = point;
    }
    for (int j = n_internal_vertex; j < n_vertices; j++) {
        VertexHandle boundary_point = our_mesh -> vertex_handle(relabel[j]);
        auto position = our_mesh -> point(boundary_point);
        PolyMesh::Point point (position[0], position[1], position[2]);
        rst[relabel[j]] = point;
    }
    std::cout << "solve completed" << std::endl;
    return rst;
}

vector<VertexHandle> Laplace_solver::get_boundary() {
    return boundary_vertices;
}
double Laplace_solver::get_boundary_length(){
    return boundary_length;
}

void Laplace_solver::boundary_detect(shared_ptr<PolyMesh> our_mesh) {
    n_vertices = our_mesh -> n_vertices();
    size_t maxItr(our_mesh -> n_halfedges());
    vector<SmartHalfedgeHandle> boundary;
    boundary_length = 0.0;
    vector<bool> is_boundary_vertex(our_mesh -> n_vertices(), false);
    // iterate through all half edges
    for (const auto& he_itr : our_mesh -> halfedges()) 
    {
        if (false == our_mesh -> is_boundary(he_itr))
            continue;
        // boundry found                
        boundary.push_back(he_itr); 
        size_t counter = 1;
        auto next_he = our_mesh->next_halfedge_handle(he_itr);
        while ( (next_he != he_itr) && counter < maxItr)
        {
            assert(our_mesh->is_boundary(next_he));
            next_he = our_mesh->next_halfedge_handle(next_he);
            counter++;
        }
    }
    // set is_boundary_vertex via boundary
    for (SmartHalfedgeHandle edge : boundary)
    {
        is_boundary_vertex[edge.to().idx()] = true;
        is_boundary_vertex[edge.from().idx()] = true;
    }
    // set the relabel map
    int counter = 0;
    for (const auto& vertex_handle : our_mesh -> vertices()) {
        if (!is_boundary_vertex[vertex_handle.idx()]) {
            relabel[counter] = vertex_handle.idx();
            counter++;
        }
    }
    // now 0 - counter correspond to internal vertices
    n_internal_vertex = counter;
    for (const auto& vertex_handle : our_mesh -> vertices()) {
        if (is_boundary_vertex[vertex_handle.idx()]) {
            relabel[counter] = vertex_handle.idx();
            counter++;
            //boundary_vertices.push_back(vertex_handle);
        }
    }
    // inverse the relabel map
    for (auto pair : relabel)
    {
        inv_relabel[pair.second] = pair.first;
    }
    // compute the boundary_length 
    for (auto edge : boundary)
    {
        const auto& v1 = edge.to();
        const auto& v2 = edge.from();
        boundary_length += (our_mesh -> point(v1) - our_mesh -> point(v2)).norm();
    }
    // obtain the boundary_vertices IN ORDER
    bool boundary_traverse_completed = false;
    for (const auto& he_itr : our_mesh -> halfedges()) 
    {
        if (false == our_mesh -> is_boundary(he_itr))
            continue;
        // boundry found    
        if (!boundary_traverse_completed) {
            boundary_vertices.push_back(he_itr.from());  
            size_t counter = 1;
            auto next_he = our_mesh->next_halfedge_handle(he_itr);
            while ( (next_he != he_itr) && counter < maxItr)
            {
                assert(our_mesh->is_boundary(next_he));
                boundary_vertices.push_back(next_he.from());
                next_he = our_mesh->next_halfedge_handle(next_he);
                counter++;
            }
            boundary_traverse_completed = true;
        }       
    }
}

std::vector<int> Laplace_solver::neighbors_relabelled(int i){
    vector<int> rst;
    VertexHandle vh = our_mesh -> vertex_handle(relabel[i]);
    for (auto it = our_mesh -> vv_iter(vh); it.is_valid(); it++)
    {
        int x = it -> idx();
        rst.push_back(inv_relabel[x]);
    }
    return rst;
}

// a few different ways to set weights
void Laplace_solver::set_uniform_weights() {
    weights.resize(n_internal_vertex, std::vector<double> (n_vertices));
    for (size_t i = 0; i < n_internal_vertex; i++)
    {
        vector<int> neighbors = neighbors_relabelled(i);
        float n_neighbors = neighbors.size();
        for (int index : neighbors)
        {
            weights[i][index] = (1.0) / (n_neighbors + 1e-7);
        }
    }
}
void Laplace_solver::set_cotangent_weights() {
    weights.resize(n_internal_vertex, std::vector<double> (n_vertices));
    for (const auto& vertex_handle : our_mesh -> vertices())
    {
        int i = inv_relabel[vertex_handle.idx()];
        if (i < n_internal_vertex) {
            for (const auto& halfedge_handle : vertex_handle.outgoing_halfedges()) {
                const auto& v = halfedge_handle.to();
                const auto& v1 = halfedge_handle.prev().opp().to();
                const auto& v2 = halfedge_handle.opp().next().to();
                int j = inv_relabel[v.idx()];
                const auto& vec11 = our_mesh -> point(v) - our_mesh -> point(v1);
                const auto& vec12 = our_mesh -> point(vertex_handle) - our_mesh -> point(v1);
                double cos1 = vec11.dot(vec12) / (vec11.norm() * vec12.norm());
                double theta1 = acosf(cos1);
                double sin1 = sin(theta1);
                double cot1 = cos1 / (sin1 + 1e-6);
                const auto& vec21 = our_mesh -> point(v) - our_mesh -> point(v2);
                const auto& vec22 = our_mesh -> point(vertex_handle) - our_mesh -> point(v2);
                double cos2 = vec21.dot(vec22) / (vec21.norm() * vec22.norm());
                double theta2 = acosf(cos2);
                double sin2 = sin(theta2);
                double cot2 = cos2 / (sin2 + 1e-6);
                weights[i][j] = cot1 + cot2;
            }
        }
    }
    
}
// obtain the matrix and target vectors
void Laplace_solver::assemble_A() {
    A.resize(n_internal_vertex, n_internal_vertex);
    A.setZero();
    for (size_t i = 0; i < n_internal_vertex; i++) {
        A.insert(i, i) = 1.0;
        vector<int> neighbors = neighbors_relabelled(i);
        int n_neighbors = neighbors.size();
        for (int index : neighbors)
        {
            if (index < n_internal_vertex && index != i) {
                A.insert(i, index) = -weights[i][index];
            }
        }
    }
}

VectorXd Laplace_solver::assemble_target_at(int coord) {
    VectorXd b(n_internal_vertex);
    for (int i = 0; i < n_internal_vertex; i++)
    {
        double sum = 0;
        for (int j = n_internal_vertex; j < n_vertices; j++)
        {
            VertexHandle boundary_point = our_mesh -> vertex_handle(relabel[j]);
            auto position = our_mesh -> point(boundary_point);
            double u_j = position[coord];
            sum += weights[i][j] * u_j;
        }
        b[i] = sum;
    }
    return b;
}
void Laplace_solver::assemble_targets() {
    for (int coord = 0; coord < dimension; coord++)
    {
        targets.push_back(assemble_target_at(coord));
    } 
}
}
