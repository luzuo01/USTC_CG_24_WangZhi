#include "parameterizer_arap.h"
#include <iostream>
#include <cmath>
#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "util_laplacian_solver.h"

namespace USTC_CG{

void Parameterizer_ARAP::initialize(shared_ptr<PolyMesh> our_mesh) {
    this -> our_mesh = our_mesh;
    backup_mesh = our_mesh;
    n_vertices = our_mesh -> n_vertices();
    // cot matrix doesn't rely on the boundary map
    // so we call assemble_A FIRST
    assemble_A();
    // init previous pos first
    for (const auto& vertex_handle : our_mesh -> vertices())
    {
        auto position = our_mesh -> point(vertex_handle);
        Vector2d vec = Vector2d(position[0], position[1]);
        previous_position[vertex_handle.idx()] = vec;
    }
    // init current_position with some 
    // kind of boundary mapping
    Laplace_solver solver;
    solver.detect(backup_mesh);
    vector<VertexHandle> boundary_vertices = solver.get_boundary();
    double length = solver.get_boundary_length();

    int numPoints = boundary_vertices.size();
    double angleIncrement = 2.0 * M_PI / numPoints;
    double radius = length / (2 * M_PI);
    double x_center = 0;
    double y_center = 0;

    for (int i = 0; i < numPoints; ++i) {
        double angle = i * angleIncrement;
        double x = x_center + radius * cos(angle);
        double y = y_center + radius * sin(angle);
        PolyMesh::Point point(x, y, 0);
        backup_mesh -> set_point(boundary_vertices[i], point);
    }
    solver.initialize(backup_mesh, 0);
    std::map<int, PolyMesh::Point> result = solver.solve();
    // update the internal points
    for (int i = 0; i < backup_mesh -> n_vertices(); i++)
    {
        Vector2d vec = Vector2d(result[i][0], result[i][1]);
        current_position[i] = vec;
    }
    // L will be inited at the first local_phase()
}


void Parameterizer_ARAP::local_phase() {
    // we build the cross variance matrix here 
    // then assign them to L
    // leave it to SVD() to further update L
    for (const auto& halfedge_handle : our_mesh -> halfedges()) {
        if (!our_mesh -> is_boundary(halfedge_handle)) {
            int i = halfedge_handle.from().idx();
            int j = halfedge_handle.to().idx();
            Matrix2d result;
            result.setZero();
            const auto& face_handle = halfedge_handle.face();
            const auto& start = face_handle.halfedge();
            SmartHalfedgeHandle he = start;
            do {
                int id1 = he.from().idx();
                int id2 = he.to().idx();
                Vector2d vecu = current_position[id1] - current_position[id2];
                Vector2d vecx = previous_position[id1] - previous_position[id2];
                int id3 = he.next().to().idx();
                double cot = cot_matrix[id2][id3];
                result += cot * vecu * vecx.transpose();
                he = he.next();
            } while (he != start);
            L[{i, j}] = result;
        } else {
            int i = halfedge_handle.from().idx();
            int j = halfedge_handle.to().idx();
            Matrix2d result;
            result.setZero();
            L[{i, j}] = result;
        }
    }
    SVD(); 
}

    
void Parameterizer_ARAP::global_phase(int phase_count) {
    if (solver.info() != Success) {
        throw std::runtime_error ("solver error");
    }
    vector<VectorXd> solutions(2);
    assemble_targets();
    solutions[0] = solver.solve(targets[0]);
    solutions[1] = solver.solve(targets[1]);
    for (int i = 0; i < n_vertices; i++)
    {
        previous_position[i] = current_position[i];
    }
    for (int i = 0; i < n_vertices; i++)
    {
        if (phase_count == 0) {
            Vector2d vec = Vector2d(solutions[0][i], solutions[1][i]);
            current_position[i] = vec;
        } else if (phase_count > 0 && i != 1 && i != n_vertices - 1) {
            Vector2d vec = Vector2d(solutions[0][i], solutions[1][i]);
            current_position[i] = vec;
        } else {
            Vector2d vec = Vector2d(solutions[0][i], solutions[1][i]);
            current_position[i] = 0.8 * vec + 0.2 * current_position[i];
        }
    }
}
void Parameterizer_ARAP::iterate(int num_iteration) {
    for (int i = 0; i < num_iteration; i++)
    {
        local_phase();
        std::cout << "local phase completed, round " << i + 1<< std::endl;
        global_phase(i);
        std::cout << "global phase completed, round " << i + 1<< std::endl;
        std::cout << "Phase " << i + 1<< " completed" << std::endl; 
        std::cout << " " << std::endl;
    }  
}
std::map<int, Vector2d> Parameterizer_ARAP::get_parameterization() {
    return current_position;
}

void Parameterizer_ARAP::assemble_A() {
    //build the cot matrix here
    cot_matrix.resize(our_mesh -> n_vertices());
    for (const auto& halfedge_handle : our_mesh -> halfedges()) {
        // first consider when the halfegde is not on the boundary
        if (!our_mesh -> is_boundary(halfedge_handle)) {
            const auto& hed1 = halfedge_handle.prev();
            const auto& hed2 = halfedge_handle.next().opp();
            const auto& vec1 = our_mesh -> point(hed1.to()) - our_mesh -> point(hed1.from());
            const auto& vec2 = our_mesh -> point(hed2.to()) - our_mesh -> point(hed2.from());
            double cos = vec1.dot(vec2) / (vec1.norm() * vec2.norm());
            double theta = acosf(cos);
            double cot = 1.0 / tan(theta);
            cot_matrix[halfedge_handle.from().idx()][halfedge_handle.to().idx()] = cot;
        } else {
            const auto& hed1 = halfedge_handle.opp().prev();
            const auto& hed2 = halfedge_handle.opp().next().opp();
            const auto& vec1 = our_mesh -> point(hed1.to()) - our_mesh -> point(hed1.from());
            const auto& vec2 = our_mesh -> point(hed2.to()) - our_mesh -> point(hed2.from());
            double cos = vec1.dot(vec2) / (vec1.norm() * vec2.norm());
            double theta = acosf(cos);
            double cot = 1.0 / tan(theta);
            cot_matrix[halfedge_handle.from().idx()][halfedge_handle.to().idx()] = cot;
        }
    }
    A.resize(our_mesh -> n_vertices(), our_mesh -> n_vertices());
    A.setZero();
    for (int i = 0; i < our_mesh -> n_vertices(); i++)
    {
        double sum = 0;
        for (int j : neighbors(i))
        {
            sum += cot_matrix[i][j] + cot_matrix[j][i];
            A.insert(i, j) = (-1) * (cot_matrix[i][j] + cot_matrix[j][i]);
        }
        A.insert(i, i) = sum;
    }
    solver.compute(A);
}
void Parameterizer_ARAP::assemble_targets() {
    targets.resize(2);
    VectorXd vec0 = assemble_target_at(0);
    VectorXd vec1 = assemble_target_at(1);
    targets[0] = vec0;
    targets[1] = vec1;
}
void Parameterizer_ARAP::SVD() {
    for (const auto& halfedge_handle : our_mesh -> halfedges()) {
        int id1 = halfedge_handle.from().idx();
        int id2 = halfedge_handle.to().idx();
        JacobiSVD<Matrix2d> svd(L[{id1, id2}], ComputeFullU | ComputeFullV);
        L[{id1, id2}] = svd.matrixU() * svd.matrixV();
    }
}

VectorXd Parameterizer_ARAP::assemble_target_at(int coord) {
    VectorXd rst(n_vertices);
    for (int i = 0; i < n_vertices; i++)
    {
        double sum = 0;
        for (int j : neighbors(i))
        {
            Vector2d vec = (cot_matrix[i][j] * L[{i, j}] + cot_matrix[j][i] * L[{j, i}]) * (previous_position[i] - previous_position[j]);
            sum += vec[coord];
        }
        rst[i] = sum;
    }
    return rst;
}

std::vector<int> Parameterizer_ARAP::neighbors(int i){
    vector<int> rst;
    VertexHandle vh = our_mesh -> vertex_handle(i);
    for (auto it = our_mesh -> vv_iter(vh); it.is_valid(); it++)
    {
        int x = it -> idx();
        rst.push_back(x);
    }
    return rst;
}
void Parameterizer_ARAP::initialize() {}
void Parameterizer_ARAP::iterate() {}
void Parameterizer_ARAP::global_phase() {};
}