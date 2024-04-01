#include "GCore/Components/MeshOperand.h"

#include "util_openmesh_bind.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/LU>

using namespace Eigen;
using namespace std;
using namespace OpenMesh;

namespace USTC_CG {
class Laplace_solver {
    public:
        Laplace_solver() = default;
        void initialize(shared_ptr<PolyMesh> our_mesh, int weight_type);
        void detect(shared_ptr<PolyMesh> our_mesh);
        std::map<int, PolyMesh::Point> solve();
        vector<VertexHandle> get_boundary();
        double get_boundary_length();
    private:
        std::map<int, int> relabel;
        std::map<int, int> inv_relabel;
        int n_internal_vertex;
        int n_vertices;
        std::vector<std::vector<double>> weights;
        SparseMatrix<double> A;
        std::vector<VectorXd> targets;
        shared_ptr<PolyMesh> our_mesh;
        int dimension = 3;
        vector<VertexHandle> boundary_vertices;
        double boundary_length;
        // the following method assigns values to the two
        // maps and the number of inner points
        void boundary_detect(shared_ptr<PolyMesh> our_mesh);
        std::vector<int> neighbors_relabelled(int i);

        // a few different ways to set weights
        void set_uniform_weights();
        void set_cotangent_weights();
        void assemble_A();
        void assemble_targets();
        VectorXd assemble_target_at(int i);
   };
}

