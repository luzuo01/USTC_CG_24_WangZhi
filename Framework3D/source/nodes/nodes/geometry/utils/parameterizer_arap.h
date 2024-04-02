#include "GCore/Components/MeshOperand.h"
#include "util_openmesh_bind.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/LU>
#include "parameterizer.h"

using namespace Eigen;
using namespace std;
using namespace OpenMesh;

namespace USTC_CG {
class Parameterizer_ARAP : public Parameterizer{
    public:
    ~Parameterizer_ARAP() = default;
    Parameterizer_ARAP() = default;
    void initialize() override;
    void initialize(shared_ptr<PolyMesh> our_mesh);
    void iterate() override;
    void iterate( int num_iteration);
    std::map<int, Vector2d> get_parameterization();

    
    private:
    SparseMatrix<double> A;
    std::vector<VectorXd> targets;

    void local_phase() override;
    void global_phase() override;
    void global_phase(int phase_count);   
    void assemble_A();
    void assemble_targets();
    void SVD();
    VectorXd assemble_target_at(int i);
    std::vector<int> neighbors(int i);


    std::vector<std::map<int, double>> cot_matrix;
    std::map<int, Vector2d> current_position;
    std::map<int, Vector2d> previous_position;
    std::map<std::pair<int, int>, MatrixX2d> L;
    shared_ptr<PolyMesh> our_mesh;
    shared_ptr<PolyMesh> backup_mesh;
    SparseLU<SparseMatrix<double>> solver;
    int n_vertices = 0;
   };
}

