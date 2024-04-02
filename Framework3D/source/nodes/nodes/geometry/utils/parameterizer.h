#include "GCore/Components/MeshOperand.h"
#include "util_openmesh_bind.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/LU>

using namespace Eigen;
using namespace std;
using namespace OpenMesh;

namespace USTC_CG {
class Parameterizer {
    public:
    virtual ~Parameterizer() = default;
    Parameterizer() = default;
    virtual void initialize() = 0;
    virtual void local_phase() = 0;
    virtual void global_phase() = 0;
    virtual void iterate() = 0;
   };
}

