#include "GCore/Components/MeshOperand.h"
#include "Nodes/node.hpp"
#include "Nodes/node_declare.hpp"
#include "Nodes/node_register.h"
#include "geom_node_base.h"
#include "utils/util_openmesh_bind.h"
#include "utils/util_laplacian_solver.h"
/*
** @brief HW4_TutteParameterization
**
** This file presents the basic framework of a "node", which processes inputs
** received from the left and outputs specific variables for downstream nodes to
** use.
** - In the first function, node_declare, you can set up the node's input and
** output variables.
** - The second function, node_exec is the execution part of the node, where we
** need to implement the node's functionality.
** - The third function generates the node's registration information, which
** eventually allows placing this node in the GUI interface.
**
** Your task is to fill in the required logic at the specified locations
** within this template, especially in node_exec.
*/

namespace USTC_CG::node_min_surf {
static void node_min_surf_declare(NodeDeclarationBuilder& b)
{
    // Input-1: Original 3D mesh with boundary
    b.add_input<decl::Geometry>("Input");
    // Input-2: Weights
    /*
    ** NOTE: You can add more inputs or outputs if necessary. For example, in some cases,
    ** additional information (e.g. other mesh geometry, other parameters) is required to perform
    ** the computation.
    **
    ** Be sure that the input/outputs do not share the same name. You can add one geometry as
    **
    **                b.add_input<decl::Geometry>("Input");
    **
    ** Or maybe you need a value buffer like:
            b.add_input<decl::Float1Buffer>("Weights");
    **
    */
    b.add_input<decl::Int>("Weight Type").min(0).max(2).default_val(0);
    // Output-1: Minimal surface with fixed boundary
    b.add_output<decl::Geometry>("Output");
}

static void node_min_surf_exec(ExeParams params)
{
    // Get the input from params
    auto input = params.get_input<GOperandBase>("Input");
    int weight_type = params.get_input<int>("Weight Type");
    // (TO BE UPDATED) Avoid processing the node when there is no input
    if (!input.get_component<MeshComponent>()) {
        throw std::runtime_error("Minimal Surface: Need Geometry Input.");
    }
    //throw std::runtime_error("Not implemented");

    /* ----------------------------- Preprocess -------------------------------
    ** Create a halfedge structure (using OpenMesh) for the input mesh. The
    ** half-edge data structure is a widely used data structure in geometric
    ** processing, offering convenient operations for traversing and modifying
    ** mesh elements.
    */
    auto halfedge_mesh = operand_to_openmesh(&input);
    Laplace_solver solver;
    solver.initialize(halfedge_mesh, weight_type);
    std::map<int, PolyMesh::Point> result = solver.solve();
    // update the internal points
    for (auto vertex : halfedge_mesh -> vertices())
    {
        halfedge_mesh -> set_point (vertex, result[vertex.idx()]);
    }
    /* ----------------------------- Postprocess ------------------------------
    ** Convert the minimal surface mesh from the halfedge structure back to
    ** GOperandBase format as the node's output.
    */
    auto operand_base = openmesh_to_operand(halfedge_mesh.get());

    // Set the output of the nodes
    params.set_output("Output", std::move(*operand_base));
}

static void node_register()
{
    static NodeTypeInfo ntype;

    strcpy(ntype.ui_name, "Minimal Surface");
    strcpy_s(ntype.id_name, "geom_min_surf");

    geo_node_type_base(&ntype);
    ntype.node_execute = node_min_surf_exec;
    ntype.declare = node_min_surf_declare;
    nodeRegisterType(&ntype);
}

NOD_REGISTER_NODE(node_register)
}  // namespace USTC_CG::node_min_surf
