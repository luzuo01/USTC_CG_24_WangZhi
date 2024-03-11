#pragma once
#include "warpper.h"
#include <imgui.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/NumericalDiff>
#include <unsupported/Eigen/NonLinearOptimization>
#include <annoylib.h>

using namespace Eigen;

namespace USTC_CG
{
class IDW : public Warpper
{
public:
    IDW() = default;
    IDW(std::vector<ImVec2> sps, std::vector<ImVec2> eps): start_points_(sps), end_points_(eps), length(sps.size()) {
        matrices.resize(length);
        for (int i = 0; i < length; i++)
            matrices[i].setZero();
    }
    std::pair<int, int> warp(int x, int y, int width, int height) override;
    std::pair<int, int> transform(int x, int y);
    void initialize() override;

private:
    std::vector<ImVec2> start_points_, end_points_;
    float weight_function(int i, int x, int y);
    float local_weight(int i, int x, int y);
    std::pair<int, int> local_appro (int i, int x, int y);
    Vector2f local_appro_with_matrix (int i, int x, int y, Matrix2f M);
    int length;
    struct ObjectiveFunction;
    std::vector<Matrix2f> matrices;
    float R = 2500; //range of influence
};
}