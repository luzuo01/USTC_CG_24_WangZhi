#pragma once
#include "warpper.h"
#include <imgui.h>
#include <Eigen/Dense>


namespace USTC_CG
{
class RBF : public Warpper
{
public:
    RBF() = default;
    RBF(std::vector<ImVec2> sps, std::vector<ImVec2> eps): start_points_(sps), end_points_(eps), length(sps.size()) {}
    std::pair<int, int> warp(int x, int y, int width, int height) override;
    // initialize method is in charge of finding the matrix
    void initialize() override;

private:
    float Gaussian(float sigma, float t);
    std::vector<ImVec2> start_points_, end_points_;
    // This is the T we interpolate
    std::pair<int, int> transform(int x, int y);
    int length;
    Eigen::MatrixXf weights;
    float sigma = 200.0;
};
}