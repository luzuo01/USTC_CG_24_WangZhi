#include "RBF.h"
#include <cmath>
#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

namespace USTC_CG
{
// there are many possible choices for the actual function used here, we just call it Gaussian for simplicity
float RBF::Gaussian(float sigma, float t){
    if (t >= 0 && t <= sigma) {
        return 1 - std::pow(t / sigma, 2) * (3 - 2 * t / sigma);
    }
    return 0;
}

void RBF::initialize(){
  //initialize the matrix first
  MatrixXf Gauss_matrix(length, length);
  for (size_t row = 0; row < length; row++) {
	for (size_t col = 0; col < length; col++) {
        float distance = std::sqrt(
        std::pow(start_points_[row].x - start_points_[col].x, 2)
        + std::pow(start_points_[row].y - start_points_[col].y, 2)
        );
        Gauss_matrix(row, col) = Gaussian(sigma, distance);
    }
  }
  // the target vector, we choose the affine function to be identity for simplicity
  MatrixXf y(length, 2);
  for (size_t row = 0; row < length; row++) {
    y(row, 0) = end_points_[row].x - start_points_[row].x;
    y(row, 1) = end_points_[row].y - start_points_[row].y;
  }
  // solve the linear eq to obtain weights
   weights = Gauss_matrix.colPivHouseholderQr().solve(y);
}
std::pair<int, int> RBF::transform(int x, int y){
    float target_x = x;
    float target_y = y;
    for (int row = 0; row < length; row++)
    {
        float distance = std::sqrt(
        std::pow(x - start_points_[row].x, 2)
        + std::pow(y - start_points_[row].y, 2)
        );
        target_x += weights(row, 0) * Gaussian(sigma, distance);
        target_y += weights(row, 1) * Gaussian(sigma, distance);
    }
    int result_x = static_cast<int> (target_x);
    int result_y = static_cast<int> (target_y);
    return {result_x, result_y};
}

std::pair<int, int> RBF::warp(int x, int y, int width, int height) {
    assert( x >= 0 && x < width && y >= 0 && y < height);
    return transform(x, y);
}
}