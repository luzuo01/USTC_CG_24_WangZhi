#include "IDW.h"
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/NumericalDiff>
#include <unsupported/Eigen/NonLinearOptimization>
#include <annoylib.h>

using namespace Eigen;

namespace USTC_CG
{

float IDW::local_weight(int i, int x, int y) {
    float distance = std::sqrt(
        std::pow(start_points_[i].x - x, 2) + std::pow(start_points_[i].y - y, 2));
    float numerator = (R - distance > 0) ? R - distance : 0;
    float denom = R * distance + 1e-5;
    return std::pow(numerator / denom, 3);
}

float IDW::weight_function(int i, int x, int y) {
    float denom = 0;
    for (int j = 0; j < length; j++)
    {
        denom += local_weight(j, x, y);
    }
    return local_weight(i, x, y) / denom;
}

std::pair<int, int> IDW::local_appro (int i, int x, int y) {
    Vector2f res = local_appro_with_matrix(i, x, y, matrices[i]);
    return {res(0,0), res(1,0)};
}

Vector2f IDW::local_appro_with_matrix (int i, int x, int y, Matrix2f M) {
    Vector2f s;
    s << x, y;
    Vector2f p_i;
    p_i << start_points_[i].x, start_points_[i].y;
    Vector2f q_i;
    q_i << end_points_[i].x, end_points_[i].y;
    Vector2f res = q_i + M * (s - p_i);
    return res;
}

std::pair<int, int> IDW::transform(int x, int y) {
    float presult_x = 0;
    float presult_y = 0;
    for (int i = 0; i < length; i++) {
        presult_x += weight_function(i, x, y) * local_appro(i, x, y).first;
        presult_y += weight_function(i, x, y) * local_appro(i, x, y).second;
    }
    int result_x = static_cast<int> (presult_x);
    int result_y = static_cast<int> (presult_y);
    return {result_x, result_y};
}

// a generic functor template
template<typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
struct Functor
{

  // Information that tells the caller the numeric type (eg. double) and size (input / output dim)
  typedef _Scalar Scalar;
  enum { // Required by numerical differentiation module
      InputsAtCompileTime = NX,
      ValuesAtCompileTime = NY
  };

  // Tell the caller the matrix sizes associated with the input, output, and jacobian
  typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
  typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
  typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

  // Local copy of the number of inputs
  int m_inputs, m_values;

  // Two constructors:
  Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
  Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

  // Get methods for users to determine function input and output dimensions
  int inputs() const { return m_inputs; }
  int values() const { return m_values; }

};

//our special functor
struct IDW::ObjectiveFunction : Functor<float>
{
    IDW& current;
    ObjectiveFunction(IDW& wassup, int index): current(wassup), index(index) {} 
    int index;
    int operator()(const VectorXf& x, VectorXf& fvec) const {
        // Extract the elements of the n 2x2 matrices from the input vector
        int n = current.length;
        for (int i = 0; i < n && 4 * i + 3 < x.size(); i++) {
            current.matrices[i] << x[i * 4], x[i * 4 + 1], x[i * 4 + 2], x[i * 4 + 3];
        }
        // Energy / loss
        for (int i = 0; i < n; ++i) {
            if (i != index) {
            Vector2f eps_i {current.end_points_[i].x, current.end_points_[i].y};
            Vector2f v = current.local_appro_with_matrix(index, current.start_points_[i].x,
            current.start_points_[i].y, current.matrices[index]) - eps_i;
            float norm = v.norm();
            fvec[i] = current.weight_function(index, current.start_points_[i].x,
            current.start_points_[i].y) * norm;
            } 
        }
        return 0;
    }
    int df(const Eigen::VectorXf &x, Eigen::MatrixXf &fjac) const {
        float epsilon;
        epsilon = 1e-8f;
        for (int i = 0; i < x.size(); i++) {
            Eigen::VectorXf xPlus(x);
            xPlus(i) += epsilon;
            Eigen::VectorXf xMinus(x);
            xMinus(i) -= epsilon;

            Eigen::VectorXf fvecPlus(values());
            operator()(xPlus, fvecPlus);

            Eigen::VectorXf fvecMinus(values());
            operator()(xMinus, fvecMinus);

            Eigen::VectorXf fvecDiff(values());
            fvecDiff = (fvecPlus - fvecMinus) / (2.0f * epsilon);

            fjac.block(0, i, values(), 1) = fvecDiff;
        }
        return 0;
    }
    int values() const { return current.length;}
};

void IDW:: initialize(){
    VectorXf x(length - 1);
    x.setZero(); 
    for (int i = 0; i < length; i++) {
    ObjectiveFunction functor(*this, i);
    NumericalDiff<ObjectiveFunction> numDiff(functor);
    LevenbergMarquardt<NumericalDiff<ObjectiveFunction>, float> lm(numDiff);
    lm.parameters.maxfev = 1000; // Maximum number of iterations
    lm.parameters.xtol = 1.0e-7; // tolerance 
    // Minimize the objective function
    lm.minimize(x);
    }
    x.resize(4 * length);
}

std::pair<int, int> IDW::warp(int x, int y, int width, int height) {
    assert( x >= 0 && x < width && y >= 0 && y < height);
    return transform(x, y);
}

}