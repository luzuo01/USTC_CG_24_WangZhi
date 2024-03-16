#pragma once

#include "comp_source_image.h"
#include "view/comp_image.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
using namespace Eigen;
namespace USTC_CG
{
// reuse the pixel structure 
struct Pixel {
    int x;
    int y;
    std::vector<unsigned char> values; // Assume grayscale value for simplicity
    Pixel(int x_, int y_, const std::vector<unsigned char>& value_) : x(x_), y(y_), values(value_) {}
    Pixel() = default;
};
class CompTargetImage : public ImageEditor
{
   public:
    enum CloneType
    {
        kDefault = 0,
        kPaste = 1,
        kSeamless = 2,
        kMixed = 3,
    };

    explicit CompTargetImage(
        const std::string& label,
        const std::string& filename);
    virtual ~CompTargetImage() noexcept = default;

    void draw() override;
    // Bind the source image component
    void set_source(std::shared_ptr<CompSourceImage> source);
    // Enable realtime updating
    void set_realtime(bool flag);
    void restore();

    // HW3_TODO: Add more types of cloning
    void set_paste();
    void set_seamless();
    void set_mixed();
    // The clone function
    void clone();

   private:
    // Store the original image data
    std::shared_ptr<Image> back_up_;
    // Source image
    std::shared_ptr<CompSourceImage> source_image_;
    CloneType clone_type_ = kDefault;

    // helper functions for Possion solver
    // the following two return boundary and gradient conditions
    // a vector of Pixels is considered to be a function/graph
    std::vector<Pixel> boundary_condition(std::shared_ptr<Image>& source_ptr);

    std::vector<Pixel> gradient_condition(std::shared_ptr<Image>& source_ptr);

    std::vector<Pixel> neighbors(Pixel);

    // the followings return the target vector and the sparse matrix
    VectorXf get_target_vector_at_channel(std::vector<Pixel> boundary,  std::vector<Pixel> region, int i);
    std::vector<VectorXf> get_target_vector(std::vector<Pixel> boundary,  std::vector<Pixel> region);

    VectorXf mixed_get_target_vector_at_channel(std::vector<Pixel> boundary,  std::vector<Pixel> region, int i);
    std::vector<VectorXf> mixed_get_target_vector(std::vector<Pixel> boundary,  std::vector<Pixel> region);


    MatrixXf get_sparse_matrix(std::vector<Pixel> boundary,  std::vector<Pixel> region);

    std::pair<int, int> inverse(int x, int y);
    bool is_neighbor(Pixel p, Pixel q);

    // we maintain a map from coordinates to number
    // so that we can build the sparse matrix faster 
    // without having to iterate through region with a counter 
    // just have to look up in the map
    // this way it takes O(nlogn) to build the sparse matrix
    std::map<std::pair<int, int>, int> region_map;

    // then the two methods above are used in the Eigen sparse solver
    std::vector<Pixel> possion_solver(std::vector<VectorXf> vector, std::vector<Pixel> region);
    // accelerate the possion solver
    std::vector<Pixel> possion_solver_pre_decomposition
    (MatrixXf matrix, std::vector<VectorXf> vector);

    // keep track of the modified region and boundary
    std::vector<std::vector<bool>> region;
    std::vector<std::vector<bool>> boundary;

    //pre-decomposed solver
    SimplicialLDLT<SparseMatrix<float>> solver;

    Eigen::SparseMatrix<float> sparseMatrix;

    ImVec2 mouse_position_;
    bool edit_status_ = false;
    bool flag_realtime_updating = false;
};

}  // namespace USTC_CG