#pragma once

#include "view/comp_image.h"

namespace USTC_CG
{
// Define struct for pixel
struct Pixel {
    int x;
    int y;
    std::vector<unsigned char> values; // Assume grayscale value for simplicity
    Pixel(int x_, int y_, const std::vector<unsigned char>& value_) : x(x_), y(y_), values(value_) {}
};
// Image component for warping and other functions
class CompWarping : public ImageEditor
{
   public:
    explicit CompWarping(const std::string& label, const std::string& filename);
    virtual ~CompWarping() noexcept = default;

    void draw() override;

    // Simple edit functions
    void invert();
    void mirror(bool is_horizontal, bool is_vertical);
    void gray_scale();
    void warping(int i);
    void restore();

    // Point selecting interaction
    void enable_selecting(bool flag);
    void select_points();
    void init_selections();

   private:
    // Store the original image data
    std::shared_ptr<Image> back_up_;
    // The selected point couples for image warping
    std::vector<ImVec2> start_points_, end_points_;

    ImVec2 start_, end_;
    bool flag_enable_selecting_points_ = false;
    bool draw_status_ = false;

   private:
    // A simple "fish-eye" warping function
    std::pair<int, int> fisheye_warping(int x, int y, int width, int height);
    std::vector<Pixel> gap_filling(std::vector<Pixel>& mapped_to, std::vector<Pixel>& unmapped_to);
};

}  // namespace USTC_CG