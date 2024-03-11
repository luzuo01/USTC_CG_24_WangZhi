#include "comp_warping.h"
#include "warpper.h"
#include <cmath>
#include "fisheye.h"
#include "IDW.h"
#include "RBF.h"
#include <annoylib.h>
#include <kissrandom.h>
using namespace Annoy;
namespace USTC_CG
{
using uchar = unsigned char;

CompWarping::CompWarping(const std::string& label, const std::string& filename)
    : ImageEditor(label, filename)
{
    if (data_)
        back_up_ = std::make_shared<Image>(*data_);
}

void CompWarping::draw()
{
    // Draw the image
    ImageEditor::draw();
    // Draw the canvas
    if (flag_enable_selecting_points_)
        select_points();
}

void CompWarping::invert()
{
    for (int i = 0; i < data_->width(); ++i)
    {
        for (int j = 0; j < data_->height(); ++j)
        {
            const auto color = data_->get_pixel(i, j);
            data_->set_pixel(
                i,
                j,
                { static_cast<uchar>(255 - color[0]),
                  static_cast<uchar>(255 - color[1]),
                  static_cast<uchar>(255 - color[2]) });
        }
    }
    // After change the image, we should reload the image data to the renderer
    update();
}
void CompWarping::mirror(bool is_horizontal, bool is_vertical)
{
    Image image_tmp(*data_);
    int width = data_->width();
    int height = data_->height();

    if (is_horizontal)
    {
        if (is_vertical)
        {
            for (int i = 0; i < width; ++i)
            {
                for (int j = 0; j < height; ++j)
                {
                    data_->set_pixel(
                        i,
                        j,
                        image_tmp.get_pixel(width - 1 - i, height - 1 - j));
                }
            }
        }
        else
        {
            for (int i = 0; i < width; ++i)
            {
                for (int j = 0; j < height; ++j)
                {
                    data_->set_pixel(
                        i, j, image_tmp.get_pixel(width - 1 - i, j));
                }
            }
        }
    }
    else
    {
        if (is_vertical)
        {
            for (int i = 0; i < width; ++i)
            {
                for (int j = 0; j < height; ++j)
                {
                    data_->set_pixel(
                        i, j, image_tmp.get_pixel(i, height - 1 - j));
                }
            }
        }
    }

    // After change the image, we should reload the image data to the renderer
    update();
}
void CompWarping::gray_scale()
{
    for (int i = 0; i < data_->width(); ++i)
    {
        for (int j = 0; j < data_->height(); ++j)
        {
            const auto color = data_->get_pixel(i, j);
            uchar gray_value = (color[0] + color[1] + color[2]) / 3;
            data_->set_pixel(i, j, { gray_value, gray_value, gray_value });
        }
    }
    // After change the image, we should reload the image data to the renderer
    update();
}


std::vector<Pixel> CompWarping::gap_filling(std::vector<Pixel>& mapped_to,
 std::vector<Pixel>& unmapped_to) {
    std::vector<Pixel> result;
    // Build Annoy index
    int dimension = 2;
    AnnoyIndex<int, float, Euclidean, Kiss32Random, AnnoyIndexSingleThreadedBuildPolicy> index(dimension);
    for (int i = 0; i < mapped_to.size(); ++i) {
        Pixel pixel = mapped_to[i];
        vector<float> coord = {(float) pixel.x, (float) pixel.y};
        index.add_item(i, coord.data());
    }
    index.build(12);  // Build the index with some trees

    int count = mapped_to.size();
    for (auto& uncolored_pixel : unmapped_to)
    {
        // Find k nearest neighbors for each uncolored pixel and average their values
        const int k = 70; 
        vector<float> coord = {(float) uncolored_pixel.x, (float) uncolored_pixel.y};
        index.add_item(count, coord.data());
        count++;
        vector<int> nearest_indices;
        vector<float> distances;
        index.get_nns_by_vector(coord.data(), k, -1, &nearest_indices, &distances);  // Get nearest neighbors

        // get the total distance
        float total_distance = 0;
        for (int i = 0; i < k; ++i) {
            total_distance += exp(distances[i]);
        }
        int colored_count = 0;
        int uncolored_count = 0;
        for (int i = 0; i < k; i++)
        {
            if (nearest_indices[i] < mapped_to.size())
            {
                // this means that it's a colored point 
                colored_count ++;
            } else {
                uncolored_count ++;
            }
        }
        //initialize average value
        vector<unsigned char> average_value(uncolored_pixel.values.size(), 0); 
       //only interpolate when it's surronded by colored points
        if (uncolored_count <= 1) {
            // now we average the colors of these k neighbors by distance squared
            for (int i = 0; i < k; ++i) {
                if (nearest_indices[i] < mapped_to.size())
                {
                    for (size_t j = 0; j < mapped_to[nearest_indices[i]].values.size(); ++j) {
                        if (nearest_indices[i] < mapped_to.size())
                        {
                            average_value[j] += mapped_to[nearest_indices[i]].values[j] * 
                            exp(distances[i]) / total_distance;
                        } else {
                            average_value[j] += 
                            unmapped_to[nearest_indices[i] - mapped_to.size()].values[j] *
                            exp(distances[i]) / total_distance;
                        }
                    }
                } 
            }
        }
        Pixel pixel(uncolored_pixel.x, uncolored_pixel.y, average_value);
        result.push_back(pixel);
    }
    return result;
}

void CompWarping::warping(int i)
{
    Warpper* our_warpper;
    if (i == 1)
    {
        our_warpper = new Fisheye();
    } else if (i == 2 && start_points_.size() != 0)
    {
        our_warpper = new IDW(start_points_, end_points_);
        //the next line should optimizes the D's in IDW, but I haven't fixed the bugs yet
        //our_warpper->initialize();

    } else if ((i == 3 || i == 4) && start_points_.size() != 0)
    {
        our_warpper = new RBF(start_points_, end_points_);
        our_warpper->initialize();
    } else {
        // to prevent the program from crashing if the user hasn't seleceted points.
        our_warpper = new RBF();
    }
    // Create a new image to store the result
    Image warped_image(*data_);
    // we'll keep track of vectors of mapped_to and unmapped_to pixels
    std::vector<Pixel> mapped_to, unmapped_to;
    std::vector<std::vector<bool>> range(data_->width(), 
    std::vector<bool> (data_->height(), false));
    // Initialize the color of result image
    for (int y = 0; y < data_->height(); ++y)
    {
        for (int x = 0; x < data_->width(); ++x)
        {
            warped_image.set_pixel(x, y, { 0, 0, 0 });
        }
    }

    // For each (x, y) from the input image, the "fish-eye" warping transfer it
    // to (x', y') in the new image:
    // Note: For this transformation ("fish-eye" warping), one can also
    // calculate the inverse (x', y') -> (x, y) to fill in the "gaps".
    for (int y = 0; y < data_->height(); ++y)
    {
        for (int x = 0; x < data_->width(); ++x)
        {
            // Apply warping function to (x, y), and we can get (x', y')
            auto [new_x, new_y] =
                our_warpper -> warp(x, y, data_->width(), data_->height());
            // Copy the color from the original image to the result image
            if (new_x >= 0 && new_x < data_->width() && new_y >= 0 &&
                new_y < data_->height())
            {
                std::vector<unsigned char> pixel = data_->get_pixel(x, y);
                warped_image.set_pixel(new_x, new_y, pixel);
                // add the pixel to mapped_to and mark the coordinates in range
                Pixel p = {new_x, new_y, pixel}; 
                if (range[new_x][new_y] == false)
                {
                mapped_to.push_back(p);
                range[new_x][new_y] = true;
                }
            }
        }
    }

    // set unmapped_to
     for (int y = 0; y < data_->height(); ++y)
    {
        for (int x = 0; x < data_->width(); ++x)
        {
            if (range[x][y] == false)
            {
                std::vector<unsigned char> color = data_->get_pixel(x, y);
                Pixel p = {x, y, color};
                unmapped_to.push_back(p);
            }
        }
    }


    // with gap filling, RBF
    if (i == 4)
    {
        //gap filling 
        std::vector<Pixel> gaps = gap_filling(mapped_to, unmapped_to);
        for (auto& pt : gaps) {
            std::vector<unsigned char> pixel = warped_image.get_pixel(pt.x, pt.y);
            warped_image.set_pixel(pt.x, pt.y, pt.values);
        }
    }
    *data_ = std::move(warped_image);
    update();
}



void CompWarping::restore()
{
    *data_ = *back_up_;
    update();
}
void CompWarping::enable_selecting(bool flag)
{
    flag_enable_selecting_points_ = flag;
}
void CompWarping::select_points()
{
    /// Invisible button over the canvas to capture mouse interactions.
    ImGui::SetCursorScreenPos(position_);
    ImGui::InvisibleButton(
        label_.c_str(),
        ImVec2(
            static_cast<float>(image_width_),
            static_cast<float>(image_height_)),
        ImGuiButtonFlags_MouseButtonLeft);
    // Record the current status of the invisible button
    bool is_hovered_ = ImGui::IsItemHovered();
    // Selections
    ImGuiIO& io = ImGui::GetIO();
    if (is_hovered_ && ImGui::IsMouseClicked(ImGuiMouseButton_Left))
    {
        draw_status_ = true;
        start_ = end_ =
            ImVec2(io.MousePos.x - position_.x, io.MousePos.y - position_.y);
    }
    if (draw_status_)
    {
        end_ = ImVec2(io.MousePos.x - position_.x, io.MousePos.y - position_.y);
        if (!ImGui::IsMouseDown(ImGuiMouseButton_Left))
        {
            start_points_.push_back(start_);
            end_points_.push_back(end_);
            draw_status_ = false;
        }
    }
    // Visualization
    auto draw_list = ImGui::GetWindowDrawList();
    for (size_t i = 0; i < start_points_.size(); ++i)
    {
        ImVec2 s(
            start_points_[i].x + position_.x, start_points_[i].y + position_.y);
        ImVec2 e(
            end_points_[i].x + position_.x, end_points_[i].y + position_.y);
        draw_list->AddLine(s, e, IM_COL32(255, 0, 0, 255), 2.0f);
        draw_list->AddCircleFilled(s, 4.0f, IM_COL32(0, 0, 255, 255));
        draw_list->AddCircleFilled(e, 4.0f, IM_COL32(0, 255, 0, 255));
    }
    if (draw_status_)
    {
        ImVec2 s(start_.x + position_.x, start_.y + position_.y);
        ImVec2 e(end_.x + position_.x, end_.y + position_.y);
        draw_list->AddLine(s, e, IM_COL32(255, 0, 0, 255), 2.0f);
        draw_list->AddCircleFilled(s, 4.0f, IM_COL32(0, 0, 255, 255));
    }
}
void CompWarping::init_selections()
{
    start_points_.clear();
    end_points_.clear();
}

std::pair<int, int>
CompWarping::fisheye_warping(int x, int y, int width, int height)
{
    float center_x = width / 2.0f;
    float center_y = height / 2.0f;
    float dx = x - center_x;
    float dy = y - center_y;
    float distance = std::sqrt(dx * dx + dy * dy);

    // Simple non-linear transformation r -> r' = f(r)
    float new_distance = std::sqrt(distance) * 10;

    if (distance == 0)
    {
        return { static_cast<int>(center_x), static_cast<int>(center_y) };
    }
    // (x', y')
    float ratio = new_distance / distance;
    int new_x = static_cast<int>(center_x + dx * ratio);
    int new_y = static_cast<int>(center_y + dy * ratio);

    return { new_x, new_y };
}
}  // namespace USTC_CG