#include "comp_target_image.h"
#include <iostream>
#include <cmath>
#include <Eigen/SparseCholesky>
#include <Eigen/Dense>
namespace USTC_CG
{
using uchar = unsigned char;

CompTargetImage::CompTargetImage(
    const std::string& label,
    const std::string& filename)
    : ImageEditor(label, filename)
{
    if (data_){
        back_up_ = std::make_shared<Image>(*data_);
        region.resize(data_ -> width(), std::vector<bool>(data_ -> height(), false));
        boundary.resize(data_ -> width(), std::vector<bool>(data_ -> height(), false));
    }
}

void CompTargetImage::draw()
{
    // Draw the image
    ImageEditor::draw();
    // Invisible button for interactions
    ImGui::SetCursorScreenPos(position_);
    ImGui::InvisibleButton(
        label_.c_str(),
        ImVec2(
            static_cast<float>(image_width_),
            static_cast<float>(image_height_)),
        ImGuiButtonFlags_MouseButtonLeft);
    bool is_hovered_ = ImGui::IsItemHovered();
    // When the mouse is clicked or moving, we would adapt clone function to
    // copy the selected region to the target.
    ImGuiIO& io = ImGui::GetIO();
    if (is_hovered_ && ImGui::IsMouseClicked(ImGuiMouseButton_Left))
    {
        mouse_position_ =
            ImVec2(io.MousePos.x - position_.x, io.MousePos.y - position_.y);
        clone();
        edit_status_ = true;
    }
    if (edit_status_)
    {
        mouse_position_ =
            ImVec2(io.MousePos.x - position_.x, io.MousePos.y - position_.y);
        if (flag_realtime_updating)
            clone();
        if (!ImGui::IsMouseDown(ImGuiMouseButton_Left))
        {
            edit_status_ = false;
        }
    }
}

void CompTargetImage::set_source(std::shared_ptr<CompSourceImage> source)
{
    source_image_ = source;
}

void CompTargetImage::set_realtime(bool flag)
{
    flag_realtime_updating = flag;
}

void CompTargetImage::restore()
{
    *data_ = *back_up_;
    update();
}

void CompTargetImage::set_paste()
{
    clone_type_ = kPaste;
}

void CompTargetImage::set_seamless()
{
    clone_type_ = kSeamless;
}

void CompTargetImage::set_mixed()
{
    clone_type_ = kMixed;
}


std::vector<Pixel> CompTargetImage::neighbors(Pixel p) {
    int x = p.x;
    int y = p.y;
    std::vector<Pixel> neighbors;
    std::vector<std::pair<int, int>> index = {{0, 1} , {0, -1}, {1, 0}, {-1, 0}};
    for (auto pair : index)
    {
        if (x + pair.first >= 0 && x + pair.first < data_-> width()
        && y + pair.second >= 0 && y + pair.second < data_ -> height())
        {
            Pixel neighbor(x + pair.first, y + pair.second, 
            data_ -> get_pixel(x + pair.first, y + pair.second));
            neighbors.push_back(neighbor);
        }
    }
    return neighbors;
}


std::vector<Pixel> CompTargetImage::boundary_condition(std::shared_ptr<Image>& source_ptr) {
    std::vector<Pixel> result;
    this -> boundary.assign(data_ -> width(), std::vector<bool>(data_ -> height(), false));
    for (int i = 0; i < data_ -> width(); i++)
    {
        for (int j = 0; j < data_ -> height(); j++)
        {
            Pixel current_pixel (i, j , data_->get_pixel(i, j));
            std::vector<Pixel> nei = neighbors(current_pixel);
            // check if the current pixel is in the region
            // only when not in region do we proceed
            if (region[i][j] == false)
            {
                for (Pixel neighbor: nei)
                {
                    // if neighbor is in the region
                    // we add it to boundary
                    if (region[neighbor.x][neighbor.y])
                    {
                        result.push_back(current_pixel);
                        boundary[i][j] = true;
                        break;
                    }
                }
            }
        } 
    }
    return result;
}

std::vector<Pixel> CompTargetImage::gradient_condition(std::shared_ptr<Image>& source_ptr) {
    std::vector<Pixel> result;
    region_map.clear();
    this -> region.assign(data_ -> width(), std::vector<bool>(data_ -> height(), false));
    int count = 0;
    for (int i = 0; i < source_ptr->width(); ++i) {
        for (int j = 0; j < source_ptr->height(); ++j)
        {
            int tar_x =
                static_cast<int>(mouse_position_.x) + i -
                static_cast<int>(source_image_->get_position().x);
            int tar_y =
                static_cast<int>(mouse_position_.y) + j -
                static_cast<int>(source_image_->get_position().y);
            if (0 <= tar_x && tar_x < image_width_ && 0 <= tar_y &&
                tar_y < image_height_ && source_ptr->get_pixel(i, j)[0] > 0)
            {
                Pixel p(
                    tar_x,
                    tar_y,
                    source_image_ -> get_data() -> get_pixel(i, j));
                result.push_back(p);
                region[tar_x][tar_y] = true;
                region_map.insert({{tar_x, tar_y}, count});
                count++;
            }
        }
    }
    return result;
}

std::pair<int, int> CompTargetImage::inverse(int x, int y) {
   int tar_x = x - static_cast<int>(mouse_position_.x) + static_cast<int>(source_image_->get_position().x);
   int tar_y = y - static_cast<int>(mouse_position_.y) + static_cast<int>(source_image_->get_position().y); 
   return {tar_x, tar_y};
}

bool CompTargetImage::is_neighbor(Pixel p, Pixel q) {
    if (p.x >= 0 && p.x <= data_ -> width()
    && p.y >= 0 && p.y <= data_ -> height()
    && q.x >= 0 && q.x <= data_ -> width()
    && q.y>= 0 && q.y <= data_ -> height()
    && (pow(p.x - q.x, 2) + pow(p.y - q.y, 2) == 1)
    ) {
        return true;
    }
    return false;
}

VectorXf CompTargetImage::get_target_vector_at_channel
(std::vector<Pixel> boundary,  std::vector<Pixel> region, int i) {
    VectorXf vec(region.size());
    int count = 0;
    for (Pixel p : region)
    {
        float entry = 0;
        // this for takes care of the boundary part
        for (Pixel q : neighbors(p))
        {
            if (this->boundary[q.x][q.y])
            {
                entry += q.values[i];
            }
        }
        // this takes care of the gradient part.
        for (Pixel q : neighbors(p))
        {
            // get data of the preimage of this pixel
            // to comopute v_pq
            std::pair pre_p = inverse(p.x, p.y);
            std::pair pre_q = inverse(q.x, q.y);
            // if the inverse is invalid we set pre_q to the 
            // preimage of p, i.e. no contribution
            if (pre_q.first < 0 || pre_q.first >= source_image_ -> get_data() -> width()
            || pre_q.second < 0 || pre_q.second >= source_image_ -> get_data() -> height())
            {
                pre_q = pre_p;
            }
            float g_p = source_image_ -> get_data() -> get_pixel (pre_p.first, pre_p.second)[i];
            float g_q = source_image_ -> get_data() -> get_pixel (pre_q.first, pre_q.second)[i];
            entry += g_p - g_q;
        }
        vec[count] = entry;
        count++;
    }
    return vec;
}

VectorXf CompTargetImage::mixed_get_target_vector_at_channel
(std::vector<Pixel> boundary,  std::vector<Pixel> region, int i) {
    VectorXf vec(region.size());
    int count = 0;
    for (Pixel p : region)
    {
        float entry = 0;
        // this for takes care of the boundary part
        for (Pixel q : neighbors(p))
        {
            if (this->boundary[q.x][q.y])
            {
                entry += q.values[i];
            }
        }
        // this takes care of the gradient part.
        for (Pixel q : neighbors(p))
        {
            // get data of the preimage of this pixel
            // to comopute v_pq
            bool is_in_region = false;
            std::pair pre_p = inverse(p.x, p.y);
            std::pair pre_q = inverse(q.x, q.y);
            // if the inverse is invalid we set pre_q to the 
            // preimage of p, i.e. no contribution
            if (pre_q.first < 0 || pre_q.first >= source_image_ -> get_data() -> width()
            || pre_q.second < 0 || pre_q.second >= source_image_ -> get_data() -> height())
            {
                pre_q = pre_p;
                is_in_region = true;
            }
            float g_p = source_image_ -> get_data() -> get_pixel (pre_p.first, pre_p.second)[i];
            float g_q = source_image_ -> get_data() -> get_pixel (pre_q.first, pre_q.second)[i];
            float f_p_star = data_ -> get_pixel(p.x, p.y)[i];
            float f_q_star = data_ -> get_pixel(q.x, q.y)[i];
            if (std::abs(g_p - g_q) < std::abs(f_p_star - f_q_star)
            && !is_in_region)
            {
                entry += f_p_star - f_q_star;
            } else if (std::abs(g_p - g_q) < std::abs(f_p_star - f_q_star)
            && is_in_region)
            {
                entry += 0;
            } else {
                entry += g_p - g_q;
            }  
        }
        vec[count] = entry;
        count++;
    }
    return vec;
}

// the followings return the target vector and the sparse matrix
std::vector<VectorXf> CompTargetImage::get_target_vector(std::vector<Pixel> boundary,  std::vector<Pixel> region) {
    std::vector<VectorXf> targets;
    for (int i = 0; i < data_ -> channels(); i++)
    { 
        targets.push_back(get_target_vector_at_channel(boundary, region, i));
    }
    return targets;
}

std::vector<VectorXf> CompTargetImage::mixed_get_target_vector(std::vector<Pixel> boundary,  std::vector<Pixel> region) {
    std::vector<VectorXf> targets;
    for (int i = 0; i < data_ -> channels(); i++)
    { 
        targets.push_back(mixed_get_target_vector_at_channel(boundary, region, i));
    }
    return targets;
}

MatrixXf CompTargetImage::get_sparse_matrix(std::vector<Pixel> boundary,  std::vector<Pixel> region) {
    MatrixXf mat = MatrixXf::Zero(region.size(), region.size());
    int row_count = 0;
    for (Pixel p : region)
    {
        mat(row_count, row_count) = neighbors(p).size();
        std::vector<int> where_to_put_negative_one;
        for (Pixel q : neighbors(p))
        {
            // check if q is in the region, constant time
            if (this->region[q.x][q.y])
            {
                int index = region_map[{q.x, q.y}];
                where_to_put_negative_one.push_back(index);
            }
        }
        for (int here : where_to_put_negative_one)
        {
            mat(row_count, here) = -1;
        }
        row_count++;
    }
    return mat;
}

std::vector<Pixel> CompTargetImage::possion_solver(std::vector<VectorXf> targets, std::vector<Pixel> region) {
    std::vector<VectorXf> solution;
    std::vector<Pixel> result;

    for (int i = 0; i < targets.size(); i++)
    {
        VectorXf sol = solver.solve(targets[i]);
        solution.push_back(sol);
    }

    int count = 0;
    for (Pixel p : region)
    {
        std::vector<unsigned char> color;
        for (int i = 0; i < solution.size(); i++)
        {
            color.push_back(solution[i][count]);
        }
        Pixel new_p(p.x, p.y, color);
        result.push_back(new_p);
        count++;
    }
    return result;
}

void CompTargetImage::clone()
{
    // The implementation of different types of cloning
    // HW3_TODO: In this function, you should at least implement the "seamless"
    // cloning labeled by `clone_type_ ==kSeamless`.
    //
    // The realtime updating (update when the mouse is moving) is only available
    // when the checkboard is selected. It is required to improve the efficiency
    // of your seamless cloning to achieve realtime editing. (Use decomposition
    // of sparse matrix before solve the linear system)
    if (data_ == nullptr || source_image_ == nullptr ||
        source_image_->get_region() == nullptr)
        return;
    std::shared_ptr<Image> mask = source_image_->get_region();

    switch (clone_type_)
    {
        case USTC_CG::CompTargetImage::kDefault: break;
        case USTC_CG::CompTargetImage::kPaste:
        {
            restore();

            for (int i = 0; i < mask->width(); ++i)
            {
                for (int j = 0; j < mask->height(); ++j)
                {
                    int tar_x =
                        static_cast<int>(mouse_position_.x) + i -
                        static_cast<int>(source_image_->get_position().x);
                    int tar_y =
                        static_cast<int>(mouse_position_.y) + j -
                        static_cast<int>(source_image_->get_position().y);
                    if (0 <= tar_x && tar_x < image_width_ && 0 <= tar_y &&
                        tar_y < image_height_ && mask->get_pixel(i, j)[0] > 0)
                    {
                        data_->set_pixel(
                            tar_x,
                            tar_y,
                            source_image_->get_data()->get_pixel(i, j));
                    }
                }
            }
            break;
        }
        case USTC_CG::CompTargetImage::kSeamless:
        {
            restore();
            std::vector<Pixel> gradient = gradient_condition(mask);
            std::vector<Pixel> boundary = boundary_condition(mask);

            // the following two return the target vector and the sparse matrix
            std::vector<VectorXf> targets = get_target_vector(boundary, gradient);

            if (flag_realtime_updating == false ||
            (flag_realtime_updating == true && edit_status_ == false))
            {
                MatrixXf sparse_matrix = get_sparse_matrix(boundary, gradient);
                // Conversion to a sparse matrix
                sparseMatrix = sparse_matrix.sparseView();
                // Make the sparse matrix compressed
                sparseMatrix.makeCompressed();
                solver.compute(sparseMatrix);
            }
            std::vector<Pixel> solution = possion_solver(targets, gradient);


            int count = 0;
            for (int i = 0; i < mask->width(); ++i)
            {
                for (int j = 0; j < mask->height(); ++j)
                {
                    int tar_x =
                        static_cast<int>(mouse_position_.x) + i -
                        static_cast<int>(source_image_->get_position().x);
                    int tar_y =
                        static_cast<int>(mouse_position_.y) + j -
                        static_cast<int>(source_image_->get_position().y);
                    if (0 <= tar_x && tar_x < image_width_ && 0 <= tar_y &&
                        tar_y < image_height_ && mask->get_pixel(i, j)[0] > 0)
                    {
                        data_->set_pixel(
                            tar_x,
                            tar_y,
                            solution[count].values);
                        count++;
                    }
                }
            }

            break;
        }
        case USTC_CG::CompTargetImage::kMixed:
        {
            restore();
            std::vector<Pixel> gradient = gradient_condition(mask);
            std::vector<Pixel> boundary = boundary_condition(mask);

            // the following two return the target vector and the sparse matrix
            std::vector<VectorXf> targets = mixed_get_target_vector(boundary, gradient);

            if (flag_realtime_updating == false ||
            (flag_realtime_updating == true && edit_status_ == false))
            {
                MatrixXf sparse_matrix = get_sparse_matrix(boundary, gradient);
                // Conversion to a sparse matrix
                sparseMatrix = sparse_matrix.sparseView();
                // Make the sparse matrix compressed
                sparseMatrix.makeCompressed();
                solver.compute(sparseMatrix);
            }
            std::vector<Pixel> solution = possion_solver(targets, gradient);


            int count = 0;
            for (int i = 0; i < mask->width(); ++i)
            {
                for (int j = 0; j < mask->height(); ++j)
                {
                    int tar_x =
                        static_cast<int>(mouse_position_.x) + i -
                        static_cast<int>(source_image_->get_position().x);
                    int tar_y =
                        static_cast<int>(mouse_position_.y) + j -
                        static_cast<int>(source_image_->get_position().y);
                    if (0 <= tar_x && tar_x < image_width_ && 0 <= tar_y &&
                        tar_y < image_height_ && mask->get_pixel(i, j)[0] > 0)
                    {
                        data_->set_pixel(
                            tar_x,
                            tar_y,
                            solution[count].values);
                        count++;
                    }
                }
            }

            break;
        }
        default: break;
    }

    update();
}

}  // namespace USTC_CG