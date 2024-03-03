#include <imgui.h>
#include "view/shapes/ellipse.h"

namespace USTC_CG
{
    void Ellipse::draw(const Config& config) const
    {
        ImDrawList* draw_list = ImGui::GetWindowDrawList();

        draw_list->AddEllipse(
            ImVec2(
            config.bias[0] + (start_point_x + end_point_x) / 2,
            config.bias[1] + (start_point_y + end_point_y) / 2),
            abs( (end_point_x - start_point_x) / 2),
            abs( (end_point_y - start_point_y) / 2),
            IM_COL32(
            config.line_color[0],
            config.line_color[1],
            config.line_color[2],
            config.line_color[3])
        );
    }

    void Ellipse::update(float x, float y)
    {
        end_point_x = x;
        end_point_y = y;
    }
}  