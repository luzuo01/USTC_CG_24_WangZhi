#include "view/shapes/polygon.h"

#include <imgui.h>

namespace USTC_CG
{
    Polygon::Polygon(float x, float y)
    {
        pList.push_back(std::make_pair(x, y));
    }
    // Draw the line using ImGui
    void Polygon::draw(const Config& config) const
    {
        ImDrawList* draw_list = ImGui::GetWindowDrawList();
        for (auto p = pList.begin(); std::next(p) != pList.end(); ++p)
        {
            Line current = Line(p -> first, p -> second, std::next(p) -> first, std::next(p) -> second);
            current.draw(config);
        }
    }

    void Polygon::update(float x, float y)
    {
         pList.push_back(std::make_pair(x, y));
    }

}  // namespace USTC_CG