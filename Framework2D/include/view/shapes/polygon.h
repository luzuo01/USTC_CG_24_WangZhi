#pragma once

#include "shape.h"
#include <vector>
#include <utility>
#include "line.h"

namespace USTC_CG
{
class Polygon : public Shape
{
   public:
    Polygon() = default;
    Polygon(float s_x, float s_y);

    virtual ~Polygon() = default;

    // Overrides draw function to implement line-specific drawing logic
    void draw(const Config& config) const override;

    // Overrides Shape's update function to adjust the end point during
    // interaction
    void update(float x, float y) override;

   private:
   std::vector<std::pair<float, float>> pList;
};
}  // namespace USTC_CG
