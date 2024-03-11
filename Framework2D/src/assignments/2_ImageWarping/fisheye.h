#pragma once
#include "warpper.h"
#include <imgui.h>

namespace USTC_CG
{
class Fisheye : public Warpper
{
public:
    Fisheye() = default;
    ~Fisheye() = default;
    std::pair<int, int> warp(int x, int y, int width, int height) override;
    void initialize() override;
};
}