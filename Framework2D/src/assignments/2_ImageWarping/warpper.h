#pragma once
#include <imgui.h>
#include <vector>

namespace USTC_CG
{

class Warpper
{
private:
public:
    virtual ~Warpper() = default;
    virtual std::pair<int, int> warp(int x, int y, int width, int height) = 0;

    virtual void initialize() = 0;
};


}