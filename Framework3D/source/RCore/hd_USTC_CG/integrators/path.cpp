#include "path.h"

#include <random>

#include "surfaceInteraction.h"
USTC_CG_NAMESPACE_OPEN_SCOPE
using namespace pxr;

VtValue PathIntegrator::Li(const GfRay& ray, std::default_random_engine& random)
{
    std::uniform_real_distribution<float> uniform_dist(
        0.0f, 1.0f - std::numeric_limits<float>::epsilon());
    std::function<float()> uniform_float = std::bind(uniform_dist, random);

    auto color = EstimateOutGoingRadiance(ray, uniform_float, 0);

    return VtValue(GfVec3f(color[0], color[1], color[2]));
}

GfVec3f PathIntegrator::EstimateOutGoingRadiance(
    const GfRay& ray,
    const std::function<float()>& uniform_float,
    int recursion_depth)
{
    if (recursion_depth >= 12) {
        return {};
    }

    SurfaceInteraction si;
    if (!Intersect(ray, si)) {
        if (recursion_depth == 0) {
            return IntersectDomeLight(ray);
        }

        return GfVec3f{ 0, 0, 0 };
    }

    // This can be customized : Do we want to see the lights? (Other than dome lights?)
    if (recursion_depth == 0) {
    }

    // Flip the normal if opposite
    if (GfDot(si.shadingNormal, ray.GetDirection()) > 0) {
        si.flipNormal();
        si.PrepareTransforms();
    }

    GfVec3f color{ 0 };
    GfVec3f directLight = EstimateDirectLight(si, uniform_float);

    // HW7_TODO: Estimate global lighting here.
    GfVec3f globalLight = GfVec3f{0.f};

    float prr = 0.8;
    float ksi = uniform_float();
    if (ksi > prr) {
    } else {
        GfVec3d dir = ray.GetDirection().GetNormalized();
        GfVec3d norm = si.geometricNormal.GetNormalized();
        float cosVal = dir * norm;
        if (Intersect(ray, si)) {
            GfVec3f newRay;
            float pdf;
            si.Sample(newRay, pdf, uniform_float);
            GfVec3f starting_point = si.position + 0.00001f * si.geometricNormal;
            GfVec3f end_point = starting_point + newRay;
            GfRay nextRay(starting_point, end_point);
            globalLight = EstimateOutGoingRadiance(nextRay, uniform_float, recursion_depth + 1) / prr;
        }
    }
    color = directLight + globalLight;

    return color;
}

USTC_CG_NAMESPACE_CLOSE_SCOPE
