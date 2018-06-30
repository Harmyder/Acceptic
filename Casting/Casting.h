#pragma once

#include <vector>
#include "Common/SDK/Types.h"

namespace Common {
    namespace Dcel {
        template <class T>
        class Mesh;
    }
}


extern const double kMaxDiff;
extern const double kMaxDiffSq;

namespace Casting
{
    std::vector<Common::SDK::Plane<double>> CombineFaces(const Common::Dcel::Mesh<int>& m, std::vector<Common::SDK::Point3<double>> verticesObj);

    const std::array<Common::SDK::Point3<double>, 4> hemisphereCoverForSphere = {
        Common::SDK::Normalize(Common::SDK::Point3<double>(1.,  .2,  0.), kMaxDiff),
        Common::SDK::Normalize(Common::SDK::Point3<double>(-1.,  .2,  0.), kMaxDiff),
        // small z-axis aligned cleft left uncovered toward negative y-axis 
        Common::SDK::Normalize(Common::SDK::Point3<double>(0., -.2,  1.), kMaxDiff),
        Common::SDK::Normalize(Common::SDK::Point3<double>(0., -.2, -1.), kMaxDiff)
    };
    const double kBound = 1 << 30;
    const int kCandidatesPerHemisphere = 3;
    const int kCoverHemispheresCount = 4;
    const Common::SDK::Point3<double> upperHemisphereDirection(0., 0., 1.);

    // Next two function are reverse of each over
    Common::SDK::HalfPlane<double> ProjectHemisphereOnZUnitPlane(Common::SDK::Point3<double> hemisphereDirection);
    Common::SDK::Point3<double> ComputeHemisphereDirFromItsZUnitPlaneBoundary(Common::SDK::HalfPlane<double> boundary);

    std::array<int, kCandidatesPerHemisphere> FindCover(const std::vector<Common::SDK::HalfPlane<double>>& matrix);
}
