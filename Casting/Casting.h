#pragma once

#include <vector>
#include "Common/SDK/Types.h"

namespace Common {
    namespace Dcel {
        template <class T>
        class Mesh;
    }
}

namespace Casting
{
    std::vector<Common::SDK::Plane<double>> CombineFaces(const Common::Dcel::Mesh<int>& m, std::vector<Common::SDK::Point3<double>> verticesObj);

    const std::array<Common::SDK::Point3<double>, 4> hemisphereCoverForSphere = {
        Common::SDK::Normalize(Common::SDK::Point3<double>(1.,  .2,  0.)),
        Common::SDK::Normalize(Common::SDK::Point3<double>(-1.,  .2,  0.)),
        // small z-axis aligned cleft left uncovered toward negative y-axis 
        Common::SDK::Normalize(Common::SDK::Point3<double>(0., -.2,  1.)),
        Common::SDK::Normalize(Common::SDK::Point3<double>(0., -.2, -1.))
    };
    const double kBound = 1 << 30;
    const int kCandidatesPerHemisphere = 3;
    const int kCoverHemispheresCount = 4;
    const Common::SDK::Point3<double> upperHemisphereDirection(0., 0., 1.);

    Common::SDK::HalfPlane<double> ProjectIntersectionOnUpperHemispherePlane(Common::SDK::Point3<double> hemisphereDirection);

    std::array<int, kCandidatesPerHemisphere> FindCover(const std::vector<Common::SDK::HalfPlane<double>>& matrix);
}
