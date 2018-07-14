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
    // I assume that the smallest possible angle between facets normals is
    const double kSmallestAngle = 1e-5;

    const double kCosSmallestAngle     = std::cos(Casting::kSmallestAngle);
    const double kCosSmallestHalfAngle = std::cos(Casting::kSmallestAngle / 2.);
    const double kCosSmallestAngleSq     = kCosSmallestAngle * kCosSmallestAngle;
    const double kCosSmallestHalfAngleSq = kCosSmallestHalfAngle * kCosSmallestHalfAngle;

    const double kSinSmallestAngle = std::sin(Casting::kSmallestAngle);
    const double kSinSmallestHalfAngle = std::sin(Casting::kSmallestAngle / 2.);
    const double kSinSmallestAngleSq = kSinSmallestAngle * kSinSmallestAngle;
    const double kSinSmallestHalfAngleSq = kSinSmallestHalfAngle * kSinSmallestHalfAngle;

    const double kBound = 1. / std::tan(kSmallestAngle / 2) + 1.;

    const int kCandidatesPerHemisphere = 4;
    const int kCoverHemispheresCount = 4;
    const Common::SDK::Point3<double> upperHemisphereDirection(0., 0., 1.);

    // Next two function are reverse of each over
    Common::SDK::HalfPlane<double> ProjectHemisphereOnZUnitPlane(Common::SDK::Point3<double> hemisphereDirection);
    Common::SDK::Point3<double> ComputeHemisphereDirFromItsZUnitPlaneBoundary(Common::SDK::HalfPlane<double> boundary);

    std::vector<int> FindMaxSubcover(const std::vector<Common::SDK::HalfPlane<double>>& matrix);

    bool IsTheSameFaceNormal(const Common::SDK::Point3<double>& n1, const Common::SDK::Point3<double>& n2);

    Common::SDK::Matrix<double, 3> ProjectFacesOnUpperHemisphere(
        const Common::SDK::Point3<double>& actualNormal,
        const std::vector<Common::SDK::Plane<double>>& bigFaces,
        const int excludeIndex,
        std::vector<Common::SDK::HalfPlane<double>>& matrix,
        std::vector<int>& indices
    );

    std::vector<int> FindS2Coverage(const std::vector<Common::SDK::Plane<double>>& bigFaces);
    std::vector<std::pair<int, Common::SDK::Point3<double>>> CheckCandidates(const std::vector<Common::SDK::Plane<double>>& bigFaces, const std::vector<int>& candidates);
}
