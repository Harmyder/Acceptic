#pragma once

#include <vector>
#include "Common/SDK/Types.h"

namespace Common {
    namespace Dcel {
        template <class T>
        class Mesh;
    }
}

std::vector<Common::SDK::Plane<double>> CombineFaces(const Common::Dcel::Mesh<int>& m, std::vector<Common::SDK::Point3<double>> verticesObj);

std::array<Common::SDK::Plane<double>, 3> FindCoverForHemisphere(
    const std::vector<Common::SDK::Plane<double>>& faces,
    const Common::SDK::Point3<double>& hemi
);
