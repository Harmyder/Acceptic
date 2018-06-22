#pragma once

#include <vector>

namespace Common {
    namespace Optimization {
        struct LpTask2Dim {
            std::array<float, 2> c;
            std::vector<float> b;
            std::vector<std::array<float, 2>> A;
        };

        using LpResult2Dim = std::array<float, 2>;
    }
}
