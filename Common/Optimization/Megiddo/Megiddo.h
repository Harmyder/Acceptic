#pragma once
#include "Optimization/lp.h"

namespace Common {
    namespace Optimization {
        LpResult2Dim solve(const LpTask2Dim& task);
    }
}
