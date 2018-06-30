#pragma once
#include <limits>

namespace Common {
    namespace SDK {
#define ASSERT_NORMAL(a, maxDiff, maxDiffSq) assert(SDK::AlmostEqualRelativeAndAbs(a.LenSq(), Normalize(a, maxDiff).LenSq(), maxDiffSq) && "Vector " #a " is not normalized.");

        enum class IsOverflow {
            Yes, No
        };

        template <class T>
        struct Overflow
        {
            static IsOverflow Sum(T a, T b) {
                return a + b == std::numeric_limits<float>::infinity() ? IsOverflow::Yes : IsOverflow::No;
            }
            static IsOverflow Mul(T a, T b) {
                return a * b == std::numeric_limits<float>::infinity() ? IsOverflow::Yes : IsOverflow::No;
            }
            static IsOverflow Div(T a, T b) {
                return a * b == std::numeric_limits<float>::infinity() ? IsOverflow::Yes : IsOverflow::No;
            }
        };

        template <class T>
        std::enable_if_t<std::is_floating_point_v<T>, bool> AlmostEqualToZero(T A, T maxDiff) {
            if (std::fabs(A) <= maxDiff)
                return true;
            return false;
        }

        // https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
        template <class T>
        std::enable_if_t<std::is_floating_point_v<T>, bool> AlmostEqualRelativeAndAbs(T A, T B, T maxDiff, T maxRelDiff = std::numeric_limits<T>::epsilon()) {
            // Check if the numbers are really close -- needed when comparing numbers near zero.
            T diff = std::fabs(A - B);
            if (AlmostEqualToZero(diff, maxDiff))
                return true;

            A = std::fabs(A);
            B = std::fabs(B);
            T largest = (B > A) ? B : A;

            if (diff <= largest * maxRelDiff)
                return true;
            return false;
        }
    }
}