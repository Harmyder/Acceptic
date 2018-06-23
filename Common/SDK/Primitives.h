#pragma once
#include "SDK/Types.h"
#include "SDK/Utility.h"

namespace Common {
    namespace SDK {
        template <class T>
        Point<T, 2> Intersection(const Line<T>& a, const Line<T>& b, T maxDiff) {
            const bool aver = a.IsVertical(maxDiff);
            const bool bver = b.IsVertical(maxDiff);

            if (aver && bver) {
                return Point2<T>(kInfinity);
            }

            if (aver) {
                const T x = a.c / a.a;
                const T y = b.ComputeY(x);
                return { x,y };
            }

            if (bver) {
                const T x = b.c / b.a;
                const T y = a.ComputeY(x);
                return { x,y };
            }

            const T slopeA = a.a / a.b;
            const T slopeB = b.a / b.b;
            if (AlmostEqualRelativeAndAbs(slopeA, slopeB, maxDiff)) {
                return Point2<T>(kInfinity);
            }

            const T x = (a.c / a.b - b.c / b.b) / (a.a / a.b - b.a / b.b);
            assert(AlmostEqualRelativeAndAbs(a.ComputeY(x), b.ComputeY(x), maxDiff));
            const T y = a.ComputeY(x);
            return { x,y };
        }

        template <class T>
        bool IsPlaneCover(const HalfPlane<T>& a, const HalfPlane<T>& b, const HalfPlane<T>& c, T maxDiff) {
            const auto checkParallel = [](const Point2<T>& p, const HalfPlane<T>& hp1, const HalfPlane<T>& hp2) {
                if (p.IsInfinite()) {
                    const T d = Dot(hp1.GetNormal(), hp2.GetNormal());
                    const bool codirected = d > 0;
                    if (!codirected) {
                        return true;
                    }
                    return false;
                }
            };
            const auto gamma  = Intersection(a.hp, b.hp, maxDiff);
            if (checkParallel(gamma, a, b)) {
                return true;
            }
            const auto alpha = Intersection(b.hp, c.hp, maxDiff);
            if (checkParallel(alpha, b, c)) {
                return true;
            }
            const auto beta  = Intersection(c.hp, a.hp, maxDiff);
            if (checkParallel(beta, c, a)) {
                return true;
            }
            const bool gammaOk = c.IsInside(gamma, maxDiff);
            const bool betaOk  = b.IsInside(beta, maxDiff);
            const bool alphaOk = a.IsInside(alpha, maxDiff);
            const bool cover = gammaOk && betaOk && alphaOk;
            return cover;
        }
    }
}