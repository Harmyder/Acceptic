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
            const auto isParallelCover = [](const Point2<T>& p, const HalfPlane<T>& hp1, const HalfPlane<T>& hp2) {
                if (p.IsInfinite()) {
                    const T d = Dot(hp1.GetNormal(), hp2.GetNormal());
                    const bool codirected = d > 0;
                    if (!codirected) {
                        return true;
                    }
                }
                return false;
            };
            const auto gamma  = Intersection(a.boundary, b.boundary, maxDiff);
            if (isParallelCover(gamma, a, b)) {
                return true;
            }
            const auto alpha = Intersection(b.boundary, c.boundary, maxDiff);
            if (isParallelCover(alpha, b, c)) {
                return true;
            }
            const auto beta  = Intersection(c.boundary, a.boundary, maxDiff);
            if (isParallelCover(beta, c, a)) {
                return true;
            }
            const bool gammaOk = c.IsInside(gamma, maxDiff);
            const bool betaOk  = b.IsInside(beta, maxDiff);
            const bool alphaOk = a.IsInside(alpha, maxDiff);
            const bool cover = gammaOk && betaOk && alphaOk;
            return cover;
        }

        template <class T>
        T Infinitize(T v, T infty) {
            if (v > 0) return infty;
            if (v < 0) return -infty;
            return 0;
        }

        template <class T>
        Point2<T> ClampOnLine(Point2<T> p, const Line<T>& line, std::pair<T, T> xb, std::pair<T, T> yb) {
            auto clampOnAxis = [&p, &line](int indA, int indB, const std::pair<T, T>& boundsA, const std::pair<T, T>& boundsB) -> bool {
                const T a = std::clamp(p.data[indA], boundsA.first, boundsA.second);
                if (p.data[indA] != a) {
                    const T b = (line.*Line<T>::Compute[indB])(a);
                    if (boundsB.first <= b && b <= boundsB.second) {
                        p.data[indA] = a;
                        p.data[indB] = b;
                        return true;
                    }
                }
                return false;
            };
            if (!clampOnAxis(0, 1, xb, yb)) {
                clampOnAxis(1, 0, yb, xb);
            }
            return p;
        }
    }
}