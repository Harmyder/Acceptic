#pragma once

#include <array>
#include <optional>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <iostream>
#include "DebugPrint.h"

#include "SDK/Primitives.h"
#include "SDK/Types.h"

namespace Common {
    namespace Optimization {
        namespace details {
            template <class T>
            bool IsParallelFeasible(const SDK::HalfPlane<T>& a, const SDK::HalfPlane<T>& b, T maxDiff) {
                SDK::Line<T> perp(-a.boundary.b, a.boundary.a, 0.);
                const auto intCurr = Intersection(perp, a.boundary, maxDiff);
                const auto intPrev = Intersection(perp, b.boundary, maxDiff);
                const auto inTheMiddle = (intCurr + intPrev) / 2;
                if (!a.IsInside(inTheMiddle, maxDiff) && !b.IsInside(inTheMiddle, maxDiff)) {
                    return false;
                }
                return true;
            }

            template <class T>
            T Infinitize(T v, T infty) {
                if (v > 0) return infty;
                if (v < 0) return -infty;
                return 0;
            }
        }

        template <class T>
        std::optional<SDK::Point2<T>> solveMax(
            const SDK::Point2<T>& c,
            const std::vector<SDK::HalfPlane<T>>& matrix, // Right side incorporated
            const std::pair<T, T>& xBounds, // Natural (domain specific bounds) on x
            const std::pair<T, T>& yBounds, // Natural (domain specific bounds) on x
            const T maxDiff
        ) {
            assert(xBounds.second > xBounds.first && yBounds.second > xBounds.first);
            bool leftToRight = c.x > 0;
            bool bottomToTop = c.y > 0;
            const Point2<T> minCorner( leftToRight ? xBounds.first : xBounds.second,  bottomToTop ? yBounds.first : yBounds.second);
            const Point2<T> maxCorner(!leftToRight ? xBounds.first : xBounds.second, !bottomToTop ? yBounds.first : yBounds.second);
            auto currP = minCorner;
            T currV = Dot(c, currP);
            const auto minCornerX = HalfPlane<T>(Line<T>(1., 0., minCorner.x), !leftToRight);
            const auto minCornerY = HalfPlane<T>(Line<T>(0., 1., minCorner.y), !bottomToTop);
            const auto maxCornerX = HalfPlane<T>(Line<T>(1., 0., maxCorner.x),  leftToRight);
            const auto maxCornerY = HalfPlane<T>(Line<T>(0., 1., maxCorner.y),  bottomToTop);

            const T infty = 1.0001 * std::max(std::max(std::abs(xBounds.first), std::abs(xBounds.first)), std::max(std::abs(yBounds.first), std::abs(yBounds.first)));

            // Protect from "bad" inputs
            const auto indices = [len = matrix.size()]() {
                std::vector<int> v(len);
                std::iota(v.begin(), v.end(), 0);
                //std::random_device rd;
                //std::minstd_rand g(rd());
                //std::shuffle(v.begin(), v.end(), g);
                return v;
            }();

            auto matrixAndBoundaries = [&minCornerX, &minCornerY, &maxCornerX, &maxCornerY, &matrix, &indices](int i) -> const HalfPlane<T>& {
                switch (i)
                {
                case 0: return minCornerX;
                case 1: return minCornerY;
                case 2: return maxCornerX;
                case 3: return maxCornerY;
                default: return matrix[indices[i - 4]];
                }
            };
            for (int i = 4; i < (int)matrix.size() + 4; ++i) {
                const HalfPlane<T>& hp = matrixAndBoundaries(i);
                if (!hp.IsInside(currP, maxDiff)) {
                    // Rotate normal 90 degrees CCW to get a separator.
                    const Point2<T> separator(SDK::Normalize(Point2<T>(-hp.GetNormal().y, hp.GetNormal().x)));
                    SDK::Point2<T> bestAgainst(details::Infinitize(separator.x, infty), details::Infinitize(separator.y, infty));
                    SDK::Point2<T> bestAlong(-bestAgainst);
                    for (int j = 0; j < i; ++j) {
                        const HalfPlane<T>& prev = matrixAndBoundaries(j);
                        const T d = Dot(prev.GetNormal(), separator);
                        const bool isParallel = AlmostEqualToZero(d, maxDiff);
                        if (isParallel) {
                            if (!details::IsParallelFeasible(hp, prev, maxDiff)) {
                                return nullopt;
                            }
                            continue;
                        }
                        const Point2<T> p = Intersection(hp.boundary, prev.boundary, maxDiff);
                        bool isAlong = d > 0;
                        if (isAlong) {
                            const Point2<T> increment = p - bestAlong;
                            if (Dot(increment, separator) > 0) {
                                bestAlong = p;
                            }                            
                        } else {
                            const Point2<T> increment = p - bestAgainst;
                            if (Dot(increment, separator) < 0) {
                                bestAgainst = p;
                            }
                        }
                    }
                    const bool bestAlongInf = bestAlong.IsInfinite();
                    const bool bestAgainstInf = bestAgainst.IsInfinite();
                    assert(!bestAlongInf || !bestAgainstInf);
                    SDK::Point2<T> candidateP(kNan);
                    T candidateV;
                    if (!bestAlongInf && !bestAgainstInf) {
                        const bool isValid = Dot(bestAgainst - bestAlong, separator) > 0;
                        const T againstV = Dot(bestAgainst, c);
                        const T alongV   = Dot(bestAlong, c);
                        if (alongV > againstV) {
                            candidateV = alongV;
                            candidateP = bestAlong;
                        }
                        else {
                            candidateV = againstV;
                            candidateP = bestAgainst;
                        }
                    } else if (bestAlong.IsInfinite()) {
                        candidateV = Dot(c, bestAgainst);
                        candidateP = bestAgainst;
                    } else {
                        assert(bestAgainst.IsInfinite());
                        candidateV = Dot(c, bestAlong);
                        candidateP = bestAlong;
                    }
                    if (candidateV > currV) {
                        currV = candidateV;
                        currP = candidateP;
                        DebugPrintf("%f, %f\n", currP.x, currP.y);
                    }
                }
            }

            return std::make_optional(currP);
        }

        //extern template std::optional<SDK::Point2<double>> solveMax(
        //    const SDK::Point2<double>& c,
        //    const std::vector<SDK::HalfPlane<double>>& matrix,
        //    const std::pair<double, double>& xBounds,
        //    const std::pair<double, double>& yBounds,
        //    const double maxDiff
        //);
    }
}
