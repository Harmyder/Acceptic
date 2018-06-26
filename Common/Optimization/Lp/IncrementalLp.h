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

            template <class T>
            SDK::Point2<T> Clamp(SDK::Point2<T> p, T min, T max) {
                p.x = std::clamp(p.x, min, max);
                p.y = std::clamp(p.y, min, max);
                return p;
            }
        }

        template <class T>
        struct LpTask2dRes {
            SDK::Point2<T> point;
            std::pair<int, int> formingHalfPlanes;
        };

        //template <class T>
        //struct LpTask2dRes {
        //    enum { kNoIndex = -1 };

        //    const int kFirstTight = 0;
        //    const int kSecondTight = 1;

        //    const int kFirstInfeas = 2;
        //    const int kSecondInfeas = 1;
        //    const int kThirdInfeas = 0;

        //    LpTask2dRes CreateTightPlanes(int i, int j) {
        //        get<kFirstTight>(res._formingHalfPlanes) = i;
        //        get<kSecondTight>(res._formingHalfPlanes) = i;
        //        _point = p;
        //    }

        //    LpTask2dRes CreateParallelInfeasibility(int i, int j) {
        //        get<kFirstInfeas>(res._formingHalfPlanes) = i;
        //        get<kSecondInfeas>(res._formingHalfPlanes) = j;
        //        _point = p;
        //    }

        //    LpTask2dRes CreateTriangularlInfeasibility(Point2<T> p, int i, int j, int k) {
        //        LpTask2dRes res;
        //        get<kFirstInfeas>(res._formingHalfPlanes) = i;
        //        get<kSecondInfeas>(res._formingHalfPlanes) = j;
        //        get<kThirdInfeas>(res._formingHalfPlanes) = k;
        //        _point = p;
        //    }

        //    bool IsFeasible() { return get<0>(_formingHalfPlanes) == kNoIndex && get<2>(_formingHalfPlanes) != kNoIndex; }
        //    auto GetPoint() const { return point; }

        //    enum class Infeasibility {
        //        Parallel,
        //        Triangular,
        //        None
        //    };
        //    Infeasibility GetInfeasibilityType
        //private:
        //    SDK::Point2<T> _point;
        //    std::tuple<int, int, int> _formingHalfPlanes;
        //};

        template <class T>
        std::optional<LpTask2dRes<T>> solveMax(
            const SDK::Point2<T>& c,
            const std::vector<SDK::HalfPlane<T>>& matrix, // Right side incorporated
            const std::pair<T, T>& xBounds, // Natural (domain specific bounds) on x
            const std::pair<T, T>& yBounds, // Natural (domain specific bounds) on x
            const T maxDiff
        ) {
            assert(xBounds.second > xBounds.first && yBounds.second > xBounds.first);
            bool leftToRight = c.x > 0;
            bool bottomToTop = c.y > 0;
            const SDK::Point2<T> minCorner( leftToRight ? xBounds.first : xBounds.second,  bottomToTop ? yBounds.first : yBounds.second);
            const SDK::Point2<T> maxCorner(!leftToRight ? xBounds.first : xBounds.second, !bottomToTop ? yBounds.first : yBounds.second);
            auto currP = maxCorner;
            T currV = Dot(c, currP);
            const auto minCornerX = SDK::HalfPlane<T>(SDK::Line<T>(1., 0., minCorner.x), !leftToRight);
            const auto minCornerY = SDK::HalfPlane<T>(SDK::Line<T>(0., 1., minCorner.y), !bottomToTop);
            const auto maxCornerX = SDK::HalfPlane<T>(SDK::Line<T>(1., 0., maxCorner.x),  leftToRight);
            const auto maxCornerY = SDK::HalfPlane<T>(SDK::Line<T>(0., 1., maxCorner.y),  bottomToTop);

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

            const int kBoundaryPlanesCount = 4;
            auto matrixAndBoundaries = [&minCornerX, &minCornerY, &maxCornerX, &maxCornerY, &matrix, &indices, &kBoundaryPlanesCount](int i) -> const SDK::HalfPlane<T>& {
                switch (i)
                {
                case 0: return minCornerX;
                case 1: return minCornerY;
                case 2: return maxCornerX;
                case 3: return maxCornerY;
                default: return matrix[indices[i - kBoundaryPlanesCount]];
                }
            };
            std::pair<int, int> tightHalfPlanes;
            for (int i = kBoundaryPlanesCount; i < (int)matrix.size() + kBoundaryPlanesCount; ++i) {
                const SDK::HalfPlane<T>& hp = matrixAndBoundaries(i);
                if (!hp.IsInside(currP, maxDiff)) {
                    // Rotate normal 90 degrees CCW to get a separator.
                    const SDK::Point2<T> separator(SDK::Normalize(SDK::Point2<T>(-hp.GetNormal().y, hp.GetNormal().x)));
                    SDK::Point2<T> bestAgainst(details::Infinitize(separator.x, infty), details::Infinitize(separator.y, infty));
                    SDK::Point2<T> bestAlong(-bestAgainst);
                    int secondHalfPlane = -1;
                    for (int j = 0; j < i; ++j) {
                        const SDK::HalfPlane<T>& prev = matrixAndBoundaries(j);
                        const T d = Dot(prev.GetNormal(), separator);
                        const bool isParallel = SDK::AlmostEqualToZero(d, maxDiff * maxDiff);
                        if (isParallel) {
                            if (!details::IsParallelFeasible(hp, prev, maxDiff)) {
                                return nullopt;
                            }
                            continue;
                        }
                        const auto pRaw = Intersection(hp.boundary, prev.boundary, maxDiff);
                        // TODO: Clamping this way is wrong as we must ensure that the point is on the hp.boundary
                        const auto p = details::Clamp(pRaw, -infty, infty);
                        bool isAlong = d > 0;
                        if (isAlong) {
                            const auto increment = p - bestAlong;
                            if (Dot(increment, separator) > 0) {
                                bestAlong = p;
                                secondHalfPlane = j;
                            }                            
                        } else {
                            const auto increment = p - bestAgainst;
                            if (Dot(increment, separator) < 0) {
                                bestAgainst = p;
                                secondHalfPlane = j;
                            }
                        }
                    }
                    const bool bestAlongInf = bestAlong.IsInfinite();
                    const bool bestAgainstInf = bestAgainst.IsInfinite();
                    assert(!bestAlongInf || !bestAgainstInf);
                    SDK::Point2<T> candidateP(SDK::kNan);
                    T candidateV;
                    if (!bestAlongInf && !bestAgainstInf) {
                        const T cos = Dot(bestAgainst - bestAlong, separator);
                        const bool isFeasible = cos > -std::sqrt(maxDiff);
                        if (!isFeasible) {
                            return nullopt;
                        }
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
                    if (candidateV <= currV) {
                        currV = candidateV;
                        currP = candidateP;
                        tightHalfPlanes = std::make_pair(
                            i - kBoundaryPlanesCount,
                            secondHalfPlane - kBoundaryPlanesCount
                        );
                        DebugPrintf("%f, %f\n", currP.x, currP.y);
                    }
                }
            }

            if (tightHalfPlanes.first >= 0) {
                tightHalfPlanes.first = indices[tightHalfPlanes.first];
            }
            if (tightHalfPlanes.second >= 0) {
                tightHalfPlanes.second = indices[tightHalfPlanes.second];
            }
            return std::make_optional(LpTask2dRes<T>{ currP, tightHalfPlanes });
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
