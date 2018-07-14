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
                const T distanceSq = (intCurr - intPrev).LenSq();
                const T farEnough = distanceSq > maxDiff * maxDiff;
                if (farEnough) {
                    const auto inTheMiddle = (intCurr + intPrev) / 2;
                    if (!a.IsInside(inTheMiddle, maxDiff) && !b.IsInside(inTheMiddle, maxDiff)) {
                        return false;
                    }
                }
                return true;
            }


        }

        template <class T>
        class LpSolver {
            const struct Boundary {
                SDK::Point2<T> minCorner;
                SDK::Point2<T> maxCorner;
                std::array<SDK::HalfPlane<T>, 4> hps;
            };

        public:
            LpSolver(
                const SDK::Point2<T>& c,
                const std::vector<SDK::HalfPlane<T>>& matrix, // Right side incorporated
                const std::pair<T, T>& xBounds, // Natural (domain specific bounds) on x
                const std::pair<T, T>& yBounds, // Natural (domain specific bounds) on y
                const T maxDiff
            ) :
                kC(c),
                kBoundary(ComputeBoundary(c, xBounds, yBounds)),
                kMatrix(matrix),
                _indices(matrix.size()),
                kMaxDiff(maxDiff),
                kInfty(1.0001 * std::max(std::max(std::abs(xBounds.first), std::abs(xBounds.first)), std::max(std::abs(yBounds.first), std::abs(yBounds.first))))
            {
                assert(matrix.size() > 0);
                std::iota(_indices.begin(), _indices.end(), 0);
            }

            // Protection from "bad" inputs, UTs don't want it
            void Randomize() {
                std::random_device rd;
                std::minstd_rand g(rd());
                std::shuffle(_indices.begin(), _indices.end(), g);
            }

            struct LpTask2dRes {
                SDK::Point2<T> point;
                std::pair<int, int> formingHalfPlanes;
            };

            std::optional<LpTask2dRes> solveMax() {
                auto currP = kBoundary.maxCorner;
                T currV = Dot(kC, currP);

                std::pair<int, int> tightHalfPlanes{-1, -1};
                for (int i = (int)kBoundary.hps.size(); i < (int)kMatrix.size() + (int)kBoundary.hps.size(); ++i) {
                    const SDK::HalfPlane<T>& hp = matrixAndBoundaries(i);
                    if (!hp.IsInside(currP, kMaxDiff)) {
                        bool isFeasible = FindMaxOnHalfPlane(i, currP, currV, tightHalfPlanes);
                        if (!isFeasible) {
                            return nullopt;
                        }
                    }
                }

                if (tightHalfPlanes.first >= 0) {
                    tightHalfPlanes.first = _indices[tightHalfPlanes.first];
                }
                if (tightHalfPlanes.second >= 0) {
                    tightHalfPlanes.second = _indices[tightHalfPlanes.second];
                }
                return std::make_optional(LpTask2dRes{ currP, tightHalfPlanes });
            }

        private:
            bool FindMaxOnHalfPlane(const int curr, SDK::Point2<T>& currP, T& currV, std::pair<int, int>& tightHalfPlanes) {
                const SDK::HalfPlane<T>& hp = matrixAndBoundaries(curr);
                // Rotate normal 90 degrees CCW to get a separator.
                const SDK::Point2<T> separator(SDK::Normalize(SDK::Point2<T>(-hp.GetNormal().y, hp.GetNormal().x), kMaxDiff));
                SDK::Point2<T> bestAgainst(SDK::Infinitize(separator.x, kInfty), SDK::Infinitize(separator.y, kInfty));
                SDK::Point2<T> bestAlong(-bestAgainst);
                int bestSecondHalfPlaneAlong = std::numeric_limits<int>::min();
                int bestSecondHalfPlaneAgainst = std::numeric_limits<int>::min();
                for (int j = 0; j < curr; ++j) {
                    const SDK::HalfPlane<T>& prev = matrixAndBoundaries(j);
                    const T d = Dot(prev.GetNormal(), separator);
                    //TODO: should maxDiff be parameterized on the region bounds?
                    const bool isParallel = SDK::AlmostEqualToZero(d, kMaxDiff);
                    if (isParallel) {
                        if (!details::IsParallelFeasible(hp, prev, sqrt(kMaxDiff))) {
                            return false;
                        }
                        continue;
                    }
                    const auto pRaw = Intersection(hp.boundary, prev.boundary, kMaxDiff);
                    const auto p = ClampOnLine(pRaw, hp.boundary);
                    bool isAlong = d > 0;
                    if (isAlong) {
                        const auto increment = p - bestAlong;
                        if (Dot(increment, separator) > 0) {
                            bestAlong = p;
                            bestSecondHalfPlaneAlong = j;
                        }
                    }
                    else {
                        const auto increment = p - bestAgainst;
                        if (Dot(increment, separator) < 0) {
                            bestAgainst = p;
                            bestSecondHalfPlaneAgainst = j;
                        }
                    }
                }
                const bool bestAlongInf = bestAlong.IsInfinite();
                const bool bestAgainstInf = bestAgainst.IsInfinite();
                assert(!bestAlongInf || !bestAgainstInf);
                SDK::Point2<T> candidateP(SDK::kNan);
                T candidateV;
                int candidateSecondHalfPlane;
                if (!bestAlongInf && !bestAgainstInf) {
                    const T cos = Dot(bestAgainst - bestAlong, separator);
                    const bool isFeasible = cos > -std::sqrt(kMaxDiff);
                    if (!isFeasible) {
                        return false;
                    }
                    const T againstV = Dot(bestAgainst, kC);
                    const T alongV = Dot(bestAlong, kC);
                    if (alongV > againstV) {
                        candidateV = alongV;
                        candidateP = bestAlong;
                        candidateSecondHalfPlane = bestSecondHalfPlaneAlong;
                    }
                    else {
                        candidateV = againstV;
                        candidateP = bestAgainst;
                        candidateSecondHalfPlane = bestSecondHalfPlaneAgainst;
                    }
                }
                else if (bestAlong.IsInfinite()) {
                    candidateV = Dot(kC, bestAgainst);
                    candidateP = bestAgainst;
                    candidateSecondHalfPlane = bestSecondHalfPlaneAgainst;
                }
                else {
                    assert(bestAgainst.IsInfinite());
                    candidateV = Dot(kC, bestAlong);
                    candidateP = bestAlong;
                    candidateSecondHalfPlane = bestSecondHalfPlaneAlong;
                }
                if (candidateV <= currV) {
                    currV = candidateV;
                    currP = candidateP;
                    tightHalfPlanes = std::make_pair(
                        curr - (int)kBoundary.hps.size(),
                        candidateSecondHalfPlane - (int)kBoundary.hps.size()
                    );
                    const auto interRaw = Intersection(hp.boundary, matrixAndBoundaries(candidateSecondHalfPlane).boundary, kMaxDiff);
                    const auto inter = ClampOnLine(interRaw, hp.boundary);
                    const T distance = (inter - currP).LenSq();
                    assert(SDK::AlmostEqualToZero(distance, kMaxDiff));
                }
                return true;
            }

            SDK::Point2<T> ClampOnLine(SDK::Point2<T> p, const SDK::Line<T>& line) {
                return SDK::ClampOnLine<T>(p, line, std::minmax(kBoundary.minCorner.x, kBoundary.maxCorner.x), std::minmax(kBoundary.minCorner.y, kBoundary.maxCorner.y));
            }
            Boundary ComputeBoundary(const SDK::Point2<T>& c, const std::pair<T, T>& xBounds, const std::pair<T, T>& yBounds) {
                assert(xBounds.second > xBounds.first && yBounds.second > xBounds.first);
                bool leftToRight = c.x > 0;
                bool bottomToTop = c.y > 0;
                Boundary b {
                    { leftToRight ? xBounds.first : xBounds.second, bottomToTop ? yBounds.first : yBounds.second },
                    { !leftToRight ? xBounds.first : xBounds.second, !bottomToTop ? yBounds.first : yBounds.second },
                    {
                        SDK::HalfPlane<T>(SDK::Line<T>(1., 0., b.minCorner.x), !leftToRight),
                        SDK::HalfPlane<T>(SDK::Line<T>(0., 1., b.minCorner.y), !bottomToTop),
                        SDK::HalfPlane<T>(SDK::Line<T>(1., 0., b.maxCorner.x), leftToRight),
                        SDK::HalfPlane<T>(SDK::Line<T>(0., 1., b.maxCorner.y), bottomToTop)
                    }
                };
                return b;
            }

            const SDK::HalfPlane<T>& matrixAndBoundaries(int i) {
                if (i < kBoundary.hps.size()) {
                    return kBoundary.hps[i];
                }
                else {
                    return kMatrix[_indices[i - (int)kBoundary.hps.size()]];
                }
            };

            const SDK::Point2<T> kC;
            const int kMinCornerX = 0;
            const int kMinCornerY = 1;
            const int kMaxCornerX = 2;
            const int kMaxCornerY = 3;
            const Boundary kBoundary;
            const std::vector<SDK::HalfPlane<T>> kMatrix;
            std::vector<int> _indices;
            const T kMaxDiff;
            const T kInfty;
        };

    }
}
