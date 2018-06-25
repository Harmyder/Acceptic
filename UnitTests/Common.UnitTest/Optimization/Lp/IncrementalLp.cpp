#include "stdafx.h"
#include "CppUnitTest.h"
#include "Common/Optimization/Lp/IncrementalLp.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace Common;
using namespace Common::SDK;
using namespace Common::Optimization;

using namespace std;

namespace CommonTest
{
    TEST_CLASS(IncrementalLp_Test)
    {
        const double maxDiff   = 1e-6;
        const double maxDiffSq = 1e-12;
        const pair<double, double> bounds{ -1e3, 1e3 };
        const pair<double, double> boundsDoubled{ -2e3, 2e3 };

        //void TestIsPlaneCover(const Line<double>& l1, const Line<double>& l2, const Line<double>& l3, bool expected, bool flip1 = false, bool flip2 = false, bool flip3 = false) {
        //    const HalfPlane<double> hp1(l1, flip1);
        //    const HalfPlane<double> hp2(l2, flip2);
        //    const HalfPlane<double> hp3(l3, flip3);
        //    const bool isCover = IsPlaneCover(hp1, hp2, hp3, maxDiff);
        //    Assert::IsTrue(expected == isCover);
        //}

    public:
        TEST_METHOD(solveMaxSquare_Unbounded) {
            // Unbounded to the right
            //      -1  0  1
            // 1     *****************************
            // 0     *****************************
            //-1     *****************************
            vector<HalfPlane<double>> matrix;
            matrix.emplace_back(Line<double>(1.,  0., -1.));
            matrix.emplace_back(Line<double>(1.,  0.,  1.));
            matrix.emplace_back(Line<double>(0.,  1.,  -1.));
            matrix.emplace_back(Line<double>(0., -1.,  -1.));
            const auto best = solveMax(Point2<double>(1., 0.), matrix, bounds, boundsDoubled, maxDiff);
            Assert::IsTrue(best.has_value());
            const Point2<double> expected(bounds.second, 1.);
            const double d = (expected - best.value()).LenSq();
            Assert::IsTrue(AlmostEqualToZero(d, maxDiffSq));
        }

        TEST_METHOD(solveMaxSquare_InsideOut_Infeasible) {
            vector<HalfPlane<double>> matrix;
            matrix.emplace_back(Line<double>(-1.,  0., 1));
            matrix.emplace_back(Line<double>( 1.,  0., 1));
            const auto best = solveMax(Point2<double>(1., 0.), matrix, bounds, bounds, maxDiff);
            Assert::IsFalse(best.has_value());
        }
            
        TEST_METHOD(solveMaxSquare_Bounded) {
            vector<HalfPlane<double>> matrix;
            matrix.emplace_back(Line<double>(-1.,  0., -1));
            matrix.emplace_back(Line<double>( 1.,  0., -1));
            matrix.emplace_back(Line<double>( 0., -1., -1));
            matrix.emplace_back(Line<double>( 0.,  1., -1));
            const auto best = solveMax(Point2<double>(1., 0.), matrix, bounds, bounds, maxDiff);
            Assert::IsTrue(best.has_value());
            const Point2<double> expected(1., 1.);
            const double d = (expected - best.value()).LenSq();
            Assert::IsTrue(AlmostEqualToZero(d, maxDiffSq));
        }
    };
}
