#include "stdafx.h"
#include "CppUnitTest.h"
#include "Common/SDK/Primitives.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace Common;
using namespace Common::SDK;

namespace CommonTest
{
    TEST_CLASS(Primitives_Test)
    {
        void TestIntersection(const Line<double>& l1, const Line<double>& l2, Point2<double> expected) {
            const auto actual = Intersection(l1, l2, maxDiff);
            if (expected.x == std::numeric_limits<double>::infinity() && expected.y == std::numeric_limits<double>::infinity()) {
                Assert::IsTrue(actual.x == std::numeric_limits<double>::infinity() && actual.y == std::numeric_limits<double>::infinity());
            }
            else {
                const auto dist = (actual - expected).Len();
                Assert::IsTrue(AlmostEqualToZero(dist, maxDiff));
            }
        }

        void TestIsPlaneCover(const Line<double>& l1, const Line<double>& l2, const Line<double>& l3, bool expected, bool flip1 = false, bool flip2 = false, bool flip3 = false) {
            const HalfPlane<double> hp1(l1, flip1);
            const HalfPlane<double> hp2(l2, flip2);
            const HalfPlane<double> hp3(l3, flip3);
            const bool isCover = IsPlaneCover(hp1, hp2, hp3, maxDiff);
            Assert::IsTrue(expected == isCover);
        }

    public:
        TEST_METHOD(IntersectionVerVer) {
            const Line<double> l(1., 0., 1.);
            TestIntersection(l, l, Point2<double>(kInfinity));
        }

        TEST_METHOD(IntersectionHorHor) {
            const Line<double> l(1., 0., 1.);
            TestIntersection(l, l, Point2<double>(kInfinity));
        }

        TEST_METHOD(IntersectionVerHor) {
            const Line<double> ver(1., 0., 1.);
            const Line<double> hor(0., 1., 1.);
            const Point2<double> expected(1., 1.);
            TestIntersection(ver, hor, expected);
            TestIntersection(hor, ver, expected);
        }

        TEST_METHOD(IntersectionSloppedSelf) {
            const Line<double> l(1., 1., 0.);
            const Point2<double> expected(kInfinity);
            TestIntersection(l, l, expected);
        }

        TEST_METHOD(IntersectionSlopped) {
            const Line<double> l1(1., -1., 0.);
            const Line<double> l2(1., 1., 0.);
            const Point2<double> expected(kZero);
            TestIntersection(l1, l2, expected);
        }

        TEST_METHOD(IsPlaneCoverByTwoLines) {
            const Line<double> l1(1., 0., 1.);   // x = 1, to the right
            const Line<double> l2(1., 0., 2.);   // x = 2, to the right
            const Line<double> l3(-1., 0., -2.); // x = 2, to the left
            TestIsPlaneCover(l1, l2, l3, true);
        }
        
        TEST_METHOD(IsPlaneCoverByTriangle) {
            const Line<double> l1(1., -1., -2.); // y =  x + 2, left leg
            const Line<double> l2( -1., -1., -2.); // y = -x + 2, right leg
            const Line<double> l3( 0., 1., 0.); // base
            TestIsPlaneCover(l1, l2, l3, true);
            TestIsPlaneCover(l1, l2, l3, false, false, false, true);
            TestIsPlaneCover(l1, l2, l3, false, false, true, false);
            TestIsPlaneCover(l1, l2, l3, false, true, false, false);
        }
    };
}
