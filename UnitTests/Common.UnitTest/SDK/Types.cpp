#include "stdafx.h"
#include "CppUnitTest.h"
#include "Common/SDK/Types.h"
#include "Common/SDK/Utility.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace Common;
using namespace Common::SDK;

namespace CommonTest
{
    TEST_CLASS(Types_Test)
    {
        double maxDiff = 1e-6;

        template <class T, int Dim>
        void TestHpProjection(const HyperPlane<T, Dim>& hp, const Point<T, Dim>& p2p, const Point<T, Dim> expectedProj) {
            const auto actualProj = hp.Project(p2p, maxDiff);
            const auto dist = (expectedProj - actualProj).Len();
            Assert::IsTrue(AlmostEqualToZero(dist, maxDiff));
        }

    public:
        TEST_METHOD(LineVertical) {
            Line<double> l(1., 0., 1.); // x = 1
            Assert::IsTrue(l.IsVertical(maxDiff));
        }

        TEST_METHOD(LineHorizontal) {
            Line<double> l(0., 1., 1.); // y = 1
            Assert::IsTrue(l.IsHorizontal(maxDiff));
        }

        TEST_METHOD(HyperPlane2dDistanceVertical) {
            Line<double> l(1., 0., 1.);
            assert(l.IsVertical(maxDiff));
            const auto d = l.Distance(Point2<double>(kZero), maxDiff);
            Assert::IsTrue(AlmostEqualRelativeAndAbs(d, l.c, maxDiff));
        }

        TEST_METHOD(HyperPlane2dDistanceHorizontal) {
            Line<double> l(0., 1., 2.);
            assert(l.IsHorizontal(maxDiff));
            const auto d = l.Distance(Point2<double>(kZero), maxDiff);
            Assert::IsTrue(AlmostEqualRelativeAndAbs(d, l.c, maxDiff));
        }
        
        TEST_METHOD(HyperPlane2dDistanceSlopped) {
            Line<double> l(1., 1., 1.);
            const auto d = l.Distance(Point2<double>(-1., 0.), maxDiff);
            Assert::IsTrue(AlmostEqualRelativeAndAbs(d, sqrt(2.), maxDiff));
        }
        
        TEST_METHOD(HyperPlane2dProjectVertical) {
            const Line<double> l(1., 0., 1.);
            assert(l.IsVertical(maxDiff));
            const Point2<double> p2p(kZero);
            const Point2<double> expected(1., 0.);
            TestHpProjection(l, p2p, expected);
        }

        TEST_METHOD(HyperPlane2dProjectHorizontal) {
            const Line<double> l(0., 1., 1.);
            assert(l.IsHorizontal(maxDiff));
            const Point2<double> p2p(kZero);
            const Point2<double> expected(0., 1.);
            TestHpProjection(l, p2p, expected);
        }

        TEST_METHOD(HyperPlane2dProjectSlopped) {
            Line<double> l(1., 1., 1.);
            const Point2<double> p2p(-1., 0.);
            const Point2<double> expected(0., 1.);
            TestHpProjection(l, p2p, expected);
        }

        TEST_METHOD(HalfSpace2dIsInside) {
            const Line<double> sep(1., -1., 0.);
            HalfPlane<double> hp(sep);
            const Point2<double> x(-1., 1);
            const bool inside = hp.IsInside(x, maxDiff);
            hp.flip();
            const bool insideFlip = hp.IsInside(x, maxDiff);
            Assert::IsFalse(inside);
            Assert::IsTrue(insideFlip);
        }
    };
}
