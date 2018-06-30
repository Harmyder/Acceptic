#include "stdafx.h"
#include "CppUnitTest.h"
#include "Common/SDK/Types.h"
#include "Common/SDK/Utility.h"
#include <string>

using namespace std;

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace Common;
using namespace Common::SDK;

namespace CommonTest
{
    TEST_CLASS(Types_Test)
    {
        template <class T, int Dim>
        void TestHpProjection(const HyperPlane<T, Dim>& hp, const Point<T, Dim>& p2p, const Point<T, Dim> expectedProj) {
            const auto actualProj = hp.Project(p2p, kMaxDiff);
            const auto dist = (expectedProj - actualProj).Len();
            Assert::IsTrue(AlmostEqualToZero(dist, kMaxDiff));
        }

    public:
        TEST_METHOD(LineVertical) {
            Line<double> l(1., 0., 1.); // x = 1
            Assert::IsTrue(l.IsVertical(kMaxDiff));
        }

        TEST_METHOD(LineHorizontal) {
            Line<double> l(0., 1., 1.); // y = 1
            Assert::IsTrue(l.IsHorizontal(kMaxDiff));
        }

        TEST_METHOD(HyperPlane2dDistanceVertical) {
            Line<double> l(1., 0., 1.);
            assert(l.IsVertical(kMaxDiff));
            const auto d = l.Distance(Point2<double>(kZero), kMaxDiff);
            Assert::IsTrue(AlmostEqualRelativeAndAbs(d, l.c, kMaxDiff));
        }

        TEST_METHOD(HyperPlane2dDistanceHorizontal) {
            Line<double> l(0., 1., 2.);
            assert(l.IsHorizontal(kMaxDiff));
            const auto d = l.Distance(Point2<double>(kZero), kMaxDiff);
            Assert::IsTrue(AlmostEqualRelativeAndAbs(d, l.c, kMaxDiff));
        }
        
        TEST_METHOD(HyperPlane2dDistanceSlopped) {
            Line<double> l(1., 1., 1.);
            const auto d = l.Distance(Point2<double>(-1., 0.), kMaxDiff);
            Assert::IsTrue(AlmostEqualRelativeAndAbs(d, sqrt(2.), kMaxDiff));
        }
        
        TEST_METHOD(HyperPlane2dProjectVertical) {
            const Line<double> l(1., 0., 1.);
            assert(l.IsVertical(kMaxDiff));
            const Point2<double> p2p(kZero);
            const Point2<double> expected(1., 0.);
            TestHpProjection(l, p2p, expected);
        }

        TEST_METHOD(HyperPlane2dProjectHorizontal) {
            const Line<double> l(0., 1., 1.);
            assert(l.IsHorizontal(kMaxDiff));
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
            const bool inside = hp.IsInside(x, kMaxDiff);
            hp.flip();
            const bool insideFlip = hp.IsInside(x, kMaxDiff);
            Assert::IsFalse(inside);
            Assert::IsTrue(insideFlip);
        }

        TEST_METHOD(CreateCrossProdMatrixXYZ) {
            Point3<double> x(kXUnit);
            Point3<double> y(kYUnit);
            const auto m = CreateCrossProdMatrix(x);
            const auto actual = m * y;
            const auto expected = Cross(x, y);
            const double d = (actual - expected).LenSq();
            Assert::IsTrue(AlmostEqualToZero(d, kMaxDiffSq));
        }

        TEST_METHOD(RotationBetweenXY) {
            Point3<double> x(kXUnit);
            Point3<double> y(kYUnit);
            const auto m = RotationBetween(x, y, kMaxDiff, kMaxDiffSq);
            const auto actual = m * x;
            const double d = (actual - y).LenSq();
            Assert::IsTrue(AlmostEqualToZero(d, kMaxDiffSq));
        }

        TEST_METHOD(RotationExample) {
            Point3<double> x(kZUnit);
            Point3<double> y(0.980580675690920, 0.196116135138184, 0);
            Matrix<double, 3> expected = { Point3<double>(0.038461538461539, -0.192307692307692, 0.980580675690920),
                {-0.192307692307692,   0.961538461538461, 0.196116135138184},
                {-0.980580675690920,   -0.196116135138184,   0.}
            };
            const auto actual = RotationBetween(x, y, kMaxDiff, kMaxDiffSq);
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    Assert::IsTrue(AlmostEqualRelativeAndAbs(actual.r[i].data[j], expected.r[i].data[j], kMaxDiff, kMaxDiff), (L"i = " + to_wstring(i) + L", j = " + to_wstring(j)).c_str());
                }
            }
        }
    };
}
