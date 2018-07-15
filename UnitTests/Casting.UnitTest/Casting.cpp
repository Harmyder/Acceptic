#include "stdafx.h"
#include "CppUnitTest.h"
#include "Common/SDK/Primitives.h"
#include "Common/SDK/Hashes.h"
#include "Common/Geometry/Dcel/Tools.h"
#include "Casting/Casting.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <unordered_set>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace Common;
using namespace Common::SDK;
using namespace Casting;
using namespace std;

namespace CommonTest
{
    TEST_CLASS(Casting_Test)
    {
        auto PrepareAndRun(const vector<Point3<double>>& verticesData, const vector<int>& trianglesData) {
            const int verticesCount = (int)verticesData.size();
            const int trianglesCount = (int)trianglesData.size() / 3;
            unordered_set<pair<int, int>, pairhash> edges;
            for (int t = 0; t < trianglesCount; ++t) {
                for (int i = 0; i < 3; ++i) {
                    const int vStart = trianglesData[t * 3 + i];
                    const int vEnd = trianglesData[t * 3 + (i + 1) % 3];
                    if (edges.find(make_pair(vEnd, vStart)) == edges.end()) {
                        edges.insert(make_pair(vStart, vEnd));
                    }
                }
            }
            const auto m = Dcel::Create(trianglesData.cbegin(), trianglesData.cend(), verticesCount, (int)edges.size());
            const auto bigFaces = CombineFaces(m, verticesData);
            return bigFaces;
        }

        void twoMissed_oneEnough_Helper(double angle1, double angle2, const vector<int>& expected) {
            // Edge between up-looking and down-looking faces must reside in a plane parallel to the projection plane
            const SDK::Point3<double> hemisphereDir(0., 0., 1.);
            const SDK::Point3<double> rotationAxis(0., 1., 0.);
            SDK::Point3<double> n1(RotationAroundAxis(rotationAxis, angle1, kMaxDiff) * hemisphereDir);
            SDK::Point3<double> n2(RotationAroundAxis(rotationAxis, angle2, kMaxDiff) * hemisphereDir);
            // Let's say that one of intersection points is origin, so free-term is zero
            std::vector<Plane<double>> bigFaces{ Plane<double>(n1, 0.), Plane<double>(n2, 0.) };
            if (expected.size() > 0) {
                const auto actual = FindS2Coverage(bigFaces, &hemisphereDir, &hemisphereDir + 1);
                Assert::AreEqual(actual.size(), expected.size());
                for (int i = 0; i < (int)actual.size(); ++i) {
                    Assert::AreEqual(actual[i], expected[i]);
                }
            }
            else {
                Assert::ExpectException<logic_error>([&]() { FindS2Coverage(bigFaces, &hemisphereDir, &hemisphereDir + 1); });
            }
        }

    public:
        TEST_METHOD(CombineFaces_TwoTriangles) {
            const vector<Point3<double>> verticesData{ Point3<double>(0,0,0), {0,0,1}, {0,1,1}, {0,1,0} };
            const vector<int> trianglesData{ 0,1,2,0,2,3 };
            const auto bigFaces = PrepareAndRun(verticesData, trianglesData);
            Assert::AreEqual((int)bigFaces.size(), 1);
        }

        TEST_METHOD(CombineFaces_Cube) {
            const vector<Point3<double>> verticesData{ Point3<double>(1., -1., -1.000000),
                {  1. ,-1., -1. },
                {  1. ,-1.,  1. },
                { -1., -1.,  1. },
                { -1., -1., -1. },
                {  1.,  1., -1. },
                {  1.,  1.,  1. },
                { -1.,  1.,  1. },
                { -1.,  1., -1. }
            };
            const vector<int> trianglesData{
                2, 3, 4,
                8, 7, 6,
                1, 5, 6,
                2, 6, 7,
                7, 8, 4,
                1, 4, 8,
                1, 2, 4,
                5, 8, 6,
                2, 1, 6,
                3, 2, 7,
                3, 7, 4,
                5, 1, 8
            };
            const auto bigFaces = PrepareAndRun(verticesData, trianglesData);
            Assert::AreEqual((int)bigFaces.size(), 6);
        }

        TEST_METHOD(ProjectHemisphereOnZUnitPlane_Reverse) {
            const auto expected = Normalize(SDK::Point3<double>(1., 1., 1.), kMaxDiff);
            const auto actual = Normalize(ComputeHemisphereDirFromItsZUnitPlaneBoundary(ProjectHemisphereOnZUnitPlane(expected)), kMaxDiff);
            const double d = (expected - actual).LenSq();
            Assert::IsTrue(SDK::AlmostEqualToZero(d, kMaxDiffSq));
        }

        TEST_METHOD(FindS2Coverage_twoMissed_oneEnough_Error) {
            twoMissed_oneEnough_Helper(1., 2., {});
        }

        TEST_METHOD(FindS2Coverage_twoMissed_oneEnough_theSameLine) {
            twoMissed_oneEnough_Helper(1., 1. + M_PI, { 0,1 });
        }

        TEST_METHOD(FindS2Coverage_twoMissed_oneEnough) {
            const double margin = 0.1;
            twoMissed_oneEnough_Helper(1., 1. + M_PI + margin, { 0,1 });
        }

        TEST_METHOD(FindS2Coverage_twoMissed_oneIsnotEnough_Error) {
            const SDK::Point3<double> hemisphereDir(0., 0., 1.);
            const SDK::Point3<double> rotationAxisUpFace(0., 1., 0.);
            const SDK::Point3<double> rotationAxisDownFace1 = Normalize(SDK::Point3<double>(0., 1.,  .2), kMaxDiff);
            const SDK::Point3<double> rotationAxisDownFace2 = Normalize(SDK::Point3<double>(0., 1., -.2), kMaxDiff);
            SDK::Point3<double> n1(RotationAroundAxis(rotationAxisUpFace, -.1, kMaxDiff) * hemisphereDir);
            SDK::Point3<double> n2(RotationAroundAxis(rotationAxisDownFace1, 1. + M_PI * .75, kMaxDiff) * hemisphereDir);
            SDK::Point3<double> n3(RotationAroundAxis(rotationAxisDownFace2, 1. + M_PI * .75, kMaxDiff) * hemisphereDir);
            // Let's say that one of intersection points is origin, so free-term is zero
            std::vector<Plane<double>> bigFaces{ Plane<double>(n1, 0.), Plane<double>(n2, 0.), Plane<double>(n3, 0) };
            Assert::ExpectException<logic_error>([&]() {FindS2Coverage(bigFaces, &hemisphereDir, &hemisphereDir + 1); });
        }

        TEST_METHOD(FindS2Coverage_twoMissed_oneIsnotEnough) {
            const SDK::Point3<double> hemisphereDir(0., 0., 1.);
            const SDK::Point3<double> rotationAxisUpFace(0., 1., 0.);
            SDK::Point3<double> upFace(RotationAroundAxis(rotationAxisUpFace, 1., kMaxDiff) * hemisphereDir);
            // Move on with the rotation a bit to make the future solution feasible
            SDK::Point3<double> downFacesMiddle(RotationAroundAxis(rotationAxisUpFace, .1, kMaxDiff) * -upFace);
            // Rotate every to the sides a bit to obtain two down faces not parallel to the up face
            SDK::Point3<double> downFace1(RotationAroundAxis(hemisphereDir, .1, kMaxDiff) * downFacesMiddle);
            SDK::Point3<double> downFace2(RotationAroundAxis(hemisphereDir, -.1, kMaxDiff) * downFacesMiddle);
            // Let's say that one of intersection points is origin, so free-term is zero
            std::vector<Plane<double>> bigFaces{ Plane<double>(upFace, 0.), Plane<double>(downFace1, 0.), Plane<double>(downFace2, 0) };
            const auto actual = FindS2Coverage(bigFaces, &hemisphereDir, &hemisphereDir + 1);
            const vector<int> expected{ 0,1,2 };
            Assert::AreEqual(actual.size(), expected.size());
            for (int i = 0; i < (int)actual.size(); ++i) {
                Assert::AreEqual(actual[i], expected[i]);
            }
        }
    };
}
