#include "stdafx.h"
#include "CppUnitTest.h"
#include "Common/SDK/Primitives.h"
#include "Common/SDK/Hashes.h"
#include "Common/Geometry/Dcel/Tools.h"
#include "Casting/Casting.h"

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
            const auto actual = ComputeHemisphereDirFromItsZUnitPlaneBoundary(ProjectHemisphereOnZUnitPlane(expected));
            const double d = (expected - actual).LenSq();
            Assert::IsTrue(SDK::AlmostEqualToZero(d, kMaxDiffSq));
        }
    };
}
