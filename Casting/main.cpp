#include "stdafx.h"

#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include "Common/SDK/Types.h"
#include "Common/Geometry/Dcel/Tools.h"
#include "Common/DebugPrint.h"
#include "Casting.h"

using namespace Common;
using namespace std;

void Read(const string& fn, vector<SDK::Point3<double>>& verticesObj, vector<int>& trianglesObj) {
    ifstream infile(fn);
    assert(!infile.fail());

    string line;
    while (getline(infile, line))
    {
        istringstream iss(line);
        string start;
        iss >> start;
        if (start == "v") {
            float x, y, z;
            iss >> x >> y >> z;
            verticesObj.emplace_back(-z, y, x);
        }
        else if (start == "f") {
            int a, b, c;
            iss >> a >> b >> c;
            trianglesObj.push_back(a - 1);
            trianglesObj.push_back(b - 1);
            trianglesObj.push_back(c - 1);
        }
    }

    for (size_t i = 0; i < trianglesObj.size() / 3; ++i) {
        const int vA = trianglesObj[i * 3];
        const int vB = trianglesObj[i * 3 + 1];
        const int vC = trianglesObj[i * 3 + 2];
        const auto side1 = verticesObj[vB] - verticesObj[vA];
        const auto side2 = verticesObj[vC] - verticesObj[vA];
        auto c = SDK::Cross(side1, side2);
        assert(c.LenSq() > 1e-12);
    }
}

void Run(const string& fn) {
    DebugPrintf("%s\n\n", fn.c_str());

    vector<SDK::Point3<double>> verticesObj;
    vector<int> trianglesObj;
    Read(fn, verticesObj, trianglesObj);

    const int verticesCount = (int)verticesObj.size();
    const int trianglesCount = (int)trianglesObj.size() / 3;

    unordered_set<pair<int, int>, Common::pairhash> edges;
    for (int t = 0; t < trianglesCount; ++t) {
        for (int i = 0; i < 3; ++i) {
            const int vStart = trianglesObj[t * 3 + i];
            const int vEnd = trianglesObj[t * 3 + (i + 1) % 3];
            if (edges.find(make_pair(vEnd, vStart)) == edges.end()) {
                edges.insert(make_pair(vStart, vEnd));
            }
        }
    }
    const auto m = Dcel::Create(trianglesObj.cbegin(), trianglesObj.cend(), verticesCount, (int)edges.size());

    const auto bigFaces = Casting::CombineFaces(m, verticesObj);

    const vector<int> candidates = Casting::FindS2Coverage(bigFaces, Casting::kHemisphereCoverForSphere.data(), Casting::kHemisphereCoverForSphere.data() + Casting::kHemisphereCoverForSphere.size());
    DebugPrintf("Candidates count: %d\n\n", candidates.size());

    const auto pullouts = Casting::CheckCandidates(bigFaces, candidates);

    for (const auto& pullout : pullouts) {
        const auto faceNormal = bigFaces[pullout.first].GetNormal();
        DebugPrintf("Face normal is     (%10f, %10f, %10f)\n", faceNormal.x, faceNormal.y, faceNormal.z);
        DebugPrintf("Pull-out normal is (%10f, %10f, %10f)\n\n", pullout.second.x, pullout.second.y, pullout.second.z);
    }

    DebugPrintf("\n\n\n\n");
}

int main()
{
    Run("../Assets/cube_1x1x1.obj");
    Run("../Assets/cube_2x2x2.obj");
    Run("../Assets/cube_10x10x10.obj");
    Run("../Assets/tapered_tube.obj");
    Run("../Assets/half_tapered_cube_2x.obj");
    Run("../Assets/tapered_cube_2x.obj");

    return 0;
}


