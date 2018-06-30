#include "stdafx.h"

#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include "Common/SDK/Types.h"
#include "Common/Geometry/Dcel/Tools.h"
#include "Common/Optimization/Lp/IncrementalLp.h"
#include "Casting.h"

using namespace Common;
using namespace std;

vector<SDK::Point3<double>> verticesObj;
vector<int> trianglesObj;

void read() {
    ifstream infile("../Assets/cube_1x1x1.obj");
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
        auto c = SDK::Cross(verticesObj[trianglesObj[i * 3 + 1]] - verticesObj[trianglesObj[i * 3]], verticesObj[trianglesObj[i * 3 + 2]] - verticesObj[trianglesObj[i * 3]]);
        assert(c.LenSq() > 1e-12);
    }
}

SDK::Matrix<double, 3> ProjectFacesOnUpperHemisphere(const SDK::Point3<double>& actualNormal, const std::vector<SDK::Plane<double>>& bigFaces, const int excludeIndex, vector<SDK::HalfPlane<double>>& matrix, vector<int>& indices) {
    ASSERT_NORMAL(actualNormal, kMaxDiff, kMaxDiffSq);
    const auto rotation = SDK::RotationBetween(actualNormal, Casting::upperHemisphereDirection, kMaxDiff, kMaxDiffSq);
    for (int j = 0; j < (int)bigFaces.size(); ++j) {
        if (j != excludeIndex) {
            const auto faceNormal = Normalize(bigFaces[j].GetNormal(), kMaxDiff);
            const double cos = Dot(faceNormal, actualNormal);
            const bool isParticipateInCover = cos >= -kMaxDiff * 100;
            if (isParticipateInCover) {
                const auto rotatedFaceNormal = rotation * faceNormal;
                assert(SDK::AlmostEqualRelativeAndAbs(Dot(rotatedFaceNormal, Casting::upperHemisphereDirection), cos, kMaxDiff));
                const auto halfPlane = Casting::ProjectHemisphereOnZUnitPlane(rotatedFaceNormal);
                matrix.push_back(halfPlane);
                indices.push_back(j);
            }
        }
    }
    return rotation;
}

int main()
{
    read();

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

    vector<int> candidates(12);
    vector<SDK::HalfPlane<double>> matrix;
    vector<int> indices;
    for (int i = 0; i < Casting::kCoverHemispheresCount; ++i) {
        matrix.clear();
        indices.clear();
        const auto actualNormal = Casting::hemisphereCoverForSphere[i];
        // Check for one element coverage
        for (int j = 0; j < (int)bigFaces.size(); ++j) {
            const auto bigFaceNormal = Normalize(bigFaces[j].GetNormal(), kMaxDiff);
            if (SDK::AlmostEqualToZero((bigFaceNormal + actualNormal).LenSq(), kMaxDiffSq)) {
                candidates.push_back(j);
                continue;
            }
        }
        ProjectFacesOnUpperHemisphere(actualNormal, bigFaces, -1, matrix, indices);
        const auto coverIndices = Casting::FindCover(matrix);
        for (int ci : coverIndices) {
            candidates.push_back(indices[ci]);
        }
    }
    // There is no varanthy that candadates' list contains all sides amenable for extraction,
    // because plane covers are not unique and have overlap for different hemispheres
    sort(candidates.begin(), candidates.end());
    candidates.erase(unique(candidates.begin(), candidates.end()), candidates.end());

    for (const int candidateIndex : candidates) {
        const auto actualNormal = Normalize(bigFaces[candidateIndex].GetNormal(), kMaxDiff);
        matrix.clear();
        indices.clear();
        const auto rotation = ProjectFacesOnUpperHemisphere(actualNormal, bigFaces, candidateIndex, matrix, indices);
        SDK::Point2<double> any{ 1.,1. };
        const auto kBounds = make_pair(-Casting::kBound, Casting::kBound);
        const auto pullout1 = Optimization::solveMax(any, matrix, kBounds, kBounds, kMaxDiff);
        const auto pullout2 = Optimization::solveMax(-any, matrix, kBounds, kBounds, kMaxDiff);
        assert(pullout1.has_value() == pullout2.has_value());
        if (pullout1.has_value()) {
            const auto p = (pullout1.value().point + pullout2.value().point) / 2;
            SDK::Point3<double> pullOutDir = actualNormal;
            if (p.LenSq() > kMaxDiffSq) {
                const auto normal = SDK::Normalize(-p, kMaxDiff);
                const double freeTerm = Dot(normal, p);
                const SDK::HalfPlane<double> hemisphereBoundary(SDK::Line<double>(normal, freeTerm));
                const auto rotatedNormal = Casting::ComputeHemisphereDirFromItsZUnitPlaneBoundary(hemisphereBoundary);
                pullOutDir = Transpose(rotation) * rotatedNormal;
            }
            DebugPrintf("Face normal is     (%10f, %10f, %10f)\n", actualNormal.x, actualNormal.y, actualNormal.z);
            DebugPrintf("Pull-out normal is (%10f, %10f, %10f)\n\n", pullOutDir.x, pullOutDir.y, pullOutDir.z);
        }
    }

    return 0;
}


