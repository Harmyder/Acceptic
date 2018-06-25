#include "stdafx.h"
#include "Casting.h"

#include "Common/Geometry/Dcel/Mesh.h"
#include "Common/SDK/Utility.h"

using namespace Common;
using namespace std;

template <class T>
SDK::Plane<T> PlaneFrom2VectorsAndPoint(const SDK::Point3<T>& v1, const SDK::Point3<T>& v2, const SDK::Point3<T>& p) {
    const auto normal = Cross(v1, v2);
    const auto d = Dot(normal, p);
    const SDK::Plane<T> res(normal, d);
    return res;
}

auto PlaneFromFace(const Dcel::Mesh<int>& m, vector<SDK::Point3<double>> verticesObj, int faceId) {
    const int edgeId = m.Faces()[faceId].GetEdge();
    const auto& eA = m.Halfedges()[edgeId];
    const auto& eB = m.Halfedges()[eA.GetNext()];
    const auto& eBTwin = m.Halfedges()[m.GetTwinId(eA.GetNext())];
    const auto& vA = verticesObj[eA.GetOrigin()];
    const auto& vB = verticesObj[eB.GetOrigin()];
    const auto& vC = verticesObj[eBTwin.GetOrigin()];
    const auto plane = PlaneFrom2VectorsAndPoint(vB - vA, vC - vB, vA);
    return plane;
}

// Combine faces into bigger faces if their are adjacent and normals are the same
// ___platoA___     ___platoB___
// |          |_____|          |
// |___________________________|
// platoA and platoB will be different big faces though they are the same for casting purposes
vector<SDK::Plane<double>> CombineFaces(const Dcel::Mesh<int>& m, vector<SDK::Point3<double>> verticesObj) {
    vector<SDK::Plane<double>> bigFaces;
    vector<bool> visited(m.Faces().size());
    int curr = -1;
    while (curr < (int)m.Faces().size()) {
        while (visited[++curr]);
        visited[curr] = true;
        bigFaces.emplace_back(PlaneFromFace(m, verticesObj, curr));
        const auto reference = SDK::Normalize(bigFaces.back().GetNormal());
        auto handleFace = [&visited, &m, &reference, &verticesObj](int faceId) {
            if (!visited[faceId]) {
                const auto normal = SDK::Normalize(PlaneFromFace(m, verticesObj, faceId).GetNormal());
                if (SDK::AlmostEqualToZero((normal - reference).LenSq(), maxDiffSq)) {
                    visited[faceId] = true;
                }
            }
        };
    }
    return bigFaces;
}

//std::array<Common::SDK::Plane<double>, 3> FindCoverForHemisphere(
//    const std::vector<Common::SDK::Plane<double>>& faces,
//    const Common::SDK::Point3<double>& hemi
//) {
//    (void)hemi;
//    (void)faces;
//    return {};
//}
