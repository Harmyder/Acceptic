#include "stdafx.h"
#include "Casting.h"

#include "Common/Geometry/Dcel/Tools.h"
#include "Common/SDK/Utility.h"
#include "Common/Optimization/Lp/IncrementalLp.h"
#include "Common/DebugPrint.h"

using namespace Common;
using namespace std;

const double kMaxDiff = 1e-6;
const double kMaxDiffSq = 1e-12;

namespace
{
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
}

namespace Casting
{
    // Combine faces into bigger faces if their are adjacent and normals are the same
    // ___platoA___     ___platoB___
    // |          |_____|          |
    // |___________________________|
    // platoA and platoB will be different big faces though they are the same for casting purposes
    vector<SDK::Plane<double>> CombineFaces(const Dcel::Mesh<int>& m, vector<SDK::Point3<double>> verticesObj) {
        vector<SDK::Plane<double>> bigFaces;
        vector<bool> visited(m.Faces().size());
        vector<int> frontier;
        int curr = -1;
        const int facesCount = (int)m.Faces().size();
        while (++curr < facesCount) {
            if (visited[curr]) {
                continue;
            }
            visited[curr] = true;
            bigFaces.emplace_back(PlaneFromFace(m, verticesObj, curr));
            const auto reference = SDK::Normalize(bigFaces.back().GetNormal(), kMaxDiff);
            frontier.push_back(curr);
            auto handleFace = [&visited, &frontier, &m, &reference, &verticesObj](int faceId) {
                if (!visited[faceId]) {
                    const auto normal = SDK::Normalize(PlaneFromFace(m, verticesObj, faceId).GetNormal(), kMaxDiff);
                    // Assume that adjacent faces can't have almost oposite normals
                    if (SDK::AlmostEqualToZero((normal - reference).LenSq(), kMaxDiffSq) || SDK::AlmostEqualToZero((normal + reference).LenSq(), kMaxDiffSq)) {
                        visited[faceId] = true;
                        frontier.push_back(faceId);
                    }
                }
            };
            while (frontier.size() > 0) {
                const int i = frontier.back();
                frontier.pop_back();
                Dcel::VisitNeighbouringFaces(m, i, handleFace);
            }
        }
        for (auto& face : bigFaces) {
            face.SetNormal(Normalize(face.GetNormal(), kMaxDiff));
        }
        return bigFaces;
    }

    SDK::HalfPlane<double> ProjectHemisphereOnZUnitPlane(SDK::Point3<double> hemisphereDirection) {
        const auto inplaneNormal = SDK::Point2<double>(hemisphereDirection.x, hemisphereDirection.y);
        const double distanceToOrigin = hemisphereDirection.z / inplaneNormal.Len();
        SDK::Point2<double> pointOnLine = -inplaneNormal * distanceToOrigin;
        const double freeTerm = Dot(pointOnLine, inplaneNormal);
        SDK::Line<double> boundary(inplaneNormal, freeTerm);
        SDK::HalfPlane<double> hp(boundary);
        return hp;
    }

    SDK::Point3<double> ComputeHemisphereDirFromItsZUnitPlaneBoundary(SDK::HalfPlane<double> boundary) {
        const double freeTerm = boundary.boundary.c;
        const auto inplaneNormal = boundary.boundary.GetNormal();
        const double distanceToOrigin = -freeTerm / inplaneNormal.LenSq();
        SDK::Point3<double> dir(inplaneNormal.x, inplaneNormal.y, distanceToOrigin * inplaneNormal.Len());
        return dir;
    }

    array<int, kCandidatesPerHemisphere> FindCover(const vector<SDK::HalfPlane<double>>& matrix) {
        SDK::Point2<double> c(1., 1.);
        const auto kBounds = make_pair(-kBound, kBound);
        const auto p1 = Optimization::solveMax(c, matrix, kBounds, kBounds, kMaxDiff).value();
        const auto p2 = Optimization::solveMax(-c, matrix, kBounds, kBounds, kMaxDiff).value();

        assert(std::abs(p1.point.x) != kBound && std::abs(p1.point.y) != kBound);
        assert(std::abs(p2.point.x) != kBound && std::abs(p2.point.y) != kBound);

        vector<int> tights{
            p1.formingHalfPlanes.first, p1.formingHalfPlanes.second,
            p2.formingHalfPlanes.first, p2.formingHalfPlanes.second
        };
        sort(tights.begin(), tights.end());
        tights.erase(unique(tights.begin(), tights.end()), tights.end());

        // TODO: Fix solver to return only halfplanes from matrix
        assert(find_if(tights.cbegin(), tights.cend(), [](int i) { return i < 0; }) == tights.cend());

        if (tights.size() == 3) {
            assert(SDK::IsPlaneCover(matrix[tights[0]], matrix[tights[1]], matrix[tights[2]], kMaxDiff));
            return { tights[0], tights[1], tights[2] };
        }

        if (SDK::IsPlaneCover(matrix[tights[0]], matrix[tights[1]], matrix[tights[2]], kMaxDiff)) return { tights[0],tights[1],tights[2] };
        if (SDK::IsPlaneCover(matrix[tights[0]], matrix[tights[1]], matrix[tights[3]], kMaxDiff)) return { tights[0],tights[1],tights[3] };
        if (SDK::IsPlaneCover(matrix[tights[0]], matrix[tights[2]], matrix[tights[3]], kMaxDiff)) return { tights[0],tights[2],tights[3] };
        if (SDK::IsPlaneCover(matrix[tights[1]], matrix[tights[2]], matrix[tights[3]], kMaxDiff)) return { tights[1],tights[2],tights[3] };

        throw logic_error("");
    }

    bool IsTheSameFaceNormal(const SDK::Point3<double>& n1, const SDK::Point3<double>& n2) {
        return Dot(n1, n2) > Casting::kCosSmallestHalfAngle;
    }
}
