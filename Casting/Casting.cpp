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
        //for (auto hp : matrix) {
        //    DebugPrintf("y = %4f*x + %4f\n", -hp.boundary.a / hp.boundary.b, hp.boundary.c / hp.boundary.b);
        //}
        SDK::Point2<double> c(1., 1.);
        const auto kBounds = make_pair(-kBound, kBound);
        const auto p1 = Optimization::solveMax(c, matrix, kBounds, kBounds, kMaxDiff).value();
        const auto p2 = Optimization::solveMax(-c, matrix, kBounds, kBounds, kMaxDiff).value();
        vector<int> tights{
            p1.formingHalfPlanes.first, p1.formingHalfPlanes.second,
            p2.formingHalfPlanes.first, p2.formingHalfPlanes.second
        };
        sort(tights.begin(), tights.end());
        tights.erase(unique(tights.begin(), tights.end()), tights.end());
        const int boundaryHpCount = (int)count_if(tights.cbegin(), tights.cend(), [](int i) { return i < 0; });
        assert(boundaryHpCount <= 2);
        copy(tights.cbegin() + boundaryHpCount, tights.cend(), tights.begin());
        tights.erase(tights.end() - boundaryHpCount, tights.end());

        if (tights.size() == 3) {
            assert(SDK::IsPlaneCover(matrix[tights[0]], matrix[tights[1]], matrix[tights[2]], kMaxDiff));
            return { tights[0], tights[1], tights[2], tights[2] };
        }
        else if (tights.size() == 2) {
            assert(SDK::IsPlaneCover(matrix[tights[0]], matrix[tights[1]], matrix[tights[1]], kMaxDiff));
            return { tights[0], tights[1], tights[1], tights[1] };
        }
        else {
            assert(tights.size() == 4);
            assert (SDK::IsPlaneCover(matrix[tights[0]], matrix[tights[1]], matrix[tights[2]], kMaxDiff) ||
                    SDK::IsPlaneCover(matrix[tights[0]], matrix[tights[1]], matrix[tights[3]], kMaxDiff) ||
                    SDK::IsPlaneCover(matrix[tights[0]], matrix[tights[2]], matrix[tights[3]], kMaxDiff) ||
                    SDK::IsPlaneCover(matrix[tights[1]], matrix[tights[2]], matrix[tights[3]], kMaxDiff));
            return { tights[0], tights[1], tights[2], tights[3] };
        }
    }

    bool IsTheSameFaceNormal(const SDK::Point3<double>& n1, const SDK::Point3<double>& n2) {
        return Dot(n1, n2) > Casting::kCosSmallestHalfAngle;
    }

    SDK::Matrix<double, 3> ProjectFacesOnUpperHemisphere(
        const SDK::Point3<double>& actualNormal,
        const std::vector<SDK::Plane<double>>& bigFaces,
        const int excludeIndex,
        vector<SDK::HalfPlane<double>>& matrix,
        vector<int>& indices
    ) {
        ASSERT_NORMAL(actualNormal, kMaxDiff, kMaxDiffSq);
        const auto rotation = SDK::RotationBetween(actualNormal, Casting::upperHemisphereDirection, kMaxDiff, kMaxDiffSq);
        assert(SDK::AlmostEqualToZero((rotation * actualNormal - upperHemisphereDirection).LenSq(), kMaxDiffSq));
        for (int j = 0; j < (int)bigFaces.size(); ++j) {
            if (j != excludeIndex) {
                const auto faceNormal = bigFaces[j].GetNormal();
                assert((faceNormal - actualNormal).LenSq() > kMaxDiffSq);
                const double cos = Dot(faceNormal, actualNormal);
                const bool isParticipateInCover = (cos >= 0);
                if (isParticipateInCover) {
                    const auto rotatedFaceNormal = rotation * faceNormal;
//                    DebugPrintf("rotated [%f %f %f]\n", rotatedFaceNormal.x, rotatedFaceNormal.y, rotatedFaceNormal.z);
                    assert(SDK::AlmostEqualRelativeAndAbs(Dot(rotatedFaceNormal, Casting::upperHemisphereDirection), cos, kMaxDiff));
                    const auto halfPlane = Casting::ProjectHemisphereOnZUnitPlane(rotatedFaceNormal);
                    matrix.push_back(halfPlane);
                    indices.push_back(j);
                }
            }
        }
        return rotation;
    }

    vector<int> FindS2Coverage(const vector<SDK::Plane<double>>& bigFaces) {
        vector<int> candidates;
        vector<SDK::HalfPlane<double>> matrix;
        vector<int> indices;
        for (int i = 0; i < Casting::kCoverHemispheresCount; ++i) {
            matrix.clear();
            indices.clear();
            const auto actualNormal = Casting::hemisphereCoverForSphere[i];
            bool singleElementCoverage = false;
            for (int j = 0; j < (int)bigFaces.size(); ++j) {
                const auto bigFaceNormal = bigFaces[j].GetNormal();
                if (IsTheSameFaceNormal(bigFaceNormal, actualNormal)) {
                    if (singleElementCoverage) {
                        candidates.pop_back();
                        break;
                    }
                    else {
                        singleElementCoverage = true;
                        candidates.push_back(j);
                    }
                }
            }
            if (!singleElementCoverage) {
                ProjectFacesOnUpperHemisphere(actualNormal, bigFaces, -1, matrix, indices);

                //DebugPrintf("actualNomral [%f %f %f]\n", actualNormal.x, actualNormal.y, actualNormal.z);
                //for (int k : indices) {
                //    DebugPrintf("bigFaces [%f %f %f]\n", bigFaces[k].GetNormal().x, bigFaces[k].GetNormal().y, bigFaces[k].GetNormal().z);
                //}
                //for (auto k : matrix) {
                //    DebugPrintf("matrix [%f %f %f]\n", k.boundary.a, k.boundary.b, k.boundary.c);
                //}
                const auto coverIndices = FindCover(matrix);
                for (int ci : coverIndices) {
                    candidates.push_back(indices[ci]);
                }
            }
        }
        sort(candidates.begin(), candidates.end());
        candidates.erase(unique(candidates.begin(), candidates.end()), candidates.end());
        return candidates;
    }

    vector<pair<int, SDK::Point3<double>>> CheckCandidates(const vector<SDK::Plane<double>>& bigFaces, const vector<int>& candidates) {
        vector<pair<int, SDK::Point3<double>>> pullouts;
        vector<SDK::HalfPlane<double>> matrix;
        vector<int> indices;
        for (const int candidateIndex : candidates) {
            const auto actualNormal = bigFaces[candidateIndex].GetNormal();
            matrix.clear();
            indices.clear();
            bool exists = true;
            for (int j = 0; j < (int)bigFaces.size(); ++j) {
                if (j != candidateIndex) {
                    const auto bigFaceNormal = bigFaces[j].GetNormal();
                    const double cos = Dot(bigFaceNormal, actualNormal);
                    if (cos > kSinSmallestHalfAngle) {
                        exists = false;
                        break;
                    }
                }
            }
            if (exists) {
                const auto rotation = ProjectFacesOnUpperHemisphere(actualNormal, bigFaces, candidateIndex, matrix, indices);
                SDK::Point2<double> any{ 1.,1. };
                const auto kBounds = make_pair(-kBound, kBound);
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
                        const auto rotatedNormal = ComputeHemisphereDirFromItsZUnitPlaneBoundary(hemisphereBoundary);
                        pullOutDir = Transpose(rotation) * rotatedNormal;
                    }
                    pullouts.push_back(make_pair(candidateIndex, SDK::Normalize(pullOutDir, kMaxDiff)));
                }
            }
        }
        return pullouts;
    }

}
