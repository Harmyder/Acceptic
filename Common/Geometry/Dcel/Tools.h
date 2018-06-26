#pragma once
#include "Geometry/Dcel/Mesh.h"
#include "Container/Dynarray.h"
#include <unordered_set>
#include "SDK/Hashes.h"

namespace Common
{
    namespace Dcel
    {
        namespace details
        {
            template <class It, class F>
            void IterateEdges(It triBegin, const typename It::value_type trianglesCount, F handler) {
                using Id = typename It::value_type;
                for (Id t = 0; t < trianglesCount; ++t) {
                    for (Id i = 0; i < 3; ++i) {
                        const Id vStart = *(triBegin + (t * 3 + i));
                        const Id vEnd = *(triBegin + (t * 3 + (i + 1) % 3));
                        handler(t, i, vStart, vEnd);
                    }
                }
            }

            template <class Id>
            void HandleEdge(Mesh<Id>& mesh, Id edgeIndex, Id faceIndex, Id vStart) {
                auto& halfedge = mesh.Halfedges()[edgeIndex];
                halfedge.SetFace(faceIndex);
                halfedge.SetOrigin(vStart);
                mesh.Faces()[faceIndex].SetEdge(edgeIndex);
                mesh.Vertices()[vStart].SetLeaving(edgeIndex);
            }
        }

        template <class Id, class Func>
        void VisitNeighbouringFaces(const Mesh<Id>& m, Id faceId, Func func) {
            const auto& face = m.Faces()[faceId];
            const Id startId = face.GetEdge();
            Id currId = startId;
            do {
                const auto& edge = m.Halfedges()[m.GetTwinId(currId)];
                if (edge.GetNext() != Traits<Id>::kNoIndex) {
                    func(edge.GetFace());
                }
                currId = m.Halfedges()[currId].GetNext();
            } while (currId != startId);
        }

        template <class It>
        Mesh<typename It::value_type> Create(It triBegin, It triEnd, const typename It::value_type verticesCount, const typename It::value_type edgesCount) {
            using Id = typename It::value_type;
            assert(distance(triBegin, triEnd) % 3 == 0);
            const Id trianglesCount = static_cast<Id>(std::distance(triBegin, triEnd) / 3);
            Mesh<Id> mesh(verticesCount, edgesCount, trianglesCount);

            std::unordered_set<std::pair<Id, Id>, Common::pairhash> boundaryEdges;
            details::IterateEdges(triBegin, trianglesCount, [&boundaryEdges](const Id, const Id, const Id vStart, const Id vEnd) {
                auto it = boundaryEdges.find(make_pair(vEnd, vStart));
                if (it == end(boundaryEdges)) {
                    boundaryEdges.insert(make_pair(vStart, vEnd));
                }
                else {
                    boundaryEdges.erase(it);
                }
            }); // Now boundaryEdges contains only boundary edges

            Dynarray<Id> edgesPerVertex(verticesCount, 0);
            details::IterateEdges(triBegin, trianglesCount, [&edgesPerVertex, &boundaryEdges](const Id, const Id, const Id vStart, const Id vEnd) {
                if (vStart < vEnd && boundaryEdges.find(make_pair(vStart, vEnd)) == end(boundaryEdges)) ++edgesPerVertex[vStart];
            });

            for (const auto& edge : boundaryEdges) {
                ++edgesPerVertex[min(edge.first, edge.second)];
                // That's important to use min(..), because then writing twin we always search w/ min vertex
            }

            Dynarray<Id> edgesPerVertexAcc(verticesCount + 1);
            edgesPerVertexAcc[0] = 0;
            for (size_t v = 1; v < edgesPerVertexAcc.size(); ++v) {
                edgesPerVertexAcc[v] = edgesPerVertexAcc[v - 1] + edgesPerVertex[v - 1];
            }

            auto& availEdgesPerVertex = edgesPerVertex;
            Dynarray<Id> edgeToDcel(trianglesCount * 3, Traits<Id>::kNoIndex); // To let us lookup next location

            // I. Write down one twin for every edge
            details::IterateEdges(triBegin, trianglesCount, [&availEdgesPerVertex, &edgesPerVertexAcc, &mesh, &edgeToDcel](const Id t, const Id i, const Id vStart, const Id vEnd) {
                if (vStart > vEnd) {}
                else {
                    Id edgeIndex = edgesPerVertexAcc[vStart];
                    const Id edgesCount = edgesPerVertexAcc[vStart + 1] - edgesPerVertexAcc[vStart];
                    assert(availEdgesPerVertex[vStart] > 0);
                    const Id availEdgesCount = availEdgesPerVertex[vStart]--;
                    edgeIndex += edgesCount - availEdgesCount;
                    edgeIndex *= 2;
                    assert(mesh.Halfedges()[edgeIndex].GetFace() == Traits<Id>::kNoIndex);
                    details::HandleEdge(mesh, edgeIndex, t, vStart);
                    // It's possible that the edge is on the boundary,
                    // in this case we will not find twin for the half-edge.
                    // If we will, then we will just overwirite w/ proper values.
                    auto& halfedgeTwin = mesh.Halfedges()[edgeIndex + 1];
                    halfedgeTwin.SetFace(Traits<Id>::kOutside);
                    halfedgeTwin.SetOrigin(vEnd);
                    edgeToDcel[t * 3 + i] = edgeIndex;
                }
            });
            assert(boundaryEdges.size() != 0 || find_if(cbegin(availEdgesPerVertex), cend(availEdgesPerVertex), [](Id i) { return i > 0; }) == cend(availEdgesPerVertex));

            // II. Write down another twin for every edge
            details::IterateEdges(triBegin, trianglesCount, [&edgesPerVertexAcc, &boundaryEdges, &availEdgesPerVertex, &mesh, &edgeToDcel](const Id t, const Id i, const Id vStart, const Id vEnd) {
                if (vStart > vEnd) {
                    // Here i assume that max degree is small
                    bool found = false;
                    // Search either twin edge or empty slot
                    for (Id ei = edgesPerVertexAcc[vEnd]; ei < edgesPerVertexAcc[vEnd + 1]; ++ei) {
                        const Id edgeIndex = ei * 2 + 1;
                        if (mesh.Halfedges()[edgeIndex - 1].GetOrigin() == vEnd &&
                            mesh.Halfedges()[edgeIndex].GetOrigin() == vStart) 
                        {
                            assert(mesh.Halfedges()[edgeIndex].GetFace() == Traits<Id>::kOutside);
                            details::HandleEdge(mesh, edgeIndex, t, vStart);
                            found = true;
                        }
                        if (mesh.Halfedges()[edgeIndex - 1].GetOrigin() == Traits<Id>::kNoIndex &&
                            mesh.Halfedges()[edgeIndex].GetOrigin() == Traits<Id>::kNoIndex) {
                            assert(mesh.Halfedges()[edgeIndex].GetFace() == Traits<Id>::kNoIndex);
                            details::HandleEdge(mesh, edgeIndex, t, vStart);
                            auto& halfedgeTwin = mesh.Halfedges()[edgeIndex - 1];
                            halfedgeTwin.SetFace(Traits<Id>::kOutside);
                            halfedgeTwin.SetOrigin(vEnd);
                            assert(boundaryEdges.find(make_pair(vStart, vEnd)) != end(boundaryEdges));
                            --availEdgesPerVertex[vEnd]; // Because vEnd is less than vStart
                            found = true;
                        }
                        if (found) {
                            edgeToDcel[t * 3 + i] = edgeIndex;
                            break;
                        }
                    }
                    assert(found);
                }
            });            
            assert(boundaryEdges.size() == 0 || find_if(cbegin(availEdgesPerVertex), cend(availEdgesPerVertex), [](Id i) { return i > 0; }) == cend(availEdgesPerVertex));

            // III. Update nexts
            details::IterateEdges(triBegin, trianglesCount, [&mesh, &edgeToDcel](const Id t, const Id i, const Id, const Id) {
                const auto currEdge = (t * 3 + i);
                const auto nextEdge = (t * 3 + (i + 1) % 3);
                mesh.Halfedges()[edgeToDcel[currEdge]].SetNext(edgeToDcel[nextEdge]);
            });            
            assert(count_if(mesh.Halfedges().cbegin(), mesh.Halfedges().cend(), [](const HalfEdge<Id>& he) { return he.GetNext() == Traits<Id>::kNoIndex; }) == (Id)boundaryEdges.size());

            return mesh;
        }
    }
}
