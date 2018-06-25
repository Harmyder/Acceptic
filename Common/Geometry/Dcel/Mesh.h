#pragma once
#include "Container/Dynarray.h"
#include "Geometry/Dcel/HalfEdge.h"
#include "Geometry/Dcel/Vertex.h"
#include "Geometry/Dcel/Face.h"
#include <iterator>

namespace Common
{
    namespace Dcel
    {
        template <class Id>
        class Mesh
        {
        public:
            using indexing_type = Id;

            Mesh() {}

            Mesh(Id verticesCount, Id edgesCount, Id facesCount) :
                _vertices(verticesCount),
                _halfedges(edgesCount * 2),
                _faces(facesCount)
            {}

            auto& Halfedges() { return _halfedges; }
            auto& Vertices() { return _vertices; }
            auto& Faces() { return _faces; }
            const auto& Halfedges() const { return _halfedges; }
            const auto& Vertices()  const { return _vertices; }
            const auto& Faces()     const { return _faces; }

            Id GetTwinId(Id i) const { return i ^ 1; }

            Mesh(Mesh&& other) {
                *this = move(other);
            }

            Mesh& operator=(Mesh&& other) {
                _vertices = move(other._vertices);
                _halfedges = move(other._halfedges);
                _faces = move(other._faces);
                return *this;
            }

        private:
            Dynarray<Vertex<Id>>   _vertices;
            Dynarray<HalfEdge<Id>> _halfedges; // Twins are stored consequently in pairs odd-even.
            Dynarray<Face<Id>>     _faces;
        };
    }
}