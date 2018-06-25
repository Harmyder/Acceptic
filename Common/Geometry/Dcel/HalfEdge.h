#pragma once

#include "Geometry/Dcel/Traits.h"

namespace Common
{
    namespace Dcel
    {
        template <class Id>
        class HalfEdge
        {
        public:
            void SetOrigin(Id o) { _origin = o; }
            Id GetOrigin() const { return _origin; }
            void SetNext(Id n) { _next = n; }
            Id GetNext() const { return _next; }
            void SetFace(Id f) { _face = f; }
            Id GetFace() const { return _face; }

            bool operator==(const HalfEdge& other) const {
                return (_origin == other._origin && _next == other._next && _face == other._face);
            }

        private:
            Id _origin = Traits<Id>::kNoIndex;
            Id _next = Traits<Id>::kNoIndex;
            Id _face = Traits<Id>::kNoIndex;
        };
    }
}
