#pragma once

#include "Geometry/Dcel/Traits.h"

namespace Common
{
    namespace Dcel
    {
        template <class Id>
        class Face
        {
        public:
            void SetEdge(Id e) { _edge = e; }
            Id GetEdge() const { return _edge; }

        private:
            Id _edge = Traits<Id>::kNoIndex;
        };
    }
}
