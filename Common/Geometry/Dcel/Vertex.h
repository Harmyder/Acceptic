#pragma once

#include "Geometry/Dcel/Traits.h"

namespace Common
{
    namespace Dcel
    {
        template <class Id>
        class Vertex
        {
        public:
            void SetLeaving(Id l) { _leaving = l; }
            Id GetLeaving() const { return _leaving; }

        private:
            Id _leaving = Traits<Id>::kNoIndex;
        };
    }
}
