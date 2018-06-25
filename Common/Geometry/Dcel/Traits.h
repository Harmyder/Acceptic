#pragma once

namespace Common
{
    namespace Dcel
    {
        template <class T>
        struct Traits
        {
            static constexpr T kNoIndex = static_cast<T>(-1);
            static constexpr T kOutside = static_cast<T>(-2);

            bool IsInner(T t) { return (t != kNoIndex) && (t != kOutside); }
        };
    }
}