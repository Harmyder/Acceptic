#pragma once
#include <limits>
#include <array>

#include "SDK/Utility.h"

namespace Common {
    namespace SDK {
        enum ZeroTag { kZero };
        enum InfinityTag { kInfinity };
        enum XUnitTag { kXUnit };
        enum YUnitTag { kYUnit };
        enum ZUnitTag { kZUnit };
        enum WUnitTag { kWUnit };

        // ************************************************************************************
        // Point
        // ************************************************************************************

        template <class T, int Dim> struct Point;

        template <class T>
        struct Point<T, 2> {
            using value_type = T;
            union {
                struct { T x, y; };
                std::array<T, 2> data;
            };

            Point(T in_x, T in_y) : x(in_x), y(in_y) {}
            explicit Point(const T* data) : Point(data[0], data[1]) {}
            explicit Point(ZeroTag)     : Point(0, 0) {}
            explicit Point(InfinityTag) : Point(std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity()) {}
            explicit Point(XUnitTag)    : Point(1, 0) {}
            explicit Point(YUnitTag)    : Point(0, 1) {}

            bool IsInfinite() const {
                const auto it = std::find_if(data.cbegin(), data.cend(), [](const T t) { return std::isinf(t); });
                return it != data.cend();
            }

            Point operator-() const { return { -x, -y }; }

            Point operator+(const Point& o) const {
                assert(Overflow<T>::Sum(x, o.x) == IsOverflow::No && Overflow<T>::Sum(y, o.y) == IsOverflow::No);
                return { x + o.x, y + o.y };
            }

            Point operator-(const Point& o) const {
                assert(Overflow<T>::Sum(x, -o.x) == IsOverflow::No && Overflow<T>::Sum(y, -o.y) == IsOverflow::No);
                return { x - o.x, y - o.y };
            }

            Point operator*(const T s) const {
                assert(Overflow<T>::Mul(x, s) == IsOverflow::No && Overflow<T>::Mul(y, s) == IsOverflow::No);
                return { x * s, y * s };
            }

            Point operator/(const T s) const {
                assert(Overflow<T>::Div(x, s) == IsOverflow::No && Overflow<T>::Div(y, s) == IsOverflow::No);
                return { x / s, y / s };
            }

            Point& operator+= (Point p) { *this = *this + p; return *this; }
            Point& operator-= (Point p) { *this = *this - p; return *this; }
            Point& operator*= (Point p) { *this = *this * p; return *this; }
            Point& operator/= (Point p) { *this = *this / p; return *this; }

            T Len() { return std::sqrt(x * x + y * y); }
            T LenSqr() { return x * x + y * y; }
        };

        template <class T> using Point2 = Point<T, 2>;

        template <class T>
        struct Point<T, 3> {
            using value_type = T;
            T x, y, z;

            Point(T in_x, T in_y, T in_z) : x(in_x), y(in_y), z(in_z) {}
            explicit Point(const T* data) : Point(data[0], data[1], data[2]) {}
            explicit Point(ZeroTag)     : Point(0, 0, 0) {}
            explicit Point(InfinityTag) : Point(std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity()) {}
            explicit Point(XUnitTag)    : Point(1, 0, 0) {}
            explicit Point(YUnitTag)    : Point(0, 1, 0) {}
            explicit Point(ZUnitTag)    : Point(0, 0, 1) {}

            Point operator-() const { return { -x, -y, -z }; }

            Point operator+(const Point& o) const {
                assert(Overflow<T>::Sum(x, o.x) == IsOverflow::No && Overflow<T>::Sum(y, o.y) == IsOverflow::No && Overflow<T>::Sum(z, o.z) == IsOverflow::No);
                return { x + o.x, y + o.y, z + o.z };
            }

            Point operator-(const Point& o) const {
                assert(Overflow<T>::Sum(x, -o.x) == IsOverflow::No && Overflow<T>::Sum(y, -o.y) == IsOverflow::No && Overflow<T>::Sum(z, -o.z) == IsOverflow::No);
                return { x - o.x, y - o.y, z - o.z };
            }

            Point operator*(const T s) const {
                assert(Overflow<T>::Mul(x, s) == IsOverflow::No && Overflow<T>::Mul(y, s) == IsOverflow::No && Overflow<T>::Mul(z, s) == IsOverflow::No);
                return { x * s, y * s, z * s };
            }

            Point operator/(const T s) const {
                assert(Overflow<T>::Div(x, s) == IsOverflow::No && Overflow<T>::Div(y, s) == IsOverflow::No && Overflow<T>::Div(z, s) == IsOverflow::No);
                return { x / s, y / s, z / s };
            }

            Point& operator+= (Point p) { *this = *this + p; return *this; }
            Point& operator-= (Point p) { *this = *this - p; return *this; }
            Point& operator*= (Point p) { *this = *this * p; return *this; }
            Point& operator/= (Point p) { *this = *this / p; return *this; }

            T LenSqr() { return x * x + y * y + z * z; }
        };

        template <class T> using Point3 = Point<T, 3>;

        template <class T> T Dot(const Point2<T>& v1, const Point2<T>& v2) { return v1.x * v2.x + v1.y * v2.y; }
        template <class T> T Dot(const Point3<T>& v1, const Point3<T>& v2) { return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z; }
            
        // ************************************************************************************
        // HyperPlane
        // ************************************************************************************
        template <class T, int Dim> struct HyperPlane;

        template <class T>
        struct HyperPlane<T, 2> {
            // a*x + b*y = c
            HyperPlane(T a, T b, T c) : a(a), b(b), c(c) {}

            using value_type = T;
            union {
                struct { T a, b, c; };
                std::array<T, 3> data;
            };

            T ComputeX(T y) const { return (c - b * y) / a; }
            T ComputeY(T x) const { return (c - a * x) / b; }

            bool IsHorizontal(T maxDiff) const { return AlmostEqualToZero(a, maxDiff); }
            bool IsVertical(T maxDiff) const { return AlmostEqualToZero(b, maxDiff); }

            Point2<T> Project(const Point2<T>& p, T maxDiff) const {
                if (IsHorizontal(maxDiff)) {
                    return { p.x, c / b };
                }
                const T y = (p.y * a * a - b * a * p.x + b * c) / (b * b + a * a);
                const T x = ComputeX(y);
                return { x, y };
            }

            T Distance(const Point2<T>& p, T maxDiff) const {
                if (IsVertical(maxDiff)) {
                    return std::fabs(p.x - c / a);
                }

                const T dist = std::fabs(a*p.x + b * p.y - c) / std::sqrt(a*a + b * b);
                return dist;
            }
        };

        template <class T> using Line = HyperPlane<T, 2>;

        template <class T>
        struct HyperPlane<T, 3> {
            // a*x + b*y + c*z = d
            HyperPlane(T a, T b, T c, T d) : a(a), b(b), c(c), d(d) {}

            using value_type = T;
            union {
                struct { T a, b, c, d; };
                std::array<T, 4> data;
            };

            T ComputeX(T y, T z) const { return (c - b * y - c * z) / a; }
            T ComputeY(T x, T z) const { return (c - a * x - c * z) / b; }
            T ComputeZ(T x, T y) const { return (c - a * x - b * y) / c; }
        };

        template <class T> using Plane = HyperPlane<T, 3>;

        // ************************************************************************************
        // HalfSpace
        // ************************************************************************************

        template <class T, int Dim>
        struct HalfSpace
        {
            // Hyper planes (1, 1, 1) and (-1, -1, -1) define the same plane, but different half spaces
            explicit HalfSpace(const HyperPlane<T, Dim>& hp, bool needsFlip = false) : hp(hp) {
                if (needsFlip) {
                    flip();
                }
            }

            HyperPlane<T, Dim> hp;

            void flip() {
                std::transform(hp.data.cbegin(), hp.data.cend(), hp.data.begin(), [](const T& t) { return -t; });
            }

            Point<T, Dim> GetNormal() const {
                Point<T, Dim> normal(hp.data.data());
                return normal;
            }

            bool IsInside(const Point<T, Dim>& p, T maxDiff) const {
                const auto base = hp.Project(p, maxDiff);
                const auto dir = p - base;
                const T d = Dot(dir, GetNormal());
                return d >= 0;
            }
        };

        template <class T> using HalfPlane = HalfSpace<T, 2>;

   }
}
