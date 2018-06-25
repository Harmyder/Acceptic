#pragma once
#include <limits>
#include <array>

#include "SDK/Utility.h"

namespace Common {
    namespace SDK {
        enum ZeroTag { kZero };
        enum IdentityTag { kIdentity };
        enum InfinityTag { kInfinity };
        enum NanTag { kNan };
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
#pragma warning(disable : 4201)
                struct { T x, y; };
#pragma warning(default : 4201)
                std::array<T, 2> data;
            };

            Point(T in_x, T in_y) : x(in_x), y(in_y) {}
            explicit Point(const T* data) : Point(data[0], data[1]) {}
            explicit Point(ZeroTag)     : Point(0, 0) {}
            explicit Point(InfinityTag) : Point(std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity()) {}
            explicit Point(NanTag)      : Point(std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN()) {}
            explicit Point(XUnitTag)    : Point(1, 0) {}
            explicit Point(YUnitTag)    : Point(0, 1) {}
            Point(const Point& o) : x(o.x), y(o.y) {}

            Point& operator=(const Point& o) { x = o.x; y = o.y; return *this; }

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

            Point& operator+= (const Point& p) { *this = *this + p; return *this; }
            Point& operator-= (const Point& p) { *this = *this - p; return *this; }
            Point& operator*= (T p) { *this = *this * p; return *this; }
            Point& operator/= (T p) { *this = *this / p; return *this; }

            T Len() const { return std::sqrt(LenSq()); }
            T LenSq() const { return x * x + y * y; }
        };

        template <class T> using Point2 = Point<T, 2>;

        template <class T>
        struct Point<T, 3> {
            using value_type = T;
            union {
#pragma warning(disable : 4201)
                struct { T x, y, z; };
#pragma warning(default : 4201)
                std::array<T, 3> data;
            };

            Point(T in_x, T in_y, T in_z) : x(in_x), y(in_y), z(in_z) {}
            explicit Point(const T* data) : Point(data[0], data[1], data[2]) {}
            explicit Point(ZeroTag)     : Point(0, 0, 0) {}
            explicit Point(InfinityTag) : Point(std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity()) {}
            explicit Point(XUnitTag)    : Point(1, 0, 0) {}
            explicit Point(YUnitTag)    : Point(0, 1, 0) {}
            explicit Point(ZUnitTag)    : Point(0, 0, 1) {}
            Point(const Point& o) : x(o.x), y(o.y), z(o.z) {}

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

            Point& operator+= (const Point& p) { *this = *this + p; return *this; }
            Point& operator-= (const Point& p) { *this = *this - p; return *this; }
            Point& operator*= (T p) { *this = *this * p; return *this; }
            Point& operator/= (T p) { *this = *this / p; return *this; }

            T Len() const { return std::sqrt(LenSq()); }
            T LenSq() const { return x * x + y * y + z * z; }
        };

        template <class T> using Point3 = Point<T, 3>;

        template <class T> T Dot(const Point2<T>& v1, const Point2<T>& v2) { return v1.x * v2.x + v1.y * v2.y; }
        template <class T> T Dot(const Point3<T>& v1, const Point3<T>& v2) { return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z; }
            
        template <class T> Point3<T> Cross(const Point3<T>& v1, const Point3<T>& v2) { 
            Point3<T> res(
                v1.y * v2.z - v1.z * v2.y,
                v1.z * v2.x - v1.x * v2.z,
                v1.x * v2.y - v1.y * v2.x
            );
            return res;
        }

        template <class T, int Dim>
        Point<T, Dim> Normalize(Point<T, Dim> p) {
            p /= p.Len();
            return p;
        }

        // ************************************************************************************
        // Matrix
        // ************************************************************************************
        template <class T, int Dim>
        struct Matrix;
        
        template <class T>
        struct Matrix<T, 3> {
            using value_type = T;
            union {
                std::array<Point3<T>, 3> r;
                std::array<T, 9> data;
            };
            Matrix(T _11, T _12, T _13, T _21, T _22, T _23, T _31, T _32, T _33) : r{Point3<T>(_11, _12, _13), Point3<T>(_21, _22, _23), Point3<T>(_31, _32, _33) } {}
            Matrix(const Point3<T>& r1, const Point3<T>& r2, const Point3<T>& r3) : r{r1, r2, r3} {}
            Matrix(const T* d) : Matrix(d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7], d[8]) {}
            explicit Matrix(IdentityTag) : Matrix(Point3<T>(kXUnit), Point3<T>(kYUnit), Point3<T>(kZUnit)) {}
            Matrix(const Matrix& m) : r{ m.r } {}

            Matrix operator+(const Matrix& m) const { Matrix res(*this); res += m; return res; }
            Matrix operator*(const Matrix& m) const { Matrix res(*this); res *= m; return res; }

            Matrix& operator+=(const Matrix& m) {
                for (int i = 0; i < 9; ++i) {
                    data[i] += m.data[i];
                }
                return *this;
            }
            Matrix& operator*=(const Matrix& m) {
                Matrix tmp(*this);
                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                        data[i * 3 + j] = 0;
                        for (int k = 0; k < 3; ++k)
                            data[i * 3 + j] += tmp.r[i].data[k] * m.r[k].data[j];
                    }
                }
                return *this;
            }

            Point3<T> operator*(const Point3<T>& p) const { return { Dot(r[0], p), Dot(r[1], p), Dot(r[2], p) }; }

            Matrix operator*(const T s) const { return { r[0] * s, r[1] * s, r[2] * s }; }
        };

        template <class T>
        Matrix<T, 3> CreateCrossProdMatrix(const Point3<T>& p) {
            Matrix<T, 3> res(0, -p.z, -p.y, p.z, 0, -p.x, -p.y, p.x, 0);
            return res;
        }

        template <class T>
        Matrix<T, 3> RotationBetween(const Point3<T>& a, const Point3<T>& b, T maxDiffSq) {
            assert(AlmostEqualRelativeAndAbs(a.LenSq(), Normalize(a).LenSq(), maxDiffSq));
            assert(AlmostEqualRelativeAndAbs(b.LenSq(), Normalize(b).LenSq(), maxDiffSq));
            assert(!AlmostEqualToZero(Cross(a, b).LenSq(), maxDiffSq));
            const auto k = Cross(a, b);
            const T sin = k.Len();
            const T cos = Dot(a, b);
            Matrix<T, 3> res(kIdentity);
            Matrix<T, 3> kx = CreateCrossProdMatrix(Normalize(k));
            res += kx * sin;
            res += kx * kx * (1 - cos);
            return res;
        }

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
#pragma warning(disable : 4201)
                struct { T a, b, c; };
#pragma warning(default : 4201)
                std::array<T, 3> data;
            };

            Point<T, 2> GetNormal() const {
                Point<T, 2> normal(data.data());
                return normal;
            }

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

            Point2<T> ProjectAlong(const Point2<T>& p, Point2<T> dir, T maxDiff) const {
                const auto l = Create(dir, p);
                const auto res = Intersection(l, *this, maxDiff);
                return res;
            }

            T Distance(const Point2<T>& p, T maxDiff) const {
                if (IsVertical(maxDiff)) {
                    return std::fabs(p.x - c / a);
                }
                const T dist = std::fabs(a*p.x + b * p.y - c) / std::sqrt(a*a + b * b);
                return dist;
            }

            static HyperPlane Create(Point2<T> dir, const Point2<T>& p) {
                swap(dir.x, dir.y);
                dir.x = -dir.x;
                const T intersept = Dot(dir, p);
                return HyperPlane(dir.x, dir.y, intersept);
            }
        };

        template <class T> using Line = HyperPlane<T, 2>;

        template <class T>
        struct HyperPlane<T, 3> {
            // a*x + b*y + c*z = d
            HyperPlane(T a, T b, T c, T d) : a(a), b(b), c(c), d(d) {}
            HyperPlane(Point3<T> p, T d) : HyperPlane(p.x, p.y, p.z, d) {}

            using value_type = T;
            union {
#pragma warning(disable : 4201)
                struct { T a, b, c, d; };
#pragma warning(default : 4201)
                std::array<T, 4> data;
            };

            Point<T, 3> GetNormal() const {
                Point<T, 3> normal(data.data());
                return normal;
            }

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
            explicit HalfSpace(const HyperPlane<T, Dim>& boundary, bool needsFlip = false) : boundary(boundary) {
                if (needsFlip) {
                    flip();
                }
            }

            HyperPlane<T, Dim> boundary;

            void flip() {
                std::transform(boundary.data.cbegin(), boundary.data.cend(), boundary.data.begin(), [](const T& t) { return -t; });
            }

            Point<T, Dim> GetNormal() const {
                return boundary.GetNormal();
            }

            bool IsInside(const Point<T, Dim>& p, T maxDiff) const {
                const auto base = boundary.Project(p, maxDiff);
                const auto dir = p - base;
                const T d = Dot(dir, GetNormal());
                return d >= 0;
            }
        };

        template <class T> using HalfPlane = HalfSpace<T, 2>;

   }
}
