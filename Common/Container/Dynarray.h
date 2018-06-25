#pragma once
// Implementation base on proposition for C++ Dynamic Arrays
// http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2013/n3662.html
#include <memory>

namespace Common
{
    template<class T>
    class Dynarray
    {
    public:
        using value_type = T;
        using reference = T&;
        using const_reference = const T&;
        using pointer = T*;
        using const_pointer = const T*;
        using iterator = pointer;
        using const_iterator = const_pointer;
        using reverse_iterator = std::reverse_iterator<iterator>;
        using const_reverse_iterator = std::reverse_iterator<const_iterator>;
        using size_type = size_t;
        using difference_type = ptrdiff_t;

        // Can be used only with types which have default constructor
        Dynarray() : data_(nullptr), size_(0) {}
        explicit Dynarray(size_type s) : data_(std::make_unique<T[]>(s)), size_(s) {}
        Dynarray(size_type s, const T& v) : data_(std::make_unique<T[]>(s)), size_(s) { fill(v); }
        Dynarray(const Dynarray& d) : data_(std::make_unique<T[]>(d.size())), size_(d.size()) { memcpy(data_.get(), d.data_.get(), size_ * sizeof(T)); }
        Dynarray(Dynarray&& d) : data_(move(d.data_)), size_(d.size()) { d.size_ = 0; }
        Dynarray(std::initializer_list<T> il) : data_(make_unique<T[]>(il.size())), size_(il.size()) { for (size_type i = 0; i < size_; ++i) data_[i] = *(std::cbegin(il) + i); }
        Dynarray(std::unique_ptr<T[]> data, size_type size) : data_(move(data)), size_(size) {}
        Dynarray& operator=(const Dynarray& d) {
            if (size_ != d.size_) { size_ = d.size_; data_ = make_unique<T[]>(size_); }
            memcpy(data_.get(), d.data_.get(), size_ * sizeof(T));
            return *this;
        }
        Dynarray& operator=(Dynarray&& d) { data_ = move(d.data_); size_ = d.size_; d.size_ = 0; return *this; }

#define CHECK_NOT_EMPTY assert(size_ > 0)
        iterator       begin()        noexcept { CHECK_NOT_EMPTY; return &data_[0]; }
        const_iterator begin()  const noexcept { CHECK_NOT_EMPTY; return &data_[0]; }
        const_iterator cbegin() const noexcept { CHECK_NOT_EMPTY; return &data_[0]; }
        iterator       end()          noexcept { CHECK_NOT_EMPTY; return &data_[size_]; }
        const_iterator end()    const noexcept { CHECK_NOT_EMPTY; return &data_[size_]; }
        const_iterator cend()   const noexcept { CHECK_NOT_EMPTY; return &data_[size_]; }

        reverse_iterator       rbegin()        noexcept { return reverse_iterator(end()); }
        const_reverse_iterator rbegin()  const noexcept { return reverse_iterator(end()); }
        const_reverse_iterator crbegin() const noexcept { return reverse_iterator(cend()); }
        reverse_iterator       rend()          noexcept { return reverse_iterator(begin()); }
        const_reverse_iterator rend()    const noexcept { return reverse_iterator(begin()); }
        const_reverse_iterator crend()   const noexcept { return reverse_iterator(cbegin()); }

        size_type size()     const noexcept { return size_; }

        reference       operator[](size_type n) { assert(n < size_); return data_[n]; }
        const_reference operator[](size_type n) const { assert(n < size_); return data_[n]; }

        reference       front()       { CHECK_NOT_EMPTY; return data_[0]; }
        const_reference front() const { CHECK_NOT_EMPTY; return data_[0]; }
        reference       back()        { CHECK_NOT_EMPTY; return data_[size_ - 1]; }
        const_reference back()  const { CHECK_NOT_EMPTY; return data_[size_ - 1]; }

        const_reference at(size_type n) const { if (size_ <= n) throw std::out_of_range("Index is out of range"); return data_[n]; }
        reference       at(size_type n)       { if (size_ <= n) throw std::out_of_range("Index is out of range"); return data_[n]; }
#undef CHECK_NOT_EMPTY

        T*       data()       noexcept { return data_.get(); }
        const T* data() const noexcept { return data_.get(); }

        void fill(const T& v) { for (size_type i = 0; i < size_; ++i) data_[i] = v; }

    private:
        std::unique_ptr<T[]> data_;
        size_type size_;
    };

}
