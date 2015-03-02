// Authors: Unknown. Please, if you are the author of this file, or if you 
// know who are the authors of this file, let us know, so we can give the 
// adequate credits and/or get the adequate authorizations.

#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <cassert>

namespace libNumerics {

// Forward declaration, definition below
template <typename T> class vector;
template <typename T> class matrix;

template <typename T> matrix<T> cat(const matrix<T>&, const matrix<T>&);
template <typename T> void swap(matrix<T>&, matrix<T>&);

/// Matrix class
template <typename T>
class matrix
{
public:
    static matrix<T> zeros(int m) { return zeros(m,m); }
    static matrix<T> zeros(int m, int n);
    static matrix<T> ones(int m) { return ones(m,m); }
    static matrix<T> ones(int m, int n);
    static matrix<T> eye(int n); ///< Identity matrix.

public:
    matrix(int m, int n);
    matrix(const matrix<T>& m);
    virtual ~matrix();
    matrix<T>& operator=(const matrix<T>& m);

    int nrow() const { return m_rows; } ///< The number of rows.
    int ncol() const { return m_cols; } ///< The number of columns.
    T  operator() (int i, int j) const;
    T& operator() (int i, int j);
    T  operator() (int i) const;
    T& operator() (int i);

    void operator=(T a);
    matrix<T> operator*(T a) const;
    matrix<T> operator/(T a) const;
    void operator*=(T a);
    void operator/=(T a);
    /// Product by scalar.
    friend matrix<T> operator*(T a, const matrix<T>& m)
    { return m * a; }

    matrix<T> operator+(const matrix<T>& m) const;
    matrix<T> operator-(const matrix<T>& m) const;
    matrix<T> operator-() const; ///< Matrix opposite.
    matrix<T> operator*(const matrix<T>& m) const;
    vector<T> operator*(const vector<T>& m) const;

    void operator+=(const matrix<T>& m);
    void operator-=(const matrix<T>& m);

    matrix<T> t() const; ///< Transpose.
    vector<T> diag() const; ///< Diagonal of matrix.
    T tr() const;
    T det() const;
    matrix<T> inv() const;

    void symUpper();
    void symLower();

    matrix<T> copy(int i0, int i1, int j0, int j1) const;
    matrix<T> copyCols(int j0, int j1) const;
    matrix<T> copyRows(int i0, int i1) const;
    void paste(int i0, int j0, const matrix<T>& block);
    friend matrix<T> cat<T>(const matrix<T>& left, const matrix<T>& right);
    vector<T> col(int j) const; ///< Copy column.
    matrix<T> row(int i) const; ///< Copy row.
    int lastCol() const {return m_cols-1;} ///< Index of last column.
    int lastRow() const {return m_rows-1;} ///< Index of last row.

    friend void swap<T>(matrix<T>&, matrix<T>&);
    void swapRows(int i0, int i1);
    void swapCols(int j0, int j1);

    template <typename U>
    void read(const U* v);
    void read(const matrix<T>& v);
    void write(T* vect) const;

protected:
    int m_rows; ///< Number of rows.
    int m_cols; ///< Number of columns.
    T* p; ///< 1-D array of coefficients.

    void alloc(int m, int n); ///< Allocate the array value.
    void free(); ///< Free the array value.
    int nElements() const; ///< Number of elements in the matrix.
    matrix<T>& sub(matrix<T>& s, int i, int j) const;
}; // class matrix

/// Column vector class (template)
template <typename T>
class vector : public matrix<T>
{
public:
    explicit vector(int m);
    vector(T x);
    vector(T x, T y);
    vector(T x, T y, T z);
    vector(const vector<T>& v);
    virtual ~vector() {}
    using matrix<T>::operator=;
    vector<T>& operator=(const vector<T>& v);
    //    void operator=(T a);

    vector<T> operator*(T a) const;
    vector<T> operator/(T a) const;
    /// Product of a vector by a scalar.
    friend vector<T> operator*(T a, const vector<T>& v)
    { return v * a; }

    vector<T> operator+(const vector<T>& v) const;
    vector<T> operator-(const vector<T>& v) const;
    vector<T> operator-() const; ///< Vector opposite.

    matrix<T> operator*(const matrix<T>& m) const;
    matrix<T> diag() const;
    T qnorm() const;

    vector<T> copy(int i0, int i1) const;
    void paste(int i0, const vector<T>& v);
};

} // namespace libNumerics

/// Output matrix coefficients.
template <typename T>
inline std::ostream& operator<<(std::ostream& out,
                                const libNumerics::matrix<T>& m)
{
    for(int i = 0; i < m.nrow(); ++i) {
        out << ((i==0)? "[": ";");
        for (int j = 0; j < m.ncol(); ++j)
            out << " " << m(i,j);
    }
    out << " ]";
    return out;
}

/// Input matrix. Need to know the dimensions in advance...
template <class T>
inline std::istream& operator>>(std::istream& in,
                                libNumerics::matrix<T>& m)
{
    char c;
    for(int i=0; i < m.nrow(); ++i) {
        in >> c;
        for(int j=0; j < m.ncol(); ++j)
            in >> m(i,j);
    }
    in >> c;
    return in;
}

template <typename T>
T dot(const libNumerics::vector<T>& u, const libNumerics::vector<T>& v);
template <typename T>
libNumerics::vector<T> cross(const libNumerics::vector<T>& u,
                             const libNumerics::vector<T>& v);

// Need to see definitions for templates...
#include "matrix.cpp"
#include "vector.cpp"

#endif
