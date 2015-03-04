// Authors: Unknown. Please, if you are the author of this file, or if you 
// know who are the authors of this file, let us know, so we can give the 
// adequate credits and/or get the adequate authorizations.

#ifdef MATRIX_H // Do nothing if not included from matrix.h

#define INDEX(i,j) ((i) * m_cols + (j))

namespace libNumerics {

/// Constructor for \a m*\a n matrix.
/// \param m number of rows.
/// \param n number of columns.
template <typename T>
matrix<T>::matrix(int m, int n)
{
    alloc(m, n);
}

/// Copy constructor.
template <typename T>
matrix<T>::matrix(const matrix<T>& m)
{
    alloc(m.m_rows, m.m_cols);
    for(int i = nElements()-1; i >= 0; i--)
        p[i] = m.p[i];
}

/// Destructor.
template <typename T>
matrix<T>::~matrix()
{
    free();
}

/// Assignment operator.
template <typename T>
matrix<T>& matrix<T>::operator=(const matrix<T>& m)
{
    if(&m == this) return *this;
    if(m.nElements() != nElements()){
        free();
        alloc(m.m_rows, m.m_cols);
    } else {
        m_rows = m.m_rows;
        m_cols = m.m_cols;
    }
    for(int i = nElements()-1; i >= 0; i--)
        p[i] = m.p[i];
    return *this;
}

/// Access the coefficient on the \a i-th row, \a j-th column.
template <typename T>
inline T matrix<T>::operator() (int i, int j) const
{
    assert(i >= 0 && i < m_rows && j >= 0 && j < m_cols);
    return p[INDEX(i,j)];
}

/// Access the coefficient on the \a i-th row, \a j-th column.
template <typename T>
inline T& matrix<T>::operator() (int i, int j)
{
    assert(i >= 0 && i < m_rows && j >= 0 && j < m_cols);
    return p[INDEX(i,j)];
}

template <typename T>
inline T matrix<T>::operator() (int i) const
{
    assert(i >= 0 && i < nElements());
    return p[i];
}

template <typename T>
inline T& matrix<T>::operator() (int i)
{
    assert(i >= 0 && i < nElements());
    return p[i];
}

/// Set matrix at constant value.
///
/// Assign all coefficients to the value \a a.
template <typename T>
inline void matrix<T>::operator=(T a)
{
    for(int i = nElements()-1; i >= 0; i--)
        p[i] = a;
}

/// Multiply a matrix by scalar.
/// \param a a scalar.
template <typename T>
matrix<T> matrix<T>::operator*(T a) const
{
    matrix<T> prod(m_rows, m_cols);
    for(int i = nElements()-1; i >= 0; i--)
        prod.p[i] = a * p[i];
    return prod;
}

/// Multiply a matrix by scalar.
/// \param a a scalar.
template <typename T>
void matrix<T>::operator*=(T a)
{
    for(int i = nElements()-1; i >= 0; i--)
        p[i] *= a;
}

/// Divide a matrix by scalar.
/// \param a a scalar.
template <typename T>
matrix<T> matrix<T>::operator/(T a) const
{
    return (*this) * ((T)1/a);
}

/// Divide a matrix by scalar.
/// \param a a scalar.
template <typename T>
void matrix<T>::operator/=(T a)
{
    *this *= (T)1 / a;
}

/// Matrix sum.
template <typename T>
matrix<T> matrix<T>::operator+(const matrix<T>& m) const
{
    assert(m.m_rows == m_rows && m.m_cols == m_cols);
    matrix<T> sum(m_rows,m_cols);
    for(int i = nElements()-1; i >= 0; i--)
        sum.p[i] = p[i] + m.p[i];
    return sum;
}

/// Matrix sum.
template <typename T>
void matrix<T>::operator+=(const matrix<T>& m)
{
    assert(m.m_rows == m_rows && m.m_cols == m_cols);
    for(int i = nElements()-1; i >= 0; i--)
        p[i] += m.p[i];
}

/// Matrix subtraction.
template <typename T>
matrix<T> matrix<T>::operator-(const matrix<T>& m) const
{
    assert(m.m_rows == m_rows && m.m_cols == m_cols);
    matrix<T> sub(m_rows,m_cols);
    for(int i = nElements()-1; i >= 0; i--)
        sub.p[i] = p[i] - m.p[i];
    return sub;
}

/// Matrix subtraction.
template <typename T>
void matrix<T>::operator-=(const matrix<T>& m)
{
    assert(m.m_rows == m_rows && m.m_cols == m_cols);
    for(int i = nElements()-1; i >= 0; i--)
        p[i] -= m.p[i];
}

template <typename T>
matrix<T> matrix<T>::operator-() const
{
    matrix<T> opp(m_rows, m_cols);
    for(int i = nElements()-1; i >= 0; i--)
        opp.p[i] = -p[i];
    return opp;
}

/// Matrix multiplication.
template <typename T>
matrix<T> matrix<T>::operator*(const matrix<T>& m) const
{
    assert(m_cols == m.m_rows);
    matrix<T> prod(m_rows, m.m_cols);
    T* out = prod.p;
    for(int i = 0; i < prod.m_rows; i++) {
        const T* left = p + i*m_cols;
        for(int j = 0; j < prod.m_cols; j++, out++) {
            const T* right = m.p + j;
            *out = 0;
            for(int k = 0; k < m_cols; k++) {
                *out += left[k] * *right;
                right += m.m_cols;
            }
        }
    }
    return prod;
}

/// Matrix-vector multiplication.
template <typename T>
vector<T> matrix<T>::operator*(const vector<T>& m) const
{
    assert(m_cols == m.m_rows);
    vector<T> prod(m_rows);
    T* out = prod.p;
    for(int i = 0; i < prod.m_rows; i++, out++) {
        const T* left = p + i*m_cols;
        const T* right = m.p;
        *out = 0;
        for(int k = 0; k < m_cols; k++)
            *out += left[k] * right[k];
    }
    return prod;
}

/// Tranposed of matrix.
template <typename T>
matrix<T> matrix<T>::t() const
{
    matrix<T> t(ncol(), nrow());
    T* out = t.p;
    for(int i = 0; i < t.nrow(); i++) {
        const T* in = p + i;
        for(int j = 0; j < t.ncol(); j++) {
            *out++ = *in;
            in += ncol();
        }
    }
    return t;
}

/// Symmetrize upper part of matrix.
template <typename T>
void matrix<T>::symUpper()
{
    assert(m_rows == m_cols);
    for(int i = 1; i < m_rows; i++) {
        const T* in = p + i;
        T* out = p + m_cols*i;
        for(int j = 0; j < i; j++) {
            *out++ = *in;
            in += m_cols;
        }
    }
}

/// Symmetrize lower part of matrix.
template <typename T>
void matrix<T>::symLower()
{
    assert(m_rows == m_cols);
    for(int i = 1; i < m_rows; i++) {
        const T* in = p + m_cols*i;
        T* out = p + i;
        for(int j = 0; j < i; j++) {
            *out = *in++;
            out += m_cols;
        }
    }
}

template <typename T>
vector<T> matrix<T>::diag() const
{
    assert(m_rows == m_cols);
    vector<T> t(m_rows);
    for(int i = 0; i < m_rows; i++)
        t.p[i] = p[i*(m_cols+1)];
    return t;
}

/// Matrix made of zeros.
template <typename T>
matrix<T> matrix<T>::zeros(int m, int n)
{
    matrix<T> M(m,n);
    for(int i = M.nElements()-1; i >= 0; i--)
        M.p[i] = (T)0;
    return M;
}

/// Matrix made of ones.
template <typename T>
matrix<T> matrix<T>::ones(int m, int n)
{
    matrix<T> M(m,n);
    for(int i = M.nElements()-1; i >= 0; i--)
        M.p[i] = (T)1;
    return M;
}

/// Identity matrix.
template <typename T>
matrix<T> matrix<T>::eye(int n)
{
    matrix<T> M(n,n);
    for(int i = M.nElements()-1; i >= 0; i--)
        M.p[i] = (T)0;
    for(int i = n-1; i >= 0; i--)
        M.p[i*(n+1)] = (T)1;
    return M;
}

/// Extract the submatrix [i0,i1]x[j0,j1].
/// \param i0 first row
/// \param i1 last row
/// \param j0 first column
/// \param j1 last column
template <typename T>
matrix<T> matrix<T>::copy(int i0, int i1, int j0, int j1) const 
{
    assert(0 <= i0 && i0 <= i1 && i1 <= m_rows &&
           0 <= j0 && j0 <= j1 && j1 <= m_cols);
    matrix<T> sub(i1-i0+1,j1-j0+1);
    T* out = sub.p;
    for(int i = i0; i <= i1; i++) {
        const T* in = p + INDEX(i, j0);
        for(int j = j0; j <= j1; j++)
            *out++ = *in++;
    }
    return sub;
}

/// Extract the columns of index in [j0,j1].
/// \param j0 first column
/// \param j1 last column
template <typename T>
matrix<T> matrix<T>::copyCols(int j0, int j1) const 
{
    return copy(0, lastRow(), j0, j1);
}

/// Extract the rows of index in [i0,i1].
/// \param i0 first row
/// \param i1 last row
template <typename T>
matrix<T> matrix<T>::copyRows(int i0, int i1) const 
{
    return copy(i0, i1, 0, lastCol());
}

/// Paste a matrix in another one, at position (\a i0,\a j0)
/// \param i0 first row where to paste in
/// \param j0 first column where to paste in
/// \param matrix to paste
template <typename T>
void matrix<T>::paste(int i0, int j0, const matrix<T>& m)
{
    assert(i0 >= 0 && i0+m.m_rows <= m_rows &&
           j0 >= 0 && j0+m.m_cols <= m_cols);
    const T* in = m.p;
    for(int i = 0; i < m.m_rows; i++) {
        T* out = p + INDEX(i0+i, j0);
        for(int j = 0; j < m.m_cols; j++)
            *out++ = *in++;
    }
}

/// Concatenate matrices.
template <typename T>
matrix<T> cat(const matrix<T>& m1, const matrix<T>& m2)
{
    assert(m1.m_rows == m2.m_rows);
    matrix<T> m(m1.m_rows, m1.m_cols+m2.m_cols);
    m.paste(0, 0, m1);
    m.paste(0, m1.m_cols, m2);
    return m;
}

/// Copy column number \a j.
template <typename T>
vector<T> matrix<T>::col(int j) const
{
    assert(j >= 0 && j < m_cols);
    vector<T> c(m_rows);
    const T* in = p + j;
    for(int i = 0; i < m_rows; i++) {
        c(i) = *in;
        in += m_cols;
    }
    return c;
}

/// Copy row number \a i.
template <typename T>
inline matrix<T> matrix<T>::row(int i) const
{
    return copy(i, i, 0, lastCol());
}

template <class T>
void swap(matrix<T>& A, matrix<T>& B)
{
    int i=A.m_rows;
    A.m_rows = B.m_rows;
    B.m_rows = i;
    i = A.m_cols;
    A.m_cols = B.m_cols;
    B.m_cols = i;
    T* p = A.p;
    A.p = B.p;
    B.p = p;
}

template <typename T>
void matrix<T>::swapRows(int i0, int i1)
{
    assert(0 <= i0 && i0 < m_rows &&
           0 <= i1 && i1 < m_rows);
    T* row0 = p + i0*m_cols;
    T* row1 = p + i1*m_cols;
    for(int j = m_cols-1; j >= 0; j--) {
        T tmp = *row0; *row0++ = *row1; *row1++ = tmp;
    }
}

template <typename T>
void matrix<T>::swapCols(int j0, int j1)
{
    assert(0 <= j0 && j0 < m_cols &&
           0 <= j1 && j1 < m_cols);
    T* col0 = p + j0;
    T* col1 = p + j1;
    for(int i = m_rows-1; i >= 0; i--) {
        T tmp = *col0; *col0 = *col1; *col1 = tmp;
        col0 += m_cols;
        col1 += m_cols;
    }
}

/// Copy the array values in a matrix, row by row.
/// \param m number of rows
/// \param n number of columns
/// \param v an array of scalar of size m*n
template <typename T> template <typename U> 
void matrix<T>::read(const U* v)
{
    for(int i = nElements()-1; i >= 0; i--)
        p[i] = (T)v[i];
}

/// Read the coefficients from \a m.
template <typename T> 
inline void matrix<T>::read(const matrix<T>& m)
{
    assert(m.nElements() == nElements());
    read(m.p);
}

/// Copy the matrix coefficients in an array.
///
/// The matrix is scanned row by row. 
template <typename T>
void matrix<T>::write(T* vect) const
{
    for(int i = nElements()-1; i >= 0; i--)
        vect[i] = p[i];
}

template <typename T>
void matrix<T>::alloc(int m, int n)
{
    assert(m > 0 && n > 0);  
    m_rows = m;
    m_cols = n;
    p = new T[m*n];
}

template <typename T>
inline void matrix<T>::free()
{
    delete [] p;
    p = NULL;
}

template <typename T>
inline int matrix<T>::nElements() const
{
    return m_rows*m_cols;
}

/// Submatrix without row \a i0 and col \a j0.
template <typename T>
matrix<T>& matrix<T>::sub(matrix<T>& s, int i0, int j0) const
{
    const T* in = p;
    T* out = s.p;
    for(int i = 0; i < i0; i++) {
        for(int j = 0; j < j0; j++)
            *out++ = *in++;
        ++in; // Skip col j0
        for(int j = j0+1; j < m_cols; j++)
            *out++ = *in++;
    }
    in += m_cols; // Skip row i0
    for(int i = i0+1; i < m_rows; i++) {
        for(int j = 0; j < j0; j++)
            *out++ = *in++;
        ++in; // Skip col j0
        for(int j = j0+1; j < m_cols; j++)
            *out++ = *in++;
    }
    return s;
}

/// Trace.
template <typename T>
T matrix<T>::tr() const
{
    assert(m_rows == m_cols);
    T res = (T)0;
    for(int i = 0; i < m_rows; i++)
        res += p[i*(m_cols+1)];
    return res;
}

/// Determinant. Slow, use only for small matrices.
template <typename T>
T matrix<T>::det() const
{
    assert(m_rows == m_cols);
    if(m_rows == 1)
        return p[0];
    if(m_rows == 2)
        return (p[0]*p[3]-p[1]*p[2]);
    T res = (T)0;
    T sign = (T)1;
    matrix<T> s(m_rows-1, m_cols-1);
    for(int j = 0; j < m_cols; j++) {
        res += sign*p[j]*sub(s,0,j).det();
        sign = -sign;
    }
    return res;
}

/// Inverse. Slow, use only for small matrices.
template <typename T>
matrix<T> matrix<T>::inv() const
{
    assert(m_rows == m_cols);
    matrix<T> res(m_rows, m_cols);
    if(m_rows == 1)
        res.p[0] = (T)1/p[0];
    else {
        T d = (T)1 / det();
        T signi = (T)1;
        T* out = res.p;
        matrix<T> s(m_rows-1, m_cols-1);
        for(int i = 0; i < m_rows; i++) {
            T signj = signi;
            for(int j = 0; j < m_cols; j++) {
                *out++ = signj*d*sub(s,j,i).det();
                signj = -signj;
            }
            signi = -signi;
        }
    }
    return res;
}

} // namespace libNumerics

#undef INDEX

#endif // MATRIX_H
