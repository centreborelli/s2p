// Authors: Unknown. Please, if you are the author of this file, or if you 
// know who are the authors of this file, let us know, so we can give the 
// adequate credits and/or get the adequate authorizations.

#include "numerics.h"
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>

namespace libNumerics {

const flnum MinLM::DEFAULT_RELATIVE_TOL = 1E-3;
const flnum MinLM::DEFAULT_LAMBDA_INIT  = 1E-3;
const flnum MinLM::DEFAULT_LAMBDA_FACT  = 10.0;
const flnum MinLM::EPSILON_KERNEL       = 1E-9;

inline flnum ABS(flnum x)
{ return (x >= 0)? x: -x; }

/// Resolution by LU decomposition with pivot.
bool solveLU(const matrix<flnum>& A, const vector<flnum>& B, vector<flnum>& X)
{
    X = B;
    return solveLU(A, X);
}

/// Replace X by A^{-1}X, by LU solver.
bool solveLU(matrix<flnum> A, vector<flnum>& X)
{
    assert(A.nrow() == A.ncol());
    int	n = A.nrow();
    vector<flnum> rowscale(n); // Implicit scaling of each row
    std::vector<int> permut(n,0); // Permutation of rows

    // Get the implicit scaling information of each row
    for(int i=0; i< n; i++) {
        flnum max = 0.0;
        for(int j=0; j< n; j++) {
            flnum tmp = ABS(A(i,j));
            if (tmp> max)
                max = tmp;
        }
        if(max == 0.0)
            return false;
        rowscale(i) = 1.0/max;
    }
    // Perform the decomposition
    for(int k=0; k < n; k++) {
        // Search for largest pivot element
        flnum max = rowscale(k)*ABS(A(k,k));
        int imax = k;
        for(int i=k+1; i < n; i++) {
            flnum tmp = rowscale(i)*ABS(A(i,k));
            if(tmp > max) {
                max = tmp;
                imax = i;
            }
        }
        if(max == 0.0)
            return false;

        // Interchange rows if needed
        if(k != imax) {
            A.swapRows(k, imax);
            rowscale(imax) = rowscale(k); // Scale of row k no longer needed
        }
        permut[k] = imax; // permut(k) was not initialized before
        flnum Akk = 1/A(k,k);
        for(int i=k+1; i < n; i++) {
            flnum tmp = A(i,k) *= Akk; // Divide by pivot
            for (int j=k+1;j < n; j++) // Reduce the row
                A(i,j) -= tmp*A(k,j);
        }
    }
    // Forward substitution
    for (int k=0; k < n; k++) {
        flnum sum = X(permut[k]);
        X(permut[k]) = X(k);
        for(int j = 0; j < k; j++)
            sum -= A(k,j)*X(j);
        X(k) = sum;
    }
    // Backward substitution
    for(int k=n-1; k >= 0; k--) {
        flnum sum = X(k);
        for(int j=k+1; j < n; j++)
            sum -= A(k,j)*X(j);
        X(k) = sum/A(k,k);
    }
    return true;
}

/// Decompose A into U diag(W) V^T with U(m,n) and V(n,n) having orthonormal
/// vectors.
SVD::SVD(const matrix<flnum>& A)
: m_U(A), m_V(A.ncol(),A.ncol()), m_W(A.ncol())
{
    compute();
    sort();
}

/// SVD computation. Initial matrix stored in m_U as input.
void SVD::compute()
{
    const flnum	EPSILON = std::numeric_limits<flnum>::epsilon();
    const int SVD_MAX_ITS = 30;

    int rows = m_U.nrow();
    int cols = m_U.ncol();
    flnum g, scale, anorm;
    vector<flnum> RV1(cols);

    // Householder reduction to bidiagonal form:
    anorm = g = scale = 0.0;
    for (int i=0; i< cols; i++) {
        int l = i + 1;
        RV1(i) = scale*g;
        g = scale = 0.0;
        if(i< rows) {
            for (int k=i; k< rows; k++)
                scale += ABS(m_U(k,i));
            if (scale != 0.0) {
                flnum invScale=1.0/scale, s=0.0;
                for (int k=i; k< rows; k++) {
                    m_U(k,i) *= invScale;
                    s += m_U(k,i) * m_U(k,i);
                }
                flnum f = m_U(i,i);
                g = - withSignOf(std::sqrt(s),f);
                flnum h = 1.0 / (f*g - s);
                m_U(i,i) = f - g;
                for (int j=l; j< cols; j++) {
                    s = 0.0;
                    for (int k=i; k< rows; k++)
                        s += m_U(k,i) * m_U(k,j);
                    f = s * h;
                    for (int k=i; k< rows; k++)
                        m_U(k,j) += f * m_U(k,i);
                }
                for (int k=i; k< rows; k++)
                    m_U(k,i) *= scale;
            }
        }
        m_W(i) = scale * g;
        g = scale = 0.0;
        if ( i< rows && i< cols-1 ) {
            for (int k=l; k< cols; k++)
                scale += ABS(m_U(i,k));
            if (scale != 0.0) {
                flnum invScale=1.0/scale, s=0.0;
                for (int k=l; k< cols; k++) {
                    m_U(i,k) *= invScale;
                    s += m_U(i,k) * m_U(i,k);
                }
                flnum f = m_U(i,l);
                g = - withSignOf(std::sqrt(s),f);
                flnum h = 1.0 / (f*g - s);
                m_U(i,l) = f - g;
                for (int k=l; k< cols; k++)
                    RV1(k) = m_U(i,k) * h;
                for (int j=l; j< rows; j++) {
                    s = 0.0;
                    for (int k=l; k< cols; k++)
                        s += m_U(j,k) * m_U(i,k);
                    for (int k=l; k< cols; k++)
                        m_U(j,k) += s * RV1(k);
                }
                for (int k=l; k< cols; k++)
                    m_U(i,k) *= scale;
            }
        }
        anorm = std::max(anorm, ABS(m_W(i)) + ABS(RV1(i)) );
    }

    // Accumulation of right-hand transformations:
    m_V(cols-1,cols-1) = 1.0;
    for (int i= cols-2; i>=0; i--) {
        m_V(i,i) = 1.0;
        int l = i+1;
        g = RV1(l);
        if (g != 0.0) {
            flnum invgUil = 1.0 / (m_U(i,l)*g);
            for (int j=l; j< cols; j++)
                m_V(j,i) = m_U(i,j) * invgUil;
            for (int j=l; j< cols; j++){
                flnum s = 0.0;
                for (int k=l; k< cols; k++)
                    s += m_U(i,k) * m_V(k,j);
                for (int k=l; k< cols; k++)
                    m_V(k,j) += s * m_V(k,i);
            }
        }
        for (int j=l; j< cols; j++)
            m_V(i,j) = m_V(j,i) = 0.0;
    }

    // Accumulation of left-hand transformations:
    for (int i=std::min(rows,cols)-1; i>=0; i--) {
        int l = i+1;
        g = m_W(i);
        for (int j=l; j< cols; j++)
            m_U(i,j) = 0.0;
        if (g != 0.0) {
            g = 1.0 / g;
            flnum invUii = 1.0 / m_U(i,i);
            for (int j=l; j< cols; j++) {
                flnum s = 0.0;
                for (int k=l; k< rows; k++)
                    s += m_U(k,i) * m_U(k,j);
                flnum f = (s * invUii) * g;
                for (int k=i; k< rows; k++)
                    m_U(k,j) += f * m_U(k,i);
            }
            for (int j=i; j< rows; j++)
                m_U(j,i) *= g;
        } else
            for (int j=i; j< rows; j++)
                m_U(j,i) = 0.0;
        m_U(i,i) = m_U(i,i) + 1.0;
    }

    // Diagonalization of the bidiagonal form:
    for (int k=cols-1; k>=0; k--) { // Loop over singular values
        for (int its=1; its<=SVD_MAX_ITS; its++) {
            bool flag = false;
            int l  = k;
            int nm = k-1;
            while(l>0 && ABS(RV1(l)) > EPSILON*anorm) { // Test for splitting
                if(ABS(m_W(nm)) <= EPSILON*anorm) {
                    flag = true;
                    break;
                }
                l--;
                nm--;
            }
            if (flag) {	// Cancellation of RV1(l), if l > 0
                flnum c=0.0, s=1.0;
                for (int i=l; i< k+1; i++) {
                    flnum f = s * RV1(i);
                    RV1(i) = c * RV1(i);
                    if (ABS(f)<=EPSILON*anorm)
                        break;
                    g = m_W(i);
                    flnum h = SVD::hypot(f,g);
                    m_W(i) = h;
                    h = 1.0 / h;
                    c = g * h;
                    s = - f * h;
                    for (int j=0; j< rows; j++)
                        rotate(m_U(j,nm),m_U(j,i), c,s); 
                }
            }
            flnum z = m_W(k);
            if (l==k) {		// Convergence of the singular value
                if (z< 0.0) {	// Singular value is made nonnegative
                    m_W(k) = -z;
                    for (int j=0; j< cols; j++)
                        m_V(j,k) = - m_V(j,k);
                }
                break;
            }
            // Exception if convergence to the singular value not reached:
            if(its==SVD_MAX_ITS) throw SvdConvergenceError();
            flnum x = m_W(l); // Get QR shift value from bottom 2x2 minor
            nm = k-1;
            flnum y = m_W(nm);
            g = RV1(nm);
            flnum h = RV1(k);
            flnum f = ( (y-z)*(y+z) + (g-h)*(g+h) ) / ( 2.0*h*y );
            g = SVD::hypot(f,1.0);
            f = ( (x-z)*(x+z) + h*(y/(f+withSignOf(g,f)) - h) ) / x;
            // Next QR transformation (through Givens reflections)
            flnum c=1.0, s=1.0;
            for (int j=l; j<=nm; j++) {
                int i = j+1;
                g = RV1(i);
                y = m_W(i);
                h = s * g;
                g = c * g;
                z = SVD::hypot(f,h);
                RV1(j) = z;
                z = 1.0 / z;
                c = f * z;
                s = h * z;
                f = x*c + g*s;
                g = g*c - x*s;
                h = y * s;
                y *= c;
                for(int jj=0; jj < cols; jj++)
                    rotate(m_V(jj,j),m_V(jj,i), c,s);
                z = SVD::hypot(f,h);
                m_W(j) = z;
                if (z!=0.0) { // Rotation can be arbitrary if z = 0.0
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = c*g + s*y;
                x = c*y - s*g;
                for(int jj=0; jj < rows; jj++)
                    rotate(m_U(jj,j),m_U(jj,i), c,s);
            }
            RV1(l) = 0.0;
            RV1(k) = f;
            m_W(k) = x;
        }
    }
}

/// Recompose from SVD. This should be the initial matrix.
matrix<flnum> SVD::compose() const
{
    return m_U * m_W.diag() * m_V.t();
}

flnum SVD::withSignOf(flnum a, flnum b)
{ return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a); }

/// Replace hypot of math.h by robust numeric implementation.
flnum SVD::hypot(flnum a, flnum b)
{
    a = ABS(a);
    b = ABS(b);
    if(a > b) {
        b /= a;
        return a*std::sqrt(1.0 + b*b);
    } else if(b) {
        a /= b;
        return b*std::sqrt(1.0 + a*a);
    }
    return 0.0;
}

/// Utility function used while computing SVD.
void SVD::rotate(flnum& a, flnum& b, flnum c, flnum s)
{
    flnum d = a;
    a = +d*c +b*s;
    b = -d*s +b*c;
}

class SVDElement {
public:
    SVDElement(const vector<flnum>& W, int i)
    : m_val(W(i)), m_i(i) {}
    bool operator<(const SVDElement& e) const
    { return (m_val>e.m_val); }

    flnum m_val;
    int m_i;
};

/// Sort SVD by decreasing order of singular value.
void SVD::sort()
{
    std::vector<SVDElement> vec;
    for(int i=0; i < m_U.ncol(); i++)
        vec.push_back( SVDElement(m_W, i) );
    std::sort(vec.begin(), vec.end());
    // Apply permutation
    for(int i=m_U.ncol()-1; i >=0; i--)
        if(vec[i].m_i != i) { // Find cycle of i
            const vector<flnum> colU = m_U.col(i);
            const vector<flnum> colV = m_V.col(i);
            const flnum w = m_W(i);
            int j = i;
            while(vec[j].m_i != i) {
                m_U.paste(0,j, m_U.col(vec[j].m_i));
                m_V.paste(0,j, m_V.col(vec[j].m_i));
                m_W(j) = m_W(vec[j].m_i);
                std::swap(j,vec[j].m_i);
            }
            vec[j].m_i = j;
            m_U.paste(0,j, colU);
            m_V.paste(0,j, colV);
            m_W(j) = w;
        }
}

/// Constructor.
MinLM::MinLM()
: iterations(0), relativeTol(DEFAULT_RELATIVE_TOL),
  lambdaInit(DEFAULT_LAMBDA_INIT), lambdaFact(DEFAULT_LAMBDA_FACT)
{}

/// In equation JtJ X = B, remove columns of J close to 0, so that JtJ can be
/// invertible
void MinLM::compress(matrix<flnum>& JtJ, vector<flnum>& B)
{
    flnum max=0;
    for(int i=0; i < JtJ.nrow(); i++)
        if(JtJ(i,i) > max)
            max = JtJ(i,i);
    max *= EPSILON_KERNEL;
    m_nullCols.clear();
    for(int i=0; i < JtJ.nrow(); i++)
        if(JtJ(i,i) <= max)
            m_nullCols.push_back(i);
    if( m_nullCols.empty() )
        return;
    int n=(int)m_nullCols.size();
    matrix<flnum> JtJ2(JtJ.nrow()-m_nullCols.size(),
                       JtJ.ncol()-m_nullCols.size());
    vector<flnum> B2(B.nrow()-(int)m_nullCols.size());
    for(int i=0,i2=0; i < JtJ.nrow(); i++) {
        if(i-i2 < n && m_nullCols[i-i2]==i)
            continue;
        for(int j=0,j2=0; j < JtJ.ncol(); j++) {
            if(j-j2 < n && m_nullCols[j-j2]==j)
                continue;
            JtJ2(i2,j2) = JtJ(i,j);
            j2++;
        }
        B2(i2) = B(i);
        i2++;
    }
    swap(JtJ,JtJ2);
    swap(B,B2);
}

/// Insert 0 in rows of B that were removed by \c compress()
void MinLM::uncompress(vector<flnum>& B)
{
    if(m_nullCols.empty())
        return;
    int n=(int)m_nullCols.size();
    vector<flnum> B2(B.nrow()+(int)m_nullCols.size());
    for(int i=0,i2=0; i2 < B2.nrow(); i2++)
        if(i2-i < n && m_nullCols[i2-i]==i2)
            B2(i2)=0;
        else
            B2(i2) = B(i++);
    swap(B,B2);
}

/// Perform minimization.
/// \a targetRMSE is the root mean square error aimed at.
/// Return the reached RMSE. Since the class does not know the dimension, the
/// real RMSE should be this value multiplied by sqrt(dim). For example, for 2-D
/// points this would be sqrt(2) times the returned value.
flnum MinLM::minimize(vector<flnum>& P, const vector<flnum>& yData,
                      flnum targetRMSE, int maxIters)
{
    flnum errorMax = targetRMSE*targetRMSE*yData.nrow();
    vector<flnum> yModel( yData.nrow() );
    modelData(P, yModel);
    vector<flnum> E( yData-yModel );
    flnum error = E.qnorm();
    matrix<flnum> J( yData.nrow(), P.nrow() );
    modelJacobian(P, J);
    matrix<flnum> Jt = J.t();
    matrix<flnum> JtJ = Jt*J;
    vector<flnum> B = Jt*E;
    compress(JtJ, B);

    flnum lambda = lambdaInit;
    for(iterations=0; iterations < maxIters && error > errorMax; iterations++) {
        matrix<flnum> H(JtJ);
        for(int i = 0; i < H.nrow(); i++)
            H(i,i) *= 1+lambda;
        vector<flnum> dP( P.nrow() );
        solveLU(H, B, dP);
        uncompress(dP);
        vector<flnum> tryP = P + dP;
        modelData(tryP, yModel);
        E = yData - yModel;
        flnum tryError = E.qnorm();
        if(ABS(tryError-error) <= relativeTol*error)
            break;
        if(tryError > error)
            lambda *= lambdaFact;
        else {
            lambda /= lambdaFact;
            error = tryError;
            P = tryP;
            modelJacobian(P, J);
            Jt = J.t();
            JtJ = Jt*J;
            B = Jt*E;
            compress(JtJ, B);
        }
    }
    return sqrt(error/yData.nrow());
}

} // namespace libNumerics
