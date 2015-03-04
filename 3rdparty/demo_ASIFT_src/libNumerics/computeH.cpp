// Authors: Unknown. Please, if you are the author of this file, or if you 
// know who are the authors of this file, let us know, so we can give the 
// adequate credits and/or get the adequate authorizations.

#include "homography.h"
#include "numerics.h"

#include <algorithm>
#include <math.h> /* For sqrt */
#include <string.h>

static const float minEigenValue = 1e-3f; // For regular matrix

namespace libNumerics {

/// Constructor. Field `b' used only for error computation.
ComputeH::ComputeH(Type type)
: _type(type), n( size(type) ), b(0)
{
    clear();
}

// Destructor
ComputeH::~ComputeH()
{}

// Dimension of matrix w.r.t. type
int ComputeH::size(Type type)
{
    switch(type) {
    case Translation:
        return 2;
    case Zoom:
        return 3;
    case Rotation: // In fact 3, but nonlinear system
    case GeneralZoom:
    case Similarity:
        return 4;
    case Affine:
        return 6;
    case Projective:
        return 8;
    }
    return 8;
}

// Return less general motion
ComputeH::Type ComputeH::restrict(Type t)
{
    switch(t) {
    case Translation:
        return Translation; // Should return identity
    case Rotation:
    case Zoom:
        return Translation;
    case Similarity:
        return Zoom; // Rotation also correct. Arbitrary choice.
    case GeneralZoom:
        return Zoom;
    case Affine:
        return Similarity;
    case Projective:
        return Affine;
    }
    return Affine;
}

// Reinitialize
void ComputeH::clear()
{
    memset(Ann, 0, n*n*sizeof(double));
    memset(Bn,  0, n*sizeof(double));
    b = 0;
}

// Add two corresponding points
void ComputeH::add(float x, float y, float X, float Y, float w)
{
    if(_type <= Similarity) { // Separate for readability
        add_4parameters(x, y, X, Y, w);
        return;
    }
    double x2 = x*x, y2 = y*y, xy = x*y;
    double xX = x*X, yX = y*X, xY = x*Y, yY = y*Y;
    double *A = Ann, *B = Bn;

    *A++ += w* x2; // Equation 1
    *A++ += w* xy;
    A += 2;
    *A++ += w* x;
    A++;
    if(_type == Projective) {
        *A++ -= w* x*xX;
        *A++ -= w* x*yX;
    }
    *B++ += w* x*X;
	
    A++; // Equation 2
    *A++ += w* y2;
    A += 2;
    *A++ += w* y;
    A++;
    if(_type == Projective) {
        *A++ -= w* y*xX;
        *A++ -= w* y*yX;
    }
    *B++ += w* y*X;
	
    A +=2; // Equation 3
    *A++ += w* x2;
    *A++ += w* xy;
    A++;
    *A++ += w* x;
    if(_type == Projective) {
        *A++ -= w* x*xY;
        *A++ -= w* x*yY;
    }
    *B++ += w* x*Y;
	
    A +=3; // Equation 4
    *A++ += w* y2;
    A++;
    *A++ += w* y;
    if(_type == Projective) {
        *A++ -= w* y*xY;
        *A++ -= w* y*yY;
    }
    *B++ += w* y*Y;
	 
    A+= 4; // Equation 5
    *A++ += w;
    A++;
    if(_type == Projective) {
        *A++ -= w* xX;
        *A++ -= w* yX;
    }
    *B++ += w* X;
	
    A += 5; // Equation 6
    *A++ += w;
    *B++ += w* Y;
    if(_type == Projective) {
        *A++ -= w* xY;
        *A++ -= w* yY;
	
        A += 6; // Equation 7
        *A++ += w* (xX*xX + xY*xY);
        *A++ += w* (xX*yX + xY*yY);
        *B++ -= w* (xX*X  + xY*Y);
	
        A+= 7; // Equation 8
        *A++ += w* (yX*yX + yY*yY);
        *B++ -= w* (yX*X  + yY*Y);
    }
    b += w* (X*X + Y*Y);
}

// Add two corresponding points, type involving at most 4 parameters
void ComputeH::add_4parameters(float x, float y, float X, float Y, float w)
{
    double *A = Ann, *B = Bn;
	
    if(_type == Translation) {
        A[0] += w;
        A[3] += w;
        B[0] += w* (X - x);
        B[1] += w* (Y - y);
        b += w* ((X-x)*(X-x) + (Y-y)*(Y-y));
        return;
    }
    b += w* (X*X + Y*Y);
    if(_type == GeneralZoom) {
        A[0] += w* x*x;
        A[2] += w* x;
        B[0] += w* x*X;

        A[5] += w* y*y;
        A[7] += w* y;
        B[1] += w* y*Y;

        A[10]+= w;
        B[2] += w* X;

        A[15]+= w;
        B[3] += w* Y;
        return;
    }

    *A++ += w* (x*x + y*y); // Equation 1
    if(_type != Zoom) // Similarity or Rotation
        A++;
    *A++ += w* x;
    *A++ += w* y;
    *B++ += w* (x*X + y*Y);

    if(_type != Zoom) { // Similarity or Rotation
        A++; // Equation 2
        *A++ += w* (x*x + y*y);
        *A++ += w* y;
        *A++ -= w* x;
        *B++ += w* (y*X - x*Y);
        A++; // Prepare for next line
    }

    A++; // Equation 3
    *A++ += w;
    A++;
    *B++ += w* X;

    A += n-1; // Equation 4
    *A++ += w;
    *B++ += w* Y;
}

// Add corresponding lines of equation ux + by + x = 0
void ComputeH::add(float x, float y, float z, float X, float Y, float Z,
                    float w)
{
    float s = 1.0f / (float)sqrt(x*x + y*y);
    x *= s;
    y *= s;
    z *= s;
    s = 1.0f / (float)sqrt(X*X + Y*Y);
    X *= s;
    Y *= s;
    Z *= s;
    if(_type <= Similarity) { // Separate for readability
        add_4parameters(x, y, z, X, Y, Z, w);
        return;
    }

    double x2 = x*x, y2 = y*y, z2 = z*z, xy = x*y, xz = x*z, yz = y*z;
    double X2 = X*X, Y2 = Y*Y, Z2 = Z*Z, XY = X*Y, XZ = X*Z, YZ = Y*Z;
    double *A = Ann, *B = Bn;

    *A++ += w* (y2+z2) * X2; // Equation 1
    *A++ -= w* xy * X2;
    *A++ += w* (y2+z2) * XY;
    *A++ -= w* xy * XY;
    *A++ -= w* xz * X2;
    *A++ -= w* xz * XY;
    if(_type == Projective) {
        *A++ += w* (y2+z2) * XZ;
        *A++ -= w* xy * XZ;
    }
    *B++ += w* xz * XZ;

    A++; // Equation 2
    *A++ += w* (x2+z2) * X2;
    *A++ -= w* xy * XY;
    *A++ += w* (x2+z2) * XY;
    *A++ -= w* yz * X2;
    *A++ -= w* yz * XY;
    if(_type == Projective) {
        *A++ -= w* xy * XZ;
        *A++ += w* (x2+z2) * XZ;
    }
    *B++ -= w* yz * XZ;
	
    A += 2; // Equation 3
    *A++ += w* (y2+z2) * Y2;
    *A++ -= w* xy * Y2;
    *A++ -= w* xz * XY;
    *A++ -= w* xz * Y2;
    if(_type == Projective) {
        *A++ += w* (y2+z2) * YZ;
        *A++ -= w* xy * YZ;
    }
    *B++ += w* xz * YZ;

    A += 3; // Equation 4
    *A++ += w* (x2+z2) * Y2;
    *A++ -= w* yz * XY;
    *A++ -= w* yz * Y2;
    if(_type == Projective) {
        *A++ -= w* xy * YZ;
        *A++ += w* (x2+z2) * YZ;
    }
    *B++ += w* yz * YZ;

    A += 4; // Equation 5
    *A++ += w* X2; // *(x2+y2=1)
    *A++ += w* XY; // *(x2+y2=1)
    if(_type == Projective) {
        *A++ -= w* xz * XZ;
        *A++ -= w* yz * XZ;
    }
    *B++ -= w* XZ; // *(x2+y2=1)

    A += 5; // Equation 6
    *A++ += w* Y2; // *(x2+y2=1)
    *B++ -= w* YZ; // *(x2+y2=1)
    if(_type == Projective) {
        *A++ -= w* xz * YZ;
        *A++ -= w* yz * YZ;

        A += 6; // Equation 7
        *A++ += w* (y2+z2) * Y2;
        *A++ -= w* xy * Z2;
        *B++ += w* xz * Z2;

        A += 7; // Equation 8
        *A++ += w* (x2+z2) * Z2;
        *B++ += w* yz * Z2;
    }
    b += w* Z2; // *(x2+y2=1) 
}

// Add two corresponding lines, type involving at most 4 parameters
void ComputeH::add_4parameters(float x, float y, float z,
                                float X, float Y, float Z, float w)
{
    double x2 = x*x, y2 = y*y, z2 = z*z, xy = x*y, xz = x*z, yz = y*z;
    double X2 = X*X, Y2 = Y*Y, Z2 = Z*Z, XY = X*Y, XZ = X*Z, YZ = Y*Z;
    double *A = Ann, *B = Bn;
    if(_type == Translation) {
        *A++ += w* X2; // *(x2+y2=1)
        *A++ += w* XY; // *(x2+y2=1)
        *B++ += w* (yz*XY + xz*X2 - XZ/* *(x2+y2=1) */);

        A++;
        *A++ += w* Y2; // *(x2+y2=1)
        *B++ += w* (xz*XY + yz*Y2 - YZ/* *(x2+y2=1) */);

        b += w* (z2 + Z2 + y2*X2 + x2*Z2 - 2*(xz*XZ + yz*YZ + xy*XZ));
        return;
    }
    b += w* Z2; // *(x2+y2=1)
    if(_type == GeneralZoom) {
        *A++ += w* (y2+z2) * X2;
        *A++ -= w* xy * XY;
        *A++ -= w* xz * X2;
        *A++ -= w* xz * XY;
        *B++ += w* xz * XZ;

        A++;
        *A++ += w* (x2+z2) * Y2;
        *A++ -= w* yz * XY;
        *A++ -= w* yz * Y2;
        *B++ += w* yz * YZ;

        A += 2;
        *A++ += w* X2; // *(x2+y2=1)
        *A++ += w* XY; // *(x2+y2=1)
        *B++ -= w* XZ; // *(x2+y2=1)

        A += 3;
        *A++ += w* Y2; // *(x2+y2=1)
        *B++ -= w* YZ; // *(x2+y2=1)
        return;
    }
	
    if(_type == Zoom) {
        *A++ += w* (z2/* *(X2+Y2=1)*/ + y2*X2 + x2*Y2 - 2*xy*XY);
        *A++ -= w* (yz*XY + xz*X2);
        *A++ -= w* (yz*Y2 + xz*XY);
        *B++ += w* (yz*YZ + xz*X2);
    } else { // Similarity or Rotation
        *A++ += w* (1 /* =x2+y2*/+ 2*(z2 - xy)) * X2;
        *A++ += w* (x2 - y2) * XY;
        *A++ -= w* (xz + yz) * X2;
        *A++ -= w* (xz + yz) * XY;
        *B++ += w* (xz + yz) * XZ;

        A++;
        *A++ += w* (1 /* =x2+y2*/+ 2*(z2 + xy)) * Y2;
        *A++ += w* (xz - yz) * XY;
        *A++ += w* (xz - yz) * Y2;
        *B++ += w* (yz - xz) * YZ;
        A++; // Prepare for next line
    }

    A++;
    *A++ += w* X2; // *(x2+y2=1)
    *A++ += w* XY; // *(x2+y2=1)
    *B++ -= w* XZ; // *(x2+y2=1)

    A += n-1;
    *A++ += w* Y2; // *(x2+y2=1)
    *B++ -= w* YZ; // *(x2+y2=1)
}	

// Wrap vector of unknowns `v' into structure `map'
void ComputeH::wrap(Homography& h, const vector<double>& v) const
{
    int i = 0;
    h.mat()(0,0) = (_type==Translation)? 1.0f: v(i++);
    h.mat()(0,1) = (_type==Translation || _type==Zoom || _type==GeneralZoom) ?
        0: v(i++);
    if(n >= 6) {
        h.mat()(1,0) = v(i++);
        h.mat()(1,1) = v(i++);
    } else {
        h.mat()(1,0) = -h.mat()(0,1);
        h.mat()(1,1) = (_type==GeneralZoom)? v(i++): h.mat()(0,0);
    }
    h.mat()(0,2) = v(i++);
    h.mat()(1,2) = v(i++);
    if(_type == Projective) {
        h.mat()(2,0) = v(i++);
        h.mat()(2,1) = v(i++);
    } else
        h.mat()(2,0) = h.mat()(2,1) = 0;
    h.mat()(2,2) = 1.0;
}

/// Unwrap parameters in \a h into vector of unknowns \a v.
void ComputeH::unwrap(const Homography& h, vector<double>& v) const
{
    int i = 0;
    if(_type != Translation) {
        v(i++) = h.mat()(0,0);
        if(_type != Zoom) {
            if(_type != GeneralZoom) {
                v(i++) = h.mat()(0,1); // Rotation or Similarity or...
                if(n >= 6) // Affine or Projective
                    v(i++) = h.mat()(1,0);
            }
            if(_type==GeneralZoom || _type==Affine || _type==Projective)
                v(i++) = h.mat()(1,1);
        }
    }
    v(i++) = h.mat()(0,2);
    v(i++) = h.mat()(1,2);
    if(_type == Projective) {
        v(i++) = h.mat()(2,0);
        v(i++) = h.mat()(2,1);
    }
}

// Sum of weights (=#correspondences)
float ComputeH::weight() const
{
    // Diagonal coefficient affecting translation
    int i = (_type == Projective) ? 6 : n;
    return static_cast<float>(Ann[(i-1)*(n+1)]); // Element (i-1,i-1)
}


// Return quadratic error when mapping with `motion'
float ComputeH::q_error(const Homography& map) const
{
    vector<double> v(n);
    unwrap(map, v);
    return q_error(v);
}

// Idem, with arguments in a vector
float ComputeH::q_error(const vector<double>& v) const
{
    double e = b;
    // Diagonal terms
    const double* A = Ann + n*n-1;
    for(int i = n-1; i >= 0; i--, A -= n+1)
        e += *A * v(i) * v(i);
    // Cross terms
    A = Ann + (n-1)*n; // Last row
    for(int i = n-1; i >= 0; i--, A -= n) {
        double vi = v(i);
        e -= 2.0 * Bn[i] * vi;
        for(int j = n-1; j > i; j--)
            e += 2.0 * A[j] * vi * v(j);
    }
    return static_cast<float>(e);
}

// LSE for rotation: solve linear system under quadratic constraint
bool ComputeH::compute_rotation(vector<double>& B) const
{
    if(Ann[15] <= 0) // No point added or absurd value
        return false;
    B(0) = Ann[15] * Bn[0] - Ann[2] * Bn[2] - Ann[3] * Bn[3];
    B(1) = Ann[15] * Bn[1] - Ann[3] * Bn[2] + Ann[2] * Bn[3];
    double root = sqrt(B(0)*B(0) + B(1)*B(1));
    if(root < minEigenValue)
        return false;
    // Test first solution
    double lambda1 = (Ann[2]*Ann[2] + Ann[3]*Ann[3] + root) / Ann[15];
    B(0) /= root;
    B(1) /= root;
    B(2) = (-Ann[2]*Bn[0] - Ann[3]*Bn[1] + lambda1 * Bn[2]) / root;
    B(3) = (-Ann[3]*Bn[0] + Ann[2]*Bn[1] + lambda1 * Bn[3]) / root;
    float v1 = q_error(B);
    // Test second solution
    vector<double> C(4);
    double lambda2 = (Ann[2]*Ann[2] + Ann[3]*Ann[3] - root) / Ann[15];
    C(0) = -B(0);
    C(1) = -B(1);
    C(2) = -(-Ann[2]*Bn[0] - Ann[3]*Bn[1] + lambda2 * Bn[2]) / root;
    C(3) = -(-Ann[3]*Bn[0] + Ann[2]*Bn[1] + lambda2 * Bn[3]) / root;
    if(v1 > q_error(C)) // Keep second solution
        B = C;
    return true;
}

// Return LSE motion and the sum of weights
float ComputeH::compute(Homography& map) const
{
    vector<double> B(n);
    B.read(Bn);

    if(_type == Rotation) {
        if(! compute_rotation(B))
            return 0;
    } else {
        matrix<double> A(n,n);
        A.read(Ann);
        Normalization left, right;
        if(_type == Projective && !normalize(left, A, B, right))
            return 0;
        A.symUpper();

        vector<double> oldB(B);
        if(! solveLU(A, B))
            return 0;

        if(_type == Projective && ! de_normalize(left, B, right))
            return 0;
    }

    wrap(map, B);
    return weight();
}

// Normalize independently original and final points so that the new
// origin is their centroid and their mean square distance (to it) is 2
bool ComputeH::normalize(Normalization& left,
                         matrix<double>& A, vector<double>& B,
                         Normalization& right) const
{
    double w = A(5,5); // Total weight
    if(w < minEigenValue)
        return false;
    double invW = 1.0 / w;

    // Find normalizations (zoom-translation)
    right.s = (A(0,0) + A(1,1)) - (A(0,4)*A(0,4) + A(1,4)*A(1,4))*invW;
    if(right.s < minEigenValue)
        return false;
    right.s = sqrt(2.0*w / right.s);
    right.x = - invW * right.s * A(0,4);
    right.y = - invW * right.s * A(1,4);

    left.s = b - (B(4)*B(4) + B(5)*B(5))*invW;
    if(left.s < minEigenValue)
        return false;
    left.s = sqrt(2.0*w / left.s);
    left.x = - invW * left.s * B(4);
    left.y = - invW * left.s * B(5);
    double norm = left.x*left.x + left.y*left.y;

    double s2 = right.s*right.s, sS = right.s*left.s, S2 = left.s*left.s;

    // Normalization of vector B
    double b0 = B(0), b1 = B(1), b2 = B(2), b3 = B(3);
    B(0) = sS * B(0) - w*right.x*left.x;
    B(1) = sS * B(1) - w*right.y*left.x;
    B(2) = sS * B(2) - w*right.x*left.y;
    B(3) = sS * B(3) - w*right.y*left.y;
    B(4) = B(5) = 0;
    B(6) = sS*(left.s*B(6) - 2*(left.x*b0 + left.y*b2)) +
        w*right.x*(norm - 2.0);
    B(7) = sS*(left.s*B(7) - 2*(left.x*b1 + left.y*b3)) +
        w*right.y*(norm - 2.0);

    // Normalization of matrix A
    double a0 = A(0,0), a1 = A(0,1), a6 = A(0,6), a7 = A(0,7), a9 = A(1,1);
    double a15 = A(1,7), a22 = A(2,6), a23 = A(2,7), a31 = A(3,7);

    A(0,0) = s2 * A(0,0) - w*right.x*right.x;
    A(0,1) = s2 * A(0,1) - w*right.x*right.y;
    A(0,4) = 0;
    A(0,6) = right.s*(sS*A(0,6) - right.s*left.x*a0 - left.s*right.x*b0) +
        w*right.x*left.x*right.x - right.x * B(0);
    A(0,7) = right.s*(sS*A(0,7) - right.s*left.x*a1 - left.s*right.x*b1) +
        w*right.x*left.x*right.y - right.y * B(0);

    A(1,1) = s2 * A(1,1) - w*right.y*right.y;
    A(1,4) = 0;
    A(1,6) = A(0,7);
    A(1,7) = right.s*(sS*A(1,7) - right.s*left.x*a9 - left.s*right.y*b1) +
        w*right.y*left.x*right.y - right.y * B(1);

    A(2,2) = A(0,0);
    A(2,3) = A(0,1);
    A(2,5) = 0;
    A(2,6) = right.s*(sS*A(2,6) - right.s*left.y*a0 - left.s*right.x*b2) +
        w*right.x*left.y*right.x - right.x * B(2);
    A(2,7) = right.s*(sS*A(2,7) - right.s*left.y*a1 - left.s*right.x*b3) +
        w*right.x*left.y*right.y - right.y * B(2);

    A(3,3) = A(1,1);
    A(3,5) = 0;
    A(3,6) = A(2,7);
    A(3,7) = right.s*(sS*A(3,7) - right.s*left.y*a9 - left.s*right.y*b3) +
        w*right.y*left.y*right.y - right.y * B(3);

    A(4,6) = -B(0);
    A(4,7) = -B(1);

    A(5,6) = -B(2);
    A(5,7) = -B(3);

    A(6,6) = s2*(S2*A(6,6) - 2*left.s*(left.x*a6+left.y*a22) + a0*norm) -
        2*right.x*(B(6) + w*right.x);
    A(6,7) = s2*(S2*A(6,7) - 2*left.s*(left.x*a7+left.y*a23) + a1*norm) -
        right.x*B(7) - right.y*B(6) - 2*w*right.x*right.y;

    A(7,7) = s2*(S2*A(7,7) - 2*left.s*(left.x*a15+left.y*a31) + a9*norm) -
        2*right.y*(B(7) + w*right.y);
    return true;
}

// `l' (left) and 'r' (right) representing zoom-translation normalizations,
// and `B' the parameters of a projective motion,
// compute l^-1 B r
bool ComputeH::de_normalize(const Normalization& l,
                            vector<double>& B,
                            const Normalization& r)
{
    // B := B r
    B(4) += r.x * B(0) + r.y * B(1); // Line 1
    B(0) *= r.s;
    B(1) *= r.s;

    B(5) += r.x * B(2) + r.y * B(3); // Line 2
    B(2) *= r.s;
    B(3) *= r.s;

    double f = r.x * B(6) + r.y * B(7) + 1.0; // Line 3
    if(-minEigenValue < f && f < minEigenValue)
        return false; // Origin of right normalization on line at infinity
    B(6) *= r.s;
    B(7) *= r.s;

	// B := l^-1 B
    double s = 1.0 / (l.s * f);
    B(0) = (B(0) - l.x*B(6)) * s; // Line 1
    B(1) = (B(1) - l.x*B(7)) * s;
    B(4) = (B(4) - l.x* f  ) * s;

    B(2) = (B(2) - l.y*B(6)) * s; // Line 2
    B(3) = (B(3) - l.y*B(7)) * s;
    B(5) = (B(5) - l.y* f  ) * s;

    B(6) /= f; // Line 3
    B(7) /= f;
    return true;
}

} // libNumerics
