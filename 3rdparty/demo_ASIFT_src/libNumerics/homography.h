#ifndef HOMOGRAPHY_H
#define HOMOGRAPHY_H

#include "matrix.h"

namespace libNumerics {

/// 2-D homography transform.
class Homography {
public:
    Homography();

    void setId();
    void setTrans(double dx, double dy);
    void setZoom(double zx, double zy);

    matrix<double>& mat() { return m_H; }
    const matrix<double>& mat() const { return m_H; }

    void operator()(double& x, double& y) const;
    Homography operator*(const Homography& rhs) const;
    Homography inverse() const;
private:
    matrix<double> m_H;
    void normalize();
};

/// Homography (and more restricted transforms) estimation.
class ComputeH {
public:
    enum Type { Translation, //                              (2 parameters)
                Rotation,    // Rotation/Translation         (3 parameters)
                Zoom,	     // Zoom/Translation             (3 parameters)
                GeneralZoom, // Non uniform zoom/Translation (4 parameters)
                Similarity,  // Zoom/Rotation/Translation    (4 parameters)
                Affine,      //                              (6 parameters)
                Projective   //                              (8 parameters)
    };
    static Type restrict(Type t); // Return less general motion
public:
    ComputeH(Type type);
    ~ComputeH();

    Type type() const { return _type; }
    void clear();

    /// Add corresponding points (x1,y1) and (x2,y2)
    void add(float x1, float y1, float x2, float y2, float w = 1.0f);
    /// Add corresponding lines of equation u x + v y + w = 0
    void add(float a1, float b1, float c1,
             float a2, float b2, float c2, float w = 1.0f);

    float weight() const; ///< Sum of weights (=#correspondences)
    float q_error(const Homography& map) const; ///< Quadratic error
    float compute(Homography& map) const; ///< LSE motion, return support weight
private:
    Type _type;
    int n; ///< Dimension of matrix = # unknown parameters
    double Ann[64], Bn[8], b; // Min (X 1) (A B) (X 1)^T is X^T = Ann^-1 Bn

    static int size(Type type);
    void add_4parameters(float x1, float y1, float x2, float y2, float w);
    void add_4parameters(float a1, float b1, float c1,
                         float a2, float b2, float c2, float w);
    void wrap(Homography& map, const vector<double>& v) const;
    void unwrap(const Homography& map, vector<double>& v) const;
    float q_error(const vector<double>& v) const; // Quadratic error

    bool compute_rotation(vector<double>& B) const;

    /// For Projective, data normalization is required
    class Normalization { public: double x, y, s; };
    bool normalize(Normalization& left,
                   matrix<double>& A, vector<double>& B,
                   Normalization& right) const;
    static bool de_normalize(const Normalization& left,
                             vector<double>& B,
                             const Normalization& right);
};

} // libNumerics

#endif
