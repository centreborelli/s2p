// Authors: Unknown. Please, if you are the author of this file, or if you 
// know who are the authors of this file, let us know, so we can give the 
// adequate credits and/or get the adequate authorizations.

#ifndef NUMERICS_H
#define NUMERICS_H

#include "matrix.h"
#include <vector>

namespace libNumerics {
    class NumericsException {};
    class SvdConvergenceError : public NumericsException {};
    typedef double flnum;

    /// Solve system AX = B.
    bool solveLU(const matrix<flnum>& A, const vector<flnum>& B,
                 vector<flnum>& X);
    bool solveLU(matrix<flnum> A, vector<flnum>& B);

    /// Singular Value Decomposition
    class SVD {
    public:
        SVD(const matrix<flnum>& A);
        matrix<flnum>& U() { return m_U; }
        vector<flnum>& W() { return m_W; }
        matrix<flnum>& V() { return m_V; }
        matrix<flnum> compose() const;

    private:
        matrix<flnum> m_U, m_V;
        vector<flnum> m_W;
        static flnum withSignOf(flnum a, flnum b);
        static flnum hypot(flnum a, flnum b);
        static void rotate(flnum& a, flnum& b, flnum c, flnum s);
        void compute();
        void sort();
    };

    /// Levenberg-Marquardt minimization.
    class MinLM {
        static const flnum DEFAULT_RELATIVE_TOL;
        static const flnum DEFAULT_LAMBDA_INIT;
        static const flnum DEFAULT_LAMBDA_FACT;
        static const flnum EPSILON_KERNEL;
    public:
        MinLM();
        flnum minimize(vector<flnum>& P, const vector<flnum>& ydata,
                       flnum targetRMSE=0.1, int maxIters=300);
        virtual void modelData(const vector<flnum>& P,
                               vector<flnum>& ymodel) const = 0;
        virtual void modelJacobian(const vector<flnum>& P,
                                   matrix<flnum>& J) const = 0;
        int iterations;
        flnum relativeTol;
        flnum lambdaInit;
        flnum lambdaFact;
    private:
        std::vector<int> m_nullCols;
        void compress(matrix<flnum>& JtJ, vector<flnum>& B);
        void uncompress(vector<flnum>& B);
    };

} // namespace libNumerics

#endif
