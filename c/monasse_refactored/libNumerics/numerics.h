/* Linear system solving, SVD, Levenberg-Marquardt minimization.
    Copyright (C) 2010 Pascal Monasse <monasse@imagine.enpc.fr>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

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
        virtual ~MinLM() {}
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
