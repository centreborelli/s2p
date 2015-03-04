// Authors: Unknown. Please, if you are the author of this file, or if you 
// know who are the authors of this file, let us know, so we can give the 
// adequate credits and/or get the adequate authorizations.

#ifdef RODRIGUES_H

namespace libNumerics {

template <class T>
matrix<T> skew(const vector<T>& v)
{
    assert(v.nrow() == 3);
    matrix<T> M(3,3);
    M(0,0) = M(1,1) = M(2,2) = 0;
    M(1,2) = -(M(2,1)=v(0));
    M(2,0) = -(M(0,2)=v(1));
    M(0,1) = -(M(1,0)=v(2));
    return M;
}

template <class T>
matrix<T> rotation(vector<T> w)
{
    assert(w.nrow() == 3);
    T n = sqrt(w.qnorm());
    T c = cos(n);
    matrix<T> R = c*matrix<T>::eye(3);
    if(n) {
        w /= n;
        R += skew(sin(n)*w);
        R += (1-c)*w*w.t();
    }   
    return R;
}

template <class T>
vector<T> rotationAxis(const matrix<T>& R)
{
    assert(R.nrow() == 3 && R.ncol() == 3);
    vector<T> w(3);
    T n = acos(0.5*(R.tr()-1));
    if(n == 0)
        w = 0;
    else {
        w(0) = R(2,1)-R(1,2);
        w(1) = R(0,2)-R(2,0);
        w(2) = R(1,0)-R(0,1);
        w *= n/(2*sin(n));
    }
    return w;
}

} // libNumerics

#endif // RODRIGUES_H
