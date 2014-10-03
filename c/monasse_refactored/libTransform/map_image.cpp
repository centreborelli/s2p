/* Apply homography to an image.
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

#include "map_image.h"
#include "spline.h"
#include "gauss_convol.h"
#include "TransformSize.h"
#include "libNumerics/numerics.h"
#include <cmath>

static std::pair<int,int> boundingBox(const libNumerics::Homography& map,
                                      int& w, int& h)
{
    TransformSize rect;
    rect.add(map, w, h);
    w = rect.w();
    h = rect.h();
    return std::make_pair(rect.x(), rect.y());
}

/// Min singular value of the Jacobian of \a H at point (\a x,\a y).
static double minSVJacob(const libNumerics::matrix<double>& H,
                         double x, double y)
{
    libNumerics::matrix<double> J(2,2);
    libNumerics::vector<double> X(3);
    X(0) = x; X(1) = y; X(2) = 1.0;
    X = H*X;
    double d = 1.0 / X(2);
    x = X(0)*d;
    y = X(1)*d;
    J(0,0) = d*(H(0,0)-H(2,0)*x);
    J(0,1) = d*(H(0,1)-H(2,1)*x);
    J(1,0) = d*(H(1,0)-H(2,0)*y);
    J(1,1) = d*(H(1,1)-H(2,1)*y);
    libNumerics::SVD svd(J);
    return svd.W()(1);
}

/// Compute min linear compression of \a H in rectangle [0,w]x[0,h].
static double minZoomOut(const libNumerics::matrix<double>& H, int w, int h)
{
    double z=1.0;
    z = std::min(z, minSVJacob(H,0.0,0.0));
    z = std::min(z, minSVJacob(H,w,0.0));
    z = std::min(z, minSVJacob(H,0.0,h));
    z = std::min(z, minSVJacob(H,w,h));
    return z;
}

/// Apply geometric transform to image.
///
/// The transformation \a map is applied to the image \a in and the result
/// stored in \a im. If \a adjustSize is \c true, \a im will be sized so that
/// it contains all the transformed rectangle, otherwise it stays at original
/// size.
///
/// The returned pair of integers is the offset of the returned image \a im
/// with respect to original image \a in. If \a adjustSize is \c false, this is
/// (0,0), otherwise the location of upper-left corner of \a im in pixel
/// coordinates of \a in.
///
/// Interpolation is done by spline. Anti-aliasing filter is optional.
///
/// \a vOut is the background value to put at pixels outside image.
std::pair<int,int> map_image(LWImage<float> in,
                             libNumerics::Homography map,
                             LWImage<float>& im,
                             int order, bool adjustSize,
                             bool antiAlias, float vOut)
{
    int w = in.w, h = in.h;
    float zoomOut = antiAlias?
        static_cast<float>( minZoomOut(map.mat(), w, h) ): 1.0f;
    const libNumerics::Homography oriMap(map);
    const int oriW=w, oriH=h;
    std::pair<int,int> offset(0,0);
    if(adjustSize) {
        offset = boundingBox(map, w, h);
        free(im.data);
        im = alloc_image<float>(w, h, in.comps);
    }
    if(zoomOut < 1.0f) {
        float zoomIn = 1.0f / zoomOut;
        // GF: added some extra space
        int wZoom=(int)std::ceil(w*zoomIn*1.5), hZoom=(int)std::ceil(h*zoomIn*1.5);
        LWImage<float> imZoom = alloc_image<float>(wZoom,hZoom,in.comps);
        libNumerics::matrix<double> mapZ(3,3);
        mapZ = 0.0;
        mapZ(0,0) = zoomIn;
        mapZ(1,1) = zoomIn;
        mapZ(2,2) = 1.0;
        map.mat() = mapZ*map.mat();
        map_image(in, map, imZoom, order, false, false, vOut);
        float sigma = 0.8f*sqrt(zoomIn*zoomIn-1.0f);
        gauss_convol(imZoom, sigma);
        map.mat() = 0.0;
        map.mat()(0,0)=zoomOut;
        map.mat()(1,1)=zoomOut;
        map.mat()(2,2)=1.0;
        in = imZoom;
    }
    LWImage<float> tmp = alloc_image(in);
    if( prepare_spline(tmp,order) ) {
        libNumerics::Homography inv = map.inverse();
        const int stepComp = im.stepComp();
        float* out = new float[im.comps];
        float* pixOut = im.data;
        for(int i = 0; i < im.h; i++)
            for(int j = 0; j < im.w; j++) {
                double x=j+offset.first, y=i+offset.second;
                inv(x,y);
                for(int k=0; k < im.comps; k++)
                    out[k] = vOut;
                interpolate_spline(tmp, order,
                                   static_cast<float>(x+.5),
								   static_cast<float>(y+.5), out);
                for(int k=0; k < im.comps; k++)
                    pixOut[k*stepComp] = out[k];
                pixOut += im.step();
            }
        delete [] out;
    }
    free(tmp.data);
    if(zoomOut < 1.0f) {
        free(in.data); // Was allocated above
        if(! is_number(vOut)) { // Put back mask
            libNumerics::Homography inv = oriMap.inverse();
            const int stepComp = im.stepComp();
            float* pixOut = im.data;
            for(int i = 0; i < im.h; i++)
                for(int j = 0; j < im.w; j++) {
                    double x=j+offset.first, y=i+offset.second;
                    inv(x,y);
                    if(x<0 || x>=oriW || y<0 || y>=oriH)
                        for(int k=0; k < im.comps; k++)
                            pixOut[k*stepComp] = NaN;
                    pixOut += im.step();
                }
        }
    }
    return offset;
}
