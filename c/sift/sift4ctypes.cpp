// Copyright (C) 2019, Julien Michel (CNES) <julien.michel@cnes.fr>

#include <stdlib.h>
#include <math.h>

#include "Utilities/Parameters.h"
#include "LibImages/LibImages.h"
#include "LibSift/LibSift.h"
#include "linalg.c"


// Compute the SQUARE Euclidean distance.
float euclidean_distance_square(const float* x, const float* y, int length)
{
    float d = 0.0;
    for (int i = 0; i < length; i++) {
        float t = (x[i] - y[i]);
        d += t*t;
    }
    return d;
}

float distance_epipolar(const float* k1, const float* k2, int length, const float* s1, const float * s2,
                        const float epi_thresh, const unsigned int offset_desc)
{
    float a = s1[0];
    float b = s1[1];
    float c = s1[2];
    float d = s2[0];
    float e = s2[1];
    float f = s2[2];

    // Naming follows Ives' convention in which x is the row index
    float x1 = k1[1];
    float y1 = k1[0];
    float x2 = k2[1];
    float y2 = k2[0];

    // rectified x coordinates (in Ives' conventions x is the row index)
    float xx1 = b * y1 + a * x1 + c; // scalar product of (b, a, c) and (x1, y1, 1)
    float xx2 = e * y2 + d * x2 + f; // scalar product of (e, d, f) and (x2, y2, 1)

    // points satisfy the epipolar constraint when the rectified x are equal
    if (fabs(xx1 - xx2) < epi_thresh)
        return euclidean_distance_square(&k1[offset_desc], &k2[offset_desc], length);
    else
        return INFINITY;
}

float distance(const float* k1, const float* k2, int length, const float* s1, const float * s2,
               const float epi_thresh, const unsigned int offset_desc)
{
    // parameters used in distance_epipolar
    (void) s1;
    (void) s2;
    (void) epi_thresh;

    return euclidean_distance_square(&k1[offset_desc], &k2[offset_desc], length);

}

extern "C"{
  /**
   * This function is meant to be mapped to python using ctypes.
   *
   * It computes sifts points of input_buffer which is interpreted as a w x h image.
   * Keypoints are returned as a linear buffer of float of size recordSize * nbRecords.
   *
   * This buffer is the responsibiliy of the caller and should be freed by her.
   */
  float * sift(const float * input_buffer, const size_t w, const size_t h,
	       const float thresh_dog,
	       const unsigned int ss_noct,
	       const unsigned int ss_nspo,
	       unsigned int & recordSize,
	       unsigned int & nbRecords) {

    // Derive Image from buffer
    Image im(input_buffer, (const size_t) w, (const size_t) h, 1);

    // prepare params object
    Parameters params;
    params.setDefaultValues();
    params.set_thresh_dog(thresh_dog);
    params.set_noct(ss_noct);
    params.set_nspo(ss_nspo);

    // run sift
    Sift sift(params);
    sift.computeKeyPoints(im);

    // Compute the number of records
    nbRecords = sift.m_keyPoints->size();

    // Compute the record length
    size_t descriptorSize = 0;
    if(nbRecords > 0)
    {
        const KeyPoint * firstPoint = sift.m_keyPoints->front();
        descriptorSize = firstPoint->getNbOri() * firstPoint->getNbHist() * firstPoint->getNbHist();
    }
    recordSize = descriptorSize + 4;

    // Allocate output buffer
    float * out = new float[recordSize*nbRecords];

    // Fill output buffer with keypoints
    std::list<KeyPoint*>::iterator key = sift.m_keyPoints->begin();
    size_t currentPoint  = 0;
    for(;key != sift.m_keyPoints->end();++key,++currentPoint)
    {
        size_t currentIndex = recordSize*currentPoint;
        out[currentIndex] = (*key)->getY();
        out[currentIndex+1] = (*key)->getX();
	out[currentIndex+2] = (*key)->getSigma();
	out[currentIndex+3] = (*key)->getTheta();
        for(unsigned int i = 0; i < descriptorSize;++i)
        {
            out[currentIndex+4+i] = (*key)->getPtrDescr()[i];
        }
    }
    return out;
  }

  float * matching(const float * k1,
                const float * k2,
                const unsigned int length_desc,
                const unsigned int offset_desc,
                const unsigned int nb_sift_k1,
                const unsigned int nb_sift_k2,
                const float sift_thresh,
                const float epi_thresh,
                double * fund_mat,
                const bool use_fundamental_mat,
                const bool use_relative_method,
                unsigned int & nb_match){
    // Structure of k1 and k2 is supposed to be the following one :
    // 4 first numbers : pos_y pos_x scale orientation
    // length_desc floats representing the descriptors
    nb_match = 0;
    const float sift_thresh_square = sift_thresh*sift_thresh;
    float * matches = new float[nb_sift_k1 * 4];

    float (*keypoints_distance_overloaded)(const float *,
                                           const float *,
                                           int,
                                           const float*, const float*,
                                           const float,
                                           const unsigned int);

    // Use fundamental matrix if given
    float s1[3]; float s2[3];
    if (use_fundamental_mat){
        rectifying_similarities_from_affine_fundamental_matrix(s1, s2, fund_mat);
        keypoints_distance_overloaded = distance_epipolar;
    }
    else{
        keypoints_distance_overloaded = distance;
    }

    // Compute the Euclidian distance for every pairs of points
    for (int i = 0; i < nb_sift_k1; i++) {
        float distA = INFINITY;
        float distB = INFINITY;
        int indexA = -1;

        const float * curr_k1_desc = &k1[i*(length_desc+offset_desc)];
        for (int j = 0; j < nb_sift_k2; j++) {
            const float * curr_k2_desc = &k2[j*(length_desc+offset_desc)];
            float dist = keypoints_distance_overloaded(curr_k1_desc, curr_k2_desc, length_desc,
                         s1, s2, epi_thresh, offset_desc);
            // find_the_two_nearest_keys
            if (dist < distA) {
                distB = distA;
                distA = dist;
                indexA = j;
            } else if (dist < distB)
                distB = dist;
        }

        float val = distA;
        if (use_relative_method){
            val = distA / distB;
        }
        if (val < sift_thresh_square) {
            matches[nb_match*4] = k1[i*(length_desc+offset_desc)];
            matches[nb_match*4+1] = k1[i*(length_desc+offset_desc)+1];
            matches[nb_match*4+2] = k2[indexA*(length_desc+offset_desc)];
            matches[nb_match*4+3] = k2[indexA*(length_desc+offset_desc)+1];
            nb_match += 1;
        }
    }
    return matches;

  }

 /**
   * This function allows to free the float buffer from python.
  */
  void delete_buffer(float * buffer)
  {
    if (buffer != NULL)
      delete [] buffer;
  }

}
