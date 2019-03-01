// Copyright (C) 2019, Julien Michel (CNES) <julien.michel@cnes.fr>

#include <stdlib.h>
#include <math.h>

#include "Utilities/Parameters.h"
#include "LibImages/LibImages.h"
#include "LibSift/LibSift.h"
#include "sift_anatomy_20141201/lib_util.h"

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
                const bool use_relative_method,
                unsigned int & nb_match){
    // Structure of k1 and k2 is supposed to be the following one :
    // 4 first numbers : pos_y pos_x scale orientation
    // length_desc floats representing the descriptors
    nb_match = 0;
    const float sift_thresh_square = sift_thresh*sift_thresh;
    float * matches = new float[nb_sift_k1 * 4];
    // TODO : add the fundamental matrix
    // Compute the Euclidian distance for every pairs of points
    for (int i = 0; i < nb_sift_k1; i++) {
        float distA = INFINITY;
        float distB = INFINITY;
        int indexA = -1;

        const float * curr_k1_desc = &k1[i*(length_desc+offset_desc) + offset_desc];
        for (int j = 0; j < nb_sift_k2; j++) {
            const float * curr_k2_desc = &k2[j*(length_desc+offset_desc) + offset_desc];
            //float dist = euclidean_distance_square(curr_k1_desc, curr_k2_desc, length_desc);
            float dist = 0.0;
            for (int ll = 0; ll < length_desc; ll++) {
                float t = (curr_k1_desc[ll] - curr_k2_desc[ll]);
                dist += t*t;
            }
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
