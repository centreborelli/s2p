// Copyright (C) 2019, Julien Michel (CNES) <julien.michel@cnes.fr>

#include <stdlib.h>

#include "Utilities/Parameters.h"
#include "LibImages/LibImages.h"
#include "LibSift/LibSift.h"

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
  /**
   * This function allows to free the float buffer from python.
  */
  void delete_buffer(float * buffer)
  {
    if (buffer != NULL)
      delete [] buffer;
  }

}
