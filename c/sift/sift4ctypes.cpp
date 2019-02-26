#include <stdlib.h>

#include "Utilities/Parameters.h"
#include "LibImages/LibImages.h"
#include "LibSift/LibSift.h"

extern "C"{

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
    recordSize = descriptorSize + 2;
    
    // Allocate output buffer
    float * out = new float[recordSize*nbRecords];

    // Fill output buffer with keypoints
    std::list<KeyPoint*>::iterator key = sift.m_keyPoints->begin();
    size_t currentPoint  = 0;
    for(;key != sift.m_keyPoints->end();++key,++currentPoint)
    {
        size_t currentIndex = recordSize*currentPoint;
        out[currentIndex] = (*key)->getX();
        out[currentIndex+1] = (*key)->getY();

        for(unsigned int i = 0; i < descriptorSize;++i)
        {
            out[currentIndex+2+i] = (*key)->getPtrDescr()[i];
        } 
    }   
    return out;
  }

  void delete_buffer(float * buffer)
  {
    if (buffer != NULL)
      delete [] buffer;
  }

}
