// WARNING: 
// This file implements an algorithm possibly linked to the patent
//
// David Lowe  "Method and apparatus for identifying scale invariant 
// features in an image and use of same for locating an object in an 
// image",  U.S. Patent 6,711,293.
//
// This file is made available for the exclusive aim of serving as
// scientific tool to verify of the soundness and
// completeness of the algorithm description. Compilation,
// execution and redistribution of this file may violate exclusive
// patents rights in certain countries.
// The situation being different for every country and changing
// over time, it is your responsibility to determine which patent
// rights restrictions apply to you before you compile, use,
// modify, or redistribute this file. A patent lawyer is qualified
// to make this determination.
// If and only if they don't conflict with any patent terms, you
// can benefit from the following license terms attached to this
// file.
//
// This program is provided for scientific and educational only:
// you can use and/or modify it for these purposes, but you are
// not allowed to redistribute this work or derivative works in
// source or executable form. A license must be obtained from the
// patent right holders for any other use.

#include "demo_lib_sift.h"
#include "library.h"

#include "libIO/io_png.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

/// Load image from PNG file.
bool load(char* nameFile, float*& pix, int& w, int& h)
{
   size_t sw,sh;
   pix = read_png_f32_gray(nameFile, &sw,&sh);
   if(! pix) {
      std::cerr << "Unable to load image file " << nameFile << std::endl;
      return false;
   }

   w = static_cast<int>(sw);
   h = static_cast<int>(sh);
   return true;
}



// Writes the keypoints keys1 into filename in text format
// Length is the lenght of the descriptor (128 is the usual), 
// zoom is a multiplier for x,y,scale to be applied before writing
void write_keypoints( char* filename, keypointslist  keys1, int Length, float zoom1)
{
   int num_keys1 = keys1.size();
   // Write all the keypoints (row, col, scale, orientation, desciptor (128 integers)) to 
   // the file argv[3] (so that the users can match the keypoints with their own matching algorithm if they wish to)
   // keypoints in the 1st image
   FILE*file_key1= fopen(filename,"w");
   if ( file_key1  )
   {
      // Follow the same convention of David Lowe: 
      // the first line contains the number of keypoints and the length of the desciptors (128)
      fprintf(file_key1,"%d %d \n", num_keys1, Length);

      keypointslist::iterator ptr = keys1.begin();
      for(int i=0; i < (int) keys1.size(); i++, ptr++)	
      {
         // x y scale angle tilt longitude_angle_theta  <descriptor>
         fprintf(file_key1,"%.8g %.8g %.8g %.8g ", zoom1*ptr->x,zoom1*ptr->y,zoom1*ptr->scale,ptr->angle);

         for (int ii = 0; ii < (int) Length; ii++)
         {
            fprintf(file_key1,"%g  ", ptr->vec[ii]);
         }

         fprintf(file_key1,"\n");
      }
   }
   else 
   {
      std::cerr << "Unable to open the file " << filename; 
   }

   fclose(file_key1);
}


// Reads the keypoint file with format, one line
// > num_points descriptor_len
// followed by num_points lines
// > zoom*ptr->x zoom*ptr->y zoom*ptr->scale  ptr->angle  <descriptor>
// where  <descriptor> indicates the descriptor
int read_keypoints( char* filename,  keypointslist  &keys1, float zoom=1) {
   int num_keys1=0;
	std::ifstream file(filename);
	if (file.is_open())
   {
      char line[4096];
      int Length;

      // extra indiates that the keypoints contain the tilt and rotation of the asift keypoints
      file.getline(line,4096);
      sscanf(line ,"%d %d", &num_keys1, &Length );
       // printf("%d %d %d\n", num_keys1, Length, extra);

      // just check that we are talking about the same vectors
      if(Length != VecLength) printf("HORROR!!! the Lenght of the descriptor is wrong? \n");

      for(int t=0;t<num_keys1;t++) {
         keypoint newkeypoint;
         file.getline(line,4096);
         char* pline = line;
         int readchars;

         // read the beginning
         sscanf(pline,"%g %g %g %g%n", &newkeypoint.x,&newkeypoint.y,&newkeypoint.scale,&newkeypoint.angle,&readchars);
         pline+=readchars;

         // fix zoom and unused values
         newkeypoint.x     *=zoom;
         newkeypoint.y     *=zoom;
         newkeypoint.scale *=zoom;

         // read the descriptor
         for(int z=0;z<Length;z++) {
            sscanf(pline,"%g%n", &(newkeypoint.vec[z]), &readchars);
            pline+=readchars;
         }

         // store in the keys vector
         keys1.push_back(newkeypoint);
      }

   }
   else 
	{
		std::cerr << "Unable to open the file keypoints."; 
	}

	file.close();

   return num_keys1;
}

void setpixel(float* u, int w, int h, float x, float y, float val){
   int xx = (int)x;
   int yy = (int)y;
   if(xx < w && yy < h && xx>=0 && yy>=0)
      u[xx + yy*w] = val;
}



/// Usage: sift imgIn fileOut [imgOut]
/// Output in text file fileOut are the SIFT keypoints of imgIn 
/// the image imgIn (PNG format). If imgOut is in the argument list,
/// this is an image file where keypoints are shown (PNG format).
int main(int argc, char **argv)
{	
   if(argc != 3 && argc != 4) {
      std::cerr << "Usage: " << argv[0] << " imgIn keyFile [imgOut]"
         << std::endl;
      return 1;
   }

   //////////////////////////////////////////////// Input
   int w1, h1;
   float* ipixels1;
   if(! load(argv[1], ipixels1, w1, h1))
      return 1;

   ///////////////////////////////////////// Applying Sift
   siftPar siftparameters;
   default_sift_parameters(siftparameters);
   siftparameters.DoubleImSize=0;

   keypointslist keyp1;
   compute_sift_keypoints(ipixels1,keyp1,w1,h1,siftparameters);
   std::cout<< "sift:: 1st image: " << keyp1.size() << " keypoints"<<std::endl;

   //////////////////////////////////////////////////////////////// Save file with keypoints
	// save keypoints of the 1st image
   write_keypoints( argv[2], keyp1, 128, 1 );

   //////////////////////////////////////////////// Output image containing keypoints
   if(argc > 3) {
      int wo = w1;
      int ho = h1;

      float *opixels = new float[wo*ho];
      for(int j = 0; j < h1; j++)
         for(int i = 0; i < w1; i++)  opixels[j*wo+i] = ipixels1[j*w1+i];

      //////////////////////////////////////////////////////////////////// Draw matches
      keypointslist::iterator ptr = keyp1.begin();
      for(int i=0; i < (int) keyp1.size(); i++, ptr++) {
         setpixel(opixels, wo, ho,  ptr->x, ptr->y,  255.0);
         setpixel(opixels, wo, ho,  ptr->x+1, ptr->y, .0);
         setpixel(opixels, wo, ho,  ptr->x-1, ptr->y, .0);
         setpixel(opixels, wo, ho,  ptr->x, ptr->y+1, .0);
         setpixel(opixels, wo, ho,  ptr->x, ptr->y-1, .0);
      }


      ///////////////////////////////////////////////////////////////// Save imgOut	
      write_png_f32(argv[3], opixels, (size_t)wo, (size_t)ho, 1);
      delete[] opixels;
   }

   /////////////////////////////////////////////////////////////// Delete memory
   free(ipixels1);
   return 0;
}
