/*
 * Copyright (c) 2011-2012, Peter Abeles. All Rights Reserved.
 *
 * This file is part of BoofCV (http://boofcv.org).
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

//package boofcv.examples;

import boofcv.abst.feature.disparity.StereoDisparity;
import boofcv.alg.distort.ImageDistort;
import boofcv.alg.geo.PerspectiveOps;
import boofcv.alg.geo.RectifyImageOps;
import boofcv.alg.geo.rectify.RectifyCalibrated;
import boofcv.core.image.ConvertBufferedImage;
import boofcv.factory.feature.disparity.DisparityAlgorithms;
import boofcv.factory.feature.disparity.FactoryStereoDisparity;
//import boofcv.gui.image.ShowImages;
import boofcv.gui.image.VisualizeImageData;
import boofcv.io.image.UtilImageIO;
import boofcv.misc.BoofMiscOps;
import boofcv.struct.calib.StereoParameters;
import boofcv.struct.image.ImageUInt8;
import boofcv.struct.image.ImageFloat32;
import georegression.struct.se.Se3_F64;
import org.ejml.data.DenseMatrix64F;



import boofcv.alg.filter.blur.GBlurImageOps;
import boofcv.abst.filter.blur.BlurFilter;
import boofcv.abst.filter.derivative.ImageHessianDirect;
import boofcv.factory.filter.derivative.FactoryDerivative;
import boofcv.alg.misc.GPixelMath;



import java.awt.image.BufferedImage;
import java.io.*;

/**
 * The disparity between two stereo images is used to estimate the range of objects inside
 * the camera's view.  Disparity is the difference in position between the viewed location
 * of a point in the left and right stereo images.  Because input images are rectified,
 * corresponding points can be found by only searching along image rows.
 *
 * Values in the disparity image specify how different the two images are.  A value of X indicates
 * that the corresponding point in the right image from the left is at "x' = x - X - minDisparity",
 * where x' and x are the locations in the right and left images respectively.  An invalid value
 * with no correspondence is set to a value more than (max - min) disparity.
 *
 * @author Peter Abeles
 */
public class pix {

	/**
	 * Computes the dense disparity between between two stereo images.  The input images
	 * must be rectified with lens distortion removed to work!  Floating point images
	 * are also supported.
	 *
	 * @param rectLeft Rectified left camera image
	 * @param rectRight Rectified right camera image
	 * @param minDisparity Minimum disparity that is considered
	 * @param maxDisparity Maximum disparity that is considered
	 * @return Disparity image
	 */
	public static ImageUInt8 denseDisparity( ImageUInt8 rectLeft , ImageUInt8 rectRight ,
											 int minDisparity , int maxDisparity ,
                                  int regionRadiusX, int regionRadiusY, double maxPerPixelError, int validateRtoL, double texture)
	{
		// A slower but more accuracy algorithm is selected
		// All of these parameters should be turned
      // regionWta(DisparityAlgorithms whichAlg, int minDisparity, int maxDisparity, 
      // int regionRadiusX, int regionRadiusY, double maxPerPixelError, 
      // int validateRtoL, double texture, java.lang.Class<T> imageType) 
		StereoDisparity<ImageUInt8,ImageUInt8> disparityAlg =
				FactoryStereoDisparity.regionWta(DisparityAlgorithms.RECT_FIVE,
						minDisparity, maxDisparity, regionRadiusX, regionRadiusY, maxPerPixelError, validateRtoL, texture, ImageUInt8.class);
                //minDisparity, maxDisparity, 3, 3, 20, 1, 0.2, ImageUInt8.class);

		// process and return the results
		disparityAlg.process(rectLeft,rectRight);


      // recover disparity
      ImageUInt8 disp = disparityAlg.getDisparity();

      // fix problems with the output disparities
      final byte[] data = disp.data;

      final int width  = disp.getWidth();
      final int height = disp.getHeight();
      final int stride = disp.stride;

      // HACK TO MANIPULATE UNSIGNED INTS
      for (int i = 0; i < height*width; i++) 
         if(data[i] < maxDisparity+1) 
            data[i] = (byte) ((int) data[i] + (int) minDisparity);

		return disp;
	}


	public static ImageFloat32 LoG( ImageFloat32 input )
	{
      ImageHessianDirect<ImageFloat32, ImageFloat32> hessian;

      ImageFloat32 ginput;
      ImageFloat32 derivXX;
      ImageFloat32 derivYY;
      ImageFloat32 derivXY;
      ImageFloat32 derivL;

      ginput  = new ImageFloat32(input.width,input.height);
      derivXX = new ImageFloat32(input.width,input.height);
      derivYY = new ImageFloat32(input.width,input.height);
      derivXY = new ImageFloat32(input.width,input.height);

      derivL  = new ImageFloat32(input.width,input.height);
      
      GBlurImageOps.gaussian(input, ginput, 1, -1, null);

      // 5-tap second order derivative correspond to the central derivatives of the image
      // http://boofcv.org/javadoc/boofcv/alg/filter/derivative/HessianThree.html
      // alternative
      // http://boofcv.org/javadoc/boofcv/alg/filter/derivative/HessianSobel.html
      hessian = FactoryDerivative.hessianDirectThree_F32();
      
      hessian.process(ginput, derivXX, derivYY,derivXY);

      GPixelMath.add(derivXX, derivYY, derivL);

		return derivL;
	}




	/**
	 * Rectified the input images using known calibration.
	 */
	public static void rectify( ImageUInt8 origLeft , ImageUInt8 origRight ,
								StereoParameters param ,
								ImageUInt8 rectLeft , ImageUInt8 rectRight )
	{
		// Compute rectification
		RectifyCalibrated rectifyAlg = RectifyImageOps.createCalibrated();
		Se3_F64 leftToRight = param.getRightToLeft().invert(null);

		// original camera calibration matrices
		DenseMatrix64F K1 = PerspectiveOps.calibrationMatrix(param.getLeft(), null);
		DenseMatrix64F K2 = PerspectiveOps.calibrationMatrix(param.getRight(), null);

		rectifyAlg.process(K1,new Se3_F64(),K2,leftToRight);

		// rectification matrix for each image
		DenseMatrix64F rect1 = rectifyAlg.getRect1();
		DenseMatrix64F rect2 = rectifyAlg.getRect2();
		// New calibration matrix,
		DenseMatrix64F rectK = rectifyAlg.getCalibrationMatrix();

		// Adjust the rectification to make the view area more useful
		RectifyImageOps.fullViewLeft(param.left, rect1, rect2, rectK);

		// undistorted and rectify images
		ImageDistort<ImageUInt8> imageDistortLeft =
				RectifyImageOps.rectifyImage(param.getLeft(), rect1, ImageUInt8.class);
		ImageDistort<ImageUInt8> imageDistortRight =
				RectifyImageOps.rectifyImage(param.getRight(), rect2, ImageUInt8.class);

		imageDistortLeft.apply(origLeft, rectLeft);
		imageDistortRight.apply(origRight, rectRight);
	}

	public static void main( String args[] ) {
//		String calibDir = "../data/applet/calibration/stereo/Bumblebee2_Chess/";
//		String imageDir = "../data/applet/stereo/";
//		String imageDir = "./";
//      System.out.print(args[0]);
//      System.out.print(args[1]);
//      System.out.print(args[2]);
      if (args.length < 3) {
         System.out.print("Usage:");
         System.out.print("   pix   im1.png im2.png out.pgm [mindisp(0) maxdisp(60) regionRadius(3) maxPerPixelError(20) validateRtoL(1) texture(0.2)]");
         return;
      }

//		StereoParameters param = BoofMiscOps.loadXML(calibDir + "stereo.xml");

		// load and convert images into a BoofCV format
//		BufferedImage origLeft = UtilImageIO.loadImage(imageDir + "a.png");
//		BufferedImage origRight = UtilImageIO.loadImage(imageDir + "b.png");
		BufferedImage origLeft = UtilImageIO.loadImage(args[0]);
		BufferedImage origRight = UtilImageIO.loadImage(args[1]);
      String outfile =  args[2];
      int mindisp=0; int maxdisp=60; int regionRadius=3; double maxPerPixelError=20; int validateRtoL=1; double texture=0.2;
      if (args.length > 3) mindisp           = Integer.parseInt(args[3]);
      if (args.length > 4) maxdisp           = Integer.parseInt(args[4]);
      if (args.length > 5) regionRadius      = Integer.parseInt(args[5]);
      if (args.length > 6) maxPerPixelError  = Double.parseDouble(args[6]);
      if (args.length > 7) validateRtoL      = Integer.parseInt(args[7]);
      if (args.length > 8) texture           = Double.parseDouble(args[8]);

		ImageUInt8 distLeft = ConvertBufferedImage.convertFrom(origLeft,(ImageUInt8)null);
		ImageUInt8 distRight = ConvertBufferedImage.convertFrom(origRight,(ImageUInt8)null);

		// rectify images
//		ImageUInt8 rectLeft = new ImageUInt8(distLeft.width,distLeft.height);
//		ImageUInt8 rectRight = new ImageUInt8(distRight.width,distRight.height);

//		rectify(distLeft,distRight,param,rectLeft,rectRight);

		// compute disparity
		ImageUInt8 disparity = denseDisparity(distLeft,distRight,mindisp,maxdisp,regionRadius,regionRadius,maxPerPixelError,validateRtoL,texture);

		// show results
		BufferedImage visualized = VisualizeImageData.disparity(disparity, null,mindisp,maxdisp,0);
		UtilImageIO.saveImage(visualized, outfile + ".png");

//		ShowImages.showWindow(rectLeft,"Rectified");
//		ShowImages.showWindow(visualized,"Disparity");


      try {
         BufferedWriter out = new BufferedWriter(new FileWriter(outfile));

         out.write("P2\n");
         out.write(String.format("%d %d\n", disparity.width, disparity.height));
         out.write("255\n");
         /* write data */
         int n=0;
         for(int y=0; y<disparity.height; y++) {
            for(int x=0; x<disparity.width; x++, n++) {
               out.write(String.format("%.8g ",(float)disparity.get(x,y)));
               if(n==16) { out.write("\n"); n = 0; }
            }
         }
         out.close();

      } catch (IOException e) {}


	}
}
