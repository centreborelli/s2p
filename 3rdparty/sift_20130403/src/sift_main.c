/*
IPOL SIFT
Copyright (C) 2013, Ives Rey Otero, CMLA ENS Cachan
<ives.rey-otero@cmla.ens-cachan.fr>

Version 20130318 (March 18, 2013)

== Patent Warning and Licence =================================================


This program is free software: you can use, modify and/or
redistribute it under the terms of the simplified BSD
License. You should have received a copy of this license along
this program. If not, see
<http://www.opensource.org/licenses/bsd-license.html>.



The SIFT method is patented 

    [3] "Method and apparatus for identifying scale invariant features
      in an image."
        David G. Lowe
        Patent number: 6711293
        Filing date: Mar 6, 2000
        Issue date: Mar 23, 2004
        Application number: 09/519,89
  
 These source codes are made available for the exclusive aim of serving as
 scientific tool to verify the soundness and completeness of the algorithm
 description. Compilation, execution and redistribution of this file may
 violate patents rights in certain countries. The situation being different
 for every country and changing over time, it is your responsibility to
 determine which patent rights restrictions apply to you before you compile,
 use, modify, or redistribute this file. A patent lawyer is qualified to make
 this determination. If and only if they don't conflict with any patent terms,
 you can benefit from the following license terms attached to this file.




*/
/**
 * @file sift_main.c
 * @brief The SIFT algorithmic chain
 *
 *
 * @author Ives Rey Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */


#include "sift_main.h"


/**  @brief Apply the SIFT transform to an image
 * 
 * \input    : image, of width X height pixels.
 * \output   : keys, list of extracted SIFT keypoints.
 * 
 * \param label  : name used to label the produced files.
 * 
 * \param flag   : verbosity flag
 *               if flag >= 0  : outputs one list of described keypoint
 *               if flag >= 1  : outputs lists of keypoints at each step
 *               if flag >= 2  : outputs lists and scale-space
 * 
 * \param  n_oct,
 *         n_spo,
 *         sigma_min,
 *         delta_min,
 *         sigma_in,
 *         C_DoG,
 *         C_edge,
 *         n_bins,
 *         lambda_ori,
 *         t,
 *         n_hist,
 *         n_ori,         
 *         lambda_descr          : SIFT method parameters
 */
void sift_transform(double* image, int width, int height,
                    struct kypt_lst* keysOUT, char* label,
                    int flag,
                    int n_oct,
                    int n_spo,
                    double sigma_min,
                    double delta_min,
                    double sigma_in,
                    double C_DoG,
                    double C_edge,
                    int n_bins,
                    double lambda_ori,
                    double t,
                    double n_hist,
                    double n_ori,
                    double lambda_descr){
        
    /** number octaves, limited by the maximum number of subsamplings */
    int noct = MIN(n_oct, (int)(log(MIN(width,height)/delta_min)/log(2)));
//     double sigma_min2 = sigma_min/delta_min; // 1.6 
    
    /** MEMORY ALLOCATION **/
    /** scale-space structure */
    //struct scsp* scalespace    = malloc_scsp_lowe_scalespace_bis(noct,n_spo,width,height,delta_min,sigma_min2);     //TODO sigma_min2
    struct scsp* scalespace    = malloc_scsp_lowe_scalespace_bis(noct,n_spo,width,height,delta_min,sigma_min);
    struct scsp* dog           = malloc_scsp_dog_from_scalespace(scalespace);
    struct scsp* dx_scalespace = malloc_scsp_from_model(scalespace);
    struct scsp* dy_scalespace = malloc_scsp_from_model(scalespace);
    /** list-of-keypoints */
    struct kypt_lst* keys                   = malloc_kypt_lst();   /* 3D (discrete) extrema   */
    struct kypt_lst* keysAcceptNoise        = malloc_kypt_lst();   /* passing the threshold on DoG  */
    struct kypt_lst* keysRejectNoise        = malloc_kypt_lst();   /* discarded by threshold on DoG */
    struct kypt_lst* keysAcceptNoiseSOFT    = malloc_kypt_lst();   /* passing the threshold on DoG  */
    struct kypt_lst* keysRejectNoiseSOFT    = malloc_kypt_lst();   /* discarded by threshold on DoG */
    struct kypt_lst* keysInterpol           = malloc_kypt_lst();   /* interpolated 3D extrema (continuous) */
    struct kypt_lst* keysRejectInterp       = malloc_kypt_lst();   /* discarded by interpolation */
    struct kypt_lst* keysAcceptEdge         = malloc_kypt_lst();   /* passing OnEdge filter */
    struct kypt_lst* keysRejectEdge         = malloc_kypt_lst();   /* discarded (on edge) */
    struct kypt_lst* cp_keys_multiori       = malloc_kypt_lst();   /* subset of keypoints with multiple orientation */
      
    
    /** KEYPOINT DETECTION ***************************************************/
    build_lowe_scalespace_ter(scalespace, image, width, height, sigma_in);  /* Builds Lowe's scale-space */    
    
    
    build_dog(scalespace,dog);
    find_candidate_keypoints(dog,keys,n_ori,n_hist,n_bins);
 
    discard_detections_duetonoise(keys, keysAcceptNoiseSOFT, keysRejectNoiseSOFT, 0.8*C_DoG);

    interpolate_position_of_candidate_keypoints(dog,keysAcceptNoiseSOFT, keysInterpol, keysRejectInterp);

    discard_detections_duetonoise(keysInterpol, keysAcceptNoise, keysRejectNoise, C_DoG); 

    compute_HarrisStephen_edge_response(dog,keysAcceptNoise); 
    discard_detections_onedge(keysAcceptNoise, keysAcceptEdge, keysRejectEdge, (C_edge+1)*(C_edge+1)/C_edge); 

    /** KEYPOINT DESCRIPTION *************************************************/
        
    build_scalespace_gradient(scalespace,dx_scalespace,dy_scalespace); /* Pre-computes gradient scale-space */
        
    attribute_orientation_to_keypoints(dx_scalespace,dy_scalespace,keysAcceptEdge,keysOUT,cp_keys_multiori,n_bins,lambda_ori,t); 
        
    attribute_featurevector_to_keypoints(dx_scalespace,dy_scalespace,keysOUT,n_hist,n_ori,lambda_descr);    
    
    /** OUTPUT ***************************************************************/
    char nom[256];
    if(flag==0){
        print_kypt_lst(keysOUT,n_hist*n_hist*n_ori);
    }
    if(flag>0){
        print_kypt_lst_extra(keysOUT,n_hist*n_hist*n_ori,n_bins);
        sprintf(nom,"extra_NES_%s.txt",label);              save_kypt_lst_detection(keys                ,nom);
        sprintf(nom,"extra_DoGSoftThresh_%s.txt",label);    save_kypt_lst_detection(keysAcceptNoiseSOFT ,nom);
        sprintf(nom,"extra_ExtrInterp_%s.txt",label);       save_kypt_lst_detection(keysInterpol        ,nom);
        sprintf(nom,"extra_ExtrInterpREJ_%s.txt",label);    save_kypt_lst_detection(keysRejectInterp    ,nom);
        sprintf(nom,"extra_DoGThresh_%s.txt",label);        save_kypt_lst_detection(keysAcceptNoise     ,nom);     
        sprintf(nom,"extra_OnEdgeResp_%s.txt",label);       save_kypt_lst_detection(keysAcceptEdge      ,nom);
        sprintf(nom,"extra_OnEdgeRespREJ_%s.txt",label);    save_kypt_lst_detection(keysRejectEdge      ,nom);
        sprintf(nom,"extra_OriAssignedMULT_%s.txt",label);  save_kypt_lst_detection(cp_keys_multiori    ,nom);
    }
    if(flag>1){
        sprintf(nom,"scalespace_%s",label);                 print_scsp_gray_nearestneighor(scalespace,nom);
        sprintf(nom,"DoG_%s",label);                        print_scsp_rgb(dog,nom);
    }
        
    /** FREEING MEMORY *******************************************************/
    /** list of keys*/
    free_kypt_lst(keys);
    free_kypt_lst(keysAcceptNoise);
    free_kypt_lst(keysRejectNoise);
    free_kypt_lst(keysAcceptNoiseSOFT);
    free_kypt_lst(keysRejectNoiseSOFT);   
    free_kypt_lst(keysInterpol);
    free_kypt_lst(keysRejectInterp);
    free_kypt_lst(keysAcceptEdge);
    free_kypt_lst(keysRejectEdge);
    free_kypt_lst(cp_keys_multiori);
    /** scalespace structures*/
    free_scsp(scalespace);
    free_scsp(dx_scalespace);
    free_scsp(dy_scalespace);
    free_scsp(dog);
    
}



/** @brief Compute the Gaussian scalespace for SIFT
 * 
 * 
 *  in @param input image  of im_width X im_height
 *     @param sigma0 : assumed level of blur in the image
 *  
 *  
 *  Other construction parameters are already stored in scalespace
 * 
 * 
 */
void build_lowe_scalespace(struct scsp* scalespace,
                           double* image, int im_width, int im_height,
                           double sigma0){
    
    /** Do we apply a zoom in by factor 2 first */
    bool dozoomin =0;
    if(im_width == scalespace->octaves[0]->width){dozoomin = 0;}
    else if (2*im_width == scalespace->octaves[0]->width){dozoomin = 1;}
    else{printf("error\n scalespace structure not valid\n");}
   
    int nOct = scalespace->nOct;  /* # of octaves in the scalespace */
    int nSca, width, height;      /* current octave dimensions */
    double delta;                 /* current octave  inter-sample distance */
    double sig_prev,sig_next,sigma_extra;
    
    /* for dereferencing */
    struct octa* octave;
    struct octa* octave_prev;
    int width_prev, height_prev;
    double* im_prev;
    double* im_next;
    
    for(int o=0;o<nOct;o++){
        octave = scalespace->octaves[o];
        nSca   = octave->nSca;
        width  = octave->width;
        height = octave->height;
        delta  = octave->delta;  /* intersample distance */
        
        
        /** first image in the stack */
        if(o==0){ /* from input image */
            if(dozoomin){
                double* tmp = (double*)malloc(2*im_width*2*im_height*sizeof(double));
                oversample_by2_bilin(image, tmp, im_width, im_height);
                sigma_extra = sqrt(octave->sigmas[0]*octave->sigmas[0] - sigma0*sigma0)/delta; /* delta = 0.5 if zoomin is applied */
                add_gaussian_blur(tmp, octave->imStack, 2*im_width,2*im_height,sigma_extra);
                free(tmp);
            }
            else{
                sigma_extra = sqrt(octave->sigmas[0]*octave->sigmas[0] - sigma0*sigma0)/delta; /* but delta = 1 if zoomin is not applied */
                add_gaussian_blur(image, octave->imStack, im_width,im_height,sigma_extra);
            }
        }        
        else{ /* from previous octave */
            octave_prev = scalespace->octaves[o-1];
            width_prev  = octave_prev->width;
            height_prev = octave_prev->height;
            subsample_by2(&octave_prev->imStack[3*width_prev*height_prev], octave->imStack, width_prev, height_prev);  /* HARDCODED 3 for the source image*/
        }
        
        /** The rest of the image stack*/ 
        for(int s=1;s<nSca;s++){ /*add blur to previous image in the stack*/
            im_prev = &octave->imStack[(s-1)*width*height];
            im_next = &octave->imStack[s*width*height];
            
            sig_prev = octave->sigmas[s-1];
            sig_next = octave->sigmas[s];
            sigma_extra = sqrt(sig_next*sig_next- sig_prev*sig_prev)/delta;                
            add_gaussian_blur(im_prev,im_next,width,height,sigma_extra);
        }    
    }
}



/** @brief Compute the Gaussian scalespace for SIFT
 * 
 * 
 *  in @param input image  of im_width X im_height
 *     @param sigma_in : assumed level of blur in the image
 *  
 *  
 *  The construction parameters are already stored in scalespace
 */
void build_lowe_scalespace_ter(struct scsp* scalespace,
                               double* image,
                               int im_width,
                               int im_height,
                               double sigma_in){
    
    
    /* the characteristic of the sead */
    double delta_min = scalespace->octaves[0]->delta;
    double sigma_min = scalespace->octaves[0]->sigmas[0];
    int width_min    = scalespace->octaves[0]->width;
    int height_min   = scalespace->octaves[0]->height;
    
    /* checking scale-space definition consistance*/
    assert(width_min == (int)(im_width/delta_min));
    assert(height_min == (int)(im_height/delta_min));
    
    int nOct = scalespace->nOct;  /* # of octaves in the scalespace */
    int nSca, width, height;      /* current octave dimensions */
    double delta;                 /* current octave  inter-sample distance */
    double sig_prev,sig_next,sigma_extra;
    
    /* for dereferencing */
    struct octa* octave;
    struct octa* octave_prev;
    int width_prev, height_prev;
    double* im_prev;
    double* im_next;
    
    for(int o=0;o<nOct;o++){
        
        octave = scalespace->octaves[o];
        nSca   = octave->nSca;
        width  = octave->width;
        height = octave->height;
        delta  = octave->delta;  /* intersample distance */
        
        
        /** first image in the stack */
        if(o==0){ /* from input image */
            assert(sigma_min>=sigma_in);
            sigma_extra = sqrt(sigma_min*sigma_min - sigma_in*sigma_in)/delta_min;
            if(delta_min < 1){
                double* tmp = (double*)malloc(width* height * sizeof(double));
                oversample_bilin(image,im_width,im_height,tmp,width,height,delta_min);
                add_gaussian_blur(tmp,octave->imStack,width,height,sigma_extra);
                free(tmp);
            }else{ /* ie delta_min = 1, width_min = in_width... */
                add_gaussian_blur(image,octave->imStack,width,height,sigma_extra);
            }
        }        
        else{ /* from previous octave */
            octave_prev = scalespace->octaves[o-1];
            width_prev  = octave_prev->width;
            height_prev = octave_prev->height;
            subsample_by2(&octave_prev->imStack[3*width_prev*height_prev], octave->imStack, width_prev, height_prev);  /* HARDCODED 3 for the source image*/
        }
        
        
        
        /** The rest of the image stack*/ 
        for(int s=1;s<nSca;s++){ /*add blur to previous image in the stack*/
            im_prev = &octave->imStack[(s-1)*width*height];
            im_next = &octave->imStack[s*width*height];
            
            sig_prev = octave->sigmas[s-1];
            sig_next = octave->sigmas[s];
            sigma_extra = sqrt(sig_next*sig_next- sig_prev*sig_prev)/delta;                
            add_gaussian_blur(im_prev,im_next,width,height,sigma_extra);
        }    
    }
}


/** @brief Compute the Gaussian scalespace for SIFT
 * 
 *  Trick : subsample ALL auxiliary images to compute the three first images of the next octave.
 * 
 */
void build_lowe_scalespace_bis(struct scsp* scalespace, double* image, int im_width, int im_height, double sigma0){
    
    /** Do we apply a zoom in by factor 2 first */
    bool dozoomin = 0;
    if(im_width == scalespace->octaves[0]->width){dozoomin = 0;}
    else if(2*im_width == scalespace->octaves[0]->width){dozoomin = 1;}
    else{printf("error\n scalespace structure not valid\n");}
    
   
    int nOct = scalespace->nOct;  /* # of octaves in the scalespace */
    int nSca, width, height;      /* current octave dimensions */
    double delta;                 /* current octave  inter-sample distance */
    double sigma_extra;
    
    /* for dereferencing */
    struct octa* octave;
    struct octa* octave_prev;
    double* imStack;
    int width_prev, height_prev;
    
    for(int o=0;o<nOct;o++){
        octave = scalespace->octaves[o];
        nSca   = octave->nSca;
        width  = octave->width;
        height = octave->height;
        delta  = octave->delta;  /* intersample distance */
        imStack = octave->imStack;
        
        /** first octave */
        if(o==0){
            /** first image */
            if(dozoomin){
                double* tmp = (double*)malloc(2*im_width*2*im_height*sizeof(double));
                oversample_by2_bilin(image, tmp, im_width, im_height);
                sigma_extra = sqrt(octave->sigmas[0]*octave->sigmas[0] - sigma0*sigma0)/delta; /* delta = 0.5 if zoomin is applied */
                add_gaussian_blur(tmp, imStack, width, height,sigma_extra);
                free(tmp);
            }
            else{
                sigma_extra = sqrt(octave->sigmas[0]*octave->sigmas[0] - sigma0*sigma0)/delta; /* but delta = 1 if zoomin is not applied */
                add_gaussian_blur(image, imStack, width, height,sigma_extra);
            }
            /** subsequent images in the first octave */
            for (int s=1;s<nSca;s++){
                sigma_extra = sqrt(octave->sigmas[s]*octave->sigmas[s] - octave->sigmas[s-1]*octave->sigmas[s-1])/delta;
                add_gaussian_blur(&imStack[(s-1)*width*height],&imStack[s*width*height], width, height, sigma_extra);
            }
        }    
        
        /** other octaves */
        else{ /* from previous octave */
            octave_prev = scalespace->octaves[o-1];
            width_prev  = octave_prev->width;
            height_prev = octave_prev->height;
            /** Subsampling of already computed scales. */
            subsample_by2(&octave_prev->imStack[3*width_prev*height_prev], imStack, width_prev, height_prev);  /* HARDCODED 3 for the source image*/
            subsample_by2(&octave_prev->imStack[4*width_prev*height_prev], &imStack[width*height], width_prev, height_prev);  /* HARDCODED 3 for the source image*/
            subsample_by2(&octave_prev->imStack[5*width_prev*height_prev], &imStack[2*width*height], width_prev, height_prev);  /* HARDCODED 3 for the source image*/
            /** Remaining scales are computed via Gaussian convolution. */
            for (int s=3;s<nSca;s++){ //nSca= 5 (that is 6 scales, and 3 for detection)
                sigma_extra = sqrt(octave->sigmas[s]*octave->sigmas[s] - octave->sigmas[s-1]*octave->sigmas[s-1])/delta;
                add_gaussian_blur(&imStack[(s-1)*width*height],&imStack[s*width*height], width, height, sigma_extra);
            }
        }
    }
}


/** @brief difference of Gaussians 
 * 
 */
void build_dog(struct scsp *scalespace,
               struct scsp *dog){
    
    int nOct = dog->nOct;
    int nSca;
    int width;
    int height;
    
    /* dereferencing */
    struct octa* dog_oct;
    struct octa* ss_oct;
    double* diff;  /* dog at scale s */
    double* im_P;  /* image at scale (s+1) */
    double* im_M;  /* image at scale s */
    
    for(int o=0;o<nOct;o++){
        dog_oct = dog->octaves[o]; 
        ss_oct  = scalespace->octaves[o];
        
        nSca   = dog_oct->nSca;
        width  = dog_oct->width;
        height = dog_oct->height;
        
        for(int s=0;s<nSca;s++){
            /* deferencing */
            diff = &dog_oct->imStack[s*width*height];
            im_P = &ss_oct->imStack[(s+1)*width*height];
            im_M = &ss_oct->imStack[s*width*height];
            
            for(int p=0;p<width*height;p++)
                diff[p]=im_P[p]-im_M[p];
        }
    }
}


/** @brief Compute the 2d gradient of each image in the scale-space
 * 
 *  @in scalespace
 *  @out dx_scalespace,
 *       scalespace structure storing the gradient x-component of each image
 * 
 *  @out dy_scalespace,
 *       scalespace structure storing the gradient y-component of each image  
 * 
 * 
 * The gradients of the auxiliary images are not computed.
 * 
 */
void build_scalespace_gradient(struct scsp* scalespace,
                               struct scsp* dx_scalespace,
                               struct scsp* dy_scalespace){
    
    int nOct = scalespace->nOct;
    int nSca;
    int width;
    int height;
    
    /* for readability */
    double* im;    /* image */
    double* dx_im;  /* image at scale (s+1) */
    double* dy_im;  /* image at scale s */
    
    for(int o=0;o<nOct;o++){
        nSca   = scalespace->octaves[o]->nSca;
        width  = scalespace->octaves[o]->width;
        height = scalespace->octaves[o]->height;
        
        for(int s=0;s<nSca;s++){
            im    =    &scalespace->octaves[o]->imStack[s*width*height];
            dx_im = &dx_scalespace->octaves[o]->imStack[s*width*height];
            dy_im = &dy_scalespace->octaves[o]->imStack[s*width*height];
            compute_gradient(im, dx_im, dy_im, width, height);                 
        }
    }
}


                
/** @brief Scan a stack of image and save the coordinates of samples 
 * 
 * 
 * @param f    : function returning bool predicate
 *               bool(*f)(int,int,int, int,int) in prototype means that f is a pointer
 *                       to a function returning a boolean value and taking 
 *                       5 integer parameters.  
 *               REFERENCE : Kernighan-Ritchie, Chapter 5 - Section 12 -Complicated declaration . 
 * 
 * 
 * @param imStack  : stack of nSca of images with width*height pixels.
 * 
 */
void scan_imStack(double* imStack, int nSca, int width, int height, struct intTrio_lst *coords){//, bool(*f)(int,int,int, int, int)){
    
    struct intTrio *coord;
    bool isMax,isMin;
    int tmpCC=0;
    
    double* ctr;     /* pointing to the center  */ 
    double *n1,*n2,*n3,*n4,*n5,*n6,*n7,*n8,*n9,
           *n10,*n11,*n12,*n13,     *n14,*n15,*n16,*n17,
           *n18,*n19,*n20,*n21,*n22,*n23,*n24,*n25,*n26;
    
    /** Loop through the samples of the image stack (one octave) */
    for(int s=1;s<nSca-1;s++){
        for(int i=1;i<height-1;i++){
            for(int j=1;j<width-1;j++){
                
                isMax=1;
                isMin=1;
                ctr = &imStack[s*width*height+i*width+j];
                
                /** Compare to the 26 neighbors ******************************/ 
                /*  ---- image below ----------- */     /* ---- image middle -----*/   /* ---- image up -------------  */
                n1  = ctr -width*height-width-1;        n10 = ctr -width-1;           n18 = ctr +width*height-width-1; 
                n2  = ctr -width*height-width;          n11 = ctr -width;             n19 = ctr +width*height-width;   
                n3  = ctr -width*height-width+1;        n12 = ctr -width+1;           n20 = ctr +width*height-width+1; 
                n4  = ctr -width*height-1;              n13 = ctr -1;                 n21 = ctr +width*height-1;       
                n5  = ctr -width*height;                                              n22 = ctr +width*height;         
                n6  = ctr -width*height+1;              n14 = ctr +1;                 n23 = ctr +width*height+1;       
                n7  = ctr -width*height+width-1;        n15 = ctr +width-1;           n24 = ctr +width*height+width-1; 
                n8  = ctr -width*height+width;          n16 = ctr +width;             n25 = ctr +width*height+width;   
                n9  = ctr -width*height+width+1;        n17 = ctr +width+1;           n26 = ctr +width*height+width+1; 
                                                    

                isMin =   (*ctr<*n1 )&&(*ctr<*n2 )&&(*ctr<*n3 )
                        &&(*ctr<*n4 )&&(*ctr<*n5 )&&(*ctr<*n6 )
                        &&(*ctr<*n7 )&&(*ctr<*n8 )&&(*ctr<*n9 )
                        &&(*ctr<*n10)&&(*ctr<*n11)&&(*ctr<*n12)
                        &&(*ctr<*n13)             &&(*ctr<*n14)
                        &&(*ctr<*n15)&&(*ctr<*n16)&&(*ctr<*n17)
                        &&(*ctr<*n18)&&(*ctr<*n19)&&(*ctr<*n20)
                        &&(*ctr<*n21)&&(*ctr<*n22)&&(*ctr<*n23)
                        &&(*ctr<*n24)&&(*ctr<*n25)&&(*ctr<*n26);
            
                isMax =   (*ctr>*n1 )&&(*ctr>*n2 )&&(*ctr>*n3 )
                        &&(*ctr>*n4 )&&(*ctr>*n5 )&&(*ctr>*n6 )
                        &&(*ctr>*n7 )&&(*ctr>*n8 )&&(*ctr>*n9 )
                        &&(*ctr>*n10)&&(*ctr>*n11)&&(*ctr>*n12)
                        &&(*ctr>*n13)             &&(*ctr>*n14)
                        &&(*ctr>*n15)&&(*ctr>*n16)&&(*ctr>*n17)
                        &&(*ctr>*n18)&&(*ctr>*n19)&&(*ctr>*n20)
                        &&(*ctr>*n21)&&(*ctr>*n22)&&(*ctr>*n23)
                        &&(*ctr>*n24)&&(*ctr>*n25)&&(*ctr>*n26);
                
           
                if(isMax||isMin){ /*if 3d discrete extrema, save a candidate keypoint*/
                    tmpCC +=1;
                    coord = (struct intTrio*)malloc(sizeof(struct intTrio));
                    coord->i=i;
                    coord->j=j;
                    coord->s=s;
                    add_intTrio_to_lst(coord,coords);
                }
            }
        }
    }
}



/** @brief Extract discrete extrema from DoG scale-space.
 * 
 * in  @param dog  : Difference of Gausssian (scale-space structure)
 * 
 * out @param keys : List of keypoints.
 * 
 * 
 * 
 *    The following parameters are for memory allocation only.
 * 
 * @param n_bins : Number of bins in the histogram
 *                 used to attribute principal orientation.
 * 
 * @param n_hist
 * @param n_ori  : The SIFT descriptor is an array of
 *                 n_hist X n_hist weighted orientation
 *                 with n_ori bins each.
 * 
 */
void find_candidate_keypoints(struct scsp* dog,
                              struct kypt_lst* keys,
                              int n_ori,
                              int n_hist,
                              int n_bins){
    
    int nOct = dog->nOct;      /* # of octaves in the scale-space */
    int nSca, width, height;   /* dimension of the image stack in
                                  the current octave */
    double delta;              /* inter-sample for the current octave */
    
    /* Dereferencing pointers */
    struct octa* octave;
    double* imStack;
    struct intTrio_lst* coords;
    struct intTrio *coord;
    struct kypt *key;
    
    for(int o=0;o<nOct;o++){
        octave  = dog->octaves[o];                                                                        
        nSca    = octave->nSca;
        width   = octave->width;
        height  = octave->height;
        delta   = octave->delta;  /* intersample distance */
        imStack = octave->imStack;
        coords = malloc_intTrio_lst();
        
        /** Scan the image stack for 3d discrete extrema *********************/
        scan_imStack(imStack, nSca, width, height, coords);
        
        /** Build keypoints struct for each discrete extrema *****************/   
        for(int k=0;k<coords->card;k++){
            /* A discrete extremum ... */
            coord = coords->lst[k];
            /* ...produce a keypoint candidate ... */
            key = malloc_kypt_extra(n_ori, n_hist, n_bins, 5);
            
            /*  ...with the following position ...*/
            key->i = coord->i;
            key->j = coord->j;
            key->s = coord->s;
            key->o = o;
            key->x = delta*key->i;
            key->y = delta*key->j;
            key->sigma = octave->sigmas[coord->s];
            /*  ... and DoG value. */
            key->val   = imStack[coord->s*width*height+coord->i*width+coord->j];
            add_kypt_to_lst(key,keys);
        }
        free_intTrio_lst(coords);
    }
}




/** @brief Filter a list of keypoints
 *
 * in  @param keysIn : list of keys to be filtered.
 * 
 * out @param keysAccept   list of keys with DoG absolute value beyond threshold.
 * 
 *  @param thresh    :constant threshold over the entire scale-space.
 * 
 */
void discard_detections_duetonoise(struct kypt_lst *keysIn,
                                   struct kypt_lst *keysAccept,
                                   struct kypt_lst *keysReject,
                                   double thresh){
 
    struct kypt *key;
    struct kypt *keycp;
    
    for(int k=0;k<keysIn->card;k++){
        key = keysIn->lst[k];
        keycp = malloc_kypt_from_model_and_copy(key);
        if(fabs(key->val)> thresh){
            add_kypt_to_lst(keycp,keysAccept);
        }
        else{
            add_kypt_to_lst(keycp,keysReject);
        }
    }
}





/** @brief Execute one interpolation (to refine extrema position
 * 
 * in  @param imStck      : a stack of images in the DoG space
 *                    width X height X nscales samples.
 *                   
 * in  @param i , j ,  s  : 3d discrete extrema
 * 
 * out @param di , dj , ds : offset in each spatial direction
 * out @param dVal         , the extremum value is : imStck[i,j,s]+dVal
 * 
 * 
 * Computes the 3D Hessian of the DoG space in one point
 * 
 * 
 */
void inverse_3D_Taylor_second_order_expansion(double *imStck,
                                              int width, int height, int nscales,
                                              int i, int j, int s,
                                              double *di, double *dj, double *ds, double *dVal){
    
    double hXX,hXY,hXS,hYY,hYS,hSS;
    double det,aa,ab,ac,bb,bc,cc;
    double gX,gY,gS;
    double ofstX, ofstY, ofstS, ofstVal;
     
    
    /** Compute the 3d Hessian at pixel (i,j,s)  Finite difference scheme  *****/
    hXX = imStck[s*width*height+(i-1)*width+j] + imStck[s*width*height+(i+1)*width+j] - 2*imStck[s*width*height+i*width+j];
    hYY = imStck[s*width*height+i*width+(j+1)] + imStck[s*width*height+i*width+(j-1)] - 2*imStck[s*width*height+i*width+j];
    hSS = imStck[(s+1)*width*height+i*width+j] + imStck[(s-1)*width*height+i*width+j] - 2*imStck[s*width*height+i*width+j];
    hXY = 0.25*(  (imStck[s*width*height+(i+1)*width+(j+1)] - imStck[s*width*height+(i+1)*width+(j-1)])
                - (imStck[s*width*height+(i-1)*width+(j+1)] - imStck[s*width*height+(i-1)*width+(j-1)]) );
    hXS = 0.25*(  (imStck[(s+1)*width*height+(i+1)*width+j] - imStck[(s+1)*width*height+(i-1)*width+j])
                - (imStck[(s-1)*width*height+(i+1)*width+j] - imStck[(s-1)*width*height+(i-1)*width+j]) );
    hYS = 0.25*(  (imStck[(s+1)*width*height+i*width+(j+1)] - imStck[(s+1)*width*height+i*width+(j-1)])
                - (imStck[(s-1)*width*height+i*width+(j+1)] - imStck[(s-1)*width*height+i*width+(j-1)]) );
    
    
    /** Compute the 3d gradient at pixel (i,j,s) */
    gX = 0.5*( imStck[s*width*height+(i+1)*width+j] - imStck[s*width*height+(i-1)*width+j] );
    gY = 0.5*( imStck[s*width*height+i*width+(j+1)] - imStck[s*width*height+i*width+(j-1)] );
    gS = 0.5*( imStck[(s+1)*width*height+i*width+j] - imStck[(s-1)*width*height+i*width+j] );
    
    /** Inverse the Hessian - Fitting a quadric function */
    det = hXX*hYY*hSS - hXX*hYS*hYS - hXY*hXY*hSS + 2*hXY*hXS*hYS - hXS*hXS*hYY ; 
    aa = (hYY*hSS - hYS*hYS)/det;
    ab = (hXS*hYS - hXY*hSS)/det;
    ac = (hXY*hYS - hXS*hYY)/det;
    bb = (hXX*hSS - hXS*hXS)/det;
    bc = (hXY*hXS - hXX*hYS)/det;
    cc = (hXX*hYY - hXY*hXY)/det;
    
    /** Compute position offset */
    ofstX = -aa*gX-ab*gY-ac*gS;
    ofstY = -ab*gX-bb*gY-bc*gS;
    ofstS = -ac*gX-bc*gY-cc*gS;
    /** Compute the DoG value offset */
    ofstVal =  0.5*(gX*ofstX+gY*ofstY+gS*ofstS);

    /** output */
    *di = ofstX;
    *dj = ofstY;
    *ds = ofstS;
    *dVal = ofstVal;
}




/** @brief Refine the position of candidate keypoints
 *         
 *  in  @param dog : difference of Gaussian
 *  in  @param keys : list candidate keypoints (3d discrete extrema)
 * 
 *  out @param keysInterpol : list of interpolated keypoints.
 *  out @param keysReject   : list of discarded keypoints.
 * 
 *  The interpolation model consists in a local 3d quadratic model of the DoG space.
 *  
 *  An interpolation is successful if the interpolated extremum lays inside the pixel area
 *  
 *  Iterative process 
 * 
 * 
 */
void interpolate_position_of_candidate_keypoints(struct scsp *dog,
                                                 struct kypt_lst *keys,
                                                 struct kypt_lst *keysInterpol,
                                                 struct kypt_lst *keysReject){
    
    
    int nMaxIntrp = 5; /* Maximum number of consecutive unsuccessful interpolation */
    int nIntrp;
    
    int o,s;
    int i,j;
    int ic,jc;  /* current pixel index - at each interpolation */
    // TODO sc; /* current scale index */
    int width,height,nscales;
    double delta;
    struct kypt *key;
    struct kypt *keycp;
    struct octa *octave;
    double* imStck;
    double ofstX, ofstY,ofstS;
    double val;
    double sigmaratio;    /* Ratio between two consecutive scales in the scalespace */
    
    bool isConv;
 
    sigmaratio = dog->octaves[0]->sigmas[1]/dog->octaves[0]->sigmas[0]; /* assuming the ratio is constant over all scales and over all octaves TEST*/ 
    
    for(int k=0;k<keys->card;k++){

        
        /* Loading keypoint and associated octave */
        key = keys->lst[k];
        o = key->o;
        s = key->s;
        i = key->i;
        j = key->j;
        octave = dog->octaves[o];
        width  = octave->width;
        height = octave->height;
        nscales = octave->nSca;
        delta  = octave->delta;
        imStck = octave->imStack;
        
        val = key->val;
        
        ic=i;   /* current value of i coordinate */
        jc=j;
        nIntrp=0;
        isConv=0;
        
        /** Discard detection too close to image borders */
        
        while((nIntrp<nMaxIntrp)){
            
            /** Extrema interpolation via a quadratic function */
            /*   only if the detection is not too close to the border (so the discrete 3D Hessian is well defined) */
              //printf(" %i Avant inverse 3D Taylor expansion (ic,jc)=(%i,%i)  (height,width)=(%i,%i)     - nombre Interpolation %i \n",k,ic,jc,height,width,nIntrp);
            if((0<ic)&&(ic<height-1)&&(0<jc)&&(jc<width-1)){
                inverse_3D_Taylor_second_order_expansion(imStck, width, height, nscales,ic, jc, s, &ofstX, &ofstY, &ofstS, &val);
            }else{
                isConv = 0;
                ofstX = 5.0;
                ofstY = 5.0;
            }
            
            /** Test if the quadric model is consistant */
            if(fabs(ofstX)<0.6 && fabs(ofstY)<0.6){
                isConv=1;
                break;
            }else{
                if((ofstX>+0.6)&&(ic+1<height-1)) {ic +=1;}
                if((ofstX<-0.6)&&(ic-1>0))        {ic -=1;}
                if((ofstY>+0.6)&&(jc+1<width-1))  {jc +=1;}
                if((ofstY<-0.6)&&(jc-1<0))        {jc -=1;}              
            }
            nIntrp +=1;
        }
        
        /** Create key and save in corresponding kypt structure */
        keycp = malloc_kypt_from_model_and_copy(key);
        keycp->x = (ic+ofstX)*delta;
        keycp->y = (jc+ofstY)*delta;
        keycp->i = ic;
        keycp->j = jc;
        keycp->sigma = octave->sigmas[s]*pow(sigmaratio,ofstS); /* logarithmic scale */
        
        if(isConv){
            add_kypt_to_lst(keycp,keysInterpol);
        }else{
            add_kypt_to_lst(keycp,keysReject);
        }
    }   
}





/** @brief  Compute Edge response 
 *    i.e.  Compute the ratio of principal curvatures
 *    i.e.  Compute the ratio (hXX + hYY)*(hXX + hYY)/(hXX*hYY - hXY*hXY);
 *          
 * 
 *     The 2D hessian of the DoG operator is computed via finite difference schemes.
 * 
 *    in  @param dog  : difference of Gaussians
 *    out @param keys : candidate keypoints
 * 
 * Note:
 *  - No keypoint is discarded here
 *
 */
void compute_HarrisStephen_edge_response(struct scsp *dog, struct kypt_lst *keys){
    
    int o,s,i,j, width, height;
    struct kypt *key;
    struct octa *octave;
    double* image;
    double hXX,hXY,hYY;
    double hsEdgeResp;  /* Harris and Stephen computed on the DoG operator */
    
    for(int k=0;k<keys->card;k++){
        
        /* Loading keypoint and associated octave */
        key = keys->lst[k];
        o = key->o;
        s = key->s;
        i = key->i;
        j = key->j;
        octave = dog->octaves[o];
        width  = octave->width;
        height = octave->height;
        image  = &octave->imStack[s*width*height];
 
        /* Compute the 2d Hessian at pixel (i,j) */
        hXX = image[(i-1)*width+j]+image[(i+1)*width+j]-2*image[i*width+j];
        hYY = image[i*width+(j+1)]+image[i*width+(j-1)]-2*image[i*width+j] ;
        hXY = 1./4*((image[(i+1)*width+(j+1)]-image[(i+1)*width+(j-1)])-(image[(i-1)*width+(j+1)]-image[(i-1)*width+(j-1)]));
        
        /* Harris and Stephen Edge response */
        hsEdgeResp = (hXX + hYY)*(hXX + hYY)/(hXX*hYY - hXY*hXY);
        key->hsEdgeResp = hsEdgeResp;
    
    }    
}




/** @brief Discard keys with edge response
 * 
 *  in  @param keysIn : list of keypoints, the edge response is stored
 * 
 *  out @param keysAccepted : passing
 *  out @param keysRejected : failling
 * 
 *  @param threshold on (hXX + hYY)*(hXX + hYY)/(hXX*hYY - hXY*hXY)
 *                   on the ratio of principal curvatures.
 * 
 * 
 */
void discard_detections_onedge(struct kypt_lst *keysIn,
                               struct kypt_lst *keysAccept,
                               struct kypt_lst *keysReject,
                               double thresh){
    
    struct kypt *key;
    struct kypt *keycp;
    for(int k=0;k<keysIn->card;k++){
        key = keysIn->lst[k];
        keycp = malloc_kypt_from_model_and_copy(key);
        if(fabs(key->hsEdgeResp)<=thresh){
            add_kypt_to_lst(keycp,keysAccept);
        }
        else{
            add_kypt_to_lst(keycp,keysReject);
        }
    }    
}



/** @brief Attribute a reference orientation to each keypoint of a list
 * 
 *  
 *   in  @param dx_scalespace, x component (up-bottom)
 *   in  @param dy_scalespace, y component (left-right) 
 * 
 *   in  @param keysIn     list of keypoints (pointer)
 * 
 *   out @param keysOut   list of oriented keypoints
 *
 *            card(keysIn) <= card(keysOut)          
 *
 *   @param t  : threshold over a local maxima to constitute a reference orientation 
 * 
 *   @param lambda_ori : size parameter for the Gaussian window
 * 
 * 
 * 
 */
void attribute_orientation_to_keypoints(struct scsp *dx_scalespace,
                                        struct scsp *dy_scalespace, /*gradient scalespaces*/
                                        struct kypt_lst *keysIn,
                                        struct kypt_lst *keysOut,
                                        struct kypt_lst *cp_keys_multiori,
                                        int n_bins,
                                        double lambda_ori,double t
                                        ){
    
    
    struct kypt *key;
    struct kypt *keycp;
    struct kypt *keycpcp;
    struct octa *dx_octave;
    struct octa *dy_octave;
    double* dx_image;
    double* dy_image;
    
    int width, height;
    double delta;
    int o,s;
    double x,y,sigma;
    
    int count;
    int n_prOri;  
    
    for(int k=0;k<keysIn->card;k++){
        
        count=0;
        
        /** Loading keypoint gradient scalespaces */
        key = keysIn->lst[k];
        x = key->x;
        y = key->y;
        o = key->o;
        s = key->s;
        sigma = key->sigma;
        dx_octave = dx_scalespace->octaves[o];
        dy_octave = dy_scalespace->octaves[o];
        width  = dx_octave->width;
        height = dx_octave->height;
        delta  = dx_octave->delta;
        dx_image = &dx_octave->imStack[s*width*height];
        dy_image = &dy_octave->imStack[s*width*height];
       
        /** Accumulate gradient orientation histogram */       
        count = accumulate_orientation_histogram_bis(x/delta,y/delta,sigma/delta, dx_image, dy_image, width, height, n_bins,lambda_ori,key->hist_prOri);
        
        /** Extract principal orientation */
        n_prOri = extract_principal_orientations(key->hist_prOri,key->n_bins, t, 1, key->ori_prOri); /*0.8 definition of a secondary orrientation */
        
        /** Updating keypoints and save in new list */
        for(int n=0;n<n_prOri;n++){
            keycp = malloc_kypt_from_model_and_copy(key);
            keycp->theta = key->ori_prOri[n];
            keycp->nsamples_ori = count;

            add_kypt_to_lst(keycp, keysOut);
            if(n_prOri>1){ /* a copy for keypoints with multiple orientations */
                keycpcp = malloc_kypt_from_model_and_copy(keycp);
                add_kypt_to_lst(keycpcp,cp_keys_multiori);
            }
        }
    }
}


/** @brief Attribute a feature vector to each keypoint of a list
 * 
 *  
 *   in  @param dx_scalespace, x component (up to bottom)
 *   in  @param dy_scalespace, y component (left to right) 
 * 
 *   in  @param keys     list of keypoints to be described
 *         
 *   @param n_hist   , the descriptor is an array of nhist^2 weighted histograms.
 *   @param n_ori    , each weighted histogram has n_ori bins.
 * 
 *   @param lambda_descr : size parameter for the Gaussian window
 *                         - the Gaussian window has a std dev of lambda_descr X sigma
 *                         - the patch considered has a width of 
 *                         2*(1+1/n_hist)*lambda_descr X sigma
 */
void attribute_featurevector_to_keypoints(struct scsp *dx_scalespace,
                                          struct scsp *dy_scalespace, /*gradient scalespaces*/
                                          struct kypt_lst *keys,
                                          int n_hist,
                                          int n_ori,
                                          double lambda_descr ){
    

    
    struct kypt *key;
    struct octa *dx_octave;
    struct octa *dy_octave;
    double* dx_image;
    double* dy_image;
    int width, height;
    double delta, theta;
    int o,s;
    double x,y,sigma;
    
    int count;
    
    for(int k=0;k<keys->card;k++){
        
        
        /** Loading keypoint gradient scalespaces */
        key = keys->lst[k];
        x = key->x;
        y = key->y;
        o = key->o;
        s = key->s;
        sigma = key->sigma;
        theta = key->theta;
        dx_octave = dx_scalespace->octaves[o];
        dy_octave = dy_scalespace->octaves[o];
        width  = dx_octave->width;
        height = dx_octave->height;
        delta  = dx_octave->delta;
        dx_image = &dx_octave->imStack[s*width*height];
        dy_image = &dy_octave->imStack[s*width*height];
        
        
        /** Compute descriptor representation */
        count = extract_sift_feature_vector_bis(x/delta,
                                                y/delta,         /* relative coordinates */
                                                sigma/delta,
                                                theta,
                                                dx_image,
                                                dy_image,        /* pointers to precomputed gradient images */
                                                width, height,       
                                                n_hist,
                                                n_ori,
                                                lambda_descr,    /* design parameters */
                                                key->descr);     /* output feature vector */ 

        
        
        /** Threshold and quantization of the descriptor */
        threshold_and_quantize_sift_descriptor(key->descr,
                                               n_hist,
                                               n_ori,
                                               1,      /* bool (Do we apply threshold ?) */
                                               0.2);   /* threshold value      */
        key->nsamples_descr = count;

    }
}
