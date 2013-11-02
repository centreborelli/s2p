/*
 *  libraryBasic.cpp
 *  iipSources
 *
 *  Created by Antoni Buades on 28/02/2012.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "libraryBasic.h"




namespace libIIPStable {
	
	
    
    
    
    //
    //! Value operations
    //
    
    void fpClear(float *fpI,float fValue, int iLength)
    {
        assert(iLength > 0);
        for (int ii=0; ii < iLength; ii++) fpI[ii] = fValue;	
    }
    
    
    void fpCopy(float *fpI,float *fpO, int iLength)
    {
        assert(iLength > 0);
        if (fpI != fpO)  memcpy((void *) fpO, (const void *) fpI, iLength * sizeof(float));
    }
    
    
    void dpClear(double *dpI,double dValue, int iLength)
    {
        assert(iLength > 0);
        for (int ii=0; ii < iLength; ii++) dpI[ii] = dValue;	
    }
    
    
    void dpCopy(float *dpI,float *dpO, int iLength)
    {
        assert(iLength > 0);
        if (dpI != dpO)  memcpy((void *) dpO, (const void *) dpI, iLength * sizeof(double));
    }
    
    
    
    float fpMax(float *u,int *pos, int size)
    {  
        assert(size > 0);
        float max=u[0];
        if (pos) *pos=0;  
        for(int i=1; i<size; i++)  if(u[i]>max){ max=u[i]; if (pos) *pos=i; } 
        return max;
    }
    
    
    float fpMin(float *u,int *pos,int size) 
    {
        assert(size > 0);
        float min=u[0];  
        if (pos) *pos=0;
        for(int i=1;i<size;i++)	if(u[i]<min){ min=u[i]; if (pos) *pos=i; }
        return min;
    }
    
    
    
    double dpMax(double *u,int *pos, int size)
    {  
        double max=u[0];
        if (pos) *pos=0;  
        for(int i=1; i<size; i++)  if(u[i]>max){ max=u[i]; if (pos) *pos=i; } 
        return max;
    }
    
    
    double dpMin(double *u,int *pos,int size) 
    {
        assert(size > 0);
        double min=u[0];  
        if (pos) *pos=0;
        for(int i=1;i<size;i++)	if(u[i]<min){ min=u[i]; if (pos) *pos=i; }
        return min;
    }
    
    float  fpMean(float *u,int size)
    {
        assert(size > 0);
        float *ptru=&u[0];
        float mean=0.0;
        for(int i=0; i<size; i++,ptru++)  mean +=  *ptru;
        mean/=(float) size;
        return mean;
    }

    
    float  fpMedian(float *u,int size)
    {
        assert(size > 0);
        if (size == 1) return u[0];
        
        fpQuickSort(u, size);
        
        int dsize = (size - 1) / 2;
        if (size % 2 == 1) return u[dsize];
        else return 0.5 * (u[dsize]+u[dsize+1]);
    }
    
    
    
    double  dpMean(double *u,int size)
    {
        assert(size > 0);
        double *ptru=&u[0];
        double mean=0.0;
        for(int i=0; i<size; i++,ptru++)  mean +=  *ptru;
        mean/=(double) size;
        return mean;
    }
    
    
    float  fpVar(float *u,int size)
    {
        
        assert(size > 0);
        float *ptru=&u[0];
        float mean=0.0;
        float mean2 = 0.0;
        for(int i=0;i<size;i++,ptru++) { mean +=  *ptru; mean2 +=  *ptru *  (*ptru);}
        
        mean/=(float) size; mean2/=(float) size;
        float var = mean2- mean*mean;
        
        var = fabsf(var);
        return var;
    }
    
    
    
    double  dpVar(double *u,int size)
    {
        
        assert(size > 0);
        double *ptru=&u[0];
        double mean=0.0;
        double mean2 = 0.0;
        for(int i=0;i<size;i++,ptru++) { mean +=  *ptru; mean2 +=  *ptru *  (*ptru);}
        
        mean/=(float) size; mean2/=(float) size;
        double var = mean2- mean*mean;
        
        var = fabs(var);
        return var;
    }
    
    
    
    void fpCombine(float *u,float a,float *v,float b, float *w,  int size)
    {
        for(int i=0;i<size ;i++)   w[i]= (a *  u[i] + b* v[i]);  
        
    }
    
    
    
    void dpCombine(double *u,double a,double *v,double b, double *w,  int size)
    {
        for(int i=0;i<size ;i++)   w[i]= (a *  u[i] + b* v[i]);  
        
    }
    
    
    
    
    void fiImageDrawCircle(float *igray, int pi,int pj, double radius, float value, int width, int height)	
    {
        
        int mark = (int) rint(radius);
        double radius2 = radius * radius;
        
        for(int s = -mark ; s <= mark ;s++)
            for(int r = -mark ; r <= mark ;r++)
                if (pj+s>=0 && pi+r>= 0 && pj+s < height && pi+r < width && (double) (r*r + s*s) < radius2)
                {	
                    
                    int l = (pj+s)*width+pi+r;
                    igray[l] = value;
                }
        
    }
    
    
    void fiImageDrawLine(float *igray, int a0, int b0, int a1, int b1, float value, int width, int height)
    {
        
        int bdx,bdy;
        int sx,sy,dx,dy,x,y,z,l;
        
        bdx = width;
        bdy = height;
        
        if (a0 < 0) a0=0; 
        else if (a0>=bdx) a0=bdx-1;
        
        if (a1<0)  a1=0; 
        else  if (a1>=bdx)   a1=bdx-1;
        
        if (b0<0) b0=0; 
        else if (b0>=bdy) b0=bdy-1;
        
        if (b1<0) 	b1=0; 
        else if (b1>=bdy) b1=bdy-1; 
        
        if (a0<a1) { sx = 1; dx = a1-a0; } else { sx = -1; dx = a0-a1; }
        if (b0<b1) { sy = 1; dy = b1-b0; } else { sy = -1; dy = b0-b1; }
        x=0; y=0;
        
        if (dx>=dy) 
        {
            z = (-dx) / 2;
            while (abs(x) <= dx) 
            {
                
                l =  (y+b0)*bdx+x+a0; 
                
                igray[l] = value;
                
                x+=sx;
                z+=dy;
                if (z>0) { y+=sy; z-=dx; }
                
            } 
            
        }
        else 
        {
            z = (-dy) / 2;
            while (abs(y) <= dy) {
                
                l = (y+b0)*bdx+x+a0;
                igray[l] = value;
                
                y+=sy;
                z+=dx;
                if (z>0) { x+=sx; z-=dy; }
            }
        }
        
    }
    
    
    
    
    void fpBinarize(float *u, float *v, float value, int inverse, int size)
	{ 
		for(int i=0;i<size;i++){
			
			if (u[i] >= value && !inverse) 	v[i]= 1.0;
			else if (u[i] <= value && inverse)  v[i]= 1.0;
			else v[i]= 0.0;
		}
	}
    
    
    
    
    float fpDistLp(float *u, float *v, int p, int size)
    {
        
        float fDist = 0.0f;
        
        for (int ii=0; ii < size; ii++)
        {
            
            float dif = fabsf(u[ii] - v[ii]);
            
            
            if (p == 0) fDist = MAX(fDist, dif);
            else if (p == 1)  fDist += dif;
            else if (p == 2)
                fDist += dif*dif;
            else
            {
                fDist += powf(dif, (float) p);
                
            }
            
        }
        
        fDist /= (float) size;
        
        if (p>0)
            fDist = powf(fDist, 1.0 / (float) p);
        
        
        return fDist;
        
    }
    
    
    
    float fpDistLp(float *u, float *v, float *m, int p, int size)
    {
        
        float fDist = 0.0f;
        int iCount = 0;
        
        for (int ii=0; ii < size; ii++)
            if (m[ii] > 0.0f)
            {
                
                float dif = fabsf(u[ii] - v[ii]);
                
                
                if (p == 0) fDist = MAX(fDist, dif);
                else if (p == 1)  fDist += dif;
                else if (p == 2)
                    fDist += dif*dif;
                else
                {
                    fDist += powf(dif, (float) p);
                    
                }
                
                iCount++;
            }
        
        fDist /= (float) iCount;
        
        if (p>0)
            fDist = powf(fDist, 1.0 / (float) p);
        
        
        return fDist;
        
    }
    

    
    
    
    
    //
    //! Gradient based
    //
    
    void fiComputeImageGradient(float * fpI,float *fpXgrad, float *fpYgrad, float *fpGrad, float * fpOri, int iWidth, int iHeight, char cType)
    {
        
        assert(fpI != NULL);
        assert(cType == 'c' || cType == 'f');
        
        int iC, iN, iS, iW, iE;
        float xgrad, ygrad;
        
        for (int ih = 0; ih < iHeight; ih++)
            for (int iw = 0; iw < iWidth; iw++) 
            {
                
                //! Indexos
                iC = ih * iWidth + iw;   
                iN = iC - iWidth;   				
                iS = iC + iWidth; 
                iW = iC - 1;  
                iE = iC + 1;  
                
                //! Boundary
                if (ih == 0) iN = iC;
                if (ih == iHeight-1) iS = iC;		
                if (iw == 0) iW = iC;		
                if (iw == iWidth-1) iE = iC;
                
                
                
                if (cType == 'f') 	
                {	
                    //! Forward
                    xgrad = fpI[iE] - fpI[iC]; 
                    
                    ygrad = fpI[iS] - fpI[iC] ; 
                    
                } else
                {
                    
                    //! Centered
                    xgrad = fpI[iE] - fpI[iW];
                    
                    ygrad = fpI[iS] - fpI[iN];
                    
                }
                
                
                if (fpXgrad) fpXgrad[ih * iWidth + iw] =  xgrad;
                
                if (fpYgrad) fpYgrad[ih * iWidth + iw] =  ygrad;
                
                if (fpGrad) fpGrad[ih * iWidth + iw] =  sqrtf(xgrad * xgrad + ygrad * ygrad);
                
                if (fpOri) fpOri[ih * iWidth + iw] =  atan2f(-ygrad,xgrad);
                
            }
        
        
    }
    
    
    
    
    void fiComputeImageGradient(float * fpI, float *fpGrad, float * fpOri, int iWidth, int iHeight, char cType)
    {
        fiComputeImageGradient( fpI, NULL, NULL, fpGrad,  fpOri, iWidth, iHeight, cType);		
    }
    
    void fiComputeImageGradient(float * fpI, float *fpGrad, int iWidth, int iHeight, char cType)
    {
        fiComputeImageGradient( fpI, NULL, NULL, fpGrad, NULL , iWidth, iHeight, cType);
        
    }
    
    
    
    
    
    
    //
    //! Noise
    //
    
    void fpAddNoiseGaussian(float *u, float *v, float std, long int randinit, int size)   
    {
        
        srand48( (long int) time (NULL) + (long int) getpid()  + (long int) randinit);
        
        for (int i=0; i< size; i++) 
        {
            
            double a = drand48();		
            double b = drand48();
            double z = (double)(std)*sqrt(-2.0*log(a))*cos(2.0*M_PI*b);
            
            v[i] =  u[i] + (float) z;
            
        }		
        
    }
    
    
    
    void fpAddNoiseGaussianAfine(float *u,float *v, float a,float b,long int randinit, int size)
    {
        
        srand48( (long int) time (NULL) + (long int) getpid()  + (long int) randinit);
        
        // Gaussian noise  
        for (int i=0; i< size; i++) 
        {
            
            double std = (double) (a + b * u[i]);
            std = sqrt(std);
            
            double a0 = drand48();
            double b0 = drand48();
            double z = (double)(std)*sqrt(-2.0*log(a0))*cos(2.0*M_PI*b0);
            
            v[i] = u[i] + (float) z;
            
        }		
        
    }
    
    
    
    
    //
    //! Float pointer ordering
    //
    
    int order_float_increasing(const void *a, const void *b)
    {
        if ( *(float*)a  > *(float*)b ) return 1;
        else if ( *(float*)a  < *(float*)b ) return -1;
        
        return 0;
    }
    
    
    
    
    int order_float_decreasing(const void *a, const void *b)
    {
        if ( *(float*)a  > *(float*)b ) return -1;
        else if ( *(float*)a  < *(float*)b ) return 1;
        
        return 0;
    }
    
    
    
    
    
    
    void fpQuickSort(float *fpI, int iLength, int inverse)
    {
        
        if (inverse)
            qsort(fpI, iLength, sizeof(float), order_float_decreasing);
        else
            qsort(fpI, iLength, sizeof(float), order_float_increasing);
        
    }
    
    
    
    
    
    struct stf_qsort
    {
        float value;
        float index;
    };
    
    
    
    int order_stf_qsort_increasing(const void *pVoid1, const void *pVoid2)
    {
        struct stf_qsort *p1, *p2;
        
        p1=(struct stf_qsort *) pVoid1;
        p2=(struct stf_qsort *) pVoid2;
        
        if (p1->value < p2->value) return -1;
        if (p1->value > p2->value) return  1;
        
        return 0;
    }
    
    
    
    
    
    int order_stf_qsort_decreasing(const void *pVoid1, const void *pVoid2)
    {
        struct stf_qsort *p1, *p2;
        
        p1=(struct stf_qsort *) pVoid1;
        p2=(struct stf_qsort *) pVoid2;
        
        if (p1->value < p2->value) return 1;
        if (p1->value > p2->value) return  -1;
        
        return 0;
        
    }
    
    
    
    
    
    void fpQuickSort(float *fpI, float *fpO, int iLength, int inverse)
    {
        
        struct stf_qsort *vector = new stf_qsort[iLength];
        
        for (int i=0; i < iLength; i++)
        {
            vector[i].value = fpI[i];
            vector[i].index = fpO[i];
            
        }
        
        
        if (inverse)
            qsort(vector, iLength, sizeof(stf_qsort), order_stf_qsort_decreasing);
        else
            qsort(vector, iLength, sizeof(stf_qsort), order_stf_qsort_increasing);
        
        
        for (int i=0; i < iLength; i++)
        {
            fpI[i] = vector[i].value;
            fpO[i] = vector[i].index;
            
        }
        
        
        delete[] vector;
    }
    
    
    
    
    
    
    
    
    //
    //! Histogram related
    //
    
    
    
    float* fpHisto(float* input, float *iminim, float *imaxim, int *n, float *s, int size, char flag)
    {
        
        assert(flag == 's' || flag == 'n');
        
        float minim;
        if (iminim) minim = *iminim;
        else minim = fpMin(input, NULL, size);
        
        float maxim;
        if (imaxim) maxim = *imaxim;
        else maxim = fpMax(input, NULL, size);
        
        int num;
        float step;
        if (flag == 'n')
        {
            num = *n;
            step = (maxim-minim)/ (float)num;
            *s = step;
            
        } else
        {
            step = *s;
            num = (int)(0.5+(maxim-minim)/step);
            *n = num;
        }
        
        float *histo = new float[num];
        fpClear(histo,0.0,num);
        
        
        for(int i=0; i < size; i++)
        {
            
            int cell = (int) floorf((input[i]-minim) / step);
            
            if (cell < 0) cell = 0;
            if (cell >= num) cell = num - 1;
            
            histo[cell]++;
        }
        
        
        return histo;
        
    }
    
    
    
    
    void fk_apply_histo(float *Histo,float *HistoInverse,int iRangeMax, float *src,float *srt, int width, int height)
    {
        
        
        
        int icounter=0;
        for (int adr=0; adr < width * height ; adr++) if ( (int) rintf(src[adr]) >= 0) icounter++;
        
        
        for (int adr=0; adr < width * height ; adr++)
        {
            
            int it = (int) rintf(src[adr]);
            //if (it < 0) it=0;
            
            if (it>=0)
            {
                
                if (it >= iRangeMax) it = iRangeMax -1;
                
                
                float x = Histo[it] * HSP_NB_CASES / (float) icounter; 
                int k= (int) rintf(x);
                
                if (k == 0) k=1;
                if (k == HSP_NB_CASES) k = HSP_NB_CASES - 1;
                
                srt[adr]=HistoInverse[k];
                
            } else srt[adr] = src[adr];
            
            
            
        }
        
    }
    
    
    
    
    
    void fk_fill_histo(float *Histo,float *HistoInverse, int iRangeMax, float *in, int width, int height)
    {
        
        
        // clear data
        for (int i=0; i < iRangeMax; i++) Histo[i]=0.0f;  
        
        
        // fill histogram
        int icounter = 0;
        for (int i=0; i < width * height; i++)
        {
            int ivalue = rintf(in[i]);
            //if (ivalue < 0) ivalue = 0;
            
            if (ivalue >= 0)
            {
                if (ivalue >= iRangeMax) ivalue = iRangeMax - 1;
                Histo[ivalue]++;
                icounter++;
            }
            
            
        }
        
        
        // accumulate histogram
        for (int i=1; i < iRangeMax;i++)   
            Histo[i]+=Histo[i-1]; 
        
        
        //for (int i=0; i < iRangeMax;i++)   
        //    printf("%d: %f\n",i,Histo[i]);
        
        
        // compute the inverse histogram
        HistoInverse[0]=0.0;
        
        
        float xp=0.0;
        
        
        for (int ii=1; ii <= HSP_NB_CASES;ii++)
        {
            
            int x=-1;
            
            float aux = (float) ii * (float) icounter / (float) (HSP_NB_CASES + 1); 
            
            for (int k=0; k <  iRangeMax && Histo[k] < aux; k++) {x=k;} 
            
            if (x==-1)  HistoInverse[ii]=0.0;
            else 
            {
                float dif = Histo[x+1]-Histo[x];
                
                if( fabs(dif) < 0.01)  HistoInverse[ii] =    (float) x;
                else HistoInverse[ii] =    (float) x + ((aux - Histo[x]) / (dif));
                
                HistoInverse[ii] =    MAX( (float) xp,  HistoInverse[ii] );
            }
            
            xp=HistoInverse[ii];
        } 
        
        
        
    }
    
    
    
    
    void fk_histogram_midway(float *in1, float *in2, float *out1, float *out2, int width1, int height1, int width2, int height2)
    {
        
        float fRangeMax;
        fRangeMax =   fpMax(in1, NULL, width1 * height1);
        fRangeMax =   MAX(fRangeMax, fpMax(in2, NULL, width1 * height1));
        
        
        // memory
        int iRangeMax = (int) rintf(fRangeMax) + 1;
        
        float *histo1 = new float[iRangeMax];
        float *histo2 = new float[iRangeMax];
        
        float *histoInverse1 = new float[HSP_NB_CASES + 1];
        float *histoInverse2 = new float[HSP_NB_CASES + 1];
        
        
        
        //! compute histograms
        fk_fill_histo(histo1,histoInverse1, iRangeMax, in1, width1, height1);
        fk_fill_histo(histo2,histoInverse2, iRangeMax, in2, width2, height2);
        
        
        
        
        //! compute average histogram
        float *histoAverage = new float[HSP_NB_CASES+1];
        for (int i=0; i < HSP_NB_CASES + 1; i++) histoAverage[i] = 0.0f;
        
        
        for (int k=0; k < HSP_NB_CASES + 1; k++)  
        {
            histoAverage[k] = 0.5f * (histoInverse1[k] + histoInverse2[k]);
        }
        
        
        
        
        // specificate histograms
        
        for (int k=0; k < HSP_NB_CASES + 1; k++)  
        {
            histoInverse1[k] =  histoInverse2[k] = histoAverage[k];
        }
        
        
        fk_apply_histo(histo1,histoInverse1,iRangeMax,in1, out1, width1, height1);
        fk_apply_histo(histo2,histoInverse2,iRangeMax,in2, out2, width2, height2);
        
        
        
    }
    
    
    
    void fk_histogram_midway_sequence(float **in, float **out, int nImages, int width, int height)
    {
        
        
        float fRangeMax;
        fRangeMax =   fpMax(in[0], NULL, width * height);		
        for (int ii=1; ii < nImages ; ii++)
        {
            fRangeMax =   MAX(fRangeMax, fpMax(in[ii], NULL, width * height));
        }
        
        
        
        //! memory
        int iRangeMax = (int) rintf(fRangeMax) + 1;
        
        
        float **histo = new float*[nImages];
        float **histoInverse = new float*[nImages];
        
        for (int ii=0; ii < nImages; ii++)
        {
            
            histo[ii] = new float[iRangeMax];	
            histoInverse[ii] = new float[HSP_NB_CASES + 1];
            
        }
        
        
        //! compute histograms
        for (int ii=0; ii < nImages; ii++)
        {
            
            fk_fill_histo(histo[ii],histoInverse[ii], iRangeMax, in[ii], width, height);
            
        }
        
        
        
        //! compute average histogram
        float *histoAverage = new float[HSP_NB_CASES+1];
        for (int i=0; i < HSP_NB_CASES + 1; i++) histoAverage[i] = 0.0f;
        
        
        for (int k=0; k < HSP_NB_CASES + 1; k++)  
        {
            for (int ii=0; ii < nImages; ii++)
            {
                
                histoAverage[k]  += histoInverse[ii][k];
            }
            
            
            histoAverage[k] /= (float) nImages;
            
        }
        
        
        
        
        //! specificate histograms
        for (int k=0; k < HSP_NB_CASES + 1; k++)  
        {
            for (int ii=0; ii < nImages; ii++)
            {
                
                histoInverse[ii][k] =  histoAverage[k];
            }
        }
        
        
        for (int ii=0; ii < nImages; ii++)
        {
            
            fk_apply_histo(histo[ii],histoInverse[ii],iRangeMax, in[ii], out[ii], width, height);
            
        }
        
        
        
    }
    
    
    
    void fk_histogram_specification(float *in1, float *in2, float *out, int width1, int height1, int width2, int height2)
    {
        
        assert(in1 !=NULL && in2 != NULL);
        
        float fRangeMax;
        fRangeMax =   fpMax(in1, NULL, width1 * height1);
        fRangeMax =   MAX(fRangeMax, fpMax(in2, NULL, width2 * height2));
        
        
        // memory
        int iRangeMax = (int) rintf(fRangeMax) + 1;
        
        
        
        assert(iRangeMax > 0);
        float *histo1 = new float[iRangeMax];
        float *histo2 = new float[iRangeMax];
        
        float *histoInverse1 = new float[HSP_NB_CASES + 1];
        float *histoInverse2 = new float[HSP_NB_CASES + 1];
        
        
        
        
        
        // compute histograms
        fk_fill_histo(histo1,histoInverse1, iRangeMax, in1, width1, height1);
        fk_fill_histo(histo2,histoInverse2, iRangeMax, in2, width2, height2);
        
        
        
        
        
        // specificate histogram
        fk_apply_histo(histo2,histoInverse1,iRangeMax,in2, out, width2, height2);
        
        delete[] histo1;
        delete[] histo2;
        
        delete[] histoInverse1;
        delete[] histoInverse2;
        
    }
    
    
    
    
    
    
    
    //
    //! Image convolutions
    //
    
    float*  fiFloatGaussKernel(float std, int & size)
    {
        
        
        
        int n = 4 * ceilf(std) + 1; 
        size = n;
        
        
        float* u = new float[n];
        
        
        if (n==1)  u[0]=1.0;
        else
        {
            
            int ishift = (n-1) / 2;
            
            for (int i=ishift; i < n; i++) 
            {
                
                float v = (float)(i - ishift) / std;
                
                u[i] = u[n-1-i] = (float) exp(-0.5*v*v); 
                
            }
            
        }	
        
        
        // normalize
        float fSum = 0.0f;
        for (int i=0; i < n; i++) fSum += u[i];	
        for (int i=0; i < n; i++)  u[i] /= fSum;
        
        
        return u;
        
    }
    
    
    
    float * fiFloatDirectionalGaussKernel(float xsigma, float ysigma, float angle, float *kernel, int kwidth, int kheight)
    {
        
        int ksize = kwidth;
        
        float xsigma2 = xsigma*xsigma;
        float ysigma2 = ysigma*ysigma;
        
        int l2 = ksize/2;
        for(int y = -l2; y <= l2; y++)
            for(int x = -l2; x <= l2; x++)
            {
                
                float a = (float) angle * PI / 180.0f;
                float sina = sin(a);
                float cosa = cos(a);
                
                float ax = (float) x * cosa + (float) y * sina;
                float ay = -(float) x * sina + (float) y * cosa;
                kernel[(y+l2) * ksize + x + l2] =  exp(-(ax*ax)/(2.0f*xsigma2)  - (ay*ay)/(2.0f*ysigma2) );  
                
            }
        
        
        float sum=0.0;
        for(int i=0; i < ksize*ksize; i++) sum += kernel[i];
        for(int i=0; i < ksize*ksize; i++) kernel[i] /= sum;
        
        kwidth = ksize;
        kheight = ksize;
        
        return kernel;
    }
    
    
    
    
    
    
    
    void fiFloatBufferConvolution(float *buffer,float *kernel,int size,int ksize)
    {
        
        for (int i = 0; i < size; i++) {
            
            float sum = 0.0;
            float *bp = &buffer[i];
            float *kp = &kernel[0];
            
            
            
            int k=0;
            for(;k + 4 < ksize;  bp += 5, kp += 5, k += 5) 
                sum += bp[0] * kp[0] +  bp[1] * kp[1] + bp[2] * kp[2] +
                bp[3] * kp[3] +  bp[4] * kp[4];
            
            
            for(; k < ksize; bp++ , kp++, k++)  sum += *bp * (*kp);
            
            buffer[i] = sum;
        }
    }
    
    
    
    
    void fiFloatHorizontalConvolution(float *u, float *v, int width, int height, float *kernel, int ksize, int boundary)
    {
        
        int halfsize = ksize / 2;
        int buffersize = width + ksize;
        float *buffer = new float[buffersize];
        
        for (int r = 0; r < height; r++) {
            
            /// symmetry
            int l = r*width;
            if (boundary == 1)
                for (int i = 0; i < halfsize; i++)
                    buffer[i] = u[l + halfsize - 1 - i ];
            else
                for (int i = 0; i < halfsize; i++)
                    buffer[i] = 0.0;
            
            
            for (int i = 0; i < width; i++)
                buffer[halfsize + i] = u[l + i];
            
            
            if (boundary == 1)
                for (int i = 0; i <  halfsize; i++)
                    buffer[i + width + halfsize] = u[l + width - 1 - i];
            else 
                for (int i = 0; i <  halfsize; i++)
                    buffer[i + width + halfsize] = 0.0;
            
            fiFloatBufferConvolution(buffer, kernel, width, ksize);
            for (int c = 0; c < width; c++)
                v[r*width+c] = buffer[c];
        }
        
        
        delete[] buffer;
        
    }
    
    
    
    void fiFloatVerticalConvolution(float *u, float *v, int width, int height, float *kernel,int ksize, int boundary)
    {
        int halfsize = ksize / 2;
        int buffersize = height + ksize;
        float *buffer = new float[buffersize];
        
        for (int c = 0; c < width; c++) {
            
            if (boundary == 1)
                for (int i = 0; i < halfsize; i++)
                    buffer[i] = u[(halfsize-i-1)*width + c];
            else 
                for (int i = 0; i < halfsize; i++)
                    buffer[i] = 0.0f;
            
            for (int i = 0; i < height; i++)
                buffer[halfsize + i] = u[i*width + c];
            
            if (boundary == 1)
                for (int i = 0; i < halfsize; i++)
                    buffer[halfsize + height + i] = u[(height - i - 1)*width+c];
            else
                for (int i = 0; i < halfsize; i++)
                    buffer[halfsize + height + i] = 0.0f;
            
            fiFloatBufferConvolution(buffer, kernel, height, ksize);
            
            for (int r = 0; r < height; r++)
                v[r*width+c] = buffer[r];
            
        }
        
        delete[] buffer;
    }
    
    
    
    void fiConvol(float *u,float *v,int width,int height,float *kernel,int kwidth,int kheight)
    {
        
        int K2 = kwidth / 2;
        int L2 = kheight / 2;
        
        for(int y=0 ; y < height; y++) 
            for (int x=0 ; x < width; x++) 
            {
                
                float S = 0.0;
                
                for (int l = -L2; l <= L2; l++) 
                    for (int k = -K2 ; k<= K2; k++)
                    { 
                        int px=x+k;
                        int py=y+l;
                        
                        if (px>=0 && px < width && py>=0 && py<height)
                            S += u[width*py + px] * kernel[kwidth*(l+L2) + k+K2];
                    }
                
                v[y*width+x] = (float) S;
                
            }
    }
    
    
    
    void fiSepConvol(float *u,float *v,int width,int height,float *xkernel, int xksize, float *ykernel, int yksize)
    {
        
        
        int boundary = 1;
        
        if (u != v) memcpy(v, u, width*height*sizeof(float));
        
        fiFloatHorizontalConvolution(v, v, width, height, xkernel, xksize, boundary);
        fiFloatVerticalConvolution(v, v, width, height,  ykernel,  yksize, boundary);
        
    }
    
    
    
    
    void fiGaussianConvol(float *u, float *v, int width, int height, float sigma)
    {
        
        int ksize;	
        float *kernel;
        kernel = fiFloatGaussKernel(sigma,ksize);
        
        
        int boundary = 1;
        
        if (u != v) memcpy(v, u, width*height*sizeof(float));
        
        fiFloatHorizontalConvolution(v, v, width, height, kernel, ksize, boundary);
        fiFloatVerticalConvolution(v, v, width, height,  kernel,  ksize, boundary);
        
        delete[] kernel;
        
    }
    
    
    
    
    
    
    
    
    void fiImageSample(float *igray,float *ogray, int factor, int width, int height)
    {
        
        int swidth = (int) floor( (float) width / (float) factor);
        int sheight = (int) floor( (float) height / (float) factor);
        
        for(int j=0; j < sheight; j++)
            for(int i=0; i < swidth; i++)
                ogray[j*swidth + i] = igray[ j * factor * width +  i*factor ];
        
    }
 
    

    
    
    void fiImageSampleCenter(float *igray,float *ogray, int factor, int width, int height)
    {
        
        int swidth = (int) floor( (float) width / (float) factor);
        int sheight = (int) floor( (float) height / (float) factor);
        int halfSample = (factor - 1) / 2; 
        
        
        for(int j=0; j < sheight; j++)
            for(int i=0; i < swidth; i++)
                if (factor % 2 == 0 && (j * factor + halfSample + 2) < height &&  (i * factor + halfSample + 2) < width)  
                {
                    float aux=0.0;
                    for (int r=0; r < 2; r++)
                        for (int s=0; s < 2; s++)
                        {
                            aux += igray[ (j * factor + halfSample +r ) * width +  i * factor + halfSample + s];
                        
                        }
                    
                    aux /= 4.0;
                    ogray[j*swidth + i] = aux;
                    
                }
                else
                    ogray[j*swidth + i] = igray[  (j * factor + halfSample)* width +  i * factor + halfSample];
        
    }
    
    
    
    
    
    
    void fiImageSampleAglomeration(float *igray,float *ogray, int factor, int width, int height)
    {
        
        int swidth = (int) floor( (float) width / (float) factor);
        int sheight = (int) floor( (float) height / (float) factor);
        
        float fFactor2 = (float) factor * (float) factor;
        
        for(int j=0; j < sheight; j++)
            for(int i=0; i < swidth; i++)
            {
                
                int pi = i * factor;
                int pj = j * factor;
                
                float fSum = 0.0f;
                for (int r=0; r < factor; r++)
                    for (int s=0; s < factor; s++)
                    {
                        fSum += igray[ (pj+s) * width + (pi+r)];
                    }
                
                
                ogray[j*swidth + i] = fSum / fFactor2;
                
            }
        
        
        
    }
    
    
    
    
    //
    //! Patch Processing
    //
    
    
    void fiPatchStatistics(float *fpIn, float *fpMinV, float *fpMaxV, float *fpMeanV, float *fpVarV, float *fpMedianV, float fRadius, int iWidth, int iHeight)
    {
        
        //! Parameters
        int iRadius = (int)(fRadius+1.0);
        int iNeigSize = (2*iRadius+1)*(2*iRadius+1);
        float fRadiusSqr = fRadius * fRadius;
        
        
        float * vector = new float[iNeigSize];
        
        
        
        //! For each pixel
        for(int x=0;x < iWidth;x++)
            for(int y=0;y< iHeight;y++)
            {
                
                int iCount=0;
                float fMin = fLarge;
                float fMax = -fLarge;
                
                
                for(int i=-iRadius;i<=iRadius;i++)
                    for(int j=-iRadius;j<=iRadius;j++)
                        if ((float) (i*i + j*j) <= fRadiusSqr)
                        {
                            
                            int x0=x+i;
                            int y0=y+j;
                            
                            if (x0 >= 0 && y0 >= 0 && x0 < iWidth && y0 < iHeight) 
                            { 
                                float fValue= fpIn[y0*iWidth+x0];
                                
                                if (fValue > fMax) fMax = fValue;
                                if (fValue < fMin) fMin = fValue;
                                
                                vector[iCount] = fValue;
                                iCount++; 
                                
                            }
                            
                        } 
                
                int l = y*iWidth+x;
                if (fpMinV) fpMinV[l] = fMin;		  
                if (fpMaxV) fpMaxV[l] = fMax;	
                
                if (fpMeanV) fpMeanV[l] = fpMean(vector, iCount);
                if (fpVarV)  fpVarV[l] = fpVar(vector, iCount);
                
                
                if (fpMedianV)
                {
                    fpQuickSort(vector, iCount, 0);
                    fpMedianV[l] = vector[iCount / 2];
                }
                
                
                
            }
        
        
    }
    
    
    
    
    void fiComputeIntegralImage(float *in, float *out, int width, int height)
    {
        
        
        
        float * tmp_im = new float[width * height];
        
        // recursivity in the rows
        
        for(int jj=0; jj < height; jj++)
            for(int ii=0; ii < width; ii++)
            {
                if  (ii==0) tmp_im[width*jj+ii] = in[width*jj+ii];
                else tmp_im[width*jj+ii] = tmp_im[width*jj+ii-1] + in[width*jj+ii];
            }
        
        // recursivity in the columns
        for(int ii=0; ii < width; ii++)
            for(int jj=0; jj < height; jj++)
            {
                if(jj==0) out[width*jj+ii] = tmp_im[width*jj+ii];
                else out[width*jj+ii] = out[width*(jj-1)+ii] + tmp_im[width*jj+ii];
            }
        
        
        delete[] tmp_im;
        
    }
    
    
    
    
    
    void fiPatchMin(float *fpIn, float *fpMinV, float fRadius, int iWidth, int iHeight)
    {
        fiPatchStatistics( fpIn, fpMinV, NULL,  NULL,  NULL, NULL,  fRadius,  iWidth, iHeight);
        
    }
    
    
    void fiPatchMax(float *fpIn, float *fpMaxV, float fRadius, int iWidth, int iHeight)
    {
        fiPatchStatistics( fpIn, NULL, fpMaxV,  NULL,  NULL, NULL,  fRadius,  iWidth, iHeight);
        
    }
    
    
    void fiPatchMean(float *fpIn, float *fpMeanV, float fRadius, int iWidth, int iHeight)
    {
        fiPatchStatistics( fpIn, NULL, NULL,  fpMeanV,  NULL, NULL,  fRadius,  iWidth, iHeight);
        
    }
    
    
    void fiPatchVar(float *fpIn, float *fpVarV, float fRadius, int iWidth, int iHeight)
    {
        fiPatchStatistics( fpIn, NULL, NULL,  NULL,  fpVarV, NULL,  fRadius,  iWidth, iHeight);
        
    }    
    
    
    void fiPatchMedian(float *fpIn, float *fpMedianV, float fRadius, int iWidth, int iHeight)
    {
        fiPatchStatistics( fpIn, NULL, NULL,  NULL,  NULL, fpMedianV,  fRadius,  iWidth, iHeight);
    }    
    
    
    
    
    
    
	/**
	 * \brief   Tabulates exp(-x) function
	 *
	 *
	 * @param[in]  lut	vector
	 * @param[in]  size	length of the vector
	 *
	 */
	
	void  wxFillExpLut(float *lut, int size)
	{
		for(int i=0; i< size;i++)   lut[i]=   expf( - (float) i / LUTPRECISION);
	}
	
	
	
	
	/**
	 * \brief   Computes exp(-x) using lut table
	 *
	 *
	 * @param[in]  dif	value
	 * @param[in]  lut	lookup table
	 *
	 */
	float wxSLUT(float dif, float *lut)
	{
		
		if (dif >= (float) LUTMAXM1) return 0.0;
		
		int  x= (int) floor( (double) dif * (float) LUTPRECISION);
		
		float y1=lut[x];
		float y2=lut[x+1];
		
		return y1 + (y2-y1)*(dif*LUTPRECISION -  x); 
	}
	

    
    
    
    
    
    
    
    
    
    //
    //! Image conversion
    //
    
    void fiRgb2Yuv(float *r,float *g,float *b,float *y,float *u,float *v,int width,int height) 
    {
        int size=height*width;
        
        for(int i=0;i<size;i++){
            y[i] = ( COEFF_YR *  r[i] + COEFF_YG * g[i] + COEFF_YB * b[i]);
            u[i] =  ( r[i] - y[i]);
            v[i] =  ( b[i] - y[i]);
        }
        
    }
    
    
    
    
    
    void fiYuv2Rgb(float *r,float *g,float *b,float *y,float *u,float *v, int width,int height)  
    {
        
        
        int iwh=height*width;
        
        for(int i=0;i<iwh;i++){
            
            g[i] =  ( y[i] - COEFF_YR * (u[i] + y[i]) - COEFF_YB * (v[i] +  y[i]) ) / COEFF_YG;
            r[i] =  ( u[i] + y[i]);
            b[i] =  ( v[i] + y[i]);
            
        }
        
    }
    
    
    void fiRgb2YuvO(float *r,float *g,float *b,float *y,float *u,float *v,int width,int height) 
    {
        int size=height*width;
        
        for(int i=0;i<size;i++)
        {
            y[i] =    0.577350 * r[i] + 0.577350 * g[i] + 0.577350  * b[i];
            u[i] =    0.707106 * r[i]					- 0.707106  * b[i];
            v[i] =    0.408248 * r[i] - 0.816496 * g[i]	+ 0.408248  * b[i];
        }
        
    }
    
    
    
    
    void fiYuvO2Rgb(float *r,float *g,float *b,float *y,float *u,float *v, int width,int height)  
    {
        
        
        int iwh=height*width;
        
        for(int i=0;i<iwh;i++)
        {
            
            r[i] =  0.577350 * y[i]	+ 0.707106 * u[i]	+  0.408248 * v[i];
            g[i] =  0.577350 * y[i]						-  0.816496 * v[i];
            b[i] =  0.577350 * y[i] - 0.707106 * u[i]	+  0.408248 * v[i];
            
        }
        
    }
    
    
    
    
    
    //
    //! Patch distances
    //
    
    
    
    float fiL2FloatDist(float *u0,float *u1,int i0,int j0,int i1,int j1,int xradius, int yradius, int width0, int width1)
    {
        
        
        
        float dist=0.0;       
        for (int s=-yradius; s<= yradius; s++){
            
            int l = (j0+s)*width0 + (i0-xradius);
            float *ptr0 = &u0[l];
            
            l = (j1+s)*width1 + (i1-xradius);
            float *ptr1 = &u1[l];
            
            for(int r=-xradius;r<=xradius;r++,ptr0++,ptr1++){	float dif = (*ptr0 - *ptr1); dist += (dif*dif); }
            
        }
        
        return dist;
    }
    
    
    
    
    float fiL2FloatDist(float *u0,float *u1,int i0,int j0,int i1,int j1,int xradius, int width0, int width1)
    {
        
        
        return fiL2FloatDist(u0,u1,i0,j0,i1,j1, xradius, xradius,  width0,  width1);
    }
    
    
    
    
    float fiL2FloatDist(float **u0,float **u1,int i0,int j0,int i1,int j1,int xradius, int yradius, int channels, int width0, int width1)
    {
        
        float dif = 0.0f;
        
        
        for (int ii=0; ii < channels; ii++) {
            
            dif += fiL2FloatDist(u0[ii],u1[ii],i0,j0,i1,j1,xradius,yradius,width0,width1);
            
        }
        
        dif /= (float) channels;
        return dif;
        
    }
    
    
    
    
    float fiL2FloatWDist ( float * u0, float *u1, int i0, int j0,int i1,int j1,int xradius, int yradius, float *kernel, int width0, int width1)
    {
        
        float *ptrk=&kernel[0];
        float dist=0.0;       
        for (int s=-yradius; s<= yradius; s++)
        {
            
            int l = (j0+s)*width0 + (i0-xradius);
            float *ptr0 = &u0[l];
            
            l = (j1+s)*width1 + (i1-xradius);
            float *ptr1 = &u1[l];
            
            for(int r=-xradius;r<=xradius;r++,ptr0++,ptr1++,ptrk++){ float dif = (*ptr0 - *ptr1); dist += *ptrk*(dif*dif); }
            
            
        }
        
        return dist;
    }
    
    
    
    
    
    float fiL2FloatWDist ( float ** u0, float **u1,int i0,int j0,int i1,int j1,int xradius, int yradius, float *kernel, int channels, int width0, int width1)
    {
        
        float dif = 0.0f;
        for (int ii=0; ii < channels; ii++) {
            
            dif += fiL2FloatWDist(u0[ii],u1[ii],i0,j0,i1,j1,xradius, yradius, kernel,width0,width1);
            
        }
        
        dif /= (float) channels;
        
        return dif;
    }
    
    
    
    
    
    
    //
    //! FFT STUFF
    //
#ifdef USE_FFTW
    
#define c_re(c) ((c)[0])
#define c_im(c) ((c)[1])
	
	
	void FFT_init (void) {
		FILE *wisdom_file;
		int res;
		
		wisdom_file = fopen("/tmp/.fftwisdom","r");
		if (wisdom_file == NULL) 
			fprintf(stderr, "Generating wisdom file...\n");
		else {
			res = fftwf_import_wisdom_from_file(wisdom_file);
			fclose (wisdom_file);
		}
		
	}
	
	void FFT_end (void) {
		FILE *wisdom_file;
		
		// export wisdom for later use
		wisdom_file = fopen("/tmp/.fftwisdom","w");
		fftwf_export_wisdom_to_file(wisdom_file);
		fclose (wisdom_file);
	}
	
    
	// compute the inverse FFT and normalizes.
	void FFT_2d_inv (fftwf_complex* in, fftwf_complex* out, int nx, int ny)
	{
		fftwf_plan p;
		int i;
		int numpix = nx*ny;
		
      #pragma omp critical (make_plan) 
		p = fftwf_plan_dft_2d(ny,nx, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
		fftwf_execute(p);
		fftwf_destroy_plan(p);
		// normalize
		for(i=0; i<numpix;i++) { c_re(out[i]) /= numpix; c_im(out[i]) /= numpix; }
		
	}
	
	// compute the forward FFT.
	void FFT_2d_fwd (fftwf_complex* in,fftwf_complex* out, int nx, int ny)
	{
		fftwf_plan p;
		
      #pragma omp critical (make_plan) 
		p = fftwf_plan_dft_2d(ny,nx, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftwf_execute(p);
		fftwf_destroy_plan(p);
	}
	
	// compute the inverse FFT and normalizes.
	void FFT_1d_inv (fftwf_complex* in, fftwf_complex* out, int nx)
	{
		fftwf_plan p;
		int i;
		int numpix = nx;
		
      #pragma omp critical (make_plan) 
		p = fftwf_plan_dft_1d(nx, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
		fftwf_execute(p);
		fftwf_destroy_plan(p);
		// normalize
		for(i=0; i<numpix;i++) { c_re(out[i]) /= numpix; c_im(out[i]) /= numpix; }
		
	}
	
	// compute the forward FFT.
	void FFT_1d_fwd (fftwf_complex* in,fftwf_complex* out, int nx)
	{
		fftwf_plan p;
		
      #pragma omp critical (make_plan) 
		p = fftwf_plan_dft_1d(nx, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftwf_execute(p);
		fftwf_destroy_plan(p);
	}
	
    
	// 1d fftw
	void fft1d(float *in_re,float *in_im,float *out_re,float *out_im, int i_flag, int nx) {
		int i;
		fftwf_complex *in,*out;
		
		if ((!out_re) && (!out_im)) 
		{
			fprintf(stderr,"At least one output needed\n");
			return;
		}
		
		
		/***** allocate images *****/
		in  = new fftwf_complex[nx];
		out = new fftwf_complex[nx];
		
		/***** copy from input *****/
		if (!in_re) {
			fprintf(stderr,"No real part!?\n");
			return;
		}
		
		/* real part */
		for(i=0; i<nx;i++) {
			c_re(in[i]) = in_re[i];
			c_im(in[i]) = 0;
		}
		/* imaginary part (if present) */
		if ( in_im ) 
			for(i=0; i<nx;i++) {
				c_im(in[i]) = in_im[i];
			}
		
		/***** transform *****/
		//   FFT_init();
		
		if(!i_flag)
			FFT_1d_fwd (in, out, nx); /* F in */
		else 
			FFT_1d_inv (in, out, nx); /* F^-1 in */
		
		//   FFT_end();
		
		/***** allocations *****/
		
		if(out_re)
			for(i=0; i<nx;i++) 
				out_re[i] = c_re(out[i]);
		if(out_im)
			for(i=0; i<nx;i++) 
				out_im[i] = c_im(out[i]);
		
		
		/***** free temporary *****/
		delete[] in;
		delete[] out;
		
	}
	
	
	
	// NEW FFT2d USES FFTW
	// 2d fftw
	//in_re, out_re and out_img has memory
	void fft2d(float *in_re, float *in_im,float *out_re,float *out_im,int i_flag, int width, int height)
	{
		int nx,ny,i;
		fftwf_complex *in,*out;
		
		nx = width;
		ny = height;
		
		if ((!out_re) && (!out_im)) 
		{
			fprintf(stderr,"At least one output needed\n");
			return;
		}
		
		
		/***** allocate images *****/
		in  = new fftwf_complex[nx*ny];
		out = new fftwf_complex[nx*ny];
		
		/***** copy from input *****/
		if (!in_re) {
			fprintf(stderr,"No real part!?\n");
			return;
		}
		
		/* real part */
		for(i=0; i<nx*ny;i++) {
			c_re(in[i]) = in_re[i];
			c_im(in[i]) = 0;
		}
		/* imaginary part (if present) */
		if ( in_im ) 
			for(i=0; i<nx*ny;i++) {
				c_im(in[i]) = in_im[i];
			}
		
		/***** transform *****/
		//   FFT_init();
		
		if(!i_flag)
			FFT_2d_fwd (in, out, nx, ny); /* F in */
		else 
			FFT_2d_inv (in, out, nx, ny); /* F^-1 in */
		
		//   FFT_end();
		
		/***** allocations *****/
		
		if(out_re)
			for(i=0; i<nx*ny;i++) 
				out_re[i] = c_re(out[i]);
		if(out_im)
			for(i=0; i<nx*ny;i++) 
				out_im[i] = c_im(out[i]);
		
		
		/***** free temporary *****/
		delete[] in;
		delete[] out;
		
	}
    
    
    
    
	void fiFFTZoom(float *in, float *out, float zoomx, float zoomy, int width, int height)
	{
		
		
		
		assert((zoomx >= 1. && zoomy >= 1.) || (zoomx < 1. && zoomy < 1.) );
		
		int swidth = width / 2 + 1;
		
		
		// Perform FFT
		// FFTW are the left half image of the mw coefficients
		// MW coefficients are placed 
		//								x y
		//								y x
		
		
		fftwf_complex  *fft_out;
		fftwf_plan p;
		
		fft_out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * (width/2+1) * height);
		
      #pragma omp critical (make_plan) 
		p =  fftwf_plan_dft_r2c_2d(height, width, in, fft_out, FFTW_ESTIMATE);
		
		fftwf_execute(p);
		
		
		
		// Padding of FFT output
		int widthz = (int)rintf( (float) width *  (float) zoomx);
		int heightz = (int)rintf( (float) height *  (float) zoomy);
		
		int swidthz = (widthz / 2 + 1);
		
		
		fftwf_complex *fft_pout;
		fft_pout = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * (widthz/2+1) * heightz);
		
		
		
		for (int i=0; i < swidthz * heightz ; i++)
		{
			fft_pout[i][0] = 0.0f;
			fft_pout[i][1] = 0.0f;
		}
		
		
		
		if (zoomx<1. && zoomy<1.) {
			
			for (int x= 0 ; x < swidthz ;x++)
				for (int y=-heightz/2 ; y<= heightz/2;y++) 
				{
					
					int adr  = x + ((y+ height )%height ) * swidth ;
					int adrz = x + ((y+ heightz)%heightz) * swidthz;
					
					float factor = zoomx * zoomy;
					
					fft_pout[adrz][0] = fft_out[adr][0] * factor;
					fft_pout[adrz][1] = fft_out[adr][1] * factor;
					
					
				}
			
			
		} else if (zoomx >= 1.0 && zoomy >= 1.0){
			
			
			int adr, adrz;
			for (int x=0; x < swidth ; x++)
				for (int y = -height/2; y <= height/2; y++) 
				{
					adr = x + ( (y + height)% height ) * swidth ; 
					adrz = x + ( (y + heightz) % heightz ) * swidthz ; 
					
					float factor = zoomx * zoomy 
					* ((2*y==height || 2*y==-height)?0.5:1.0);
					
					fft_pout[adrz][0] = fft_out[adr][0] * factor;
					fft_pout[adrz][1] = fft_out[adr][1] * factor;
					
				}
			
			
			
		} 
		
		
		
		
		//   Backward fft	
      #pragma omp critical (make_plan) 
		p =  fftwf_plan_dft_c2r_2d(heightz, widthz, fft_pout, out, FFTW_ESTIMATE);
		
		fftwf_execute(p);
		
		
		for (int i=0; i < widthz * heightz; i++) out[i] /= (float) (widthz * heightz);
		
		
		
		fftwf_destroy_plan(p);
		fftwf_free(fft_out);	 
		fftwf_free(fft_pout);
		
	}
	
    
    // NEW SHEAR USES FFTW
	void fiFFTShearPascal(double dAmountOfShear, float *pFloatImageInput,float *pFloatImageOutput, int iAxis, float fDelta, int PleaseDoNotCenter, int width, int height)
	{
		int i, j, iSize, iHalfSize, iOtherSize, iHalfOtherSize;
		fftwf_complex *pmxsignal, *pmxsignalTransformed;
		double dTranslation;
		float fCos, fSin, fReal;
		int odd,oddOther;
		
		
		// iSize is the size of the image in the direction of the axis of shear. 
		// iOtherSize the size of the image in the orthogonal direction 
		if(iAxis == 0) // HORIZONTAL
		{
			iHalfSize = (iSize = width) >> 1;
			iHalfOtherSize = (iOtherSize = height) >> 1;
		}
		else
		{
			iHalfSize = (iSize = width) >> 1;
			iHalfOtherSize = (iOtherSize = height) >> 1;
		}
		
		if ((iOtherSize & 1)!= 0) oddOther=1; else oddOther=0;
		if ((iSize & 1) !=0) odd=1; else odd=0;
		
		// Create temporary mxsignals to compute the fourier transform of a line 
		// (or a column) of the image 
		pmxsignal = new fftwf_complex[iSize]; // mw_change_fmxsignal(0, iSize);
		pmxsignalTransformed  = new fftwf_complex[iSize]; //= mw_new_fmxsignal();
		
		// THERE IS A SMALL INPROVEMENT IN RECYCLING THE PLANS
		// Naming of forward and backward fft's are reversed wrt fftw
      #pragma omp critical (make_plan) 
		fftwf_plan pfwd = fftwf_plan_dft_1d(iSize, pmxsignal, pmxsignalTransformed, FFTW_BACKWARD, FFTW_ESTIMATE);
      #pragma omp critical (make_plan) 
		fftwf_plan pbwd = fftwf_plan_dft_1d(iSize, pmxsignalTransformed, pmxsignal, FFTW_FORWARD, FFTW_ESTIMATE);
		
		
		for(i = 0; i < iOtherSize; i++)
		{
			memset(pmxsignal, 0, iSize * sizeof(fftwf_complex));
			
			if(iAxis == 0) {// HORIZONTAL
				for(j = 0; j < iSize; j++)
					c_re(pmxsignal[j]) = pFloatImageInput[i * iSize + j];
			}
			else {
				for(j = 0; j < iSize; j++)
					c_re(pmxsignal[j]) = pFloatImageInput[j * iOtherSize + i];
			}
			
			// Compute the FFT of the current line (or column) of the image 
			//      FFT_1d_fwd (pmxsignal,pmxsignalTransformed, iSize);
			fftwf_execute(pfwd);          // RECYCLING PLANS
			
			
			if (PleaseDoNotCenter) 
				dTranslation = - (i * dAmountOfShear + fDelta) * 2. * M_PI;
			else 
				if (oddOther)
					dTranslation = - ((i - iHalfOtherSize ) * dAmountOfShear + fDelta) * 2. * M_PI;
				else
					dTranslation = - ((i - iHalfOtherSize + .5) * dAmountOfShear + fDelta) * 2. * M_PI;
			
			for(j = 1; j < iHalfSize+odd; j++)
			{
				fCos = (float) cos(j * dTranslation / iSize);
				fSin = (float) sin(j * dTranslation / iSize);
				fReal = c_re(pmxsignalTransformed[j]); 
				c_re(pmxsignalTransformed[j]) = fCos * fReal - fSin * c_im(pmxsignalTransformed[j]);
				c_im(pmxsignalTransformed[j]) = fSin * fReal + fCos * c_im(pmxsignalTransformed[j]);
				c_re(pmxsignalTransformed[iSize - j]) = c_re(pmxsignalTransformed[j]);
				c_im(pmxsignalTransformed[iSize - j]) = -c_im(pmxsignalTransformed[j]);
			}
			
			if (odd == 0) {
				c_re(pmxsignalTransformed[iHalfSize]) = cos(dTranslation * .5) * c_re(pmxsignalTransformed[iHalfSize]) - sin(dTranslation * .5) * c_im(pmxsignalTransformed[iHalfSize]);
				c_im(pmxsignalTransformed[iHalfSize]) = 0.;
			}
			
			// Compute the inverse FFT of the current line (or column) 
			//FFT_1d_inv (pmxsignalTransformed, pmxsignal,iSize);
			fftwf_execute(pbwd); // RECYCLING PLANS
			for(j=0; j<iSize;j++) { c_re(pmxsignal[j]) /= iSize; c_im(pmxsignal[j]) /= iSize; }
			
			if(iAxis == 0) {// HORIZONTAL
				for(j = 0; j < iSize; j++)
					pFloatImageOutput[i * iSize + j] = c_re(pmxsignal[j]);
			}
			else {
				for(j = 0; j < iSize; j++)
					pFloatImageOutput[j * iOtherSize + i] = c_re(pmxsignal[j]);
			}
		}
		
		// Delete the previously allocated temporary mxsignals 
		fftwf_destroy_plan(pfwd); // RECYCLING PLANS
		fftwf_destroy_plan(pbwd);
		delete[] pmxsignalTransformed;
		delete[] pmxsignal;
	}
	
    
    
    
#endif	
	
    
    
    
    
    ////////////////////////////////////////////////////////////////////////////////
	//! Begin Spline Interpolation MEGAWAVE
	////////////////////////////////////////////////////////////////////////////////
	
	
	
	//! extract image value (even outside image domain) 
	float vMW(float *in,int x,int y,float bg, int width, int height)
	{
		if (x<0 || x>=width || y<0 || y>=height)
			return bg; 
        else return(in[y*width+x]);
	}
	
	
	
	
	
	//! c[] = values of interpolation function at ...,t-2,t-1,t,t+1,... 
	//! coefficients for cubic interpolant (Keys' function)
	void keysMW(float *c,float t,float a)
	{
		float t2,at;
		
		t2 = t*t;
		at = a*t;
		c[0] = a*t2*(1.0-t);
		c[1] = (2.0*a+3.0 - (a+2.0)*t)*t2 - at;
		c[2] = ((a+2.0)*t - a-3.0)*t2 + 1.0;
		c[3] = a*(t-2.0)*t2 + at;
	}
	
    
    
	//! Coefficients for cubic spline
	void spline3MW(float *c,float t)
	{
		float tmp;
		
		tmp = 1.-t;
		c[0] = 0.1666666666*t*t*t;
		c[1] = 0.6666666666-0.5*tmp*tmp*(1.+t);
		c[2] = 0.6666666666-0.5*t*t*(2.-t);
		c[3] = 0.1666666666*tmp*tmp*tmp;
	}
    
	
	
	
	
	
	//! Init spline n
	void init_splinenMW(float *a,int n)
	{
		int k;
		
		a[0] = 1.;
		for (k=2;k<=n;k++) a[0]/=(float)k;
		for (k=1;k<=n+1;k++)
			a[k] = - a[k-1] *(float)(n+2-k)/(float)k;
	}
	
    
	//! Fast integral power function
	float ipowMW(float x,int n)
	{
		float res;
		
		for (res=1.;n;n>>=1) {
			if (n&1) res*=x;
			x*=x;
		}
		return(res);
	}
	
    
	//! coefficients for spline of order >3 
	void splinenMW(float *c,float t,float *a,int n)
	{
		int i,k;
		float xn;
		
		memset((void *)c,0,(n+1)*sizeof(float));
		for (k=0;k<=n+1;k++) { 
			xn = ipowMW(t+(float)k,n);
			for (i=k;i<=n;i++) 
				c[i] += a[i-k]*xn;
		}
	}
    
	
	
	
	//! called by invspline1d
	//! takes c as input and output some kind of sum combining c and z
	double initcausalMW(double *c,int n,double z)
	{
		double zk,z2k,iz,sum;
		int k;
		
		zk = z; iz = 1./z;
		z2k = pow(z,(double)n-1.);
		sum = c[0] + z2k * c[n-1];
		z2k = z2k*z2k*iz;
		for (k=1;k<=n-2;k++) {
			sum += (zk+z2k)*c[k];
			zk *= z;
			z2k *= iz;
		}
		return (sum/(1.-zk*zk));
	}
	
	
	
	//! called by invspline1d
	//! takes c as input and output some kind of sum combining c and z
	double initanticausalMW(double *c,int n,double z)
	{
		return((z/(z*z-1.))*(z*c[n-2]+c[n-1]));
	}
	
	
	
	//! called by finvspline
	void invspline1DMW(double *c,int size,double *z,int npoles)
	{
		double lambda;
		int n,k;
		
		/* normalization */
		for (k=npoles,lambda=1.;k--;) lambda *= (1.-z[k])*(1.-1./z[k]);
		for (n=size;n--;) c[n] *= lambda;
		
		/*----- Loop on poles -----*/
		for (k=0;k<npoles;k++) {
			
			/* forward recursion */
			c[0] = initcausalMW(c,size,z[k]);
			for (n=1;n<size;n++) 
				c[n] += z[k]*c[n-1];
			
			/* backwards recursion */
			c[size-1] = initanticausalMW(c,size,z[k]);
			for (n=size-1;n--;) 
				c[n] = z[k]*(c[n+1]-c[n]);
			
		}
	}
	
	
	
	//! Main function filling coefficients
	void finvsplineMW(float *in,int order,float *out, int width, int height)
	{
		double *c,*d,z[5];
		int npoles,nx,ny,x,y;
		
		ny = height; nx = width;
		
		/* initialize poles of associated z-filter */
		switch (order) 
		{
			case 2: z[0]=-0.17157288;  /* sqrt(8)-3 */
				break;
				
			case 3: z[0]=-0.26794919;  /* sqrt(3)-2 */ 
				break;
				
			case 4: z[0]=-0.361341; z[1]=-0.0137254;
				break;
				
			case 5: z[0]=-0.430575; z[1]=-0.0430963;
				break;
				
			case 6: z[0]=-0.488295; z[1]=-0.0816793; z[2]=-0.00141415;
				break;
				
			case 7: z[0]=-0.53528; z[1]=-0.122555; z[2]=-0.00914869;
				break;
				
			case 8: z[0]=-0.574687; z[1]=-0.163035; z[2]=-0.0236323; z[3]=-0.000153821;
				break;
				
			case 9: z[0]=-0.607997; z[1]=-0.201751; z[2]=-0.0432226; z[3]=-0.00212131;
				break;
				
			case 10: z[0]=-0.636551; z[1]=-0.238183; z[2]=-0.065727; z[3]=-0.00752819;
				z[4]=-0.0000169828;
				break;
				
			case 11: z[0]=-0.661266; z[1]=-0.27218; z[2]=-0.0897596; z[3]=-0.0166696; 
				z[4]=-0.000510558;
				break;
				
			default:
				printf("finvspline: order should be in 2..11.\n");
				exit(-1);
		}
		
		npoles = order/2;
		
		/* initialize double array containing image */
		c = (double *)malloc(nx*ny*sizeof(double));
		d = (double *)malloc(nx*ny*sizeof(double));
		for (x=nx*ny;x--;) 
			c[x] = (double)in[x];
		
		/* apply filter on lines */
		for (y=0;y<ny;y++) 
			invspline1DMW(c+y*nx,nx,z,npoles);
		
		/* transpose */
		for (x=0;x<nx;x++)
			for (y=0;y<ny;y++) 
				d[x*ny+y] = c[y*nx+x];
		
		/* apply filter on columns */
		for (x=0;x<nx;x++) 
			invspline1DMW(d+x*ny,ny,z,npoles);
		
		/* transpose directy into image */
		for (x=0;x<nx;x++)
			for (y=0;y<ny;y++) 
				out[y*nx+x] = (float)(d[x*ny+y]);
		
		/* free array */
		free(d);
		free(c);
	}
	
	
    
    float evaluate_splineMW(float *input, float *ref, float xp, float yp, float *ak, int order, float bg, int width, int height)
    {
        
        float  cx[12],cy[12];
        float p=-0.5;
        
        int xi = (int)floor((double)xp); 
        int yi = (int)floor((double)yp);
        
        float res;
        int n1, n2;
        
        if (order == 0) 
        { 
            if (xi<0 || xi>=width || yi<0 || yi>=height)        res = bg; 
            else res = input[yi*width+xi];
            
        } else 
        { 
            
            
            if (xp<0. || xp>=(float)width || yp<0. || yp>=(float)height) res=bg; 
            else {
                
                //xp -= 0.5; yp -= 0.5;
                
                float ux = xp-(float)xi;
                float uy = yp-(float)yi;
                
                switch (order) 
                {
                    case 1:  
                        n2 = 1;
                        cx[0]=ux; cx[1]=1.-ux;
                        cy[0]=uy; cy[1]=1.-uy;
                        break;
                        
                    case -3:  
                        n2 = 2;
                        keysMW(cx,ux,p);
                        keysMW(cy,uy,p);
                        break;
                        
                    case 3:  
                        n2 = 2;
                        spline3MW(cx,ux);
                        spline3MW(cy,uy);
                        break;
                        
                    default:  
                        n2 = (1+order)/2;
                        splinenMW(cx,ux,ak,order);
                        splinenMW(cy,uy,ak,order);
                        break;
                }
                
                res = 0.; n1 = 1-n2;
                if (xi+n1>=0 && xi+n2<width && yi+n1>=0 && yi+n2<height) {
                    
                    int adr = yi*width+xi; 
                    for (int dy=n1;dy<=n2;dy++) 
                        for (int dx=n1;dx<=n2;dx++) 
                            res += cy[n2-dy]*cx[n2-dx]*ref[adr+width*dy+dx];
                } else 
                    
                    for (int dy=n1;dy<=n2;dy++)
                        for (int dx=n1;dx<=n2;dx++) 
                            res += cy[n2-dy]*cx[n2-dx]*vMW(ref,xi+dx,yi+dy,bg,width,height);
                
            }
            
            
        }	
        
        
        return res;
        
        
    }
	
    
	////////////////////////////////////////////////////////////////////////////////
	//! End Spline Interpolation
	////////////////////////////////////////////////////////////////////////////////
	
    
    
    
    
    
    
    
    
    
    
    
    ////////////////////////////////////////////////////////////////////////////////
	//! BEGIN Apply transformations
	////////////////////////////////////////////////////////////////////////////////
	
    
    
    
    
    void compute_planar_homography_n_points(float *x0, float *y0, float *x1, float *y1, int n, laMatrix &H)
    {
        
        ////////////////// Compute Baricenter
        float lx = 0.0,  ly = 0.0;
        float rx = 0.0,  ry = 0.0;
        
        for (int i=0; i< n; i++) {
            
            lx += x0[i]; 
            ly += y0[i]; 
            
            rx += x1[i]; 
            ry += y1[i]; 
        }
        
        lx /= (float) n; 
        ly /= (float) n; 
        rx /= (float) n; 
        ry /= (float) n; 
        
        
        /////////////// Normalize points without modifying original vectors
        
        float *px0 = new float[n];
        float *py0 = new float[n];
        
        float *px1 = new float[n];
        float *py1 = new float[n];
        
        float spl= 0.0f, spr = 0.0f;
        for(int i=0; i < n; i++)
        {
            
            px0[i] = x0[i] - lx;
            py0[i] = y0[i] - ly;
            
            px1[i] = x1[i] - rx;
            py1[i] = y1[i] - ry;
            
            spl += sqrtf(px0[i] * px0[i] + py0[i]*py0[i]);
            spr += sqrtf(px1[i] * px1[i] + py1[i]*py1[i]);
            
        }
        
        
        spl = sqrtf(2.0f) / spl;
        spr = sqrtf(2.0f) / spr;
        
        for (int i=0; i< n; i++) {
            
            px0[i] *= spl;
            py0[i] *= spl;
            
            px1[i] *= spr;
            py1[i] *= spr;
            
        }
        
        
        /////////////////////////////////////////////////////// Minimization problem || Ah ||
        
        
        laMatrix Tpl(3,3), Timinv(3,3);
        
        // similarity transformation of the plane 
        Tpl[0][0] = spl; Tpl[0][1] = 0.0; Tpl[0][2] = -spl*lx;
        Tpl[1][0] = 0.0; Tpl[1][1] = spl; Tpl[1][2] = -spl*ly;
        Tpl[2][0] = 0.0; Tpl[2][1] = 0.0; Tpl[2][2] = 1.0;
        
        // inverse similarity transformation of the image
        Timinv[0][0] = 1.0/spr; Timinv[0][1] =   0.0  ; Timinv[0][2] = rx;
        Timinv[1][0] =   0.0  ; Timinv[1][1] = 1.0/spr; Timinv[1][2] = ry;
        Timinv[2][0] =   0.0  ; Timinv[2][1] =   0.0  ; Timinv[2][2] = 1.0;
        
        
        ///////////////// System matrix
        
        laMatrix A(2*n,9);
        for(int i=0, eq=0; i < n; i++, eq++) {
            
            float	xpl = px0[i], ypl = py0[i],
            xim = px1[i], yim = py1[i];
            
            A[eq][0] = A[eq][1] = A[eq][2] = 0.0;
            A[eq][3] = -xpl;
            A[eq][4] = -ypl;
            A[eq][5] = -1.0;
            A[eq][6] =  yim * xpl;
            A[eq][7] =  yim * ypl;
            A[eq][8] =  yim;
            
            eq++;
            
            A[eq][0] =  xpl;
            A[eq][1] =  ypl;
            A[eq][2] =  1.0;
            A[eq][3] = A[eq][4] = A[eq][5] = 0.0;
            A[eq][6] = -xim * xpl;
            A[eq][7] = -xim * ypl;
            A[eq][8] = -xim;
        }
        
        
        ///////////////// SVD
        laMatrix U(2*n,9), V(9,9);
        laVector W(9);
        
        compute_svd(A,U,V,W);
        
        
        // Find the index of the least singular value
        int imin = 0;
        for (int i=1; i < 9; i++) 
            if ( W[i] < W[imin] ) imin = i;
        
        
        ////////////////// Denormalize H = Timinv * V.col(imin)* Tpl;
        laMatrix matrix(3,3), result(3,3);
        
        
        int k=0;
        for(int i=0; i < 3; i++)
            for(int j=0; j < 3; j++, k++)
                matrix[i][j] = V[k][imin];
        
        
        result = Timinv * matrix;
        H = result * Tpl;
        
        
        
        delete[] px0;
        delete[] py0;
        delete[] px1;
        delete[] py1;
        
        
    }
    
    
    
    
    
    /// Compute homography using svd + Ransac
    int compute_ransac_planar_homography_n_points(float *x0, float *y0, float *x1, float *y1, int n, int niter, float tolerance, laMatrix &H, int *accorded)
    {
        
        
        // Initialize seed		
        srand48( (long int) time (NULL) + (long int) getpid() );	
        float tolerance2 = tolerance * tolerance;
        
        
        H.create(3,3);
        laMatrix Haux(3,3);
        
        
        
        int *paccorded = new int[n];	
        int naccorded = 0;
        int pnaccorded = 0;
        
        for(int iter = 0; iter < niter; iter++)
        {
            
            
            // Choose 4 indexos from 1..n without repeated values
            int indexos[4];
            int acceptable = 0;
            while (!acceptable)
            {
                acceptable = 1;
                for(int i=0; i < 4; i++) indexos[i] = (int)  floor(drand48() * (double) n);
                
                // Check if indexos are repeated
                for(int i=0; i < 4; i++)
                    for(int j=i+1; j < 4; j++)
                        if (indexos[i] == indexos[j]) acceptable = 0; 
            }
            
            
            
            
            
            // Store selected matches 
            float px0[4] , py0[4], px1[4] , py1[4];
            for(int i=0; i < 4; i++)
            {
                px0[i] = x0[indexos[i]];
                py0[i] = y0[indexos[i]];
                
                px1[i] = x1[indexos[i]];
                py1[i] = y1[indexos[i]];
            }
            
            
            // Compute planar homography
            compute_planar_homography_n_points(px0, py0, px1, py1, 4, Haux);
            
            
            // Which matches are according to this transformation
            pnaccorded = 0;
            for(int i=0; i < n; i++)
            {
                
                laVector vec(3), res(3);
                vec[0] = x0[i];
                vec[1] = y0[i];
                vec[2] = 1.0f;
                
                res = Haux * vec;
                
                if (res[2] != 0.0f) {
                    
                    res[0] /= res[2]; res[1] /= res[2];
                    
                    float dif = (res[0] - x1[i]) * (res[0] - x1[i]) + (res[1] - y1[i]) * (res[1] - y1[i]);
                    
                    if (dif < tolerance2) {  paccorded[pnaccorded] = i; pnaccorded++; }
                    
                }
                
            }
            
            
            
            // if more according points --> save positions
            if (pnaccorded > naccorded)
            {
                
                naccorded = pnaccorded;
                
                if (accorded != NULL)
                    for(int i=0; i < naccorded; i++) accorded[i] = paccorded[i];
                
                H = Haux;
                
                
            }
            
        }
        
        
        
        return naccorded;
        
    }
    
    
    
    int compute_ransac_planar_homography_n_points(float *x0, float *y0, float *x1, float *y1, int n, int niter, float tolerance, laMatrix &H)
    {
        
        return  compute_ransac_planar_homography_n_points(x0, y0, x1, y1, n, niter, tolerance, H, NULL);
        
    }
    
    
    
    
    
    void apply_planar_homography(float *input, int width, int height, laMatrix &H, float bg, int order, float *out, float x0, float y0, int nwidth, int nheight)
    {
        
        
        //! Inverse matrix
        laMatrix V(3,3);        
        luinv(H, V);
        
        
        //! Spline interpolation
        float *coeffs;
        float *ref;
        
        float  ak[13];
        
        if (order!=0 && order!=1 && order!=-3 && 
            order!=3 && order!=5 && order!=7 && order!=9 && order!=11)
        {	
            printf("unrecognized interpolation order.\n");
            exit(-1);
        }
        
        if (order>=3) {
            
            coeffs = new float[width*height];
            finvsplineMW(input,order,coeffs,width,height);
            ref = coeffs;
            if (order>3) init_splinenMW(ak,order);
            
        } else 
        {
            coeffs = NULL;
            ref = input;
        }
        
        
        //! For each point in new image we compute its anti image and interpolate the new value
        for(int i=0; i < nwidth; i++)
            for(int j=0; j < nheight; j++)
            {
                
                
                //! Compute transformed vector
                laVector vec(3), vres(3);
                vec[0] = (float) i + x0;
                vec[1] = (float) j + y0;
                vec[2] = 1.0f;
                
                vres = V * vec;
                
                
                
                
                if (vres[2] != 0.0f) 
                {
                    
                    vres[0] /= vres[2]; vres[1] /= vres[2];
                    
                    float xp =  (float) vres[0];
                    float yp =  (float) vres[1];
                    
                    float res =  evaluate_splineMW(input,  ref,   xp,   yp, &ak[0],   order,   bg,   width,   height); 
                    
                    out[j*nwidth+i] = res;
                    
                } else
                {
                    out[j*nwidth+i] = bg;
                }    
                
                
            }
        
        if (coeffs != NULL )  delete[] coeffs;
        
    }
    
   
    
    
    
    
    void apply_planar_homography_zoom(float *input, int width, int height, laMatrix &H, float bg, int order, float *out, float x0, float y0, int nwidth, int nheight, float fZoom)
    {
        
        
        //! Inverse matrix
        laMatrix V(3,3);        
        luinv(H, V);
        
        
        //! Spline interpolation
        float *coeffs;
        float *ref;
        
        float  ak[13];
        
        if (order!=0 && order!=1 && order!=-3 && 
            order!=3 && order!=5 && order!=7 && order!=9 && order!=11)
        {	
            printf("unrecognized interpolation order.\n");
            exit(-1);
        }
        
        if (order>=3) {
            
            coeffs = new float[width*height];
            finvsplineMW(input,order,coeffs,width,height);
            ref = coeffs;
            if (order>3) init_splinenMW(ak,order);
            
        } else 
        {
            coeffs = NULL;
            ref = input;
        }
        
        
        //! For each point in new image we compute its anti image and interpolate the new value
        for(int i=0; i < nwidth; i++)
            for(int j=0; j < nheight; j++)
            {
                
                
                
                //! Compute transformed vector
                laVector vec(3), vres(3);
                vec[0] = (float) i / fZoom + x0;
                vec[1] = (float) j / fZoom + y0;
                vec[2] = 1.0f;
                
                vres = V * vec;

                
                if (vres[2] != 0.0f) 
                {
                    
                    vres[0] /= vres[2]; vres[1] /= vres[2];
                    
                    float xp =  (float) vres[0];
                    float yp =  (float) vres[1];
                    
                    float res =  evaluate_splineMW(input,  ref,   xp,   yp, &ak[0],   order,   bg,   width,   height); 
                    
                    out[j*nwidth+i] = res;
                    
                } else
                {
                    out[j*nwidth+i] = bg;
                }    
                
                
            }
        
        if (coeffs != NULL) delete[] coeffs;
        
        
    }
    
    
    
    
    
    void apply_zoom(float *input, float *out, float zoom, int order, int width, int height)
    {
        
        int nwidth = (int)rintf( zoom * (float) width);
        int nheight = (int)rintf( zoom * (float) height);
        
    
        
        laMatrix H(3,3);
        H[0][0] = (double) zoom; H[0][1] = 0.0; H[0][2] = 0.0;
        H[1][0] = 0.0; H[1][1] = (double) zoom; H[1][2] = 0.0;
        H[2][0] = 0.0; H[2][1] = 0.0; H[2][2] = 1.0;
        
        
        float bg=-1;
        apply_planar_homography(input,  width,  height,  H,  bg,  order, out, 0.0, 0.0,  nwidth,  nheight);
        
        
    }
    
    
    
    
    
    
    
    
    
    
    
    ////////////////////////////////////////////////////////////////////////////////
	//! END Apply transformations
	////////////////////////////////////////////////////////////////////////////////
    
    
    
    
    
    
    
    
    
    
    flData3D::flData3D() : npoints(0), dx(NULL), dy(NULL), dz(NULL), dr(NULL), dg(NULL), db(NULL),  ok(false), color(false)
    {
        
    }
	
    
    
    flData3D::flData3D(flData3D &inData) 
        : npoints(0), dx(NULL), dy(NULL), dz(NULL), dr(NULL), dg(NULL), db(NULL),  ok(false), color(false)
    {
        
        
        
        if (inData.npoints > 0){
            
            npoints = inData.npoints;
            allocate_coordinates(npoints);
            
            color = inData.color;
            if (color) 
                allocate_color(npoints);
            else{
                dr = dg = db = NULL;
            }
            
            for(int i = 0; i < npoints; i++)
            {
                
                dx[i] = inData.dx[i]; dy[i] = inData.dy[i]; dz[i] = inData.dz[i];
                
                if (color){dr[i] = inData.dr[i]; dg[i] = inData.dg[i]; db[i] = inData.db[i];}
            }
            
            ok = true;
        }
        
    }
    
    
    
    flData3D::flData3D(int inpoints, float *idx, float *idy, float *idz, float *idr, float *idg, float *idb) 
        :  npoints(0), dx(NULL), dy(NULL), dz(NULL), dr(NULL), dg(NULL), db(NULL),ok(false), color(false)
    {
        
        
        if (inpoints > 0)
        {
            
            npoints = inpoints;
            allocate_coordinates(npoints);
            
            if (idr && idg && idb) color = true;
            else color = false;
            
            if (color) 
                allocate_color(npoints);
            else{
                dr = dg = db = NULL;
            }
            
            for(int i = 0; i < npoints; i++)
            {
                
                dx[i] = idx[i]; dy[i] = idy[i]; dz[i] = idz[i];
                
                if (color){dr[i] = idr[i]; dg[i] = idg[i]; db[i] = idb[i];}
            }
            
            ok = true;
        }
        
    }
    
    
    
    
    
    
    flData3D::flData3D(const char * filename)
            :  npoints(0), dx(NULL), dy(NULL), dz(NULL), dr(NULL), dg(NULL), db(NULL),ok(false), color(false)
    {
        npoints = 0;
        dx = dy = dz = NULL;
        dr = dg = db = NULL;
        ok = false;
        color = false;
        
        loadFile(filename);
        
    }
    
    
    
    void flData3D::allocate_coordinates(int n)
    {
        
        if (!dx) dx = new float[n];
        if (!dy) dy = new float[n];
        if (!dz) dz = new float[n];
        
    }
    
    void flData3D::allocate_color(int n)
    {
        
        if (!dr) dr = new float[n];
        if (!dg) dg = new float[n];
        if (!db) db = new float[n];
        
    }
    
    void flData3D::desallocate_coordinates()
    {
        
        if (!dx) {delete[] dx; dx = NULL;}
        if (!dy) {delete[] dy; dy = NULL;}
        if (!dz) {delete[] dz; dz = NULL;}
    }
    
    
    
    void flData3D::desallocate_color()
    {
        
        if (!dr) {delete[] dr; dr = NULL;}
        if (!dg) {delete[] dg; dg = NULL;}
        if (!db) {delete[] db; db = NULL;}
    }
    
    
    
    
    
    
    
	
    
    flData3D::~flData3D()
    {
        
        if (dx) delete[] dx;
        if (dy) delete[] dy;
        if (dz) delete[] dz;
        
        if (dr) delete[] dr;
        if (dg) delete[] dg;
        if (db) delete[] db;
        
    }
    
    
    
    flData3D& flData3D::operator= (const flData3D &in)
    {
        
        if (npoints != in.npoints)
        {
            desallocate_color();
            desallocate_coordinates();
            
        }
        
        npoints = in.npoints;
        
        if (color && !in.color) desallocate_color();
        
        allocate_coordinates(npoints); // It only allocates if pointers are NULL
        if (color) allocate_color(npoints);
        
        for(int i = 0; i < npoints; i++)
        {
            dx[i] = in.dx[i]; dy[i] = in.dy[i]; dz[i] = in.dz[i];
            if (color){dr[i] = in.dr[i]; dg[i] = in.dg[i]; db[i] = in.db[i];}	
        }
        
        return *this;
        
        
    }
    
    
    
    
    void flData3D::loadFile(const char * filename)
    {
        
        std::ifstream file(filename);		
        if (!file.is_open())
        {
            printf("Exit(flData3D::loadData(wxString filename)): Surface file not found or impossible to open\n");
            ok = false;
            return;
            
        }
        
		
        char  iinfo[256];
        file >> iinfo;
		
        if( strcmp(iinfo, "flData3D") == 0)
        {
            
            if (dx) delete[] dx;
            if (dy) delete[] dy;
            if (dz) delete[] dz;
            
            if (dr) delete[] dr;
            if (dg) delete[] dg;
            if (db) delete[] db;
            
            
            file >> npoints;
			
            dx = new float[npoints];
            dy = new float[npoints];
            dz = new float[npoints];
            
            file >>	iinfo;
            
            if (strcmp(iinfo,"color") != 0)
            {
                printf("Warning (flData3D::loadData(wxString filename)): Color flag not found, not a correct 3ddata file.\n");
                ok = false;
                return;
            }
            
            int colorflag;
            file >> colorflag;
            
            if (colorflag) color = true;
            else color = false;
            
            if (color)
            {		
                dr = new float[npoints];
                dg = new float[npoints];
                db = new float[npoints];
            } else
            {
                dr = dg = db = NULL;
            }
			
			
            for(int i=0; i < npoints; i++)
            {
                file >> dx[i];
                file >> dy[i];
                file >> dz[i];
                
                
                if (color)
                {
                    file >> dr[i];
                    file >> dg[i];
                    file >> db[i];
                    
                }		
            }
            
            ok = true;
			
        } else
        {
			
			printf("Error: not a correct flData3D file\n");
			ok = false;
			return ;
			
        }
		
        
    }
    
    
    
    
    
    
    
    
    void flData3D::SaveFile(const char* filename)
    {
        
        std::ofstream outfile(filename);
        
        if (!outfile.is_open())
        {
            printf("Data file not writtable\n");
            exit(-1); 
        }
        
        outfile << "flData3D" << std::endl;
        outfile << npoints << std::endl;
        outfile << "color" << std::endl;
        if (color) outfile << "1" << std::endl;
        else outfile << "0" << std::endl;
        
        
        for(int i=0; i < npoints; i++)
        {
            
            if (color)
                outfile << dx[i] << "  " << dy[i] << " " << dz[i] << " " << dr[i]  << "  " << dg[i] << " " << db[i] << std::endl;
            else
                outfile << dx[i] << "  " << dy[i] << " " << dz[i] << std::endl;
            
        }
        
    }
    
    
    
    
    
    
    void flData3D::normalize()
    {
        
        
        // Move to baricenter
        float bx, by,bz;
        
        bx=by=bz=0.0;
        for(int i=0; i < npoints; i++)
        {
            bx += dx[i];
            by += dy[i];
            bz += dz[i];
            
        }	
        
        bx /= npoints;
        by /= npoints;
        bz /= npoints;
        
        
        for(int i=0; i < npoints; i++)
        {
            
            dx[i] -= bx;
            dy[i] -= by;
            dz[i] -= bz;
            
        }
        
        
        
        // Changing range to [-0.5,0.5]
        float max = fabsf(dx[0]);
        for(int i=0; i < npoints; i++)
        {
            
            if (fabsf(dx[i]) > max)  max = fabsf(dx[i]);			
            if (fabsf(dy[i]) > max)  max = fabsf(dy[i]);			
            if (fabsf(dz[i]) > max)  max = fabsf(dz[i]);			
            
        }
        
        max = 2.0*max;
        
        for(int i=0; i < npoints; i++)
        {
            
            dx[i] /= max;
            dy[i] /= max;
            dz[i] /= max;
            
        }
        
        
    }
    
    /*
     
     
     void flData3D::getNormParam(float *ibx, float *iby,float *ibz, float *iamp)
     {
     
     // Move to baricenter
     float bx, by,bz;
     
     bx=by=bz=0.0;
     for(int i=0; i < npoints; i++)
     {
     bx += dx[i];
     by += dy[i];
     bz += dz[i];
     
     }	
     
     bx /= npoints;
     by /= npoints;
     bz /= npoints;
	 
     *ibx = bx;
     *iby = by;
     *ibz = bz;
     
     
     // Changing range to [-0.5,0.5]
     float max = fabsf(dx[0] - bx);
     for(int i=0; i < npoints; i++)
     {
     
     if (fabsf(dx[i] - bx) > max)  max = fabsf(dx[i] - bx);			
     if (fabsf(dy[i] - by) > max)  max = fabsf(dy[i] - by);			
     if (fabsf(dz[i] - bz) > max)  max = fabsf(dz[i] - bz);			
     
     }
     
     *iamp = max;
     
     }
     
     
     
     void flData3D::normalize(float bx, float by, float bz, float amp)
     {
     
     for(int i=0; i < npoints; i++)
     {
     
     dx[i] -= bx;
     dy[i] -= by;
     dz[i] -= bz;
     
     }
     
     
     amp = 2.0*amp;
     
     for(int i=0; i < npoints; i++)
     {
     
     dx[i] /= amp;
     dy[i] /= amp;
     dz[i] /= amp;
     
     }
     
     }
     
     
     */
    
	
	
    
    
    
    
    ////////////////////////////////////////////////////////////////////////////////
	//! BEGIN Numerics
	////////////////////////////////////////////////////////////////////////////////
	
    
    
    
    
    
    
    laVector::laVector() : d_n(0), d_v(0) {}
    
    
    
    laVector::laVector(int n) : d_n(n), d_v(new double[n]) {}
    
    
    laVector::laVector(const double& a, int n) : d_n(n), d_v(new double[n])
    {
        for(int i=0; i<n; i++)
            d_v[i] = a;
    }
    
    
    laVector::laVector(const double *a, int n) : d_n(n), d_v(new double[n])
    {
        for(int i=0; i<n; i++)
            d_v[i] = *a++;
    }
    
    
    laVector::laVector(const laVector &rhs) : d_n(rhs.d_n), d_v(new double[d_n])
    {
        for(int i=0; i<d_n; i++)
            d_v[i] = rhs[i];
    }
    
    
    laVector & laVector::operator=(const laVector &rhs)
    {
        if (this != &rhs)
        {
            if (d_n != rhs.d_n) {
                if (d_v != 0) delete [] d_v;
                d_n =rhs.d_n;
                d_v = new double[d_n];
            }
			
            for (int i=0; i<d_n; i++)
                d_v[i]=rhs[i];
        }
        return *this;
    }
    
    
    laVector & laVector::operator=(const double &a)	//assign a to every element
    {
        for (int i=0; i<d_n; i++)
            d_v[i]=a;
        return *this;
    }
    
    
	double & laVector::operator[](const int i)	//subscripting
    {
        return d_v[i];
    }
    
    
	const double & laVector::operator[](const int i) const	//subscripting
    {
        return d_v[i];
    }
    
    
	double * laVector::v() {return d_v;}
	
    
    
	int laVector::size() const
    {
        return d_n;
    }
    
    
    laVector::~laVector()
    {
        if (d_v != 0)
            delete[] d_v;
    }
    
    
    
    
    
    
    laMatrix::laMatrix() : d_n(0), d_m(0), d_v(0) {}
    
    
    laMatrix::laMatrix(int n, int m) : d_n(n), d_m(m), d_v(new double*[n])
    {
        d_v[0] = new double[m*n];
        for (int i=1; i< n; i++)
            d_v[i] = d_v[i-1] + m;
    }
    
    
    
    void laMatrix::create(int n, int m)
    {
        
        if (d_v != 0) {
            delete[] d_v[0];
            delete[] d_v;
        }
        
        d_n = n;
        d_m = m;
        
        d_v = new double*[n];
        d_v[0] = new double[m*n];
        for (int i=1; i< n; i++)
            d_v[i] = d_v[i-1] + m;
    }
    
    
    
    
    
    laMatrix::laMatrix(const double &a, int n, int m) : d_n(n), d_m(m), d_v(new double*[n])
    {
        int i,j;
        d_v[0] = new double[m*n];
        for (i=1; i< n; i++)
            d_v[i] = d_v[i-1] + m;
        for (i=0; i< n; i++)
            for (j=0; j<m; j++)
                d_v[i][j] = a;
    }
    
    
    laMatrix::laMatrix(const double *a, int n, int m) : d_n(n), d_m(m), d_v(new double*[n])
    {
        int i,j;
        d_v[0] = new double[m*n];
        for (i=1; i< n; i++)
            d_v[i] = d_v[i-1] + m;
        for (i=0; i< n; i++)
            for (j=0; j<m; j++)
                d_v[i][j] = *a++;
    }
    
    
    laMatrix::laMatrix(const laMatrix &rhs) : d_n(rhs.d_n), d_m(rhs.d_m), d_v(new double*[rhs.d_n])
    {
        assert(d_m>0 && d_n>0);
        
        d_v[0] = new double[d_m * d_n];
        
        for (int i=1; i< d_n; i++)
            d_v[i] = d_v[i-1] + d_m;
        
        for (int i=0; i< d_n; i++)
            for (int j=0; j<d_m; j++)
                d_v[i][j] = rhs[i][j];
    }
    
    
    laMatrix & laMatrix::operator=(const laMatrix &rhs)
    {
        if (this != &rhs) {
            int i,j;
            if (d_n != rhs.d_n || d_m != rhs.d_m) {
                if (d_v != 0) {
                    delete[] (d_v[0]);
                    delete[] (d_v);
                }
                
                d_n=rhs.d_n;
                d_m=rhs.d_m;
                d_v = new double*[d_n];
                d_v[0] = new double[d_m*d_n];
                
                
            }
            
            
            for (i=1; i< d_n; i++)
                d_v[i] = d_v[i-1] + d_m;
            for (i=0; i< d_n; i++)
                for (j=0; j<d_m; j++)
                    d_v[i][j] = rhs[i][j];
        }
        return *this;
    }
    
    
    laMatrix & laMatrix::operator=(const double &a)	//assign a to every element
    {
        for (int i=0; i< d_n; i++)
            for (int j=0; j<d_m; j++)
                d_v[i][j] = a;
        return *this;
    }
    
    
	
	double* laMatrix::operator[](const int i)	//subscripting: pointer to row i
    {
        return d_v[i];
    }
    
    
	const double* laMatrix::operator[](const int i) const
    {
        return d_v[i];
    }
    
    
	double ** laMatrix::v() {return d_v;}
    
	
    
    
    
    int laMatrix::nrows() const
    {
        return d_n;
    }
    
    
    int laMatrix::ncols() const
    {
        return d_m;
    }
    
    
    laMatrix::~laMatrix()
    {
        if (d_v != 0) {
            delete[] (d_v[0]);
            delete[] (d_v);
        }
    }
    
    
    
    
    
    
    
    laMatrix operator*(double a, const laMatrix& rhs)
    {
        
        laMatrix r(rhs.d_n, rhs.d_m);
        
        for (int k=r.d_n * r.d_m - 1; k>=0; k--)
            r.d_v[0][k] = a * rhs.d_v[0][k];
        
        return r;
    }
    
    
    
    laMatrix operator/(const laMatrix& lhs, double a)
    {
        laMatrix r(lhs.d_n, lhs.d_m);
        
        for (int k=r.d_n * r.d_m - 1; k>=0; k--)
            r.d_v[0][k] = lhs.d_v[0][k] / a;
        return r;
        
    }
    
    
    
    laMatrix operator+(const laMatrix& lhs, const laMatrix& rhs)
    {
        laMatrix r(lhs.d_n, lhs.d_m);
        for (int k=r.d_m * r.d_n -1; k>=0; k--)
            r.d_v[0][k] = lhs.d_v[0][k] + rhs.d_v[0][k];
        return r;
    }
    
    
    
    laMatrix operator-(const laMatrix& lhs, const laMatrix& rhs)
    {
        laMatrix r(lhs.d_m, lhs.d_n);
        for (int k=r.d_m * r.d_n - 1; k>=0; k--)
            r.d_v[0][k] = lhs.d_v[0][k] - rhs.d_v[0][k];
        return r;
    }
    
    
    
    laMatrix operator*(const laMatrix& lhs, const laMatrix& rhs)
    {
        double aux;
        
        laMatrix r(lhs.d_n, rhs.d_m);
        for (int i=0; i< r.d_n; i++)
            for (int j=0; j< r.d_m; j++) 
            {
                aux = 0.0;
                for (int k=0; k< lhs.d_m; k++)
                    aux += lhs.d_v[i][k] * rhs.d_v[k][j];
                
                r.d_v[i][j] = aux;
            }
        
        return r;
    }
    
    
    
    laVector operator*(const laMatrix& lhs, const laVector& rhs)
    {
        
		
        laVector r(lhs.d_n);
        for (int i=0; i < r.size(); i++) 
        {
            
            r[i] = 0;
            for (int k=0; k< rhs.size(); k++)
            {   r[i] += lhs.d_v[i][k] * rhs[k];}
        }
        
        return r;
    }
    
	
	
	laMatrix laMatrix::transposed()
	{
		laMatrix r(d_m, d_n);
		
		for (int ii=0; ii < d_m; ii++)
			for (int jj=0; jj < d_n; jj++)
            {
                r.d_v[ii][jj] = d_v[jj][ii];
            }
        
		return r;
	}
	
    
    laMatrix laMatrix::copyBlock(int i0, int j0, int rowb, int colb) 
    {
        
        laMatrix block(rowb, colb);
        
        for (int i=0; i < rowb; i++)
            for (int j=0; j < colb; j++)
                block.d_v[i][j] = d_v[i0+i][j0+j];
        
        return block;
    }
    
    
    
    void ludcmp(laMatrix &a, laVector &indx, double &d)
    {
        const double TINY=1.0e-20;
        int i,imax,j,k;
        double big,dum,sum,temp;
        
        int n=a.nrows();
        laVector vv(n);
        d=1.0;
        
        for (i=0;i<n;i++) {
            big=0.0;
            for (j=0;j<n;j++)
                if ((temp=fabs(a[i][j])) > big) big=temp;
            if (big == 0.0) {printf("Singular matrix in routine ludcmp"); exit(-1);}
            vv[i]=1.0/big;
        }
        for (j=0;j<n;j++) {
            for (i=0;i<j;i++) {
                sum=a[i][j];
                for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
                a[i][j]=sum;
            }
            big=0.0;
            for (i=j;i<n;i++) {
                sum=a[i][j];
                for (k=0;k<j;k++) sum -= a[i][k]*a[k][j];
                a[i][j]=sum;
                if ((dum=vv[i]*fabs(sum)) >= big) {
                    big=dum;
                    imax=i;
                }
            }
            if (j != imax) {
                for (k=0;k<n;k++) {
                    dum=a[imax][k];
                    a[imax][k]=a[j][k];
                    a[j][k]=dum;
                }
                d = -d;
                vv[imax]=vv[j];
            }
            indx[j]=(double) imax;
            if (a[j][j] == 0.0) a[j][j]=TINY;
            if (j != n-1) {
                dum=1.0/(a[j][j]);
                for (i=j+1;i<n;i++) a[i][j] *= dum;
            }
        }
    }
    
    
    /* Solves the set of n linear equations Ax=b. Here a[0..n-1][0..n-1] as input, not as the matrix A but rather as its LU decomposition,*/
	/* determined by the routine ludcmp. indx[0..n-1] is input as the permutation vector returned by ludcmp. b[0..n-1] is input as the    */
	/* right hand side vector and returns with the solution vector x. */
	void lubksb(laMatrix &a, laVector &indx, laVector &b)
    {
        int i,ii=0,ip,j;
        double sum;
        
        int n=a.nrows();
        for (i=0;i<n;i++) {
            ip = (int) indx[i];
            sum=b[ip];
            b[ip]=b[i];
            if (ii != 0)
                for (j=ii-1;j<i;j++) sum -= a[i][j]*b[j];
            else if (sum != 0.0)
                ii=i+1;
            b[i]=sum;
        }
        for (i=n-1;i>=0;i--) {
            sum=b[i];
            for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
            b[i]=sum/a[i][i];
        }
    }
    
    
    
    // Solves Ax=b by using lu decomposition 
	void lusolve(laMatrix &a, laVector &x, laVector &b)
	{
		
        int n = b.size();
        
        laVector indx(n);
		laMatrix Aux = a;
        
        double d;
		ludcmp(Aux,indx,d);
        
        
        x=b;
        lubksb(Aux,indx,x);
        
	}
    
    
    
    
    // Computes the inverse of a column by column using the LU decomposition
	void luinv(laMatrix &a, laMatrix &inv)
	{
        
        
        int n = a.nrows();
		laVector col(n), indx(n);
        
        laMatrix invaux=a;
        inv=a;
        
        double d;
		ludcmp(invaux,indx,d);
        
        
        for(int j=0; j<n; j++)
        {
            
            for(int i=0; i<n; i++)  col[i]=0.0;
            col[j]=1.0;
            
            lubksb(invaux,indx,col);
            
            for(int i=0;i<n;i++) {inv[i][j] = col[i];}
			
			
		} 
        
	}
    
    
    
    void invCC(laMatrix &inv)
    {
        
        double z,*p,*q,*r,*s,*t; 
        int n= inv.nrows();
        int j,k;
        double *v = inv[0];
        
        for(j=0,p=v; j<n ;++j,p+=n+1){
            for(q=v+j*n; q<p ;++q) *p-= *q* *q;
            
            if(*p<=0.) { printf("error :: non defite positive matrix\n"); return;}
            
            *p=sqrt(*p);
            for(k=j+1,q=p+n; k<n ;++k,q+=n){
                for(r=v+j*n,s=v+k*n,z=0.; r<p ;) z+= *r++ * *s++;
                *q-=z; *q/= *p;
            }
            
        }
        
        
        inv = inv.transposed();
        
        for(j=0,p=v; j<n ;++j,p+=n+1){ *p=1./ *p;
            for(q=v+j,t=v; q<p ;t+=n+1,q+=n){
                for(s=q,r=t,z=0.; s<p ;s+=n) z-= *s * *r++;
                *q=z* *p; }
        }
        for(j=0,p=v; j<n ;++j,p+=n+1){
            for(q=v+j,t=p-j; q<=p ;q+=n){
                for(k=j,r=p,s=q,z=0.; k<n ;++k) z += *r++ * *s++;
                *t++ =(*q=z); }
        }
        
    }
    
    
    
    
    double withSignOf(double a, double b)
	{ return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a); }
	
	double svdhypot(double a, double b)
	{
		a = fabsf(a);
		b = fabsf(b);
		if(a > b) {
			b /= a;
			return a*sqrt(1.0 + b*b);
		} else if(b) {
			a /= b;
			return b*sqrt(1.0 + a*a);
		}
		return 0.0;
	}
	
    
	void svdrotate(double& a, double& b, double c, double s)
	{
		double d = a;
		a = +d*c +b*s;
		b = -d*s +b*c;
	}
    
    
    
    
    
    
    
	
    
    
    
    void compute_svd(laMatrix &A, laMatrix &m_U, laMatrix &m_V, laVector &m_W)
	{
        
        int rows = A.nrows();
        int cols = A.ncols();
        
		const double	EPSILON = 0.00001;
		const int SVD_MAX_ITS = 100;
		
		double g, scale, anorm;
		double * RV1 = new double[cols];
		
		
		for(int i=0; i < rows; i++)
			for(int j=0; j < cols; j++)
				m_U[i][j] = A[i][j];		
		
		// Householder reduction to bidiagonal form:
		anorm = g = scale = 0.0;
		for (int i=0; i< cols; i++) 
		{
			int l = i + 1;
			RV1[i] = scale*g;
			g = scale = 0.0;
			if(i < rows) 
			{
				for (int k=i; k< rows; k++)
					scale += fabsf(m_U[k][i]);
				if (scale != 0.0) {
					double invScale=1.0/scale, s=0.0;
					for (int k=i; k< rows; k++) {
						m_U[k][i] *= invScale;
						s += m_U[k][i] * m_U[k][i];
					}
					double f = m_U[i][i];
					g = - withSignOf(sqrt(s),f);
					double h = 1.0 / (f*g - s);
					m_U[i][i] = f - g;
					for (int j=l; j< cols; j++) {
						s = 0.0;
						for (int k=i; k< rows; k++)
							s += m_U[k][i] * m_U[k][j];
						f = s * h;
						for (int k=i; k< rows; k++)
							m_U[k][j] += f * m_U[k][i];
					}
					for (int k=i; k< rows; k++)
						m_U[k][i] *= scale;
				}
			}
			
			
			m_W[i] = scale * g;
			g = scale = 0.0;
			if ( i< rows && i< cols-1 ) {
				for (int k=l; k< cols; k++)
					scale += fabsf(m_U[i][k]);
				if (scale != 0.0) {
					double invScale=1.0/scale, s=0.0;
					for (int k=l; k< cols; k++) {
						m_U[i][k] *= invScale;
						s += m_U[i][k] * m_U[i][k];
					}
					double f = m_U[i][l];
					g = - withSignOf(sqrt(s),f);
					double h = 1.0 / (f*g - s);
					m_U[i][l] = f - g;
					for (int k=l; k< cols; k++)
						RV1[k] = m_U[i][k] * h;
					for (int j=l; j< rows; j++) {
						s = 0.0;
						for (int k=l; k< cols; k++)
							s += m_U[j][k] * m_U[i][k];
						for (int k=l; k< cols; k++)
							m_U[j][k] += s * RV1[k];
					}
					for (int k=l; k< cols; k++)
						m_U[i][k] *= scale;
				}
			}
			anorm = MAX(anorm, fabsf(m_W[i]) + fabsf(RV1[i]) );
		}
		
		// Accumulation of right-hand transformations:
		m_V[cols-1][cols-1] = 1.0;
		for (int i= cols-2; i>=0; i--) {
			m_V[i][i] = 1.0;
			int l = i+1;
			g = RV1[l];
			if (g != 0.0) {
				double invgUil = 1.0 / (m_U[i][l]*g);
				for (int j=l; j< cols; j++)
					m_V[j][i] = m_U[i][j] * invgUil;
				for (int j=l; j< cols; j++){
					double s = 0.0;
					for (int k=l; k< cols; k++)
						s += m_U[i][k] * m_V[k][j];
					for (int k=l; k< cols; k++)
						m_V[k][j] += s * m_V[k][i];
				}
			}
			for (int j=l; j< cols; j++)
				m_V[i][j] = m_V[j][i] = 0.0;
		}
		
		// Accumulation of left-hand transformations:
		for (int i=MIN(rows,cols)-1; i>=0; i--) {
			int l = i+1;
			g = m_W[i];
			for (int j=l; j< cols; j++)
				m_U[i][j] = 0.0;
			if (g != 0.0) {
				g = 1.0 / g;
				double invUii = 1.0 / m_U[i][i];
				for (int j=l; j< cols; j++) {
					double s = 0.0;
					for (int k=l; k< rows; k++)
						s += m_U[k][i] * m_U[k][j];
					double f = (s * invUii) * g;
					for (int k=i; k< rows; k++)
						m_U[k][j] += f * m_U[k][i];
				}
				for (int j=i; j< rows; j++)
					m_U[j][i] *= g;
			} else
				for (int j=i; j< rows; j++)
					m_U[j][i] = 0.0;
			m_U[i][i] = m_U[i][i] + 1.0;
		}
		
		// Diagonalization of the bidiagonal form:
		for (int k=cols-1; k>=0; k--) { // Loop over singular values
			for (int its=1; its<=SVD_MAX_ITS; its++) {
				bool flag = false;
				int l  = k;
				int nm = k-1;
				while(l>0 && fabsf(RV1[l]) > EPSILON*anorm) { // Test for splitting
					if(fabsf(m_W[nm]) <= EPSILON*anorm) {
						flag = true;
						break;
					}
					l--;
					nm--;
				}
				if (flag) {	// Cancellation of RV1[l], if l > 0
					double c=0.0, s=1.0;
					for (int i=l; i< k+1; i++) {
						double f = s * RV1[i];
						RV1[i] = c * RV1[i];
						if (fabsf(f)<=EPSILON*anorm)
							break;
						g = m_W[i];
						double h = svdhypot(f,g);
						m_W[i] = h;
						h = 1.0 / h;
						c = g * h;
						s = - f * h;
						for (int j=0; j< rows; j++)
							svdrotate(m_U[j][nm],m_U[j][i], c,s); 
					}
				}
				double z = m_W[k];
				if (l==k) {		// Convergence of the singular value
					if (z< 0.0) {	// Singular value is made nonnegative
						m_W[k] = -z;
						for (int j=0; j< cols; j++)
							m_V[j][k] = - m_V[j][k];
					}
					break;
				}
				// Exception if convergence to the singular value not reached:
				if(its==SVD_MAX_ITS) {printf("svd::convergence_error\n"); exit(-1);}
				double x = m_W[l]; // Get QR shift value from bottom 2x2 minor
				nm = k-1;
				double y = m_W[nm];
				g = RV1[nm];
				double h = RV1[k];
				double f = ( (y-z)*(y+z) + (g-h)*(g+h) ) / ( 2.0*h*y );
				g = svdhypot(f,1.0);
				f = ( (x-z)*(x+z) + h*(y/(f+withSignOf(g,f)) - h) ) / x;
				// Next QR transformation (through Givens reflections)
				double c=1.0, s=1.0;
				for (int j=l; j<=nm; j++) {
					int i = j+1;
					g = RV1[i];
					y = m_W[i];
					h = s * g;
					g = c * g;
					z = svdhypot(f,h);
					RV1[j] = z;
					z = 1.0 / z;
					c = f * z;
					s = h * z;
					f = x*c + g*s;
					g = g*c - x*s;
					h = y * s;
					y *= c;
					for(int jj=0; jj < cols; jj++)
						svdrotate(m_V[jj][j],m_V[jj][i], c,s);
					z = svdhypot(f,h);
					m_W[j] = z;
					if (z!=0.0) { // Rotation can be arbitrary if z = 0.0
						z = 1.0 / z;
						c = f * z;
						s = h * z;
					}
					f = c*g + s*y;
					x = c*y - s*g;
					for(int jj=0; jj < rows; jj++)
						svdrotate(m_U[jj][j],m_U[jj][i], c,s);
				}
				RV1[l] = 0.0;
				RV1[k] = f;
				m_W[k] = x;
			}
		}
	}
	
    
    
    void compute_pca_svd(laMatrix &X, laVector &S, laMatrix &V, laMatrix &U)
	{
		
        int n = X.nrows();
        int p = X.ncols();
        
		compute_svd(X,U,V,S);
		
		
		// U contain new coefficients, which must be normalized by eigenvalues
		for(int i=0; i < n; i++)
			for(int j=0; j < p; j++)
				U[i][j] *= S[j];
		
		// Normalize eigenvalues
		double norm = (double) (n-1);
		for(int i=0; i < p; i++)
			S[i] = S[i] * S[i] / norm;
		
		// If n < p, principal component should be zero from n to p-1			
		// Coefficients of these principal components should be zero
		if (n < p)
		{
			for(int i=n-1; i < p; i++) S[i] = 0.0f;
			
			for(int j=0; j < n; j++)
				for(int i=n-1; i < p; i++)
					U[j][i] = 0.0f;
		}
		
		
	}
    
    
    
    
    ////////////////////////////////////////////////////////////////////////////////
	//! END Numerics
	////////////////////////////////////////////////////////////////////////////////
    
    
    
    
    
    
    
    
    
    
    
    
    
	
	//////////////////////////////////////////////////////////////////
	//! Begin Cflimage
	//////////////////////////////////////////////////////////////////
	
	//! Constructors
	cflimage::cflimage() : d_c(0), d_w(0), d_h(0), d_wh(0), d_whc(0), d_v(0),  visuMin(0.0f), visuMax(255.f)
	{
	}	
	
	
	
	cflimage::cflimage(int w, int h, int c) : d_c(c), d_w(w), d_h(h),d_wh(w*h), d_whc(c*w*h), d_v(new float[c*w*h]),  visuMin(0.0f), visuMax(255.f)
	{
		for (int ii=0; ii < d_whc; ii++) d_v[ii] = 0.0f;
		
	}	
	
	
	
	cflimage::cflimage(int w, int h, float *igray) : d_c(1), d_w(w), d_h(h),d_wh(w*h), d_whc(w*h), d_v(new float[w*h]), visuMin(0.0f), visuMax(255.f)
	{
		memcpy(d_v, igray, w * h * sizeof(float));
	}
	
	
	cflimage::cflimage(int w, int h, float *ired, float *igreen, float *iblue) : d_c(3), d_w(w), d_h(h),d_wh(w*h), d_whc(3*w*h), d_v(new float[3*w*h]),  visuMin(0.0f), visuMax(255.f)
	{
		memcpy(d_v, ired, w * h * sizeof(float));
		memcpy(d_v + w*h, igreen, w * h * sizeof(float));
		memcpy(d_v + 2*w*h, iblue, w * h * sizeof(float));
	}
	
	
	
	cflimage::cflimage(flimage &red, flimage &green, flimage &blue) : d_c(0), d_w(0), d_h(0), d_wh(0), d_whc(0), d_v(0),  visuMin(0.0f), visuMax(255.f)
	{
		
		assert(red.d_w == green.d_w && green.d_w == blue.d_w);
		assert(red.d_h == green.d_h && green.d_h == blue.d_h);
		
		d_w = red.d_w;
		d_h = red.d_h;
		d_c = 3;
		
		d_wh = d_w * d_h;
		d_whc = d_wh * d_c;
		
		d_v = new float[3 * d_wh];
		memcpy(d_v, red.d_v, d_wh * sizeof(float));
		memcpy(d_v + d_wh, green.d_v, d_wh * sizeof(float));
		memcpy(d_v + 2 * d_wh, blue.d_v, d_wh * sizeof(float));
		
		
	}
	
	
	cflimage::cflimage(const cflimage& im) : d_c(im.d_c), d_w(im.d_w), d_h(im.d_h),d_wh(im.d_wh), d_whc(im.d_whc), d_v(0),  visuMin(0.0f), visuMax(255.f)
	{
		
		if (d_whc > 0)
		{
			d_v = new float[d_whc];
			memcpy((void *) d_v, (const void *) im.d_v, d_whc * sizeof(float));
			
			//for (int ii=0; ii < d_whc; ii++) d_v[ii] = im.d_v[ii];
		}
		
	}
	
	
	
	void cflimage::create(int w, int h, int c)
	{	
		erase();
		d_c = c; d_w = w; d_h=h; d_wh = w*h; d_whc = c*w*h;
		d_v = new float[d_whc];
		visuMin=0.0f;
		visuMax=255.0f;
		
		for (int ii=0; ii < d_whc; ii++) d_v[ii] = 0.0f;	
	}
	
	
	
	void cflimage::erase() 
	{
		d_w = d_h = d_wh = d_whc = 0; 
		if (d_v) delete[] d_v; 
		d_v=0;
	} 
	
	
	
	
	cflimage::~cflimage()
	{
		erase();
	}
	
	
	
    
	
	
	
	//////////////////////////////////////////////////////////////////
	//! Begin Operators
	//////////////////////////////////////////////////////////////////
	
	
	
	
	cflimage&  cflimage::operator= (const cflimage& im)
	{
		if (&im == this)
		{ 
			return *this;
		}
		
		
		if (d_c != im.d_c || d_w != im.d_w || d_h != im.d_h)
		{ 			
			erase();
			d_c = im.d_c; d_w = im.d_w; d_h=im.d_h; d_wh = d_w * d_h; d_whc=d_c * d_w * d_h;
			d_v = new float[d_whc];
		}
		
		
		memcpy((void *) d_v, (const void *) im.d_v, d_whc * sizeof(float));
		
		visuMin = im.visuMin;
		visuMax = im.visuMax;
		
		return *this;	
		
	}
	
	
	
	
    
	
	
	//////////////////////////////////////////////////////////////////
	//! Begin value operations
	//////////////////////////////////////////////////////////////////
	
	void  cflimage::addGaussianNoise(float std)
	{
		assert(d_v != NULL);
		fpAddNoiseGaussian(d_v, d_v, std, 0, d_whc);
	}
	
    
    void  cflimage::addGaussianNoiseSD(float std, float sdstd)
	{
		assert(d_v != NULL);
		fpAddNoiseGaussianAfine(d_v, d_v, std, sdstd, 0, d_whc);
	}
	
    
    
    
	
	cflimage  cflimage::gradient(char cType)
	{
		assert(d_v != NULL);
        
		cflimage grad(d_w, d_h, d_c);
		
		for (int ii=0; ii < d_c; ii++)
		{
			fiComputeImageGradient(&d_v[ii * d_wh], &grad.d_v[ii * d_wh], NULL, d_w, d_h,cType);
		}
        
		return grad;
	}
    

	
    cflimage  cflimage::xgradient(char cType)
	{
		assert(d_v != NULL);
        
		cflimage grad(d_w, d_h, d_c);
		
		for (int ii=0; ii < d_c; ii++)
		{
            fiComputeImageGradient(&d_v[ii * d_wh], &grad.d_v[ii * d_wh], NULL, NULL, NULL, d_w, d_h,cType);
        }
        
		return grad;
	}
    
    
    
    cflimage  cflimage::ygradient(char cType)
	{
		assert(d_v != NULL);
        
		cflimage grad(d_w, d_h, d_c);
		
		for (int ii=0; ii < d_c; ii++)
		{
            fiComputeImageGradient(&d_v[ii * d_wh], NULL, &grad.d_v[ii * d_wh], NULL, NULL, d_w, d_h,cType);
		}
        
		return grad;
	}
    
    
    
    
	cflimage  cflimage::gradient(cflimage &orientation, char cType)
	{
		assert(d_v != NULL);
		
		cflimage grad(d_w, d_h, d_c);
		if (!orientation.isSameSize(grad)) {orientation.erase(); orientation.create(d_w, d_h, d_c);}
		
		for (int ii=0; ii < d_c; ii++)
		{
			fiComputeImageGradient(&d_v[ii * d_wh], &grad.d_v[ii * d_wh], &orientation.d_v[ii * d_wh], d_w, d_h,cType);
		}
		
		return grad;
	}
	
	
    
	cflimage  cflimage::gradient(cflimage &xgrad, cflimage &ygrad, cflimage &orientation, char cType)
	{
		assert(d_v != NULL);
		
		cflimage grad(d_w, d_h, d_c);
		if (!orientation.isSameSize(grad)) {orientation.erase(); orientation.create(d_w, d_h, d_c);}
		if (!xgrad.isSameSize(grad)) {xgrad.erase(); xgrad.create(d_w, d_h, d_c);}
		if (!ygrad.isSameSize(grad)) {ygrad.erase(); ygrad.create(d_w, d_h, d_c);}
		
		for (int ii=0; ii < d_c; ii++)
		{
			fiComputeImageGradient(&d_v[ii * d_wh], &xgrad.d_v[ii * d_wh], &ygrad.d_v[ii * d_wh], &grad.d_v[ii * d_wh], &orientation.d_v[ii * d_wh], d_w, d_h,cType);
		}
		
		return grad;
	}
    
	
	
	
    
	//////////////////////////////////////////////////////////////////
	//! Begin Load/Save
	//////////////////////////////////////////////////////////////////
	
	
	
	void cflimage::load(const char *filename)
	{
		
		
		
		// erase current image
		erase();
		
		
		
		// look if it exists
        std::ifstream fileImage(filename);
        if (!fileImage.good())
		{
            std::cout << "... failed to read image " << filename << std::endl;
            exit(-1);
            
        } else fileImage.close();
		
        
        
		// the sure bet
		d_v = iio_read_image_float_split(filename, &d_w, &d_h, &d_c);
		
		if (d_v) 
		{
            if (d_c == 2) d_c = 1;
            if (d_c > 3)  d_c = 3;
            
			d_wh = d_w * d_h;
			d_whc = d_c * d_w * d_h;
            
			return;
		}
		
		
		
		
		std::cout << "... failed to read image " << filename << std::endl;
		exit(-1);
		
		
	}
	
	
	
	void cflimage::save(const char *filename)
	{
		
		if (!d_v)
		{
			std::cout << "... failed to save image " << filename << std::endl;
			exit(-1);
		}
		
		iio_save_image_float_split((char*)filename, d_v, d_w, d_h, d_c);
			
		
	}
	
	
    
    
    
    
	//////////////////////////////////////////////////////////////////
	//! Begin Get Basic Data	
	//////////////////////////////////////////////////////////////////
	
	
	
	
	cflimage::operator  flimage()
	{
		return getGray(); 
	}
	
	
	
	
	
	int  cflimage::isSameSize(cflimage &inIm)
	{
		
		if (d_c != inIm.d_c || d_w != inIm.d_w || d_h != inIm.d_h) return 0;
		else return 1;
		
	}
	
	
	
	
	
	
	flimage cflimage::getChannel(int i)
	{
		
		assert(i < d_c);
		
		flimage image(d_w,d_h);
		
		for (int jj=0; jj < d_wh; jj++) image.d_v[jj] = d_v[ i * d_wh + jj];
		
		return image;
	}
	
	
	
	
	
	
	
	//////////////////////////////////////////////////////////////////
	//! Begin Math	
	//////////////////////////////////////////////////////////////////
	
	
	
	cflimage& cflimage::operator= (float a)
	{	
		if (d_v)  for (int j=0; j<d_whc ; j++) d_v[j] = a;	
		
		return *this;
	}
	
	
	void cflimage::operator-= (float a)
	{	
		if (d_v)  for (int j=0; j<d_whc ; j++) d_v[j] -= a;
		
	}
	
	
	void cflimage::operator+= (float a)
	{	
		if (d_v)  for (int j=0; j<d_whc ; j++) d_v[j] += a;
		
	}
	
	
	void cflimage::operator*= (float a)
	{	
		if (d_v)  for (int j=0; j<d_whc ; j++) d_v[j] *= a;
		
	}
	
	
	
	
	float  cflimage::max ()
	{		
		
		assert(d_v != NULL);
		
		float fmax = d_v[0];
		for (int j=0; j < d_whc ; j++)  if (d_v[j] > fmax)  fmax = d_v[j];
		return fmax;
		
		
	}
	
	
	
	float  cflimage::min ()
	{		
		assert(d_v != NULL);
		
		float fmin = d_v[0];
		for (int j=0; j < d_whc ; j++)  if (d_v[j] < fmin)  fmin = d_v[j];
		return fmin;
	}
	
	
	
	float  cflimage::min_channel (int i)
	{		
		assert(d_v != NULL);
		assert(i>= 0 && i < d_c);
		
		float *ptr = &d_v[i * d_wh];
		float fmin = *ptr;
		
		for (int j=0; j < d_wh ; j++, ptr++)  if (*ptr < fmin)  fmin = *ptr;
		return fmin;
	}
	
	
	
	float  cflimage::max_channel (int i)
	{		
		assert(d_v != NULL);
		assert(i>= 0 && i < d_c);
		
		float *ptr = &d_v[i * d_wh];
		float fmax = *ptr;
		
		for (int j=0; j < d_wh ; j++, ptr++)  if (*ptr > fmax)  fmax = *ptr;
		return fmax;
	}
	
	
	
	void  cflimage::min (float m)
	{		
		assert(d_v != NULL);
		for (int j=0; j < d_whc ; j++)  if (d_v[j] > m)  d_v[j] = m;
	}
	
	
	
	void  cflimage::max (float M)
	{		
		assert(d_v != NULL);
		for (int j=0; j < d_whc ; j++)  if (d_v[j] < M)  d_v[j] = M;
	}
	
	
	
	
	void cflimage::thre(float m, float M)
	{
		
		assert(d_v != NULL);	
		for (int ii=0; ii < d_whc; ii++)
		{
			if (d_v[ii] >= M) 	d_v[ii]= M;
			else if (d_v[ii] <= m)  d_v[ii]= m;			
		}
		
	}
	
	
	
	void cflimage::normalizeL1()
	{
		
		float fSum = 0.0f;
		for (int j=0; j < d_whc ; j++)    fSum += d_v[j];
		
		assert(fSum != 0.0f);
		float dfSum = 1.0 / fSum;
		for (int j=0; j < d_whc ; j++)    d_v[j] *= dfSum;
		
	}
	
	
	
	void cflimage::rint()
	{
		
		assert(d_v != NULL);	
		for (int ii=0; ii < d_whc; ii++)		d_v[ii]= rintf(d_v[ii]);			
		
	}
    
	
	void cflimage::abs()
	{
		
		assert(d_v != NULL);	
		for (int ii=0; ii < d_whc; ii++)		d_v[ii]= fabsf(d_v[ii]);			
		
	}
	
	
	
	
	//////////////////////////////////////////////////////////////////
	//! Begin Block Operations	
	//////////////////////////////////////////////////////////////////
	
	
	
	cflimage cflimage::padding(int w, int h, float fValue)
	{
		
		assert(w >= d_w  && h >= d_h);
		
		cflimage image(w,h,d_c);
		image=fValue;
		
		for (int ii=0; ii < d_c; ii++)
		{
			
			for(int j=0; j < d_h; j++)
				for(int i=0; i < d_w; i++)	
					image.d_v[ii * image.d_wh + j * image.d_w + i] = d_v[ ii * d_wh + j * d_w + i];		
			
		}
		
		return image;
	}
	
	
	
	
	
	
	cflimage cflimage::copy(int ipx, int ipy, int iw, int ih)
	{
		
		assert(iw>0 && ih>0);
		assert(ipx>=0 && ipy>=0);
		assert(ipx + iw - 1 < d_w && ipy + ih - 1 < d_h);
		
		cflimage image(iw, ih, d_c);
		
		int nn=0;
		for (int ii=0; ii < d_c; ii++)
		{
			
			int l = ii * d_wh +  ipy * d_w + ipx;
			
			for (int jj = 0; jj < ih; jj++)
			{
				
				for (int kk = 0; kk < iw; kk++,nn++,l++)
					image[nn] = d_v[l];
				
				l += d_w - iw;
			}
			
			
		}
		
		return image;	
	}
	
	
	
	
	void cflimage::paste(const cflimage &im, int ipx, int ipy)
	{
		
		assert(ipx>=0 && ipy>=0);
		assert(ipx + im.d_w - 1 < d_w && ipy + im.d_h - 1 < d_h);
		assert(d_c == im.d_c);
		
		
		for (int ii=0; ii < d_c; ii++)
		{
			
			
			int ll = ii * im.d_wh;
			int nn = ii * d_wh +  ipy * d_w + ipx;
			
			
			for (int jj = 0; jj < im.d_h; jj++)
			{
				
				for (int kk = 0; kk < im.d_w; ll++, nn++, kk++)
					d_v[nn] = im.d_v[ll];
				
				nn += d_w - im.d_w;
			}
			
			
		}
		
	}
	
	
	
	
	
	
	cflimage cflimage::append(cflimage &imIn, int extension)
	{
		
		assert(d_c == imIn.d_c);
		
		//! create image
		cflimage image;
		if (extension == iipHorizontal)
			image.create(d_w + imIn.d_w, MAX(d_h, imIn.d_h), d_c);
		else
			image.create(MAX(d_w, imIn.d_w), d_h + imIn.d_h, d_c);
		
		
		image = 0.0f;
		image.paste(*this, 0, 0);
		if (extension == iipHorizontal)
			image.paste(imIn, d_w, 0);
		else
			image.paste(imIn, 0, d_h);
        
		return image;
	}
	
	
    
    
    
	
	
	///////////////////////////////////////////////
	//! Begin Geometrical transforms
	///////////////////////////////////////////////
	
	cflimage cflimage::mirror(int Orientation)
	{
		
		cflimage image;
		if (Orientation == iipHorizontal)  
		{
			image.create(2 * d_w, d_h, d_c);
			
			for (int ic = 0; ic < d_c; ic++)
			{
				float *iptr = &d_v[ic * d_wh];
				float *optr = &image.d_v[ic * image.d_wh];				
				
				for (int ij = 0; ij < d_h; ij++) 
					for (int ii = 0; ii < d_w; ii++) 
					{
						optr[ ij * 2 * d_w + ii] =  iptr[ ij * d_w + ii];
						optr[ ij * 2 * d_w + 2 * d_w - 1 -ii] =  iptr[ ij * d_w + ii];
					}
			}	
			
			
		}else
		{
			image.create(d_w, 2 * d_h, d_c);
			
			for (int ic = 0; ic < d_c; ic++)
			{
				float *iptr = &d_v[ic * d_wh];
				float *optr = &image.d_v[ic * image.d_wh];	
				
				for (int ij = 0; ij < d_h; ij++) 
					for (int ii = 0; ii < d_w; ii++) 
					{
						optr[ ij * d_w + ii] =  iptr[ ij * d_w + ii];
						optr[ (2 * d_h -1 - ij) * d_w + ii] =  iptr[ ij * d_w + ii];
					}
			}			
			
		}
		
		
		return image;
		
	}
	
	
	
	
	cflimage cflimage::fftUpSample(float fFactorx, float fFactory)
	{
		
		int zd_w = (int) rintf( fFactorx * (float) d_w);
		int zd_h = (int) rintf( fFactory * (float) d_h);
		
		cflimage image(zd_w, zd_h, d_c);
		
		for(int i=0; i < d_c ;i++)
		{
			fiFFTZoom(&d_v[i * d_wh], &image.d_v[i * image.d_wh], fFactorx, fFactory, d_w, d_h);
		}	
		
		return image;
	}
	
    
	
	cflimage cflimage::fftUpSample(float fFactorx)
	{
		cflimage image = (*this).fftUpSample(fFactorx, fFactorx);
		return image;
	}
	
	
	
	
	void cflimage::fftRot(float angle, float xtrans , float ytrans, int flagNoCenter, int flagSymmetric)
	{
		
		angle = -(angle * M_PI / 180.);
		
		// The rotation is decomposed into three shears, two horizontal and one vertical
		(*this).fftShear(  tan(angle * .5),  iipHorizontal, 0.0,   flagNoCenter, flagSymmetric);			
		(*this).fftShear(- sin(angle     ),  iipVertical, -(ytrans),  flagNoCenter, flagSymmetric);		
		(*this).fftShear(  tan(angle * .5),  iipHorizontal, -(xtrans),  flagNoCenter, flagSymmetric);
		
	}
	
	
	
	void cflimage::fftTrans(float xtrans,float ytrans, int flagNoCenter, int flagSymmetric)
	{
		
		if (ytrans != 0.0f)
			(*this).fftShear(0.0, iipVertical,  -ytrans, flagNoCenter, flagSymmetric);			 
		
		if (xtrans != 0.0f)
			(*this).fftShear(0.0, iipHorizontal,  -(xtrans), flagNoCenter, flagSymmetric);			 
		
	}
	
	
	
	
	void cflimage::fftShear(float dAmountOfShear, int iOrientation, float fDelta, int flagNoCenter, int flagSymmetric)
	{
		
		cflimage extended;
		if (!flagSymmetric) 
			extended = *this;
		else
			extended = (*this).mirror(iOrientation);
		
		for (int i=0; i < d_c; i++)
		{
			fiFFTShearPascal((double) dAmountOfShear, &extended.d_v[ i * extended.d_wh],&extended.d_v[i*extended.d_wh], iOrientation,  fDelta,  flagNoCenter,extended.d_w, extended.d_h);
		}
		
		
		if (!flagSymmetric) (*this)=extended;
		else (*this)=extended.copy(0, 0, d_w, d_h);
		
	}
    
    
    
    
    ///////////////////////////////////////////////
	//! Begin Zooming
	///////////////////////////////////////////////
	
    cflimage cflimage::upSampleSplines(float fFactor, int order)
    {
        
        
		int zd_w = (int)rintf( fFactor * (float) d_w);
        int zd_h = (int)rintf( fFactor * (float) d_h);
        
        cflimage image(zd_w, zd_h, d_c);

        

        for(int i=0; i < d_c ;i++)
            apply_zoom(&d_v[i * d_wh], &image.d_v[i * image.d_wh], fFactor, order, d_w, d_h);
        
        
        return image;
		
        
    }
	
    
    
    
    
	
	///////////////////////////////////////////////
	//! Begin Color Conversion
	///////////////////////////////////////////////
	
	
	flimage cflimage::getGray()
	{
		
		flimage image(d_w, d_h);	image=0.0f;
		
		for (int i=0; i < d_whc; i++)
		{
			image.d_v[i % d_wh] += d_v[i];
		}
		
		for (int i=0; i < d_wh; i++)
			image.d_v[i] /= (float) d_c;
		
		return image;
	}
	
	
	
	flimage cflimage::getGray(float wr, float wg, float wb)
	{
		assert(d_c == 1  || d_c == 3);
		
		flimage image(d_w, d_h);	image=0.0f;
        
		if (d_c == 1)  image = (flimage) (*this);
		else
		{
			for (int i=0; i < d_wh; i++) image.d_v[i] = wr * d_v[i] + wg * d_v[d_wh + i] + wb * d_v[2*d_wh+i];
		}
		
		
		return image;
	}
	
	
	
	
	
	
	
    
	cflimage cflimage::binarize(float value, int inverse)
	{
		assert(d_v != NULL);	
		cflimage binary(d_w,d_h,d_c);
		
		fpBinarize(d_v, binary.d_v,  value, inverse, d_whc);
        
		return binary;
	}
	
    
    
    
	void cflimage::Rgb2Yuv(int iflagOrto)
	{
		
		assert(d_c==3);
		cflimage image = *this;
		
		if (iflagOrto)
			fiRgb2YuvO(image.d_v, &image.d_v[d_wh], &image.d_v[2*d_wh], d_v, &d_v[d_wh], &d_v[2*d_wh], d_w, d_h);
		else
			fiRgb2Yuv( image.d_v, &image.d_v[d_wh], &image.d_v[2*d_wh], d_v, &d_v[d_wh], &d_v[2*d_wh], d_w, d_h);
		
	}
	
	
	
	void  cflimage::Yuv2Rgb(int iflagOrto)
	{
		
		assert(d_c==3);
		cflimage image = *this;
		
		
		if (iflagOrto)
			fiYuvO2Rgb(d_v, &d_v[d_wh], &d_v[2*d_wh], image.d_v, &image.d_v[d_wh], &image.d_v[2*d_wh], d_w, d_h);
		else
			fiYuv2Rgb(d_v, &d_v[d_wh], &d_v[2*d_wh], image.d_v, &image.d_v[d_wh], &image.d_v[2*d_wh], d_w, d_h);
	}
	
    
    
    
	
	
    
	///////////////////////////////////////////////
	//! Begin Sampling and Convolution
	///////////////////////////////////////////////
	
	
    
    
    
	
	cflimage cflimage::subSample(int fFactor, int flagCenter)
	{
		
		assert(d_v != NULL);
		int sd_w = (int) floor((float) d_w / (float) fFactor);
		int sd_h = (int) floor((float) d_h / (float) fFactor);
        
		cflimage image(sd_w, sd_h, d_c);
		
		
		for (int ii=0; ii < d_c; ii++)
		{
            if (!flagCenter)
                fiImageSample(&d_v[ii*d_wh], &image.d_v[ii* image.d_wh], fFactor, d_w, d_h);
            else
                fiImageSampleCenter(&d_v[ii*d_wh], &image.d_v[ii* image.d_wh], fFactor, d_w, d_h);    
		}
		
		
		return image;
		
	}
	
	
 
	
	
	cflimage cflimage::subSampleAglomeration(int iFactor)
	{
		
		assert(d_v != NULL);
        
        
		int sd_w = (int) floor((float) d_w / (float) iFactor);
		int sd_h = (int) floor((float) d_h / (float) iFactor);
		
        
		cflimage image(sd_w, sd_h, d_c);
		
		
		for (int ii=0; ii < d_c; ii++)
		{
            fiImageSampleAglomeration(&d_v[ii*d_wh], &image.d_v[ii*image.d_wh], iFactor, d_w, d_h);
		}
		
		
		return image;
		
	}
	
	
	
	
	
	cflimage cflimage::subSample(int fFactor, float fSigma, int flagCenter)
	{
		
        cflimage convolved =  (*this).convolveGauss(fSigma);
        
		cflimage image = convolved.subSample(fFactor,  flagCenter);

		return image;	
	}
	
	
    
	
	
	cflimage cflimage::convolveGauss(float fSigma)
	{
		
		assert(d_v != NULL);
		cflimage image= *this;
		
		
		int ksize;	
		float *kernel;
		kernel = fiFloatGaussKernel(fSigma,ksize);
		
		
		for (int i=0; i < d_c; i++)
		{
			fiFloatHorizontalConvolution( &image.d_v[i*d_wh], &image.d_v[i*d_wh], d_w, d_h, kernel, ksize, 1);
			fiFloatVerticalConvolution( &image.d_v[i*d_wh], &image.d_v[i*d_wh], d_w, d_h, kernel,  ksize, 1);
		}
		
		
		delete[] kernel;
		
		return image;
		
	}
	
	
	
	
	cflimage cflimage::convolve(flimage &kernel)
	{
		assert(d_v != NULL);
		cflimage image(d_w, d_h, d_c);
		
		for (int i=0; i < d_c; i++)
		{
			fiConvol( &d_v[i*d_wh], &image.d_v[i*d_wh], d_w, d_h, kernel.d_v, kernel.d_w, kernel.d_h);
			
		}
		
		return image;	
	}
	
	
	
    
    
    
    cflimage  cflimage::patchMean(float fRadius)
    {
        cflimage image(d_w, d_h, d_c);
		image = *this;
        
        for (int i=0; i < d_c; i++)
		{
			fiPatchMean( &d_v[i*d_wh], &image.d_v[i*d_wh], fRadius, d_w, d_h);
			
		}
        
        return image;
    }
    
    
    cflimage  cflimage::patchVar(float fRadius)
    {
        cflimage image(d_w, d_h, d_c);
		image = *this;
        
        for (int i=0; i < d_c; i++)
		{
			fiPatchVar( &d_v[i*d_wh], &image.d_v[i*d_wh], fRadius, d_w, d_h);
			
		}
        
        return image;
    }
    
    
    cflimage  cflimage::patchMin(float fRadius)
    {
        cflimage image(d_w, d_h, d_c);
		image = *this;
        
        for (int i=0; i < d_c; i++)
		{
			fiPatchMin( &d_v[i*d_wh], &image.d_v[i*d_wh], fRadius, d_w, d_h);
			
		}
        
        return image;
        
    }
    
    
    
    cflimage  cflimage::patchMax(float fRadius)
    {
        cflimage image(d_w, d_h, d_c);
		image = *this;
        
        for (int i=0; i < d_c; i++)
		{
			fiPatchMax( &d_v[i*d_wh], &image.d_v[i*d_wh], fRadius, d_w, d_h);
			
		}
        
        return image;
    }
    
    
    
    cflimage  cflimage::patchMedian(float fRadius)
    {
        cflimage image(d_w, d_h, d_c);
        image = *this;
        
        for (int i=0; i < d_c; i++)
        {
            fiPatchMedian( &d_v[i*d_wh], &image.d_v[i*d_wh], fRadius, d_w, d_h);
            
        }
        return image;
    }
    
    
    
    
    
    
    
    
	cflimage  cflimage::patchMean(flimage &kernel)
	{
		
		cflimage image(d_w, d_h, d_c);
		image = *this;
		
		int spd_w = kernel.d_w / 2;
		int spd_h = kernel.d_h / 2;
		int boundary =  MAX(kernel.d_h, kernel.d_w) + 1;
		
		
		for (int ipx = 0; ipx < d_w - boundary; ipx++) 
			for (int ipy = 0; ipy < d_h - boundary ; ipy++) 
			{
				
				for (int iC=0; iC < d_c; iC++)
				{
					
					float fMean = 0.0f; 
					float *ptr = &d_v[ iC * d_wh + ipy *  d_w + ipx];
					float *ptrk = kernel.d_v;			
					
					for (int s = 0; s < kernel.d_h; s++)
					{
						
						for(int r = 0 ; r < kernel.d_w; r++, ptr++, ptrk++)
						{
							fMean += *ptrk * (*ptr); 
						}
						
						ptr += d_w - kernel.d_w;
					}
					
					
					image[iC * d_wh + (ipy + spd_h) * d_w + ipx + spd_w] = fMean;
					
				}
				
			}
		
		return image;
	}
	
	
	
	
    // ATTENTION: THIS IS A LIST-KERNEL. PLACEHOLDER. 
	cflimage  cflimage::patchListMean(flimage &kernel)
	{
		
		cflimage image(d_w, d_h, d_c);
		image = *this;

      // this kernel is represented as a list of offset & values : offx, offy, val.
		int k_len = kernel.list_len;
		assert(kernel.list_len&& kernel.offx && kernel.offy && kernel.offval); 
		
      // the kernel should also be stored as an image so we can acces its support
		int spd_w = kernel.d_w / 2;
		int spd_h = kernel.d_h / 2;
		int boundary =  MAX(kernel.d_h, kernel.d_w) + 1;
		
		for (int ipx = 0; ipx < d_w - boundary; ipx++) 
			for (int ipy = 0; ipy < d_h - boundary ; ipy++) 
			{
				
				for (int iC=0; iC < d_c; iC++)
				{

               int   *ptrk_offx   = kernel.offx;
               int   *ptrk_offy   = kernel.offy;
               float *ptrk_val    = kernel.offval;

					float fMean = 0.0f; 
					float *ptr1 = d_v + iC * d_wh;

               for (int jj = 0; jj < k_len; jj++)
               {
                  float *ptr = ptr1 + ((ipy + *ptrk_offy) * d_w + ipx + *ptrk_offx);

						fMean += *ptrk_val * (*ptr); 

                  ptrk_offx++;
                  ptrk_offy++;
                  ptrk_val++;

               }

					image[iC * d_wh + (ipy + spd_h) * d_w + ipx + spd_w] = fMean;
					
				}
				
			}
		
		return image;
		
	}
	
	
	
	
	
	
	
	
	
	
	cflimage  cflimage::patchVar(flimage &kernel)
	{
		
		cflimage image(d_w, d_h, d_c);
		image = 0.0f;
		
		
		int spd_w = (kernel.d_w - 1) / 2;
		int spd_h = (kernel.d_h - 1) / 2;
		int boundary = MAX(spd_w, spd_h) + 1;
		
		
		for (int ipx = boundary; ipx < d_w - boundary; ipx++) 
			for (int ipy = boundary; ipy < d_h - boundary ; ipy++) 
			{
				
				for (int iC=0; iC < d_c; iC++)
				{
					
					float fMean = 0.0f; 
					float fMean2 = 0.0f;
					
					float *ptr = &d_v[ iC * d_wh + (ipy - spd_h ) *  d_w + (ipx - spd_w)];
					float *ptrk = kernel.d_v;			
					
					for (int s = 0; s < kernel.d_h; s++)
					{
						
						for(int r = 0 ; r < kernel.d_w; r++, ptr++, ptrk++)
						{
							fMean += *ptrk * (*ptr); 
							fMean2 += *ptrk * (*ptr) * (*ptr); 
						}
						
						
						ptr += d_w - kernel.d_w;
						
					}
					
					image[iC * d_wh + ipy * d_w + ipx] = fMean2 - fMean*fMean;
					
				}
				
			}
		
		return image;
	}
	

    
    
    ///////////////////////////////////////////////
	//! Begin Patch Image Distances
	///////////////////////////////////////////////

    float distanceL2(cflimage &input1, cflimage &input2, int ipx, int ipy, int iqx, int iqy)
	{
		
		assert(ipx >=0 && ipy>=0 && iqx>=0 && iqy>=0);
		assert(ipx < input1.w() && ipy  < input1.h() && iqx  < input2.w() && iqy  < input2.h() );
		
        float fDif = 0.0f;
		float fDist = 0.0f;
		for (int ii=0; ii < input1.d_c; ii++)
		{
			
			float *ptr1 = &input1.d_v[ii * input1.d_wh +  ipy  * input1.d_w + ipx];
			float *ptr2 = &input2.d_v[ii * input2.d_wh +  iqy  * input2.d_w + iqx];
            
            fDif = *ptr1 - *ptr2;
            fDist += fDif * fDif;
			
		}
		
		return fDist;
	}
	
    
    
    
	float distancePatchL2(cflimage &input1, cflimage &input2, int ipx, int ipy, int iqx, int iqy, int r_w, int r_h)
	{
		
		assert(ipx >=0 && ipy>=0 && iqx>=0 && iqy>=0);
		assert(ipx + r_w -1 < input1.w() && ipy + r_h - 1 < input1.h() && iqx + r_w - 1 < input2.w() && iqy + r_h - 1 < input2.h() );
		
		
		
		float fDist = 0.0f;
		float dif = 0.0f;
		float fSum = 1.0 / ((float) input1.d_c * (float) r_w * (float) r_h);
        
		for (int ii=0; ii < input1.d_c; ii++)
		{
			
			float *ptr1 = &input1.d_v[ii * input1.d_wh +  ipy  * input1.d_w + ipx];
			float *ptr2 = &input2.d_v[ii * input2.d_wh +  iqy  * input2.d_w + iqx];
			
			for (int jj = 0; jj < r_h; jj++)
			{
				
				for (int kk = 0; kk < r_w; kk++,ptr1++, ptr2++)
				{
					dif = *ptr1 - *ptr2;
					fDist += dif * dif;
				}
				
				
				ptr1 += input1.d_w - r_w;
				ptr2 += input2.d_w - r_w;
				
			}
			
			
			
		}
		
		
		return fSum * fDist;
	}
	
	
	
	
	
	float distancePatchWL2(cflimage &input1, cflimage &input2, int ipx, int ipy, int iqx, int iqy, flimage &kernel)
	{
		
		int pd_w = kernel.d_w ;
		int pd_h = kernel.d_h ;
		
		assert(ipx >=0 && ipy>=0 && iqx>=0 && iqy>=0);
		assert(ipx + pd_w < input1.w() && ipy + pd_h < input1.h() && iqx + pd_w < input2.w() && iqy + pd_h < input2.h() );
		
		
		
		float fDist = 0.0f;
		float dif = 0.0f;
        
        
		for (int ii=0; ii < input1.d_c; ii++)
		{
			
			float *ptr1 = &input1.d_v[ii * input1.d_wh +  ipy  * input1.d_w + ipx ];
			float *ptr2 = &input2.d_v[ii * input2.d_wh +  iqy  * input2.d_w + iqx ];
			
			float *ptrk = kernel.d_v;
			
			for (int jj = 0; jj < kernel.d_h; jj++)
			{
				
				for (int kk = 0; kk < kernel.d_w; kk++,ptr1++, ptr2++, ptrk++)
				{
					dif = *ptr1 - *ptr2;
					fDist += *ptrk * dif * dif;
				}
				
				
				ptr1 += input1.d_w - kernel.d_w;
				ptr2 += input2.d_w - kernel.d_w;
				
			}
			
			
			
		}
		
		
		return  fDist / (float) input1.d_c;
	}
	

	
	// ATTENTION: there is no check if the patch goes outside the domain of the image
	float distancePatchListWL2(cflimage &input1, cflimage &input2, int ipx, int ipy, int iqx, int iqy, flimage &kernel)
	{
      // this kernel is represented as a list of offset & values : offx, offy, val.
		int k_len = kernel.list_len;
		assert(kernel.list_len&& kernel.offx && kernel.offy && kernel.offval); 
		
      // the kernel should also be stored as an image so we can acces its support
		int pd_w = kernel.d_w ;
		int pd_h = kernel.d_h ;
		int sph = pd_h / 2 ;
		int spw = pd_w / 2 ;
		assert(ipx >=0 && ipy>=0 && iqx>=0 && iqy>=0);
		assert(ipx + pd_w < input1.w() && ipy + pd_h < input1.h() && 
            iqx + pd_w < input2.w() && iqy + pd_h < input2.h() );

	   // short names	
		int in1_w = input1.w();
		int in2_w = input2.w();

		
		float fDist = 0.0f;
		float dif = 0.0f;
        
		for (int ii=0; ii < input1.d_c; ii++)
		{
			float *ptr1 = &input1.d_v[ii * input1.d_wh +  ipy  * in1_w + ipx ];
			float *ptr2 = &input2.d_v[ii * input2.d_wh +  iqy  * in2_w + iqx ];
		
			int   *ptrk_offx   = kernel.offx;
			int   *ptrk_offy   = kernel.offy;
			float *ptrk_val    = kernel.offval;
			
			for (int jj = 0; jj < k_len; jj++)
			{
				float *tptr1 = ptr1 + ((*ptrk_offy) * in1_w + *ptrk_offx);
				float *tptr2 = ptr2 + ((*ptrk_offy) * in2_w + *ptrk_offx);

				dif = *tptr1 - *tptr2;
				fDist += *ptrk_val * dif * dif;
				
				ptrk_offx++;
				ptrk_offy++;
				ptrk_val++;
				
			}
						
		}
		
		return  fDist / (float) input1.d_c;
	}
	
	
	
	
	float distancePatchL2M(cflimage &input1, cflimage &input2, int ipx, int ipy, int iqx, int iqy, int r_w, int r_h, cflimage &mean1, cflimage &mean2)
	{
		
		assert(ipx >=0 && ipy>=0 && iqx>=0 && iqy>=0);
		assert(ipx + r_w < input1.w() && ipy + r_h < input1.h() && iqx + r_w < input2.w() && iqy + r_h < input2.h() );
		
		int sph = r_h / 2;
		int spw = r_w / 2;
		
		float fDist = 0.0f;
		float dif = 0.0f;
        float fSum = 1.0 / ((float) input1.d_c * (float) r_w * (float) r_h);
        
        
		for (int ii=0; ii < input1.d_c; ii++)
		{
			
			float *ptr1 = &input1.d_v[ii * input1.d_wh +  ipy  * input1.d_w + ipx];
			float *ptr2 = &input2.d_v[ii * input2.d_wh +  iqy  * input2.d_w + iqx];
			
			float fMean1 = mean1.d_v[ii * input1.d_wh +  (ipy + sph) * input1.d_w + ipx + spw];
			float fMean2 = mean2.d_v[ii * input2.d_wh +  (iqy + sph) * input2.d_w + iqx + spw];
			float fMean = -fMean1 + fMean2;
			
			
			
			
			for (int jj = 0; jj < r_h; jj++)
			{
				
				for (int kk = 0; kk < r_w; kk++,ptr1++, ptr2++)
				{
					dif = *ptr1 - *ptr2 + fMean;
					fDist += dif * dif;
				}
				
				
				ptr1 += input1.d_w - r_w;
				ptr2 += input2.d_w - r_w;
				
			}
			
			
			
		}
		
		
		return fSum * fDist;
	}
	
	
	
	
	
	float distancePatchWL2M(cflimage &input1, cflimage &input2, int ipx, int ipy, int iqx, int iqy, flimage &kernel,  cflimage &mean1, cflimage &mean2)
	{
		
		int pd_w = kernel.d_w ;
		int pd_h = kernel.d_h ;
		
		int sph = pd_h / 2 ;
		int spw = pd_w / 2 ;
		
        
        if (ipx < 0 || ipy < 0 || iqx < 0 || iqy < 0) printf("%d %d %d %d\n", ipx,ipy,iqx,iqy);
        
		assert(ipx>=0 && ipy>=0 && iqx>=0 && iqy>=0);
		assert(ipx + pd_w < input1.w() && ipy + pd_h < input1.h() && iqx + pd_w < input2.w() && iqy + pd_h < input2.h() );
		
		
		
		float fDist = 0.0f;
		float dif = 0.0f;
		
		for (int ii=0; ii < input1.d_c; ii++)
		{
			
			float *ptr1 = &input1.d_v[ii * input1.d_wh +  ipy  * input1.d_w + ipx ];
			float *ptr2 = &input2.d_v[ii * input2.d_wh +  iqy  * input2.d_w + iqx ];
			
			float fMean1 = mean1.d_v[ii * input1.d_wh +  (ipy + sph) * input1.d_w + ipx + spw];
			float fMean2 = mean2.d_v[ii * input2.d_wh +  (iqy + sph) * input2.d_w + iqx + spw];
			float fMean = -fMean1 + fMean2;
			
			
			
			float *ptrk = kernel.d_v;
			
			for (int jj = 0; jj < kernel.d_h; jj++)
			{
				
				for (int kk = 0; kk < kernel.d_w; kk++,ptr1++, ptr2++, ptrk++)
				{
					dif = *ptr1 - *ptr2 + fMean;
					fDist += *ptrk * dif * dif;
				}
				
				
				ptr1 += input1.d_w - kernel.d_w;
				ptr2 += input2.d_w - kernel.d_w;
				
			}
			
			
			
		}
		
		
		return  fDist / (float) input1.d_c;

	}
	
	
	
	
	
	// ATTENTION: there is no check if the patch goes outside the domain of the image
	// ATTENTION: the patch mean must be pre-computed using the window sph spw are not good here
	float distancePatchListWL2M(cflimage &input1, cflimage &input2, int ipx, int ipy, int iqx, int iqy, flimage &kernel,  cflimage &mean1, cflimage &mean2)
	{
      // this kernel is represented as a list of offset & values : offx, offy, val.
		int k_len = kernel.list_len;
		assert(kernel.list_len&& kernel.offx && kernel.offy && kernel.offval); 

      // the kernel should also be stored as an image so we can acces its support
		int pd_w = kernel.d_w ;
		int pd_h = kernel.d_h ;
		int sph = pd_h / 2 ;
		int spw = pd_w / 2 ;
		assert(ipx >=0 && ipy>=0 && iqx>=0 && iqy>=0);
		assert(ipx + pd_w < input1.w() && ipy + pd_h < input1.h() && 
            iqx + pd_w < input2.w() && iqy + pd_h < input2.h() );

	   // short names	
		int in1_w = input1.w();
		int in2_w = input2.w();
		
		
		float fDist = 0.0f;
		float dif = 0.0f;
        
        
		for (int ii=0; ii < input1.d_c; ii++)
		{
			float *ptr1 = &input1.d_v[ii * input1.d_wh +  ipy  * in1_w + ipx ];
			float *ptr2 = &input2.d_v[ii * input2.d_wh +  iqy  * in2_w + iqx ];
			
			float fMean1 = mean1.d_v[ii * input1.d_wh +  (ipy + sph)  * in1_w + ipx + spw ];
			float fMean2 = mean2.d_v[ii * input2.d_wh +  (iqy + sph)  * in2_w + iqx + spw ];
			float fMean = -fMean1 + fMean2;

			
			int   *ptrk_offx   = kernel.offx;
			int   *ptrk_offy   = kernel.offy;
			float *ptrk_val    = kernel.offval;
			
			for (int jj = 0; jj < k_len; jj++)
			{
				float *tptr1 = ptr1 + ((*ptrk_offy) * in1_w + *ptrk_offx);
				float *tptr2 = ptr2 + ((*ptrk_offy) * in2_w + *ptrk_offx);
				
				dif = *tptr1 - *tptr2 + fMean;
				fDist += *ptrk_val * dif * dif;
				
				ptrk_offx++;
				ptrk_offy++;
				ptrk_val++;
				
			}
			
		}
		
		return  fDist / (float) input1.d_c;
	}

	
	
	
	
	
	float distancePatchL1(cflimage &input1, cflimage &input2, int ipx, int ipy, int iqx, int iqy, int r_w, int r_h)
	{
		
		assert(ipx >=0 && ipy>=0 && iqx>=0 && iqy>=0);
		assert(ipx + r_w < input1.w() && ipy + r_h < input1.h() && iqx + r_w < input2.w() && iqy + r_h < input2.h() );
		
		
		
		float fDist = 0.0f;
		float dif = 0.0f;
		
		for (int ii=0; ii < input1.d_c; ii++)
		{
			
			float *ptr1 = &input1.d_v[ii * input1.d_wh +  ipy  * input1.d_w + ipx];
			float *ptr2 = &input2.d_v[ii * input2.d_wh +  iqy  * input2.d_w + iqx];
			
			for (int jj = 0; jj < r_h; jj++)
			{
				
				for (int kk = 0; kk < r_w; kk++,ptr1++, ptr2++)
				{
					dif = *ptr1 - *ptr2;
					fDist += fabsf(dif);
				}
				
				
				ptr1 += input1.d_w - r_w;
				ptr2 += input2.d_w - r_w;
				
			}
			
			
			
		}
		
		
		return fDist;
	}
	
	
	
	
	
	float distancePatchWL1(cflimage &input1, cflimage &input2, int ipx, int ipy, int iqx, int iqy, flimage &kernel)
	{
		
		int pd_w = kernel.d_w ;
		int pd_h = kernel.d_h ;
		
		assert(ipx >=0 && ipy>=0 && iqx>=0 && iqy>=0);
		assert(ipx + pd_w < input1.w() && ipy + pd_h < input1.h() && iqx + pd_w < input2.w() && iqy + pd_h < input2.h() );
		
		
		
		float fDist = 0.0f;
		float dif = 0.0f;
		
		for (int ii=0; ii < input1.d_c; ii++)
		{
			
			float *ptr1 = &input1.d_v[ii * input1.d_wh +  ipy  * input1.d_w + ipx ];
			float *ptr2 = &input2.d_v[ii * input2.d_wh +  iqy  * input2.d_w + iqx ];
			
			float *ptrk = kernel.d_v;
			
			for (int jj = 0; jj < kernel.d_h; jj++)
			{
				
				for (int kk = 0; kk < kernel.d_w; kk++,ptr1++, ptr2++, ptrk++)
				{
					dif = *ptr1 - *ptr2;
					fDist += *ptrk * fabsf(dif);
				}
				
				
				ptr1 += input1.d_w - kernel.d_w;
				ptr2 += input2.d_w - kernel.d_w;
				
			}
			
			
			
		}
		
		
		return fDist;
	}
	
	
	
	
	
	///////////////////////////////////////////////
	//! End Patch Image Distances
	///////////////////////////////////////////////
	
	
	

    
    
    
    
    
    
    
    
    
	/// Class flimage
	
	flimage::flimage() : cflimage()
	{
	}	
	
	
	
	
	flimage::flimage(int w, int h) : cflimage(w, h, 1)
	{
	}	
	
	
	flimage::flimage(int w, int h, float *ptr) : cflimage(w, h, 1)
	{
		memcpy(this->d_v, ptr, w * h * sizeof(float));
	}	
	
	
	
	flimage::flimage(const flimage& im)
	: cflimage(im)
	{
      this->list_len=im.list_len;
      this->offx=im.offx;
      this->offy=im.offy;
      this->offval=im.offval;
  }
	
	
	
	
	flimage& flimage::operator=(const flimage & im)
	{
		cflimage::operator=(im);

      this->list_len=im.list_len;
      this->offx=im.offx;
      this->offy=im.offy;
      this->offval=im.offval;
		return *this;
	}
	
	
	
	
	void flimage::create(int w, int h)
	{	
		cflimage::create(w,h,1);
	}
    
    
    
    
    
    
    cflmovie::cflmovie():
    d_n(0), filename(NULL), fileadress(NULL), writable(0), strings(NULL), pos(-1), outfile(0)
    {}
    
    
    
    
    cflmovie::cflmovie(const char* ifilename,  int inframes) 
    :  d_n(inframes), filename(NULL), fileadress(NULL), writable(0), strings(NULL), pos(-1), outfile(0)
    {
        
        filename = new char[128];
        strcpy(filename,ifilename); 
        
        
        outfile.open(filename);
        if (!outfile.is_open())
        {
            printf("cflmovie file not writtable\n");
            exit(-1); 
        }
        
        outfile << "cflmovie" << std::endl;
        outfile << d_n << std::endl;
        
        
        writable = true;
        pos = -1;
        
    }
    
    
    
    cflmovie::cflmovie(const cflmovie & imovie)
	:  d_n(imovie.d_n), filename(NULL), fileadress(NULL), writable(imovie.writable), strings(NULL), pos(imovie.pos), outfile(0)
    {
		
        strcpy(filename, imovie.filename);
        strcpy(fileadress, imovie.fileadress);
        
        if (imovie.strings)
        {
            strings = new char*[d_n];
            for (int i=0; i < d_n; i++)
            {
                strcpy(strings[i], imovie.strings[i]);
            }
            
        }
        
        
    }
    
    
    
    cflmovie::cflmovie(const char * ifilename) 
	:  d_n(0), filename(NULL), fileadress(NULL), writable(0), strings(NULL), pos(-1), outfile(0) 
    {
        
        /// Reading 
        writable = false; pos = -1;
        fileadress = new char[128];
        filename = new char[128];
        
        strcpy(filename,ifilename); 
        
        
        
        std::ifstream file(ifilename);
        if (!file.is_open())
        {
            printf("cflmovie file not found or impossible to open\n");
            exit(-1); 
        }
        
        /// Reading format
        char  iinfo[128];
        file >> iinfo;
        
        if( strcmp(iinfo, "cflmovie") == 0)
        {
            
			file >> d_n;
			
			
			strings = new char*[d_n];
			for(int i=0; i < d_n; i++) strings[i] = new char[128];
            
			for(int i = 0; i < d_n ; i++ )
			{
				file >>iinfo;
				strcpy(strings[i],iinfo);
				
			}
			
            
			
        }else{
			
			printf("Error: not a correct cflmovie list file\n");
			exit(-1);
			
        }
		
    }
    
    
    void cflmovie::load(const char * ifilename)
    {
        
        /// Reading 
        writable = false; pos = -1;
        fileadress = new char[128];
        filename = new char[128];
        
        strcpy(filename,ifilename); 
        
        
        
        std::ifstream file(ifilename);
        if (!file.is_open())
        {
            printf("cflmovie file not found or impossible to open\n");
            exit(-1); 
        }
        
        /// Reading format
        char  iinfo[128];
        file >> iinfo;
        
        if( strcmp(iinfo, "cflmovie") == 0)
        {
            
            file >> d_n;
            
            
            
            strings = new char*[d_n];
            for(int i=0; i < d_n; i++) strings[i] = new char[128];
            
            for(int i = 0; i < d_n ; i++ )
            {
                file >>iinfo;
                strcpy(strings[i],iinfo);
                
            }
            
            
            
        }else{
            
            printf("Error: not a correct cflmovie list file\n");
            exit(-1);
            
        }
        
    }
    
    
    
    
    
    cflmovie& cflmovie::operator= (const cflmovie& im)
    {
        printf("warning :: using cflmovie operator = which is not defined properly\n");
        
        if (&im == this)
        { 
            return *this;
        }	
        
        return *this;	
    }
    
    
    
    
    cflmovie::~cflmovie()
    {
        
        if (outfile.is_open()) outfile.close();	
        
    }
    
    
    
    
    
    cflimage cflmovie::getframe(int fpos)
    {
        
        /*	
         if (strcmp(fileadress,"") != 0)
         {
         
         char *imname = new char[128];
         strcpy(imname, fileadress);
         strcat(imname, "/");
         strcat(imname, strings[fpos]);
         
         cflimage image;
         image.load(imname);
         
         if (DEBUG) printf("......END get frame %d\n", fpos);
         return image;
         
         } else
         */
        //{
		
		cflimage image;
		
		
		image.load(strings[fpos]);
		return image;
		
        //}		
        
        
    }
    
    
    
    
    
    void cflmovie::write(cflimage &frame)
    {
        
        if (writable)
        {
            
            pos++;
            
            char* imfilename = new char[128];
            strcpy(imfilename,filename); 
            
            char buf[128];
            sprintf(buf, "%d", pos);
            
            strcat(imfilename, "_");
            strcat(imfilename, buf);
            strcat(imfilename, ".png");
            
			
            frame.save(imfilename);
            
            outfile <<  imfilename << std::endl;
            
        }
        
    }
    
    
    
    
}    






void  printusage(char *pname,
				 char *gp,
				 std::vector<OptStruct*>  &opt,
				 std::vector<ParStruct*>  &par)
{
	
	//int nopt = opt.size();
	int npar = par.size();
	
	///// USAGE
	printf("\nusage: %s ", pname);
	for(int i=0; i < (int) strlen(gp); i++)
		if (gp[i] != ':')
		{
			printf("[-%c",gp[i]);
			
			if (i+1 < (int) strlen(gp) && gp[i+1] ==  ':') printf(" %c] ", gp[i]);
			else printf("] ");
			
		}
	
	for(int i=0; i < npar; i++)
		printf(" %s ", par[i]->name);
	
	printf("\n");
	//// PARAMETERS
	
	int j=0;
	for(int i=0; i < (int) strlen(gp); i++)
		if (gp[i] != ':')
		{
			printf("\t-%c",gp[i]);
			
			if (i+1 < (int) strlen(gp) && gp[i+1] ==  ':') {
				
				printf("  %c\t %s ", gp[i], opt[j]->comment);
				if (opt[j]->defvalue != NULL) printf("(Default: %s)",opt[j]->defvalue);
				
				printf("\n");
				
			}		
			else printf("\t %s \n", opt[j]->comment);
			
			j++;
		}
	
	
	for(int i=0; i < npar; i++)
	{
		printf("\t%s",par[i]->name);
		printf("\t %s\n", par[i]->comment);
	}
	
	
}




int parsecmdline(char *pname,
				 char *function,
				 int argc, char **argv,
                 std::vector <OptStruct*> & opt,
				 std::vector <ParStruct*> & par)
{
	int nopt = opt.size();
	int npar = par.size();
	
	
	char *gp = new char[2*nopt+1];
	gp[0]='\0';
	
	
	for(int i=0; i < nopt; i++) { opt[i]->flag = 0; opt[i]->value=NULL; strcat(gp, opt[i]->gp);}
	for(int i=0; i < npar; i++) { par[i]->value = NULL;}
	
	opterr = 0;	// No messages by getopt
	
	int c;
	while ((c = getopt (argc, argv, gp)) != -1)
	{
		
		int j=0;
		for(unsigned int i=0; i < strlen(gp); i++)
			if (c == gp[i])
			{
				
				opt[j]->flag = 1;
				/*				if (optarg != NULL && optarg[0] == '-')
				 {	 
				 printf("\n%s: %s\n", pname, function);
				 fprintf (stderr, "\nerror: option -%c requires an argument.\n", c);
				 printusage(pname,function, gp, argc, argv, opt, nopt, par, npar);
				 return 0;
				 
				 }*/
				
				opt[j]->value = optarg;
				break;
				
			} else if (gp[i] != ':') j++;
		
		
		
		if (c == '?')
		{	
			
			unsigned int i = 0;
			for(i=0; i < strlen(gp); i++)
				if (optopt == gp[i])
				{
					printf("\n%s: %s\n", pname, function);
					fprintf (stderr, "\nerror: option -%c requires an argument.\n", optopt);
					break;	
				}
			
			if (i == strlen(gp)) { 	printf("\n%s: %s\n", pname, function);
				fprintf (stderr, "\nerror: unknown option `-%c'.\n", optopt);
			}
			
			printusage(pname, gp,  opt,  par);
			return 0;
			
		}
		
	}
	
	
	//// Setting default values for non selected options
	for(int j=0; j < nopt; j++)
		if (opt[j]->flag == 0 && opt[j]->defvalue != NULL) opt[j]->value =  opt[j]->defvalue;
	
	
	if (argc - optind != npar) {
		printf("\n%s: %s\n", pname, function);
		fprintf (stderr, "\nerror: incorrect number of parameters\n");
		printusage(pname, gp,  opt,par);
		return 0;
	}
	
	int i=0;
	for (int index = optind; index < argc ; index++, i++){
		par[i]->value = argv[index];
	}
	
	return 1;
	
	
}




