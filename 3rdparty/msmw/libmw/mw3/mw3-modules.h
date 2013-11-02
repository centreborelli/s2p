/**
 * @file mw3-modules.h
 *
 * @author Nicolas Limare (2009)
 *
 * mw3-modules API header
 */

#ifndef _MW3_MODULES_
#define _MW3_MODULES_

#include "mw3-x11.h"

/* src/wave/biowave2.c */
void biowave2(int *NumRec, int *Haar, int *Edge, int *FilterNorm,
              Fimage Image, Wtrans2d Output, Fsignal Ri1, Fsignal Ri2);
/* src/wave/dybiowave2.c */
void dybiowave2(int *NumRec, int *StopDecim, int *Ortho, int *Edge,
                int *FilterNorm, Fimage Image, Wtrans2d Output, Fsignal Ri1,
                Fsignal Ri2);
/* src/wave/ibiowave2.c */
void ibiowave2(int *NumRec, int *Haar, int *Edge, int *FilterNorm,
               Wtrans2d Wtrans, Fimage Output, Fsignal Ri1, Fsignal Ri2);
/* src/wave/ibiowave1.c */
void ibiowave1(int *NumRec, int *Edge, int *FilterNorm, Wtrans1d Wtrans,
               Fsignal Output, Fsignal Ri1, Fsignal Ri2);
/* src/wave/biowave1.c */
void biowave1(int *NumRec, int *Edge, int *FilterNorm, Fsignal Signal,
              Wtrans1d Output, Fsignal Ri1, Fsignal Ri2);
/* src/wave/packets/wp2drecomp.c */
void wp2drecomp(Wpack2d pack, Fimage F);
/* src/wave/packets/wp2dfreqorder.c */
void wp2dfreqorder(int level, Fsignal order, char *inverse_flag);
/* src/wave/packets/wp2dchangepack.c */
void wp2dchangepack(Wpack2d pack_in, Wpack2d pack_out, Cimage tree_out);
/* src/wave/packets/wp2doperate.c */
void wp2doperate(Fimage input, Fsignal Ri, Fsignal Ri_biortho, Cimage tree,
                 Fimage output, float *threshold_hard, float *threshold_soft,
                 int *translation, float *track_noise_hard,
                 float *track_noise_soft, int *track_noise_level,
                 Fimage pfilter, int *convolution_level);
/* src/wave/packets/wp2dview.c */
void wp2dview(Fimage A, Fsignal Ri, Fsignal Ri_biortho, Cimage tree,
              char *do_not_reorder_flag, char *do_not_rescale_flag,
              Fimage toDisplay, Wpack2d input_pack);
/* src/wave/packets/wp2dchecktree.c */
int wp2dchecktree(Cimage tree, char *display_on_flag);
/* src/wave/packets/wp2deigenval.c */
void wp2deigenval(Fimage filter, Cimage tree, Fimage pfilter, float *limit);
/* src/wave/packets/wp2dmktree.c */
void wp2dmktree(int level, Cimage tree, char *sinc_flag, char *quinc_flag,
                char *wave_flag, char *full_flag, char *mirror_flag);
/* src/wave/packets/wp2ddecomp.c */
void wp2ddecomp(Fimage A, Fsignal Ri, Fsignal Ri_biortho, Wpack2d pack,
                Cimage tree);
/* src/wave/packets/wp2dchangetree.c */
void wp2dchangetree(Cimage treeIn, Cimage treeOut, char *up_tree,
                    char *down_tree, int *new_tree_size, char *prune_tree,
                    Cimage tree_for_max, Cimage tree_for_min);
/* src/wave/packets/wpsconvolve.c */
void wpsconvolve(Fsignal Signal, Fsignal Output, Fsignal Ri, char *upSample,
                 char *band, char *oddSize);
/* src/wave/owave2.c */
void owave2(int *NumRec, int *Haar, int *Edge, int *Precond, int *Inverse,
            int *FilterNorm, Fimage Image, Wtrans2d Output, Fsignal Ri,
            Fimage Edge_Ri);
/* src/wave/precond1d.c */
void precond1d(int *Inverse, Fsignal Signal, Fsignal Output, Fimage Edge_Ri);
/* src/wave/iowave1.c */
void iowave1(int *NumRec, int *Haar, int *Edge, int *Precond, int *Inverse,
             int *FilterNorm, Wtrans1d Wtrans, Fsignal Output, Fsignal Ri,
             Fimage Edge_Ri);
/* src/wave/iowave2.c */
void iowave2(int *NumRec, int *Haar, int *Edge, int *Precond, int *Inverse,
             int *FilterNorm, Wtrans2d Wtrans, Fimage Output, Fsignal Ri,
             Fimage Edge_Ri);
/* src/wave/owave1.c */
void owave1(int *NumRec, int *Haar, int *Edge, int *Precond, int *Inverse,
            int *FilterNorm, Fsignal Signal, Wtrans1d Output, Fsignal Ri,
            Fimage Edge_Ri);
/* src/wave/sconvolve.c */
void sconvolve(Fsignal Signal, Fsignal Output, int *DownRate, int *UpRate,
               int *ReflIR, int *Band, int *Edge, int *Prolong, Fsignal Ri,
               Fimage Edge_Ri);
/* src/wave/owtrans_fimage.c */
void owtrans_fimage(Wtrans2d Wtrans, Fimage Output, int *NumRec,
                    double *Contrast);
/* src/wave/ridgelet/ridgrecpol.c */
void ridgrecpol(Fimage in_re, Fimage in_im, Fimage out_re, Fimage out_im);
/* src/wave/ridgelet/istkwave1.c */
void istkwave1(int np, Fsignal in, Fsignal out);
/* src/wave/ridgelet/ridgelet.c */
void ridgelet(Fimage in_re, Fimage in_im, int np, Fimage out_re,
              Fimage out_im);
/* src/wave/ridgelet/iridgelet.c */
void iridgelet(Fimage in_re, Fimage in_im, int np, Fimage out_re,
               Fimage out_im);
/* src/wave/ridgelet/stkwave1.c */
void stkwave1(int np, Fsignal in, Fsignal out);
/* src/wave/ridgelet/ridgthres.c */
void ridgthres(Fimage in_re, Fimage in_im, int np, Fimage out_re,
               Fimage out_im, int pourcent);
/* src/wave/ridgelet/ridgpolrec.c */
void ridgpolrec(Fimage out_re, Fimage out_im, Fimage in_re, Fimage in_im);
/* src/wave/dyowave2.c */
void dyowave2(int *NumRec, int *StopDecim, int *Ortho, int *Edge,
              int *Precond, int *FilterNorm, Fimage Image, Wtrans2d Output,
              Fsignal Ri, Fimage Edge_Ri);
/* src/wave/precond2d.c */
void precond2d(int *Inverse, Fimage Image, Fimage Output, Fimage Edge_Ri);
/* src/curve/io/kplot.c */
void kplot(Cimage A, char *line, Curves curves, Cimage B, int *d);
/* src/curve/io/kplotasc.c */
void kplotasc(Curve cv);
/* src/curve/io/fkreadasc.c */
Fcurves fkreadasc(void);
/* src/curve/io/readpoly.c */
Polygons readpoly(Cimage image, char *window);
/* src/curve/io/flreadasc.c */
Flists flreadasc(int dim);
/* src/curve/io/fkview.c */
void fkview(Flists in, Ccimage * out, int *sx, int *sy, Flists ref,
            Fimage bg, int *i, char *a, char *s, char *e, int *d, int *g,
            int *c, int *C, char *n, char *window, int *x_0, int *y_0,
            Flist * curve);
/* src/curve/io/fkprintfig.c */
void fkprintfig(Fcurves in, int *d, char *e_flag, char *s_flag, float *m,
                float *r);
/* src/curve/io/kreadfig.c */
Curves kreadfig(void);
/* src/curve/io/fkprintasc.c */
void fkprintasc(Fcurves fcvs);
/* src/curve/io/fkplot.c */
Cimage fkplot(Fcurves cs, Cimage out, char *s_flag);
/* src/curve/io/dkinfo.c */
void dkinfo(Dlists in);
/* src/curve/io/flprintasc.c */
void flprintasc(Flists in, char *v, int *n);
/* src/curve/fkzrt.c */
Fcurves fkzrt(Fcurves cs, float zoom, float angle, float x, float y,
              char *sym);
/* src/curve/flconcat.c */
Flists flconcat(Flists l1, Flists l2, Flists out);
/* src/curve/dsplit_convex.c */
Dlists dsplit_convex(Dlist in, Dlists out, int *ncc, double *eps);
/* src/curve/fkbox.c */
void fkbox(Fcurves cs, float *xmin, float *ymin, float *xmax, float *ymax,
           float *z, Fcurve box);
/* src/curve/fillpolys.c */
void fillpolys(int *dx, int *dy, Polygons polys, Cimage bitmap);
/* src/curve/perimeter.c */
double perimeter(Dlist in, double *min, double *max);
/* src/curve/extract_connex.c */
Fcurves extract_connex(Cimage in, Fcurves curves, int *g);
/* src/curve/flscale.c */
Flists flscale(Flists in, Fsignal s);
/* src/curve/fkcrop.c */
Fcurves fkcrop(float X1, float Y1, float X2, float Y2, Fcurves cs,
               Fcurve box);
/* src/curve/ksplines.c */
void ksplines(int *C, float *Step, Curves ctrl_pts, Curves splines);
/* src/curve/fsplit_convex.c */
Flists fsplit_convex(Flist in, Flists out, int *ncc, double *eps);
/* src/curve/fkcenter.c */
void fkcenter(Fcurves cs, float *xg, float *yg);
/* src/curve/area.c */
double area(Dlist in);
/* src/curve/matching/km_prematchings.c */
Flists km_prematchings(float maxError, Flists dict1, Flists dict2,
                       Flists matchings);
/* src/curve/matching/km_codecurve_ai.c */
Flists km_codecurve_ai(Flist curve, Flist curve_IP, Flist curve_FP,
                       Flist curve_BP, Flists dict, int NC, int NN, float FN);
/* src/curve/matching/km_match_ai.c */
void km_match_ai(float maxError1, float maxError2, float minLength,
                 float minComplex, Flists levlines1, Flists levlines2,
                 Flists dict1, Flists dict2, Flist matchings,
                 Flist matching_pieces);
/* src/curve/matching/km_match_si.c */
void km_match_si(float maxError1, float maxError2, float minLength,
                 float minComplex, Flists levlines1, Flists levlines2,
                 Flists dict1, Flists dict2, Flist matchings,
                 Flist matching_pieces);
/* src/curve/matching/km_flatpoints.c */
Flist km_flatpoints(Flist curve, Flist curve_IP, Flist curve_FP, float angle,
                    float dist);
/* src/curve/matching/km_inflexionpoints.c */
Flist km_inflexionpoints(Flist curve, Flist curve_IP);
/* src/curve/matching/km_createdict_ai.c */
void km_createdict_ai(float *FNorm, int *NNorm, Flists list_curves,
                      Flists dict);
/* src/curve/matching/km_savematchings.c */
void km_savematchings(Flist matching_pieces, Flists levlines1,
                      Flists levlines2, Flists aux1, Flists aux2,
                      Flists pieceaux1, Flists pieceaux2);
/* src/curve/matching/km_createdict_si.c */
void km_createdict_si(float *FNorm, int *NNorm, Flists list_curves,
                      Flists dict);
/* src/curve/matching/km_bitangents.c */
Flist km_bitangents(Flist curve, Flist curve_IP, Flist curve_BP);
/* src/curve/matching/km_codecurve_si.c */
Flists km_codecurve_si(Flist curve, Flist curve_IP, Flist curve_FP,
                       Flist curve_BP, Flists dict, int NC, int NN, float FN);
/* src/curve/smooth/fksmooth.c */
Flist fksmooth(Flist in, int *n, float *std, float *t, char *P);
/* src/curve/smooth/gass.c */
Dlists gass(Dlists in, Dlists out, double *first, double *last, double *eps,
            double *step, int *n, double *r, char *v);
/* src/curve/smooth/iter_fksmooth.c */
Flists iter_fksmooth(Flist in, Flists out, int *niter, int *n, float *std,
                     float *t, char *P);
/* src/curve/smooth/iter_gass.c */
Dlists iter_gass(Dlist in, Dlists out, double *scale, int *niter, double *e,
                 double *s, int *n);
/* src/curve/smooth/iter_gcsf.c */
Dlists iter_gcsf(Dlist in, Dlists out, double *gam, double *last,
                 double *er_area, double *eps, int *n, int *N, char *v);
/* src/curve/smooth/gcsf.c */
Dlists gcsf(Dlists in, Dlists out, double *gam, double *first, double *last,
            double *eps, double *area_sz, int *n, double *r, char *v,
            int *iter, char *conv);
/* src/curve/disc.c */
Curve disc(float r, float *i);
/* src/curve/kspline.c */
void kspline(int *C, float *Step, Curve P, Curve spline);
/* src/curve/circle.c */
Dlist circle(Dlist out, double *r, int *n);
/* src/curve/fillpoly.c */
void fillpoly(int *dx, int *dy, Polygon poly, Cimage bitmap);
/* src/compression/scalar/fiscalq.c */
void fiscalq(int *Print, int *NRow, int *NCol, Cimage Compress,
             Fimage Result, double *Rate);
/* src/compression/scalar/fscalq.c */
void fscalq(int *PrintSNR, int *SmallHeader, int *NStep, float *SStep,
            int *Center, Cimage Compress, Fimage Image, Fimage Result,
            double *MSE, double *SNR, double *Ent, double *RateAr);
/* src/compression/lossless/cvsencode.c */
double cvsencode(int *L, Curves O, Curves C, unsigned int *N, double *B);
/* src/compression/lossless/arencode2.c */
void arencode2(int *Print, long int *Size, int *NSymb, long int *Cap_Histo,
               int *Predic, Fsignal Histo, int *Header, Fimage Input,
               double *Rate, Cimage Output);
/* src/compression/lossless/cvsfrecode.c */
double cvsfrecode(Curves C, unsigned int *N, double *B);
/* src/compression/lossless/cvsorgcode.c */
double cvsorgcode(Curves C, unsigned int *N, double *B);
/* src/compression/lossless/fencode.c */
double fencode(Fimage U);
/* src/compression/lossless/ardecode2.c */
void ardecode2(int *Print, int *NRow, int *NSymb, long int *Cap_Histo,
               int *Predic, Fsignal Histo, Cimage Input, double *Rate,
               Fimage Output);
/* src/compression/ezwave/iezw.c */
void iezw(int *PrintFull, float *WeightFac, Cimage Compress, Wtrans2d Output);
/* src/compression/ezwave/cfiezw.c */
void cfiezw(Fimage Edge_Ri, Fsignal Ri2, float *WeightFac, int *Conv,
            Cimage Compress, Fsignal Ri, Cfimage Output);
/* src/compression/ezwave/ezw.c */
void ezw(int *PrintFull, int *NumRec, float *WeightFac, float *Thres,
         int *Max_Count_AC, int *DistRate, float *Rate, float *PSNR,
         Polygons SelectedArea, Cimage Compress, Wtrans2d Wtrans,
         Wtrans2d Output, char *PtrDRC);
/* src/compression/ezwave/fezw.c */
void fezw(int *NumRec, Fimage Edge_Ri, Fsignal Ri2, int *FilterNorm,
          float *WeightFac, int *DistRate, float *Rate, float *PSNR,
          Polygons SelectedArea, Cimage Output, Fimage Image, Fsignal Ri,
          Fimage QImage, char *PtrDRC);
/* src/compression/ezwave/fiezw.c */
void fiezw(Fimage Edge_Ri, Fsignal Ri2, float *WeightFac, Cimage Compress,
           Fsignal Ri, Fimage Output);
/* src/compression/ezwave/cfezw.c */
void cfezw(int *NumRec, Fimage Edge_Ri, Fsignal Ri2, int *FilterNorm,
           float *WeightFac, int *DistRate, float *Rate, float *GRate,
           float *BRate, float *PSNR, int *Conv, Polygons SelectedArea,
           Cimage Output, Cfimage Image, Fsignal Ri, Cfimage QImage);
/* src/compression/vector/mk_trainset.c */
void mk_trainset(int *Width, int *Height, int *Lap, int *Decim, int *Edge,
                 float *ThresVal1, float *ThresVal2, float *ThresVal3,
                 int *SizeCB, Fimage Image2, Fimage Image3, Fimage Image4,
                 Fimage Image5, Fimage Image6, Fimage Image7, Fimage Image8,
                 Fimage Result2, Fimage Result3, Fimage Result4,
                 Fimage Image, Fimage Result);
/* src/compression/vector/mk_codebook.c */
void mk_codebook(int *Normal, double *Mean, double *Variance, int *Size,
                 int *BlockSize, Fimage CodeBook);
/* src/compression/vector/flbg_adap.c */
void flbg_adap(int *Size, int *Width, int *Height, int *Lap, int *Decim,
               int *Edge, float *ThresVal1, float *ThresVal2,
               float *ThresVal3, Fsignal Weight, int *MultiCB, int *PrintSNR,
               int *Size2, int *Size3, int *Size4, Fimage Image2,
               Fimage Image3, Fimage Image4, Fimage Image5, Fimage Image6,
               Fimage Image7, Fimage Image8, Fimage Output2, Fimage Output3,
               Fimage Output4, Fimage Image1, Fimage Output);
/* src/compression/vector/flbg.c */
void flbg(int *Size, int *Width, int *Height, int *Lap, int *Decim,
          int *Edge, Fsignal Weight, int *MultiCB, int *PrintSNR,
          Fimage Image2, Fimage Image3, Fimage Image4, Fimage Image5,
          Fimage Image6, Fimage Image7, Fimage Image8, int *InitRandCB,
          int *RepeatRand, int *NResCB, int *NResResCB, Fimage ResCodeBook,
          Fimage ResResCodeBook, Fimage Image1, Fimage CodeBook, float *MSE);
/* src/compression/vector/flbg_train.c */
void flbg_train(int *Size, Fsignal Weight, int *MultiCB, Fimage InitCodeBook,
                int *NResCB, int *NResResCB, Fimage ResCodeBook,
                Fimage ResResCodeBook, int *PrintSNR, Fimage TrainSet,
                float *MSE, Fimage CodeBook);
/* src/compression/vector/fvq.c */
void fvq(int *PrintSNR, int *SmallHeader, int *BitMapCode, int *RateDist,
         int *NCB1, int *NCB2, int *NCB3, int *NCB4, Fimage CodeBook2,
         Fimage CodeBook3, Fimage CodeBook4, int *NResCB1, int *NResCB2,
         int *NResCB3, int *NResCB4, Fimage ResCodeBook1,
         Fimage ResCodeBook2, Fimage ResCodeBook3, Fimage ResCodeBook4,
         int *NResResCB1, int *NResResCB2, Fimage ResResCodeBook1,
         Fimage ResResCodeBook2, Cimage Compress, Fimage Image,
         Fimage CodeBook1, Fimage Result, double *MSE, double *SNR,
         double *Entropy, double *RateAr);
/* src/compression/vector/fivq.c */
void fivq(int *Print, int *NRow, int *NCol, Fimage CodeBook2,
          Fimage CodeBook3, Fimage CodeBook4, Fimage ResCodeBook1,
          Fimage ResCodeBook2, Fimage ResCodeBook3, Fimage ResCodeBook4,
          Fimage ResResCodeBook1, Fimage ResResCodeBook2, Cimage Compress,
          Fimage CodeBook1, Fimage Result, double *Rate);
/* src/compression/vqwave/fwlbg_adap.c */
void fwlbg_adap(int *NumRecMax, int *Level, int *Orient, Fimage Edge_Ri,
                Fsignal Ri2, int *FilterNorm, int *StopDecim, int *Width,
                int *Height, int *MultiCB, int *Lap, int *Sizec1,
                int *Sizec2, int *Sizec3, float *ThresVal1, float *ThresVal2,
                float *ThresVal3, Fimage OldCodeBook,
                Fimage OldAdapCodeBook2, Fimage OldAdapCodeBook3,
                Fimage Output2, Fimage Output3, Fimage Image2, Fimage Image3,
                Fimage Image4, Fimage ResCodeBook, Fimage ResResCodeBook,
                Fimage Image1, Fsignal Ri, Fimage Output1);
/* src/compression/vqwave/fwvq.c */
void fwvq(int *NumRec, Fimage Edge_Ri, Fsignal Ri2, int *FilterNorm,
          float *WeightFac, int *NumRecScal, int *NStep, int *MultiCB,
          Fimage CodeBook2, Fimage CodeBook3, Fimage ResCodeBook1,
          Fimage ResCodeBook2, Fimage ResCodeBook3, Fimage ResResCodeBook1,
          Fimage ResResCodeBook2, int *DistRate, float *TargRate,
          Cimage Output, Fimage Image, Fimage CodeBook1, Fsignal Ri,
          Fimage QImage);
/* src/compression/vqwave/wlbg_adap.c */
void wlbg_adap(int *NumRecMax, int *Level, int *Orient, int *Dyadic,
               int *Lap, int *Edge, int *Width, int *Height, int *MultiCB,
               int *Sizec1, int *Sizec2, int *Sizec3, float *ThresVal1,
               float *ThresVal2, float *ThresVal3, Wtrans2d OldCodeBook,
               Wtrans2d OldAdapCodeBook2, Wtrans2d OldAdapCodeBook3,
               Wtrans2d * Output2, Wtrans2d * Output3, Wtrans2d TrainWtrans2,
               Wtrans2d TrainWtrans3, Wtrans2d TrainWtrans4,
               Wtrans2d ResCodeBook, Wtrans2d ResResCodeBook,
               Wtrans2d TrainWtrans1, Wtrans2d * Output1);
/* src/compression/vqwave/fwivq.c */
void fwivq(Fimage Edge_Ri, Fsignal Ri2, int *FilterNorm, float *WeightFac,
           Fimage CodeBook2, Fimage CodeBook3, Fimage ResCodeBook1,
           Fimage ResCodeBook2, Fimage ResCodeBook3, Fimage ResResCodeBook1,
           Fimage ResResCodeBook2, Cimage Compress, Fimage CodeBook1,
           Fsignal Ri, Fimage Output);
/* src/image/io/ccview.c */
void ccview(Ccimage image, int *x0, int *y0, float *zoom, int *order,
            int *no_refresh, char *window);
/* src/image/io/cccopy.c */
void cccopy(Ccimage Input, Ccimage * Output);
/* src/image/io/ccmview.c */
void ccmview(Ccmovie input, int *x0, int *y0, float *zoom, int *order,
             int *loop, int *pause, char *window);
/* src/image/io/ccopy.c */
void ccopy(Cimage Input, Cimage * Output);
/* src/image/io/fline_extract.c */
void fline_extract(char *cflag, Fimage Image, Fsignal Signal, long int Index);
/* src/image/io/fprintasc.c */
void fprintasc(Fimage u, int *x1, int *y1, int *x2, int *y2, char *verbose);
/* src/image/io/cmview.c */
void cmview(Cmovie input, int *x0, int *y0, float *zoom, int *order,
            int *loop, int *pause, char *window);
/* src/image/io/cview.c */
void cview(Cimage input, int *x0, int *y0, float *zoom, int *order,
           int *no_refresh, char *window);
/* src/image/io/cfgetchannels.c */
void cfgetchannels(Cfimage Image, Fimage RImage, Fimage GImage,
                   Fimage BImage);
/* src/image/io/creadasc.c */
void creadasc(Cimage u, int dx, int dy);
/* src/image/io/fview.c */
void fview(Fimage input, int *x0, int *y0, float *zoom, int *order,
           int *no_refresh, Wframe * window, char *linear, float *m,
           float *M, float *d);
/* src/image/io/cfchgchannels.c */
void cfchgchannels(int *Conv, int *Inverse, int *Norm, Cfimage Image,
                   Cfimage TImage);
/* src/image/io/cline_extract.c */
void cline_extract(char *cflag, Cimage Image, Fsignal Signal, long int Index,
                   Cimage OutImage);
/* src/image/io/freadasc.c */
void freadasc(Fimage u, int *Dx, int *Dy);
/* src/image/io/flip.c */
void flip(Fimage in1, Fimage in2, float *z, int *o);
/* src/image/io/cprintasc.c */
void cprintasc(Cimage u, int *x1, int *y1, int *x2, int *y2, char *verbose);
/* src/image/io/cfputchannels.c */
void cfputchannels(Fimage RImage, Fimage GImage, Fimage BImage,
                   Cfimage Image);
/* src/image/io/ccputstring.c */
Ccimage ccputstring(Ccimage in, int x, int y, int *c, int *C, float *r,
                    char *str);
/* src/image/io/fcopy.c */
void fcopy(Fimage Input, Fimage * Output);
/* src/image/level_lines/flstb_dualchain.c */
void flstb_dualchain(Shapes pTree, Shape pShape, Flist pBoundary,
                     char *ctabtabSaddleValues);
/* src/image/level_lines/ml_fml.c */
Fmorpho_line ml_fml(Morpho_line lline);
/* src/image/level_lines/ml_draw.c */
Cimage ml_draw(Mimage m_image, Cimage bimage, char *border, Cmovie movie);
/* src/image/level_lines/fml_ml.c */
Morpho_line fml_ml(Fmorpho_line flline);
/* src/image/level_lines/flstb_dual.c */
void flstb_dual(Shapes pTree, Shapes pDualTree);
/* src/image/level_lines/ml_decompose.c */
Mimage ml_decompose(Mimage m_image_in, int *ml_opt, char *m_flag,
                    Fimage image_in);
/* src/image/level_lines/mscarea.c */
int mscarea(char *connex8, Cimage U, Cimage O, int *stoparea, int a, int b,
            int x0, int y0);
/* src/image/level_lines/flst_boundary.c */
Flist flst_boundary(Shapes pTree, Shape pShape, Flist pBoundary);
/* src/image/level_lines/tjmap.c */
int tjmap(char *connex8, char *values, int *tarea, int *tquant, Cimage U,
          Cimage J);
/* src/image/level_lines/ll_distance.c */
void ll_distance(Fimage in, Fimage out, float *level, int *MaxDist,
                 Fimage nearest);
/* src/image/level_lines/ll_remove.c */
Mimage ll_remove(Mimage mimage, int *L);
/* src/image/level_lines/llmap.c */
Cimage llmap(short int *ls, char *tmap, Cimage input, Cimage output);
/* src/image/level_lines/flstb_boundary.c */
Flist flstb_boundary(int *pPrecision, Fimage pImage, Shapes pTree,
                     Shape pShape, Flist pDualchain, Flist pBoundary,
                     char *ctabtabSaddleValues);
/* src/image/level_lines/ml_reconstruct.c */
Fimage ml_reconstruct(char *v_flag, Mimage m_image);
/* src/image/level_lines/cml_draw.c */
Ccimage cml_draw(Cmimage cmimage, Ccimage bimage, char *border,
                 Ccmovie movie);
/* src/image/level_lines/flstb_quantize.c */
void flstb_quantize(Fsignal pLevels, float *pOffset, float *pStep,
                    Shapes pTree, Shapes pQuantizedTree);
/* src/image/level_lines/cml_reconstruct.c */
Cfimage cml_reconstruct(char *v_flag, Cmimage m_image);
/* src/image/level_lines/flst_bilinear.c */
void flst_bilinear(int *pMinArea, Fimage pImageInput, Shapes pTree);
/* src/image/level_lines/cml_decompose.c */
Cmimage cml_decompose(Cmimage cmimage_in, int *ml_opt, int *L,
                      Cfimage image_in);
/* src/image/level_lines/llremove.c */
Fimage llremove(Fimage Input, int *Lmin, int *Lmax, char *invert);
/* src/image/level_lines/fsaddles.c */
int fsaddles(Fimage pImage, Fimage pSaddlesImage);
/* src/image/level_lines/flst_pixels.c */
void flst_pixels(Shapes pTree);
/* src/image/level_lines/cll_remove.c */
Cmimage cll_remove(Cmimage cmimage, int *L);
/* src/image/level_lines/tjpoint.c */
int tjpoint(char *connex8, int *tarea, int *tquant, Cimage U, int x0,
            int y0, int *lambda, int *mu, int *xlambda, int *ylambda,
            int *xmu, int *ymu, unsigned char *M, int *P);
/* src/image/level_lines/ml_extract.c */
void ml_extract(float *level, Fsignal levels, int *opt, Cimage c_out,
                char *m_flag, Fimage image_org, Mimage m_image);
/* src/image/level_lines/flst_reconstruct.c */
void flst_reconstruct(Shapes pTree, Fimage pFloatImageOutput);
/* src/image/level_lines/ll_extract.c */
Flists ll_extract(Fimage in, Fsignal levels, float *offset, float *step,
                  int *prec, int *min_area, Shapes tree, char *z);
/* src/image/level_lines/flst.c */
void flst(int *pMinArea, Fimage pImageInput, Shapes pTree);
/* src/image/level_lines/flstb_tv.c */
void flstb_tv(float *pScale, float *pQuantizationLevel,
              int *pQuantizationCurve, Shapes pTree, Fimage pImage);
/* src/image/level_lines/llview.c */
void llview(Fimage in, float *z, int *s, int *l, int *d, int *i, int *b,
            Ccimage * out, char *n);
/* src/image/shape_recognition/sr_genhypo.c */
Fcurves sr_genhypo(Fimage sg, Cimage img, Fimage hypo_list);
/* src/image/shape_recognition/sr_normalize.c */
Fcurves sr_normalize(Fcurves in);
/* src/image/shape_recognition/sr_signature.c */
Fimage sr_signature(Fcurves in, int *n, Fimage base);
/* src/image/shape_recognition/sr_distance.c */
float sr_distance(Fcurves Shape1, Fcurves Shape2);
/* src/image/misc/ccdisocclusion.c */
void ccdisocclusion(Ccimage Input, Ccimage Output, Fimage Holes,
                    int *energy_type, char *angle);
/* src/image/misc/lsnakes_demo.c */
void lsnakes_demo(Cimage u, Cmovie out, int *Niter, int *Nframes,
                  float *thre, float *force);
/* src/image/misc/drawocclusion.c */
void drawocclusion(Cimage image, Fimage labelimage, Cimage * hole_image,
                   char *window, int *zoom, int *gray);
/* src/image/misc/lsnakes.c */
Fimage lsnakes(Fimage in, Fimage ref, int *Niter, float *thre, float *force);
/* src/image/misc/mac_snakes.c */
Dlists mac_snakes(Fimage u, Dlists in, int *niter, double *step,
                  double *power, char *v, float *V);
/* src/image/misc/skeleton.c */
Cmovie skeleton(int *iteration, int *infsup_iter, int *extremite,
                char *average, Fmovie cmovie, Cimage image, Cmovie output);
/* src/image/misc/thinning.c */
Cimage thinning(int *niter, char *inv, Cmovie M, Cimage imageI,
                Cimage imageO);
/* src/image/misc/cdisc.c */
Cimage cdisc(Cimage out, int nx, int ny, float *x, float *y, float *r);
/* src/image/misc/disocclusion.c */
void disocclusion(Cimage Input, Cimage Output, Fimage Holes,
                  int *energy_type, char *angle);
/* src/image/misc/emptypoly.c */
Cimage emptypoly(Cimage A, Cimage B);
/* src/image/values/frank.c */
void frank(Fimage u, Fimage rank, Flist g, float *w, int *c);
/* src/image/values/bicontrast.c */
void bicontrast(Fimage u, Fimage v, char *verb, Flist r, Flist g, Fimage out);
/* src/image/values/amle_init.c */
void amle_init(Fimage in, float delta, Fimage out);
/* src/image/values/amle3d_init.c */
void amle3d_init(Fmovie in, float delta, Fmovie out);
/* src/image/values/cmnoise.c */
Cmovie cmnoise(Cmovie in, float *std, float *p, float *q);
/* src/image/values/ccontrast.c */
void ccontrast(Cimage in, Cimage out, float *g);
/* src/image/values/cfquant.c */
void cfquant(Cfimage A, Cfimage Q, int Mr, int Mg, int Mb, char *left);
/* src/image/values/amle3d.c */
void amle3d(int *num, Fmovie init, Fmovie in, Fmovie out);
/* src/image/values/binarize.c */
void binarize(Fimage in, Cimage out, float *t, char *i);
/* src/image/values/fthre.c */
Fimage fthre(Fimage in, Fimage out, char *norm, char *N, float *min,
             float *max, float *p, float *q, float *d, char *aff,
             char *linear);
/* src/image/values/flgamma.c */
Flist flgamma(Flist out, float *g, float *s, int *n, float *f);
/* src/image/values/fhisto.c */
Fsignal fhisto(Fimage in, Fsignal out, float *l, float *r, int *n, float *s,
               char *t, Fimage w);
/* src/image/values/fquant.c */
float fquant(Fimage A, Fimage Q, int M, char *left, float *min, float *max);
/* src/image/values/fnoise.c */
void fnoise(Fimage u, Fimage v, float *std, float *p, float *q, char *n_flag);
/* src/image/values/amle.c */
void amle(Cimage Init, Cimage Input, Fimage Output, float *omega, int *n,
          float *ht, float *mse);
/* src/image/values/fcontrast.c */
Fimage fcontrast(Fimage in, Flist g);
/* src/image/values/ccontrast_local.c */
void ccontrast_local(Cimage in, Cimage out, int *d, char *n_flag);
/* src/image/values/fvalues.c */
Fsignal fvalues(char *i_flag, Fsignal multiplicity, Fimage image_rank,
                Fimage image_in);
/* src/image/values/cnoise.c */
void cnoise(Cimage u, Cimage v, float *std, float *p, float *q, char *n_flag);
/* src/image/values/chisto.c */
Fsignal chisto(Cimage A, Fsignal S, char *i_flag);
/* src/image/domain/ccmcollect.c */
Ccimage ccmcollect(Ccmovie in, Ccimage out, int *n, int *i, int *o, int *c,
                   int *x, int *y);
/* src/image/domain/cfunzoom.c */
Cfimage cfunzoom(Cfimage in, Cfimage out, float *z, int *o, float *tx,
                 float *ty);
/* src/image/domain/fzrt.c */
void fzrt(Fimage in, Fimage out, float zoom, float angle, float x, float y,
          int *o, float *p, float *b);
/* src/image/domain/fproj.c */
void fproj(Fimage in, Fimage out, int *sx, int *sy, float *bg, int *o,
           float *p, char *i, float X1, float Y1, float X2, float Y2,
           float X3, float Y3, float *x4, float *y4);
/* src/image/domain/cextcenter.c */
Cimage cextcenter(int *Fact, Cimage Image);
/* src/image/domain/cmcollect.c */
Cimage cmcollect(Cmovie in, Cimage out, int *n, int *i, int *o, int *c,
                 int *x, int *y);
/* src/image/domain/fmaskrot.c */
Fimage fmaskrot(Fimage in, Fimage out, float *bg, float *s);
/* src/image/domain/fshift.c */
Fimage fshift(Fimage in, Fimage out, int *x, int *y, int *h, int *i);
/* src/image/domain/frot.c */
void frot(Fimage in, Fimage out, float *a, float *b, char *k_flag);
/* src/image/domain/cmzoom.c */
void cmzoom(Cmovie Input, Cmovie Output, char *x_flg, char *y_flg,
            float *factor, int *o, char *i_flg);
/* src/image/domain/csample.c */
Cimage csample(Cimage in, Cimage out, double step);
/* src/image/domain/fextract.c */
Fimage fextract(float *b, Fimage in, Fimage bg, Fimage out, int X1, int Y1,
                int X2, int Y2, int *Xc, int *Yc, char *r);
/* src/image/domain/flocal_zoom.c */
Fimage flocal_zoom(Fimage Input, int *X, int *Y, int *W, int *factor);
/* src/image/domain/ccmzoom.c */
void ccmzoom(Ccmovie Input, Ccmovie Output, char *x_flg, char *y_flg,
             float *factor, int *o, char *i_flg);
/* src/image/domain/fsample.c */
Fimage fsample(Fimage in, Fimage out, double step, double *delta, int *norm);
/* src/image/domain/finvspline.c */
void finvspline(Fimage in, int order, Fimage out);
/* src/image/domain/ccextract.c */
Ccimage ccextract(int *b, Ccimage in, Ccimage bg, Ccimage out, int X1, int Y1,
                  int X2, int Y2, int *Xc, int *Yc, char *r);
/* src/image/domain/czoom.c */
Cimage czoom(Cimage in, Cimage out, char *x_flg, char *y_flg, float *zoom,
             int *o, char *i_flg);
/* src/image/domain/fdirspline.c */
Fimage fdirspline(Fimage in, int n, Fimage out);
/* src/image/domain/fzoom.c */
Fimage fzoom(Fimage in, Fimage out, char *x_flg, char *y_flg, float *zoom,
             int *o, float *p, char *i_flg);
/* src/image/domain/cextract.c */
Cimage cextract(int *b, Cimage in, Cimage bg, Cimage out, int X1, int Y1,
                int X2, int Y2, int *Xc, int *Yc, char *r);
/* src/image/domain/cmextract.c */
Cmovie cmextract(int *b, Cmovie in, Cmovie bg, int X1, int Y1, int T1, int X2,
                 int Y2, int T2, int *Xc, int *Yc, int *Tc, char *r);
/* src/image/domain/funzoom.c */
Fimage funzoom(Fimage in, Fimage out, float *z, int *o, float *tx, float *ty);
/* src/image/domain/fcrop.c */
Fimage fcrop(Fimage in, Fimage out, float *sx, float *sy, float *z, float *bg,
             int *o, float *p, float X1, float Y1, float X2, float Y2);
/* src/image/domain/clocal_zoom.c */
Cimage clocal_zoom(Cimage Input, int *X, int *Y, int *W, int *factor);
/* src/image/domain/cczoom.c */
Ccimage cczoom(Ccimage in, Ccimage out, char *x_flg, char *y_flg, float *zoom,
               int *o, char *i_flg);
/* src/image/domain/cccrop.c */
Ccimage cccrop(Ccimage in, Ccimage out, float *sx, float *sy, float *z,
               unsigned char *bg, int *o, float *p, float X1, float Y1,
               float X2, float Y2);
/* src/image/domain/cclocal_zoom.c */
Ccimage cclocal_zoom(Ccimage Input, int *X, int *Y, int *W, int *factor);
/* src/image/domain/cmparitysep.c */
Cmovie cmparitysep(Cmovie u, char *e, char *l);
/* src/image/domain/cfextcenter.c */
Cfimage cfextcenter(int *Fact, Cfimage Image);
/* src/image/motion/ofdraw.c */
Cmovie ofdraw(Fmovie U, Fmovie V, int *a, float *m, int *p, int *z, float *h);
/* src/image/motion/ws_flow.c */
void ws_flow(float *percent, int *n, float *tau, float *lambda, float *eps,
             float *alpha, Fmovie norm, Fmovie movie, Fmovie wsU, Fmovie wsV);
/* src/image/motion/hs_flow.c */
void hs_flow(int *niter, float *alpha, Fmovie in, Fmovie xflow, Fmovie yflow);
/* src/image/motion/motionseg.c */
void motionseg(int *n, float *prec, float *e, float *alphac, float *alphabr,
               float *alphacr, float *seu, Fmovie N, Fmovie C, Fimage B);
/* src/image/detection/ll_boundaries.c */
Flists ll_boundaries(Fimage in, Shapes tree, float *eps, char *all,
                     float *step, int *precision, char *z, char *weak);
/* src/image/detection/ll_edges.c */
Flists ll_edges(Fimage in, Shapes tree, float *eps, float *step,
                int *precision, char *z);
/* src/image/detection/falign.c */
Flist falign(Fimage u, int *d, int *nl, double *eps, float *g, char *all,
             Flists crv);
/* src/image/detection/harris.c */
Flist harris(Fimage in, float *k, float *g, int *size, double *t);
/* src/image/detection/canny.c */
Fimage canny(float *alpha, Fimage IN, Fimage OUT);
/* src/image/detection/ll_boundaries2.c */
Flists ll_boundaries2(Fimage in, float *eps, Shapes tree, float *step,
                      int *prec, float *std, float *hstep, char *all,
                      int *visit, char *loc, Fimage image_out,
                      Shapes keep_tree);
/* src/image/detection/falign_mdl.c */
Fimage falign_mdl(Fimage u, int *d, int *nd, int *no_mdl, int *nl,
                  double *eps, float *g, Flists crv);
/* src/image/detection/vpsegplot.c */
void vpsegplot(Fimage image, Fimage allsegs, Flist vpoints, Flists vsegs,
               int n, Flists crv, char *lines);
/* src/image/detection/vpoint.c */
int vpoint(Fimage imagein, Fimage allsegs, Flist output, Flists segs,
           double *eps, char *all, char *masked, char *verbose,
           int *maskedVPs);
/* src/image/seg/msegct.c */
Cimage msegct(Fsignal weight, int *sgrid, int *nb_of_regions, float *lambda,
              Curves curves, Fmovie u, int *f_nb_of_regions, float *f_lambda,
              Fmovie orig_data);
/* src/image/seg/mschannel.c */
void mschannel(int *N, int *S, int *W, int *p, Fimage in, Fmovie mov);
/* src/image/seg/segct.c */
Cimage segct(int *sgrid, int *nb_of_regions, float *lambda, Curves curves,
             Fimage u, int *f_nb_of_regions, float *f_lambda,
             Fimage image_org);
/* src/image/seg/one_levelset.c */
void one_levelset(float *level, Cimage cb, Fpolygons pb, Fimage fu,
                  Cimage cu, Fimage image_org);
/* src/image/seg/segtxt.c */
Cimage segtxt(int *N, int *S, int *W, int *p, int *nr, Fimage in, Fmovie mov);
/* src/image/filter/fsharpen.c */
Fimage fsharpen(Fimage A, Fimage B, float *p);
/* src/image/filter/cfsharpen.c */
Cfimage cfsharpen(Cfimage A, Cfimage B, float *p, char *LumOnly);
/* src/image/filter/tvdenoise.c */
Fimage tvdenoise(Fimage in, Fimage out, double *s, char *c, char *v,
                 double *e, int *n, double *W, Fimage ref, double *eps,
                 double *p);
/* src/image/filter/resthline.c */
void resthline(Cimage u, Cimage v);
/* src/image/filter/osamss.c */
void osamss(char *isotrop, char *power, float *Step, float *MinGrad,
            float *firstScale, float *lastScale, Fimage input, Fimage output);
/* src/image/filter/prolate.c */
float prolate(int s, float d, int *m, Fimage ker);
/* src/image/filter/cfdiffuse.c */
void cfdiffuse(float *deltat, float *epsilon, Cfimage in, Cfimage out,
               Fsignal MDiag0, Fsignal MDiag1, Fsignal U0, Cfimage Yimage,
               Cfimage Vimage, Fimage L2h, Fimage L2v);
/* src/image/filter/tvdeblur.c */
Fimage tvdeblur(Fimage in, Fimage out, Fimage ker, double *s, char *c,
                char *v, double *e, int *n, double *W, Fimage ref,
                double *eps, double *p);
/* src/image/filter/amss.c */
void amss(char *isotrop, char *power, float *Step, float *MinGrad,
          float *outputStep, float *firstScale, float *lastScale,
          Fimage image, Fimage * imageD, Fimage * imageG, Fimage * imageC,
          Cmovie cmovieD, Cmovie cmovieG, Cmovie cmovieC, char *no_norm);
/* src/image/filter/prolatef.c */
float prolatef(int s, float d, int *n, Fimage ker, int *p, int *c);
/* src/image/filter/flipschitz.c */
Fimage flipschitz(Fimage in, float lip, Fimage out, float *r, Curve s, int *n,
                  char *i);
/* src/image/filter/mam.c */
void mam(Cmovie in, Cmovie out, float *ptime, float *ppower, int *n_iter,
         short int *pMAXvit, short int *pMINvit, short int *pfmxa);
/* src/image/filter/opening.c */
Cimage opening(Cimage u, Cimage v, float *r, Curve s, int *n, char *i);
/* src/image/filter/rotaffin.c */
void rotaffin(int *Nrota, int *Naffi, int *Size, int *Type, float *Area,
              int *Definition, double *OptSym, float *Factor, Cimage cimage,
              Cmovie cmovie);
/* src/image/filter/cfmdiffuse.c */
void cfmdiffuse(float *deltat, int *N, float *epsilon, Cfimage in,
                Cfmovie out);
/* src/image/filter/infsup.c */
Cimage infsup(int *Niter, float *deginf, float *degsup, char *average,
              Cimage image, Fmovie fmovie, Cimage output);
/* src/image/filter/heat.c */
Fimage heat(Fimage in, Fimage out, int *n, float *s);
/* src/image/filter/nlmeans.c */
Fimage nlmeans(Fimage in, Fimage out, double *h, double *a, int *s, int *d,
               double *c);
/* src/image/filter/fgrain.c */
void fgrain(int *pMinArea, Fimage pFloatImageInput, Fimage pFloatImageOutput);
/* src/image/filter/median.c */
Cimage median(Cimage u, Cimage v, float *r, Curve s, int *n);
/* src/image/filter/fconvol.c */
void fconvol(Fimage in, Fimage filtre, Fimage out);
/* src/image/filter/forder.c */
void forder(int *e, int *N, Fimage in, Fimage out);
/* src/image/filter/shock.c */
Fimage shock(Fimage in, int *n, float *s);
/* src/image/filter/tvdenoise2.c */
Fimage tvdenoise2(Fimage in, Fimage out, double *s, int *v, int *n, double *r,
                  double *W, int *V);
/* src/image/filter/fsepconvol.c */
Fimage fsepconvol(Fimage in, Fimage out, Fsignal xker, Fsignal yker,
                  int *width, float *std, int *b);
/* src/image/filter/ll_sharp.c */
void ll_sharp(float *pPercentIncreaseArea, Fimage pFloatImageInput,
              Fimage pFloatImageOutput);
/* src/image/filter/erosion.c */
Cimage erosion(Cimage u, Cimage v, float *r, Curve s, int *n, char *i);
/* src/image/filter/fsmooth.c */
void fsmooth(int *S, int *W, Fimage in, Fimage out);
/* src/image/operations/fmse.c */
void fmse(Fimage Img1, Fimage Img2, int *Norm, char *PsnrFlg, double *SNR,
          double *PSNR, double *MSE, double *MRD);
/* src/image/operations/fop.c */
void fop(Fimage B, Fimage out, Fimage A, float *a, char *plus, char *minus,
         char *times, char *divide, char *dist, char *norm, char *inf,
         char *sup, char *less, char *greater, char *equal);
/* src/image/operations/fpsnr255.c */
double fpsnr255(int *Norm, Fimage image);
/* src/image/operations/finfo.c */
void finfo(Fimage u);
/* src/image/operations/faxpb.c */
void faxpb(Fimage in, Fimage out, float *a, float *s, float *b, float *m,
           char *k, float *M, float *N);
/* src/image/operations/fsize.c */
void fsize(Fimage in);
/* src/image/operations/fdiff.c */
void fdiff(char *absd, Fimage A, Fimage B, Fimage O);
/* src/image/operations/fabso.c */
void fabso(Fimage A, Fimage O);
/* src/image/operations/frthre.c */
void frthre(Fimage A, Fimage B, float *noise);
/* src/image/operations/fconst.c */
Fimage fconst(float g, int x, int y);
/* src/image/operations/fmask.c */
void fmask(Fimage mask, Fimage A, Fimage B, Fimage out, char *i_flag, int *v,
           float *c);
/* src/image/operations/fmean.c */
float fmean(Fimage A);
/* src/image/operations/fentropy.c */
float fentropy(Fimage A);
/* src/image/operations/fnorm.c */
float fnorm(Fimage in, Fimage ref, float *p, char *s, char *v, int *b,
            char *n, float *t);
/* src/image/operations/fvar.c */
float fvar(Fimage A, char *e, char *s);
/* src/image/operations/fadd.c */
Fimage fadd(Fimage A, Fimage B, Fimage C, float *min, float *max, char *a);
/* src/image/operations/cfdiff.c */
void cfdiff(char *absd, Cfimage A, Cfimage B, Cfimage O);
/* src/image/operations/fpset.c */
Fimage fpset(Fimage in, int x, int y, float g);
/* src/image/operations/fderiv.c */
void fderiv(Fimage in, Fimage curv, Fimage anti, Fimage canny_img,
            Fimage laplacian, Fimage gradx, Fimage grady, Fimage gradn,
            Fimage gradp, float *MinGrad, int *nsize);
/* src/image/fourier/fftrot.c */
void fftrot(Fimage in, Fimage out, float *a, float *x, float *y,
            char *o_flag);
/* src/image/fourier/fft2dpol.c */
void fft2dpol(Fimage in_re, Fimage in_im, Fimage out1, Fimage out2,
              char *i_flag);
/* src/image/fourier/frandphase.c */
void frandphase(Fimage in, Fimage out, char *i_flag);
/* src/image/fourier/fftconvol.c */
Fimage fftconvol(Fimage in, Fimage filter, Fimage out);
/* src/image/fourier/fft2dshrink.c */
Fimage fft2dshrink(Fimage in, Fimage out, float *p, char *v);
/* src/image/fourier/fftgrad.c */
void fftgrad(Fimage in, Fimage gradx, Fimage grady, Fimage gradn,
             Fimage gradp);
/* src/image/fourier/fkeepphase.c */
void fkeepphase(Fimage in, Fimage mod, Fimage out);
/* src/image/fourier/fshrink2.c */
void fshrink2(Fimage in, Fimage out);
/* src/image/fourier/wiener.c */
Fimage wiener(Fimage in, Fimage out, Fimage kernel, Fsignal rad, float *g,
              float *lambda);
/* src/image/fourier/fftzoom.c */
void fftzoom(Fimage in, Fimage out, float *z, char *i_flag);
/* src/image/fourier/fhamming.c */
void fhamming(Fimage in, Fimage out);
/* src/image/fourier/fft2dview.c */
void fft2dview(int *type, char *h_flag, Fimage in, Fimage out, char *i_flag,
               float *d);
/* src/image/fourier/fft2drad.c */
Fsignal fft2drad(Fimage in_re, Fimage in_im, Fsignal out, char *l_flag,
                 int *size);
/* src/image/fourier/fsym2.c */
void fsym2(Fimage in, Fimage out, char *i);
/* src/image/fourier/fft2d.c */
void fft2d(Fimage in_re, Fimage in_im, Fimage out_re, Fimage out_im,
           char *i_flag);
/* src/signal/smse.c */
void smse(Fsignal Sig1, Fsignal Sig2, int *Norm, double *SNR, double *PSNR,
          double *MSE, double *MRD);
/* src/signal/stvrestore.c */
Fsignal stvrestore(char *relax, int *Niter, double *alpha, Fsignal O,
                   Fimage I, Fsignal V, Fsignal M, Fsignal u);
/* src/signal/w1threshold.c */
void w1threshold(float *P, float *T, float *D, char *s, char *o, Fsignal M,
                 Wtrans1d in, Wtrans1d * out);
/* src/signal/sshrink2.c */
Fsignal sshrink2(Fsignal in, Fsignal out);
/* src/signal/sreadasc.c */
Fsignal sreadasc(Fsignal out, int *n, float *s, float *t, float *g);
/* src/signal/splot.c */
void splot(Fsignal in, int *x_0, int *y_0, int *sx, int *sy, int *no_refresh,
           char *window, Ccimage * out, char *n);
/* src/signal/entropy.c */
void entropy(Fsignal Histo, double *Entropy);
/* src/signal/sintegral.c */
Fsignal sintegral(Fsignal in, char *n, char *r);
/* src/signal/sop.c */
void sop(Fsignal B, Fsignal out, Fsignal A, float *a, char *plus,
         char *minus, char *times, char *divide, char *dist, char *norm,
         char *inf, char *sup, char *less, char *greater, char *equal);
/* src/signal/sconst.c */
void sconst(int *size, float *A, Fsignal signal);
/* src/signal/sinfo.c */
void sinfo(Fsignal s);
/* src/signal/saxpb.c */
void saxpb(Fsignal in, Fsignal out, float *a, float *s, float *b, float *m,
           char *k, float *M);
/* src/signal/ssinus.c */
void ssinus(int *size, float *A, float *D, Fsignal signal);
/* src/signal/fft1dshrink.c */
Fsignal fft1dshrink(Fsignal in, Fsignal out, float *p, char *v);
/* src/signal/snorm.c */
float snorm(Fsignal in, Fsignal ref, float *p, char *s, char *v, int *b,
            char *n, float *t);
/* src/signal/snoise.c */
Fsignal snoise(Fsignal in, Fsignal out, float *std, float *t, char *n_flag);
/* src/signal/sderiv.c */
void sderiv(Fsignal A, Fsignal B);
/* src/signal/fft1d.c */
void fft1d(Fsignal Xr, Fsignal Xi, Fsignal Yr, Fsignal Yi, char *inverse);
/* src/signal/sdirac.c */
void sdirac(int *size, float *A, Fsignal signal);
/* src/signal/fct1d.c */
void fct1d(Fsignal X, Fsignal Y, char *inverse);
/* src/signal/sprintasc.c */
void sprintasc(Fsignal s, int *i1, int *i2, char *verbose);
/* src/signal/sgauss.c */
Fsignal sgauss(float std, Fsignal out, int *size, float *prec);

#endif                          /* !_MW3_MODULES_ */
