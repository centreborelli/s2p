#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif


// gcc -o anatomy2lowe anatomy2lowe.c -std=c99




int main(int argc, char** argv)
{
    if ((argc != 2)&&(argc != 5)){
        fprintf(stderr, "usage: anatomy2lowe keypoints [nhist(4) nori(8)]\n");
        return EXIT_FAILURE;
    }
    // length of the feature vector
    int nhist;
    int nori;
    if ( argc == 4){
        nhist = atoi(argv[argc-2]);
        nori = atoi(argv[argc-1]);
    }
    else{
        nhist = 4;
        nori = 8;
    }

    // How many keypoints in the file
    char text[4096];
    FILE* f =  fopen(argv[1],"r");
    int nkeys = 0;
    while(fgets(text, 4096, f) != NULL)
        nkeys += 1;
    fclose(f);

    // memory allocation
    int ndescr = nhist*nhist*nori;
    int length = 4+ndescr;
    float* l = malloc( nkeys*(length)*sizeof(*l));

    // read keypoints
    f =  fopen(argv[1],"r");
    int pos, read;
    for(int n = 0; n < nkeys; n++){
        fgets(text, 4096, f);
        pos = 0; read = 0;
        sscanf(text+pos,"%f %f %f %f %n", &l[n*length]
                                        , &l[n*length+1]
                                        , &l[n*length+2]
                                        , &l[n*length+3]
                                        , &read);
        pos += read;
        for(int m = 0; m < ndescr; m++){
            sscanf(text+pos,"%f %n", &l[n*length+4+m],&read);
            pos += read;
        }
    }
    fclose(f);

    // write keypoints in standard output
    for(int n = 0; n < nkeys; n++){
        // conversion orientation
        float ori = l[n*length+3] - M_PI/2+2*M_PI;
        if (ori>M_PI)
            ori -= 2*M_PI;
        fprintf(stdout, "%f %f %f %f ", l[n*length]
                                      , l[n*length+1]
                                      , l[n*length+2]
                                      , ori);
        float* descr = &l[n*length+4];
        for(int i = 0; i < nhist; i++){
            for(int j = 0; j < nhist; j++){
                // conversion descriptor 
                int iA = j;
                int jA = nhist-1-i;
                for(int k = 0; k < nori; k++){
                    fprintf(stdout, "%i ", (int)descr[iA*nhist*nori+jA*nori+k]);
                }
            }
        }
        fprintf(stdout,"\n");
    }
    
}
