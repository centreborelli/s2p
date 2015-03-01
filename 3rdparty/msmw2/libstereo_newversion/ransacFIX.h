#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


int DEBUG() {
return 0;
}


#include "vvector.h"

	/*   dim=3
    * | a x + b y + c    -    d |
    * */
float test_model(int dim, float *data, float *model)  {
	return fabs (model[0]*data[0] + model[1]*data[1] + model[2]   -   data[2]);
}


float generate_model ( int dim, int numel, float *data, int *indices, int n, float *result) {

   double A[3][3]={{0,0,0},{0,0,0},{0,0,0}};
   double b[3]={0,0,0};

//   assert(n==3);
//   for(int t=0;t<3;t++) {
//      int i = indices[t];
////      printf("%d\t", i);
//
//      assert (i<numel);
//      A[t][0] = data[dim*i+0];
//      A[t][1] = data[dim*i+1];
//      A[t][2] = 1;
//      b[t]    = data[dim*i+2];
////      printf("%f %f %f    %f\n", A[t][0],A[t][1],A[t][2],b[t]);
//   }
////      printf("\n");
      
   for(int t=0;t<n;t++) {
      int ii = indices[t];
      for (int i=0;i<3;i++) {
         double d[4] = {data[dim*ii+0], data[dim*ii+1], 1, data[dim*ii+2]};
         for (int j=0;j<3;j++) {
            A[i][j] += d[i]*d[j];
         }
         b[i] += d[i]*d[3];
      }
   }

   double iA[3][3];
   double detA;
   INVERT_3X3(iA, detA, A);

   if(fabs(detA)>0.00001) {
      MAT_DOT_VEC_3X3(result, iA,b);
   } else {
		result[0]=0;
		result[1]=0;
		result[2]=INFINITY;
	}

//   for(int i=0;i<3;i++){
//      printf("%f \n", test_model(dim, data, result));
//   }
//   printf("%f %f %f\n", result[0], result[1], result[2]);

}




// /RANDOM UTILS
void swapi(int *a, int *b) {
   int swap=*a;
   *a=*b;
   *b=swap;
}

static int randombounds(int a, int b, unsigned int *seed)
{
   if (b < a) {
      printf(" wrong use of randombounds\n");
      exit(1);
   }
   if (b == a) return b;
   return a + rand_r(seed)%(b - a + 1);
}
//fisher-yates
void shuffle(int *t, int n)
{
	unsigned int seed = time(NULL);
   for (int i = 0; i < n-1; i++)
      swapi(&t[i], &t[randombounds(i, n-1, &seed)]);
}





// dim=3; // dimension of the problem
// numel=1000; // number of elements in data
//	n=4;   // estimate the model with n random samples of data
//	k=100; // maximum ransac iterations 
//	t=0.5; // error threshold
//	d=8;   // consider only models that verify at least d data points 
int ransacFIX(int dim, int numel, float *data, int n, int k, float t, int d, float *best_model, int *num_best_consensus, int *best_consensus_set_idx, float *best_error, int num_FIXED, int *FIXED_inliers_idx) 
{
	float better_model[dim];
	float maybe_model[dim];

	*best_error=INFINITY;
	*num_best_consensus=0;

   int random_set_idx[numel];

   // reset random_set_idx
   for(int i=0;i<numel;i++){
      random_set_idx[i] = i;
   }
   
   // set the fixed data points
   assert(num_FIXED<n);
   for(int i=0;i<num_FIXED;i++){
      swapi(&random_set_idx[i], &random_set_idx[FIXED_inliers_idx[i]]);
   }

	unsigned int seed = time(NULL);
	for (int iterations=0; iterations<k; iterations++) {
		if(DEBUG())	printf("%d %f\n", iterations, *best_error);

      // complete shuffle of the remaining indices //fisher-yates
      // this can be done several times on the randomized index vector
      for(int i=num_FIXED;i<n;i++){
         swapi(&random_set_idx[i], &random_set_idx[randombounds(i, numel-1,&seed)]);
      }

		/* estimate the model from the maybe inliers fit */
      // CHECK THAT num_inliers is sufficient to estimate any model
      generate_model ( dim, numel, data, random_set_idx, n, maybe_model);
		if(DEBUG())	printf("%.2f %.2f %.2f\n", maybe_model[0], maybe_model[1], maybe_model[2]);

		/* expand consensus set */
	   int consensus_set_idx[numel];
      int num_consensus = 0;
      float this_error = 0;
		for(int i=0;i<numel;i++) {
		   float err = test_model(dim, &(data[i*dim]), maybe_model); 
		   if(DEBUG())	printf("%f\t", err);

		   if(err < t ) { 
		   	consensus_set_idx[num_consensus]=i; 
		   	num_consensus++;
            this_error += err; // accumulate error 
		   } 
		}
      this_error /= num_consensus;
		if(DEBUG())	printf("num_consensus %d \n",num_consensus);


		/* solution candidate */
		if(num_consensus > d ) {
//       // fitting of a better model using the consensus set
//			generate_model ( dim, numel, data, consensus_set_idx, num_consensus, better_model );
//       this_error = 0;
//       num_consensus =0;
//		   for(int i=0;i<numel;i++) {
//		      float err = test_model(dim, &(data[i*dim]), better_model); 
//          if(err < t) {
//             this_error+=err;
//             num_consensus++;
//          }
//          this_error/=num_consensus;
//       }


			//if(  num_consensus  >= *num_best_consensus )  // use this if the number of inliers is critic
         {
            if (this_error < *best_error) {
   				*best_error = this_error;
               // copy the model
               for(int p=0;p<dim;p++)    best_model[p] = maybe_model[p];  
               // copy the consensus set
   				*num_best_consensus=num_consensus;
               if ( best_consensus_set_idx ) {
                  for(int i=0;i<num_consensus;i++) best_consensus_set_idx[i] = consensus_set_idx[i];
               }
   				if(DEBUG())	printf("iteration %d error %f num %d\n",iterations, *best_error, *num_best_consensus);
            }
			}
      }

	}

	return(*num_best_consensus!=0);
}



#if 0

test_ransac() {
	int dim=3;
	int numel = 25;
	double model[3];
	double data[numel*dim];	
	int set[numel];
	int num_set;
	double err;
	int k,d,n,i;
	double t;
	double x,y,z;

	/* generate data */
	srand(time(NULL));
	double mm[]=  {1,1,0};
	for(i=0;i<10;i++) {
		x = rand() % 20;
		y = rand() % 20;
		z = x*mm[0] + y*mm[1] + mm[2];
		data[i*dim+0]=x;
		data[i*dim+1]=y;
		data[i*dim+2]=z+rand()/10000000000.;
		printf("%f ", z+rand()/10000000000.);
	}	

	mm[1]=0;
	for(i=10;i<25;i++) {
		x = rand() % 20;
		y = rand() % 20;
		z = x*mm[0] + y*mm[1] + mm[2];
		data[i*dim+0]=x;
		data[i*dim+1]=y;
		data[i*dim+2]=z+rand()/10000000000.;
	}	


	n=5;   // estimate the model with 5 elements
	k=100; // maximum iterations 
	t=0.1; // error threshold variance?
	d=8;   // accept a model with 8 elements within the band
	/* call ransac */
	for(i=1;i<10;i++) {
		if (ransac(dim, numel, data, n,k,t,d,model, &num_set, set, &err)) 
			/* display result */
			printf("%d) %f %f %f  %f,%d\n", i, model[0], model[1], model[2], err, num_set);
		else { 
			printf("RANSAC RETOURNED 0 , MUST RE-RUN \n");
			break;
		}
	}

}


int main() {
	test_ransac();
	test_ransac2();
}
#endif
