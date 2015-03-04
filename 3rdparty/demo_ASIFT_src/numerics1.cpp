// Authors: Unknown. Please, if you are the author of this file, or if you 
// know who are the authors of this file, let us know, so we can give the 
// adequate credits and/or get the adequate authorizations.


#include "numerics1.h"

float **allocate_float_matrix(int nrows, int ncols)
{
	float ** matrix;
	matrix = new float*[nrows];
	for(int i=0; i < nrows; i++) matrix[i] = new float[ncols];
	return matrix;
}


void desallocate_float_matrix(float **matrix, int nrows, int ncols)
{
	if (matrix == NULL) return;

	for(int i=0; i < nrows; i++) { delete[] matrix[i]; matrix[i] = 0;}
	matrix = 0;
	
	ncols = ncols; // to remove the warning unused parameter ¡®ncols¡¯
}


// **********************************************
//  LU based algorithms
// **********************************************



// Solves Ax=b by using lu decomposition 
int lusolve(float **a, float *x, float *b, int n)
{

	float d;
	int *indx = new int[n];

	if (ludcmp(a,n,indx,&d))
	{

		for(int i=0; i < n; i++) x[i] = b[i];
		lubksb(a,n,indx,x);
		
		delete[] indx; /*memcheck*/
		return 1;
	} else
	{

		printf("lusolve::lu decomposition failed\n");
		delete[] indx; /*memcheck*/
		return 0;
	}	
	delete[] indx; // Guoshen Yu
}



int ludcmp(float **a, int n, int *indx, float *d)
{

	int i,imax=0,j,k,aux;
	float big,dum,sum,temp;
	float *vv;


	vv=(float *) malloc(n*sizeof(float));
	*d=1.0;


	for(i=0;i<n;i++) indx[i]=i;


	/**** Look for the largest value of every line and store 1/this value****/
	for(i=0;i<n;i++){

		big=0.0;
		for(j=0;j<n;j++) 
			if ( (temp=fabs(a[i][j]))>big ) big=temp;

			if (big==0.0) { return 0; printf("LU Decomposition failed\n");}
			
		vv[i]=1.0/big;
	}


  

	for(j=0;j<n;j++){

		for(i=0;i<j;i++){
      
			sum=a[i][j];
			for(k=0;k<i;k++) sum-= a[i][k]*a[k][j];
			a[i][j]=sum; 

		}


		big=0.0;
		for(i=j;i<n;i++){
      
			sum=a[i][j];
			for(k=0;k<j;k++)   
				sum -=a[i][k]*a[k][j];
			a[i][j]=sum;

			if ( (dum=vv[i]*fabs(sum))>=big){
				big=dum;
				imax=i;
			}
		}

    
		if (j != imax){

			for(k=0;k<n;k++){

				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}

			*d=-(*d);
			vv[imax]=vv[j];

			aux=indx[j];
			indx[j]=indx[imax];
			indx[imax]=aux;
		}

		if (a[j][j]==0.0)  a[j][j]=NRTINY;
    

		if (j!=n-1 ){

			dum=1.0 / a[j][j];
			for(i=j+1;i<n;i++) a[i][j]*=dum;
		}
    

	}

	free(vv);
	return 1;
}



/* Solves the set of n linear equations Ax=b. Here a[0..n-1][0..n-1] as input, not as the matrix A but rather as its LU decomposition,*/
/* determined by the routine ludcmp. indx[0..n-1] is input as the permutation vector returned by ludcmp. b[0..n-1] is input as the    */
/* right hand side vector and returns with the solution vector x. */
void lubksb(float **a, int n, int *indx, float *b)
{

	int i,ii=0,j;
	float sum,*aux;

	aux= (float *) malloc(n*sizeof(float));
	
	
	for(i=0;i<n;i++)
		aux[i]=b[indx[i]];


	for(i=0;i<n;i++){
		sum=aux[i];

		if (ii)
			for(j=ii-1;j<i;j++) sum-=a[i][j]*aux[j];
		else if (sum) ii=i+1;
		aux[i]=sum;
	}

	for(i=n-1;i>=0;i--){
		sum=aux[i];
		for(j=i+1;j<n;j++) sum-=a[i][j]*aux[j];
		aux[i]=sum/a[i][i];
		b[i]=sum/a[i][i];
	}

	free(aux);
}