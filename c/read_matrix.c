#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static void *freadwhole_f (FILE *f, long *on)
{
   int r, n = 0, nt = 0;
   char *t = NULL;
   while(1) {
      if (n >= nt)
      {
         nt = 1000+2 * (nt + 1);
         t = realloc(t, nt);
      }
      r = fgetc(f);
      if (r == EOF)
         break;
      t[n] = r;
      n += 1;
   }
   *on = n;
   return t;
}

int read_matrix_3x4(double H[3][4], char* name)
{
   // read the file
   long number_of_bytes;
   FILE* f = fopen(name,"r");
   char* tmpstr = (char*) freadwhole_f(f, &number_of_bytes);
   fclose(f);

   // all of this just to add a wimzy space in front of the string
   char str[number_of_bytes+2];
   strcpy(str," ");
   strncat(str, tmpstr, number_of_bytes+1);
   free(tmpstr);

   // read the 12 fileds of the matrix
   char* ff=str;
   int rr=0, r=0, nc=0;
   double* HH = (double*) H;
   while(rr++ < 12) {
      // %*[][ \n\t\r,:;] excludes: '[' ']' ':' ';' '.' 'space' '\t' '\n' '\r'
      // %n: number of read characters
      r += sscanf(ff,"%*[][ (){}\n\t\r,:;]%lf%n", &HH[r],&nc); ff += nc;
   }

   // fail if the number entries is not 12
   if(r!=12) {
      fprintf(stderr, "Expecting 12 numbers in \"%s\", %d found\n", name, r);
      abort();
      return 0;
   }

   return 1;
}

int read_matrix(double H[3][3], char* name)
{
   // read the file
   long number_of_bytes;
   FILE* f = fopen(name,"r");
   char* tmpstr = (char*) freadwhole_f(f, &number_of_bytes);
   fclose(f);

   // all of this just to add a wimzy space in front of the string
   char str[number_of_bytes+2];
   strcpy(str," ");
   strncat(str, tmpstr, number_of_bytes+1);
   free(tmpstr);

   // read the 9 fileds of the matrix
   char* ff=str;
   int rr=0, r=0, nc=0;
   double* HH = (double*) H;
   while(rr++ < 9) {
      // %*[][ \n\t\r,:;] excludes: '[' ']' ':' ';' '.' 'space' '\t' '\n' '\r'
      // %n: number of read characters
      r += sscanf(ff,"%*[][ (){}\n\t\r,:;]%lf%n", &HH[r],&nc); ff += nc;
   }

   // fail if the number entries is not 9
   if(r!=9) {
      fprintf(stderr, "Expecting 9 numbers in \"%s\", %d found\n", name, r);
      abort();
      return 0;
   }

   return 1;
}


#ifdef USE_TEST_MAIN
int main (int c, char **v)
{
   if (c<2) return 1;
   double H[3][3];
   read_matrix(H, v[1]);
   printf("%f %f %f\n %f %f %f\n %f %f %f\n",
         H[0][0], H[0][1],H[0][2],
         H[1][0], H[1][1],H[1][2],
         H[2][0], H[2][1],H[2][2]);
   return 0;
}
#endif
