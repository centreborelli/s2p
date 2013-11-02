#include <stdio.h>
#include <stdbool.h>

/* Original QNM's version from E. Meinhardt-Llopis
 * Megawave port by G. Facciolo*/

/*/ a smart parameter is just like a regular parameter, but it can be
// re-defined at the shell-environment.  Instead of
// 
// #define NUMBER 42
// ...
// printf("%g", NUMBER);
// 
// do
// SMART_PARAMETER(NUMBER,42)
// ...
// printf("%g", NUMBER());
// 
// Notice that the environment only gets queried once, at the first use.
// 
*/
#ifndef VERBOSE_SMART_PARAMETER
#define VERBOSE_SMART_PARAMETER 0
#endif

#define SMART_PARAMETER(n,v) static double n(void)\
{\
   static int smapa_known_ ## n = false;\
   static double smapa_value_ ## n = v;\
   if (!smapa_known_ ## n)\
   {\
      int r;\
      char *sv;\
      double y;\
      sv = getenv(#n);\
      if(VERBOSE_SMART_PARAMETER) fprintf(stderr,"scanning the environment for \"%s\"... ", #n);\
      if (sv)\
         r = sscanf(sv, "%lf", &y);\
      if (sv && r == 1)\
      {\
         if(VERBOSE_SMART_PARAMETER) fprintf(stderr, "got value %g\n", y);\
         smapa_value_ ## n = y;\
      } else {\
         if(VERBOSE_SMART_PARAMETER) fprintf(stderr, "kept default value %g\n",\
               smapa_value_ ## n);\
      }\
      smapa_known_ ## n = true;\
   }\
   return smapa_value_ ## n;\
}


#define SMART_PARAMETER_INT(n,v) static int n(void)\
{\
   static int smapa_known_ ## n = false;\
   static int smapa_value_ ## n = v;\
   if (!smapa_known_ ## n)\
   {\
      int r;\
      int y;\
      char *sv;\
      sv = getenv(#n);\
      if(VERBOSE_SMART_PARAMETER)fprintf(stderr,"scanning the environment for \"%s\"... ", #n);\
      if (sv)\
         r = sscanf(sv, "%d", &y);\
      if (sv && r == 1)\
      {\
         if(VERBOSE_SMART_PARAMETER)fprintf(stderr, "got int value %d\n", y);\
         smapa_value_ ## n = y;\
      } else {\
         if(VERBOSE_SMART_PARAMETER)fprintf(stderr, "kept default int value %d\n",\
               smapa_value_ ## n);\
      }\
      smapa_known_ ## n = true;\
   }\
   return smapa_value_ ## n;\
}


