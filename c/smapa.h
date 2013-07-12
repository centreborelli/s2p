
// a smart parameter is just like a regular parameter, but it can be
// re-defined at the shell-environment.  Instead of
//
// 	#define NUMBER 42
// 	...
// 	printf("%g", NUMBER);
//
// do
// 	SMART_PARAMETER(NUMBER,42)
// 	...
// 	printf("%g", NUMBER());
//
// Notice that the environment only gets queried once, at the first use.
//
#define SMART_PARAMETER(n,v) static double n(void)\
{\
	static bool smapa_known_ ## n = false;\
	static double smapa_value_ ## n = v;\
	if (!smapa_known_ ## n)\
	{\
		fprintf(stderr,"scanning the environment for \"%s\"... ", #n);\
		int r;\
		char *sv = getenv(#n);\
		double y;\
		if (sv)\
			r = sscanf(sv, "%lf", &y);\
		if (sv && r == 1)\
		{\
			fprintf(stderr, "got value %g\n", y);\
			smapa_value_ ## n = y;\
		} else {\
			fprintf(stderr, "kept default value %g\n",\
					smapa_value_ ## n);\
		}\
		smapa_known_ ## n = true;\
	}\
	return smapa_value_ ## n;\
}

#define SMART_PARAMETER_SILENT(n,v) static double n(void)\
{\
	static bool smapa_known_ ## n = false;\
	static double smapa_value_ ## n = v;\
	if (!smapa_known_ ## n)\
	{\
		int r;\
		char *sv = getenv(#n);\
		double y;\
		if (sv)\
			r = sscanf(sv, "%lf", &y);\
		if (sv && r == 1)\
			smapa_value_ ## n = y;\
		smapa_known_ ## n = true;\
	}\
	return smapa_value_ ## n;\
}
