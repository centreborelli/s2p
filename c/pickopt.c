#include <string.h>
#include <stdbool.h>

// @c pointer to original argc
// @v pointer to original argv
// @o option name (after hyphen)
// @d default value
char * pick_option(int *c, char ***v, char *o, char *d)
{
	int argc = *c;
	char **argv = *v;
	int id = d ? 1 : 0;
	for (int i = 0; i < argc - id; i++)
		if (argv[i][0] == '-' && 0 == strcmp(argv[i]+1, o))
		{
			char *r = argv[i+id]+1-id;
			*c -= id+1;
			for (int j = i; j < argc - id; j++)
				(*v)[j] = (*v)[j+id+1];
			return r;
		}
	return d;
}

// char *oval = pick_option(&argc, &argv, "o", "37");
// returns "37" or the value of the option, removes 0 or 2 arguments
//
// bool o = pick_option(&argc, &argv, "o", NULL);
// returns NULL or true, removes 0 or 1 arguments
//

#ifdef MAIN_PICKOPT
#include <stdio.h>

static void print_args(int c, char **v)
{
	for (int i = 0; i <= c; i++)
		printf("ARG[%d/%d] = \"%s\"\n", i, c, v[i]);
}

int main(int c, char **v)
{
	printf("arguments before processing:\n");
	print_args(c, v);
//	char *o = pick_option(&c, &v, "o", "42");
//	printf("o = \"%s\"\n", o);
	printf("pick_option 'o' is: \"%d\"\n", (bool) pick_option(&c, &v, "o", NULL));
	printf("arguments after processing:\n");
	print_args(c, v);
	return 0;
}
#endif//MAIN_PICKOPT
