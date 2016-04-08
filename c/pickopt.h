// @c pointer to original argc
// @v pointer to original argv
// @o option name (after hyphen)
// @d default value
//
// char *oval = pick_option(&argc, &argv, "o", "37");
// returns "37" or the value of the option, removes 0 or 2 arguments
//
// bool o = pick_option(&argc, &argv, "o", NULL);
// returns NULL or true, removes 0 or 1 arguments
const char * pick_option(int *c, char ***v, const char *o, const char *d);
