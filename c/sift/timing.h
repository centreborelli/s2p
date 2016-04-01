#include <time.h>

void portable_gettime(struct timespec *ts);
void print_elapsed_time(struct timespec *ts, char *title, int title_width);
