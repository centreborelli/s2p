#include<stdio.h>
#include<stdlib.h>
#include<string.h>


struct list {
    char *current;
    struct list *next;
};

struct list *push(struct list *l, char *s)
{
    if (l == NULL) {
        l = (struct list *) malloc(sizeof(struct list));
        l->current = strdup(s);
        l->next = NULL;
    } else {
        l->next = push(l->next, s);
    }
    return l;
}


void list_print(struct list *l)
{
    if (l != NULL) {
        printf("%s\n", l->current);
        list_print(l->next);
    }
}


#ifdef MAIN_LISTS
int main(int c, char *v[])
{
	// process input arguments
	if (c != 1) {
		fprintf(stderr, "usage:\n\t"
				"echo strings | %s\n", *v);
		return 1;
	}

	// add each string from stdin to the list
    struct list *l = NULL;
	char s[FILENAME_MAX];

	while (fgets(s, FILENAME_MAX, stdin))
	{
		strtok(s, "\n");
        l = push(l, s);
	}

    // print list
    list_print(l);
    return 0;
}
#endif
