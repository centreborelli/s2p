#ifndef _ABSTRACT_DSF_H
#define _ABSTRACT_DSF_H


// implementation of a DSF with path compression but without union by rank

//void adsf_assert_consistency(int *t, int n);
void adsf_begin(int *t, int n);
int adsf_union(int *t, int n, int a, int b);
int adsf_find(int *t, int n, int a);
int adsf_number_of_classes(int *t, int n);


// implementation of a DSF with path compression and union by rank

//void adsfr_assert_consistency(int *t, int *r, int n);
//void adsfr_begin(int *t, int *r, int n);
//int adsfr_union(int *t, int *r, int n, int a, int b);
//int adsfr_find(int *t, int *r, int n, int a);
//int adsfr_number_of_classes(int *t, int *r, int n);



#endif /* _ABSTRACT_DSF_H */
