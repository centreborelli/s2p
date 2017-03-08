#ifndef _RPC_H
#define _RPC_H

// rational polynomial coefficient stuff

#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>


// rational polynomial coefficients (and those of the inverse model)
struct rpc {
	double numx[20];
	double denx[20];
	double numy[20];
	double deny[20];
	double scale[3], offset[3];

	double inumx[20];
	double idenx[20];
	double inumy[20];
	double ideny[20];
	double iscale[3], ioffset[3];

	double dmval[4];
	double imval[4];
};


// read an XML file specifying an RPC model
void read_rpc_file_xml(struct rpc *p, char *filename);

void print_rpc(FILE *f, struct rpc *p, char *n);

// evaluate the direct rpc model
void eval_rpc(double *result,
		struct rpc *p, double x, double y, double z);

// evaluate the inverse rpc model
void eval_rpci(double *result,
		struct rpc *p, double x, double y, double z);

// evaluate an epipolar correspondence
static void eval_rpc_pair(double xprime[2],
		struct rpc *a, struct rpc *b,
		double x, double y, double z);

// compute the height of a point given its location inside two images
double rpc_height(struct rpc *rpca, struct rpc *rpcb,
		double x, double y, double xp, double yp, double *outerr);

#endif // _RPC_H
