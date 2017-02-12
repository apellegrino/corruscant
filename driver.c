#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include "kdtree.h"
#include "mpi.h"

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

#define MASTER 0

struct kdtree tree_construct(int size, double x[], double y[], double z[]) {

	struct kdtree t;
	t.size = size;
	t.x = x;
	t.y = y;
	t.z = z;

	// Argsort the inputs
	int *x_arg, *y_arg, *z_arg;
	x_arg = argsort(x, size);
	y_arg = argsort(y, size);
	z_arg = argsort(z, size);

	t.root = build(x,y,z,x_arg,y_arg,z_arg,0,size-1,X);
	//verify(t.root,0);

	return t;
}

long long two_point_correlation(struct kdtree tree, double x[], double y[],
									double z[], int n, double r, MPI_Comm comm) {
    int mpi_size, mpi_rank;
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);


	int i;
	long long result;
	result = 0;
	for (i=mpi_rank; i<n; i += mpi_size) {
		result += radius(tree.root, 0, x[i], y[i], z[i], r);
	}

    MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_INTEGER, MPI_SUM, comm);

	return result;
}

double landy_szalay(struct kdtree data, struct kdtree random, double rad, MPI_Comm comm) {


    int mpi_size, mpi_rank;
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);

	int total;
	int dd, dr, rr;
	dd = 0;
	dr = 0;
	rr = 0;
	double f;
	f = (double)random.size/(double)data.size; //N_random / N_data

	printf("rsize = %d, dsize = %d, rad = %2.20f\n",random.size,data.size,rad);

	int i;
	for (i=mpi_rank; i<data.size; i += mpi_size) {
		dd += radius(data.root,0,data.x[i],data.y[i],data.z[i],rad);
		dr += radius(data.root,0,random.x[i],random.y[i],random.z[i],rad);
	}
	for (i=mpi_rank; i<random.size; i += mpi_size) {
		rr += radius(random.root,0,random.x[i],random.y[i],random.z[i],rad);
	}

    MPI_Allreduce(MPI_IN_PLACE, &dd, 1, MPI_INTEGER, MPI_SUM, comm);
    MPI_Allreduce(MPI_IN_PLACE, &dr, 1, MPI_INTEGER, MPI_SUM, comm);
    MPI_Allreduce(MPI_IN_PLACE, &rr, 1, MPI_INTEGER, MPI_SUM, comm);

	double ls;
	//if(rr == 0) return 0.0;

	ls = (f*f*(double)dd - 2*f*(double)dr + (double) rr)/(double)rr;
	printf("For r = %f:\n", rad);
	printf("DD = %d\n", dd);
	printf("DR = %d\n", dr);
	printf("RR = %d\n", rr);
	printf("Landy-Szalay = %f\n", ls);

	/*
	// correct for self-counting as well as double pair counting
	//dd = (dd - data.size)/2;
	//rr = (rr - random.size)/2;
	dd = dd/2;
	rr = rr/2;
	ls = (f*f*(double)dd - 2*f*(double)dr + (double) rr)/(double)rr;
	
	printf("For r = %f:\n", rad);
	printf("DD = %d\n", dd);
	printf("DR = %d\n", dr);
	printf("RR = %d\n", rr);
	printf("Landy-Szalay = %f\n", ls);
	*/

	/*
	destroy(data_t.root);
	destroy(rand_t.root);
	free(x);
	free(y);
	free(z);
	free(data_x);
	free(data_y);
	free(data_z);
	*/
	return ls;
}
