/**
 * File: main.c
 * Author: Cui Donghang
 * Last Update: 2022/04/24
 *
 * This is an implementation of the Distributed Pivoter.
 * It requires the library OpenMPI/MPICH and OpenMP to run this program.
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "graph.h"
#include "dkc.h"


int main(int argc, char **argv) {
	
	MPI_Init(NULL, NULL);

	dkc(argv[1], atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]));

	MPI_Finalize();

	return 0;
}
