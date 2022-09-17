#include <stdio.h>
#include <mpi.h>
#include <iostream>

using namespace std;

int main(int* argc, char** argv) {
	int ProcNum, ProcRank, Nring = 2, SendMessage = 0, RecvNring;
	MPI_Init(argc, &argv);
	MPI_Status Status;
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	if (ProcRank == 0) {
		MPI_Send(&SendMessage, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
	}
	for (int i = 0; i < Nring; i++) {
		MPI_Recv(&RecvNring, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
		RecvNring++;
		if (ProcRank + 1 != ProcNum) {
			MPI_Send(&RecvNring, 1, MPI_INT, ProcRank + 1, 0, MPI_COMM_WORLD);
		}
		else {
			MPI_Send(&RecvNring, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
	}
	if (ProcRank == 1) {
		MPI_Recv(&RecvNring, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
		RecvNring++;
	}
	printf("PROC: %d message: %d\n", ProcRank, RecvNring);
	MPI_Finalize();
}