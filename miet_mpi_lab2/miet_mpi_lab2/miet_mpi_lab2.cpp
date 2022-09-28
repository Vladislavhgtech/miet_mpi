// Starting script for manual call in unix systems
// ./start.sh

/*
Задание. Модифицировать программу, написанную на Л.Р. No1, так чтобы она работала на
основе коллективной передачи сообщений. Результаты работы сравнить и занести в отчет.
*/


#include <stdio.h>
#include <mpi.h>
#include <iostream>

using namespace std;

void DataInitialization(int* array, int size) {
	for (int i = 0; i < size; i++) {
		array[i] = 0;
	}
}

void PrintData(int* array, int size) {
	for (int i = 0; i < size; i++) {
		printf("message: %d\n", array[i]);
	}
}

int main(int* argc, char** argv) {
	int ProcNum, ProcRank, ProcPrevios = 0, Nring = 2, SendMessage = 0, * RecvNring;
	MPI_Init(argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	int DataSize = Nring * ProcNum + 1;
	RecvNring = new int[DataSize];
	if (ProcRank == 0) {
		DataInitialization(RecvNring, DataSize);
	}
	for (int i = 0; i < DataSize - 1; i++) {
		MPI_Bcast(RecvNring, DataSize, MPI_INT, ProcPrevios, MPI_COMM_WORLD);
		if ((ProcRank == ProcPrevios + 1) || (ProcPrevios + 1 == ProcNum && ProcRank == 0)) {
			RecvNring[ProcRank + ((i + 1) / ProcNum) * ProcNum] = RecvNring[ProcRank + ((i + 1) / ProcNum) * ProcNum - 1] + 1;
		}
		ProcPrevios = (ProcPrevios + 1 == ProcNum) ? 0 : ProcPrevios + 1;
		MPI_Barrier(MPI_COMM_WORLD);
	}
	MPI_Bcast(RecvNring, DataSize, MPI_INT, ProcPrevios, MPI_COMM_WORLD);
	printf("Proc: %d\n", ProcRank);
	PrintData(RecvNring, DataSize);
	MPI_Finalize();
	delete[] RecvNring;
}
