/*
Задание. Модифицировать программу, написанную на Л.Р. №1 так чтобы она работала на
основе обобщенной передачи сообщений. Результаты работы сравнить (с результатами
полученными в Л.Р. №2) и занести в отчет.

*/

#include <stdio.h>
#include <mpi.h>
#include <iostream>

using namespace std;

void DataInitialization(int* array, int size, int fillNumber) {
	for (int i = 0; i < size; i++) {
		array[i] = fillNumber;
	}
}

void PrintData(int* array, int size) {
	for (int i = 0; i < size; i++) {
		printf("message: %d\n", array[i]);
	}
}

int main(int* argc, char** argv) {
	int ProcNum, ProcRank, ProcPrevios = 0, Nring = 1, * Msg, RecvMsg;
	MPI_Init(argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	Msg = new int[ProcNum];
	DataInitialization(Msg, ProcNum, 0);
	for (int i = 0; i < Nring * ProcNum - 1; i++) {
		MPI_Scatter(Msg, 1, MPI_INT, &RecvMsg, 1, MPI_INT, ProcPrevios, MPI_COMM_WORLD);
		if ((ProcRank == ProcPrevios + 1) || (ProcPrevios + 1 == ProcNum && ProcRank == 0)) {
			DataInitialization(Msg, ProcNum, RecvMsg + ProcRank);
		}
		ProcPrevios = (ProcPrevios + 1 == ProcNum) ? 0 : ProcPrevios + 1;
		MPI_Barrier(MPI_COMM_WORLD);
	}
	MPI_Scatter(Msg, 1, MPI_INT, &RecvMsg, 1, MPI_INT, ProcPrevios, MPI_COMM_WORLD);
	if (ProcRank == 0) {
		DataInitialization(Msg, ProcNum, RecvMsg + ProcRank);
		printf("Proc: %d\n", ProcRank);
		PrintData(Msg, ProcNum);
	}
	MPI_Finalize();
	delete[] Msg;
}