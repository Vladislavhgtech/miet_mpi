#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include "mpi.h"

const int CIRCLE1 = 2;
const int CIRCLE2 = 3;
const int SIZE = 16;
const int EXIT_SIZE = 4;

struct ComplexMatrix
{
	int real[16][16];
	int imaginary[16][16];
	int size;
};

void printMatrix(ComplexMatrix cm) {
	for (int i = 0; i < cm.size; i++)
	{
		for (int j = 0; j < cm.size; j++)
		{
			printf("%d+%di ", cm.real[i][j], cm.imaginary[i][j]);
		}
		printf("\n");
	}
}

ComplexMatrix Sum(ComplexMatrix cm1, ComplexMatrix cm2)
{
	struct ComplexMatrix ans;
	ans.size = cm1.size;

	for (int i = 0; i < ans.size; i++)
	{
		for (int j = 0; j < ans.size; j++)
		{
			ans.real[i][j] = cm1.real[i][j] + cm2.real[i][j];
			ans.imaginary[i][j] = cm1.imaginary[i][j] + cm2.imaginary[i][j];
		}
	}

	return ans;
};

ComplexMatrix Substract(ComplexMatrix cm1, ComplexMatrix cm2)
{
	struct ComplexMatrix ans;
	ans.size = cm1.size;

	for (int i = 0; i < ans.size; i++)
	{
		for (int j = 0; j < ans.size; j++)
		{

			ans.real[i][j] = cm1.real[i][j] - cm2.real[i][j];
			ans.imaginary[i][j] = cm1.imaginary[i][j] - cm2.imaginary[i][j];
		}
	}

	return ans;
};

ComplexMatrix Multiply(ComplexMatrix cm1, ComplexMatrix cm2)
{
	struct ComplexMatrix ans;
	ans.size = cm1.size;
	for (int i = 0; i < ans.size; i++) {
		for (int j = 0; j < ans.size; j++)
		{
			ans.real[i][j] = 0;
			ans.imaginary[i][j] = 0;
		}
	}

	if (ans.size <= EXIT_SIZE)
	{
		for (int i = 0; i < ans.size; i++)
		{
			for (int j = 0; j < ans.size; j++)
			{
				for (int r = 0; r < ans.size; r++)
				{
					ans.real[i][j] += cm1.real[i][r] * cm2.real[r][j];
					ans.real[i][j] += -1 * cm1.imaginary[i][r] * cm2.imaginary[r][j];
					ans.imaginary[i][j] += cm1.real[i][r] * cm2.imaginary[r][j];
					ans.imaginary[i][j] += cm1.imaginary[i][r] * cm2.real[r][j];
				}
			}
		}
		return ans;
	}
	else
	{
		ComplexMatrix** a = new ComplexMatrix * [2];
		ComplexMatrix** b = new ComplexMatrix * [2];
		ComplexMatrix** c = new ComplexMatrix * [2];
		for (int i = 0; i < 2; i++) {
			a[i] = new ComplexMatrix[2];
			b[i] = new ComplexMatrix[2];
			c[i] = new ComplexMatrix[2];
			for (int j = 0; j < 2; j++) {
				a[i][j].size = cm1.size / 2;
				b[i][j].size = cm1.size / 2;
				c[i][j].size = cm1.size / 2;
				for (int k = 0; k < a[i][j].size; k++) {
					for (int l = 0; l < a[i][j].size; l++) {
						a[i][j].real[k][l] = cm1.real[i * cm1.size / 2 + k][j * cm1.size / 2 + l];
						b[i][j].real[k][l] = cm2.real[i * cm2.size / 2 + k][j * cm2.size / 2 + l];
						a[i][j].imaginary[k][l] = cm1.imaginary[i * cm1.size / 2 + k][j * cm1.size / 2 + l];
						b[i][j].imaginary[k][l] = cm2.imaginary[i * cm2.size / 2 + k][j * cm2.size / 2 + l];
					}
				}
			}
		}

		ComplexMatrix p1 = Multiply(Sum(a[0][0], a[1][1]), Sum(b[0][0], b[1][1]));
		ComplexMatrix p2 = Multiply(Sum(a[1][0], a[1][1]), b[0][0]);
		ComplexMatrix p3 = Multiply(a[0][0], Substract(b[0][1], b[1][1]));
		ComplexMatrix p4 = Multiply(a[1][1], Substract(b[1][0], b[0][0]));
		ComplexMatrix p5 = Multiply(Sum(a[0][0], a[0][1]), b[1][1]);
		ComplexMatrix p6 = Multiply(Substract(a[1][0], a[0][0]), Sum(b[0][0], b[0][1]));
		ComplexMatrix p7 = Multiply(Substract(a[0][1], a[1][1]), Sum(b[1][0], b[1][1]));

		c[0][0] = Substract(Sum(Sum(p1, p4), p7), p5);
		c[0][1] = Sum(p3, p5);
		c[1][0] = Sum(p2, p4);
		c[1][1] = Substract(Sum(Sum(p1, p3), p6), p2);

		for (int i = 0; i < c[0][0].size; i++) {
			for (int j = 0; j < c[0][0].size; j++) {
				ans.real[i][j] = c[0][0].real[i][j];
				ans.imaginary[i][j] = c[0][0].imaginary[i][j];
				ans.real[i + ans.size / 2][j] = c[1][0].real[i][j];
				ans.imaginary[i + ans.size / 2][j] = c[1][0].imaginary[i][j];
				ans.real[i][j + ans.size / 2] = c[0][1].real[i][j];
				ans.imaginary[i][j + ans.size / 2] = c[0][1].imaginary[i][j];
				ans.real[i + ans.size / 2][j + ans.size / 2] = c[1][1].real[i][j];
				ans.imaginary[i + ans.size / 2][j + ans.size / 2] = c[1][1].imaginary[i][j];
			}
		}

		return ans;
	}

};

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);

	MPI_Datatype MPI_COMPLEX_MATRIX;
	MPI_Type_contiguous(513, MPI_INT, &MPI_COMPLEX_MATRIX);
	MPI_Type_commit(&MPI_COMPLEX_MATRIX);

	int ProcNum, ProcRank;
	MPI_Status Status;
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	MPI_Group WorldGr1;

	MPI_Group ZeroToOne;
	MPI_Comm ZtoComm;

	MPI_Comm_group(MPI_COMM_WORLD, &WorldGr1);

	int* ztoRank = new int[3] { 0, 1, 2 };

	MPI_Group_incl(WorldGr1, 3, ztoRank, &ZeroToOne);
	MPI_Comm_create(MPI_COMM_WORLD, ZeroToOne, &ZtoComm);

	ComplexMatrix complexMatrix;
	complexMatrix.size = SIZE;

	if (ProcRank == 0) {
		for (int i = 0; i < complexMatrix.size; i++)
		{
			for (int j = 0; j < complexMatrix.size; j++)
			{
				complexMatrix.real[i][j] = i % CIRCLE1 + j % CIRCLE2;
				complexMatrix.imaginary[i][j] = i % CIRCLE2 + j % CIRCLE1;
			}
		}
		MPI_Send(&complexMatrix, 1, MPI_COMPLEX_MATRIX, 1, 1, ZtoComm);
	}

	if (ProcRank == 1) {
		MPI_Recv(&complexMatrix, 1, MPI_COMPLEX_MATRIX, 0, 1, ZtoComm, &Status);
		printf("data:\n");
		printMatrix(complexMatrix);
		complexMatrix = Multiply(complexMatrix, complexMatrix);
		MPI_Send(&complexMatrix, 1, MPI_COMPLEX_MATRIX, 2, 0, ZtoComm);
	}

	if (ProcRank == 2) {
		MPI_Recv(&complexMatrix, 1, MPI_COMPLEX_MATRIX, 1, 0, ZtoComm, &Status);
		printf("result:\n");
		printMatrix(complexMatrix);
	}

	MPI_Finalize();
	return 0;
}