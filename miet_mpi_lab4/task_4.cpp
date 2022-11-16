#include "mpi.h"
#include <iostream>
#include <queue>

using std::cout;
using std::endl;
using std::queue;
using std::min;

#define sizeMatrix 3
#define countMatrix 3

struct Complex {
    int real;
    int imaginary;
};
Complex** initializationMatrix(int n) {
    Complex* data = new Complex[n * n];
    Complex** array = new Complex * [n];
    for (int i = 0; i < n; i++) {
        array[i] = &(data[n * i]);
    }
    return array;
}
void fillZeroMatrix(Complex** array, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            array[i][j].real = 0;
            array[i][j].imaginary = 0;
        }
    }
}
void fillMatrix(Complex** array, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            array[i][j].real = i + j;
            array[i][j].imaginary = i + j;
        }
    }
}
void printMatrix(Complex** array, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << "(" << array[i][j].real << ", " << array[i][j].imaginary << "i" << ") ";
        }
        cout << endl;
    }
}

int main(int argc, char* argv[]) {

    int procNum, procRank;

    MPI_Init(&argc, &argv);
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

    //регистрируем тип
    MPI_Datatype complex_number;
    MPI_Datatype type[2] = { MPI_INT, MPI_INT };
    int blocklen[2] = { 1, 1 };
    MPI_Aint disp[2] = { 0, 1 };
    MPI_Type_create_struct(2, blocklen, disp, type, &complex_number);
    MPI_Type_commit(&complex_number);

    queue<Complex**> queueMatrix;
    bool isEnd;

    //Нулевой процесс
    if (procRank == 0) {
        //Инициализируем все матрицы и заполняем
        for (int i = 0; i < countMatrix; i++) {
            Complex** matrix;
            matrix = initializationMatrix(sizeMatrix);
            fillMatrix(matrix, sizeMatrix);
            queueMatrix.push(matrix);
        }

        //Пока не отсанется результируещий матрицы отправляем для перемножения
        while (queueMatrix.size() != 1) {
            isEnd = false;
            int size = (int)queueMatrix.size();

            for (int i = 0; i < min((int)(size / 2), procNum - 1); i++) {

                MPI_Send(&isEnd, 1, MPI_INT, i + 1, 0, MPI_COMM_WORLD);

                Complex** matrix1;
                matrix1 = initializationMatrix(sizeMatrix);
                matrix1 = queueMatrix.front();
                queueMatrix.pop();

                Complex** matrix2;
                matrix2 = initializationMatrix(sizeMatrix);
                matrix2 = queueMatrix.front();
                queueMatrix.pop();

                MPI_Send(&(matrix1[0][0]), 9, complex_number, i + 1, 0, MPI_COMM_WORLD);
                MPI_Send(&(matrix2[0][0]), 9, complex_number, i + 1, 0, MPI_COMM_WORLD);
            }

            queue<Complex**> tmpQueueMatrix;
            for (int i = 1; i <= min((int)(size / 2), procNum - 1); i++) {
                Complex** multiplyMatixResult;
                multiplyMatixResult = initializationMatrix(sizeMatrix);

                //Получаем результат перемножения
                MPI_Recv(&(multiplyMatixResult[0][0]), sizeMatrix * sizeMatrix, complex_number, i, 1, MPI_COMM_WORLD, &status);
                tmpQueueMatrix.push(multiplyMatixResult);
            }

            while (queueMatrix.size()) {
                tmpQueueMatrix.push(queueMatrix.front());
                queueMatrix.pop();
            }
            queueMatrix = tmpQueueMatrix;
        }

        for (int i = 1; i < procNum; i++) {
            isEnd = true;
            MPI_Send(&isEnd, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }

        //Выводим результат
        printMatrix(queueMatrix.front(), sizeMatrix);
    }

    //Все процессы кроме 0
    else {
        while (true) {

            MPI_Recv(&isEnd, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

            if (isEnd) {
                break;
            }

            Complex** matrix1;
            Complex** matrix2;
            Complex** multiplyMatixResult;
            matrix1 = initializationMatrix(sizeMatrix);
            matrix2 = initializationMatrix(sizeMatrix);
            multiplyMatixResult = initializationMatrix(sizeMatrix);

            MPI_Recv(&(matrix1[0][0]), sizeMatrix * sizeMatrix, complex_number, 0, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&(matrix2[0][0]), sizeMatrix * sizeMatrix, complex_number, 0, 0, MPI_COMM_WORLD, &status);

            fillZeroMatrix(multiplyMatixResult, sizeMatrix);

            //Перемножаем матрицы
            for (int i = 0; i < sizeMatrix; i++) {
                for (int j = 0; j < sizeMatrix; j++) {
                    for (int k = 0; k < sizeMatrix; k++)
                    {
                        multiplyMatixResult[i][j].real += matrix1[i][k].real * matrix2[k][j].real;
                        multiplyMatixResult[i][j].imaginary += matrix1[i][k].imaginary * matrix2[k][j].imaginary;
                    }
                }
            }
            //Отправляем результат перемножения 0 процессу
            MPI_Send(&(multiplyMatixResult[0][0]), sizeMatrix * sizeMatrix, complex_number, 0, 1, MPI_COMM_WORLD);
        }
    }

    //удалаем свой тип
    MPI_Type_free(&complex_number);

    MPI_Finalize();

    return 0;
}