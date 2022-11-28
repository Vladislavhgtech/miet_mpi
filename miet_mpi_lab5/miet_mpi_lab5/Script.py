import numpy as np

CIRCLE1 = 2
CIRCLE2 = 3
SIZE = 16


def printMatrix(matrix):
    for i in range (SIZE):
        print(matrix[i])


if __name__ == '__main__':
    cm = []
    for i in range (SIZE):
        row = []
        for j in range(SIZE):
            row.append(complex(i % CIRCLE1 + j % CIRCLE2, i % CIRCLE2 + j % CIRCLE1))
        cm.append(row)
    printMatrix(cm)
    print("--------------------------------------------------------------")
    matrix = np.array(cm)
    printMatrix(matrix.dot(matrix))