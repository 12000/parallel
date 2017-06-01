#ifndef DECOMPOSITION_H
#define DECOMPOSITION_H
#include "constants.h"
#include "math.h"

using namespace std;
//Возвращает матрицу matrix без row-ой строки и col-того столбца, результат в newMatrix
void getMatrixWithoutRowAndCol(double **matrix, int size, int row, int col, double **newMatrix) {
    int offsetRow = 0; //Смещение индекса строки в матрице
    int offsetCol = 0; //Смещение индекса столбца в матрице
	int i, j;
	
	//#pragma omp parallel for private(i, j, offsetRow, offsetCol) shared (newMatrix, matrix) //работает, но делает хуже по времени
    for( i = 0; i < size-1; i++) {
        //Пропустить row-ую строку
        if(i == row) {
            offsetRow = 1; //Как только встретили строку, которую надо пропустить, делаем смещение для исходной матрицы
        }
 
        offsetCol = 0; //Обнулить смещение столбца
        for(j = 0; j < size-1; j++) {
            //Пропустить col-ый столбец
            if(j == col) {
                offsetCol = 1; //Встретили нужный столбец, проускаем его смещением
            }
 
            newMatrix[i][j] = matrix[i + offsetRow][j + offsetCol];
        }
    }
}

//разложение по первой строке
double matrixDet(double **matrix, int size) {
    double det = 0.0;
    int degree = 1;
 
    //Условие выхода из рекурсии
    if(size == 1) {
        return matrix[0][0];
    }
    //Условие выхода из рекурсии
    else if(size == 2) {
        return matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0];
    }
    else {
        //Матрица без строки и столбца
        double **newMatrix = new double*[size-1];
        for(int i = 0; i < size-1; i++) {
            newMatrix[i] = new double[size-1];
        }
		
		for(int j = 0; j < size; j++) {
			degree = pow((double)-1.0, (double)(j));
            getMatrixWithoutRowAndCol(matrix, size, 0, j, newMatrix);
            det += (degree * matrix[0][j] * matrixDet(newMatrix, size-1));			
        }

        for(int i = 0; i < size-1; i++) {
            delete [] newMatrix[i];
        }
        delete [] newMatrix;
    }
	
    return det;
}
int decomposition(vector<double> M){
	cout << endl << "-----------------------------------------" << endl;
	cout << "decompositions" << endl;
	cout << endl;
	double det = 1;

	double **Matrix = new double* [n];
	for(int j=0; j<n; ++j){
			Matrix[j] = new double [n];
	}

	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			Matrix[i][j] = M[i*n + j];
		}
	}

	__int64 start, end, tps;
	QueryPerformanceFrequency((LARGE_INTEGER *)&tps);
	QueryPerformanceCounter((LARGE_INTEGER *)&start);

	det = matrixDet(Matrix, n);
	QueryPerformanceCounter((LARGE_INTEGER *)&end);
    cout << "Determinant: " << det << endl;
	double time = ((double)(end - start) / tps) * 100.;
	cout << "Time:" << time << " ms" << endl;
        for(int i = 0; i < n; i++) {
            delete [] Matrix[i];
        }
        delete [] Matrix;
	return (int) floor(time + 0.5);
}



#endif