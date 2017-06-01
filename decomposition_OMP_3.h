#ifndef DECOMPOSITION_OMP_3_H
#define DECOMPOSITION_OMP_3_H
#include "constants.h"
#include "math.h"

using namespace std;

double start_decomposition_OMP_3(double **matrix, int size) {
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
		
		//3-мерная матрица, для каждого потока своя 2-мерная матрица
		double ***newMatrix = new double**[NUM_THREADS];
		for(int i = 0; i < NUM_THREADS; i++) {
			newMatrix[i] = new double*[size-1];
		}
		for(int i = 0; i < NUM_THREADS; i++) {
			for(int j = 0; j < size-1; j++)
				newMatrix[i][j] = new double[size-1];
		}
		
		//для каждого потока свой детерминант
		double det_m[NUM_THREADS];
		for(int i=0; i<NUM_THREADS; ++i)
			det_m[i] = 0;
		
		int j;
		#pragma omp parallel for private(degree, j) shared(matrix, det_m, newMatrix) //schedule(dynamic) //firstprivate(newMatrix) //reduction(+:det)  //schedule(static)// работает правильно
		for(j = 0; j < size; j++) {
			int th_i = omp_get_thread_num(); //номер потока для индекса
			degree = pow((double)-1.0, (double)(j));
            getMatrixWithoutRowAndCol(matrix, size, 0, j, newMatrix[th_i]);
			det_m[th_i] += (degree * matrix[0][j] * matrixDet(newMatrix[th_i], size-1)); //создать массив с детерминантами, куда будут суммироваться, номер потока - номер ячейки
        }
		//подсчет итогового детерминанта
		for(int i=0; i<NUM_THREADS; ++i)
			det += det_m[i];

		//удаление 3-мерной матрицы
		for(int i = 0; i < NUM_THREADS; i++) {
			for(int j = 0; j < size-1; j++)
				delete [] newMatrix[i][j];
		}
		for(int i=0; i<NUM_THREADS; ++i)
			delete [] newMatrix[i];
		delete [] newMatrix;
    }
	
    return det;
}
int decomposition_OMP_3(vector<double> M){
	cout << endl << "-----------------------------------------" << endl;
	cout << "decompositions_OMP" << endl;
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

	omp_set_num_threads(NUM_THREADS);
	det = start_decomposition_OMP_3(Matrix, n);

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