#ifndef DECOMPOSITION_THREADS_H
#define DECOMPOSITION_THREADS_H
#include "constants.h"
#include "math.h"
#include "conio.h"

using namespace std;

void getMatrixWithoutRowAndCol_threads(double **matrix, int size, int row, int col, double **newMatrix);
double matrixDet_th(double **matrix, int size);

//аргумент для функции в поток
typedef struct funcArg{
	double **matrix;
	int size;
	double *det;
	int j;
	int end;
	int thread_num;
}  funcArg_t;

void *decomp_func(void *args){
	funcArg_t *arg = (funcArg_t*) args;
	int size = arg->size;
	int degree = 1;
	int j = arg->j;
	int end = arg->end;
	//cout << "Start thread # " << arg->thread_num  << endl;
	//cout << ", j = " << j << ", end = " << end << endl;

	//Матрица без строки и столбца
    double **newMatrix = new double*[size-1];
    for(int i = 0; i < size-1; i++) {
		newMatrix[i] = new double[size-1];
    }

	double tmp = 0.0;
	for(j; j < end; j++) {
		degree = pow((double)-1.0, (double)(j));
        getMatrixWithoutRowAndCol(arg->matrix, size, 0, j, newMatrix);
		tmp += (degree * arg->matrix[0][j] * matrixDet(newMatrix, size-1));
		*arg->det = tmp;
    }
	//cout << "Stop thread # " << arg->thread_num << endl;
	for(int i = 0; i < size-1; i++) {
		delete [] newMatrix[i];
	}
    delete [] newMatrix;
	return NULL;
}

/*
void getMatrixWithoutRowAndCol_threads(double **matrix, int size, int row, int col, double **newMatrix) {
    int offsetRow = 0; //Смещение индекса строки в матрице
    int offsetCol = 0; //Смещение индекса столбца в матрице
	int i, j;
	
    for(i = 0; i < size-1; i++) {
        //Пропустить row-ую строку
		//offsetRow = 1;
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
			//newMatrix[i][j] = 0;
        }
    }
}

//Вычисление определителя матрицы разложение по первой строке
double matrixDet_th(double **matrix, int size) {
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
        //Раскладываем по 0-ой строке, цикл бежит по столбцам
		for(int j = 0; j < size; j++) {
			degree = pow((double)-1.0, (double)(j));
            getMatrixWithoutRowAndCol(matrix, size, 0, j, newMatrix);
            det += (degree * matrix[0][j] * matrixDet_th(newMatrix, size-1));
        }
		for(int i = 0; i < size-1; i++) {
			delete [] newMatrix[i];
		}
		delete [] newMatrix;
    }
    return det;
}
*/

double start_decomposition_threads(double **matrix, int size) {
	double det = 0.0;
	//детерминант для каждого потока
    double dets[NUM_THREADS];
	for (int i=0; i < NUM_THREADS; i++){
		dets[i] = 0.0;
	}

	//потоки
	int status_addr[NUM_THREADS];
	pthread_t thread[NUM_THREADS];

	//аргументы для потоков
	funcArg_t args[NUM_THREADS];

    int degree = 1; // (-1)^(1+j) из формулы определителя
	
    if(size == 1) {
        return matrix[0][0];
    }
    else if(size == 2) {
        return matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0];
    }
    else {
		//заполнение содержимого аргументов
		//cout << "Size > 2 MAIN THREAD START" << endl;
		for(int th_i=0; th_i<NUM_THREADS; th_i++){
			args[th_i].matrix = matrix;
			args[th_i].size = size;
			//cout << th_i << " argument, size " << args[th_i].matrix.size() << " x " <<args[th_i].matrix[0].size() << endl;
			args[th_i].det = &dets[th_i];
			args[th_i].thread_num = th_i;
			if(th_i==0) args[th_i].j = 0;
			else args[th_i].j = (size/NUM_THREADS)*th_i;

			if(th_i==NUM_THREADS-1) args[th_i].end = size;
			else args[th_i].end = (size/NUM_THREADS)*(th_i+1);
		}
		//cout << "Threads creating ... " << endl;
		//создание потоков
		for(int th_i=0; th_i < NUM_THREADS; th_i++)
			pthread_create(&thread[th_i], NULL, decomp_func, (void*) &args[th_i]); 
		//cout << "... Done. Threads executing" << endl;
		for(int th_i=0; th_i < NUM_THREADS; th_i++)
			int status = pthread_join(thread[th_i], (void**)&status_addr[th_i]);
		
		//сумма
		//cout << "Add determinant" << endl;
		for(int th_i=0; th_i < NUM_THREADS; th_i++)
			det += *args[th_i].det;
    }
    return det;
}

int decomposition_threads(vector<double> M){
	cout << endl << "-----------------------------------------" << endl;
	cout << "decomposition_threads" << endl;
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
	
	det = start_decomposition_threads(Matrix, n);
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