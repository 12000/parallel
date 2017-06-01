#include "decomposition.h"
#include "decomposition_threads.h"
#include "decomposition_OMP_3.h"
#include "constants.h"

using namespace std;
int main()
{	
	//генерация данных от 0 до 9
	/*ofstream out("m.txt");
	srand( time(0) );
	for(int i=0; i < n*n; ++i){
		int num = rand()%10;
		out << num << endl;
	}*/

	vector<double> M(n*n);
	ifstream in("m.txt");
	int num = 0;
	for(int i=0; i < n*n; ++i){
			in >> num;
			M[i] = (double)num;
	}

	cout <<"Start" << endl << "n = " << n << endl << "NUM_THREADS " << NUM_THREADS << endl << "-----------------------------------" << endl;

	
	decomposition(M);
	decomposition_threads(M);
	decomposition_OMP_3(M);
	
    system("pause");
    
    return 0;
}