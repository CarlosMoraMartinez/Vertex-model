#include <cmath>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cstdlib>
#include <random>
#include <chrono> 


using namespace std::chrono; 

const int RANDOM_SEED = 1234;
const int N_ITERS =10000000;
const int SHORT_LEN = 6;
const int LONG_LEN = 1000;
const int SHORT_REPS = 10;
typedef std::map<int, double> mapa;

double sumArray(double *shortarr, int len){
	double a = 0;
	for(int i = 0; i<len; i++) a+= pow(shortarr[i], 2) - shortarr[i] + pow(shortarr[i], 3);
	return a;
}

inline double sumArrayInline(double *shortarr, int len){
	double a = 0;
	for(int i = 0; i<len; i++) a+= pow(shortarr[i], 2) - shortarr[i] + pow(shortarr[i], 3);
	return a;
}

double sumArraySimple(double *shortarr, int len){
	double a = 0;
	for(int i = 0; i<len; i++) a+= shortarr[i];
	return a;
}

inline double sumArraySimpleInline(double *shortarr, int len){
	double a = 0;
	for(int i = 0; i<len; i++) a+= shortarr[i];
	return a;
}


void testVectorInit(){
	std::vector<double> vec;
	for(int i = 0; i < LONG_LEN; i++){
		vec.push_back(i);
	}
	std::cout << "Final number: " << vec[LONG_LEN - 1] << std::endl;
}

void testArrayInit(){
	double *longarr1, *longarr2;
	longarr1 = new double[1];
	longarr1[0] = 0.0;
	for(int i = 2; i <= LONG_LEN; i++){
		longarr2 = new double [i];
		for(int j = 0; j < i - 1; j++){
			longarr2[j] = longarr1[j];
		}
		longarr2[i - 1] = i*0.1;
		longarr1 = longarr2;
	}	
	std::cout << "Final number: " << longarr1[LONG_LEN - 1] << std::endl;
}


double distance(double x1, double x2, double y1, double y2){
	return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}
inline double distanceInline(double x1, double x2, double y1, double y2){
	return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}


int main(int argc, char *argv[]){


	srand(RANDOM_SEED);
	std::default_random_engine generator;
	std::uniform_real_distribution<double> unif;
	

	double a;
	int d;
	double x1 = 0.3245, x2 = 543.55, y1= 55442, y2 = -0.4325;
	double shortarr[6]= { 16, 2, 77, 40, 12071, 435 };
	std::vector<double> vec{ 16, 2, 77, 40, 12071, 435 };
	mapa res;
	res.insert(std::pair<int, double >(0, 16));
	res.insert(std::pair<int, double >(1, 2));
	res.insert(std::pair<int, double >(2, 77));
	res.insert(std::pair<int, double >(3, 40));
	res.insert(std::pair<int, double >(4, 12071));
	res.insert(std::pair<int, double >(5, 435));	
	

	//Test random generators
	std::cout << "Testing random generators - unif()" << std::endl;

	auto start = high_resolution_clock::now(); 
	for(int i = 0; i < N_ITERS; i++){
		a = unif(generator);
	}
	auto stop = high_resolution_clock::now();
	auto duration1 = duration_cast<microseconds>(stop - start); 



	std::cout << "Testing random generators - rand()" << std::endl;
	start = high_resolution_clock::now(); 
	for(int i = 0; i < N_ITERS; i++){
		d = std::rand();
	}
	stop = high_resolution_clock::now();
	auto duration2 = duration_cast<microseconds>(stop - start); 



	std::cout << "Testing sum in vector" << std::endl;	
	//Test array sum
	start = high_resolution_clock::now(); 
	for(int i = 0; i < N_ITERS; i++){
		a = 0;
		for(int j = 0; j<SHORT_LEN; j++) a+= pow(vec[j], 2) - vec[j] + pow(vec[j], 3);
	}
	stop = high_resolution_clock::now();
	auto duration3 = duration_cast<microseconds>(stop - start); 
	std::cout << "Result: " << a << std::endl;



	std::cout << "Testing sum in array" << std::endl;	
	//Test array sum
	start = high_resolution_clock::now();
	for(int i = 0; i < N_ITERS; i++){
		a = 0;
		for(int j = 0; j<SHORT_LEN; j++) a+= pow(shortarr[j], 2) - shortarr[j] + pow(shortarr[j], 3);
	}
	stop = high_resolution_clock::now();
	auto duration4 = duration_cast<microseconds>(stop - start); 
	std::cout << "Result: " << a << std::endl;



	std::cout << "Testing sum in map" << std::endl;
	//Test map sum
	start = high_resolution_clock::now(); 
	for(int i = 0; i < N_ITERS; i++){
		a = 0;
		for(int j = 0; j<SHORT_LEN; j++) a+= pow(res[j], 2) - res[j] + pow(res[j], 3);
	}
	stop = high_resolution_clock::now();
	auto duration5 = duration_cast<microseconds>(stop - start); 
	std::cout << "Result: " << a << std::endl;



	std::cout << "Testing sum in array - function" << std::endl;	
	//Test array sum in function
	start = high_resolution_clock::now(); 
	for(int i = 0; i < N_ITERS; i++){
		a = sumArray(shortarr, SHORT_LEN);
	}
	stop = high_resolution_clock::now();
	auto duration6 = duration_cast<microseconds>(stop - start); 
	std::cout << "Result: " << a << std::endl;


	std::cout << "Testing sum in array - inline function" << std::endl;	
	//Test array sum in Inline function
	start = high_resolution_clock::now(); 
	for(int i = 0; i < N_ITERS; i++){
		a = sumArrayInline(shortarr, SHORT_LEN);
	}
	stop = high_resolution_clock::now();
	auto duration7 = duration_cast<microseconds>(stop - start); 
	std::cout << "Result: " << a << std::endl;



	std::cout << "Testing realocating vector" << std::endl;	
	//Test realocating to vector
	start = high_resolution_clock::now(); 
	for(int i = 0; i < SHORT_REPS; i++){
		testVectorInit();
	}
	stop = high_resolution_clock::now();
	auto duration8 = duration_cast<microseconds>(stop - start); 


	std::cout << "Testing realocating array" << std::endl;		
	//Test realocating to array
	start = high_resolution_clock::now(); 
	for(int i = 0; i < SHORT_REPS; i++){
		testArrayInit();
	}
	stop = high_resolution_clock::now();
	auto duration9 = duration_cast<microseconds>(stop - start); 
	
	
	std::cout << "Testing mathematical function (distance) - main" << std::endl;	
	//Test array sum in Inline function
	start = high_resolution_clock::now(); 
	for(int i = 0; i < N_ITERS; i++){
		a = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
	}
	stop = high_resolution_clock::now();
	auto duration10 = duration_cast<microseconds>(stop - start); 
	std::cout << "Result: " << a << std::endl;


	std::cout << "Testing mathematical function (distance) - function" << std::endl;	
	//Test array sum in Inline function
	start = high_resolution_clock::now(); 
	for(int i = 0; i < N_ITERS; i++){
		a = distance(x1, x2, y1, y2);
	}
	stop = high_resolution_clock::now();
	auto duration11 = duration_cast<microseconds>(stop - start); 
	std::cout << "Result: " << a << std::endl;


	std::cout << "Testing mathematical function (distance) - function Inline" << std::endl;	
	//Test array sum in Inline function
	start = high_resolution_clock::now(); 
	for(int i = 0; i < N_ITERS; i++){
		a = distanceInline(x1, x2, y1, y2);
	}
	stop = high_resolution_clock::now();
	auto duration12 = duration_cast<microseconds>(stop - start); 
	std::cout << "Result: " << a << std::endl;




	//Now testing sum array again but with a simpler function

	std::cout << "Testing sum in vector" << std::endl;	
	//Test array sum
	start = high_resolution_clock::now(); 
	for(int i = 0; i < N_ITERS; i++){
		a = 0;
		for(int j = 0; j<SHORT_LEN; j++) a+= vec[j];
	}
	stop = high_resolution_clock::now();
	auto duration13 = duration_cast<microseconds>(stop - start); 
	std::cout << "Result: " << a << std::endl;



	std::cout << "Testing sum in array" << std::endl;	
	//Test array sum
	start = high_resolution_clock::now();
	for(int i = 0; i < N_ITERS; i++){
		a = 0;
		for(int j = 0; j<SHORT_LEN; j++) a+= shortarr[j];
	}
	stop = high_resolution_clock::now();
	auto duration14 = duration_cast<microseconds>(stop - start); 
	std::cout << "Result: " << a << std::endl;



	std::cout << "Testing sum in map" << std::endl;
	//Test map sum
	start = high_resolution_clock::now(); 
	for(int i = 0; i < N_ITERS; i++){
		a = 0;
		for(int j = 0; j<SHORT_LEN; j++) a+= res[j];
	}
	stop = high_resolution_clock::now();
	auto duration15 = duration_cast<microseconds>(stop - start); 
	std::cout << "Result: " << a << std::endl;



	std::cout << "Testing sum in array - function" << std::endl;	
	//Test array sum in function
	start = high_resolution_clock::now(); 
	for(int i = 0; i < N_ITERS; i++){
		a = sumArraySimple(shortarr, SHORT_LEN);
	}
	stop = high_resolution_clock::now();
	auto duration16 = duration_cast<microseconds>(stop - start); 
	std::cout << "Result: " << a << std::endl;


	std::cout << "Testing sum in array - inline function" << std::endl;	
	//Test array sum in Inline function
	start = high_resolution_clock::now(); 
	for(int i = 0; i < N_ITERS; i++){
		a = sumArraySimpleInline(shortarr, SHORT_LEN);
	}
	stop = high_resolution_clock::now();
	auto duration17 = duration_cast<microseconds>(stop - start); 
	std::cout << "Result: " << a << std::endl;







	std::cout << "Finished" << std::endl << "***" << std::endl;
	std::cout << "gen unif():" << duration1.count()/static_cast<double>(N_ITERS) << std::endl;
	std::cout << "gen rand():" << duration2.count()/static_cast<double>(N_ITERS) << std::endl;
	std::cout << "sum Vector '+=x² - x + x³':" << duration3.count()/static_cast<double>(N_ITERS) << std::endl;
	std::cout << "sum Array '+=x² - x + x³':" << duration4.count()/static_cast<double>(N_ITERS) << std::endl;
	std::cout << "sum Map() '+=x² - x + x³':" << duration5.count()/static_cast<double>(N_ITERS) << std::endl;
	std::cout << "sum Array func() '+=x² - x + x³':" << duration6.count()/static_cast<double>(N_ITERS) << std::endl;
	std::cout << "sum Array inline() '+=x² - x + x³':" << duration7.count()/static_cast<double>(N_ITERS) << std::endl;

	std::cout << "sum Vector '+=x'" << duration13.count()/static_cast<double>(N_ITERS) << std::endl;
	std::cout << "sum Array '+=x':" << duration14.count()/static_cast<double>(N_ITERS) << std::endl;
	std::cout << "sum Map() '+=x':" << duration15.count()/static_cast<double>(N_ITERS) << std::endl;
	std::cout << "sum Array func() '+=x':" << duration16.count()/static_cast<double>(N_ITERS) << std::endl;
	std::cout << "sum Array inline() '+=x':" << duration17.count()/static_cast<double>(N_ITERS) << std::endl;


	std::cout << "alloc Vector():" << duration8.count()/static_cast<double>(SHORT_REPS) << std::endl;
	std::cout << "alloc Array():" << duration9.count()/static_cast<double>(SHORT_REPS) << std::endl;
	std::cout << "distance main:" << duration10.count()/static_cast<double>(N_ITERS) << std::endl;
	std::cout << "distance func:" << duration11.count()/static_cast<double>(N_ITERS) << std::endl;
	std::cout << "distance func inline:" << duration12.count()/static_cast<double>(N_ITERS) << std::endl;
	std::cout << "***" << std::endl;

	exit(0);
	
}
