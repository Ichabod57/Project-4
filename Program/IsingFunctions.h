#include <iostream>
#include <cmath>
#include <fstream>
#include <random>
#include <ctime>
using namespace std;

inline int periodic(int i, int limit, int add){
	return (i+limit+add) % (limit);
}

void SpinLattice1(double **, int, double &, double &);
void SpinLattice2(double **, int, double &, double &);
void Equate(double **, double **, int);
void Ising(double **, int, double &, double &, double *, int &);
void InitialConditions(double &, double &, double &, double &, int &, int &, int &);
void ThermoQuantities(double *, double *, double);
void Probability(double *, int, double, int );