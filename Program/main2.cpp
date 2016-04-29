#include <iostream>
#include <cmath>
#include <fstream>
#include <random>
#include <ctime>
#include <cstdlib>
using namespace std;
#include "IsingFunctions.h"

void main(){
	srand(time(NULL));


///////////////define variables//////////////////
double E, M, De[5],T, Values[5], J, T2, T_step, Results[6], *Energies;
int size, MC_Cycle, accept, Start;
accept=0;
ofstream IsingData;
IsingData.open("80x80_FineScale2.txt");



////define initial conditions
InitialConditions(T,T2,T_step,J,MC_Cycle,size, Start);

IsingData<<"MC_Cycles "<<MC_Cycle<<endl;
IsingData<<"MC_Cycyle"<<" "<<"accept"<<" "<<"Temp"<<" "<<"Energy"<<" "<<"Specific_Heat"<<" "<<"Mag"<<" "<<"Sucept"<<" "<<"Sucept2"<<" "<<"Mag2"<<endl;


for (int i=0; i<5; i++){
	Values[i]=0;}
int dE;
M=0; E=0;
double **SpinMat;

SpinMat = new double *[size];
for (int i=0; i<5; i++){
	SpinMat[i]=0;
}

for (int i = 0; i<size; i++){
	SpinMat[i] = new double [size];
}

////////////adjusting the number of MC cycles to detemine stability///////////////
//for (int q=MC_Cycle; q<=50000; q+=1000){
	//MC_Cycle=q;
	Energies = new double [MC_Cycle];
///////////////////////////////


///////////////////initalize the particle lattice///////////////////
SpinLattice2(SpinMat, size, M, E);

int u=0;
while (T<T2){
	u=u+1;
	cout<<u<<endl;
//////////////create possible energy differences/////////////////////
for (dE=-2; dE <=2; dE++){De[dE+2]=exp(- J*4*(double) dE/T);}

////////////////////Preform Monte Carlo Tests///////////////////////

for (int i=0; i<MC_Cycle; i++){
	
Ising(SpinMat,size,M,E,De, accept);

if (i>=Start){
Values[0]+=E; Values[1]+=(E*E); Values[2]+=M; Values[3]+=M*M;
Values[4]+=fabs(M);
}
Energies[i]=E;


}
//////calculate the average values for energy and such///////////////
for (int i=0; i<5; i++){
	Values[i] = Values[i]/(MC_Cycle-Start);
}
///////////////calculate the thermodynamic quantities of interest///////////
ThermoQuantities(Values,Results,T);


/////////////calculate the probability of having a certain energy////////////////
//Probability(Energies,MC_Cycle,Values[0], Start);


IsingData<<MC_Cycle<<" "<<accept <<" "<<T<<" "<<Results[0]/(size*size)<<" "<<Results[1]/(size*size)<<" "<<Results[2]/(size*size)<<" "<<Results[3]/(size*size)<<" "<<Results[4]/(size*size)<<" "<<
	Results[5]/(size*size)<<endl;


////////////end of temp changing/////////////
T=T+T_step;
}
///////////////////////////////////////


	for (int i=0; i<5; i++){
		Values[i]=0;
	}
	M=0;
	E=0;

/////////////end of loop for changing number of MC cycles
	accept=0;
//}
//////////////////////////////////////////////


	cout << Values[0] <<" "<<Values[2]<<" "<<Values[1]-(Values[0]*Values[0])<<" "<<Values[4]<<endl;
cout << Results[0] << " "<<Results[1]<<" "<<Results[2]<<" "<<Results[3]<<endl;

cout <<"/////////////"<<endl;
for (int v=0; v<size; v++){
	for (int q=0; q<size; q++){
		cout << SpinMat[v][q]<<" ";
	}
	cout << endl;
}

cin.get();
}