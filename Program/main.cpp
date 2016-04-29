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
double E, M, De[5],T, Values[5], J, T2, T_step, Results[6];
int size, MC_Cycle,e;
ofstream IsingData;
e=0;
IsingData.open("junk");


int abc=0;
////define initial conditions
InitialConditions(T,T2,T_step,J,MC_Cycle,size,abc);

IsingData<<"MC_Cycles "<<MC_Cycle<<endl;
IsingData<<"MC_Cycyle"<<" "<<"Temp"<<" "<<"Energy"<<" "<<"Specific_Heat"<<" "<<"Mag"<<" "<<"Sucept"<<" "<<"Sucept2"<<" "<<"Mag2"<<endl;


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
//for (int q=MC_Cycle; q<=16000; q+=200){
	//MC_Cycle=q;
///////////////////////////////


///////////////////initalize the particle lattice///////////////////
SpinLattice2(SpinMat, size, M, E);


//////////////create possible energy differences/////////////////////
for (dE=-2; dE <=2; dE++){De[dE+2]=exp(- J*4*(double) dE/T);}

////////////////////Preform Monte Carlo Tests///////////////////////

for (int i=0; i<MC_Cycle; i++){
	
Ising(SpinMat,size,M,E,De,e);
Values[0]+=E; Values[1]+=(E*E); Values[2]+=M; Values[3]+=M*M;
Values[4]+=fabs(M);

//cout << E << endl;
}
//calculate the average values for energy and such
for (int i=0; i<5; i++){
	Values[i] = Values[i]/(MC_Cycle+1);
}






ThermoQuantities(Values,Results,T);

IsingData<<MC_Cycle<<" "<<T<<" "<<Results[0]<<" "<<Results[1]<<" "<<Results[2]<<" "<<Results[3]<<" "<<Results[4]<<" "<<
	Results[5]<<endl;

	for (int i=0; i<5; i++){
		Values[i]=0;
	}
	M=0;
	E=0;

/////////////end of loop for changing number of MC cycles	
//}
//////////////////////////////////////////////


	cout << Values[0] <<" "<<Values[2]<<" "<<Values[1]-(Values[0]*Values[0])<<" "<<Values[4]<<endl;
cout << Results[0] << " "<<Results[1]<<" "<<Results[2]<<" "<<Results[4]<<endl;

cout <<"/////////////"<<endl;
for (int v=0; v<size; v++){
	for (int q=0; q<size; q++){
		cout << SpinMat[v][q]<<" ";
	}
	cout << endl;
}

cin.get();
}