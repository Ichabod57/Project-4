
#include "IsingFunctions.h"



///////////////////create the initial state for the Ising model////////////////////
void SpinLattice1(double **SpinMat, int size, double& M, double &E){
	srand(time(NULL));
int v;

//fill the spin matrix randomly with 
for (int i=0; i<size; i++){
	for (int j=0; j<size; j++){
		SpinMat[i][j]=-1;
		M += SpinMat[i][j];
		
	} 
}

for (int i=0; i<size; i++){
	for (int j=0; j<size; j++){
		E -= SpinMat[i][j]*(SpinMat[periodic(i,size,-1)][j] + SpinMat[i][periodic(j,size,-1)]);
	} 
}


}


void SpinLattice2(double **SpinMat, int size, double &M, double &E){
	srand(time(NULL));
int v;

//fill the spin matrix randomly with 
for (int i=0; i<size; i++){
	for (int j=0; j<size; j++){
		
		
		v = rand() % 100;
		if (v<50){
			SpinMat[i][j]=1;
		}
		else{
			SpinMat[i][j] = -1;
		}
       M += SpinMat[i][j];
	} 
}

for (int i=0; i<size; i++){
	for (int j=0; j<size; j++){
		E -= SpinMat[i][j]*(SpinMat[periodic(i,size,-1)][j] + SpinMat[i][periodic(j,size,-1)]);
	} 
}

}

////////////////set two matrices equal to eachother////////////
void Equate(double **Mat1, double **Mat2, int size){
	for (int i=0; i<size; i++){
		for (int j=0; j<size; j++){
			Mat1[i][j]=Mat2[i][j];
		} 
	}
}




///////////////Metropolis to solve Ising Model/////////////////
void Ising(double **SpinMat, int size, double &M, double &E, double De[5], int &e){
	
	int DeltaE=0;
	int a, b, c, d;
	int x, y;
	double r=0;
	
	//////////generate a random place in the matrix//////////////
	for (int i=0; i<size*size; i++){
        x=rand() % size;
        y=rand() % size;
		a=SpinMat[x][periodic(y,size,-1)];
		b=SpinMat[periodic(x,size,-1)][y];
		c=SpinMat[x][periodic(y,size,1)];
		d=SpinMat[periodic(x,size,1)][y];
		//cout <<a+b+c+d<<" "<<a<<" "<<b<<" "<<c<<" "<<d<<endl;
	DeltaE = 2*SpinMat[x][y]*(a+b+c+d);
	
	r = ((double) rand() / (RAND_MAX));

	DeltaE=DeltaE/4;
	
		
	////////Metropolis testing//////////
	
	   if (r <= De[DeltaE+2]){
		SpinMat[x][y] *= -1;
	   M += 2*SpinMat[x][y];
	   E += DeltaE*4;
	   e++;
	  //cout << E << endl;
	  }
	
}
}



//////////////////define initial conditions for the system////////////////////
void InitialConditions(double & Start_Temp, double & End_Temp, double & Temp_Step, double & J, int & MC_Cycles, int & size, int & Start){

cout <<"How many particles do you want in the lattice row? ";
cin >> size;
cin.get();

cout <<"How many Monte-Carlo Cycles do you want to do? ";
cin >> MC_Cycles;
cin.get();

cout <<"What is the constant J? ";
cin >> J;
cin.get();

cout <<"What is the Starting Temperature? ";
cin >> Start_Temp;
cin.get();

cout <<"What is the Ending Temperature? ";
cin >> End_Temp;
cin.get();

cout <<"What is the Step size of the temperature? ";
cin >> Temp_Step;
cin.get();

cout <<"After how many MC Cycles is the solution stable? ";
cin >> Start;
cin.get();
}


/////////////////calculate all relevent thermodynamic quantities///////////////////

void ThermoQuantities(double * Values,double * Results, double T){
	//Values[0]=E, Values[1]=E*E, Values[2]=M, Values[3]=M*M, Values[4]=|M|
    //Results[0]=E, Results[1]=Specific Heat, Results[2]=|M|, Results[3]=Susceptability
	//Results[4]=Susceptability2

	Results[0]=Values[0];
	Results[1]=(1/(T*T))*(Values[1]-(Values[0]*Values[0]));
	Results[2]=Values[4];
	Results[3]=(1/(T))*(Values[3]-(Values[2]*Values[2]));
	Results[4]=(1/(T))*(Values[3]-(Values[4]*Values[4]));
	
	Results[5]=Values[2];
}


////////////////Probability/////////////////
void Probability(double *Energies, int MC_Cycles, double Energy, int Start){
	
	int EnergyRange[2][600];
	int En = (int) Energy;

	/////////Fill the array energy range with + to - 200 around the last energy value in the mc cycle//////
	for (int i=0; i<600; i++){
		EnergyRange[1][i]=0;
		EnergyRange[0][i]=En-300 + i;
		//cout << EnergyRange[0][i] << endl;
	}


	////////count number of times eneries appear///////////
	for (int i=Start; i<MC_Cycles; i++){
		for (int j=0; j<600; j++){
          if (EnergyRange[0][j]==(int) Energies[i]){
	           EnergyRange[1][j] = EnergyRange[1][j] + 1;
           }

		}
	}



	ofstream Probability;
	Probability.open("junk");
	for (int i=0; i<600; i++){
		if (EnergyRange[1][i] !=0){
		Probability<<EnergyRange[0][i]<<" "<<EnergyRange[1][i]<<" "<<(double) EnergyRange[1][i]/(MC_Cycles-Start)<<endl;
	}
	}
}



