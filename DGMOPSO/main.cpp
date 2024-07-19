
/*
//---------------- Author: Francesco Castellini* and Annalisa Riccardi** ---------------//
//------- Affiliation: -----------------------------------------------------------------//
//--------------* Dipartimento di Ingegneria Aerospaziale, Politecnico di Milano -------//
//--------------**Centrum for Industrial Mathematics, University of Bremen -------------//
//--------E-mail: castellfr@gmail.com, nina1983@gmail.com ------------------------------//

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

#include "pso.h"
#include "dgmopso.h"

USING_PART_OF_NAMESPACE_EIGEN
using namespace std;

//SINGLE OBJECTIVE OPTIMIZATION PROBLEMS
Eigen::VectorXd rastragin(Eigen::VectorXd& x);
Eigen::VectorXd schwefel(Eigen::VectorXd& x);
Eigen::VectorXd griewangk(Eigen::VectorXd& x);
Eigen::VectorXd griewangk_constrained(Eigen::VectorXd& x);
Eigen::VectorXd f_1disolated(Eigen::VectorXd& x);
Eigen::VectorXd f_1disolated2(Eigen::VectorXd& x);

//BI-OBJECTIVE OPTIMIZATION PROBLEMS
Eigen::VectorXd t1(Eigen::VectorXd& x);
Eigen::VectorXd t2(Eigen::VectorXd& x);
Eigen::VectorXd t3(Eigen::VectorXd& x);
Eigen::VectorXd t4(Eigen::VectorXd& x);
Eigen::VectorXd t5(Eigen::VectorXd& x);
Eigen::VectorXd t6(Eigen::VectorXd& x);
Eigen::VectorXd t7(Eigen::VectorXd& x);
Eigen::VectorXd t3c1(Eigen::VectorXd& x);
Eigen::VectorXd t3c2(Eigen::VectorXd& x);
Eigen::VectorXd t3c3(Eigen::VectorXd& x);


//Example main to show a single-objective PSO run and a multi-objective DGMOPSO run on mathematical test problems
//The model must be a function taking as input an Eigen VectorXd of nvar elements with the values of the optimization 
//variables and returning an Eigen VectorXd of nobj+ncstr elements with the evaluated optimization objectives and constraints
//All optimization options and functionalities are commented below
//NOTE: Problems griewangk and t1 are tested, with a Sleep(10 ms) for each evalation to show parallel computing speed-up
void main ()
{
	clock_t start, end;	double CPUtime;

	//PSO test: Griewangk function (multimodal)
	int seed=time(0);	//initialize random seed with current clock time
	srand(seed); cout << "Random seed: " << seed << endl;		
	start = clock();	//initialize cpu time count
	int nvar=2;
	int ncstr=0;
	int niter=100;
	int niterRevComm=0;	//ask for user action every niterRevComm iterations. Do not ask for user actions if =0
	int nparticles=100;
    Eigen::VectorXd LOWERBOUND, UPPERBOUND;
	LOWERBOUND = Eigen::VectorXd::Constant(nvar,-600);
	UPPERBOUND = Eigen::VectorXd::Constant(nvar,600);
	PSO pso(nvar,ncstr,niter,nparticles,LOWERBOUND,UPPERBOUND);
	pso.setModel(griewangk);
	pso.setProblemName("griewangk");	//for output files naming
	pso.setPrintLevel(1);				//amount of information printed to screen, from 0 to 3
	pso.setMultiThread(false);			//activate multi-threading --> NOTE: at the moment only multiThreadType=1 is implemented
	pso.setMultiThreadType(1);			//1: use OpenMP implementation for shared-memory machines (only if multiThread=true)
	pso.setNIterReverseComm(niterRevComm);
	pso.optimize();
	end=clock(); CPUtime=((double)(end - start))/CLOCKS_PER_SEC;
	cout << "PSO CPU time: " << CPUtime << " seconds." << endl;

	//DGMOPSO test: 
	seed=time(0);		//initialize random seed with current clock time
	srand(seed); cout << "Random seed: " << seed << endl;		
	start = clock();	//initialize cpu time count
	int nobj=2;
	nvar=10;
	ncstr=0;
	niter=100;
	niterRevComm=0;					//ask for user action every niterRevComm iterations. Do not ask for user actions if =0
	nparticles=100;
	int archivesize=100;
	LOWERBOUND = Eigen::VectorXd::Constant(nvar,-1);
	UPPERBOUND = Eigen::VectorXd::Constant(nvar,1);
	DGMOPSO dgmopso(nvar,nobj,ncstr,niter,nparticles,archivesize,LOWERBOUND,UPPERBOUND);
	dgmopso.setModel(t1);
	dgmopso.setProblemName("t1");	//for output files naming
	dgmopso.setPrintLevel(1);		//amount of information printed to screen, from 0 to 3
	dgmopso.setMultiThread(false);	//activate multi-threading --> NOTE: at the moment only multiThreadType=1 is implemented
	dgmopso.setMultiThreadType(1);	//1: use OpenMP implementation for shared-memory machines (only if multiThread=true)
	dgmopso.setNIterReverseComm(niterRevComm);
	dgmopso.optimize();
	end=clock(); CPUtime=((double)(end - start))/CLOCKS_PER_SEC;
	cout << "DGMOPSO CPU time: " << CPUtime << " seconds." << endl;
}