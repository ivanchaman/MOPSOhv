
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

#include "optimizer.h"

#include <conio.h>
#include <Windows.h>
#include <Winbase.h>

#ifdef _DEBUG
#ifdef WIN32
int omp_get_thread_num() { return 0; };
int omp_get_num_threads() { return 1; };
#endif
#endif

USING_PART_OF_NAMESPACE_EIGEN
using namespace std;

Optimizer::Optimizer(int NVAR, int NCONS, Eigen::VectorXd LOWERBOUND, Eigen::VectorXd UPPERBOUND, Eigen::VectorXd GLOWERBOUND, Eigen::VectorXd GUPPERBOUND): m_bestFeasible(NVAR){
	m_name = "O";
	m_numberModelEval = 0;
	m_NVAR = NVAR;
	m_NCONS = NCONS;
	m_LOWERBOUND = LOWERBOUND;
	m_UPPERBOUND = UPPERBOUND;
	m_GLOWERBOUND = GLOWERBOUND;
	m_GUPPERBOUND = GUPPERBOUND;
	opt_print_level = 2;
	m_feasibleOnly = false;	
	m_nIterReverseComm = 0;
	m_multiThread = false;	
	m_multiThreadType = 0;	
	m_initialization = 0;
	m_indexInit = 0;
	m_filenameInit = "";
}

Optimizer::~Optimizer() {
	// TODO Auto-generated destructor stub
}

void Optimizer::setModel(Eigen::VectorXd (*model)(Eigen::VectorXd& x)){
	m_model = model;
}

void Optimizer::setProblemName(string name){
	m_name = name;
}

void Optimizer::setPrintLevel(int i){
	opt_print_level = i;
}

void Optimizer::setFeasibleOnly(){
	m_feasibleOnly = true;
}

double Optimizer::constrViolation(double constrValue, int constrIndex){
	double violation = 0;

	if(constrValue>m_GUPPERBOUND(constrIndex))
			violation += constrValue - m_GUPPERBOUND(constrIndex);
	if(constrValue<m_GLOWERBOUND(constrIndex))
			violation += m_GLOWERBOUND(constrIndex) - constrValue;

	return violation;
}

//This function loads from file to matrix the optimization solutions to Eigen matrix of nsol x nvar double variables in [0;1])
//the procedure reads the number of lines in the file nsol (number of solutions loaded)
// but NOTE: no check is performed on the input file except for all variables being in [0;1] --> all lines MUST have 
//				the same number of columns, equal to nvar (one column per each variable)
//INPUTS:
//	nvar: number of variables in each solution (number of optimization variables)
//	name: string with the name of the input file containing the solutions matrix, including extension. 
//			NOTE: The file will be searched in data/optimization folder and all its subfolders
//OUTPUTS:
//	sol: Eigen matrix of nsol x nvar elements containing the requested solutions (nsol rows for the nsol solutions in the file)
//	flagread: true if the file has been found and the loaded solution is correct (i.e. all variables in [0;1])
bool Optimizer::loadFirstGuessSolutionFromFile(int nvar, string name, Eigen::MatrixXd& sol)
{
	//Look for input file containing first guess solution(s) in all optimization subfolders: 
	bool flagread=false; string filename; ifstream filestream;
	if(!flagread) { filename="optimization/"+name; 
	filestream.open(filename.data(), ios::in); if(filestream.good()) flagread=true; else filestream.clear(); }
	if(!flagread) { filename="optimization/DGMOPSO/"+name;
	filestream.open(filename.data(), ios::in); if(filestream.good()) flagread=true; else filestream.clear(); }
	if(!flagread) { filename="optimization/HybridGO/"+name;
	filestream.open(filename.data(), ios::in); if(filestream.good()) flagread=true; else filestream.clear(); }
	if(!flagread) { filename="optimization/MADS/"+name;
	filestream.open(filename.data(), ios::in); if(filestream.good()) flagread=true; else filestream.clear(); }
	if(!flagread) { filename="optimization/MOACOr/"+name;
	filestream.open(filename.data(), ios::in); if(filestream.good()) flagread=true; else filestream.clear(); }
	if(!flagread) { filename="optimization/NSGA2/"+name;
	filestream.open(filename.data(), ios::in); if(filestream.good()) flagread=true; else filestream.clear(); }
	if(!flagread) { filename="optimization/PSO/"+name;
	filestream.open(filename.data(), ios::in); if(filestream.good()) flagread=true; else filestream.clear(); }
	if(!flagread) { filename="optimization/WORHP/"+name;
	filestream.open(filename.data(), ios::in); if(filestream.good()) flagread=true; else filestream.clear(); }
	if(!flagread)
		cout << "Error opening file " << name << ", no such file found in any of the optimization subfolders" << endl;
	else		//if the given file has been found, read first guess solution data
		{
		int nsol = std::count(istreambuf_iterator<char>(filestream), istreambuf_iterator<char>(),'\n');
		filestream.close();
		filestream.open(filename.data(), ios::in);
		sol=Eigen::MatrixXd(nsol,nvar);
		//Read sol-th solution from solutions matrix in the file:
		for(int j=0;j<nsol;j++)
			{
			if(!flagread)
				break;
			for(int k=0;k<nvar;k++)
				if(filestream.eof())
					{
					flagread=false;
					cout << "Error loading from file variable "<<k<<" of solution "<<j<<": end of file reached" << endl;
					break;
					}
				else
					{
					filestream >> sol(j,k);				//load k-th variable of j-th solution
					if(sol(j,k)<0 || sol(j,k)>1)		//check if k-th variable of j-th solution is within [0;1]
						{
						flagread=false;
						cout << "Error loading from file variable "<<k<<" of solution "<<j<<": variable value "<<sol(j,k)<<" is not within [0;1]" << endl;
						break;
						}
					}
			}
		}
	return flagread;
}



// Function used by PSO and DGMOPSO to perform the model evaluation for all swarm particles in a given iteration, also 
// implementing the parallel computing specifications and printing outputs to video if requested by the opt_print_level
//INPUTS:
//  * iternumber: number of current iteration, just for screen printing
//  * useOpenMP: boolean for using OpenMP parallel implementation for shared memory machines 
//			--> NOTE: OMP_NUM_THREADS environment variables MUST be set to the number of processors you want to use!
//  * useMPI:    boolean for using MPI parallel implementation for distributed memory machines 
//  * useOpenCL: boolean for using OpenCL GPU computing implementation for graphic cards (NVIDIA and ATI supported)
//  * swarm: matrix of (nparticles x nvars) elements with the current position of the swarm, to be evalauted
//OUTPUTS:
//  * Jswarm: matrix of (nparticles x nobj) elements with the computed objective function values for the current swarm positions
//  * Cswarm: matrix of (nparticles x ncstrs) elements with the computed constraints values for the current swarm positions
//	* quitOptimization: false --> continue optimization with following iterations; true --> quit optimization upon user's request
bool Optimizer::evaluateParticles(const int& iternumber, const bool& useOpenMP, const bool& useMPI, const bool& useOpenCL, int swarmsize, double **swarm, double **Jswarm, double **Cswarm)
{
	if((useOpenMP+useOpenCL+useMPI)>1)
		{ cerr << "Attempting to use more than one parallel computing paradigm in PSO::evaluateParticles(), exiting program." << endl; }	//at most one of the three parallel implementations can be used
	Eigen::VectorXd var(m_NVAR), objconst(m_NOBJ+m_NCONS);
	int chunk=1; double penalty=0;
	#pragma omp parallel default(shared) private(var,objconst,penalty) if(useOpenMP)
		{
		var.resize(m_NVAR); objconst.resize(1+m_NCONS);
		int nthreads=omp_get_num_threads(); 
		#pragma omp for schedule(dynamic,chunk)
		for(int i=0;i<swarmsize;i++)
			{
			for(int j=0;j<m_NVAR;j++)
				var(j)=swarm[i][j];									//copy current swarm to temp. vector
			if(opt_print_level>1)									//print particle output if requested
				{
				if(iternumber==0)
					{
					if(nthreads==1) cout << "Initialization, particle " << i << endl;
					else cout << "Initialization, particle " << i << ", processor " << omp_get_thread_num()+1 << " of " << nthreads << endl;
					}
				else
					{
					if(nthreads==1) cout << "Iteration " << iternumber << ", particle " << i << endl;
					else cout << "Iteration " << iternumber << ", particle " << i << ", processor " << omp_get_thread_num()+1 << " of " << nthreads << endl;
					}
				}
			//cout << "Thread " << omp_get_thread_num() << ", PsoBeforeEval, mem:" << (int)(this) << endl;
			objconst = m_model(var); m_numberModelEval++;			//evaluate model
			//cout << "Thread " << omp_get_thread_num() << ", PsoAfterEval, mem:" << (int)(this) << endl;
			penalty=0;
			if(opt_print_level>2)
				{
				for(int j=0;j<m_NOBJ;j++)
					cout << " " << objconst[j];
				cout << "  | ";
				}
			for(int j=0;j<m_NCONS;j++)								//copy constraints from temp. vector
				{
				Cswarm[i][j]=constrViolation(objconst[j+m_NOBJ],j);	//the function returns the constraint violation from the constraint value, according to the defined boundaries
				if(opt_print_level>2) cout << " " << Cswarm[i][j];
				penalty+=Cswarm[i][j];								//adding current constraint violation to penalty
				}
			for(int j=0;j<m_NOBJ;j++)
				Jswarm[i][j]=objconst[j]+penalty*hugeCstrVio();		//adding a huge value if the constraint is violated
			if(opt_print_level>2)
				{
				cout << "  |  " << penalty << "  | ";
				for(int j=0;j<m_NOBJ;j++)
					cout << " " << Jswarm[i][j];
				 cout << endl;
				}
			}
		}	//CLOSE PARALLEL cycle on solutions to be evaluated at the current iteration

	//REVERSE COMMUNICATION FUNCTIONALITY for global optimization algorithms: ask for user action every n iterations:
	bool quitOptimization=false;	//by default, continue optimization (user can modify this below)
	double nSecondsWait=60;			//deafult wait value [s], beyond this the optimization resumes
	double nSecondsRefresh=1;		//refresh rate [s] for optimization resume count-down
	if( m_nIterReverseComm!=0 && iternumber!=0 && iternumber%m_nIterReverseComm==0 )
		{
		cout << "Waiting for user decision, input a single integer number followed by enter (optimization will automatically continue with iteration " << iternumber+1 << " in " << nSecondsWait << " seconds):\n   Input 0 to exit the optimization at the current iteration\n   Input 1 to continue optimization process\n   Input n>1 to continue waiting other n seconds" << endl;
		cout << "Remaining time [s] (" << (int)(nSecondsRefresh*1000) << " ms refresh time): " << endl;
		clock_t endwait;
		endwait = clock () + nSecondsWait * CLOCKS_PER_SEC;
		int useraction; 
		bool continueWait=true;
		while (clock()<endwait && continueWait)
			{
			Sleep((int)(nSecondsRefresh*1000));		//sleep for refresh time seconds
			if( _kbhit() )				//if user has input any key in this wait
				{
				int tmpflag = _cscanf("%d", &useraction); _getch();
				cout << "Requested user action: " << useraction << endl;
				if(useraction<0 || tmpflag!=1)
					cout << "Unrecognized user action, input a single integer number followed by enter: 0 to exit optimization, 1 to continue optimization or any n>1 to continue waiting other n seconds" << endl;
				else if(useraction==0)	//user requests to terminate optimization process
					{
					continueWait=false;
					quitOptimization=true;
					}
				else if(useraction==1)	//user requests to continue optimization process with next iteration(s)
					continueWait=false;
				else					//user requests to continue waiting for useraction seconds
					{
					endwait = clock () + useraction * CLOCKS_PER_SEC;
					continueWait=true;
					cout << "Remaining time [s] (" << (int)(nSecondsRefresh*1000) << " ms refresh time): " << endl;
					}
				}
			else
				{
				double timetogo=((double)(endwait - clock ())) / CLOCKS_PER_SEC;
				cout << " " << (ROUND(timetogo));
				}
			}
		}
	return quitOptimization;
}


// Function used by PSO and DGMOPSO to print to file the optimization variables of the current archive of 1 best 
//  solution for single-objective PSO or of N non dominated solutions for multi-objective DGMOPSO
//INPUTS:
//  * fp: ofstream for the file to write on (MUST HAVE BEEN OPENED AND MUST BE CLOSED OUTSIDE THIS FUNCTION)
//  * iternumber: number of current iteration, to be printed in the first column
//  * size: current size of the archive (1 for single-objective, variable or fixed to maxsize for multi-objective)
//  * archive: matrix of (nparticles x nvars) elements with the optimization variables of the current solutions 
//			contained in the archive (DIMENSIONS AND MEMORY ALLOCATION IS NOT CHECKED HERE, MUST BE TAKEN CARE OF OUTSIDE!)
//OUTPUTS: None, the variables are printed on a row of the file fp

void Optimizer::printArchiveVars(ofstream& fp, const int& iternumber, const int& size, double **archive)
{
	for(int i=0;i<size;i++)
		{
		fp << iternumber << " ";
		for(int j=0;j<m_NVAR;j++)
			fp << archive[i][j] << " ";
		fp << endl;
		}
}

// Function used by PSO and DGMOPSO to print to file the optimization variables of the current archive of 1 best 
//  solution for single-objective PSO or of N non dominated solutions for multi-objective DGMOPSO
//INPUTS:
//  * fp: ofstream for the file to write on (MUST HAVE BEEN OPENED AND MUST BE CLOSED OUTSIDE THIS FUNCTION)
//  * iternumber: number of current iteration, to be printed in the first column
//  * size: current size of the archive (1 for single-objective, variable or fixed to maxsize for multi-objective)
//  * Jarchive: matrix of (nparticles x nobj) elements with the optimization objectives of the current solutions 
//			contained in the archive (DIMENSIONS AND MEMORY ALLOCATION IS NOT CHECKED HERE, MUST BE TAKEN CARE OF OUTSIDE!)
//  * Carchive: matrix of (nparticles x ncstr) elements with the optimization constraints of the current solutions 
//			contained in the archive (DIMENSIONS AND MEMORY ALLOCATION IS NOT CHECKED HERE, MUST BE TAKEN CARE OF OUTSIDE!)
//OUTPUTS: None, the variables are printed on a row of the file fp

void Optimizer::printArchiveObjCstr(ofstream& fp, const int& iternumber, const int& size, double **Jarchive, double **Carchive)
{
	for(int i=0;i<size;i++)
		{
		fp << iternumber << " ";
		for(int j=0;j<m_NOBJ;j++)
			fp << Jarchive[i][j] << " ";
		for(int j=0;j<m_NCONS;j++)
			fp << Carchive[i][j] << " ";
		fp << endl;
		}
}

