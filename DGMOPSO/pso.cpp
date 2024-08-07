
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

USING_PART_OF_NAMESPACE_EIGEN
using namespace std;

#ifdef _DEBUG
#ifdef WIN32
int omp_get_thread_num();
int omp_get_num_threads();
#endif
#endif

PSO::PSO(int NVAR, int NCONS, int MAXITER, int SWARMSIZE, Eigen::VectorXd LOWERBOUND, Eigen::VectorXd UPPERBOUND, Eigen::VectorXd GLOWERBOUND, Eigen::VectorXd GUPPERBOUND) : Optimizer(NVAR, NCONS, LOWERBOUND, UPPERBOUND, GLOWERBOUND, GUPPERBOUND)
{
	m_NOBJ = 1;				//single-objective only
	m_MAXITER = MAXITER;	//set number of iterations
	//default values
	m_swarmsize=SWARMSIZE;	//number of particles in the swarm
	m_inertiainit=0.5;		//initial inertia parameter (linear)
	m_inertiafinal=0.5;		//final inertia parameter (linear)
	m_selfconf=1;			//self confidence parameter (constant)
	m_swarmconf=2;			//swarm confidence parameter (constant)
	m_mutprob=0.01;			//mutation probability (constant)
	//double proportialParamScaling=1.0;
	//m_inertiainit*=proportialParamScaling;
	//m_inertiafinal*=proportialParamScaling;
	//m_selfconf*=proportialParamScaling;
	//m_swarmconf*=proportialParamScaling;

	m_optimalSolutions=Eigen::MatrixXd(1,m_NVAR);
	m_optimalValues=Eigen::MatrixXd(1,m_NCONS+m_NOBJ);
}

PSO::~PSO() {}


//============ Function implementing PSO algorithm =============//
/*Implementation of the basic PSO single objective optimizer: constraints, linear inertia weight, constant parameters,
mutation on all variables of the mutated particle, only gbest and pbest (no lbest), no variables discretization, no mailboxes

Reference paper for original single-objective PSO:
- Eberhart J., Kennedy J., Particle Swarm Optimization, 1995
*/

void PSO::iterationoutput(){
		cout << "=========================================Global Optimization========================================="<<endl;
		cout << "================================================PSO=================================================="<<endl;
		cout << "Algorithm: Eberhart J., Kennedy J., Particle Swarm Optimization, 1995"<<endl;
		cout << endl;
		cout << "-----------------------------------Problem definition------------------------------------"<< endl;
		cout << "Dimension: "<<setw(22)<<m_NVAR<<endl;
		cout << "Number of constraints: " <<setw(10)<<m_NCONS<<endl;
		cout << setprecision(5);
		cout << "---------------------------------------Parameters----------------------------------------"<<endl;
		cout << "Maximum number of iterations: " <<setw(19)<<m_MAXITER<<endl;
		cout << "Number of particle in the swarm: " <<setw(16)<<m_swarmsize<<endl;
		cout << "Initial inertia parameter (linear): " <<setw(13)<<m_inertiainit<<endl;
		cout << "Final inertia parameter (linear): " <<setw(15)<<m_inertiafinal<<endl;
		cout << "Self confidence parameter (constant): " <<setw(11)<<m_selfconf<<endl;
		cout << "Swarm confidence parameter (constant): " <<setw(10)<<m_swarmconf<<endl;
		cout << "Mutation probability (constant): "<<setw(16)<<m_mutprob<<endl;
		cout << setprecision(10);
		cout << endl;
}

void PSO::optimize(){

	//Multi-threaded (parallel computing) options definition:
	bool useOpenMP=false; bool useMPI=false; bool useOpenCL=false;
	if(m_multiThread && m_multiThreadType==1)
		{ useOpenMP=true; if(opt_print_level>0) cout << "Parallel PSO execution with OpenMP selected (multi-processor, shared memory)" << endl; }
	if(opt_print_level>0 && useOpenMP)
		{
		#pragma omp parallel
			{ if(!omp_get_thread_num()) cout << "Total number of threads: " << omp_get_num_threads() << endl; }
		}
	if(m_multiThread && m_multiThreadType==2)
		{ useMPI   =true; if(opt_print_level>0) cout << "Parallel PSO execution with MPI selected (multi-processor, distributed memory)" << endl; }
	if(m_multiThread && m_multiThreadType==3)
		{ useOpenCL=true; if(opt_print_level>0) cout << "Parallel PSO execution with OpenCL selected (GPU computing, NVIDIA or ATI graphics card required)" << endl; }
	if(opt_print_level>0 && !useOpenMP && !useMPI && !useOpenCL) cout << "Sequential PSO execution selected (single-processor)" << endl;

	//Print PSO initialization message:
	if(opt_print_level>0){
		iterationoutput();
	}
	
	//Open optimization output files:
	string file="PSO_optVAR_"+m_name+".dat"; ofstream fpSolutions(file.data());
	file="PSO_optOBJandCONS_"+m_name+".dat"; ofstream fpOptValues(file.data());
	if(!fpSolutions || !fpOptValues) { cout << "Cannot open optimization output files.\n"; return; }
	fpSolutions.setf(ios::scientific); fpSolutions.precision(16);
	fpOptValues.setf(ios::scientific); fpOptValues.precision(16);
	//fpSolutions<<"This file contains the optimal solution found by the single-objective PSO algorithm"<<endl;
	//fpOptValues<<"This file contains the optimal value and the constraints violations (<0 violated, 0 satisfied) found by the single-objective PSO algorithm"<<endl;
	
	//Local variables:
    double **vel,**swarm,**pbest,**gbest,**Jswarm,*Jpbest,Jgbest,**Cswarm,**Cpbest,**Cgbest,
                            tempselfconf,tempswarmconf,tempinertia,temprand,penalty,**tmpJgbest;

    //memory allocations:
    swarm = dmatrix_allocation(m_swarmsize,m_NVAR);	//particles positions (position=optimization variables)
    pbest = dmatrix_allocation(m_swarmsize,m_NVAR);	//personal best positions
    gbest =  dmatrix_allocation(1,m_NVAR);			//global best position
    Jswarm = dmatrix_allocation(m_swarmsize,m_NOBJ);//objective value for the swarm
    Jpbest = new double [m_swarmsize];				//personal best function value
    Cswarm = dmatrix_allocation(m_swarmsize,m_NCONS);//constraint values for the swarm
    Cpbest = dmatrix_allocation(m_swarmsize,m_NCONS);//constraint values for the pbest
    Cgbest = dmatrix_allocation(1,m_NCONS);			//gbest constraints values
    vel = dmatrix_allocation(m_swarmsize,m_NVAR);	//particles velocities
    tmpJgbest = dmatrix_allocation(1,m_NOBJ);		//gbest objective value
	//NOTE: Jswarm, gbest and Cgbest could be vectors instead of matrixes, and tmpJgbest could be a scalar, but they are all 
	//		defined AS MATRIXES to be consistent with the definitions in DGMOPSO --> to use the same common functions 
	//		Optimizer::evaluateParticles(), Optimizer::printArchiveVars() and Optimizer::printArchiveObjCstr()


    /* ---------------------------------------------------------------- MAIN ------------------------------------------------------------------------*/

    /* ----------------------------------------------------------------------------------------------------------------------------------------------*/
    /* --------------------------------------------------------- 1. INITIALIZATION -------------------------------------------------------------------*/
    /* ----------------------------------------------------------------------------------------------------------------------------------------------*/
	if(opt_print_level>0)
		cout << "-Random initialization of the swarm------------------------------------------------------"<< endl;
    Jgbest=hugeCstrVio()*hugeCstrVio();	//initial high number for Jgbest
	Eigen::MatrixXd initMatrix;
	bool flagInit;
    // 1a. Initialize particles either randomly (case 0), from stored previous solutions (case 1) or from file (case 2):
	switch(m_initialization){
	case 0:
		for(int i=0;i<m_swarmsize;i++)
		{
			for(int j=0;j<m_NVAR;j++)
			{
				temprand=((double)rand()/RAND_MAX);			//define random number in [0;1]
				swarm[i][j]=m_LOWERBOUND[j]+(m_UPPERBOUND[j]-m_LOWERBOUND[j])*temprand;	//positions random initialization
				//swarm[i][j]=0.5;	//positions random initialization
				pbest[i][j]=swarm[i][j];					//initial pbest=initial swarm
				vel[i][j]=0;								//null initial velocity
			}
		}
		break;
	case 1:
		for(int i=0;i<m_swarmsize;i++)
		{
			for(int j=0;j<m_NVAR;j++)
			{
				if(i<m_optimalSolutions.rows())
					swarm[i][j]=m_optimalSolutions(i,j);
				else{
					temprand=((double)rand()/RAND_MAX);			//define random number in [0;1]
					swarm[i][j]=m_LOWERBOUND[j]+(m_UPPERBOUND[j]-m_LOWERBOUND[j])*temprand;	//positions random initialization
				}
				pbest[i][j]=swarm[i][j];					//initial pbest=initial swarm
				vel[i][j]=0;								//null initial velocity
			}
		}
		break;
	case 2:
		flagInit = loadFirstGuessSolutionFromFile(m_NVAR,m_filenameInit,initMatrix);
		for(int i=0;i<m_swarmsize;i++)
		{
			for(int j=0;j<m_NVAR;j++)
			{
				if(i<initMatrix.rows())
					swarm[i][j]=initMatrix(i,j);	//positions random initialization
				else{ 
					temprand=((double)rand()/RAND_MAX);	
					swarm[i][j]=m_LOWERBOUND[j]+(m_UPPERBOUND[j]-m_LOWERBOUND[j])*temprand;
				}
				pbest[i][j]=swarm[i][j];					//initial pbest=initial swarm
				vel[i][j]=0;								//null initial velocity
			}
		}
		break;
	default: cout<<"Warning: Wrong initialization parameters for PSO optimization process"<<endl; break; 
	}

    // 1b. Perform first model evaluation for all particles, with parallel paradigms if required:
	evaluateParticles(0,useOpenMP,useMPI,useOpenCL, m_swarmsize, swarm,Jswarm,Cswarm);	// swarm --> m_model --> Jswarm,Cswarm

    // 1c. Update (sequentially) Cpbest,Jpbest,Cgbest,Jgbest:
	for(int i=0;i<m_swarmsize;i++)
		{
		//pbest initialized with initialization solution:
		for(int j=0;j<m_NCONS;j++)
			Cpbest[i][j]=Cswarm[i][j];
		Jpbest[i]=Jswarm[i][0];
		if(Jswarm[i][0]<Jgbest)									//Updating gbest if the function value is lower than previous gbest
			{
			Jgbest=Jswarm[i][0];
			for(int j=0;j<m_NCONS;j++)
				Cgbest[0][j]=Cswarm[i][j];							//initial Cgbest values
			for(int j=0;j<m_NVAR;j++)
				gbest[0][j]=swarm[i][j];							//initial gbest values
			}
		}
		
	// 1d. Write initial objectives and constraints values to file:
    penalty=0;
    for(int j=0;j<m_NCONS;j++)
	    penalty+=Cgbest[0][j];
    tmpJgbest[0][0]=Jgbest-penalty*hugeCstrVio();
	printArchiveVars(fpSolutions,0,1,gbest);				//print to file opt. vars of gbest at initialization
	printArchiveObjCstr(fpOptValues,0,1,tmpJgbest,Cgbest);	//print to file opt. objs and cstrs of gbest at initialization

	//Print to screen:
	if(opt_print_level>0) 
		{
		cout << "-Swarm initialization completed ---------------------------------------------------------"<< endl;
		cout <<setw(5)<< "Iter"<<setw(15)<<"Touched LB "<<setw(12)<<"Touched UB "<<setw(18)<<"Jbest"<< setw(34)<<"Constraints violation"<<endl; 
		cout << "[init]"<<setw(11)<<0<<" | "<<setw(9) <<0<< " | "<<setw(17)<<tmpJgbest[0][0]<< " | "<<setw(16);
		if(m_NCONS>0)
			{
			for(int j=0;j<m_NCONS;j++)
				cout<<Cgbest[0][j]<<" ";
			cout<<endl;
			}
		else cout<<"-----"<<endl;
		}

	/* ----------------------------------------------------------------------------------------------------------------------------------------------*/
    /* ---------------------------------------------------------- 2. MAIN CYCLE ---------------------------------------------------------------------*/
    /* ----------------------------------------------------------------------------------------------------------------------------------------------*/
    //MAIN CYCLE ON k: maxiter iterations of the PSO algorithm (no other stop criterion)
	long int countUB=0,countLB=0;
	double bestObj = hugeCstrVio();
    // 2a: Update particles position and velocity, including mutation:
    bool quitOptimizationUponUserRequest=false;
	for(int k=1;k<m_MAXITER+1;k++)
        {
        tempinertia=(m_inertiafinal-m_inertiainit)/m_MAXITER*k+m_inertiainit;
        for(int i=0;i<m_swarmsize;i++)				//for each particle i at iteration k
            {
			// Update velocity and swarm according to PSO equation:
            temprand=((double)rand()/RAND_MAX);
            tempselfconf=temprand*m_selfconf;		//randomizing the value of self confidence for the particle i at iteration k
            temprand=((double)rand()/RAND_MAX);
            tempswarmconf=temprand*m_swarmconf;		//randomizing the value of swarm confidence for the particle i at iteration k
            for(int j=0;j<m_NVAR;j++)	//for each decision variable j of particle i at iteration k
                {
				//if(j==0) cout<< " Before: " << vel[i][j] << " "  << swarm[i][j];
                vel[i][j]=tempinertia*vel[i][j]+tempselfconf*(pbest[i][j]-swarm[i][j])+tempswarmconf*(gbest[0][j]-swarm[i][j]);	//velocity update
                swarm[i][j]=swarm[i][j]+vel[i][j];			//position update
				//if(j==0) cout<< " --> After: " << vel[i][j] << " " << swarm[i][j] << " "  << pbest[i][j] << " "  << gbest[j] << endl;
                //Restraint the position within the boundaries:
                if (swarm[i][j]<m_LOWERBOUND[j])
					{
					if (j==0) countLB++;
					swarm[i][j]=m_LOWERBOUND[j];		//set variable to lower limit
					vel[i][j]=-vel[i][j];	//invert velocity that leads out of the boudaries
                    }
				else if (swarm[i][j]>m_UPPERBOUND[j])
					{
					if (j==0) countUB++;
					swarm[i][j]=m_UPPERBOUND[j];	//set variable to upper limit
					vel[i][j]=-vel[i][j];			//invert velocity that leads out of the boudaries
					}
                }

			// Mutate particles position:
            temprand=((double)rand()/RAND_MAX);
            if (temprand<m_mutprob)					//decide if mutate or not the particle
				{
				for(int j=0;j<m_NVAR;j++)			//for each decision variable j of particle i at iteration k
                    {
                    temprand=((double)rand()/RAND_MAX);
                    swarm[i][j]=m_LOWERBOUND[j]+(m_UPPERBOUND[j]-m_LOWERBOUND[j])*temprand;	//position random reset is used as mutation
                    }
				}
			}
		
		// 2b. Perform model evaluation for all particles, with parallel paradigms if required:
		quitOptimizationUponUserRequest=evaluateParticles(k,useOpenMP,useMPI,useOpenCL, m_swarmsize, swarm,Jswarm,Cswarm);	// swarm --> m_model --> Jswarm,Cswarm

	    // 2c. Check feasible only option and update Cpbest,Jpbest,Cgbest,Jgbest:
        for(int i=0;i<m_swarmsize;i++)				//for each particle i at iteration k
			{
			//Check feasible only option:
			if(m_feasibleOnly)
				{
				int counter = 0;
				for(int j=0;j<m_NCONS;j++)
					if(Cswarm[i][j]==0) counter++; 
				if(counter==m_NCONS && bestObj>Jswarm[i][0])//if solution i is feasible and it's the best feasible found so far
					{
					bestObj = Jswarm[i][0];					//update bestObj
		            for(int j=0;j<m_NVAR;j++)
						m_bestFeasible(j)=swarm[i][j];		//update m_bestFeasible
					}
				}
			//Personal best update:
            if (Jswarm[i][0]<Jpbest[i])						//if current J value is lower than previous Jpbest
                {
                Jpbest[i]=Jswarm[i][0];						//Jpbest value update
                for(int j=0;j<m_NVAR;j++)
					pbest[i][j]=swarm[i][j];				//pbest update
                for(int j=0;j<m_NCONS;j++)
					Cpbest[i][j]=Cswarm[i][j];				//Cpbest update
                }
            //Global best update:
            if (Jswarm[i][0]<Jgbest)							//if current J value is lower than previous Jgbest
                {
                Jgbest=Jswarm[i][0];							//Jgbest update
                for(int j=0;j<m_NVAR;j++)
					gbest[0][j]=swarm[i][j];					//gbest update
                for(int j=0;j<m_NCONS;j++)
					Cgbest[0][j]=Cswarm[i][j];					//Cgbest update
                }
			}	//close for cycle on the particles

        // 2d. Write the objectives, constraints and variables values of current gbest to output files, print to screen:
        penalty=0;
        for(int j=0;j<m_NCONS;j++)
			penalty+=Cgbest[0][j];
        tmpJgbest[0][0]=Jgbest-penalty*hugeCstrVio();
		printArchiveVars(fpSolutions,k,1,gbest);				//print to file opt. vars of gbest at k-th iteration
		printArchiveObjCstr(fpOptValues,k,1,tmpJgbest,Cgbest);	//print to file opt. obj and cstrs of gbest at k-th iteration

		if(opt_print_level>0)
			{
			cout << endl << "["<<setw(4)<<k<<"]"<<setw(11)<<countLB<<" | "<<setw(9) <<countUB << " | "<<setw(17)<<tmpJgbest[0][0]<< " | "<<setw(16);
			if(m_NCONS>0)
				{
				for(int j=0;j<m_NCONS;j++)
					cout<<Cgbest[0][j]<<" ";
				cout<<endl;
				}
			else cout<<"-----"<<endl;
		}

        // 2e. Break iterations cycle if feasible only option is active and a feasible solution has been found:
		if(m_feasibleOnly && bestObj < hugeCstrVio())
			{ 
			if(opt_print_level>0)
				{
				cout << endl;
				cout << "-Feasibilty reached ---------------------------------------------------------------------"<< endl;
				cout << "-Number of Model evaluation: "<<m_numberModelEval<< endl;
				cout << endl;
				}	
			break;
			}
		// or if the user has requested to stop the optimization process
		if(quitOptimizationUponUserRequest)
			{
			if(opt_print_level>0)
				{
				cout << endl;
				cout << "-User requested termination -------------------------------------------------------------"<< endl;
				cout << "-Number of Model evaluation: "<<m_numberModelEval<< endl;
				cout << endl;
				}
			break;
			}
        }	//close for cycle on the iterations

	/* ----------------------------------------------------------------------------------------------------------------------------------------------*/
    /* --------------------------------------------------------- 3. POSTPROCESS ---------------------------------------------------------------------*/
    /* ----------------------------------------------------------------------------------------------------------------------------------------------*/
	if(!m_feasibleOnly) // operations done just when  m_feasibleOnly is off
		{
	    //3a. COPY FINAL RESULTS FOR THE NEXT OPTIMIZATION RUN OR FOR THE FINALIZATION:
		for(int j=0;j<m_NVAR;j++)
			m_optimalSolutions(0,j)=gbest[0][j];		//copy gbest variables
		m_optimalValues(0,0)=Jgbest;				//copy gbest objective
		for(int j=0;j<m_NCONS;j++)
            m_optimalValues(0,j+1)=Cgbest[0][j];		//copy gbest constraints
	
		//3b. WRITE OUTPUT ON STD OUTPUT:
		if(opt_print_level>0)
			{
			cout << endl;
			cout << "-Max number of iteration reached --------------------------------------------------------"<< endl;
			cout << "-Number of Model evaluation: " << m_numberModelEval << endl;
			cout << endl;
			}
	}//close the  !m_feasibleOnly if statement
	m_numberModelEval = 0;
	if(opt_print_level>0)
		{
		cout << "SOLUTION FOUND: " << endl;
		for(int j=0;j<m_NVAR;j++)
			cout<<gbest[0][j]<<" ";
		cout << endl;
		cout << "Objective value: " << Jgbest << endl;
		if(m_NCONS>0)
			{
			cout<<"Constraints violation: ";
			for(int j=0;j<m_NCONS;j++)
				cout << Cgbest[0][j] << " ";
			}
		cout << endl;
		}

	if(opt_print_level>0)
		{
		cout << "Writing to files: PSO_solutions_"+m_name+".dat"<< endl;
		cout << "                  PSO_paretofront_"+m_name+".dat"<<endl;
		cout << "Number of Model evaluation: "<<m_numberModelEval<< endl;
		}

	//3c. CLOSE OUTPUT FILES and FREE MEMORIES:
	fpSolutions.close();
	fpOptValues.close();

    dmatrix_free(vel,   m_swarmsize,m_NVAR);
    dmatrix_free(swarm, m_swarmsize,m_NVAR);
    dmatrix_free(Cswarm,m_swarmsize,m_NCONS);
    dmatrix_free(Jswarm,m_swarmsize,m_NOBJ);
    dmatrix_free(pbest, m_swarmsize,m_NVAR);
    dmatrix_free(Cpbest,m_swarmsize,m_NCONS);
    delete[] Jpbest;
    dmatrix_free(gbest,    1,m_NVAR);
    dmatrix_free(Cgbest,   1,m_NCONS);
    dmatrix_free(tmpJgbest,1,m_NOBJ);

}
