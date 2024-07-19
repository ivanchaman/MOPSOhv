
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

#ifndef DGMOPSO_H_
#define DGMOPSO_H_

#include <search.h>
#include <stdlib.h>
#include <iomanip>

#include "optimizer.h"


class DGMOPSO: public Optimizer {
public:
	DGMOPSO(int NVAR, int NOBJ, int NCONS, int MAXITER, int SWARMSIZE, int NsolPareto, Eigen::VectorXd LOWERBOUND, Eigen::VectorXd UPPERBOUND, Eigen::VectorXd GLOWERBOUND = Eigen::VectorXd::Zero(1), Eigen::VectorXd GUPPERBOUND = Eigen::VectorXd::Zero(1));
	virtual ~DGMOPSO();

	void iterationoutput();
	void optimize();

	//setting parameters
	void setArchiveSize(int archivesize){m_archivesize = archivesize; 
	m_swarmsize = archivesize;
	m_optimalSolutions.resize(archivesize, m_NVAR);
	m_paretoFront.resize(archivesize, m_NOBJ+m_NCONS);
	m_oldSwarm.resize(m_swarmsize, m_NVAR);
	m_oldPbest.resize(m_swarmsize,m_NVAR);
	m_oldVel.resize(m_swarmsize,m_NVAR); 
	m_oldJpbest.resize(m_swarmsize,m_NOBJ);
	m_oldJswarm.resize(m_swarmsize,m_NOBJ);}
	void setGridBisectionOut(int gridbisect_out){m_gridbisect_out = gridbisect_out;}
	void setGridBisectionIn(int gridbisect_in){m_gridbisect_in = gridbisect_in;}
	void setEps(double eps){m_eps = eps;}
	void setInitialInertial(double inertiainit){m_inertiainit = inertiainit;}
	void setFinalInertial(double inertiafinal){m_inertiafinal = inertiafinal;}
	void setInitialSelfConfidence(double selfconfinit){m_selfconfinit = selfconfinit;}
	void setFinalSelfConfidence(double selfconffinal){m_selfconffinal = selfconffinal;}
	void setInitialSwarmConfidence(double swarmconfinit){m_swarmconfinit = swarmconfinit;}
	void setFinalSwarmConfidence(double swarmconffinal){m_swarmconffinal = swarmconffinal;}
	void setInitialMutationProbability(double mutprobinit){m_mutprobinit = mutprobinit;}
	void setFinalMutationProbability(double mutprobfinal){m_mutprobfinal = mutprobfinal;}
	void setInitialMutationDistribution(double mutdistrinit){m_mutdistrinit = mutdistrinit;}
	void setFinalMutationDistribution(double mutdistrfinal){m_mutdistrfinal = mutdistrfinal;}

	const Eigen::MatrixXd& getOptimalSolutions(){return m_optimalSolutions;}
	const Eigen::MatrixXd& getParetoFront(){return m_paretoFront;}

	void setOptimalSolutions(Eigen::MatrixXd& sol){m_optimalSolutions = sol;}

private:
	int* integer2binary(int n, int nbits);
	int binary2integer(int *binary, int nbits);
	void check_boundaries_real(double **position, double **velocity, int i, int nvar, double *LB, double *UB);
	int checkgridbound(double **Jarchive,int *erase,int counterase,double **gridbound,double *newsol,int pos,int nobj,int *minB,int *maxB);
	void define_fitness_in(int Nocc_out, int *gridocc_out, int *Nocc_in, int **gridocc_in, int **gridNO_in, double **fitness_in);
	void define_fitness_out(int Nocc_out, int *gridNO_out, int *gridocc_out, double *fitness_out);
	int roulette_wheel_selection(double *fitness, int size);
	int dgchoose_leader(int solgrid_out,double *fitness_out,double **fitness_in,int *archivegrid_out,int *archivegrid_in,int **gridID_out,int ***gridID_in,int *gridNO_out,int *gridocc_out,int *gridoccref_out,int **gridocc_in,int Nocc_out,int *Nocc_in);
	void erase_sol(double **archive,double **Jarchive,double **Carchive,int j_erase,int l_erase,int sol_erase,int size,int *Nocc_out,int *Nocc_in,int *archivegrid_out,int *archivegrid_in,int **gridID_out,int ***gridID_in,int *gridNO_out,int **gridNO_in,int *gridocc_out,int **gridocc_in,int *gridoccref_out,int **gridoccref_in,int **gridIDref_out,int ***gridIDref_in,double *fitness_out,double **fitness_in,int nvar,int nobj,int ncons,int *miB, int *maB);
	int findgrid(double *sol, double **gridbound, int gridbisect, int nobj);
	void findgridbound(double **J, double **gridbound, int size, int nobj, int *minB, int *maxB);
	void findinnerbound(double **outer,double ***inner, int bisect, int nobj);
	int findnondom(double **archive, double **swarm, double **Jarchive, double **Jswarm, double **Carchive, double **Cswarm, int swarmsize, int nvars, int nobj, int ncons);
	void init_real(double **swarm,double **pbest,double **vel, double *LB, double *UB, int swarmsize, int nvar);
	void init_fromPreviousIteration(double **swarm,double **pbest,double **vel, double *LB, double *UB, int swarmsize, int nvar);
	void init_fromFile(double **swarm,double **pbest,double **vel, double *LB, double *UB, int swarmsize, int nvar);
	int isequal(double* fun1, double* fun2, int nobj);
	void move_particle_real(double **position, double **velocity, double **archive, double **pbest, int i, int k, int nvar, int archivesize, double swarmconf, double selfconf, double inertiainit, double inertiafinal,  int maxiter, int leader);
	void mutate_nsga2(double **position, int i, int nvar, double *LB, double *UB, double mutdistr);
	int pareto2tournament(double* fun1, double* fun2, int nobj);
	void reset_grid(int *archivegrid_out,int *archivegrid_in, int ngrids_out, int ngrids_in, int size, int nobj, int **gridID_out, int ***gridID_in, int *gridNO_out, int **gridNO_in, int *gridocc_out, int **gridocc_in, int *Nocc_out, int *Nocc_in, int *gridoccref_out, int **gridoccref_in, int **gridIDref_out, int ***gridIDref_in);
	int update_gbestandgrid(double **archive,double **Jarchive,double **Carchive,double *sol,double *Jsol,double *Csol,int size,int *Nocc_out,int *Nocc_in,int *agr_out,int *agr_in,int **gridID_out,int ***gridID_in,int *gridNO_out,int **gridNO_in,int *gridocc_out,int **gridocc_in,int *gridoccref_out,int **gridoccref_in,int **gridIDref_out,int ***gridIDref_in,int solgrid_out,int solgrid_in,int nvar,int nobj,int ncons,double *fitness_out,double **fitness_in,double **ngrb,int *miB,int *maB);
	int update_pbest_real_fitness(double **Jswarm,double **Jpbest,double **swarm,double **pbest,int i,int nvar,int nobj,double **gridbound_out,double ***gridbound_in,int gridbisect_out,int gridbisect_in,int *gridNO_out,int **gridNO_in,int *solgrid_out,int *solgrid_in,double *pbest_fitness_out,double *pbest_fitness_in);

	int m_NsolPareto, m_MAXITER;
	int m_swarmsize;
	int m_archivesize;			//non-dominated archive size
	int m_gridbisect_out;		//number of successive bisections for each dimension to create the outer grid
	int m_gridbisect_in;		//number of further successive bisections for each dimension to create the inner grid
	double m_eps;				//percentage error on the external boundaries over which the grid is resetted
	double m_inertiainit;		//initial inertia parameter (linear)
	double m_inertiafinal;		//final inertia parameter (linear)
	double m_selfconfinit;		//initial self confidence parameter (linear)
	double m_selfconffinal;		//final self confidence parameter (linear)
	double m_swarmconfinit;		//initial swarmconfidence parameter (linear)
	double m_swarmconffinal;	//final swarmconfidence parameter (linear)
	double m_mutprobinit;		//initial mutation probability (linear)
	double m_mutprobfinal;		//final mutation probability (linear)
	double m_mutdistrinit;		//initial NSGA-II mutation distribution (linear)
	double m_mutdistrfinal;		//final NSGA-II mutation distribution (linear)

	Eigen::MatrixXd m_optimalSolutions, m_paretoFront;
	Eigen::MatrixXd m_oldSwarm, m_oldPbest, m_oldVel, m_oldJpbest, m_oldJswarm;
};

#endif /* DGMOPSO_H_ */
