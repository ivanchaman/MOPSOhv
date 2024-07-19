
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

#ifndef PSO_H_
#define PSO_H_

#include <iomanip>
#include <stdio.h>
#include <math.h>

#include "optimizer.h"

class PSO: public Optimizer {
public:
	PSO(int NVAR, int NCONS, int MAXITER, int SWARMSIZE, Eigen::VectorXd LOWERBOUND, Eigen::VectorXd UPPERBOUND, Eigen::VectorXd GLOWERBOUND = Eigen::VectorXd::Zero(1), Eigen::VectorXd GUPPERBOUND = Eigen::VectorXd::Zero(1));
	virtual ~PSO();

	void iterationoutput();
	void optimize();

	//setting parameters
	void setSwarmSize(int swarmsize){m_swarmsize = swarmsize;}
	void setInitialInertia(double inertiainit) {m_inertiainit = inertiainit;}
	void setFinalInertia(double inertiafinal) {m_inertiafinal = inertiafinal;}
	void setSelfConfidence(double selfconf){m_selfconf = selfconf;}
	void setSwarmConfidence(double swarmconf){m_swarmconf = swarmconf;}
	void setMutationProbability(double mutprob){m_mutprob = mutprob;}

	//get optimization results
	const Eigen::MatrixXd& getOptimalSolutions(){return m_optimalSolutions;}
	const Eigen::MatrixXd& getOptimalValues(){return m_optimalValues;}

private:
	int m_MAXITER;
	int m_swarmsize;
	double m_inertiainit, m_inertiafinal, m_selfconf, m_swarmconf, m_mutprob;
	Eigen::MatrixXd m_optimalSolutions, m_optimalValues;

};

#endif /* PSO_H_ */
