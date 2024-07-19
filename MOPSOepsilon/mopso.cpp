/*
 MOPSO con dominancia epsilon y clusters
*/
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <time.h>
#include "efile.h"
using namespace std;

MOPSO::~MOPSO(){
}

MOPSO::MOPSO(int _dims, int _objs,int _parts,int _nclusters,int _gmax,double *_lb,double *_ub ){
    srand(time(0));
    ndimensions=_dims;
    nobjectives=_objs;
    nparticles=_parts;
    nclusters=_nclusters;
    gmax=_gmax;
    gen=0;
    lb.resize(ndimensions);
    ub.resize(ndimensions);
    particles.resize(nparticles);

    for(int _j(0);_j<ndimensions;_j++){
        lb[_j]=_lb[_j];
        ub[_j]=_ub[_j];
    }
    archive.init(ndimensions, nobjectives, nclusters,100);
}

void MOPSO::initialize(){
    C1=2.0; C2=2.0;  W=0.4;
    for(int _i(0);_i<nparticles;_i++){
        for(int _j(0);_j<ndimensions;_j++){
            particles[_i].vel[_j] = 0;
            particles[_i].x[_j]=rnd(lb[_j],ub[_j]);
            particles[_i].xpbest[_j]=particles[_i].x[_j];
        }
        function(_i);
        for(int _j(0);_j<nobjectives;_j++){
            particles[_i].fxpbest[_j]=particles[_i].fx[_j];
        }
        archive.add(particles[_i],-1);
    }
}

void MOPSO::execute(){
    initialize();
    for(gen=0;gen < gmax;gen++){
        //printf("%d",gen);
        archive.updatematrix=true;
        archive.hierarchicalClustering();
        flight();
        if(gen%10==0)
            cout<<gen<<""<<endl;
    }
    archive.output();
}

void MOPSO::flight(){
    for(int _i(0);_i<nparticles;_i++){
        int _whichcluster=(int)_i/(nparticles/nclusters);
        int _gbestselected;
        _gbestselected=archive.selectClusteredRandomSolution(_whichcluster);
        for(int _k(0);_k<5;_k++){
            for(int _j(0);_j<ndimensions;_j++){
                particles[_i].vel[_j]=W*particles[_i].vel[_j]+C1*rnd(0,1)*(archive.solutions[_gbestselected].x[_j]-particles[_i].x[_j])+C2*rnd(0,1)*(particles[_i].xpbest[_j]-particles[_i].x[_j]);
                particles[_i].x[_j]+=particles[_i].vel[_j];
                if(particles[_i].x[_j]<lb[_j])
                    particles[_i].x[_j]=lb[_j];
                if(particles[_i].x[_j]>ub[_j])
                    particles[_i].x[_j]=ub[_j];
            }
            function(_i);
            int _tmp=archive.domine1(particles[_i].fx,particles[_i].fxpbest);
            if(_tmp==11||_tmp==1){
                copy(particles[_i].fxpbest,particles[_i].fx);
                copy(particles[_i].xpbest,particles[_i].x);
                archive.add(particles[_i],(int)_i/(nparticles/nclusters));
            }
        }
    }
}

void MOPSO::perturbation(int _whichparticle){
    int _dimension=0;
    double _lb,_ub,_rango;
    double _gt=(double)gen/gmax;
    double  _pM=pow(_gt,1.7)-2.0*_gt+1.0;
    int _flag(0);
    if((rnd(0.0,1.0)>_pM)&&_flag<=ndimensions){
        _dimension= (int) rnd(0,ndimensions);
        _rango=(ub[_dimension]-lb[_dimension])*_pM/2.0;//totGen
        if(particles[_whichparticle].x[_dimension]-_rango<lb[_dimension])
            _lb=lb[_dimension];
        else
            _lb=particles[_whichparticle].x[_dimension]-_rango;
        if(particles[_whichparticle].x[_dimension]+_rango>ub[_dimension])
            _ub=ub[_dimension];
        else
            _ub=particles[_whichparticle].x[_dimension]+_rango;
        particles[_whichparticle].x[_dimension]=rnd(_lb,_ub);
        _flag++;
    }
}

void MOPSO::function(int _w){

  /*
  double sum1=0.0,sum2=0.0;

  for(int i=0;i<ndimensions;i++){
    if(i<ndimensions-1)
      sum1+=-10.0*exp((-0.2)*sqrt( pow(particles[_w].x[i],2.0)+pow(particles[_w].x[i+1],2.0)));
    sum2+=pow(fabs(particles[_w].x[i]),0.8)+5*sin(pow(particles[_w].x[i],3.0));
  }
  particles[_w].fx[0]=sum1;
  particles[_w].fx[1]=sum2;
  //*/
  //     }
  /*
    double f;
    double g;
    double h;
    double x1=particles[_w].x[0];
    double x2=particles[_w].x[1];

    f=x1;
    g=11+pow(x2,2)-10*cos(2*3.1415926*x2);
    if(f<=g)
    h=1-pow((f/g),(double)(1.0/2.0));
    else
    h=0;
     particles[_w].fx[0]=f;
    particles[_w].fx[1]=g*h;

    //    */

  //*zdt1

    double g(0),h(0);
    particles[_w].fx[0]=particles[_w].x[0];
    for(int _i(1);_i<ndimensions;_i++)
        g+=(double)particles[_w].x[_i]/((double)ndimensions-1.0);
    g*=9.0;
    g+=1.0;
    h=1.0-sqrt((double)particles[_w].fx[0]/g);
    particles[_w].fx[1]=g*h;
  //*/
}

double MOPSO::rnd(double _min,double _max){
    return((double)(_min + ((double)(_max-_min)*rand()/(double)(RAND_MAX+_min))));
}

void MOPSO::copyx(vector<double> &_a,vector<double> &_b){
    for(int _i(0);_i<ndimensions;_i++){
        _a[_i]=_b[_i];
    }
}

void MOPSO::copyfx(vector<double> &_a,vector<double> &_b){
    for(int _i(0);_i<nobjectives;_i++){
        _a[_i]=_b[_i];
    }
}

void MOPSO::copy(vector<double> &_a,vector<double> &_b){
    for(int _i(0);_i<_a.size();_i++){
        _a[_i]=_b[_i];
    }
}

int main()
{
  /*

  int _ndims=3;
  int _nobjs=2;
  int _parts=40;
  int _nclusters=4;
  int _gmax=40;
  double *_lbound=new(double[_ndims]);//=0.1;
  double *_ubound=new(double[_ndims]);//=1;
  _lbound[0]=-5.0;
  _lbound[1]=-5.0;
  _lbound[2]=-5.0;
  _ubound[0]=5.0;
  _ubound[1]=5.0;
  _ubound[2]=5.0;
  //*/

  /*//dtlz6
  int _ndims=22;
  int _nobjs=3;
  int _parts=40;
  int _nclusters=4;
  int _gmax=100;
  double *_lbound=new(double[_ndims]);//=0.1;
  double *_ubound=new(double[_ndims]);//=1;
  for(int _i(0);_i<_ndims;_i++){
    _lbound[_i]=0.0;
    _ubound[_i]=1.0;
  }
    //*/
  //*//zdt1
    int _ndims=30;
    int _nobjs=2;
    int _parts=40;
    int _nclusters=4;
    int _gmax=100;
    double *_lbound=new(double[_ndims]);//=0.1;
    double *_ubound=new(double[_ndims]);//=1;
    for(int _i(0);_i<_ndims;_i++){
        _lbound[_i]=0.0;
        _ubound[_i]=1.0;
    }
    //*/
    MOPSO  pso(_ndims,_nobjs,_parts,_nclusters , _gmax,_lbound,_ubound);
    pso.execute();
    delete [] _lbound;
    delete [] _ubound;
}
