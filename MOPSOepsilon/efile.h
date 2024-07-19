#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <math.h>
using namespace std;
#define DIMENSIONS 30
#define OBJS 2
#define NPARTICLES 40
#define MAXGBESTS 200

class EPSs{
    public:
        vector<double> eps;
};

class Particle{
    public:
        vector<double> vel;
        vector<double> x;
        vector<double> nx;
        vector<double> fx;
        int subswarm;
        vector<double> xpbest;
        vector<double> fxpbest;
        vector<EPSs> epss;
  //Particle::
        Particle(){
            vel.resize(DIMENSIONS);
            fx.resize(OBJS);
            fxpbest.resize(OBJS);
            x.resize(DIMENSIONS);
            nx.resize(DIMENSIONS);
            xpbest.resize(DIMENSIONS);
        }
  //Particle::
        Particle(int _dimensions, int _objectives){
            vel.resize(_dimensions);
            fx.resize(_objectives);
            fxpbest.resize(_objectives);
            x.resize(_dimensions);
            nx.resize(_dimensions);
            xpbest.resize(_dimensions);
        }
};

class EFILE{
    vector<double> tlb,tub;
    int ndimensions;
    int nobjectives;
    int nsolutions;
    int nclusters;
    int allnsolutions;
    int maxsolutions;
    bool hasentered;
    bool firsttime;
    vector<double> xtlb,xtub;
    vector<int> nsolutionsbyswarm;
    bool hasclusteredfirsttime;
    int pmsize;
    double proximitymatrix[MAXGBESTS][MAXGBESTS];
    int proximitymatrixindex[800][800];
    public:
        bool updatematrix;
        vector<double> EPS;
        vector<Particle> solutions;
        vector<Particle> allsolutions;
        void init(int _ndimensions, int _nobjectives, int _nclusters,int _maxsolutions);
        EFILE();
        int nClusteredSolutions(int _whichcluster);
        int selectClusteredSolution(int _whichcluster,int _whichsolution);
  //  void add(vector<Particle> &_sol);
  //void add(Particle _s);
        void output();
        bool add(Particle &_s,int _whichcluster=-1);
        bool addAll(Particle &_s);
        void update();
        bool reAdd();
        void finalSolutions();
        void deleteExcedent();
        int falseReAdd();
        void printSolutions();
        int  domine(vector<double> &_a,vector<double> &_b);
        int  domine1(vector<double> &_a,vector<double> &_b);
        int  domine2(vector<double> &_a,vector<double> &_b);
        double  euclideanDistance(vector<double> &_a,vector<double> &_b);
        void normalizeDistances1();
        void normalizeDistances2();
        int deleteFromCluster(int _whichsolution);
        void add2Cluster(int _whichcluster);
        int substractFromTo(int _from, int _to, int _quantity);
        int updateClusters(vector<int> &_whichparticles);
        void normalizeDistances();
        void hierarchicalClustering();
        void createProximityMatrix();
        void copy2Column(int _i);
        void copy2Row(int _i);
        void findLowestValue(int *_i,int *_j);
        void deleteRow(int _j);
        void deleteColumn(int _i);
        void joinColumns(int _i,int _j);
        int selectRandomSolution(int _whichcluster);
        int selectClusteredRandomSolution(int _whichcluster);
        double rnd(double _min,double _max);
        int nSolutions();
};

class MOPSO{
    vector<double> lb,ub;
    int ndimensions;
    int nobjectives;
    int nparticles;
    int gmax;
    vector<Particle> particles;
    double W,C1,C2;
    int nclusters;
    EFILE archive;
    int gen;
public:
    MOPSO(int _dims, int _objs,int _parts,int _nclusters,int _gmax,double *_lb,double *_ub );
    ~MOPSO();
    double rnd(double _min,double _max);
    void function(int _w);
    void initialize();
    int selectGBest(int _swarm);
    void flight();
    void execute();
    void copyx(vector<double> &_a,vector<double> &_b);
    void copyfx(vector<double> &_a,vector<double> &_b);
    void copy(vector<double> &_a,vector<double> &_b);
    void addClusters();
    int selectGbest(int _whichcluster);
    void perturbation(int );
};



/********************************************************************************************/


void EFILE::normalizeDistances(){

  for(int _i(0);_i<ndimensions;_i++)
    {
      xtlb[_i]=10*pow(10.0,10.0);
      xtub[_i]=-10*pow(10.0,10.0);
    }
  for(int _i(0);_i<nclusters;_i++){
    nsolutionsbyswarm[_i]=0;
  }

  for(int _j(0);_j<nsolutions;_j++){
    for(int _i(0);_i<ndimensions;_i++)
      {
	if(solutions[_j].x[_i] < xtlb[_i]) xtlb[_i]=solutions[_j].x[_i];
	if(solutions[_j].x[_i] > xtub[_i]) xtub[_i]=solutions[_j].x[_i];

      }
  }

  for(int _j(0);_j<nsolutions;_j++){
    for(int _i(0);_i<ndimensions;_i++)
      {
	solutions[_j].nx[_i]=(double)(solutions[_j].x[_i]-xtlb[_i]);
	if(abs(xtub[_i]-xtlb[_i])==0)solutions[_j].nx[_i]=0;
	else
	  solutions[_j].nx[_i]/=(double)(xtub[_i]-xtlb[_i]);
      }
  }

}




void EFILE::hierarchicalClustering()
{
  int _i,_j;
  if(nsolutions>nclusters+1){
    hasclusteredfirsttime=true;
    updatematrix=false;
    normalizeDistances();
    createProximityMatrix();
    while(pmsize>nclusters){
      findLowestValue(&_i,&_j);
      copy2Column(_i);
      copy2Column(_j);
      joinColumns(_i,_j);
      copy2Row(_i);
      deleteColumn(_j);
      deleteRow(_j);
      pmsize--;
    }
    for(_i=0;_i<nclusters;_i++){
      for(_j=0;proximitymatrixindex[_i][_j]!=-1;_j++);
      /*
	{

	for(int _k=0;_k<ndimensions;_k++){

	gbests4swarm[_i][_j].x[_k]=solutions[proximitymatrixindex[_i][_j]].x[_k];
      }

      for(int _k=0;_k<nobjectives;_k++){
      gbests4swarm[_i][_j].fx[_k]=solutions[proximitymatrixindex[_i][_j]].fx[_k];
      }
      }*/

      nsolutionsbyswarm[_i]=_j;
    }

  }
}

void EFILE::createProximityMatrix(){
  int _i(0),_j(0);
  for(int _k(0);_k<nsolutions-1;_k++,_i++){
    proximitymatrix[_i][_i]=0;
    _j=_i+1;
    for(int _l(_k+1);_l<nsolutions;_l++,_j++){
      proximitymatrix[_i][_j]=euclideanDistance(solutions[_k].nx,solutions[_l].nx);
    }
    pmsize=_i+1;
    for(int _i(0);_i<pmsize;_i++){
      proximitymatrixindex[_i][0]=_i;
      for(int _j(1);_j<pmsize;_j++){
	proximitymatrixindex[_i][_j]=-1;
      }
    }


  }
}
void EFILE::copy2Column(int _i){
  for(int _k(_i+1);_k<pmsize;_k++){
    proximitymatrix[_k][_i]=proximitymatrix[_i][_k];
  }
}

void EFILE::copy2Row(int _i){
  for(int _k(_i+1);_k<pmsize;_k++){
    proximitymatrix[_i][_k]=proximitymatrix[_k][_i];
  }
}

void EFILE::findLowestValue(int *_i,int *_j){
  double _maxvalue=10*pow(10.0,6.0);
  int _tempi=-1,_tempj=-1;
  for(int _k(0);_k<pmsize;_k++){
    for(int _l=_k+1;_l<pmsize;_l++){
      if(proximitymatrix[_k][_l]<_maxvalue){
	_tempi=_k;
	_tempj=_l;
	_maxvalue=proximitymatrix[_k][_l];
      }
    }
  }
  if(_tempi>_tempj)
    {
      int _tmp=_tempi;
      _tempi=_tempj;
      _tempj=_tmp;
    }
  *_i=_tempi;
  *_j=_tempj;
}

void EFILE::deleteRow(int _j){
  for(int _k(_j);_k+1<pmsize;_k++){
    for(int _l(0);_l<pmsize;_l++){//ojo hay que optimizar
      proximitymatrix[_k][_l]=proximitymatrix[_k+1][_l];
    }
    int _m(0);
    for(;proximitymatrixindex[_k+1][_m]!=-1;_m++)
      proximitymatrixindex[_k][_m]=proximitymatrixindex[_k+1][_m];
    proximitymatrixindex[_k][_m]=-1;
  }
}
void EFILE::deleteColumn(int _i){
  for(int _k(_i);_k+1<pmsize;_k++){
    for(int _l(0);_l<pmsize;_l++){//ojo hay que optimizar
      proximitymatrix[_l][_k]=proximitymatrix[_l][_k+1];
    }
  }
}

void EFILE::joinColumns(int _i,int _j){
  for(int _k=0;_k < pmsize;_k++)
    {
      double _ii=proximitymatrix[_k][_i];
      double _jj=proximitymatrix[_k][_j];

      proximitymatrix[_k][_i]=(_ii<_jj)?_ii:_jj;//si se toma el más cercano
    }
  int _l(0);
  for(;proximitymatrixindex[_i][_l]!=-1;_l++);
  int _k(0);
  for(;proximitymatrixindex[_j][_k]!=-1;_k++){
    proximitymatrixindex[_i][_l+_k]=proximitymatrixindex[_j][_k];
  }
  proximitymatrixindex[_i][_l+_k]=-1;
}




/********************************************************************************************/


/*

Algorithm
.- Begin
.- Hacer
.- Crear la matriz de eps.
.- Buscar el valor de eps más chico.
.- Eliminar el renglon y la columna del valor más pequenio.
.- Buscar actualizar el nuevo valor más pequenio para la dimensión eliminada.
.- Mientras el valor de ndimensions sea alcanzado.
.- End.





 */



EFILE::EFILE(){
  nsolutions=0;
  hasentered=false;
  firsttime=true;
  updatematrix=true;
}


void EFILE::init(int _ndimensions, int _nobjectives, int _nclusters, int _maxsolutions){
  nobjectives=_nobjectives;
  ndimensions=_ndimensions;
  nclusters=_nclusters;
  maxsolutions=_maxsolutions;
  updatematrix=true;
  hasclusteredfirsttime=false;
  EPS.resize(_nobjectives);

  tlb.resize(_nobjectives);
  tub.resize(_nobjectives);

  xtlb.resize(_ndimensions);
  xtub.resize(_ndimensions);

  nsolutions=0;
  for(int _i(0);_i<nobjectives;_i++){
    EPS[_i]=0.01;
  }

  for(int _i(0);_i<nobjectives;_i++)
    {
      tlb[_i]=10*pow(10.0,10.0);
      tub[_i]=-10*pow(10.0,10.0);
    }


  for(int _i(0);_i<ndimensions;_i++)
    {
      xtlb[_i]=10*pow(10.0,10.0);
      xtub[_i]=-10*pow(10.0,10.0);
    }

  nsolutionsbyswarm.resize(nclusters);
  for(int _i(0);_i<nclusters;_i++){
    nsolutionsbyswarm[_i]=0;
  }



  hasentered=false;
  firsttime=true;
}

void EFILE::update(){
    if(firsttime){
      normalizeDistances1();
      firsttime=false;
    }
    else normalizeDistances2();
    hasentered=true;
}


void EFILE::printSolutions(){
  for(int _i(0);_i<nsolutions;_i++){
    for(int _j(0);_j<nobjectives;_j++)
      cout<<solutions[_i].fx[_j]<<" ";

    cout<<endl;

  }
}

//
//  Funcion que indica si un individuo domina a otro
//  si el primero domina al segundo retorna 1
//  si el segundo domina al primero retoruna -1
//  si son iguales retorna 0
//  si ninguno domina retorna 11
//

int  EFILE::domine1(vector<double> &_a,vector<double> &_b){
  int anterior = 0, mejor;
  for(int i=0;i<nobjectives;i++){
    if(_a[i] <_b[i])	mejor = 1;
    else if(_b[i]<_a[i])mejor = -1;
    else mejor = 0;
    if(mejor!=anterior&&anterior!=0&&mejor!=0) return(11);
    if(mejor!=0) anterior = mejor;
  }
  return(anterior);
}

int EFILE::domine(vector<double> &_a,vector<double> &_b){
  int _tmp;
      if(!hasentered){
      _tmp=domine1(_a,_b);
      }
    else _tmp=domine2(_a,_b);

    return _tmp;
}


bool EFILE::add(Particle &_particle,int _whichcluster){
  bool _flag=true;

  vector<bool> erase;
  vector<int> _delparticles;

  addAll(_particle);

  for(int _i(0);_i<nsolutions;_i++){
    int _tmp=domine(solutions[_i].fx,_particle.fx);

    if(_tmp==1||_tmp==0) {
      _flag=false;
    }
    if(_tmp==-1){
      _delparticles.push_back(_i);
      erase.push_back(true);
    }
    else erase.push_back(false);
  }
  if(_delparticles.size()>0)
  updateClusters(_delparticles);

   vector<bool>::iterator _iterator=erase.begin();
   vector<Particle>::iterator _iterator2=solutions.begin();
   int _delindicator(0);
   for(_iterator=erase.begin(); _iterator != erase.end();++_iterator,_delindicator) {
     if(*_iterator==true){
       solutions.erase(_iterator2);
       nsolutions--;
     }
     else
       _iterator2++;
   }
  if(_flag==true){
    if(nsolutions>=maxsolutions+maxsolutions/10){
      update();
      add(_particle);
      return _flag;
    }
    else{
      if(_whichcluster!=-1)add2Cluster(_whichcluster);
      solutions.push_back(_particle);
      nsolutions++;
    }
  }

  if(nsolutions<=maxsolutions-maxsolutions/10){
    hasentered=false;
    firsttime=true;
  }
  //  update();

  return _flag;

}


void EFILE::add2Cluster(int _whichcluster){

  if(hasclusteredfirsttime==false)return;
  //  cout<<_whichcluster<<nsolutionsbyswarm[_whichcluster]<<endl;
  proximitymatrixindex[_whichcluster][nsolutionsbyswarm[_whichcluster]]=nsolutions;
  proximitymatrixindex[_whichcluster][nsolutionsbyswarm[_whichcluster]+1]=-1;
  nsolutionsbyswarm[_whichcluster]++;
}

int EFILE::deleteFromCluster(int _whichsolution){
  int _i;
  int _j;
  int _k;
  for(_j=0;_j<nclusters;_j++)
    for(_i=0;_i<nsolutionsbyswarm[_j];_i++)
      if(proximitymatrixindex[_j][_i]==_whichsolution){
	for(_k=_i;_k<nsolutionsbyswarm[_j];_k++){
	  proximitymatrixindex[_j][_k]=proximitymatrixindex[_j][_k+1];
	}
	proximitymatrixindex[_j][_k]=-1;
	nsolutionsbyswarm[_j]--;
	return 1;
      }
  return -1;
}

int EFILE::substractFromTo(int _from, int _to, int _quantity){
  int _i;
  int _j;
  for(_j=0;_j<nclusters;_j++){
    for(_i=0;_i<nsolutionsbyswarm[_j];_i++){
      if(proximitymatrixindex[_j][_i]>_from&&proximitymatrixindex[_j][_i]<_to){
	proximitymatrixindex[_j][_i]-=_quantity;
      }
    }
  }
  return -1;

}
int EFILE::updateClusters(vector<int> &_whichparticles){
  int _counter(0);
  if(hasclusteredfirsttime==false)return -1;
    _whichparticles.push_back(nsolutions);
    while(_whichparticles.size()>_counter+1){
      int _from(_whichparticles[_counter]);
      int _to(_whichparticles[_counter+1]);
       deleteFromCluster(_from);
       substractFromTo(_from,_to,_counter+1);
      _counter++;
    }

  return 1;
}


bool EFILE::addAll(Particle &_particle){
  bool _flag=true;
  vector<bool> erase;
  //  vector<EPSs>::iterator _iterator;

  for(int _i(0);_i<allnsolutions;_i++){
    int _tmp;
    _tmp=domine(allsolutions[_i].fx,_particle.fx);
    if(_tmp==1||_tmp==0) _flag=false;

    if(_tmp==-1){
      erase.push_back(true);
    }
    else erase.push_back(false);
  }

   vector<bool>::iterator _iterator=erase.begin();
   vector<Particle>::iterator _iterator2=allsolutions.begin();
   for(_iterator=erase.begin(); _iterator != erase.end();++_iterator) {
     if(*_iterator==true){
       allsolutions.erase(_iterator2);
       allnsolutions--;
     }
     else
       _iterator2++;
   }

  if(_flag==true){
    if(allnsolutions>maxsolutions*3)return _flag;
      allsolutions.push_back(_particle);
      allnsolutions++;
  }
  return _flag;

}


bool EFILE::reAdd(){
  bool _flag=true;
  vector<bool> erase;
  vector<int> _delparticles;
  erase.resize(nsolutions);

  //  vector<EPSs>::iterator _iterator;
  for(int _i(0);_i<nsolutions;_i++){
    erase[_i]=false;
    for(int _j(0);_j<nsolutions;_j++){
      int _tmp;
      if(_i==_j)continue;
      _tmp=domine2(solutions[_i].fx,solutions[_j].fx);

      if(_tmp==1||_tmp==0) _flag=false;

      if(_tmp==-1){
      _delparticles.push_back(_i);
	erase[_i]=true;
      }
    }
  }
  //  if(_delparticles.size()>0)
  // updateClusters(_delparticles);

   vector<bool>::iterator _iterator=erase.begin();
   vector<Particle>::iterator _iterator2=solutions.begin();
   for(_iterator=erase.begin(); _iterator != erase.end();++_iterator) {
     if(*_iterator==true){
       solutions.erase(_iterator2);
       nsolutions--;
     }
     else
       _iterator2++;
   }

  return _flag;

}



int EFILE::falseReAdd(){
  int _howmany(0);
  vector<bool> erase;
  erase.resize(nsolutions);

  //  vector<EPSs>::iterator _iterator;
  for(int _i(0);_i<nsolutions;_i++){
    erase[_i]=false;
    for(int _j(0);_j<nsolutions;_j++){
      int _tmp;
      if(_i==_j)continue;
      _tmp=domine2(solutions[_i].fx,solutions[_j].fx);
      if(_tmp==-1)
	erase[_i]=true;
    }
  }
   vector<bool>::iterator _iterator=erase.begin();
   vector<Particle>::iterator _iterator2=solutions.begin();
   for(_iterator=erase.begin(); _iterator != erase.end();++_iterator) {
     if(*_iterator==true){
       _howmany++;
     }
     else
       _iterator2++;
   }
  return _howmany;
}


int EFILE::domine2(vector<double> &_a,vector<double>&_b){
  vector<double> _box(nobjectives),_box2(nobjectives),_box3(nobjectives);
  bool _flag=true;
  //calculate the box of both particles
  for (int _i = 0; _i < nobjectives; _i++){
    _box[_i] = (int) floor ((fabs(tlb[_i]-_a[_i]) / EPS[_i]));
    _box2[_i] = (int)floor ((fabs(tlb[_i]-_b[_i]) / EPS[_i]));
    //

    //_box[_i] = (int) floor ((_a[_i] / EPS[_i]));
    //_box2[_i] = (int)floor ((_b[_i] / EPS[_i]));
    //_box3[_i] = (_box[_i]<_box2[_i])?_box[_i]*EPS[_i]:_box2[_i]*EPS[_i];//
    _box3[_i] = (int)_box[_i]*EPS[_i];
    //if they are in the same box
    if(_box[_i]!=_box2[_i])_flag=false;
  }
  if(_flag==true){//check for dominance
    int anterior = 0, mejor;

    for(int _i=0;_i<nobjectives;_i++){
      if(_a[_i] <_b[_i])	mejor = 1;
      else if(_b[_i]<_a[_i])mejor = -1;
      else mejor = 0;
      if(mejor!=anterior&&anterior!=0&&mejor!=0){
	if(euclideanDistance(_a,_box3)<euclideanDistance(_b,_box3))
	  return 1;
	else return -1;	}
      if(mejor!=0) anterior = mejor;
    }
    //      if(anterior==1) return true;
    //else return false;
    return(anterior);

  }
  int anterior = 0, mejor;
  for(int _i=0;_i<nobjectives;_i++){
    if(_box[_i] <_box2[_i])	mejor = 1;
    else if(_box2[_i]<_box[_i])mejor = -1;
    else mejor = 0;
    if(mejor!=anterior&&anterior!=0&&mejor!=0)return 11;
    if(mejor!=0) anterior = mejor;
  }
  return(anterior);


}

void EFILE::normalizeDistances2(){

  for(int _j(0);_j<nsolutions;_j++){
    for(int _i(0);_i<nobjectives;_i++)
      {
	if(solutions[_j].fx[_i] < tlb[_i]) tlb[_i]=solutions[_j].fx[_i];
	if(solutions[_j].fx[_i] > tub[_i]) tub[_i]=solutions[_j].fx[_i];

      }
  }

  double min;
  int _index;
  int _ncurrentsolutions=nsolutions;
  int _a1(0);
  int _a2(0);
  while(_ncurrentsolutions>maxsolutions+maxsolutions/20||_ncurrentsolutions<maxsolutions-maxsolutions/10){
    _ncurrentsolutions=nsolutions-falseReAdd();

    if(_ncurrentsolutions<maxsolutions-maxsolutions/10){
      for(int _i(0);_i<nobjectives;_i++){
	EPS[_i]/=1.2;
      }
	if(_a1==2&&_a2==1) break;
	_a2=_a1;
	_a1=1;

    }
    else if(_ncurrentsolutions>maxsolutions+maxsolutions/20){
      for(int _i(0);_i<nobjectives;_i++){
	EPS[_i]*=1.2;
      }
	if(_a1==1&&_a2==2) break;
	_a2=_a1;
	_a1=2;
    }


  }
  reAdd();

}


void EFILE::normalizeDistances1(){

  for(int _j(0);_j<nsolutions;_j++){
    for(int _i(0);_i<nobjectives;_i++)
      {
	if(solutions[_j].fx[_i] < tlb[_i]) tlb[_i]=solutions[_j].fx[_i];
	if(solutions[_j].fx[_i] > tub[_i]) tub[_i]=solutions[_j].fx[_i];

      }
  }
  for(int _i(0);_i<nobjectives;_i++){
    EPS[_i]=(double)(tub[_i]-tlb[_i])/nsolutions;
  }

  double min;
  int _index;
  int _a1(0),_a2(0);
  int _ncurrentsolutions=nsolutions;

  while(_ncurrentsolutions>maxsolutions+maxsolutions/20||_ncurrentsolutions<maxsolutions-maxsolutions/10){
    _ncurrentsolutions=nsolutions-falseReAdd();

    if(_ncurrentsolutions<maxsolutions-maxsolutions/10){
      for(int _i(0);_i<nobjectives;_i++){
	EPS[_i]/=1.2;
      }
	if(_a1==2&&_a2==1) break;
	_a2=_a1;
	_a1=1;
    }
    else if(_ncurrentsolutions>maxsolutions+maxsolutions/20){
      for(int _i(0);_i<nobjectives;_i++){
	EPS[_i]*=1.2;
      }
	if(_a1==1&&_a2==2) break;
	_a2=_a1;
	_a1=2;
    }

  }
  reAdd();
}






void EFILE::deleteExcedent(){

  for(int _j(0);_j<nsolutions;_j++){
    for(int _i(0);_i<nobjectives;_i++)
      {
	if(solutions[_j].fx[_i] < tlb[_i]) tlb[_i]=solutions[_j].fx[_i];
	if(solutions[_j].fx[_i] > tub[_i]) tub[_i]=solutions[_j].fx[_i];

      }
  }

  vector<double> distances;
  vector<int> _indexes;
  vector<double> distances2;
  vector<int> _indexes2
;
  while(nsolutions>maxsolutions){

    _indexes.resize(nsolutions);

    distances.resize(nsolutions);
    _indexes2.resize(nsolutions);

    distances2.resize(nsolutions);

    for(int _j(0);_j<nsolutions;_j++){
      distances[_j]=10.0*pow(10.0,10.0);
      for(int _k(0);_k<nsolutions;_k++){
	if(_j==_k)continue;
	double _sum(0);
	for(int _i(0);_i<nobjectives;_i++)
	  {
	    _sum+=pow(solutions[_j].fx[_i]-solutions[_k].fx[_i],2.0);
	  }
	double _tot=sqrt(_sum);

	if(_tot<distances[_j]){
	  distances2[_j]=distances[_j];
	  distances[_j]=_tot;
	  _indexes2[_j]=_indexes[_j];
	  _indexes[_j]=_k;
	}
      }
    }

  double min;
  int _index;
  min=10.0*pow(10.0,10.0);
  for(int _j(0);_j<nsolutions;_j++){
    if(distances[_j]+distances2[_j]<min){
      min=distances[_j]+distances2[_j];
      _index=_j;
    }
  }

   vector<Particle>::iterator _iterator2=solutions.begin();
   _iterator2+=_index;
   solutions.erase(_iterator2);
   nsolutions--;
  }
}


void EFILE::finalSolutions(){
  for(int _i(0);_i<allnsolutions;_i++){
    add(allsolutions[_i]);
  }
  //cout<<"deleteExcedent"<<endl;
  deleteExcedent();
}

double EFILE::euclideanDistance(vector<double> &_v1, vector<double> & _v2){
  double s(0);
  vector<double>::iterator _iterator1=_v1.begin();
  vector<double>::iterator _iterator2=_v2.begin();
  for(_iterator1=_v1.begin(); _iterator1 != _v1.end();++_iterator1) {
    s+= pow(*_iterator2-*_iterator1,2.0);
  }
  return sqrt(s);
}


int EFILE::nSolutions(){
  return nsolutions;
}

int EFILE::selectRandomSolution(int _whichcluster){
    return (int) rnd(0,nsolutions);
}

int EFILE::selectClusteredRandomSolution(int _whichcluster){
  if(_whichcluster>nclusters) _whichcluster=(int) rnd(0,nclusters);//ojo hay que verificar
  if(nsolutionsbyswarm[_whichcluster]==0){
    return(int) rnd(0,nsolutions);
  }
  int _tmp=(int )rnd(0,nsolutionsbyswarm[_whichcluster]);

  return proximitymatrixindex[_whichcluster][_tmp];
}

int EFILE::nClusteredSolutions(int _whichcluster){
  return nsolutionsbyswarm[_whichcluster];
}
int EFILE::selectClusteredSolution(int _whichcluster,int _whichsolution){
  return proximitymatrixindex[_whichcluster][_whichsolution];
}




double EFILE::rnd(double _min,double _max){
  return((double)(_min + ((double)(_max-_min)*rand()/(double)(RAND_MAX+_min))));
}

void EFILE::output(){
  finalSolutions();
  printSolutions();
}



