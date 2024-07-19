/***************************************************************************
                          mainlib.h  -  description
                             -------------------
    begin                : Sun Jan 28 2001
    copyright            : (C) 2001 by Max Salazar
    email                : max@maxnet.cc
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
using namespace std;

void vuelo(unsigned int fun, unsigned int M, unsigned int Gmax, unsigned int opt, unsigned int D, unsigned int NF, unsigned int MEM, unsigned int ndiv, char *cad) {

  unsigned int h, i, j, t, rest, Gbestpos, gbestpos, GT, lastPart = 0, maxCube = (unsigned int) pow(ndiv, NF);
  double Gbest, gbest, x, passed;
  double *pop, *pbests, *vel, *fitness, *fbests;
  double *noDomP, *noDomF;
  double *amp, *start;
  unsigned int  *linf, *lsup, *hyperspace, *partPos;
  unsigned int *hyperPoolL, poolC;
  double *hyperPoolF, hypsum;
  hyperPoolL = new unsigned int[MEM];
  hyperPoolF = new double[MEM];
  clock_t  now, later;
  //MPTR noDom;
  //DSPTR aux;
  //DPTR auxV;
  pop = new double [M * D];
  pbests = new double [M * D];
  vel = new double [M * D];
  fitness = new double [M * NF];
  fbests = new double [M * NF];
  noDomP = new double [MEM * D];
  noDomF = new double[MEM * NF];
  ///
  amp = new double [NF];
  start = new double [NF];
  linf = new unsigned int[NF];
  lsup = new unsigned int[NF];
  hyperspace = new unsigned int[maxCube];
  partPos = new unsigned int[MEM];
  //cambio por gtp 25 junio 2002
  unsigned int *selec= new unsigned int[80/10];
  bool state=false;
  //fin cambio
  //cout << "Inicializacion..." << endl;

  ///////////////////////////////
  //Inicializacion de Variables//
  ///////////////////////////////
  /*
    Handle the seed range erraaors
    i = First random number seed must be between 0 and 31328
    j = Second seed must have a value between 0 and 30081
  */
  srand((unsigned int)time((time_t *)NULL));
  i = (unsigned int) (31329.0 * rand() / (RAND_MAX + 1.0));
  j = (unsigned int) (30082.0 * rand() / (RAND_MAX + 1.0));
  //i = 1;
  //j = 100;
  RandomInitialise(i,j);

  now = clock();
  //Numero de restricciones
  rest = num_rest(fun);
  //Numero de Generaciones
  t = 0;
  //Se generan aleatoriamente las particulas
  initpop(pop, fun, M, D);
  //Se inicializa la velocidad de cada particula
  for(i = 0; i < M; i++)
    for(j = 0; j < D; j++)
      vel[(D * i) + j] = 0.0;
  //Se evalua cada una de las particulas
  evaluation(pop, fitness, M, D, fun, NF);
  //Copia los valores en memoria
  copy_array(fitness, fbests, M * NF);
  copy_array(pop, pbests, M * D);
  if (NF == 1) {
    //Mejor global
    gbestpos = search_opt_cons(fitness, M, D, rest, fun, pop, opt);
    gbest = fitness[gbestpos];
    Gbestpos = gbestpos;
    Gbest = fitness[Gbestpos];
    x = Gbest;
    GT = 0;
    cout << "Mejor particula al inicio: " << gbest << endl;
    cout << "Mejor particula Global inicial: " << Gbest<< endl;
    cout << "----------------------------------------------" << endl;
  }
  else {
    //Busca a todas las partículas que no sean dominadas y las inserta al repositorio.
    search_insert(noDomP, noDomF, pop, fitness, D, NF, M, MEM, opt, &lastPart,fun);
    //Se divide el espacio de búsqueda explorado hasta ahora y se investiga la posición de las partículas en el.
    makeHyper(noDomF, linf, lsup, amp, start, hyperspace, partPos, ndiv, MEM, NF, lastPart, maxCube);
  }
  //Ciclo de vuelo
  do {
    //Se obtienen las nuevas velocidades de cada particula
    if (NF == 1) {
      for(i = 0; i < M; i++)
	for(j = 0; j < D; j++)
	  vel[(D * i) + j] = velocity(0.4, vel[(D * i) + j], 2.0, 2.0, pbests[(D * i) + j], pbests[(D * Gbestpos) + j], pop[(D * i) + j]);
    }
    else {
      //h = 0;
      hyperFit(hyperspace, hyperPoolL, hyperPoolF, &hypsum, maxCube, &poolC);
      state=false;
      for(i = 0; i < M; i++) {
	if(state==false)state=true;
	//h=lead2(pop,lastPart,selec,state,i,noDomF,NF);
	h = lead(hyperspace, hyperPoolL, partPos, hyperPoolF, hypsum, lastPart, poolC);
	for(j = 0; j < D; j++)
	  vel[(D * i) + j] = velocity(0.4, vel[(D * i) + j], 1.0, 1.0, pbests[(D * i) + j], noDomP[(D * h) + j], pop[(D * i) + j]);
	//h += 1;
	//if (h > lastPart - 1)
	//	h = 0;
      }
    }
    //Se calculan las nuevas posiciones
    for(i = 0; i < M; i++)
      for(j = 0; j < D; j++)
	pop[(D * i) + j] =	pop[(D * i) + j] + vel[(D * i) + j];

    //Se mantinen las particulas dentro del espacio
    keepin(pop, vel, fun, M, D);
    //Muta cada una de las particulas
    if(t<Gmax*pM){
      //cout<<t<<endl;
      muta(pop, M, D,t,Gmax,fun);
    }
    //Se evalua cada una de las particulas
    evaluation(pop, fitness, M, D, fun, NF);

    if (NF == 1) {
      //Se actualizan las mejores posiciones de cada particula
      ifbest_interchange(fitness, pop, fbests, pbests, M, D, fun, rest, opt);

      //Se obtiene al mejor del ciclo
      gbestpos = search_opt_cons(fitness, M, D, rest, fun, pop, opt);
      gbest = fitness[gbestpos];

      //Se busca al mejor lider
      Gbestpos = search_opt_cons(fbests, M, D, rest, fun, pbests, opt);
      Gbest = fbests[Gbestpos];
      if (Gbest != x) {
	GT = t + 1;
	x = Gbest;
      }

      cout << "Ciclo #"<< t + 1 << endl;
      cout << "Mejor particula del ciclo: " << gbest << endl;
      cout << "Mejor particula del ciclo >: " << chfun(&pop[D * gbestpos], fun) << endl;
      cout << "Mejor particula Global: " << Gbest << endl;
      cout << "Mejor particula Global >: " << chfun(&pbests[D * Gbestpos], fun) << endl;
      cout << "------" << endl;
      cout << "----------------------------------------------" << endl;
    }
    else {
      search_insert2(noDomP, noDomF, pop, fitness, linf, lsup, amp, start, hyperspace, partPos, ndiv, D, NF, M, MEM, opt, maxCube, &lastPart,fun);
	ifbest_interchange_moo(fitness, pop, fbests, pbests, M, D, NF, opt,fun);
    }
    //Se actualiza el contador de generaciones
    t++;
  } while (Gmax-1 > t);
  //cout << "Hola1" << endl;
  if (NF == 1) {
    //Mejor Particula Global
    cout << endl << "Mejor particula Global: " << Gbest << endl;
    cout << "Mejor particula Global >: " << chfun(&pbests[D * Gbestpos], fun) << endl;
    cout << "Puntos: " << endl;
    for(i = 0; i < D; ++i) {
      cout << "x"<< i+1 << ": " << pbests[(D * Gbestpos) + i] << endl;
    }
    cout << "Generacion: " << GT << endl;
    cout << "Restricciones: " << endl;
    for(i = 0; i < rest; ++i) {
      cout << "g"<< i+1 << ": " << chconst(&pbests[D * Gbestpos], fun, i + 1) << endl;
    }
  }
  else {
    //cout << CONT_FUN / 2 << endl;
    mo_out(noDomF, lastPart, NF, cad);
  }
  later = clock();
  passed = ( later - now ) / (double)CLOCKS_PER_SEC;
  //  printf("%lf\n", passed);
  //cout << "Hola2" << endl;
  //delete [] noDomF;
  //delete [] noDomP;
  delete [] pop;
  delete [] pbests;
  delete [] vel;
  delete [] fitness;
  delete [] fbests;
  delete [] noDomP;
  delete [] noDomF;
  delete [] amp;
  delete [] start;
  delete [] linf;
  delete [] lsup;
  delete [] hyperspace;
  delete [] partPos;
  delete [] hyperPoolL;
  delete [] hyperPoolF;
  //cambio gtp 25 junio 2002
  delete [] selec;
  //end cambio
  //cout << "Hola3" << endl;
}
