/***************************************************************************
                          psolib.h  -  description
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

//Copia el arreglo s a d (De tipo double y longitud l)
void copy_array(double *s, double *d, unsigned int l) {
  unsigned int i;

  for(i = 0; i < l; i++)
    d[i] = s[i];
	
}

//Busca el valor maximo de un arreglo de tipo double
unsigned int search_max(double *f, unsigned int M) {
  unsigned int i, fmax = 0;

  for(i = 1; i < M; i++)
    if (f[fmax] < f[i])
      fmax = i;
	
  return fmax;
}

//Busca el valor maximo de un arreglo de tipo entero sin signo
unsigned int search_max(unsigned int *f, unsigned int M) {
  unsigned int i, fmax = 0;

  for(i = 1; i < M; i++)
    if (f[fmax] < f[i])
      fmax = i;
	
  return fmax;
}

//Busca el valor minimo de un arreglo de tipo double
unsigned int search_min(double *f, unsigned int M) {
  unsigned int i, fmin = 0;

  for(i = 1; i < M; i++)
    if (f[fmin] > f[i])
      fmin = i;
	
  return fmin;
}

//Busca el valor minimo de un arreglo de tipo entero sin signo
unsigned int search_min(unsigned int *f, unsigned int M) {
  unsigned int i, fmin = 0;

  for(i = 1; i < M; i++)
    if (f[fmin] > f[i])
      fmin = i;
	
  return fmin;
}
//*
void muta(double *pop, int M, unsigned int D,int Gen,int totGen,int fun){
  unsigned int i,j;
  int dimension=0;
  double linf[D],lsup[D],linftemp,lsuptemp,rango;
  double valtemp=0;
  ranges(fun,D,linf,lsup);
  for(i = 0;i < M;i++){
    if(flip(pow(1.0-(double)Gen/(totGen*pM),1.5))){
      dimension=RandomInt(0,D-1);
      rango=(lsup[dimension]-linf[dimension])*pow(1.0-(double)Gen/(totGen*pM),1.5)/2;//totGen
      valtemp=RandomDouble(rango,-rango);
      if(pop[dimension+(D*i)]-rango<linf[dimension]) linftemp=linf[dimension]; 
      else linftemp=pop[dimension+(D*i)]-rango;
      if(pop[dimension+(D*i)]+rango>lsup[dimension]) lsuptemp=lsup[dimension]; 
      else lsuptemp=pop[dimension+(D*i)]+rango;

      pop[dimension+(D*i)]=RandomDouble(linftemp,lsuptemp);//mutation_operator(pop[ngen+(D*i)],Gen);

    }
  }
}

//Se evaluan las particulas
void evaluation(double *pop, double *fitness, unsigned int M, unsigned int D, unsigned int fun, unsigned int NF) {
  unsigned int i, j, f;
		
  for(i = 0; i < M; i++)
    for(j = 0, f = fun; j < NF; j++, f++)
      fitness[(NF * i) + j] = chfun(&pop[D * i], f);
}

//Calcula la velocidad de las particulas
double velocity(double W, double Vi, double C1, double C2, double Pb, double Pg, double Xi) {
  return W * Vi + C1 * RandomDouble(0.0, 1.0) * (Pb - Xi)
    + C2 * RandomDouble(0.0, 1.0) * (Pg - Xi);
}

//Funcion que regresa 1 en caso de que el indiviudo (pop) a
//evaluar sea mejor que el establecido (pbests), de lo contrario
//la funcion devuelve 0.
unsigned int g_is_best(double *pop, double *pbests, unsigned int fun, unsigned int rest, double fitness, double fbests, unsigned int opt) {
  unsigned int i, j;
  double cont1, cont2, aux;
	
  i = 0;
  cont1 = 0;
  cont2 = 0;
 		
  for(j = 0; j < rest; j++)
    aux = chconst(pop, fun, j + 1);
  if (aux > 0.0)				//<<<<----Cuando las restricciones deben ser menores o iguales a 0 (>)
    cont1 += aux;
 		
  for(j = 0; j < rest; j++)
    aux = chconst(pbests, fun, j + 1);
  if (aux > 0.0)				//<<<<----Cuando las restricciones deben ser menores o iguales a 0 (>)
    cont2 += aux;
 					
  if ((cont1 == 0) && (cont2 == 0)) {
    if (((!(opt)) && (fitness < fbests)) || ((opt) && (fitness > fbests)))	// <<<------ min(<) - max(>)
      i = 1;
  }
  else {
    if ((cont1 == 0) || (cont2 == 0)) {
      if (cont1 == 0)
	i = 1;
    }
    else {
      if (cont1 < cont2)
	i = 1;
      else
	if ((cont1 == cont2) && (((!(opt)) && (fitness < fbests)) || ((opt) && (fitness > fbests))))	// <<<------ min(<) - max(>)
	  i = 1;
    }
  }
	
  return i;
}

//Intercambia los registros de los arreglos en
//caso de que se trate de una mejor particula.
void ifbest_interchange(double *fitness, double *pop, double *fbests, double *pbests, unsigned int M, unsigned int D, unsigned int fun, unsigned int rest, unsigned int opt) {
  unsigned int i, j;
	
  for(i = 0; i < M; i++) {
    if (g_is_best(&pop[D * i], &pbests[D * i], fun, rest, fitness[i], fbests[i], opt)) {
      fbests[i] = fitness[i];
      for(j = 0; j < D; j++)
	pbests[(D * i) + j] = pop[(D * i) + j];
    }
  }
	
	
}

//Busca a la mejor de las particulas, intentando
//que cumpla con todas sus restricciones
unsigned int search_opt_cons(double *f, unsigned int M, unsigned int D, unsigned int rest, unsigned int fun, double *pop, unsigned int opt) {
  unsigned int i, j, band = 0, cmin = 0, cmax = 0, it = 0, bestp;
  double aux, *cons;
  cons = new double [M];
	
  for(i = 0; i < M; i++) {
    cons[i] = 0;
    for(j = 0; j < rest; j++) {
      aux = chconst(&pop[D * i], fun, j + 1);
      if (aux > 0.0)				//<<<<----Cuando las restricciones deben ser menores o iguales a 0 (>)
	cons[i] += aux;
    }
  }
	
  cmax = search_max(cons, M);
  cmin = search_min(cons, M);
	
  bestp = cmin;
	
  if(cons[bestp] == 0.0)
    band = 1;
		
  do {
    if (band)
      if (cons[cmin] == 0.0) {
	if (((!(opt)) && (f[cmin] < f[bestp])) || ((opt) && (f[cmin] > f[bestp]))) {	// <<<------ min(<) - max(>)
	  bestp = cmin;
	}
	cons[cmin] = cons[cmax] + 1.0;
	cmin = search_min(cons, M);
      }
    band = 0;
				
    it++;
  } while ((band) && (it < M));
	
  delete [] cons;

  return bestp;	
}


///////////////////////////////////////////////////////
//MULTI-OBJETIVO
///////////////////////////////////////////////////////
/*
  unsigned int sumcompVec(double *noDomF, unsigned int NF, unsigned int h, unsigned int i, unsigned int opt) {
  unsigned int sum = 0, j;
	
  for(j = 0; j < NF; j++)
  if (((noDomF[(NF * h) + j] < noDomF[(NF * i) + j]) && (opt == 0)) || ((noDomF[(NF * h) + j] > noDomF[(NF * i) + j]) && (opt == 1)))
  //if (p->next == NULL)
  sum += 1;
  //else
  //	sum = sumcompVec(p->next, q->next, sum += 1, opt);
  //else
  //if (p->next != NULL)
  //	sum = sumcompVec(p->next, q->next, sum, opt);
				
  return sum;
  }
*/

unsigned int best_moo(double *fitness, double *fbests, unsigned int NF, unsigned int opt,unsigned int fun,double *pop) {
  unsigned int i, sum, ret;
	
  //  if(constraints(pop,fun) < constraints(fbests,fun))return 0;
  //if(constraints(pop,fun) > constraints(fbests,fun))return 1;
  //  cout<<"                        holas"<<endl;
  for(i = 0; i < NF; i++)
    if (((fitness[i] < fbests[i]) && (opt == 0)) || ((fitness[i] > fbests[i]) && (opt == 1)))
      sum += 1;
	
  if (sum == NF)
    ret = 0;
  else
    if (sum == 0)
      ret = 1;
    else
      ret = RandomInt(0, 1);
		
  return ret;
}

void ifbest_interchange_moo(double *fitness, double *pop, double *fbests, double *pbests, unsigned int M, unsigned int D, unsigned int NF, unsigned int opt,unsigned int fun) {
  unsigned int i, j;
		
  for(i = 0; i < M; i++) {
    if (best_moo(&fitness[NF * i], &fbests[NF * i], NF, opt,fun,pop+(NF*i))) {
      for(j = 0; j < NF; j++)
	fbests[(NF * i) + j] = fitness[(NF * i) + j];
      for(j = 0; j < D; j++)
	pbests[(D * i) + j] = pop[(D * i) + j];
    }
  }
}

unsigned int compVec(double *noDomF, unsigned int NF, unsigned int h, unsigned int i, unsigned int opt,unsigned int fun,double *noDomP,unsigned int D) {
  unsigned int sum = 0, ret, j;
	
  //sum = sumcompVec(noDomF, NF, h, i, opt);
    if(constraints(noDomF+(NF*h),fun)>0&&constraints(noDomF+(NF*i),fun)>0){
   return 4;
   }
  else if(constraints(noDomF+(NF*h),fun)>0)return 2;
  else if(constraints(noDomF+(NF*i),fun)>0)return 1;
  //if(constraints(noDomP+(D*h),fun) > constraints(pop+(D*i),fun))return 1;
  //if(constraints(noDomP+(D*h),fun) < constraints(pop+(D*i),fun))return 2;

  //  cout<<"                                          rakam"<<endl;

  for(j = 0; j < NF; j++)
    if (((noDomF[(NF * h) + j] < noDomF[(NF * i) + j]) && (opt == 0)) || ((noDomF[(NF * h) + j] > noDomF[(NF * i) + j]) && (opt == 1)))
      sum += 1;
	
  if (sum == NF)
    ret = 1;
  else
    if (sum == 0)
      ret = 2;
    else
      ret = 3;
			
  return ret;		
}
	
void delPartDom(double * noDomP, double *noDomF, unsigned int *lastPart, unsigned int D, unsigned int NF, unsigned int opt, unsigned int fun) {
  unsigned int h, i, j;// k;
  h = 0;
  i = h + 1;
  do {
    do {
      switch (compVec(noDomF, NF, h, i, opt,fun,noDomP,D)) {
      case 1: for(j = 0; j < D; j++)
	noDomP[(D * i) + j] = noDomP[(D * (*lastPart - 1)) + j];
      for(j = 0; j < NF; j++)
	noDomF[(NF * i) + j] = noDomF[(NF * (*lastPart - 1)) + j];/*
								    for(k = i; k < (*lastPart -1); k++) {
								    for(j = 0; j < D; j++)
								    noDomP[(D * k) + j] = noDomP[(D * (k + 1)) + j];
								    for(j = 0; j < NF; j++)
								    noDomF[(NF * k) + j] = noDomF[(NF * (k + 1)) + j];
								    }*/
      *lastPart -= 1;
      break;
      case 2:	 for(j = 0; j < D; j++)
	noDomP[(D * h) + j] = noDomP[(D * (*lastPart - 1)) + j];
      for(j = 0; j < NF; j++)
	noDomF[(NF * h) + j] = noDomF[(NF * (*lastPart - 1)) + j];/*
								    for(k = h; k < (*lastPart -1); k++) {
								    for(j = 0; j < D; j++)
								    noDomP[(D * k) + j] = noDomP[(D * (k + 1)) + j];
								    for(j = 0; j < NF; j++)
								    noDomF[(NF * k) + j] = noDomF[(NF * (k + 1)) + j];
								    }*/
      *lastPart -= 1;

      i = h + 1;
      break;
      case 3:	i += 1;
	break;
      case 4:
	noDomP[(D * i) + j] = noDomP[(D * (*lastPart - 1)) + j];
      for(j = 0; j < NF; j++)
	noDomF[(NF * i) + j] = noDomF[(NF * (*lastPart - 1)) + j];/*
								    for(k = i; k < (*lastPart -1); k++) {
								    for(j = 0; j < D; j++)
								    noDomP[(D * k) + j] = noDomP[(D * (k + 1)) + j];
								    for(j = 0; j < NF; j++)
								    noDomF[(NF * k) + j] = noDomF[(NF * (k + 1)) + j];
								    }*/
      *lastPart -= 1;

	noDomP[(D * h) + j] = noDomP[(D * (*lastPart - 1)) + j];
      for(j = 0; j < NF; j++)
	noDomF[(NF * h) + j] = noDomF[(NF * (*lastPart - 1)) + j];/*
								    for(k = h; k < (*lastPart -1); k++) {
								    for(j = 0; j < D; j++)
								    noDomP[(D * k) + j] = noDomP[(D * (k + 1)) + j];
								    for(j = 0; j < NF; j++)
								    noDomF[(NF * k) + j] = noDomF[(NF * (k + 1)) + j];
								    }*/
      *lastPart -= 1;

      i = h + 1;

	break;

      }
    } while (i < *lastPart);
    h += 1;
    i = h + 1;
  } while (h < *lastPart - 1);
}

//Imprime la salida
void mo_out(double *noDomF, unsigned int lastPart, unsigned int NF, char *cad) {
  unsigned int i, j;
  FILE *salF;
	
  salF = fopen(cad,"w");
	
  for(i = 0; i < lastPart; i++) {
    for(j = 0; j < NF; j++)
      fprintf(salF, "%.4f ", (noDomF[(NF * i) + j]));
    fprintf(salF, "\n");
  }
	
  fclose(salF);
}

void search_insert(double *noDomP, double *noDomF, double *pop, double *fitness, unsigned int D, unsigned int NF, unsigned int M, unsigned int MEM, unsigned int opt, unsigned int *lastPart,unsigned int fun) {
  unsigned int i, j, k;

  for(i = *lastPart, k = 0; k < M; i++, k++) {
    for(j = 0; j < D; j++)
      noDomP[(D * i) + j] = pop[(D * k) + j];
    for(j = 0; j < NF; j++) {
      noDomF[(NF * i) + j] = fitness[(NF * k) + j];
    }
    *lastPart += 1;
  }
	
  //elimina las partículas dominadas
  delPartDom(noDomP, noDomF, lastPart, D, NF, opt,fun);

}

unsigned int locHyper(double *vect, double *start, double *amp, unsigned int ndiv, unsigned int NF, unsigned int part) {
  unsigned int temp, loc = 0;
  unsigned int i = NF - 1, j;

  for(j = 0; j < NF ; j++ , i--) {
    temp = (long int) floor((vect[(NF * part) + j] - start[j]) / amp[j]);
    if (temp > ndiv || temp < 0) return 0;
    loc += (long int) floor(temp * pow(ndiv, i));
  }
  return loc + 1;
}

void makeHyper(double *noDomF, unsigned int *linf, unsigned int *lsup, double *amp, double *start, unsigned int *hyperspace, unsigned int *partPos, unsigned int ndiv, unsigned int MEM, unsigned int NF, unsigned int lastPart, unsigned int maxCube) {
  double *aux;
  unsigned int i, j;
  aux = new double [MEM];

  for(j = 0; j < NF; j++) {
    for(i = 0; i < lastPart; i++)
      aux[i] = noDomF[(NF * i) + j];
    linf[j] = search_min(aux, i);
    lsup[j] = search_max(aux, i);
  }

  for(i = 0; i < NF; i++) {
    amp[i] = (noDomF[(NF * lsup[i]) + i] - noDomF[(NF * linf[i]) + i]) / (double)(ndiv - 1);
    start[i] = noDomF[(NF * linf[i]) + i] - (amp[i] / 2.0);
  }

  for(i = 0; i < maxCube; i++)
    hyperspace[i] = 0;

  for(i = 0; i < lastPart; i++) {
    partPos[i] = locHyper(noDomF, start, amp, ndiv, NF, i) - 1;
    hyperspace[partPos[i]] += 1;
  }
	
  delete [] aux;

  /*
    Monumento a "Lo que hace la ignorancia":
    //Se genera el primer hipercubo
    //Punto 1
    for(i = 0; i < NF; i++)
    grid[(NF * 0) + i] = start[i];
    //Punto 2
    for(i = 0; i < NF; i++)
    grid[(NF * 1) + i] = start[i] + amp[i];

    //El número de hipercubos esta dado por el número de divisiones^número de funciones
    //Se forman los primeros h hipercubos en base al primero (o a su antesesor)
	
	
    //Viejo enfoque
    for(h = 1; h < ndiv; h++)
    for(i = 0; i < Point; i++)
    for(j = 0; j < NF; j++)
    if (j == i)
    grid[((NF * Point) * h) + (NF * i) + j] = grid[((NF * Point) * (h - 1)) + (NF * i) + j] + amp[j];
    else
    grid[((NF * Point) * h) + (NF * i) + j] = grid[((NF * Point) * (h - 1)) + (NF * i) + j];

    //Se generan los restantes hipercubos en base a los primeros h (o a sus antecesores)
    for(g = ndiv; g < ndiv * ndiv; g += ndiv)
    for(h = 0; h < ndiv; h++)
    for(i = 0; i < Point; i++)
    for(j = 0; j < NF; j++)
    if(j == i)
    grid[((NF * Point) * (g + h)) + (NF * i) + j] = grid[((NF * Point) * ((g - ndiv) + h)) + (NF * i) + j];
    else
    grid[((NF * Point) * (g + h)) + (NF * i) + j] = grid[((NF * Point) * ((g - ndiv) + h)) + (NF * i) + j] + amp[j];
  */
}

unsigned int compVec2(double *noDomF, double *fitness, unsigned int NF, unsigned int h, unsigned int i, unsigned int opt, unsigned int fun ,double *noDomP,double *pop,unsigned int D) {
  unsigned int sum = 0, ret, j;

  if(constraints(noDomP+(D*h),fun) < constraints(pop+(D*i),fun))return 1;
  if(constraints(noDomP+(D*h),fun) > constraints(pop+(D*i),fun))return 2;
  //  cout<<"      chepao"<<endl;

  //sum = sumcompVec(noDomF, NF, h, i, opt);
  for(j = 0; j < NF; j++)
    if (((noDomF[(NF * h) + j] < fitness[(NF * i) + j]) && (opt == 0)) || ((noDomF[(NF * h) + j] > fitness[(NF * i) + j]) && (opt == 1)))
      sum += 1;
	
  if (sum == NF)
    ret = 1;
  else
    if (sum == 0)
      ret = 2;
    else
      ret = 3;
			
  return ret;		
}


unsigned int verDom(double *noDomP, double *noDomF, double *fitness, unsigned int *hyperspace, unsigned int *partPos, unsigned int D, unsigned int NF, unsigned int opt, unsigned int *lastPart, unsigned int i,unsigned int fun,double *pop) {
  unsigned int h = 0, j;

  do {
    switch (compVec2(noDomF, fitness, NF, h, i, opt,fun,noDomP,pop,D)) {
      //Si fitness es dominada por un noDom.
    case 1: return 0;
                 		//Si fitness domina a un noDom.
    case 2:	//La partícula noDom dominada se elimina del repositorio
      for(j = 0; j < D; j++)
	noDomP[(D * h) + j] = noDomP[(D * (*lastPart - 1)) + j];
      for(j = 0; j < NF; j++)
	noDomF[(NF * h) + j] = noDomF[(NF * (*lastPart - 1)) + j];
      //La partícula se elimina del hyperespacio
      hyperspace[partPos[h]] -= 1;
      partPos[h] = partPos[*lastPart - 1];
      //Se disminute el número de partículas en el repositorio
      *lastPart -= 1;
      break;
    case 3:	h += 1;
      break;
    case 4:	//La partícula noDom dominada se elimina del repositorio
      for(j = 0; j < D; j++)
	noDomP[(D * h) + j] = noDomP[(D * (*lastPart - 1)) + j];
      for(j = 0; j < NF; j++)
	noDomF[(NF * h) + j] = noDomF[(NF * (*lastPart - 1)) + j];
      //La partícula se elimina del hyperespacio
      hyperspace[partPos[h]] -= 1;
      partPos[h] = partPos[*lastPart - 1];
      //Se disminute el número de partículas en el repositorio
      *lastPart -= 1;
      break;

    }
  } while(h < *lastPart);

  return 1;
}

void search_insert2(double *noDomP, double *noDomF, double *pop, double *fitness, unsigned int *linf, unsigned int *lsup, double *amp, double *start, unsigned int *hyperspace, unsigned int *partPos, unsigned int ndiv, unsigned int D, unsigned int NF, unsigned int M, unsigned int MEM, unsigned int opt, unsigned int maxCube, unsigned int *lastPart, unsigned int fun) {
  unsigned int i, j, k, l, aux_1, aux_2;

  for(k = 0; k < M; k++) {
    //Se verifica que la partícula candidata a ingresar al repositorio NO sea dominada por las demás
    //y de paso elimina a las que logre dominar (actualizando el hypercubo).
    //        if(constraints(noDomP+(NF*k),fun)>0)continue;
	

    if (verDom(noDomP, noDomF, fitness, hyperspace, partPos, D, NF, opt, lastPart, k,fun,pop) == 1) {
      //Si el repositorio aún no se ha llenado, ingresa la partícula y se actualiza el hypercubo
      if (*lastPart < MEM) {
	i = *lastPart;
	for(j = 0; j < D; j++)
	  noDomP[(D * i) + j] = pop[(D * k) + j];
	for(j = 0; j < NF; j++)
	  noDomF[(NF * i) + j] = fitness[(NF * k) + j];
	*lastPart += 1;
	partPos[i] = locHyper(noDomF, start, amp, ndiv, NF, i);
	if (partPos[i] != 0) {
	  partPos[i] -= 1;
	  hyperspace[partPos[i]] += 1;
	}
	else
	  makeHyper(noDomF, linf, lsup, amp, start, hyperspace, partPos, ndiv, MEM, NF, *lastPart, maxCube);
      } 
      //En caso de que la partícula no haya sido dominada, para ingresar en el repositorio debe de
      //encontrarse en un hypercubo no muy poblado.
      else {
	aux_1 = locHyper(fitness, start, amp, ndiv, NF, k);
	aux_2 = search_max(hyperspace, maxCube);
	//Si la partícula esta fuera del rango se recalcula el hyperespacio eliminando
	//una partícula del área más poblada
	if (aux_1 == 0) {
	  for(l = 0; l < *lastPart; l++)
	    if (partPos[l] == aux_2)
	      break;
	  for(j = 0; j < D; j++)
	    noDomP[(D * l) + j] = pop[(D * k) + j];
	  for(j = 0; j < NF; j++)
	    noDomF[(NF * l) + j] = fitness[(NF * k) + j];
	  makeHyper(noDomF, linf, lsup, amp, start, hyperspace, partPos, ndiv, MEM, NF, *lastPart, maxCube);
	}
				//Sino se busca eliminar a una partícula del área más poblada
	else {
	  	  aux_1 -= 1;
	  //	  if(aux_2<1){cout<<"cotzito"<<endl;exit(0);}

          			//Si el área de la partícula candidata es menos poblada que el área más poblada en el hyperespacio
          			//se inserta la partícula reemplazandola por una del área más poblada.
          			//Sino se descarta.
	  if (hyperspace[aux_1] < hyperspace[aux_2]) {
	    for(l = 0; l < *lastPart; l++)
	      if (partPos[l] == aux_2)
		break;
	    for(j = 0; j < D; j++)
	      noDomP[(D * l) + j] = pop[(D * k) + j];
	    for(j = 0; j < NF; j++)
	      noDomF[(NF * l) + j] = fitness[(NF * k) + j];
	    /*
	    if(partPos[l]-aux_1<0)cout<<"saludos"<<endl;
	    if(aux_2<0)cout<<"saludos"<<endl;
	    if(l<0||l>=MEM)cout<<"saludos"<<endl;
	    if(partPos[l]>=maxCube)cout<<"saludos"<<endl;
	    if(hyperspace[aux_2]<0)cout<<"saludos"<<endl;
	    cout<<partPos[l]<<" "<<aux_1<<" "<<maxCube<<endl;
	    */
	    partPos[l] = aux_1;
	    hyperspace[partPos[l]] += 1;
	    hyperspace[aux_2] -= 1;

	  /*GTP
	    // GTP*/
	  }
	}
      }
    }
  }
}

void hyperFit(unsigned int *hyperspace, unsigned int *hyperPoolL, double *hyperPoolF, double *hypsum, unsigned int maxCube, unsigned int *poolC) {
  unsigned int i;
  *poolC = 0;
  *hypsum = 0;
	
  for(i = 0; i < maxCube; i++)
    if (hyperspace[i] > 0) {
      hyperPoolF[*poolC] = 10.0 / hyperspace[i];
      hyperPoolL[*poolC] = i;
      *hypsum += hyperPoolF[*poolC];
      *poolC += 1;
    }
}

unsigned int hyperWheel(unsigned int *hyperPoolL, double *hyperPoolF, double hypsum, unsigned int poolC) {
  unsigned int i = 0;
  double rand, sum = 0;

  rand = RandomDouble(0.0, 1.0) * hypsum;
  do {
    sum = sum + hyperPoolF[i];
    i++;
  } while ((i < poolC) && (sum < rand));
  return hyperPoolL[i - 1];
}

unsigned int lead(unsigned int *hyperspace, unsigned int *hyperPoolL, unsigned int *partPos, double *hyperPoolF, double hypsum, unsigned int lastPart, unsigned int poolC) {
  unsigned int h, i, rand;

  h = hyperWheel(hyperPoolL, hyperPoolF, hypsum, poolC);
  rand = RandomInt(1, hyperspace[h]);
  //for(i = 0; i < 900;i++){
  //if(hyperspace[i]>0) cout<<hyperspace[h]<<endl;
  //}
  //  cout<<"lastPart "<<lastPart<<endl;
  //for(i=0;i<lastPart;i++)
  //cout<<partPos[i]<<endl;
  //exit(0);
  for(i = 0; i < lastPart; i++) {
    if (partPos[i] == h)
      rand -= 1;
    if (rand < 1)
      break;
  }
	
  return i;
}
unsigned int lead2(double *pop, unsigned int lastPart, unsigned int *selec, bool state,int indpob,double *noDomF, int NF) {
  unsigned int h, i, j, rand;
  /*
  h = hyperWheel(hyperPoolL, hyperPoolF, hypsum, poolC);
  rand = RandomInt(1, hyperspace[h]);
  for(i = 0; i < 900;i++){
    if(hyperspace[i]>0) cout<<hyperspace[h]<<endl;
  }
  cout<<"lastPart "<<lastPart<<endl;
  */
  //for(i=0;i<lastPart;i++)
  //cout<<partPos[i]<<endl;
  //exit(0);

  if(state==false)
    for(int i=0;i<80/10;i++)
      selec[i]=RandomInt(0,lastPart-1);

  
  double minimo=pow(10,10);
  int minimo1=-1;
  double disttemp=0;
  for( i = 0; i < 80/10; i++){
    for(j = 0;j<NF;j++)
      disttemp+=pow(pop[indpob*NF+j],2)+pow(noDomF[selec[i]*NF+j],2);

    if(disttemp < minimo){
      minimo=disttemp;
      minimo1=i;
    }
  }
  return selec[minimo1];
}
/*
  void clasifica(unsigned int *hyperspace, unsigned int *partPos, unsigned int *nivel, unsigned int maxCube, unsigned int lastPart) {
  unsigned int aux, max, min, x, y, i;
	
  max = search_max(hyperspace, maxCube);
  min = search_min(hyperspace, maxCube);
	
  aux = (unsigned int) ((max -min) / 3);
	
  x = min + aux;
  y = max - aux;
	
  for(i =0; i < lastPart; i++)
  if (hyperspace[partPos[i]] < x)
  nivel[i] = 1;
  else	if ((hyperspace[partPos[i]] > x) && (hyperspace[partPos[i]] < y))
  nivel[i] = 2;
  else
  nivel[i] = 3;
  }
*/
