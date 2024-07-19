/***************************************************************************
                          main.cpp  -  description
                             -------------------
    begin                : Sun Jan 28 14:04:34 CST 2001
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

//#ifdef HAVE_CONFIG_H
//#include <config.h>
//#endif

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
//unsigned int CONT_FUN = 0;

#include <string.h>
#include <math.h>
#include <time.h>
double pM=0.05;//10300 0.1 10400 10500 0.1
#include "randomlib.h"
#include "fun-res.h"
#include "fun-eng.h"
#include "fun-SR.h"
#include "fun-moo.h"
//#include <estructuras.h>
#include "variables.h"
#include "psolib.h"
#include "mainlib.h"

//namespace;
using namespace std;

int main(int argc, char *argv[])
{
	unsigned int funcion, particulas, ciclos, optimizacion, MEM, ndiv, i;
	char arch1[20];
	clock_t  now, later;
	double   passed=0.0;
	//FILE *time;

	//time = fopen("time3.dat","w");

	/* Funcion a optimizar
	   Funciones mono-objetivo sin restricciones:
	   9,10,11
	   Funciones mono-objetivo con restricciones:
	   1,2,3,4,5,6
	   Funciones de ingenieria:
	   7,8
	   Funciones multiobjetivo sin restricciones:
	   100,200,300,400,500,*600,700,*800,900,1000

	   10100 50 0.05//kita
	   10300 80 0.05//truss *
	   10400 100 //deb2
	   10500 40 0.05//deb
	   200 120 0.05//kursawe
	*/
	funcion = 10100;
	// Numero de particulas
	particulas = 100;
	// Numero de ciclos
	ciclos = 80;
	// En caso de minimizar = 0, en caso de maximizar = 1
	optimizacion = 0;
	//Tamaño de la memoria
	MEM = 100;
	//Divisiones del espacio
	ndiv = 30;

	sprintf(arch1,"Pareto.dat");

	now = clock();
	//PSO
	vuelo(funcion, particulas, ciclos, optimizacion, num_dim(funcion), num_fun(funcion), MEM, ndiv, arch1);
	//	cout<<"fin"<<endl;

  return EXIT_SUCCESS;
}
