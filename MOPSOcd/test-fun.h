/****** Objective functions for test problems *******/

double kita_f1(unsigned int i)
{
	return ( -(popVar[i][0] * popVar[i][0]) + popVar[i][1] );
}

double kita_f2(unsigned int i)
{
	return ( (popVar[i][0] / 2.0) + popVar[i][1] + 1 );
}

double kursawe_f1(unsigned int i)
{
	double r = 0.0;
	unsigned int j;

	for(j = 0; j < 2; j++)
		r += -10.0 * exp(-0.2 * sqrt(pow(popVar[i][j], 2) + pow(popVar[i][j + 1], 2)) );

	return r;
}

double kursawe_f2(unsigned int i)
{
	double r = 0.0;
	unsigned int j;

	for(j = 0; j < 3; j++)
		r += pow(fabs(popVar[i][j]), 0.8) + 5.0 * sin(pow(popVar[i][j], 3));
	return r;
}

double deb_f1(unsigned int i)
{
  return (popVar[i][0]);
}

double deb_f2(unsigned int i)
{
  double g = 2.0 - exp(-pow(((popVar[i][1]-0.2)/0.004),2)) - 0.8 * exp(-pow(((popVar[i][1]-0.6)/0.4),2));

  return ((double)g / popVar[i][0]);
}


double DTLZ6_f1(unsigned int i)
{
  return (popVar[i][0]);
}

double DTLZ6_f2(unsigned int i)
{
  return (popVar[i][1]);
}

double DTLZ6_f3(unsigned int i)
{
	unsigned int j = 0;
    unsigned int n = maxvar;
    unsigned int k = n - maxfun + 1;
	double g, h, s, t;

	s = 0;
	for (j = 2; j < maxvar; j++){
		s += popVar[i][j];
    }

	g = 1 + (9 / 20)  * s;

	t = 0;
	for (j = 0; j < maxfun-1; j++){
		t += (popVar[i][j]  * (1 + sin(3 * PI * popVar[i][j]))) / (1 + g);
    }

	h = 3 - t;

  return ( (1 + g) * h);
}

void deb2(unsigned int i) {
	double g;
	double h;
	double f1;

	f1 = popVar[i][0];
	g = 1.0 + 10.0 * popVar[i][1];
	h = 1.0 - pow((f1 / g), 2.0) - (f1 / g) * sin(12.0 * PI * f1);

	popFit[i][0] = f1;
	popFit[i][1] = g * h ;
	return;
}

void deb3(unsigned int i) {
	double g;
	double h;
	double f1;
	double alfa;
	double beta;
	double seno;
	//double pi = 3.1415926535;
	double argumento;
	double q;

	alfa = 10.0;
	q = 10.0;
	beta = 1.0;

	argumento = q * PI * popVar[i][0];
	seno = sin(argumento);
	f1 = 1.0 - exp(-4.0 * popVar[i][0]) * pow(seno, 4);
	g = 1.0 + popVar[i][1] * popVar[i][1];
	if (f1 <= beta * g)
	{
		h = 1.0 - pow((f1 / (beta * g)), alfa);
	}
	else
	{
		h = 0.0;
	}
	popFit[i][0] = f1;
	popFit[i][1] = g * h;
	return;
}

void fonseca2(unsigned int i)
{
	double s1, s2;
	int j;
	s1 = s2 = 0.0;
	for (j = 0; j < maxvar; j++)
	{
		s1 += pow((popVar[i][j] - (1.0 / sqrt((double) maxvar))), 2.0);
		s2 += pow((popVar[i][j] + (1.0 / sqrt((double) maxvar))), 2.0);
	}

	popFit[i][0] = 1.0 - exp(-s1);
	popFit[i][1] = 1.0 - exp(-s2);
	return;
}


/*ZDT*/

void ZDT1(unsigned int i)
{
	double f1, f2, g, h, sum;
	int j;

	f1 = popVar[i][0];
	sum = 0.0;
	for (j = 1; j < maxvar; j++)
	{
		sum += popVar[i][j];
	}
	g = 1.0 + 9.0 * sum / ((double)maxvar - 1.0);
	h = 1.0 - sqrt(f1 / g);
	f2 = g * h;
	popFit[i][0] = f1;
	popFit[i][1] = f2;
	return;
}

/* ZDT2 - MOP */
void ZDT2(unsigned int i)
{
	double f1, f2, g, h, sum;
	int j;

	f1 = popVar[i][0];
	sum = 0.0;
	for (j = 1; j < maxvar; j++)
	{
		sum += popVar[i][j];
	}
	g = 1.0 + 9.0 * sum / ((double)maxvar - 1.0);
	h = 1.0 - pow((f1 / g), 2.0);
	f2 = g * h;
	popFit[i][0] = f1;
	popFit[i][1] = f2;
	return;
}

/* ZDT3 - MOP */
void ZDT3(unsigned int i)
{
	double f1, f2, g, h, sum;
	int j;

	f1 = popVar[i][0];
	sum = 0.0;
	for (j = 1; j < maxvar; j++)
	{
		sum += popVar[i][j];
	}
	g = 1.0 + 9.0 * sum / ((double)maxvar - 1.0);
	h = 1.0 - sqrt(f1 / g) - (f1 / g) * sin(10.0 * PI * f1);
	f2 = g * h;
	popFit[i][0] = f1;
	popFit[i][1] = f2;
	return;
}

/* ZDT4 - MOP */
void ZDT4(unsigned int i)
{
	double f1, f2, g, h, sum;
	int j;

	f1 = popVar[i][0];
	sum = 0.0;
	for (j = 1; j < maxvar; j++)
	{
		sum += (pow(popVar[i][j], 2.0) - 10.0 * cos(4.0 * PI * popVar[i][j]));
	}
	g = 1.0 + 10.0 * ((double)maxvar - 1.0) + sum;
	h = 1.0 - sqrt(f1 / g);
	f2 = g * h;
	popFit[i][0] = f1;
	popFit[i][1] = f2;
	return;
}

/* ZDT6 - MOP */
void ZDT6(unsigned int i)
{
	double f1, f2, g, h, sum;
	int j;

	f1 = 1.0 - exp(-4.0 * popVar[i][0]) * pow(sin(6.0 * PI * popVar[i][0]), 6.0);
	sum = 0.0;
	for (j = 1; j < maxvar; j++)
	{
		sum += popVar[i][j];
	}
	g = 1.0 + 9.0 * pow((sum / ((double)maxvar - 1.0)), 0.25);
	h = 1.0 - pow((f1 / g), 2.0);
	f2 = g * h;
	popFit[i][0] = f1;
	popFit[i][1] = f2;
	popFit[i][1] = f2;
	return;
}


/*DTLZ*/

void DTLZ1(unsigned int i)
{
	int j;
	double sum;
	double g;
	int n = maxvar;
	int m = maxfun;
	int k = n - m + 1;

	sum = 0.0;
	for (j = m - 1; j < n; j++)
	{
		sum += pow(popVar[i][j] - 0.5, 2.0) - cos(20.0 * PI * (popVar[i][j] - 0.5));
	}
	g = 100.0 * (k + sum);

	popFit[i][0]= (0.5 * (1.0 + g) * popVar[i][0] * popVar[i][1]);
	popFit[i][1] = 0.5 * (1.0 + g) * (1.0 - popVar[i][1]) * popVar[i][0];
	popFit[i][2] = 0.5 * (1.0 + g) * (1.0 - popVar[i][0]);
	return;
}


/* DTLZ2 - MOP */
void DTLZ2(unsigned int i)
{
	int j;
	double sum;
	double g;
	int n = maxvar;
	int m = maxfun;

	sum = 0.0;
	for (j = m - 1; j < n; j++)
	{
		sum += pow(popVar[i][j] - 0.5, 2.0);
	}
	g = sum;

	popFit[i][0] = (1.0 + g) * cos(popVar[i][0] * PI * 0.5) * cos(popVar[i][1] * PI * 0.5);
	popFit[i][1] = (1.0 + g) * cos(popVar[i][0] * PI * 0.5) * sin(popVar[i][1] * PI * 0.5);
	popFit[i][2] = (1.0 + g) * sin(popVar[i][0] * PI * 0.5);
	return;
}

/* DTLZ3 - MOP */
void DTLZ3(unsigned int i)
{
	int j;
	double sum;
	double g;
	int n = maxvar;
	int m = maxfun;
	int k = n - m + 1;

	sum = 0.0;
	for (j = m - 1; j < n; j++)
	{
		sum += pow(popVar[i][j] - 0.5, 2.0) - cos(20.0 * PI * (popVar[i][j] - 0.5));
	}
	g = 100.0 * (k + sum);

	popFit[i][0] = (1.0 + g) * cos(popVar[i][0] * PI * 0.5) * cos(popVar[i][1] * PI * 0.5);
	popFit[i][1] = (1.0 + g) * cos(popVar[i][0] * PI * 0.5) * sin(popVar[i][1] * PI * 0.5);
	popFit[i][2] = (1.0 + g) * sin(popVar[i][0] * PI * 0.5);
	return;
}

/* DTLZ4 - MOP */
void DTLZ4(unsigned int i)
{
	int j;
	double alpha;
	double g, factor1, factor2;
	int n = maxvar;
	int m = maxfun;

	alpha = 100;
	g = 0.0;
	for (j = m - 1; j < n; j++)
	{
		g += pow(popVar[i][j] - 0.5, 2.0);
	}

	factor1 = pow(popVar[i][0], alpha);
	factor2 = pow(popVar[i][1], alpha);
	popFit[i][0] = (1.0 + g) * cos(factor1 * PI * 0.5) * cos(factor2 * PI * 0.5);
	popFit[i][1] = (1.0 + g) * cos(factor1 * PI * 0.5) * sin(factor2 * PI * 0.5);
	popFit[i][2] = (1.0 + g) * sin(factor1 * PI * 0.5);
	return;
}

/* DTLZ5 - MOP */
void DTLZ5(unsigned int i)
{
	int j;
	double g;
	double tetha[3];
	int n = maxvar;
	int m = maxfun;

	g = 0.0;
	for (j = m - 1; j < n; j++)
	{
		g += pow(popVar[i][j] - 0.5, 2.0);
	}

	tetha[0] = (PI * 0.5) * popVar[i][0];
	tetha[1] = (PI / (4.0 * (1.0 + g))) * (1.0 + 2.0 * g * popVar[i][1]);

	popFit[i][0] = (1.0 + g) * cos(tetha[0]) * cos(tetha[1]);
	popFit[i][1] = (1.0 + g) * cos(tetha[0]) * sin(tetha[1]);
	popFit[i][2] = (1.0 + g) * sin(tetha[0]);
	return;
}

/* DTLZ6 - MOP */
void DTLZ6(unsigned int i)
{
	int j;
	double g;
	double tetha[3];
	int n = maxvar;
	int m = maxfun;

	g = 0.0;
	for (j = m - 1; j < n; j++)
	{
		g += pow(popVar[i][j], 0.1);
	}

	tetha[0] = (PI * 0.5) * popVar[i][0];
	tetha[1] = (PI / (4.0 * (1.0 + g))) * (1.0 + 2.0 * g * popVar[i][1]);
	popFit[i][0] = (1.0 + g) * cos(tetha[0]) * cos(tetha[1]);
	popFit[i][1] = (1.0 + g) * cos(tetha[0]) * sin(tetha[1]);
	popFit[i][2] = (1.0 + g) * sin(tetha[0]);
	return;
}

/* DTLZ7 - MOP */
void DTLZ7(unsigned int i)
{
	int j;
	double sum;
	double g, h;
	int n = maxvar;
	int m = maxfun;
	int k = n - m + 1;

	sum = 0.0;
	for (j = m - 1; j < n; j++)
	{
		sum += popVar[i][j];
	}
	g = 1.0 + (9.0 / (double) k) * sum;

	popFit[i][0] = popVar[i][0];
	popFit[i][1] = popVar[i][1];
	h = m - (((popFit[i][0] / (1.0 + g)) * (1.0 + sin(3.0 * PI * popFit[i][0]))) + ((popFit[i][1]/(1.0 + g)) * (1.0 + sin(3.0 * PI * popFit[i][1]))));
	popFit[i][2] = (1.0 + g) * h;
	return;
}
