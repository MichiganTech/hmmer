//Header file for sre_random.c

double
sre_random_positive();


double
ExponentialRandom();


double
Gaussrandom(
	double mean,
	double stddev);


int   
DChoose(
	double *p,
	int N);


int   
FChoose(
	float *p,
	int N);


#define CHOOSE(a)   ((int) (sre_random() * (a)))
