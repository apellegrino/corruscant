#include <math.h>
#include <stdio.h>

/* for call by scipy.integrate */
double f(int n, double args[n]) {
	return 1.0/sqrt(0.3*args[0]+0.7*args[0]*args[0]*args[0]*args[0]);
}
