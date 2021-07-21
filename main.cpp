#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h> 
#include <sys/shm.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <pthread.h>
#include <math.h>


double next_exp(int max, double lambda);

int main(int argc, char ** argv){
	if(argc != 8){
		fprintf(stderr, "Not enough command-line arguments.\n");
		return EXIT_FAILURE;
	}
	int numProc = atoi(argv[1]);
	if(numProc < 1 || numProc > 26){
		fprintf(stderr, "Invalid number of processes.\n");
		return EXIT_FAILURE;
	}
	int seed = atoi(argv[2]);
	if(seed == 0){
		fprintf(stderr, "Invalid seed.\n");
		return EXIT_FAILURE;
	}
	double lambda = atof(argv[3]);
	if(lambda == 0.0f){
		fprintf(stderr, "Invalid lambda.\n");
		return EXIT_FAILURE;
	}
	int upBound = atoi(argv[4]);
	if(upBound == 0){
		fprintf(stderr, "Invalid upper bound.\n");
		return EXIT_FAILURE;
	}
	int csTime = atoi(argv[5]);
	if(csTime == 0){
		fprintf(stderr, "Invalid context switch time.\n");
		return EXIT_FAILURE;
	}
	double alpha = atof(argv[6]);
	if(alpha == 0.0f){
		fprintf(stderr, "Invalid alpha.\n");
		return EXIT_FAILURE;
	}
	int timeSlice = atoi(argv[7]);
	if(timeSlice == 0){
		fprintf(stderr, "Invalid time slice.\n");
		return EXIT_FAILURE;
	}
	srand48(seed);
	next_exp(upBound, lambda);
}

double next_exp(int max, double lambda){
  double x = 0; 

  while(true)
  {
    double r = drand48();   /* uniform dist [0.00,1.00) -- also see random() */

    x = -log( r ) / lambda;  /* log() is natural log */

    /* skip values that are far down the "long tail" of the distribution */
    if ( x > max ) { continue; }

    break;
  }	
  printf("%f\n", x);
  return x;
}