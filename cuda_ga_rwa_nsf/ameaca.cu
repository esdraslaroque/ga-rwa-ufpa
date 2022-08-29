#include <helper_cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "info.h"
#include "nsf.h"
#include "rwa.h"
#include "util.cuh"

__device__ int count[SIM_MAX_LOAD-1];

__global__ void ameaca ()
{
    int load = threadIdx.x + 1;

    // double until_next   = -log(1-((int)dRand(BIGN+1))/(double)((unsigned)BIGN+1))/load;
	// double holding_time = -log(1-((int)dRand(BIGN+1))/(double)((unsigned)BIGN+1));
    unsigned long long int BIGNN = 2147483647;

    if (blockIdx.x == 1) {
        if (threadIdx.x == 1) { printf ("BIGN: %d\n", BIGNN+1); }
        printf ("%d\n", (int)dRand(BIGNN+1));

    }
}

int main (void)
{
    ameaca <<<SIM_NUM_GEN,(SIM_MAX_LOAD-1)>>>();
    // printf (" [curand: %f]\n", count_bp[0]);    
    cudaDeviceReset();
    return EXIT_SUCCESS;
}