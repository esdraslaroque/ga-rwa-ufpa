#include <curand_kernel.h>
#include <cuda_runtime.h>
#include "util.cuh"

__device__ curandStateXORWOW_t state;
__device__ unsigned long long int BIGNN = 2147483647;

__device__ double dRand()
{
    //curandState state;
    //curandStateXORWOW_t state;
    int tId = (blockDim.x * blockIdx.x) + threadIdx.x;
    curand_init((unsigned long long)clock() + tId, 0, 0, &state);

    return (BIGNN+1) * curand_uniform(&state);
}

