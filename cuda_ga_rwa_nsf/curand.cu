#include <curand_kernel.h>
#include <cuda_runtime.h>
#include "util.cuh"

__device__ curandStateXORWOW_t state;

__device__ double dRand(int max)
{
    //curandState state;
    //curandStateXORWOW_t state;
    int tId = (blockDim.x * blockIdx.x) + threadIdx.x;
    curand_init((unsigned long long)clock() + tId, 0, 0, &state);

    return -(max) * curand_uniform(&state);
}

