/*
 * Programa sequencial para simular AG com RWA
 * --------------------------------------------------
 * Autor: Esdras La-Roque <esdras.laroque@gmail.com>
 * Criado: 25/05/2018
 *
 * Adapted from:
 * ----------------
 * GA: RWA with GOF
 * Genetic Algorithm: 
 * Routing and Wavelength Assignment with General Objective Function
 *
 * Copyright 2017 Universidade Federal do Pará (PPGCC UFPA)
 *
 * Authors: April 2017
 * Cassio Trindade Batista - cassio.batista.13@gmail.com
*/
#include <helper_cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "info.h"
#include "nsf.h"
#include "rwa.h"
#include "util.cuh"

__global__ void rwa (struct Vertice *d_verts, struct Edge *d_edges, double *erlangs)
{
	int myId = threadIdx.x + blockDim.x * blockIdx.x;
	int stride = blockDim.x * gridDim.x;
	__shared__ int count_bp;

	if (myId > stride) return;

	int fitness[GA_MIN_GEN];
	int load = blockIdx.x + 1;

	// Local structures for GA and NSF
	struct Vertice verts[NSF_NUM_NODES];
	struct Edge edges[NSF_NUM_EDGES];

	// Copying vertices from global memory
	for (int i = 0; i < NSF_NUM_NODES; i++)
		verts[i] = d_verts[i];
	// Copying edges from global memory
	for (int i = 0; i < NSF_NUM_EDGES; i++)
		edges[i] = d_edges[i];

	for (int i = 0; i < NSF_NUM_NODES; i++) {
	        for (int e = 0; e < verts[i].adjLen; e++) {
	                int link = verts[i].edges[e];

			if (edges[link].u < i) {
				// Making symetric
				for (int ee = 0; ee < verts[edges[link].u].adjLen; ee++) {
					int link_cp = verts[edges[link].u].edges[ee];

					if (edges[link_cp].u != i) continue;

					for (int c = 0; c < NSF_NUM_CHANNELS; c++)
						edges[link].time[c] = edges[link_cp].time[c];
				}
			} else {
	                	for (int c = 0; c < NSF_NUM_CHANNELS; c++) {
	                	        if ( !((edges[link].wave & (1<<c)) >> c) ) { 
	                	                //edges[link].time[c] = erlangs[blockIdx.x*NSF_NUM_CHANNELS+c];
	                	                edges[link].time[c] = -log(1-(dRand(BIGN+1))/(double)((unsigned)BIGN+1));
					}
	                	}
			}
	        }
	}
#ifdef DEBUG_NSF
	d_printNsf (verts, edges, myId); 
#endif
        // Initializing individuals
	struct Chromosome individual[GA_SIZE_POP];

        for (int i = 0; i < GA_SIZE_POP; i++)
                individual[i].length = 0;

	double until_next = -log(1-((int)dRand(BIGN+1))/(double)((unsigned)BIGN+1))/load;
	double holding_time = -log(1-((int)dRand(BIGN+1))/(double)((unsigned)BIGN+1));

        /* Generating individuals */
        for (int i = 0; i < GA_SIZE_POP; i++) 
                addIndividual (verts, edges, &individual[i], i, NSF_SRC_NODE, NSF_DST_NODE);
#ifdef DEBUG_POP
	printPopulation (individual, PREVIEW_AVAIL, myId);
#endif
	// ******** Processing GA ******** //
	for (int gen = 0; gen < GA_MIN_GEN; gen++) {
		evaluate (individual, verts, edges);

		int *mating_pool = selection (individual, GA_MAX_CROSS_RATE);

                if (mating_pool) 
                        crossover (individual, mating_pool);

		struct Chromosome *chromosome = NULL;
		for (int i = 0; i < GA_SIZE_POP; i++) {
			if ((int)dRand(100) > GA_MIN_MUT_RATE) continue;

			chromosome = mutate (verts, edges, &individual[i]);

                 	if (chromosome) {
                 	        individual[i].idx = i;
                 	        individual[i].length = chromosome->length;
                 	        individual[i].wl_avail = chromosome->wl_avail;
                 	        individual[i].L = chromosome->L;
                 	        individual[i].S = chromosome->S;
                 	        for (int j = 0; j < NSF_NUM_NODES; j++) {
                 	                individual[i].gene[j] = chromosome->gene[j];
                 	        }
                 	}
		} // Mutate

		insertionSort (individual, BY_LENGTH);
		insertionSort (individual, BY_WL_AVAIL);

		//if (gen == 20) printPopulation (individual, PREVIEW_AVAIL, myId);
		
		fitness[gen] = 0;

                for (int i = 0; i < GA_SIZE_POP; i++) {
                        if (individual[i].wl_avail > 0)
                                fitness[gen] += 1;
		}
	} // GA_MIN_GEN

        //if (myId == 0) {
        //        for (int i = 0; i < GA_MIN_GEN; i++)
        //                printf ("%d\n", fitness[i]);
        //}

	evaluate (individual, verts, edges);
	insertionSort (individual, BY_LENGTH);
	insertionSort (individual, BY_WL_AVAIL);

	int color;
	
	if (individual[0].wl_avail > 0) {
	        for (color = 0; color < NSF_NUM_CHANNELS; color++)
	                if (((individual[0].L & (1<<color)) >> color)) break;
	
	        for (int i = 0; i < individual[0].length-1; i++) {
	                for (int j = 0; j < verts[i].adjLen; j++) {
	                        int e = verts[i].edges[j];
	
	                        if (edges[e].u != individual[0].gene[i+1]) continue;
	
	                        edges[e].wave -= pow (2.0, color);
	                        edges[e].time[color] = holding_time;
	
	                        for (int z = 0; z < verts[i+1].adjLen; z++) {
	                                e = verts[i+1].edges[z];
	                                if (verts[e].label != individual[0].gene[i]) continue;
	                                edges[e].wave -= pow (2.0, color);
	                                edges[e].time[color] = holding_time;
	                        }
	
	                        break;
	                }
	        }
	} else atomicAdd(&count_bp, 1);

	// Updating traffic matrix
	for (int i = 0; i < NSF_NUM_NODES; i++) {
	        for (int e = 0; e < verts[i].adjLen; e++) {
	                int link = verts[i].edges[e];
	
	                for (int c = 0; c < NSF_NUM_CHANNELS; c++) {
	                        if (edges[link].time[c] > until_next) {
	                                edges[link].time[c] -= until_next;
	                        } else {
	                                edges[link].time[c] = 0.0;
	                                if ( !((edges[link].wave & (1<<c)) >> c) ) {
	                                        edges[link].wave += (int) pow (2.0, c);
	                                }
	                        }
	                }
	        }
	}

	d_printNsf (verts, edges, myId);
        //__syncthreads();
        //if (threadIdx.x == 0) printf ("Done\n");
        //if (threadIdx.x == 127) {
        //      float bp = 100.0*(count_bp/(float)SIM_NUM_GEN);
        //      //printf ("Allocated: %d\nBlocked: %d\nBlock Probability: %.1f%%\n", counts[0], counts[1], bp);
        //      printf ("Block Probability: %.1f\n", bp);
        //}
}

/***********
 * Main ()
 ***********/
int main (void)
{
	clock_t t;

        // Creating NSF
        struct Vertice verts[NSF_NUM_NODES];
        struct Edge edges[NSF_NUM_EDGES];
	int initErlangs = (NSF_NUM_NODES*NSF_NUM_CHANNELS)*(SIM_MAX_LOAD-1);

        for (int i = 0; i < NSF_NUM_NODES; i++) {
                verts[i].label = i;
                verts[i].adjLen = 0;
                verts[i].isVisited = FALSE;
        }

	// Device pointers
	struct Vertice *d_verts;
	struct Edge *d_edges;
	double *d_erlangs;

	checkCudaErrors( cudaMalloc ((void**)&d_verts, NSF_NUM_NODES * sizeof(struct Vertice)) );
	checkCudaErrors( cudaMalloc ((void**)&d_edges, NSF_NUM_EDGES * sizeof(struct Edge)) );
	checkCudaErrors( cudaMalloc ((void**)&d_erlangs, initErlangs * sizeof(double)) );

	srand ((unsigned) time(NULL));

	// 20 uplinks in GPU memory
	int E = 0;
        addLink (verts, edges,  0,  1, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
        addLink (verts, edges,  0,  2, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
        addLink (verts, edges,  0,  5, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
        addLink (verts, edges,  1,  2, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
        addLink (verts, edges,  1,  3, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
        addLink (verts, edges,  2,  8, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
        addLink (verts, edges,  3,  4, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
        addLink (verts, edges,  3,  6, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
        addLink (verts, edges,  3, 13, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
        addLink (verts, edges,  4,  9, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
        addLink (verts, edges,  5,  6, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
        addLink (verts, edges,  5, 10, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
        addLink (verts, edges,  6,  7, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
        addLink (verts, edges,  7,  8, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
        addLink (verts, edges,  8,  9, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
        addLink (verts, edges,  9, 11, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
        addLink (verts, edges,  9, 12, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
        addLink (verts, edges, 10, 11, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
        addLink (verts, edges, 10, 12, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
        addLink (verts, edges, 11, 13, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);

	double erlangs[initErlangs];

	for (int i = 0; i < initErlangs; i++)
		erlangs[i] = -log(1-(rand()%RAND_MAX+1)/(double)((unsigned)RAND_MAX+1));

	checkCudaErrors( cudaMemcpy (d_verts, verts, NSF_NUM_NODES * sizeof(struct Vertice), cudaMemcpyHostToDevice) );
	checkCudaErrors( cudaMemcpy (d_edges, edges, NSF_NUM_EDGES * sizeof(struct Edge), cudaMemcpyHostToDevice) );
	checkCudaErrors( cudaMemcpy (d_erlangs, erlangs, initErlangs * sizeof(double), cudaMemcpyHostToDevice) );

	t = clock();
	//rwa <<<(SIM_MAX_LOAD-1),SIM_NUM_GEN>>>(d_nsf, d_erlangs);
	rwa <<<(SIM_MAX_LOAD-1),SIM_NUM_GEN>>>(d_verts, d_edges, d_erlangs);
	//rwa <<<1,1>>>(d_verts, d_edges, d_erlangs);
	cudaDeviceSynchronize();


	t = clock() - t;
	double time_taken = ((double) t)/CLOCKS_PER_SEC;
	//float bp = h_counts[1]/(float)SIM_NUM_GEN;

	//printf ("Allocated: %d\nBlocked: %d\nBlock Probability: %.1f%%\n\nTime taken: %f\n", h_counts[0], h_counts[1], bp, time_taken);
	printf ("Time taken: %.2f\n", time_taken);

	cudaFree(d_verts);
	cudaFree(d_edges);
	cudaFree(d_erlangs);
	cudaDeviceReset();
	return EXIT_SUCCESS;
}

