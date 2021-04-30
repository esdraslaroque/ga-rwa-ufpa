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

__device__ int count_bp[SIM_MAX_LOAD-1];

__global__ void rwa (struct Vertice *d_verts, struct Edge *d_edges)
{
	int myId = threadIdx.x + blockDim.x * blockIdx.x;
	int stride = blockDim.x * gridDim.x;
	int offsets = threadIdx.x * NSF_NUM_EDGES;
	int offsetf = (threadIdx.x * NSF_NUM_EDGES)+NSF_NUM_EDGES;

	if (myId > stride-1) return;

	int fitness[GA_MIN_GEN];
	int load = threadIdx.x + 1;

	// Local structures for GA and NSF
	struct Vertice verts[NSF_NUM_NODES];
	struct Edge edges[NSF_NUM_EDGES];

	// Copying vertices from global memory
	for (int i = 0; i < NSF_NUM_NODES; i++)
		verts[i] = d_verts[i];

	// Copying edges from global memory
	int pv = 0;
	for (int i = offsets; i < offsetf; i++) {
		edges[pv] = d_edges[i];
		pv += 1;
	}
#ifdef DEBUG_NSF
	d_printNsf (verts, edges, 0); 
#endif
        // Initializing individuals
	struct Chromosome individual[GA_SIZE_POP];

        for (int i = 0; i < GA_SIZE_POP; i++)
                individual[i].length = 0;

	double until_next   = -log(1-((int)dRand(BIGN+1))/(double)((unsigned)BIGN+1))/load;
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

	//printPopulation (individual, PREVIEW_AVAIL, myId);
	int color;
	
	if (individual[0].wl_avail > 0) {
	        for (color = 0; color < NSF_NUM_CHANNELS; color++)
	                if (((individual[0].L & (1<<color)) >> color)) break;
	
	        for (int i = 0; i < individual[0].length-1; i++) {
	                for (int j = 0; j < verts[i].adjLen; j++) {
	                        int e = verts[i].edges[j];
	
	                        if (edges[e].u != individual[0].gene[i+1]) continue;
	
	                        edges[e].wave -= (int) pow (2.0, color);
	                        edges[e].time[color] = holding_time;
	
	                        for (int z = 0; z < verts[i+1].adjLen; z++) {
	                                e = verts[i+1].edges[z];
	                                if (verts[e].label != individual[0].gene[i]) continue;
	                                edges[e].wave -= (int) pow (2.0, color);
	                                edges[e].time[color] = holding_time;
	                        }
	
	                        break;
	                }
	        }
	} else atomicAdd(&count_bp[threadIdx.x], 1);

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

	//d_printNsf (verts, edges, 0);
	pv = 0;
	for (int i = offsets; i < offsetf; i++) {
		d_edges[i] = edges[pv];
		pv += 1;
	}
        __syncthreads();
	
	//d_printNsf (verts, edges, 0); 
	//if (myId == 0) printf ("until_next: %2f\tholding_time: %2f\n", until_next, holding_time);

        if (threadIdx.x == 0 && blockIdx.x == 127) {
		float bp = 0.0;
		for (int i = 0; i < SIM_MAX_LOAD-1; i++) {
			bp = 100.0*(count_bp[i]/(float)SIM_NUM_GEN);
              		//printf ("Allocated: %d\nBlocked: %d\nBlock Probability: %.1f%%\n", counts[0], counts[1], bp);
			//printf ("Block Probability [Erlang %2d]: %.1f%%\n", i, bp);
              		//printf ("Block Probability: %.1f\n", bp);
              		printf ("%.1f ", bp);
		}
        }

	if (myId == stride-1) printf ("\nblockDim.x: %d\ngridDim.x: %d\nstride: %d\n\n", blockDim.x, gridDim.x, stride);
}

/***********
 * Main ()
 ***********/
int main (void)
{
	clock_t t;

        // Creating NSF
        struct Vertice verts[NSF_NUM_NODES];
        struct Edge edges[NSF_NUM_EDGES*(SIM_MAX_LOAD-1)];

        for (int i = 0; i < NSF_NUM_NODES; i++) {
                verts[i].label = i;
                verts[i].adjLen = 0;
                verts[i].isVisited = FALSE;
        }

	// Device pointers
	struct Vertice *d_verts;
	struct Edge *d_edges;

	checkCudaErrors( cudaMalloc ((void**)&d_verts, NSF_NUM_NODES * sizeof(struct Vertice)) );
	checkCudaErrors( cudaMalloc ((void**)&d_edges, NSF_NUM_EDGES*(SIM_MAX_LOAD-1) * sizeof(struct Edge)) );

	srand ((unsigned) time(NULL));

	// 20 uplinks in GPU memory
	int E = 0;
        addLink (verts, edges, 0, 1, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 0, 11, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 0, 12, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 0, 13, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 0, 16, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 0, 25, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 1, 2, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 1, 3, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 3, 4, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 3, 6, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 4, 5, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 6, 7, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 7, 8, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 7, 10, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 8, 9, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 9, 10, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 10, 11, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 10, 12, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 10, 25, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 11, 23, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 12, 22, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 12, 25, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 13, 14, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 13, 16, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 15, 16, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 16, 17, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 16, 24, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 17, 18, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 18, 19, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 18, 20, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 18, 22, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 18, 25, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 19, 27, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 20, 21, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 21, 22, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 22, 23, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 24, 25, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 25, 26, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);
    addLink (verts, edges, 26, 27, (rand() % (1<<NSF_NUM_CHANNELS)), NSF_NUM_CHANNELS, &E, 2);

        /* init erlangs randomnly */
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
                                                edges[link].time[c] = -log(1-(rand()%RAND_MAX+1)/(double)((unsigned)RAND_MAX+1));
                                        }
                                }
                        }
                }
        }

	int dataset = NSF_NUM_EDGES, pv = 0; //, k = 0;
	for (int i = NSF_NUM_EDGES; i < (NSF_NUM_EDGES*(SIM_MAX_LOAD-1)); i++) {
		edges[i] = edges[i-dataset];

		pv += 1;
		//printf ("Copying to edges[%d] from edges[%d] - pv: %d\n", i, i-dataset, pv);

		if (pv == NSF_NUM_EDGES) {
			dataset = dataset+NSF_NUM_EDGES;
			pv = 0;
			//printf ("**[%d] Now, dataset is %d and pv: %d\n", k, dataset, pv);
			//k++;
		}
	}
	//printf ("Terminou de gerar edges..\n");

	//for (int i = 0; i < (NSF_NUM_EDGES*(SIM_MAX_LOAD-1)); i++) {
	//	printf ("Edges[%4d] => v: %2d\tu: %2d\twave: %2d\ttimes: (", i, edges[i].v, edges[i].u, edges[i].wave);
	//	for (int c = 0; c < NSF_NUM_CHANNELS; c++)
	//		printf (" %.2f", edges[i].time[c]);
	//	printf (")\n");
	//}

	checkCudaErrors( cudaMemcpy (d_verts, verts, NSF_NUM_NODES * sizeof(struct Vertice), cudaMemcpyHostToDevice) );
	checkCudaErrors( cudaMemcpy (d_edges, edges, (NSF_NUM_EDGES*(SIM_MAX_LOAD-1)) * sizeof(struct Edge), cudaMemcpyHostToDevice) );

	t = clock();
	//v1: rwa <<<(SIM_MAX_LOAD-1),SIM_NUM_GEN>>>(d_verts, d_edges, d_erlangs);
	rwa <<<SIM_NUM_GEN,(SIM_MAX_LOAD-1)>>>(d_verts, d_edges);
	//rwa <<<1,1>>>(d_verts, d_edges, d_erlangs);
	cudaDeviceSynchronize();
	t = clock() - t;

	double time_taken = ((double) t)/CLOCKS_PER_SEC;
	//float bp = h_counts[1]/(float)SIM_NUM_GEN;

	//printf ("Allocated: %d\nBlocked: %d\nBlock Probability: %.1f%%\n\nTime taken: %f\n", h_counts[0], h_counts[1], bp, time_taken);
	printf ("Time taken: %.2f\n", time_taken);

	cudaFree(d_verts);
	cudaFree(d_edges);
	cudaDeviceReset();
	return EXIT_SUCCESS;
}

