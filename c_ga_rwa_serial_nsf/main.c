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
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "info.h"
#include "nsf.h"
#include "rwa.h"

/***********
 * Main ()
 ***********/
int main (void)
{
	clock_t t;
	time_t seed;

	// Creating NSF
	struct Vertice verts[NSF_NUM_NODES];
	struct Edge edges[NSF_NUM_EDGES];

	for (int i = 0; i < NSF_NUM_NODES; i++) {
		verts[i].label = i;
		verts[i].adjLen = 0;
		verts[i].isVisited = FALSE;
	}

	srand ((unsigned) time(&seed));

	// 20 uplinks in CPU memory
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

	t = clock();
	int count_bp;

	for (int load = SIM_MIN_LOAD; load < SIM_MAX_LOAD; load++) {
		count_bp = 0;

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

		// print edges
		//for (int i = 0; i < NSF_NUM_EDGES; i++) {
		//	printf ("Load[%d] Edges[%4d] => v: %2d\tu: %2d\twave: %2d\ttimes: (", load, i, edges[i].v, edges[i].u, edges[i].wave);
		//	for (int c = 0; c < NSF_NUM_CHANNELS; c++)
		//		printf (" %.2f", edges[i].time[c]);
		//	printf (")\n");
		//}

#ifdef DEBUG_NSF
		printNsf (verts, edges);
#endif
		/* Processing GA */
		for (int sim = 0; sim < SIM_NUM_GEN; sim++) {
			//printf ("Simulation %3d", sim);
			struct Chromosome individual[GA_SIZE_POP];
			double until_next   = -log(1-(rand()%RAND_MAX+1)/(double)((unsigned)RAND_MAX+1))/load;
			double holding_time = -log(1-(rand()%RAND_MAX+1)/(double)((unsigned)RAND_MAX+1));
			//printf (" [until_next: %f | holding_time: %f]", until_next, holding_time);
#ifdef DEBUG_FITNESS
			int fitness[GA_MIN_GEN];
#endif
			/* Generating individuals */
		        for (int i = 0; i < GA_SIZE_POP; i++) 
                		addIndividual (verts, edges, &individual[i], i, NSF_SRC_NODE, NSF_DST_NODE);
#ifdef DEBUG_POP
        		printPopulation (individual, PREVIEW_AVAIL);
#endif
			//printf ("\n GA[");
			for (int gen = 0; gen < GA_MIN_GEN; gen++) {
				//printf (".");
				/* perform evaluate */
				evaluate (individual, verts, edges);

				/* perform selection */
				int *mating_pool = selection (individual, GA_MAX_CROSS_RATE);
				
				/* perform cross */
				if (mating_pool) {
					//printf ("s"); 
					crossover (individual, mating_pool);
				}

				/* perform mutation */
				struct Chromosome *chromosome = NULL;
				for (int i = 0; i < GA_SIZE_POP; i++) {
					if ((rand() % 100) > GA_MIN_MUT_RATE) continue;
#ifdef DEBUG_MUTATE
                                        printf ("\nMutate [%d]:\n---\n", i);
#endif
					chromosome = mutate (verts, edges, &individual[i]);
					//printf ("m");

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
#ifdef DEBUG_MUTATE
                        		else printf ("Mutation: Returned NULL!\n\n");
#endif
				} // Mutate

				/* sorting by length of route */
				insertionSort (individual, BY_LENGTH);
				/* sorting by wavelength availability */
				insertionSort (individual, BY_WL_AVAIL);
				//if (sim == 0 && load == 1 && gen == 1) printPopulation (individual, PREVIEW_AVAIL);
#ifdef DEBUG_FITNESS
			if (load == 1 && sim == 1) {
				fitness[gen] = 0;
				for (int i = 0; i < GA_SIZE_POP; i++) {
					if (individual[i].wl_avail > 0)
						fitness[gen] += 1;
				}
				//printf ("fitness[%d]: %d\n", gen, fitness[gen]);
				//printf ("%d\n", fitness[gen]);
			}
#endif
			} // GA_MIN_GEN
			//printf ("] OK\n");

			/* Final Evaluation - Fitness calculation */
			evaluate (individual, verts, edges);
			insertionSort (individual, BY_LENGTH);
			insertionSort (individual, BY_WL_AVAIL);
			//printPopulation (individual, PREVIEW_AVAIL);

			int color;

			if (individual[0].wl_avail > 0) {
					for (color = 0; color < NSF_NUM_CHANNELS; color++) 
							if (((individual[0].L & (1<<color)) >> color)) break; 

					for (int i = 0; i < individual[0].length-1; i++) {
							for (int j = 0; j < verts[i].adjLen; j++) {
									int e = verts[i].edges[j];

									if (edges[e].u != individual[0].gene[i+1]) continue;

									edges[e].wave -= pow (2, color);
									edges[e].time[color] = holding_time;

									for (int z = 0; z < verts[i+1].adjLen; z++) {
											e = verts[i+1].edges[z];
											if (verts[e].label != individual[0].gene[i]) continue;
										edges[e].wave -= pow (2, color);
										edges[e].time[color] = holding_time;
									}

									break;
							}
					}
			} else count_bp += 1;

			// Updating traffic matrix
			for (int i = 0; i < NSF_NUM_NODES; i++) {
					for (int e = 0; e < verts[i].adjLen; e++) {
							int link = verts[i].edges[e];

							if (edges[link].u < i) {
								// Making Symetric
								for (int ee = 0; ee < verts[edges[link].u].adjLen; ee++) {
									int link_cp = verts[edges[link].u].edges[ee];

									if (edges[link_cp].u != i) continue;

									for (int c = 0; c < NSF_NUM_CHANNELS; c++)
											edges[link].time[c] = edges[link_cp].time[c];
								}
							} else {
								for (int c = 0; c < NSF_NUM_CHANNELS; c++) {
										if (edges[link].time[c] > until_next) {
												edges[link].time[c] -= until_next;
										} else {
												edges[link].time[c] = 0.0;
												if ( !((edges[link].wave & (1<<c)) >> c) ) {
														edges[link].wave += (int) pow (2, c);
												}
										}
								}
							}

					}
			}
			//break;
		} // SIM_NUM_GEN
		//if (load == 1) printNsf(verts, edges);
		//break;
		float bp = 100.0 * (count_bp/(float)SIM_NUM_GEN);
		printf ("%.1f\n", bp);
		//printf ("Block Probability [Erlang %2d]: %.1f%%\n", load, bp);
		//printf ("Blocked: %d\nBlock Probability: %.1f%%\n", count_bp, bp);
	} // SIM_MAX_LOAD

	t = clock() - t;
	double time_taken = ((double) t)/CLOCKS_PER_SEC;

	printf ("Time taken: %.2f\n", time_taken);

	return EXIT_SUCCESS;
}
