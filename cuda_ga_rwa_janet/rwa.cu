
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include "nsf.h"
#include "rwa.h"
#include "info.h"
#include "util.cuh"

/* device */
__device__ void addIndividual (struct Vertice *verts, struct Edge *edges, struct Chromosome *individual, int idx, int src, int dst)
{
#ifdef DEBUG_ADDIND
	int myId = threadIdx.x + blockDim.x * blockIdx.x;

	if (myId == 0)
	printf (">> Adding individual[%d]\n", idx);
#endif
        resetNsfVisited (verts, NSF_NUM_NODES);

        int route[NSF_NUM_NODES], route_next, route_curr = src;

        route[0] = route_curr;
        individual->length = 1;
        individual->S = FALSE;
        verts[route_curr].isVisited = TRUE;
	
        for (int i = 1; i < NSF_NUM_NODES;) {
		int notVisited = 0;
		int adj[NSF_NUM_NODES];

#ifdef DEBUG_ADDIND
	if (myId == 0)
		printf ("route[%d]: %d => Adjacencies (not visited): { ", i, route[i-1]);
#endif
                for (int j = 0; j < verts[route_curr].adjLen; j++) {
			int e = edges[verts[route_curr].edges[j]].u;
			if (verts[e].isVisited == TRUE) continue;
			adj[notVisited] = verts[e].label;
#ifdef DEBUG_ADDIND
	if (myId == 0)
			printf ("%d ", adj[notVisited]);
#endif
			notVisited++;
		}

		if (notVisited == 0) {
#ifdef DEBUG_ADDIND
	if (myId == 0)
			printf ("No adjacencies left }\n");
#endif
			route_curr = route[i-2];

			if (route_curr == src) {
#ifdef DEBUG_ADDIND
	if (myId == 0)
				printf (">> Backtracked to src. Terminating!\n");
#endif
				break;
			}
#ifdef DEBUG_ADDIND
	if (myId == 0)
			printf (">> Backtracked to route[%d]: %d ..\n", i-1, route_curr);
#endif
			i--;
                        individual->length--;
			continue;
		}
#ifdef DEBUG_ADDIND
		else {
	if (myId == 0)
			printf ("}\n");
		}
#endif

                for (int j = 0; j < notVisited; j++) {
#ifdef DEBUG_ADDIND 
	if (myId == 0)
			printf (">> Try %d/%d sorting next hop.. ", j+1, notVisited);
#endif
                        route_next = adj[(int)dRand(notVisited)];
#ifdef DEBUG_ADDIND
	if (myId == 0)
                        printf ("%d\n", route_next);
#endif
                        if (verts[route_next].isVisited == TRUE) continue;

                        verts[route_next].isVisited = TRUE;
                        individual->length++;
                        route_curr = route_next;
                        route[i] = route_next;

                        break;
                }

                if (route_curr == dst) break;

		i++;
        }
#ifdef DEBUG_ADDIND
	if (myId == 0)
	printf (">> Inserting individual[%d] in population.. ", idx);
#endif
        for (int i = 0; i < individual->length; i++)
                individual->gene[i] = route[i];

        individual->L = 0;
        individual->wl_avail = 0;
#ifdef DEBUG_ADDIND
	if (myId == 0)
	printf ("Done!\n");
#endif
}

__device__ void evaluate (struct Chromosome* individual, struct Vertice* verts, struct Edge* edges)
{
        float L[NSF_NUM_CHANNELS];
        for (int i = 0; i < NSF_NUM_CHANNELS; i++)
                L[i] = 0.0;

        for (int i = 0; i < GA_SIZE_POP; ++i) {
                for (int j = 1; j < NSF_NUM_CHANNELS+1; ++j) {
                        int num = 0;

                        //printf ("Evaluating individual[%d]:\n", i);
                        for (int k = 0; k < individual[i].length-1; k++) {
                                int rcurr = individual[i].gene[k];
                                int rnext = individual[i].gene[k+1];
                                //printf ("%d -> %d [%f]\n", rcurr, rnext);

                                for (int n = 0; n < verts[rcurr].adjLen; n++) {
                                        int e = verts[rcurr].edges[n];

                                        if (edges[e].u != rnext)
                                                continue;

                                        if (edges[e].wave & (1 << (j-1)))
                                                num += j;
                                }
                        }

                        L[j-1] = (float) num/(float)(j*(individual[i].length-1));
                }

                individual[i].wl_avail = 0;
                individual[i].L = 0;

                for (int j = 0; j < NSF_NUM_CHANNELS; ++j) {
                        if (L[j] == 1.0) {
                                individual[i].wl_avail += 1;
                                individual[i].L += 1<<j;
                        }
                }
        }
}

__device__ int *selection (struct Chromosome* individual, int cross_rate)
{
        int count = 0;
        int candidate[2];
        int selected[GA_SIZE_POP];

        for (int i = 0; i < GA_SIZE_POP; i++) {
                selected[i] = 0;
                candidate[0] = (int)dRand(GA_SIZE_POP);

                if ((int)dRand(100) < cross_rate) {
                        for (int j = 0; j < 3; j++) {
                                candidate[1] = (int)dRand(GA_SIZE_POP);

                                if (candidate[0] == candidate[1]) continue;

                                if (individual[candidate[0]].wl_avail >= individual[candidate[1]].wl_avail) {
                                        selected[count] = candidate[0];
					individual[candidate[0]].S = TRUE;
                                } else {
					selected[count] = candidate[1];
					individual[candidate[1]].S = TRUE;
				}
                        }
                }

                if (count > 0) {
                        for (int j = 0; j < count; j++)
                                if (selected[j] == selected[count]) --count;
                }

                ++count;
        }

        if (count == 0) return NULL;

#ifdef DEBUG_SELECT
	int myId = threadIdx.x + blockDim.x * blockIdx.x;

	if (myId == 0) {
        	printf ("\nSelection (sum: %d):\n---\n", count);
        	for (int i = 0; i < count; i++) 
        	        printf ("Parent %d: %d\n", i, selected[i]);
	}
#endif
        return selected;
}

__device__ void crossover (struct Chromosome *individual, int *mating_pool)
{
        struct Chromosome children[2];
        int common[NSF_NUM_NODES], dad = 0, mom = 0, idx = 0, count = 0, S = 0;

	for (int i = 0; i < GA_SIZE_POP; i++)
		if (individual[i].S == TRUE) S++;

        for (int i = 0; i < S; i++) {
                dad = mating_pool[(int)dRand(S)];
                mom = mating_pool[(int)dRand(S)];

                for (int j = 0; j < 10; j++) {
                        if (dad == mom) {
                                mom = (int)dRand(S);
                                continue;
                        }
                        break;
                }

                for (int j = 1; j < individual[dad].length-1; ++j) {
                        for (int k = 1; k < individual[mom].length-1; ++k) {
                                if (individual[dad].gene[j] == individual[mom].gene[k]) {
                                        common[count] = individual[dad].gene[j];
                                        count++;
                                }
                        }
                }
                // If don't find a common node between
                // dad and mom, try again
                ///
                if (!count) continue;

                idx = common[(int)dRand(count)];

                break;
        }
#ifdef DEBUG_CROSS
	int myId = threadIdx.x + blockDim.x * blockIdx.x;

	if (myId == 0) {
        	printf ("\nCrossover:\n---\nDad: %d\nMom: %d\n\nCommon nodes: ", dad, mom);
        	for (int i = 0; i < count; i++) {
        	        if (common[i] == idx)
        	                printf (" [%d] ", common[i]);
        	        else    printf (" %d ", common[i]);
        	}
        	printf ("\n");
	}
#endif
        // Children generation 
        children[0].length = 0;
        children[0].L = 0;
        children[0].wl_avail = 0;
        children[1].length = 0;
        children[1].wl_avail = 0;
        children[1].L = 0;

        for (count = 0; individual[dad].gene[count] != idx; ++count) {
                        children[0].length++;
                        children[0].gene[count] = individual[dad].gene[count];
        }

        children[0].gene[count] = individual[dad].gene[count];
        children[0].length++;
        int idx0 = children[0].length;

        for (count = 0; individual[mom].gene[count] != idx; ++count) {
                        children[1].length++;
                        children[1].gene[count] = individual[mom].gene[count];
        }

        children[1].gene[count] = individual[mom].gene[count];
        children[1].length++;
        int idx1 = children[1].length;

        if ((idx0 + (individual[mom].length - idx1)) <= NSF_NUM_NODES) {
                int idxmom = idx1;
                for (count = idx0; count < (idx0 + (individual[mom].length - idx1)); ++count) {
                        children[0].gene[count] = individual[mom].gene[idxmom];
                        children[0].length++;
                        idxmom++;
                }
        }

        if ((idx1 + (individual[dad].length - idx0)) <= NSF_NUM_NODES) {
                int idxdad = idx0;
                for (count = idx1; count < (idx1 + (individual[dad].length - idx0)); ++count) {
                        children[1].gene[count] = individual[dad].gene[idxdad];
                        children[1].length++;
                        idxdad++;
                }
        }
        // Children generation ends 

#ifdef DEBUG_CROSS
	if (myId == 0) {
        printf ("\n");
        for (int j = 0; j < 2; ++j) {
                printf ("Child %d: ", j);
                for (int i = 0; i < NSF_NUM_NODES; ++i) {
                        if (children[j].gene[i] == NSF_DST_NODE) {
                                printf (" %d", children[j].gene[i]);
                                break;
                        }
                        printf (" %d ->", children[j].gene[i]);
                }
                printf ("\n");
        }}
#endif
        // Testing and inserting child in population 
        for (count = 0; count < 2; count++) {
                int ok = TRUE;

                if (children[count].gene[children[count].length-1] != NSF_DST_NODE) continue; 

                for (int i = 0; i < children[count].length-1; i++) {
                        for (int j = i+1; j < children[count].length; j++) {
                                if (i == (i+j)) continue;

                                if (children[count].gene[i] == children[count].gene[j]) {
#ifdef DEBUG_CROSS
	if (myId == 0) 
                                        printf ("DUP: child[%d] => %d\n", count, children[count].gene[i]);
#endif
                                        ok = FALSE;
                                        break;
                                }
                        }
                        if (ok == FALSE) break;
                }
                if (ok == FALSE) continue;
#ifdef DEBUG_CROSS
		if (myId == 0) 
                printf ("Inserting child %d in population..\n", count);
#endif
		//children[count].idx = GA_SIZE_POP - (count+1);
                individual[GA_SIZE_POP - (count+1)] = children[count];
        }
}

__device__ struct Chromosome *mutate (struct Vertice* verts, struct Edge* edges, struct Chromosome* ind)
{
        if (ind->length == 2) return NULL;

        struct Chromosome newInd;
        int route[NSF_NUM_NODES];

        //int idx = (int)dRand(ind->length-1);
        int idx = (int)dRand(ind->length);
#ifdef DEBUG_MUTATE
	int myId = threadIdx.x + blockDim.x * blockIdx.x;
	if (myId ==0) {
        	printf ("Individual:");
        	for (int i = 0; i < ind->length; ++i) {
        	        if (i == ind->length-1) {
        	                printf (" %d", ind->gene[i]);
        	                continue;
        	        }

        	        if (i == idx)
        	                printf (" [%d] ->", ind->gene[i]);
        	        else    printf (" %d ->", ind->gene[i]);
        	}
        	printf ("\n");
	}
#endif
        resetNsfVisited (verts, NSF_NUM_NODES);
        newInd.length = 0;

        for (int i = 0; i <= idx; ++i) {
#ifdef DEBUG_MUTATE
	if (myId ==0)
                printf ("Marking gene %d as visited..\n", ind->gene[i]);
#endif
                verts[ind->gene[i]].isVisited = TRUE;
                route[i] = ind->gene[i];
                newInd.length++;
        }

        int adj[NSF_NUM_NODES], notVisited, route_next, route_curr = route[newInd.length-1];

        // Like addIndividual()
        for (int i = newInd.length; i < NSF_NUM_NODES;) {
		notVisited = 0;

                for (int j = 0; j < verts[route_curr].adjLen; j++) {
                        //int e = nsf->vert[route_curr].edges[j];
                        int e = edges[verts[route_curr].edges[j]].u;
                        if (verts[e].isVisited == TRUE) continue;
                        adj[notVisited] = verts[e].label;
                        notVisited++; 
                }

		if (notVisited == 0) {
			route_curr = route[i-2];

			if (route_curr == NSF_SRC_NODE) break;

			i--;
			newInd.length--;
			continue;
		}

                for (int j = 0; j < notVisited; j++) {
                        route_next = adj[(int)dRand(notVisited)];

                        if (verts[route_next].isVisited == TRUE) continue;

                        verts[route_next].isVisited = TRUE;
                        newInd.length++;
                        route_curr = route_next;
                        route[i] = route_next;

                        break;
                }

                if (route_curr == ind->gene[ind->length-1]) break;
		i++;
        }

        for (int i = 0; i < newInd.length; i++)
                newInd.gene[i] = route[i];

        if (newInd.gene[newInd.length-1] != ind->gene[ind->length-1]) {
#ifdef DEBUG_MUTATE
	if (myId ==0)
                printf ("Mutation: Destination node can't reached. Aborting! [1]\n");
#endif
                return NULL;
        }

        // Checking if newInd is different of ind 
        if (newInd.length == ind->length) {
                int eq = 0;

                for (int i = 0; i < newInd.length; ++i)
                        for (int j = 0; j < ind->length; ++j)
                                if (newInd.gene[i] == ind->gene[j]) ++eq;

                if (eq == newInd.length) {
#ifdef DEBUG_MUTATE
	if (myId ==0)
                        printf ("Mutation: Individual mutated aren't different. Aborting! [2]\n");
#endif
                        return NULL;
                }
        }
#ifdef DEBUG_MUTATE
	if (myId ==0) {
        	printf ("Mutation:");
        	for (int i = 0; i < newInd.length; ++i) {
        	        if (newInd.gene[i] == ind->gene[ind->length-1]) {
        	                printf (" %d", newInd.gene[i]);
        	                continue;
        	        }
        	        printf (" %d ->", newInd.gene[i]);
        	}
        	printf ("\n\n");
	}
#endif
        newInd.L = 0;
        newInd.wl_avail = 0;

        return &newInd;
}

__device__ void insertionSort (struct Chromosome* individual, int field)
{
        int j;
        struct Chromosome key;

        for (int i = 1; i < GA_SIZE_POP; i++) {
                j = i-1;

                key.idx = individual[i].idx;
                key.length = individual[i].length;
                key.wl_avail = individual[i].wl_avail;
                key.L = individual[i].L;
                key.S = individual[i].S;

                for (int k = 0; k < key.length; ++k)
                        key.gene[k] = individual[i].gene[k];

                if (field == BY_LENGTH) {
                        while (j >= 0 && individual[j].length > key.length) {
                                individual[j+1].idx = individual[j].idx;
                                individual[j+1].length = individual[j].length;
                                individual[j+1].wl_avail = individual[j].wl_avail;
                                individual[j+1].L = individual[j].L;
                                individual[j+1].S = individual[j].S;
                                for (int k = 0; k < individual[j+1].length; ++k)
                                        individual[j+1].gene[k] = individual[j].gene[k];

                                j = j-1;
                        }
                }

                if (field == BY_WL_AVAIL) {
                        while (j >= 0 && individual[j].wl_avail < key.wl_avail) {
                                individual[j+1].idx = individual[j].idx;
                                individual[j+1].length = individual[j].length;
                                individual[j+1].wl_avail = individual[j].wl_avail;
                                individual[j+1].L = individual[j].L;
                                individual[j+1].S = individual[j].S;
                                for (int k = 0; k < individual[j+1].length; ++k)
                                        individual[j+1].gene[k] = individual[j].gene[k];

                                j = j-1;
                        }
                }

                individual[j+1].idx = key.idx;
                individual[j+1].length = key.length;
                individual[j+1].wl_avail = key.wl_avail;
                individual[j+1].L = key.L;
                individual[j+1].S = key.S;
                for (int k = 0; k < key.length; ++k)
                        individual[j+1].gene[k] = key.gene[k];
        }
}

__device__ void printPopulation (struct Chromosome *individual, int type, int myId)
{
	if (myId != 0) return;

	printf ("\nPopulation:\n---\n");
	for (int i = 0; i < GA_SIZE_POP; ++i){
	        printf ("[%02d/%d](%2d)(wl: %d)[%2d]:", i, GA_SIZE_POP-1, individual[i].length, individual[i].wl_avail, individual[i].L);
	        for (int j = 0; j < individual[i].length; ++j) {
	                if (individual[i].gene[j] != NSF_DST_NODE) {
	                        printf (" %d ->", individual[i].gene[j]);
	                        continue;
	                }
	                printf (" %d", individual[i].gene[j]);
	        }
		if (type == FINAL_AVAIL && i == 0)
			printf (" *Best route*");
	        printf ("\n");
	}
}
