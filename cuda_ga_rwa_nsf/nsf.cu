#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nsf.h"
#include "info.h"

__device__ void resetNsfVisited (struct Vertice* verts, int V)
{
        for (int i = 0; i < V; i++)
                verts[i].isVisited = FALSE;
}

__device__ struct Nsf *d_createNsf (int V, int E)
{
        struct Nsf* nsf = (struct Nsf*) malloc(sizeof(struct Nsf));
	
        nsf->V = V;
        nsf->E = E;

	nsf->vert = (struct Vertice*) malloc(nsf->V * sizeof(struct Vertice));
	nsf->edge = (struct Edge*) malloc(nsf->E * sizeof(struct Edge));

        for (int i = 0; i < nsf->V; i++) {
                nsf->vert[i].label = i;
		nsf->vert[i].adjLen = 0;
                nsf->vert[i].isVisited = FALSE;

		for (int j = 0; j < nsf->V; j++) {
                	nsf->vert[i].edges[j] = BIGN;
		}
        }

        return nsf;
}

void addLink (struct Vertice *verts, struct Edge *edges, int src, int dst, int wave, int nc, int *E, int ctl)
{
        //printf ("Adding nsf->edge[%d]: %d -> %d\n", nsf->E, src, dst);
        //if (ctl == 1) printf ("\n");
        int nE = *E;

        edges[nE].v = src;
        edges[nE].u = dst;
        edges[nE].wave = wave;

        // Init for traffic matrix by wavelength
        for (int i = 0; i < nc; ++i)
                edges[nE].time[i] = 0.0;

        verts[src].edges[verts[src].adjLen] = nE;
        verts[src].adjLen++;

        nE++; *E = nE;
        ctl--;

        if (ctl == 0) return;

        addLink (verts, edges, dst,  src, wave, nc, E, ctl);

}
/*
void updateTime (struct Nsf *nsf, int vert1, int vert2, double until_next)
{
	struct NodeList *node1;
	struct NodeList *node2;

	node1 = nsf->array[vert1].head;
	node2 = nsf->array[vert2].head;

	for (int i = 0; i < NSF_NUM_CHANNELS; i++) {

		if (node1->time[i] > until_next) {
			node1->time[i] -= until_next;
			node2->time[i] -= until_next;
		} else {
			node1->time[i] = 0.0;
			node2->time[i] = 0.0;

			if ( !((node1->wave & (1<<i)) >> i) ) {
				node1->wave += (int)pow(2,i);
				node2->wave += (int)pow(2,i);
			}
		}
	}
}
*/
const char *bin (int num, int nbits)
{
        static char str[9] = {0};

        for (int i = nbits-1; i >= 0; i--) {
                str[i] = (num&1)?'1':'0';
                num >>= 1;
        }
        return str;
}

// Printing function by Host 
void printNsf (struct Vertice *verts, struct Edge *edges)
{
        printf ("\n[Host] NSF - Wavelength Network:\n---\nVertices: %d\nEdges: %d\n---\n", NSF_NUM_NODES, NSF_NUM_EDGES);
        for (int v = 0; v < NSF_NUM_NODES; ++v) {

                printf ("\nAdjacency for Node: %d/%d (sum: %d)\n", verts[v].label, NSF_NUM_NODES-1, verts[v].adjLen);
                for (int e = 0; e < verts[v].adjLen; e++) {
                        //printf ("%2d -> %2d \t[%3d](",        edges[verts[v].edges[e]].v,
                        //                              edges[verts[v].edges[e]].u,
                        //                              edges[verts[v].edges[e]].wave);
                        printf ("%2d -> %2d \t[%s](",   edges[verts[v].edges[e]].v,
                                                        edges[verts[v].edges[e]].u,
                                                        bin(edges[verts[v].edges[e]].wave, NSF_NUM_CHANNELS));

                        for (int i = 0; i < NSF_NUM_CHANNELS; i++) {
                                if (i == (NSF_NUM_CHANNELS-1)) {
                                        printf ("%.2f", edges[verts[v].edges[e]].time[i]);
                                        continue;
                                }

                                printf ("%.2f ", edges[verts[v].edges[e]].time[i]);
                        }

                        printf (")\n");
                }
        }
        printf ("\n");
}

// Printing function by Device
__device__ void d_printNsf (struct Vertice *verts, struct Edge *edges, int myId)
{
	if (myId != threadIdx.x + blockDim.x * blockIdx.x) return;

        printf ("\n[Device] NSF - Wavelength Network:\n---\nVertices: %d\nEdges: %d\n---\n", NSF_NUM_NODES, NSF_NUM_EDGES);
        for (int v = 0; v < NSF_NUM_NODES; ++v) {

                printf ("\nAdjacency for Node: %d/%d (sum: %d)\n", verts[v].label, NSF_NUM_NODES-1, verts[v].adjLen);
                for (int e = 0; e < verts[v].adjLen; e++) {
                        printf ("%2d -> %2d \t[%2d](", 	edges[verts[v].edges[e]].v,
							edges[verts[v].edges[e]].u,
							edges[verts[v].edges[e]].wave);

			for (int i = 0; i < NSF_NUM_CHANNELS; i++) {
				if (i == (NSF_NUM_CHANNELS-1)) {
					printf ("%.2f", edges[verts[v].edges[e]].time[i]);
					continue;
				}

				printf ("%.2f ", edges[verts[v].edges[e]].time[i]);
			}

			printf (")\n");
                }
        }
        printf ("\n");
}


