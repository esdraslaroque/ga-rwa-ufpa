#include "info.h"
struct Edge {
        int v;
        int u;
        int wave;
        double time[NSF_NUM_CHANNELS];
};

struct Vertice {
	int label;
        int adjLen;
        int isVisited;
	int edges[NSF_NUM_NODES];
};

struct Nsf {
        int V;
        int E;
        struct Vertice *vert;
        struct Edge *edge;
};

struct NodeList *newNodeList (int node, struct Nsf *nsf);

void addLink (struct Vertice *verts, struct Edge *edges, int src, int dst, int wave, int nc, int *E, int ctl);

//void resetNsfVisited (struct Nsf* nsf);
void resetNsfVisited (struct Vertice* verts, int V);

//void updateTime (struct Nsf *nsf, int vert1, int vert2, double until_next);

const char *bin (int num, int nbits);
void printNsf (struct Vertice *verts, struct Edge *edges);

