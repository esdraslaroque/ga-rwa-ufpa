struct Chromosome {
	int idx;
        int length;
        int wl_avail;
        int L;
        int S;
        int gene[NSF_NUM_NODES];
};

__device__ void addIndividual (struct Vertice *verts, struct Edge *edges, struct Chromosome *individual, int idx, int src, int dst);

__device__ void evaluate (struct Chromosome* individual, struct Vertice* verts, struct Edge* edges);

__device__ int *selection (struct Chromosome* individual, int cross_rate);

__device__ void crossover (struct Chromosome *individual, int *mating_pool);

__device__ struct Chromosome *mutate (struct Vertice* verts, struct Edge* edges, struct Chromosome* ind);

__device__ void insertionSort (struct Chromosome* individual, int field);

__device__ void printPopulation (struct Chromosome *individual, int type, int myId);

