struct Chromosome {
	int idx;
        int length;
        int wl_avail;
        int L;
        int S;
        int gene[NSF_NUM_NODES];
};

void addIndividual (struct Vertice *verts, struct Edge *edges, struct Chromosome *individual, int idx, int src, int dst);

void evaluate (struct Chromosome* individual, struct Vertice* verts, struct Edge* edges);

int *selection (struct Chromosome* individual, int cross_rate);

void crossover (struct Chromosome *individual, int *mating_pool);

struct Chromosome *mutate (struct Vertice* verts, struct Edge* edges, struct Chromosome* ind);

void insertionSort (struct Chromosome* individual, int field);

void printPopulation (struct Chromosome *individual, int type);

