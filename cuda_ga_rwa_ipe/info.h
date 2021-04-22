/*************************
 * Arquivo de parametros
 ************************/

#define FALSE 			0
#define TRUE  			1
#define BIGN			32767	// RAND_MAX

#define BY_LENGTH 		0
#define BY_WL_AVAIL 		1
#define PREVIEW_AVAIL 		2
#define FINAL_AVAIL 		3

/* Simulation Parameters */
//#define DEBUG			TRUE
//#define DEBUG_NSF		TRUE
//#define DEBUG_POP		TRUE
//#define DEBUG_ADDIND		TRUE
//#define DEBUG_SELECT		TRUE
//#define DEBUG_CROSS		TRUE
//#define DEBUG_MUTATE		TRUE

#define SIM_NUM_GEN		150
#define SIM_MIN_LOAD		1
#define SIM_MAX_LOAD		65
#define SIM_RAND_MAX		9000

/* NSF Parameters */
#define NSF_SRC_NODE		0 	// source node
#define NSF_DST_NODE		18	// destination node
#define NSF_NUM_NODES		28 	// number of nodes on NSF network
#define NSF_NUM_EDGES		80 	// number of edges on NSF graph
#define NSF_NUM_CHANNELS	4

#define FIELDS	2
//#define NODE_VISITED		0
//#define NODE_WEIGHT		1

/* Genetic Algorithm Parameters */
#define GA_SIZE_POP		30 	// size of population of each species

#define GA_MIN_GEN		35 	// min number of generations
#define GA_MAX_GEN		80 	// max number of generations

#define GA_MIN_CROSS_RATE	0.15 	// min crossover rate
#define GA_MAX_CROSS_RATE	40 	// max crossover rate

#define GA_MIN_MUT_RATE		2 	// min mutation rate
#define GA_MAX_MUT_RATE		20 	// max mutation rate

#define GA_GEN_INTERVAL		8 	// interval to update rates
