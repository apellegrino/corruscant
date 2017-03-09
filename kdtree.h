#define FLOAT double

enum dim;

typedef struct kdtree {
	struct node* root;
	int size;
	FLOAT *x, *y, *z;
} kdtree_t;

typedef struct node {
	double x, y, z;
	//double x,y,z;
	struct node * lchild;
	struct node * rchild;
} node_t;

/*
static inline double max(double, double);
static inline double min(double, double);
static inline double max3(double, double, double);
static inline double min3(double, double, double);
*/
static inline FLOAT norm2(FLOAT, FLOAT, FLOAT);

static inline void swapDouble(FLOAT *, FLOAT *);
static inline void swapInt(int *, int *);
void quicksort(FLOAT *, int *, int, int);
int * argsort(FLOAT *, int);

void partition(FLOAT *, int *, FLOAT, int, int, int, int);
node_t * build(FLOAT *, FLOAT *, FLOAT *,
					int *, int *, int *,
					int, int, enum dim);
void verify(node_t *, enum dim);
void destroy(node_t *);
int radius(node_t *, enum dim, FLOAT, FLOAT, FLOAT, FLOAT);
kdtree_t tree_construct(int, FLOAT [], FLOAT [], FLOAT []);
/*
long long two_point_correlation(kdtree_t tree, double [], double [],
									double [], int, double, MPI_Comm);
*/
