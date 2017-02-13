typedef struct kdtree {
	//int n_dimensions;
	struct node* root;
	int size;
	double *x, *y, *z;
	//double * dim_data[10];
} kdtree_t;

typedef struct node {
	double x, y, z;
	int size;
	double xmin, xmax, ymin, ymax, zmin, zmax;
	struct node * lchild;
	struct node * rchild;
} node_t;

