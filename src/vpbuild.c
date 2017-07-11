/*
 * Andrew Pellegrino, 2017
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

#include "vptree.h"

#define ROOT 1

static int _include;
static int _use_err;

static void count_node(field_counter_t * counter, int this_field, int query_field)
{
    // index 0 represents whole dataset
    (counter->array)[0]++;

    // no error
    if(!_use_err) return;

    // field-to-field error: increment counter for this_field if it
    // coincides with query_field
    if (_include) {
        if (this_field == query_field) (counter->array[this_field])++;
    }
    else { //jackknife error: increment all but this_field and query_field
        int i;
        // needs vectorization?
        for(i=1; i<counter->size; i++) { (counter->array[i])++; }
        /* 
         * un-count this pair for the querying point field or that of the
         * point on the tree, but don't double count if they have the same
         * field
         */
        (counter->array)[this_field]--;
        if (this_field != query_field) (counter->array)[query_field]--;
    }
    return;
}

static inline void swapPoint(vpoint_t * a, vpoint_t * b)
{
    vpoint_t temp;
    temp = *a;
    *a = *b;
    *b = temp;
}

inline int left_child(int p)
{
    return 2*p;
}

inline int right_child(int p)
{
    return 2*p+1;
}

void destroy(vptree_t t)
{
    free(t.node_data);
}

/* Flag node as having a left child */
static inline void set_lchild(node_t * node)
{
    node->flags |= HAS_LCHILD;
}

/* Flag node as having a right child */
static inline void set_rchild(node_t * node)
{
    node->flags |= HAS_RCHILD;
}

static inline int has_lchild(node_t p)
{
    return p.flags & HAS_LCHILD;
}

static inline int has_rchild(node_t p)
{
    return p.flags & HAS_RCHILD;
}

double vincenty(vpoint_t a, vpoint_t b)
{
    double c, d, e;
    double dlon = a.lon - b.lon;
    c = cos(b.lat) * sin(dlon);
    d = cos(a.lat) * sin(b.lat) - sin(a.lat) * cos(b.lat) * cos(dlon);
    e = sin(a.lat) * sin(b.lat) + cos(a.lat) * cos(b.lat) * cos(dlon);
    return atan2( sqrt( c*c + d*d ), e );
}

double great_circle(vpoint_t a, vpoint_t b)
{
    return acos( sin(a.lat) * sin(b.lat)
                + cos(a.lat) * cos(b.lat) * cos(a.lon - b.lon) );
}

double distance(vpoint_t a, vpoint_t b)
{
    if (a.lat == b.lat && a.lon == b.lon) return 0.0;

    double c, d, e;
    double dlon = a.lon - b.lon;
    c = cos(b.lat) * sin(dlon);
    d = cos(a.lat) * sin(b.lat) - sin(a.lat) * cos(b.lat) * cos(dlon);
    e = sin(a.lat) * sin(b.lat) + cos(a.lat) * cos(b.lat) * cos(dlon);
    return atan2( sqrt( c*c + d*d ), e );

    // chord length method
    /*
    double dx = cos(b.lat) * cos(b.lon) - cos(a.lat) * cos(a.lon);
    double dy = cos(b.lat) * sin(b.lon) - cos(a.lat) * sin(a.lon);
    double dz = sin(b.lat) - sin(a.lat);
    return 2 * asin( sqrt(dx*dx+dy*dy+dz*dz) / 2.0 );
    */
}


/*
 * Shuffle `list` so that the nth element is the same as if `list` were sorted.
 * Sorting is based on distance to `vantage`. Uses quickselect, which acts like
 * quicksort except function only recurses to the side that n is in
 */
void nth_element(vpoint_t * list, int left, int right, vpoint_t vantage, int n)
{
    if ( right <= left ) return;

    vpoint_t pivot = list[left];
    int i = left+1;
    int j = right;

    while(1) {
        while( distance(vantage,list[i]) <= distance(vantage,pivot) && i <= right ) i++;
        while( distance(vantage,list[j]) >  distance(vantage,pivot) ) j--;

        if ( i >= j ) break;
        swapPoint(list+i, list+j);
    }
    swapPoint(&list[left], &list[j]);

    if (j > n) {
        nth_element(list, left, j-1, vantage, n);
    } else if (j < n) {
        nth_element(list, j+1, right, vantage, n);
    }
    return;
}

// round up a positive integer x to the next power of 2
static int pow2ceil(int x)
{
    int i = 1;
    while(x > 1) { x /= 2; i++; }
    while(i > 0) { x *= 2; i--; }
    return x;
}

void build(vptree_t * tree, int node_id, vpoint_t * data, int left, int right)
{
    node_t * node = &((tree->node_data)[node_id]);
    // select root -- 0th for now
    node->point = data[left];
    if (left == right) {
        node->r = 0.0;
        return;
    }

    if (right - left == 1) {
        set_rchild(node);
        node_t * rc = &((tree->node_data)[right_child(node_id)]);
        rc->point = data[right];
        rc->r = 0.0;
        node->r = distance(node->point, rc->point);
        return;
    }

    int median = ((left+1)+right+1)/2;
    nth_element(data, left+1, right, node->point, median);
    node->r = distance(node->point, data[median]);

    set_lchild(node);
    set_rchild(node);
    build(tree, left_child(node_id), data, left + 1, median - 1);
    build(tree, right_child(node_id), data, median, right);

    return;
}

/*
void verify(node_t * data, int index)
{
    
    node_t *parent, *lc, *rc;
    //char *errmsg = "node %p should not be %s child of %p\n";

    parent = data + index;

	if (parent->flags & HAS_LCHILD) {
        lc = data + left_child(index);
		verify(data,left_child(index));
	}
	if (parent->flags & HAS_RCHILD) {
        rc = data + right_child(index);
		verify(data,right_child(index));
	}

}
*/

vptree_t make_vp_tree(double * lat, double * lon, int * fields, int length)
{
    vptree_t tree;
    tree.size = length;
    tree.memsize = pow2ceil(length);

    // allocate all space for the tree in one malloc
    tree.node_data = (node_t *) calloc( tree.memsize, sizeof(node_t) );

    // put data into list of points for easy manipulation
    vpoint_t * points = (vpoint_t *) malloc( length * sizeof(vpoint_t) );

    int i;
    for(i=0; i<length; i++) {
        points[i].lat = lat[i];
        points[i].lon = lon[i];
    }

    if(fields != NULL) {
        for(i=0; i<length; i++) {
            points[i].field = fields[i];
        }
    }

    build(&tree, 1, points, 0, length-1);

    free(points);
    return tree;
}

void radq(vptree_t tree, field_counter_t * counter, int pi, vpoint_t q,
          double radius)
{
    node_t * node = tree.node_data + pi;
    vpoint_t p = node->point;

    if( distance(p, q) < radius )
        count_node(counter, p.field, q.field);

    if(!has_rchild(*node)) return;
    if(!has_lchild(*node)) {
        radq(tree, counter, right_child(pi), q, radius);
        return;
    }

    if(distance(p, q) + radius < node->r) {
        radq(tree, counter, left_child(pi), q, radius);
    } else if(distance(p, q) - radius >= node->r) {
        radq(tree, counter, right_child(pi), q, radius);
    } else {
        radq(tree, counter, left_child(pi), q, radius);
        radq(tree, counter, right_child(pi), q, radius);
    }
    return;
}

typedef struct shared_args {
    vpdata_t data;
    int data_length;
    //int node_start;
    //int node_stop; // actually one after the last index to stop at
    int num_fields;
    int num_threads;
    field_counter_t ** counters;
    double radius;
    vptree_t tree;
} shared_args_t;

typedef struct thread_args {
    int thread_rank;
    struct shared_args * ptr;
} thread_args_t;

static inline int min(int a, int b)
{
    return a < b ? a : b;
}

void assign_idx(int rank, int size, int n_array, int * start, int * stop)
{
    int n_array_local = n_array / size + ( n_array % size > rank );

    *start = rank * ( n_array / size ) + min( rank, n_array % size );
    *stop = *start + n_array_local;
}

void * pair_count_wrap(void *voidargs)
{
    // thread_args_t: Overall args
    // shared_args_t: Args which are common to all threads
    thread_args_t *targs;
    targs = (thread_args_t *) voidargs;

    int start, stop, i;
    vpoint_t q;

    shared_args_t * args;
    args = targs->ptr;

    vpdata_t data = args->data;
    vptree_t tree = args->tree;
    double radius = args->radius;
    
    // divide up data array equally among threads by getting start, stop values
    // based on rank
    int rank = targs->thread_rank;
    assign_idx(rank, args->num_threads, args->data.length, &start, &stop);

    // make sure to free later
    //field_counter_t * new_counter = malloc(sizeof(field_counter_t));
    //new_counter->array = calloc(args->num_fields+1, sizeof(long long));
    //new_counter->size = args->num_fields+1;
    //args->counters[rank] = new_counter;
    field_counter_t * counter = args->counters[rank];
   
    for(i=start; i<stop; i++) {
        q.lat = (data.lat)[i];
        q.lon = (data.lon)[i];
        if(data.fields != NULL) {
            q.field = (data.fields)[i];
        }
        radq(tree, counter, ROOT, q, radius);
    }

    return NULL;
}

long long * pair_count(vptree_t tree, vpdata_t data, double radius, int num_threads)
{
    int i, j;
    int num_fields = data.num_fields;

    long long * results;
    results = calloc(num_fields+1, sizeof(long long));

    if(num_threads > 1<<16) {
        fprintf(stderr, "More than %d threads!\n", 1<<16);
        exit(1);
    }
    
    shared_args_t ss;

    ss.counters = (field_counter_t **) malloc(num_threads * sizeof(field_counter_t *));
    field_counter_t * new_counter;

    for(i=0; i < num_threads; i++) {
        new_counter = malloc(sizeof(field_counter_t));
        new_counter->size = num_fields+1;
        new_counter->array = calloc(new_counter->size, sizeof(long long));
        ss.counters[i] = new_counter;
    }

    ss.data = data;
    ss.tree = tree;
    ss.data_length = data.length;
    ss.num_fields = data.num_fields;
    ss.radius = radius;
    ss.num_threads = num_threads;
    
    thread_args_t targs[num_threads];

    // give each thread a unique index, like an MPI worker
    for(i=0; i<num_threads; i++) {
        targs[i].thread_rank = i;
        targs[i].ptr = &ss;
    }

    pthread_t threads[num_threads];

    // create threads
    for(i=0; i<num_threads; i++) 
        pthread_create(&threads[i], NULL, pair_count_wrap, &targs[i]);

    // join threads, sum the array for each thread into one array
    for(i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
        for(j = 0; j <= num_fields; j++) {
            results[j] += ((ss.counters[i])->array)[j];
        }
        free(ss.counters[i]->array);
        free(ss.counters[i]);
    }
    free(ss.counters);
    return results;
}

long long * pair_count_jackknife(vptree_t tree, double * lat, double * lon, int * fields,
                            int length, int num_fields, double radius, int num_threads)
{
    vpdata_t data;
    data.lat = lat;
    data.lon = lon;
    data.fields = fields;
    data.length = length;
    data.num_fields = num_fields;
    // exclude only the field of the querying point
    _include = 0;
    _use_err = 1;
    return pair_count(tree, data, radius, num_threads);
}

long long * pair_count_ftf(vptree_t tree, double * lat, double * lon, int * fields,
                            int length, int num_fields, double radius, int num_threads)
{
    vpdata_t data;
    data.lat = lat;
    data.lon = lon;
    data.fields = fields;
    data.length = length;
    data.num_fields = num_fields;
    // include only the field of the querying point
    _include = 1;
    _use_err = 1;
    return pair_count(tree, data, radius, num_threads);
}

long long * pair_count_noerr(vptree_t tree, double * lat, double * lon,
                            int length, double radius, int num_threads)
{
    vpdata_t data;
    data.lat = lat;
    data.lon = lon;
    data.fields = NULL;
    data.length = length;
    data.num_fields = 1;
    _use_err = 0;

    return pair_count(tree, data, radius, num_threads);
}
