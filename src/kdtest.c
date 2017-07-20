#include <stdio.h>

#include "kdtree.h"

int nodesize(void) {
    return sizeof(node_t);
}


static inline int left_child(int p)
{
    return 2*p;
}

static inline int right_child(int p)
{
    return 2*p+1;
}

static inline unsigned long long convert(double x)
{
    /* could do return *((long long *) &x) here, but the following does not
     * raise a strict-aliasing warning from gcc */ 
    union u {
        double d;
        long long l;
    };
    return ((union u *) &x)->d;
}

unsigned long long node_hash(node_t * n)
{
    int i;
    unsigned long long sum = 0;
    for(i=0; i<NDIM; i++) {
        sum += convert(n->data.value[i]);
    }
    return sum;
}

unsigned long long _hash(node_t * data, int i)
{
    node_t n = *(data+i);
    unsigned long long sum = node_hash(data+i);

    if(n.has_lchild) {
        if(n.has_rchild) {
            sum += _hash(data,right_child(i));
        }
        sum += 2 * _hash(data,left_child(i));
    }

    return sum;
}

unsigned long long hash(kdtree_t t)
{
    return _hash(t.node_data,1);
}

static void _verify(node_t * data, int index) {
    node_t *parent, *lc, *rc;
    char *errmsg = "node %d should not be %s child of %d in dimension %d\n";

    parent = data + index;
    int dim = parent->dim;

	if (parent->has_lchild) {
        lc = data + left_child(index);
        if( (parent->data).value[dim] < (lc->data).value[dim] ) {
            fprintf(stderr,errmsg,left_child(index),"left",index,dim);
        }
		_verify(data,left_child(index));
	}
	if (parent->has_rchild) {
        rc = data + right_child(index);
        if( (parent->data).value[dim] > (rc->data).value[dim] ) {
            fprintf(stderr,errmsg,right_child(index),"right",index,dim);
        }
		_verify(data,right_child(index));
	}
}

void verify_tree(kdtree_t t)
{
    _verify(t.node_data,1);
}

static int _count(kdtree_t t, int index)
{
    node_t * p = t.node_data + index;
    if(!(p->has_rchild)) return (p->has_lchild) + 1;
    return (_count(t, left_child(index)) + _count(t, right_child(index)) + 1);
}

int count_tree(kdtree_t t)
{
    return _count(t,1);
}
