#include "kdtree.h"
#include <stdio.h>

int nodesize(void) {
    return sizeof(node_t);
}

static void verify(node_t * data, int index, enum dim d) {
	d=d%3;

    node_t *parent, *lc, *rc;
    char *errmsg = "node %p should not be %s child of %p (%lf, %lf) \n";
    int comp;

    parent = data + index;

	if (parent->flags & HAS_LCHILD) {
        lc = data + left_child(index);
		switch(d) {
        case X:
            comp = parent->x < lc->x;
            break;
        case Y:
            comp = parent->y < lc->y;
            break;
        case Z:
            comp = parent->z < lc->z;
            break;
		}
        if(comp) {
            fprintf(stderr,errmsg,left_child(index),"left",index);
        }
		verify(data,left_child(index),d+1);
	}
	if (parent->flags & HAS_RCHILD) {
        rc = data + right_child(index);
		switch(d) {
        case X:
            comp = parent->x > rc->x;
            break;
        case Y:
            comp = parent->y > rc->y;
            break;
        case Z:
            comp = parent->z > rc->z;
            break;
		}
        if(comp) {
            fprintf(stderr,errmsg,right_child(index),"right",index);
        }
		verify(data,right_child(index),d+1);
	}
}

void verify_main(kdtree_t t, enum dim d)
{
    verify(t.node_data,1,d);
}

static int count(kdtree_t t, int index)
{
    node_t * p = t.node_data + index;
    if(!(p->flags & HAS_RCHILD)) return (p->flags & HAS_LCHILD) + 1;
    return (count(t, left_child(index)) + count(t, right_child(index)) + 1);
}

int count_main(kdtree_t t)
{
    return count(t,1);
}
