#include "kdtree.h"
#include <stdio.h>

static node_t * tree_data = NULL;

int nodesize(void) {
    return sizeof(node_t);
}

static inline node_t * left_child(node_t * p)
{
    int i = p - tree_data;
    return tree_data + 2 * i;
}

static inline node_t * right_child(node_t * p)
{
    int i = p - tree_data;
    return tree_data + 2 * i + 1;
}

void verify(node_t *root, enum dim d) {
	d=d%3;
    char *errmsg = "node %p should not be %s child of %p (%lf, %lf) \n";

	if (root->flags & HAS_LCHILD) {
		switch(d) {
        case X:
            if(root->x < left_child(root)->x) {
                fprintf(stderr,errmsg,
                (void*)left_child(root),"left",(void*)root);
            }
            break;
        case Y:
            if(root->y < left_child(root)->y) {
                fprintf(stderr,errmsg,
                (void*)left_child(root),"left",(void*)root);
            }
            break;
        case Z:
            if(root->z < left_child(root)->z) {
                fprintf(stderr,errmsg,
                (void*)left_child(root),"left",(void*)root);
            }
            break;
		}
		verify(left_child(root),d+1);
	}
	if (root->flags & HAS_RCHILD) {
		switch(d) {
        case X:
            if(root->x > right_child(root)->x) {
                fprintf(stderr,errmsg,
                (void*)right_child(root),"right",(void*)root);
            }
            break;
        case Y:
            if(root->y > right_child(root)->y) {
                fprintf(stderr,errmsg,
                (void*)right_child(root),"right",(void*)root);
            }
            break;
        case Z:
            if(root->z > right_child(root)->z) {
                fprintf(stderr,errmsg,
                (void*)right_child(root),"right",(void*)root);
            }
            break;
		}
		verify(right_child(root),d+1);
    }
}

void verify_main(node_t *root, enum dim d)
{
    tree_data = root - 1;
    verify(root, d);
}

int count(node_t *p) {
    if(!(p->flags & HAS_RCHILD)) return (p->flags & HAS_LCHILD) + 1;
    return (count(left_child(p)) + count(right_child(p)) + 1);
}
