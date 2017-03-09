CC= mpicc

all: libkdtree
.PHONY: cosmology libkdtree threaded clean

# Cosmology lib does not provide much speed benefit so ignore it for now
cosmology:
	${CC} -shared -lm -o cosmology.so -fPIC cosmology.c

libkdtree:
	${CC} -O2 -shared -o libkdtree.so -fPIC kdtree.c

threaded:
	${CC} -O2 -shared -lpthread -o libkdtree.so -fPIC kdtree_thread.c

bench: bench.o kdtree_thread.o
	${CC} bench.o kdtree_thread.o -lpthread -o bench

kdtree_thread.o: kdtree_thread.c kdtree.h
	${CC} -g -c kdtree_thread.c -O2

bench.o: bench.c kdtree.h
	${CC} -g -c bench.c -O2

#%.o: %.c
#	${CC} -c *.c -O2

clean:
	rm -f *.o *.so *.pyc bench
