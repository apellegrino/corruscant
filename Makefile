CC= mpicc
BINS= bench libkdtree.so cosmology.so
all: cosmology.so libkdtree.so bench
.PHONY: clean python

# Cosmology lib does not provide much speed benefit so ignore it for now
cosmology.so: cosmology.c
	${CC} -O2 -shared -fPIC $^ -o $@ -lm

python: libkdtree.so

libkdtree.so: kdtree_build.c kdtree_query.c
	${CC} -O2 -shared -fPIC $^ -o $@ -lpthread

bench: bench.o kdtree_build.o kdtree_query.o kdtest.o
	${CC} $^ -o $@ -lpthread

kdtest: kdtest.o kdtree_build.o kdtree_query.o
	${CC} $^ -o $@

%.o: %.c
	${CC} -g -O2 -c -Wall $^

clean:
	rm -f *.o *.so *.pyc ${BINS}
