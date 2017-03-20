CC= mpicc
BINS= bench libkdtree.so cosmology.so
all: cosmology.so libkdtree.so bench
.PHONY: cosmology libkdtree threaded clean

# Cosmology lib does not provide much speed benefit so ignore it for now
cosmology.so: cosmology.c
	${CC} -O2 -shared -fPIC $^ -o $@ -lm

libkdtree.so: kdtree.c
	${CC} -O2 -shared -fPIC $^ -o $@ -lpthread

bench: bench.o kdtree.o kdtest.o
	${CC} $^ -o $@ -lpthread

kdtest: kdtest.o kdtree.o
	${CC} $^ -o $@

%.o: %.c
	${CC} -g -O2 -c -Wall $^

clean:
	rm -f *.o *.so *.pyc ${BINS}
