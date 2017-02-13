CC= mpicc

all: libkdtree
.PHONY: clean

# Cosmology lib does not provide much speed benefit so ignore it for now

cosmology:
	${CC} -shared -lm -o cosmology.so -fPIC cosmology.c

libkdtree:
	${CC} -O2 -shared -o libkdtree.so -fPIC kdtree.c

clean:
	rm *.so *.pyc
