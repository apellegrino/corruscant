CC= mpicc

all: driver

# Cosmology lib does not provide much speed benefit so ignore it for now

cosmology:
	${CC} -shared -lm -o cosmology.so -fPIC cosmology.c

driver:
	${CC} -O2 -shared -o driver.so -fPIC driver.c kdtree.c

.PHONY: clean
clean:
	rm *.so *.pyc
