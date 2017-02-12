all: driver

# Cosmology lib does not provide much speed benefit so ignore it for now

cosmology:
	gcc -shared -lm -o cosmology.so -fPIC cosmology.c

driver:
	gcc -O2 -shared -o driver.so -fPIC driver.c kdtree.c

.PHONY: clean
clean:
	rm *.so *.pyc
