CC= gcc
BINS= bench libkdtree.so cosmology.so
all: cosmology.so libkdtree.so bench
.PHONY: clean python

# Cosmology lib does not provide much speed benefit so ignore it for now
cosmology.so: cosmology.c
	${CC} -O2 -shared -fPIC $^ -o $@ -lm

python: libkdtree.so

libkdtree.so: build.c query.c kdtest.c
	${CC} -O2 -shared -fPIC $^ -o $@ -lpthread

bench: bench.o build.o query.o kdtest.o
	${CC} $^ -o $@ -lpthread -lrt

kdtest: kdtest.o build.o query.o
	${CC} $^ -o $@

%.o: %.c
	${CC} -O2 -c -Wall $^

clean:
	rm -f *.o *.so *.pyc ${BINS}
