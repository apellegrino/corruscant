CFLAGS= -Wall -O3 -funroll-loops -march=native -mtune=native
SRC= src
BIN= bin
OBJ= obj

all: python benchmark

.PHONY: python
python: mkdirs ${BIN}/libkdtree.so ${BIN}/libcoords.so

.PHONY: benchmark
benchmark: mkdirs ${BIN}/kdbench

${BIN}/libkdtree.so: ${SRC}/kdbuild.c ${SRC}/kdquery.c
	${CC} ${CFLAGS} -shared -fPIC $^ -lpthread -o $@

${BIN}/libcoords.so: ${SRC}/coordmath.c
	${CC} ${CFLAGS} -shared -fPIC $^ -o $@

${BIN}/libvptree.so: ${SRC}/vpbuild.c
	${CC} ${CFLAGS} -shared -fPIC $^ -lpthread -o $@

${BIN}/kdbench: ${OBJ}/kdbench.o ${OBJ}/kdbuild.o ${OBJ}/kdquery.o ${OBJ}/kdtest.o
	${CC} $^ -lpthread -lrt -o $@

${BIN}/kdbench_ang: ${OBJ}/kdbench_ang.o ${OBJ}/kdbuild.o ${OBJ}/kdquery.o ${OBJ}/kdtest.o
	${CC} $^ -lpthread -lrt -lm -o $@

${OBJ}/kdtest: ${OBJ}/kdtest.o ${OBJ}/kdbuild.o ${OBJ}/kdquery.o
	${CC} $^ -o $@

${OBJ}/%.o: ${SRC}/%.c
	${CC} -c ${CFLAGS} $< -o $@

.PHONY: mkdirs
mkdirs:
	mkdir -p obj bin

.PHONY: clean
clean:
	rm -f ./${OBJ}/*.o ./${BIN}/*
	find . -name '*.pyc' -delete
