CC= gcc
CFLAGS= -O3 -funroll-loops -march=native -mtune=native -Wall
#CFLAGS= -g -Wall
SRC= src
BIN= bin
OBJ= obj

all: python benchmark

.PHONY: python
python: mkdirs ${BIN}/libkdtree.so

.PHONY: benchmark
benchmark: mkdirs ${BIN}/kdbench

${BIN}/libkdtree.so: ${SRC}/kdbuild.c ${SRC}/kdquery.c
	${CC} ${CFLAGS} -shared -fPIC $^ -lpthread -o $@

${BIN}/libvptree.so: ${SRC}/vpbuild.c
	${CC} ${CFLAGS} -shared -fPIC $^ -lpthread -o $@

${BIN}/kdbench: ${OBJ}/kdbench.o ${OBJ}/kdbuild.o ${OBJ}/kdquery.o ${OBJ}/kdtest.o
	${CC} $^ -lpthread -lrt -o $@

${BIN}/kdbench_ang: ${OBJ}/kdbench_ang.o ${OBJ}/kdbuild.o ${OBJ}/kdquery.o ${OBJ}/kdtest.o
	${CC} $^ -lpthread -lrt -lm -o $@

#${BIN}/vpbench: ${OBJ}/vpbench.o ${OBJ}/vpbuild.o
#	${CC} $^ -lrt -lm -o $@

${OBJ}/kdtest: ${OBJ}/kdtest.o ${OBJ}/kdbuild.o ${OBJ}/kdquery.o
	${CC} $^ -o $@

${OBJ}/%.o: ${SRC}/%.c
	${CC} -c ${CFLAGS} $< -o $@

.PHONY: mkdirs
mkdirs:
	mkdir -p obj bin

.PHONY: clean
clean:
	rm -f ./${OBJ}/*.o ./${BIN}/* *.pyc
