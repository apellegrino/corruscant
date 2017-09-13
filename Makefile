CFLAGS= -Wall -Wextra -O3 -march=native -mtune=native
SRC= src
BIN= bin
OBJ= obj

.PHONY: python
python: mkdirs_py ${BIN}/libkdtree.so ${BIN}/libcoords.so

.PHONY: benchmark
benchmark: mkdirs_c ${BIN}/kdbench

.PHONY: all
all: python benchmark

${BIN}/libkdtree.so: ${SRC}/kdbuild.c ${SRC}/kdquery.c
	${CC} ${CFLAGS} -shared -fPIC $^ -lpthread -o $@

${BIN}/libcoords.so: ${SRC}/coordmath.c
	${CC} ${CFLAGS} -shared -fPIC $^ -o $@

${BIN}/libvptree.so: ${SRC}/vpbuild.c
	${CC} ${CFLAGS} -shared -fPIC $^ -lpthread -o $@

${BIN}/kdbench: ${OBJ}/kdbench.o ${OBJ}/kdbuild.o ${OBJ}/kdquery.o ${OBJ}/kdtest.o
	${CC} $^ -lpthread -lrt -o $@

${BIN}/analytic: ${OBJ}/analytic.o ${OBJ}/kdbuild.o ${OBJ}/kdquery.o
	${CC} $^ -lpthread -lrt -o $@

${BIN}/kdbench_ang: ${OBJ}/kdbench_ang.o ${OBJ}/kdbuild.o ${OBJ}/kdquery.o ${OBJ}/kdtest.o
	${CC} $^ -lpthread -lrt -lm -o $@

${OBJ}/kdtest: ${OBJ}/kdtest.o ${OBJ}/kdbuild.o ${OBJ}/kdquery.o
	${CC} $^ -o $@

${OBJ}/%.o: ${SRC}/%.c
	${CC} -c ${CFLAGS} $< -o $@

.PHONY: mkdirs_c
mkdirs_c:
	mkdir -p obj bin

.PHONY: mkdirs_py
mkdirs_py:
	mkdir -p bin

.PHONY: clean
clean:
	rm -f ./${OBJ}/*.o ./${BIN}/*
	find . -name '*.pyc' -delete
