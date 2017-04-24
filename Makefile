CC= gcc
CFLAGS= -O2 -Wall
#CFLAGS= -g -Wall
BIN_NAMES= bench libkdtree.so
SRC= src
BIN= bin
OBJ= obj

all: python benchmark

.PHONY: python
python: mkdirs ${BIN}/libkdtree.so

.PHONY: benchmark
benchmark: mkdirs ${BIN}/bench

${BIN}/libkdtree.so: ${SRC}/build.c ${SRC}/query.c ${SRC}/kdtest.c
	${CC} ${CFLAGS} -shared -fPIC $^ -o $@ -lpthread

${BIN}/bench: ${OBJ}/bench.o ${OBJ}/build.o ${OBJ}/query.o ${OBJ}/kdtest.o
	${CC} $^ -o $@ -lpthread -lrt

${OBJ}/kdtest: ${OBJ}/kdtest.o ${OBJ}/build.o ${OBJ}/query.o
	${CC} $^ -o $@

${OBJ}/%.o: ${SRC}/%.c
	${CC} -c ${CFLAGS} $< -o $@

.PHONY: mkdirs
mkdirs:
	mkdir -p obj bin

.PHONY: clean
clean:
	rm -f ${OBJ}/*.o ${BIN}/bench ${BIN}/*.so *.pyc ${BIN_NAMES}
