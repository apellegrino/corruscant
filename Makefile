CC= gcc
BIN_NAMES= bench libkdtree.so
SRC= src
BIN= bin
OBJ= obj

all: mkdirs python ${BIN}/bench

.PHONY: python
python: mkdirs ${BIN}/libkdtree.so

${BIN}/libkdtree.so: ${SRC}/build.c ${SRC}/query.c ${SRC}/kdtest.c
	${CC} -O2 -shared -fPIC $^ -o $@ -lpthread

${BIN}/bench: ${OBJ}/bench.o ${OBJ}/build.o ${OBJ}/query.o ${OBJ}/kdtest.o
	${CC} $^ -o $@ -lpthread -lrt

${OBJ}/kdtest: ${OBJ}/kdtest.o ${OBJ}/build.o ${OBJ}/query.o
	${CC} $^ -o $@

${OBJ}/%.o: ${SRC}/%.c
	${CC} -O2 -c -Wall $< -o $@

.PHONY: mkdirs
mkdirs:
	mkdir -p obj bin

.PHONY: clean
clean:
	rm -f ${OBJ}/*.o ${BIN}/bench ${BIN}/*.so tpcf/*.pyc ${BIN_NAMES}
