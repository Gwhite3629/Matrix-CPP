CC=g++
LINK=g++
OBJS=../vector.o ../matrix.o

DET-ECHELON=det-echelon
DET-LAPLACE=det-laplace
MATMUL=matmul
QR=qr
INVERT=invert
CHOLESKY=cholesky
LU=lu

DET-ECHELON-OBJ=det-echelon.o ${OBJS}
DET-LAPLACE-OBJ=det-laplace.o ${OBJS}
MATMUL-OBJ=matmul.o ${OBJS}
QR-OBJ=qr.o ${OBJS}
INVERT-OBJ=invert.o ${OBJS}
CHOLESKY-OBJ=cholesky.o ${OBJS}
LU-OBJ=lu.o ${OBJS}

LIBS=-lm
CFLAGS=-O3 -Wall -Wextra
LFLAGS=-O3

all: ${DET-ECHELON} ${DET-LAPLACE} ${MATMUL} ${QR} ${INVERT} ${CHOLESKY} ${LU}

${DET-ECHELON}: ${DET-ECHELON-OBJ}
	${CC} ${LFLAGS} -o $@ $^

${DET-LAPLACE}: ${DET-LAPLACE-OBJ}
	${CC} ${LFLAGS} -o $@ $^

${MATMUL}: ${MATMUL-OBJ}
	${CC} ${LFLAGS} -o $@ $^

${QR}: ${QR-OBJ}
	${CC} ${LFLAGS} -o $@ $^

${INVERT}: ${INVERT-OBJ}
	${CC} ${LFLAGS} -o $@ $^

${CHOLESKY}: ${CHOLESKY-OBJ}
	${CC} ${LFLAGS} -o $@ $^

${LU}: ${LU-OBJ}
	${CC} ${LFLAGS} -o $@ $^

det-echelon.o: det-echelon.cpp
det-laplace.o: det-laplace.cpp
matmul.o: matmul.cpp
qr.o: qr.cpp
invert.o: invert.cpp
cholesky.o: cholesky.cpp
lu.o: lu.cpp

vector.o: ../vector.cpp ../vector.h
matrix.o: ../matrix.cpp ../matrix.h

.PHONY : clean

clean:
	rm -f ${DET-ECHELON} core*
	rm -f ${DET-LAPLACE} core*
	rm -f ${MATMUL} core*
	rm -f ${QR} core*
	rm -f ${INVERT} core*
	rm -f *.o core*
	rm -f ../*.o core*