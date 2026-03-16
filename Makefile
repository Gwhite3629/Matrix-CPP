CC=g++
LINK=g++

TESTDIR:=tests
OBJDIR:=objs
SRCDIR:=.

CONVERG=${TESTDIR}/convergencetest
DECOMP=${TESTDIR}/decomptest
DETER=${TESTDIR}/determinanttest
MATRIX=${TESTDIR}/mattest
GIVENS=${TESTDIR}/givenstest
SPEED=${TESTDIR}/speedtest

LIBOBJS=${OBJDIR}/matrix.o ${OBJDIR}/vector.o
CONVERGOBJS=${OBJDIR}/convergencetest.o
DECOMPOBJS=${OBJDIR}/decomp.o
DETEROBJS=${OBJDIR}/determinanttest.o
MATRIXOBJS=${OBJDIR}/mattest.o
GIVENSOBJS=${OBJDIR}/givenstest.o
SPEEDOBJS=${OBJDIR}/speedtest.o

CFLAGS=-fPIC -Werror -Wextra -pedantic
DFLAGS=-fPIC -ggdb3 -Og -Werror -Wextra -pedantic -Wcast-align -Wcast-qual -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-include-dirs -Wredundant-decls -Wshadow -Wstrict-overflow=5 -Wundef -Wno-unused -Wno-variadic-macros -Wno-parentheses -fdiagnostics-show-option
BFLAGS=-O3
LFLAGS=-ggdb3 -Og
OFLAGS=-fPIC -O3 -rdynamic

all:${CONVERG} ${DECOMP} ${DETER} ${MATRIX} ${GIVENS} ${SPEED} ${OBJDIR}

all : CFLAGS += ${DFLAGS}

${CONVERG}: ${CONVERGOBJS} ${LIBOBJS}
	${CC}  ${CFLAGS} -o $@ $^

${DECOMP}: ${DECOMPOBJS} ${LIBOBJS}
	${CC}  ${CFLAGS} -o $@ $^

${DETER}: ${DETEROBJS} ${LIBOBJS}
	${CC}  ${CFLAGS} -o $@ $^

${MATRIX}: ${MATRIXOBJS} ${LIBOBJS}
	${CC}  ${CFLAGS} -o $@ $^

${GIVENS}: ${GIVENSOBJS} ${LIBOBJS}
	${CC}  ${CFLAGS} -o $@ $^

${SPEED}: ${SPEEDOBJS} ${LIBOBJS}
	${CC}  ${CFLSGS} -o $@ $^

${OBJDIR}/%.o: ${SRCDIR}/%.cpp
	${CC}  ${CFLAGS} -c -o $@ $<

.PHONY : clean

clean:
	rm -f ${TESTDIR}/* core*
	rm -f ${OBJDIR}/* core*