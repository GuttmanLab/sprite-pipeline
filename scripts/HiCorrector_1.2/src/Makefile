# compiler and options for ANSI C
CC=cc
CFLAGS= -Wall -O3
LIBS= -lm

# mpi compiler and options for ANSI C and mpi library
CCMPI=/usr/usc/openmpi/default/bin/mpicc
CFLAGSMPI=
LIBSMPI= -L/usr/usc/openmpi/1.6.4/lib
INCLUDE= 

# getopt functions library
getoptDIR=./my_getopt-1.5

# programs to be compiled and built
OBJ_GETOPT= my_getopt.o
OBJECTS1= ic.o matutils.o
OBJECTS2= ic_mes.o matutils.o
OBJECTS3= ic_mep.o matutils.o ${OBJ_GETOPT}
OBJECTS4= split_data.o matutils.o
OBJECTS5= split_data_parallel.o matutils.o
OBJECTS6= export_norm_data.o matutils.o
PROG1=ic
PROG2=ic_mes
PROG3=ic_mep
PROG4=split_data
PROG5=split_data_parallel
PROG6=export_norm_data

all: ${PROG1} ${PROG2} ${PROG3} ${PROG4} ${PROG5} ${PROG6}

${PROG1}: ${OBJECTS1}
	${CC} ${CFLAGS} ${OBJECTS1} ${LIBS} -o ${PROG1}
	cp ${PROG1} ../bin

${PROG2}: ${OBJECTS2}
	${CC} ${CFLAGS} ${OBJECTS2} ${LIBS} -o ${PROG2}
	cp ${PROG2} ../bin

${PROG3}: ${OBJECTS3}
	${CCMPI} ${CFLAGSMPI} ${OBJECTS3} ${LIBSMPI} ${LIBS} -o ${PROG3}
	cp ${PROG3} ../bin

${PROG4}: ${OBJECTS4}
	${CC} ${CFLAGS} ${OBJECTS4} ${LIBS} -o ${PROG4}
	cp ${PROG4} ../bin

${PROG5}: ${OBJECTS5}
	${CCMPI} ${CFLAGSMPI} ${OBJECTS5} ${LIBSMPI} ${LIBS} -o ${PROG5}
	cp ${PROG5} ../bin

${PROG6}: ${OBJECTS6}
	${CC} ${CFLAGS} ${OBJECTS6} ${LIBS} -o ${PROG6}
	cp ${PROG6} ../bin

matutils.o: matutils.c matutils.h
	${CC} ${CFLAGS} ${INCLUDE} -c $< -o $@

ic.o: ic.c matutils.h
	${CC} ${CFLAGS} ${INCLUDE} -c $< -o $@

ic_mes2.o: ic_mes2.c matutils.h
	${CC} ${CFLAGS} ${INCLUDE} -c $< -o $@

ic_mep.o: ic_mep.c matutils.h
	${CCMPI} ${CFLAGS} ${INCLUDE} -c $< -o $@

ic_mep2.o: ic_mep2.c matutils.h
	${CCMPI} ${CFLAGS} ${INCLUDE} -c $< -o $@

split_data.o: split_data.c matutils.h
	${CC} ${CFLAGS} ${INCLUDE} -c $< -o $@

split_data_parallel.o: split_data_parallel.c matutils.h
	${CCMPI} ${CFLAGS} ${INCLUDE} -c $< -o $@

export_norm_data.o: export_norm_data.c matutils.h
	${CC} ${CFLAGS} ${INCLUDE} -c $< -o $@

my_getopt.o: $(getoptDIR)/my_getopt.c $(getoptDIR)/my_getopt.h $(getoptDIR)/getopt.h
	cd $(getoptDIR)
	$(CC) $(Flags) -I$(getoptDIR) -c $(getoptDIR)/my_getopt.c -o my_getopt.o

# --- remove binary and executable files
clean:
	rm -f *.o
	rm -f ${PROG1}
	rm -f ${PROG2}
	rm -f ${PROG3}
	rm -f ${PROG4}
	rm -f ${PROG5}
	rm -f ${PROG6}
