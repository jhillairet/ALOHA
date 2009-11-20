
include ../variables

LIB     = ../lib/libquad.a
INCLUDE = -I.
CFLAGS  = ${C_FLAGS} ${INCLUDE} 
FFLAGS  = ${F77_FLAGS} ${INCLUDE}

CSRC =  

FSRC = \
quadpack.f\

COBJ  = ${CSRC:.c=.o}
FOBJ  = ${FSRC:.f=.o}

${LIB}: ${COBJ} ${FOBJ}
	${AR} ${LIB} ${COBJ} ${FOBJ}
	${RANLIB} ${LIB}

.f.o:
	${F77} ${FFLAGS} -c $<

clean:
	rm -f *.o
