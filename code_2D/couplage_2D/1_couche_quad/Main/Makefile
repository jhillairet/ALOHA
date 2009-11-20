
include ../variables

LIB     = ../lib/libMain_bis.a
INCLUDE = -I.
CFLAGS  = ${C_FLAGS} ${INCLUDE} 
FFLAGS  = ${F77_FLAGS} ${INCLUDE}

CSRC =  

FSRC =  \
couplage_plasma_2D.f

COBJ  = ${CSRC:.c=.o}
FOBJ  = ${FSRC:.f=.o}

${LIB}: ${COBJ} ${FOBJ}
	${AR} ${LIB} ${COBJ} ${FOBJ}
	${RANLIB} ${LIB}

.f.o:
	${F77} ${FFLAGS} -c $<

clean:
	rm -f *.o
