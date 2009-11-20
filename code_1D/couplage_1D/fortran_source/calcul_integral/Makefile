include ../Makefile.config


INCLUDE = -I
FFLAGS  = ${F77_FLAGS} 

FOBJ=${FSRC:.f=.o}

EXEC=calcul


all: $(EXEC)

calcul: calcul_K calcul_K_ANF calcul_K_APF calcul_K_pertes

calcul_K:
	${F77} ${FFLAGS} $(COUPL_LIBS) -o $@ calcul_K.f
	
calcul_K_ANF: 
	${F77} ${FFLAGS} $(COUPL_LIBS) -o $@ calcul_K_ANF.f
	
calcul_K_APF: 
	${F77} ${FFLAGS} $(COUPL_LIBS) -o $@ calcul_K_APF.f
	
calcul_K_pertes: 
	${F77} ${F77_FLAGS} $(COUPL_LIBS) -o $@ calcul_K_pertes.f

clean:
	rm -f *.o
	
mrproper: clean
	rm -f calcul_K_ANF calcul_K_APF calcul_K_pertes	
	
	
