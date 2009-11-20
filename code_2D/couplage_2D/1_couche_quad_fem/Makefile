
include variables

all: compile link

compile:
	@for i in ${COUPL_DIRS}; do (cd $$i && ${MAKE}); done

link:
	#${LINKER} -o bin/couplage_plasma_2D ${COUPL_LIBS}
	$(F77) $(F77_FLAGS) Main/couplage_plasma_2D.f -o bin/couplage_plasma_2D lib/libsous_prog.a lib/libairy.a lib/libwhittaker.a lib/libquadpack.a lib/libpivot_gauss.a


clean:
	for i in lib ${COUPL_DIRS}; do (cd $$i && ${MAKE} clean ); done
	rm -fr lib/*.a

depend:
	for i in ${GETDP_DIRS};\
          do (cd $$i && ${MAKE} depend "SOLVER=-DHAVE_SPARSKIT -DHAVE_GSL");\
        done

nodepend:
	for i in ${GETDP_DIRS}; do (cd $$i && ${MAKE} nodepend); done
