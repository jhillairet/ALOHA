
include ../variables

LIB     = ../lib/libpivot_gauss.a
INCLUDE = -I.
CFLAGS  = ${C_FLAGS} ${INCLUDE} 
FFLAGS  = ${F77_FLAGS} ${INCLUDE}

CSRC =  

FSRC = \
zgbsv.f \
zgbtf2.f \
zgbtrf.f \
zgbtrs.f \
zlacgv.f \
zlaswp.f \
xerbla.f\
lsame.f\
zswap.f\
zgeru.f\
ztbsv.f\
zgemv.f\
izamax.f\
zscal.f\
zcopy.f\
ztrsm.f\
zgemm.f\
dcabs1.f\
ilaenv.f\
iparmq.f\
ieeeck.f\

COBJ  = ${CSRC:.c=.o}
FOBJ  = ${FSRC:.f=.o}

${LIB}: ${COBJ} ${FOBJ}
	${AR} ${LIB} ${COBJ} ${FOBJ}
	${RANLIB} ${LIB}

.f.o:
	${F77} ${FFLAGS} -c $<

clean:
	rm -f *.o
