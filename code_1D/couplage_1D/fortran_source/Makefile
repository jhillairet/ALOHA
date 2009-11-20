include Makefile.config

# le code suivant est assez 'rustique', car les possibilites de make sont limitees avec la version disponible sur deneb.
all: LIB \
    coupl_plasma_version3 \
    coupl_plasma_version6

LIB:
	for dir in $(SUBFUNS_DIR) ; do ( cd $$dir ; $(MAKE) all ) ; done
#	#@(cd $(SUBFUNS_DIR) && $(MAKE))

coupl_plasma_version1: LIB
	$(FORTRAN) $(FORTRAN_FLAGS) $(LIBS) -o coupl_plasma_version1 coupl_plasma_version1.f S_grill/S_grill_version1.o 

coupl_plasma_version2: LIB
	$(FORTRAN) $(FORTRAN_FLAGS) $(LIBS) -o coupl_plasma_version2 coupl_plasma_version2.f S_grill/S_grill_version2.o

coupl_plasma_version3: coupl_plasma_version3.f S_grill/S_grill_version3.o
	$(FORTRAN) $(FORTRAN_FLAGS) $(LIBS) -o $@ coupl_plasma_version3.f S_grill/S_grill_version3.o

coupl_plasma_version4: LIB
	$(FORTRAN) $(FORTRAN_FLAGS) $(LIBS) -o coupl_plasma_version4 coupl_plasma_version4.f S_grill/S_grill_version4.o

coupl_plasma_version5: LIB
	$(FORTRAN) $(FORTRAN_FLAGS) $(LIBS) -o coupl_plasma_version5 coupl_plasma_version5.f S_grill/S_grill_version5.o

coupl_plasma_version6: LIB
	$(FORTRAN) $(FORTRAN_FLAGS) $(LIBS) -o coupl_plasma_version6 coupl_plasma_version6.f S_grill/S_grill_version6.o S_grill/num_int.o
#
#coupl_plasma_version7: LIB
#	$(FORTRAN) $(FORTRAN_FLAGS) $(LIBS) -o coupl_plasma_version7 coupl_plasma_version7.f S_grill/S_grill_version7.o

copy:
	@echo "copy binaries into parent directory #? (yes[y]/no[n]) :"
	#read CHOICE
	#if test "${CHOICE}" = "y" || test "${CHOICE}" = "yes" ; then
	cp coupl_plasma_version[1-7] ..
	#fi
clean:
	@(cd $(SUBFUNS_DIR) && $(MAKE) clean)
	rm -fr $(EXEC) *.o *.*~
