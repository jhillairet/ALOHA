include ../Makefile.config

# One use of the wildcard function is to get a list of all the source files in a directory, like this:
SRC= $(wildcard *.f)
# la ligne precedente ne fonctionne pas avec la version de make disponible sur deneb
SRC=	S_grill_version1.f \
	S_grill_version2.f \
	S_grill_version3.f \
	S_grill_version4.f \
	S_grill_version5.f \
	S_grill_version6.f \
	S_grill_version7.f
# There is a facility to substitute portions of a predefined variable. The syntax of this is as follows:
#    $(VAR:pattern=replacement)
# should be equivalent to
#    $(patsubst pattern,replacement,$(VAR))    
OBJ=$(SRC:.f=.o) 

# All of these is equivalent to :
# OBJ=$(patsubst %.f,%.o,$(wildcard *.f))
 	  
all: num_int.o $(OBJ)

# actual style (linux)
#.o: .f num_int.o
# old style (unix:deneb)
%.o: %.f num_int.o 
	$(FORTRAN) $(FORTRAN_FLAGS) $(LIBS) -o $@ -c $< 

num_int.o: num_int.f
	$(FORTRAN) $(FORTRAN_FLAGS) $(LIBS) -c $<

clean:
	rm -f *.o
