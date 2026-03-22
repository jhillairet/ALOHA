include ../Makefile.inc

SRC = $(wildcard *.f90)
OBJ = $(SRC:.f90=.o)
MOD = $(SRC:.f90=.mod)

all: $(OBJ) $(LIBQUADPACK)

%.o: %.f90
	$(FC) -O3 $(FCFLAGS) -c $< -o $(OBJDIR)/$@ -I$(MODDIR)

$(LIBQUADPACK):
	ar cr $(OBJDIR)/$(LIBQUADPACK) $(addprefix $(OBJDIR)/,$(OBJ))
	ranlib $(OBJDIR)/$(LIBQUADPACK)

clean:
	$(RM) $(addprefix $(OBJDIR)/,$(OBJ))
	$(RM) $(addprefix $(MODDIR)/,$(MOD))
	$(RM) $(addprefix $(OBJDIR)/,$(LIBQUADPACK))

# dependencies
quadpack.o: quadpack_constants.o
