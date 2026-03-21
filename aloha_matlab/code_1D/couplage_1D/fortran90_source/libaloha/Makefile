include ../Makefile.inc

SRC = $(wildcard *.f90)
OBJ = $(SRC:.f90=.o)
MOD = $(SRC:.f90=.mod)

all: $(OBJ) $(LIBALOHA)

#%.o: %.f90
#	$(FC) $(FCFLAGS) -c $< -o $(OBJDIR)/$@ -I$(MODDIR)

$(LIBALOHA): $(OBJ)
	ar cr $(OBJDIR)/$@ $(addprefix $(OBJDIR)/,$(OBJ))
	ranlib $(OBJDIR)/$@

clean:
	$(RM) $(addprefix $(OBJDIR)/, $(OBJ))
	$(RM) $(addprefix $(MODDIR)/, $(MOD))

# dependencies
aloha_constants.o: aloha_constants.f90
	$(FC) -O3 $(FCFLAGS) -c $< -o $(OBJDIR)/$@ -I$(MODDIR)

aloha_config.o: aloha_config.f90 aloha_constants.o
	$(FC) -O3 $(FCFLAGS) -c $< -o $(OBJDIR)/$@ -I$(MODDIR)

aloha_function.o: aloha_function.f90 aloha_config.o aloha_constants.o
	$(FC) $(FCFLAGS) -c $< -o $(OBJDIR)/$@ -I$(MODDIR)

aloha_scattering_matrix_utils.o: aloha_scattering_matrix_utils.f90 aloha_constants.o
	$(FC) -O3 $(FCFLAGS) -c $< -o $(OBJDIR)/$@ -I$(MODDIR)

aloha_sgrill.o: aloha_sgrill.f90 aloha_config.o aloha_constants.o aloha_function.o aloha_scattering_matrix_utils.o
	$(FC) -O3 $(FCFLAGS) -c $< -o $(OBJDIR)/$@ -I$(MODDIR)

aloha_integration.o: aloha_integration.f90 quadpack.o
	$(FC) -O3 $(FCFLAGS) -c $< -o $(OBJDIR)/$@ -I$(MODDIR)

