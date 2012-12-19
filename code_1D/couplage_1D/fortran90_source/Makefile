include Makefile.inc

all: $(addprefix $(OBJDIR)/,$(LIBS))
	$(FC) $(FCFLAGS) -o coupl_plasma_version3 coupl_plasma_version3.f90 $(LDFLAGS)
	$(FC) $(FCFLAGS) -o coupl_plasma_version6 coupl_plasma_version6.f90 $(LDFLAGS)

$(OBJDIR)/$(LIBBLAS):
	@(cd $(LIBBLAS_DIR) && $(MAKE))

$(OBJDIR)/$(LIBLAPACK):
	@(cd $(LIBLAPACK_DIR) && $(MAKE))

$(OBJDIR)/$(LIBBLAS):
	@(cd $(LIBBLAS_DIR) && $(MAKE))

$(OBJDIR)/$(LIBLAPACK):
	@(cd $(LIBLAPACK_DIR) && $(MAKE))

$(OBJDIR)/$(LIBLAPACK95):
	@(cd $(LIBLAPACK95_DIR) && $(MAKE) single_double_complex_dcomplex)

$(OBJDIR)/$(LIBQUADPACK):
	@(cd $(LIBQUADPACK_DIR) && $(MAKE))

$(OBJDIR)/$(LIBAIRY):
	@(cd $(LIBAIRY_DIR) && $(MAKE))

$(OBJDIR)/$(LIBALOHA): $(OBJDIR)/$(LIBAIRY) $(OBJDIR)/$(LIBQUADPACK) $(OBJDIR)/$(LIBLAPACK95)
	@(cd $(LIBALOHA_DIR) && $(MAKE))

clean:
	@(cd $(LIBBLAS_DIR) && $(MAKE) clean)
	@(cd $(LIBLAPACK_DIR) && $(MAKE) clean)
	@(cd $(LIBLAPACK95_DIR) && $(MAKE) clean)
	@(cd $(LIBQUADPACK_DIR) && $(MAKE) clean)
	@(cd $(LIBAIRY_DIR) && $(MAKE) clean)
	@(cd $(LIBALOHA_DIR) && $(MAKE) clean)
	$(RM) $(LIBLAPACK95_DIR)/lapack95_modules/*.mod 
	$(OBJDIR)/*.o | xargs $(RM)
	$(RM) $(OBJDIR)/*.a
	$(RM) $(OBJDIR)/*.o
	$(RM) $(MODDIR)/*.mod
	$(RM) coupl_plasma_version3 coupl_plasma_version6
