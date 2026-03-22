# HELP for compiling LAPACK95 with GORTRAN :
# http://gcc.gnu.org/wiki/GfortranBuild#head-093aadd510514f31d177acf4d75acb21c83b814c


#
#  -- LAPACK95 interface driver routine (version 2.0) --
#     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
#     August 5, 2000
#

# g95
#FC	 = g95 -ffree-form
#FC1      = g95 -ffixed-form

# gfortran
FC       = gfortran -ffree-form
FC1	 = gfortran -ffixed-form

# -dcfuns  Enable recognition of non-standard double
#          precision  complex intrinsic functions
# -dusty   Allows the compilation and execution of "legacy"
#          software by downgrading the category  of  common
#          errors  found  in  such software from "Error" to
# -ieee=full enables all IEEE arithmetic facilities
#          including non-stop arithmetic.
#OPTS0    = #-ieee=full -march=opteron # -u -V -dcfuns -dusty -ieee=full # g95 options

# gfortran specific options :
OPTS0    = -mtune=generic -pipe -ftree-vectorize
# -pipe : Use pipes rather than temporary files for communication between the various stages of compilation. This fails to work on some systems where the assembler is unable to read from a pipe; but the GNU assembler has no trouble.
# -ftree-vectorize : Perform loop vectorization on trees.
# -mtune=generic
# if you know the CPU on which your code will run, then you should use the corresponding mtune switch found in :
# http://gcc.gnu.org/onlinedocs/gcc-4.3.0/gcc/i386-and-x86_002d64-Options.html#i386-and-x86_002d64-Options

MODLIB   = -I./../lapack95_modules
OPTS1    = -c -O3 $(OPTS0)
OPTS3    = $(OPTS1) $(MODLIB)
OPTL     = -o 
OPTLIB   = 

LAPACK_PATH = ../obj

LAPACK95 = ../lapack95.a
LAPACK77 = $(LAPACK_PATH)/liblapack.a
TMG77    = $(LAPACK_PATH)/tmglib.a
BLAS     = $(LAPACK_PATH)/libblas.a


LIBS     = $(LAPACK95) $(TMG77) $(LAPACK77) $(BLAS)
SUF      = f90

XX = 'rm' -f $@; \
        'rm' -f $@.res; \
	$(FC) $(OPTS0) -o $@ $(MODLIB) $@.$(SUF) $(OPTLIB) $(LIBS); \
        $@ < $@.dat > $@.res; \
        'rm' -f $@

YY = $(FC) $(OPTS0) -o $@ $(MODLIB) $@.$(SUF) $(OPTLIB) $(LIBS)

.SUFFIXES: .f90 .f .o

.$(SUF).o: 
	$(FC) $(OPTS3) $<

.f.o:
	$(FC1) $(OPTS3) $<
