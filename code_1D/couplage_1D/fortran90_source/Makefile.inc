####################################################################
#  ALOHA  make include and initialisation file.       #
#
#  Last update :
#  -------------
# 04/03/2017
#                                              
####################################################################

# Path and directories configuration
# TOP path 
# the next line retrieve automaticcaly the relative path to the top directory
# thus avoiding the need to set an absolute path 
ALOHAPATH:=$(dir $(lastword $(MAKEFILE_LIST)))

LIBALOHA_DIR = $(ALOHAPATH)/libaloha
LIBALOHA_HOT_DIR = $(ALOHAPATH)/libaloha_hot

LIBBLAS_DIR=$(ALOHAPATH)/liblapack/BLAS/SRC
LIBLAPACK_DIR=$(ALOHAPATH)/liblapack/SRC

LIBLAPACK95_DIR = $(ALOHAPATH)/liblapack95/SRC

LIBQUADPACK_DIR = $(ALOHAPATH)/libquadpack

LIBAIRY_DIR = $(ALOHAPATH)/libairy/Src/Fortran90/Src/


MODDIR = $(ALOHAPATH)/mod
OBJDIR = $(ALOHAPATH)/obj

VPATH = $(MODDIR):$(OBJDIR)

# fortran compiler 
#FC = gfortran # nashira
#FC = /usr/local/pgi/linux86-64/7.1-6/bin/pgf90 # saturne
FC = g95 #saturne

# libraries
LIBALOHA=libaloha.a
LIBALOHA_HOT=libaloha_hot.a
LIBQUADPACK=libquadpack.a
LIBAIRY=libairy.a
LIBNAG=#libnag.a

LIBBLAS=libblas.a
LIBLAPACK=liblapack.a
LIBLAPACK95=liblapack95.a

# attention : l'ordre est tres important pour gfortran. 
LIBS= $(LIBALOHA) $(LIBLAPACK95) $(LIBLAPACK) $(LIBBLAS) $(LIBNAG) $(LIBQUADPACK) $(LIBAIRY)

# Compiler flags
# -g : debug
# -pg: profile with gprof
# -O : optimisation level [1,2,3] : a proscrire : donne des resultats numeriques errones a cause de libaloha

# optimisation options
OPTS0    = -O3 #g95
# OPTS0    = -mtune=generic -pipe -ftree-vectorize #gfortran

DEBUG = -g -fbounds-check -ftrace=full




# for pgf90 : -module $(MODDIR)
#FCFLAGS = -g -module $(MODDIR)

# for g95
FCFLAGS = $(OPTS0) -fmod=$(MODDIR) -L $(OBJDIR) # g95

# for gfortran
#FCFLAGS = $(OPTS0) -M $(MODDIR) -L $(OBJDIR) # gfortran

# linker flags
LDFLAGS = -I $(MODDIR)/ $(addprefix $(OBJDIR)/,$(LIBS))

# HELP GFORTRAN :
# GCC command options.
# http://gcc.gnu.org/onlinedocs/gcc-4.3.0/gcc/Option-Index.html
# -L : Add directory dir to the list of directories to be searched for -l.
# -llibrary ou -l library : Search the library named library when linking.
#
# GFORTRAN options
# http://gcc.gnu.org/onlinedocs/gcc-4.3.0/gfortran/Invoking-GNU-Fortran.html
# -Idir : Add the directory dir to the head of the list of directories to be searched for header files.
# -Mdir ou -Jdir : This option specifies where to put .mod files for compiled modules.
#                  It is also added to the list of directories to searched by an USE statement.
# g95: -fmod=dir
#

