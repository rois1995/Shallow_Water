############################################################
#                     MAKEFILE FOR F95                     #
############################################################

# Zero axiom: the file name of a makefile must be mandatorily
#
#              makefile    OR    Makefile

# To start the execution, type in terminal: make

# The names (identifictars of variables) for aliases
# are completely free!   They are usually in CAPITALS.

# FUNDAMENTAL WARNINGS
#
# The actions, typically $(CMP) and $(LNK), must be
# NECESSARILY preceded by a DAMNED! TAB character
#

############################################################
#
# Compiler and linker
#
# The following definitions are for Fortran 95 only,
# but they can extended if another language is used,
# or modified when other compilers are considered.

# LNK alias to be used as $(LNK) to compile and link main files .f90
# CMP alias to be used as $(CMP) to compile the source files .f90
# OPT alias to be used as $(OPT) to select the compilation options

# NAG Compiler
# ============
#
CMP = f95 -c
LNK = f95
OPT = -g -mcmodel=large
# OPT = -gline # Line number is indicated for errors
               # occurring in runtime (for debugging)
# OPT = -O3    # Execution optimization


# Intel Compiler
# ==============
#
# CMP = ifort  -c  -heap-arrays 100
# LNK = ifort
# OPT =
# OPT = -g  -DD  -debug -traceback  -fpe0  -check all
# OPT = -fast
# OPT = -O3


# GNU Fortran Compiler
# ===================
#
# CMP = gfortran  -c
# LNK = gfortran
# OPT = -fno-automatic  -fbounds-check  -ffpe-trap=invalid,zero,overflow  -ggdb3
# OPT = -fno-automatic  -O3

############################################################
#
# Objects: list of all objects *.o

OBJS = gnufor_mio.o LOGO.o conversione.o inversa.o problema_vero.o fluxes.o\
       Boundary_conditions.o #fix_BC.o
# OBJS alias to be used as $(OBJS)

# Alternative way for line feed
#OBJS = Riemann_solvers.o  \
#       polytropic_ideal_gas.o

############################################################
# Executable generation by the linker / the effect of command
#                                     / make starts from here

finale.exe:  SW_completo.o  $(OBJS)
	$(LNK) $(OPT)    SW_completo.o  $(OBJS)  \
              -o finale.exe  # Name chosen for .exe

############################################################
# Objects generation

SW_completo.o:     SW_completo.f90  $(OBJS)
	$(CMP) $(OPT)     SW_completo.f90

gnufor_mio.o:  gnufor_mio.f90
	$(CMP) $(OPT)     gnufor_mio.f90

LOGO.o: LOGO.f90
	$(CMP) $(OPT)  LOGO.f90

conversione.o: conversione.f90
	$(CMP) $(OPT)  conversione.f90

inversa.o: inversa.f90
	$(CMP) $(OPT)  inversa.f90

problema_vero.o: problema_vero.f90  inversa.o
	$(CMP) $(OPT)  problema_vero.f90

Boundary_conditions.o: Boundary_conditions.f90  problema_vero.o conversione.o
	$(CMP) $(OPT)  Boundary_conditions.f90

#fix_BC.o: fix_BC.f90  conversione.o
#	$(CMP) $(OPT)  fix_BC.f90

fluxes.o: fluxes.f90  problema_vero.o conversione.o Boundary_conditions.o
	$(CMP) $(OPT)  fluxes.f90



#plot.o: plot.f90 checking.o compl_Newton.o gnufor_mio.o
#	$(CMP) $(OPT)     plot.f90

#
# il carattere TAB c'e' ma non appare a causa di #
# provare a fare return appena dopo # per credere
# il TAB (con il suo fottutissimo doppio spazio) resuscita


############################################################
# Cleaning command to remove all objects *.o, files *.mod and
# the executables *.exe

clean:
	@echo cleaning objects, modules, executables
	rm  *.o  *.mod  *.exe

cleanimages:
	@echo cleaning images files
	rm  video/plot_video/*

cleanvideo:
	@echo cleaning video files
	rm  plot_video/*

# The command "make clean" deletes all the indicated files

############################################################
