SHELL           = /bin/bash
FC              = gfortran
MAKEDEPF90      = /usr/local/bin/makedepf90
NETCDF_INCLUDES = -I/usr/lib64/gfortran/modules/
NETCDF_LIBS     = -L/usr/lib64/ -lnetcdff
MKLROOT         = /home/gblumberg/intel/oneapi/mkl/latest/
FCFLAGS         = -g -pg -fbacktrace -fcheck=all -Wall $(NETCDF_INCLUDES) -m64 -I$(MKLROOT)/include
#FCFLAGS         = -Ofast -pg -fbacktrace -fcheck=all -Wall $(NETCDF_INCLUDES) -m64 -I$(MKLROOT)/include
LDFLAGS         = $(NETCDF_LIBS) -L$(MKLROOT)/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm

PROG            = ozo
SRC             = ozo.f90 mod_const.f90 mod_poisson_DFT.f90 mod_zo.f90 \
	          mod_common_subrs.f90 mod_time_step_loop.f90 \
	          mod_wrf_file.f90 mkl_poisson.f90 mkl_dfti.f90 \
	          mod_poisson_green.f90 mod_poisson_interp.f90 \
		  mod_omega.f90
OBJS            = $(SRC:.f90=.o) 

vpath %.f90 $(MKLROOT)/include src

.PHONY : all test clean realclean

all : $(PROG)

$(PROG) : $(OBJS)
	$(FC) -o $@ $(FCFLAGS) $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

test : $(PROG)
	@for test in $^; do cd test; LD_LIBRARY_PATH=$(MKLROOT)/lib/ ../$$test < namelist; done

clean :
	rm  -f *.o *.mod $(PROG)

# Create dependencies automatically
deps.mk : src/*.f90 mkl_poisson.f90 mkl_dfti.f90
	[ -x $(MAKEDEPF90) ] && $(MAKEDEPF90) -b '' $^ > $@

include deps.mk
