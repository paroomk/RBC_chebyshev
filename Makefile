FC=gfortran
FFLAGS= -c -Wall -Wno-maybe-uninitialized -O3

LIBFLAGS  = -L/usr/lib64
LIBFLAGS2 = -L/Users/dsondak/fftw/install/lib
LDFLAGS   = -I/Users/dsondak/fftw/install/include

MAIN = TC

OBJECTS = types.o fftw.o global.o write_pack.o makeICs.o chebutils.o time_integrators.o statutils.o $(MAIN).o
PROGRAMS = $(MAIN).exe

all: $(PROGRAMS)

$(PROGRAMS) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $(PROGRAMS) $(OBJECTS) $(LIBFLAGS1) -llapack -lblas $(LIBFLAGS2) -lfftw3 -lm

types.o : types.f03
	$(FC) $(FFLAGS) types.f03

fftw.o : fftw.f03
	$(FC) $(FFLAGS) fftw.f03

global.o : global.f03
	$(FC) $(FFLAGS) global.f03

write_pack.o : write_pack.f03
	$(FC) $(FFLAGS) write_pack.f03

makeICs.o : makeICs.f03
	$(FC) $(FFLAGS) makeICs.f03

chebutils.o : chebutils.f03
	$(FC) $(FFLAGS) chebutils.f03

time_integrators.o : time_integrators.f03
	$(FC) $(FFLAGS) time_integrators.f03

statutils.o : statutils.f03
	$(FC) $(FFLAGS) statutils.f03


$(MAIN).o : $(MAIN).f03
	$(FC) $(FFLAGS) $(MAIN).f03

clean :
	rm -rf *.mod $(OBJECTS) $(PROGRAMS) 

cleanall :
	rm -rf *.mod *.txt $(OBJECTS) $(PROGRAMS) 
