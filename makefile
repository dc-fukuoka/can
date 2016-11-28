SERIAL   = seri
OMP      = omp
DGEMM    = dgemm
CAN      = can
CAN_HYB  = can_hyb
VSL      = mkl_vsl.mod

FC      = ifort
MPIFC   = mpiifort
FCFLAGS = -g -O3 -mavx -mcmodel=medium -shared-intel
LFLAGS  = 
LIBS    = -mkl
OPENMP  = -fopenmp
PREP    = -fpp
SRCS    = $(MKLROOT)/include/mkl_vsl.f90 generate.f seri.f dgemm.f omp.f

.SUFFIXES: .f.o

ALL: $(SERIAL) $(DGEMM) $(OMP) $(CAN) $(CAN_HYB)

$(VSL): $(MKLROOT)/include/mkl_vsl.f90
	ifort -c $(FCFLAGS) $(MKLROOT)/include/mkl_vsl.f90

generate.o: $(VSL)
	$(FC) -c $(FCFLAGS) generate.f
seri.o: seri.f
	$(FC) -c $(FCFLAGS) seri.f
dgemm.o: dgemm.f
	$(FC) -c $(FCFLAGS) dgemm.f
omp.o: omp.f
	$(FC) -c $(OPENMP) $(FCFLAGS) omp.f
can.o: can.f
	$(MPIFC) -c $(FCFLAGS) can.f
can_hyb.o: can_hyb.f
	$(MPIFC) $(PREP) -c $(OPENMP) $(FCFLAGS) can_hyb.f

$(SERIAL): generate.o seri.o
	$(FC) $(LFLAGS) $(LIBS) seri.o generate.o -o $(SERIAL)
$(DGEMM): generate.o dgemm.o
	$(FC) $(LFLAGS) $(LIBS) dgemm.o generate.o -o $(DGEMM)
$(OMP): generate.o omp.o
	$(FC) $(OPENMP) $(LFLAGS) $(LIBS) omp.o generate.o -o $(OMP)
$(CAN): generate.o can.o
	$(MPIFC) $(LFLAGS) $(LIBS) can.o generate.o -o $(CAN)
$(CAN_HYB): generate.o can_hyb.o
	$(MPIFC) $(OPENMP) $(LFLAGS) $(LIBS) can_hyb.o generate.o -o $(CAN_HYB)
clean:
	rm -f *.o *.mod $(SERIAL) $(DGEMM) $(OMP) $(CAN) $(CAN_HYB)
