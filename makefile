CREATE_INPUT = create_input
SERIAL       = seri
OMP          = omp
DGEMM        = dgemm
CAN          = can
CAN_HYB      = can_hyb
VSL          = mkl_vsl.mod
CHECK        = check

FC      = ifort
MPIFC   = mpiifort
FCFLAGS = -g -O3 -march=core-avx2 -mcmodel=medium -shared-intel
LFLAGS  = 
LIBS    = -mkl
OPENMP  = -fopenmp
PREP    = -fpp
SRCS    = $(MKLROOT)/include/mkl_vsl.f90 generate.f seri.f dgemm.f omp.f

.SUFFIXES: .f.o

ALL: $(CREATE_INPUT) $(SERIAL) $(DGEMM) $(OMP) $(CAN) $(CAN_HYB) $(CHECK)

$(VSL): $(MKLROOT)/include/mkl_vsl.f90
	ifort -c $(FCFLAGS) $(MKLROOT)/include/mkl_vsl.f90

create_input.o: create_input.f
	$(FC) -c $(FCFLAGS) $<
generate.o: generate.f $(VSL)
	$(FC) -c $(FCFLAGS) $<
seri.o: seri.f
	$(FC) -c $(FCFLAGS) $<
dgemm.o: dgemm.f
	$(FC) -c $(FCFLAGS) $<
omp.o: omp.f
	$(FC) -c $(OPENMP) $(FCFLAGS) $<
can.o: can.f
	$(MPIFC) -c $(FCFLAGS) $<
can_hyb.o: can_hyb.f
	$(MPIFC) $(PREP) -c $(OPENMP) $(FCFLAGS) $<
check.o: check.f
	$(FC) -c $(OPENMP) $(FCFLAGS) $<

$(CREATE_INPUT): generate.o create_input.o
	$(FC) $(LFLAGS) $(LIBS) $^ -o $@
$(SERIAL): generate.o seri.o
	$(FC) $(LFLAGS) $(LIBS) $^ -o $@
$(DGEMM): generate.o dgemm.o
	$(FC) $(LFLAGS) $(LIBS) $^ -o $@
$(OMP): generate.o omp.o
	$(FC) $(OPENMP) $(LFLAGS) $(LIBS) $^ -o $@
$(CAN): generate.o can.o
	$(MPIFC) $(LFLAGS) $(LIBS) $^ -o $@
$(CAN_HYB): generate.o can_hyb.o
	$(MPIFC) $(OPENMP) $(LFLAGS) $(LIBS) $^ -o $@
$(CHECK): check.o
	$(FC) $(OPENMP) $(LFLAGS) $(LIBS) $^ -o $@
clean:
	rm -f *.o *.mod *~ $(SERIAL) $(DGEMM) $(OMP) $(CAN) $(CAN_HYB) $(CREATE_INPUT) $(CHECK) a b c.seri c.omp c.dgemm c.can c.can_hyb
