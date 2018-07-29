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
LDFLAGS  = 
LIBS    = -mkl
OPENMP  = -fopenmp
PREP    = -fpp
SRCS    = $(MKLROOT)/include/mkl_vsl.f90 generate.f seri.f dgemm.f omp.f

.SUFFIXES: .f.o

ALL: $(CREATE_INPUT) $(SERIAL) $(DGEMM) $(OMP) $(CAN) $(CAN_HYB) $(CHECK)

$(VSL): $(MKLROOT)/include/mkl_vsl.f90
	ifort -c $(FCFLAGS) $(MKLROOT)/include/mkl_vsl.f90

create_input.o: create_input.f param.f
	$(FC) -c $(FCFLAGS) $<
generate.o: generate.f $(VSL)
	$(FC) -c $(FCFLAGS) $<
seri.o: seri.f param.f
	$(FC) -c $(FCFLAGS) $<
dgemm.o: dgemm.f param.f
	$(FC) -c $(FCFLAGS) $<
omp.o: omp.f param.f
	$(FC) -c $(OPENMP) $(FCFLAGS) $<
can.o: can.f param.f
	$(MPIFC) -c $(FCFLAGS) $<
can_hyb.o: can_hyb.f param.f
	$(MPIFC) $(PREP) -c $(OPENMP) $(FCFLAGS) $<
check.o: check.f param.f
	$(FC) -c $(OPENMP) $(FCFLAGS) $<

$(CREATE_INPUT): generate.o create_input.o
	$(FC) $(LDFLAGS) $(LIBS) $^ -o $@
$(SERIAL): seri.o
	$(FC) $^ -o $@
$(DGEMM): dgemm.o
	$(FC) $(LDFLAGS) $(LIBS) $^ -o $@
$(OMP): omp.o
	$(FC) $(OPENMP) $^ -o $@
$(CAN): can.o
	$(MPIFC) $^ -o $@
$(CAN_HYB): can_hyb.o
	$(MPIFC) $(OPENMP) $^ -o $@
$(CHECK): check.o
	$(FC) $(OPENMP) $^ -o $@
clean:
	rm -f *.o *.mod *~ $(SERIAL) $(DGEMM) $(OMP) $(CAN) $(CAN_HYB) $(CREATE_INPUT) $(CHECK) a b c.seri c.omp c.dgemm c.can c.can_hyb
