CREATE_INPUT = create_input
SERI         = seri
CAN_ACC      = can_acc
VSL          = mkl_vsl.mod
CHECK        = check

FC      = pgfortran
MPIFC   = mpif90
FCFLAGS = -O4 -Mvect=simd:256 -Mfma -mcmodel=medium
LDFLAGS = -L$(MKLROOT)/lib/intel64 -L$(MKLROOT)/../compiler/lib/intel64
LIBS    = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
OPENMP  = -mp
OPENACC = -acc -ta=tesla:cc60
PREP    = -Mpreprocess
SRCS    = $(MKLROOT)/include/mkl_vsl.f90 generate.f dclock.f check.f can_acc.f param.f

.SUFFIXES: .f.o

ALL: $(CREATE_INPUT) $(SERI) $(CAN_ACC) $(CHECK)

$(VSL): $(MKLROOT)/include/mkl_vsl.f90
	$(FC) -c $(FCFLAGS) $(MKLROOT)/include/mkl_vsl.f90
dclock.o: dclock.f
	$(FC) -c $(FCFLAGS) $<
create_input.o: create_input.f param.f
	$(FC) -c $(PREP) $(FCFLAGS) $<
generate.o: generate.f $(VSL) param.f
	$(FC) -c $(FCFLAGS) $<
seri.o: seri.f param.f
	$(FC) -c $(PREP) $(FCFLAGS) $<
can_acc.o: can_acc.f param.f
	$(MPIFC) $(PREP) -c $(OPENACC) $(FCFLAGS) $<
check.o: check.f param.f
	$(FC) -c $(PREP) $(OPENMP) $(FCFLAGS) $<

$(CREATE_INPUT): generate.o create_input.o
	$(FC) $(LDFLAGS) $(LIBS) $^ -o $@
$(SERI): seri.o dclock.o
	$(FC) $^ -o $@
$(CAN_ACC): can_acc.o dclock.o
	$(MPIFC) $(OPENACC) $^ -o $@
$(CHECK): check.o
	$(FC) $(OPENMP) $^ -o $@
clean:
	rm -f *.o *.mod *~ $(SERI) $(CAN_ACC) $(CREATE_INPUT) $(CHECK) a b c.can_acc
