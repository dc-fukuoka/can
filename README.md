can - A simple matrix-matrix mutiplication benchmark.
======
there are serial, intel MKL dgemm(), OpenMP, MPI and hybrid versions.  
MPI version is based on Cannon's algorithm.  
intel compiler and intel MKL library are needed.  
Input matrix is a psudorandom number, that is generated by intel MKL Mersenne Twister(MT19937)  
  
- binary names:  
  - serial: seri  
  - OpenMP: omp  
  - intel MKL dgemm(): dgemm  
  - MPI: can  
  - hybrid: can_hyb  
  
- matrix size: imax x imax (param.f)  
  
- Some notes for MPI and hybrid version:
  - imax/sqrt(np) must be an integer.
  - sqrt(np) must be an integer.

performance comparison(matrix size: 4096x4096, Intel(R) Xeon(R) CPU E5-2680 0 @ 2.70GHz, 8 cores/socket, 2 sockets/node, 4 nodes, 4xFDR infiniband):
-------

* serial
~~~
$ ./seri
 serial time:   12.8100309371948        10.7290102690491      Gflops
 trace:   4194330.21842071
 ~~~
* MKL dgemm() (single thread)
~~~
$ MKL_NUM_THREADS=1 ./dgemm
 dgemm time:   7.46148896217346        18.4197757537077      Gflops
 trace:   4194330.21842071
 ~~~
* MKL dgemm() (16 threads)
~~~
$ MKL_NUM_THREADS=16 KMP_AFFINITY=compact ./dgemm
 dgemm time:   1.30764889717102        105.103865241914      Gflops
 trace:   4194330.21842071
~~~
* OpenMP (16 threads)
~~~
$ OMP_NUM_THREADS=16 KMP_AFFINITY=compact ./omp
 omp time:  0.919203042984009        149.519689388573      Gflops
 trace:   4194330.21842071
~~~
* MPI
~~~
$ mpirun -hosts `nodeset -e -S"," selene[71-74]` -ppn 16 ./can
 MPI time:  0.723188161849976        190.045911592938      Gflops
 trace:   4194330.21842071
~~~
* hybrid
~~~
$ OMP_NUM_THREADS=4 KMP_AFFINITY=compact I_MPI_PIN_DOMAIN=omp mpirun -hosts `nodeset -e -S"," selene[71-74]` -ppn 4 ./can_hyb
 hybrid time:  0.671340942382812        204.723032359956      Gflops
 trace:   4194330.21842071
~~~

-------
larger scale test:  
with imax = 16384, 16 nodes, 256 cores:  

* MPI
~~~
$ srun --mpi=pmi2 -N16 -n256 -c1 --cpu_bind=cores -m block:block ./can
 MPI time:   50.1612188816071        175.356445045105      Gflops
 trace:   67110044.2091656
~~~
* hybrid
~~~
$ OMP_NUM_THREADS=4 KMP_AFFINITY=compact I_MPI_PIN_DOMAIN=omp srun --mpi=pmi2 -N16 -n64 -c4 --cpu_bind=cores -m block:block ./can_hyb
 hybrid time:   20.0198850631714        439.367808279245      Gflops
 trace:   67110044.2091656
~~~
