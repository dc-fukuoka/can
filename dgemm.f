      program main
      implicit none
      include "param.f"
      integer*8::flop=imax*imax*imax*2
      real*8::a(imax,imax),b(imax,imax),c(imax,imax)
      integer::i,j,k
      real*8::t0,time
      real*8::dclock
      real*8::trace

      open(unit=ia,file="a",      form="unformatted",access="stream")
      open(unit=ib,file="b",      form="unformatted",access="stream")
      open(unit=ic,file="c.dgemm",form="unformatted",access="stream")

      read(ia) a
      read(ib) b
      c = 0.0d0

      close(ia)
      close(ib)

      t0 = dclock()
c.... dgemm
      call dgemm("n","n",imax,imax,imax,1.0d0,
     &     a,imax,b,imax,0.0d0,c,imax)

      time = dclock() - t0
      write(6,*) "dgemm time:",time,
     &     flop/time/1000000000,"Gflops"

      trace=0.0d0
c$omp parallel do reduction(+:trace)
      do i = 1,imax
         trace = trace + c(i,i)
      end do
      write(6,*) "trace:",trace
      
      write(ic) c
      close(ic)

      stop
      end program main
