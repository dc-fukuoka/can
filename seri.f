      program main
      implicit none
      include "param.f"
      integer*8::flop=imax*imax*imax*2
      real*8::a(imax,imax),b(imax,imax),c(imax,imax)
      integer::i,j,k
      real*8::t0,time,dclock
      real*8::trace

      open(unit=ia,file="a",     form="unformatted",access="stream")
      open(unit=ib,file="b",     form="unformatted",access="stream")
      open(unit=ic,file="c.seri",form="unformatted",access="stream")

      read(ia) a
      read(ib) b
      c = 0.0d0

      close(ia)
      close(ib)

      t0 = dclock()
c.... serial
      do j=1,imax
      do k=1,imax
      do i=1,imax
         c(i,j) = c(i,j) + a(i,k)*b(k,j)
      end do
      end do
      end do
      time = dclock() - t0
      write(6,*) "serial time:",time,
     &     flop/time/1000000000,"Gflops"

      trace = 0.0d0
      do i = 1, imax
         trace = trace + c(i,i)
      end do

      write(6,*) "trace:",trace

      write(ic) c
      close(ic)
      
      stop
      end program main
