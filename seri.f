      program main
      implicit none
      include "param.f"
      integer*8::flop=imax*imax*imax*2
      real*8::a(imax,imax),b(imax,imax),c(imax,imax)
      integer::i,j,k
      real*8::t0,time,dclock
      real*8::trace

      c(:,:) = 0.0d0
      call generate_randomno(a,imax,55555)
      call generate_randomno(b,imax,77777)

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

      write(55) c
c         do i = 1,imax
c            write(55,"(1024(1pe24.15))") (c(i,j),j=1,imax)
c            write(32768,"(65536(1pe14.5))") (c(i,j),j=1,imax)
c         end do
      
      stop
      end program main
