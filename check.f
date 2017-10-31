      program main
      implicit none
      include "param.f"
      character(len=16)::file(2)
      real(8),dimension(imax*imax)::c1,c2
      real(8)::maxerr
      integer::i,j

      if (iargc().ne.2) then
         write(6,*) "two arguments are needed."
         stop
      end if

      do i=1,2
         call getarg(i, file(i))
      end do

      open(10000,file=file(1),form="unformatted")
      open(10001,file=file(2),form="unformatted")
      read(10000) c1
      read(10001) c2

      maxerr = 0.0d0
c$omp parallel do reduction(max:maxerr)
      do i=1,imax*imax
         maxerr = max(abs(c1(i)-c2(i)),maxerr)
      end do

      write(6,*) "maximum error:",maxerr

      close(10000)
      close(10001)
      
      stop
      end program main
