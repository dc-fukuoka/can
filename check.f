      program main
      implicit none
      include "param.f"
      character(len=32)::file(2)
      real(8),dimension(imax*imax)::c1,c2
      real(8)::maxerr
      integer::i,j
      integer,parameter::iic(2)=(/100,200/)
      integer::iargc

      if (iargc().ne.2) then
         write(6,*) "two arguments are needed."
         stop
      end if

      do i=1,2
         call getarg(i, file(i))
      end do

      open(unit=iic(1),file=file(1),form="unformatted",access="stream")
      open(unit=iic(2),file=file(2),form="unformatted",access="stream")
      read(iic(1)) c1
      read(iic(2)) c2

      maxerr = 0.0d0
c$omp parallel do reduction(max:maxerr)
      do i=1,imax*imax
         maxerr = max(abs(c1(i)-c2(i)),maxerr)
      end do

      write(6,*) "maximum error:",maxerr

      close(iic(1))
      close(iic(2))
      
      stop
      end program main
