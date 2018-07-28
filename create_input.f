      program main
      implicit none
      include "param.f"
      real*8,dimension(imax,imax) :: a, b

      call generate_randomno(a, imax, ia)
      call generate_randomno(b, imax, ib)

      open(unit=ia, file="a", form="unformatted", access="stream")
      open(unit=ib, file="b", form="unformatted", access="stream")

      write(ia) a
      write(ib) b

      close(ia)
      close(ib)
      end program main
