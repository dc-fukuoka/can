      subroutine generate_randomno(a_in, n, seed)
      use mkl_vsl_type
      use mkl_vsl
      implicit none
      integer,intent(in) :: n
      real*8,intent(inout)::a_in(n, n)
      integer,intent(in) :: seed
      integer::ierr
      integer::brng, method
      type(vsl_stream_state)::stream

      brng   = vsl_brng_mt19937
c     for older intel MKL
c     method = vsl_method_duniform_std
      method = vsl_rng_method_uniform_std_accurate

      ierr = vslnewstream(stream, brng, seed)
      ierr = vdrnguniform(method, stream, n*n, a_in, 0.0d0, 1.0d0)
      ierr = vsldeletestream(stream)

      return
      end subroutine generate_randomno
