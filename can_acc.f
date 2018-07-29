      program main
      use openacc
      implicit none
      include "mpif.h"
      include "param.f"

      integer::iam,np,ierr
      integer*8::flop=int(imax,kind=8)**3*2
      real*8,allocatable::a_l(:,:),b_l(:,:),c_l(:,:)
      integer::imax_l
      integer::i,j,k,l
      integer::up,down,left,right,coords(2)
      integer::dims(2),prds(2)
      integer::istat(mpi_status_size)
      integer::cart_comm,cart_rank
      integer::src,dest
      integer::ireqs(4),istats(mpi_status_size,4)
      integer::index_i,index_j
      real*8::t0,time
      real*8::trace
      integer::sizes(2),subsizes(2),starts(2)
      integer::ifiletype
      character*16::file_a="a",file_b="b",file_c="c.can_acc"
      integer::fh_a,fh_b,fh_c
      integer(kind=mpi_offset_kind)::idisp
      integer::color,key,diag_comm,diag_rank
#ifdef _OPENACC
      integer::ndevs,mydev
#endif

      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,iam,ierr)
      call mpi_comm_size(mpi_comm_world,np, ierr)

      call acc_init(acc_device_nvidia)
      ndevs = acc_get_num_devices(acc_device_nvidia)
      mydev = mod(iam, ndevs)
      call acc_set_device_num(mydev, acc_device_nvidia)

      if (mod(imax*imax,np).ne.0.or.
     &     (dsqrt(dble(np))-dint(dsqrt(dble(np)))).ne.0.0d0) then
         if (iam.eq.0) then
            write(6,*) "imax/sqrt(np) must be an integer."
            write(6,*) "sqrt(np) must be an integer"
         end if
         call mpi_finalize(ierr)
         stop
      end if

      prds(1) = .true.
      prds(2) = prds(1)
      dims(1) = dsqrt(dble(np))
      dims(2) = dims(1)

      call mpi_cart_create(mpi_comm_world,2,dims,prds,.true.,
     &     cart_comm,ierr)
      call mpi_comm_rank(cart_comm,cart_rank,ierr)
      call mpi_cart_coords(cart_comm,cart_rank,2,coords,ierr)

      call mpi_cart_shift(cart_comm,0,1,up,down,ierr)
      call mpi_cart_shift(cart_comm,1,1,left,right,ierr)
      imax_l = imax/dsqrt(dble(np))

      allocate(a_l(imax_l,imax_l),b_l(imax_l,imax_l),
     &     c_l(imax_l,imax_l))
      
c.... datatype for parallel read/write
      sizes(1)    = imax
      sizes(2)    = imax
      subsizes(1) = imax_l
      subsizes(2) = imax_l
      starts(1)   = coords(1)*imax_l
      starts(2)   = coords(2)*imax_l

      call mpi_type_create_subarray(2,sizes,subsizes,starts,
     &     mpi_order_fortran,mpi_real8,ifiletype,ierr)
      call mpi_type_commit(ifiletype,ierr)

c.... read a and b
      call mpi_file_open(mpi_comm_world,file_a,mpi_mode_rdonly,
     &     mpi_info_null,fh_a,ierr)
      call mpi_file_open(mpi_comm_world,file_b,mpi_mode_rdonly,
     &     mpi_info_null,fh_b,ierr)
      idisp = 0
      call mpi_file_set_view(fh_a,idisp,mpi_real8,ifiletype,"native",
     &     mpi_info_null,ierr)
      call mpi_file_set_view(fh_b,idisp,mpi_real8,ifiletype,"native",
     &     mpi_info_null,ierr)
      call mpi_file_read_all(fh_a,a_l,imax_l*imax_l,mpi_real8,
     &     istat,ierr)
      call mpi_file_read_all(fh_b,b_l,imax_l*imax_l,mpi_real8,
     &     istat,ierr)
      call mpi_file_close(fh_a,ierr)
      call mpi_file_close(fh_b,ierr)

c.... initialize the order    process grid(C order...)
c a11 a12 a13    a11 a12 a13  |p0 p1 p2|
c a21 a22 a23 => a22 a23 a21  |p3 p4 p5|
c a31 a32 a33    a33 a31 a32  |p6 p7 p8|
c
c b11 b12 b13    b11 b22 b33  |p0 p1 p2|
c b21 b22 b23 => b21 b32 b13  |p3 p4 p5|
c b31 b32 b33    b31 b12 b23  |p6 p7 p8|
      call mpi_cart_shift(cart_comm,1,-coords(1),src,dest,ierr)
      call mpi_sendrecv_replace(a_l,imax_l*imax_l,mpi_real8,
     &     dest,0,src,0,cart_comm,istat,ierr)

      call mpi_cart_shift(cart_comm,0,-coords(2),src,dest,ierr)
      call mpi_sendrecv_replace(b_l,imax_l*imax_l,mpi_real8,
     &     dest,1,src,1,cart_comm,istat,ierr)
c.... looks ok
      call mpi_barrier(mpi_comm_world,ierr)
      t0 = mpi_wtime()

c.... main loop
c a : rotate left direction
c a11 a12 a13    a12 a13 a11    a13 a11 a12
c a22 a23 a21 => a23 a21 a22 => a21 a22 a23
c a33 a31 a32    a31 a32 a33    a32 a33 a31
c
c b : rotate up direction
c b11 b22 b33    b21 b32 b13    b31 b12 b23
c b21 b32 b13 => b31 b12 b23 => b11 b22 b33
c b31 b12 b23    b11 b22 b33    b21 b32 b13
c
c process grid(C order!)
c
c |p0 p1 p2|
c |p3 p4 p5|
c |p6 p7 p8|
c
c ex. process 0
c c11 = a11*b11 + a12*b21 + a13*a31 => OK
c process 1
c c12  = a12*b22 + a13*b32 + a11*b12
c      = a11*b12 + a12*b22 + a13*b32 => OK

c$acc enter data copyin(a_l(:,:),b_l(:,:))
c$acc& create(c_l(:,:))
      do l=1,dims(1)
c$acc host_data use_device(a_l(:,:),b_l(:,:))
         call mpi_sendrecv_replace(a_l,imax_l*imax_l,mpi_real8,
     &     left,0,right,0,cart_comm,istat,ierr)
         call mpi_sendrecv_replace(b_l,imax_l*imax_l,mpi_real8,
     &     up,1,down,1,cart_comm,istat,ierr)
c$acc end host_data
c$acc parallel private(i,j,k,l)
c$acc& present(a_l(:,:),b_l(:,:),c_l(:,:))        
c$acc loop
      do j=1,imax_l
      do k=1,imax_l
      do i=1,imax_l
         c_l(i,j) = c_l(i,j) + a_l(i,k)*b_l(k,j)
      end do
      end do
      end do
c$acc end parallel
      end do
c$acc exit data copyout(c_l(:,:))
c$acc& delete(a_l(:,:),b_l(:,:))
      call mpi_barrier(mpi_comm_world,ierr)
      time = mpi_wtime() - t0

      if (iam.eq.0) then
         write(6,*) "MPI time:",
     &        time,
     &        real(flop,8)/time/1000000000.0d0, "Gflops" 
      end if

c.... write the result
      call mpi_file_open(mpi_comm_world,file_c,
     &     mpi_mode_wronly+mpi_mode_create,mpi_info_null,fh_c,ierr)
      idisp = 0
      call mpi_file_set_view(fh_c,idisp,mpi_real8,ifiletype,"native",
     &     mpi_info_null,ierr)
      call mpi_file_write_all(fh_c,c_l,imax_l*imax_l,mpi_real8,
     &     istat,ierr)
      call mpi_file_close(fh_c,ierr)

c.... check result
      if (coords(1).eq.coords(2)) then
c.... diagonal in the process grid
         color = 0
         key   = coords(1)
      else
         color = 1
         key   = iam
      end if
      call mpi_comm_split(mpi_comm_world,color,key,diag_comm,ierr)
      call mpi_comm_rank(diag_comm,diag_rank,ierr)
      trace   = 0.0d0
c$acc parallel private(i)
c$acc loop reduction(+:trace)
      do i = 1,imax_l
         trace = trace + c_l(i,i)
      end do
c$acc end parallel
      if (diag_rank.eq.0) then
         call mpi_reduce(mpi_in_place,trace,1,mpi_real8,mpi_sum,
     &        0,diag_comm,ierr)
      else
         call mpi_reduce(trace,trace,1,mpi_real8,mpi_sum,
     &        0,diag_comm,ierr)
      end if
      
      if (iam.eq.0) write(6,*) "trace:",trace

c.... finalize
      call mpi_comm_free(cart_comm,ierr)
      call mpi_comm_free(diag_comm,ierr)
      call mpi_type_free(ifiletype,ierr)
      deallocate(a_l,b_l,c_l)

      call mpi_finalize(ierr)

      stop
      end program main
