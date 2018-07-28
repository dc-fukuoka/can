      program main
      implicit none
      include "mpif.h"
      include "param.f"

      integer::iam,np,ierr
      integer*8::flop=imax*imax*imax*2
      real*8,allocatable::a(:,:),b(:,:),c(:,:)
      real*8,allocatable::a_l(:,:),b_l(:,:),c_l(:,:)
      real*8,allocatable::a_l_tmp(:,:),b_l_tmp(:,:),c_l_tmp(:,:)
      integer::imax_l
      integer::i,j,k,l
      integer::up,down,left,right,coords(2)
      integer::dims(2),prds(2)
      integer::istat(mpi_status_size)
      integer::cart_comm,cart_rank
      integer::src,dest
      integer::newtype
      integer::ireqs(4),istats(mpi_status_size,4)
      integer::index_i,index_j
      real*8::t0,time
      real*8::trace
      real*8::dclock
#ifdef _OPENMP
      integer::ireq,iprov
#endif

#ifdef _OPENMP
      ireq=mpi_thread_serialized
      call mpi_init_thread(ireq,iprov,ierr)
      if (iprov<ireq) then
         write(6,*) "mpi_thread_serialized is required."
         call mpi_finalize(ierr)
         stop
      end if
#else
      call mpi_init(ierr)
#endif
      call mpi_comm_rank(mpi_comm_world,iam,ierr)
      call mpi_comm_size(mpi_comm_world,np, ierr)

      if (mod(imax*imax,np).ne.0.or.
     &     (dsqrt(dble(np))-dint(dsqrt(dble(np)))).ne.0.0d0) then
         if (iam.eq.0) then
            write(6,*) "imax/sqrt(np) must be an integer."
            write(6,*) "sqrt(np) must be an integer"
         end if
         call mpi_finalize(ierr)
         stop
      end if

      if (iam.eq.0) then
         allocate(a(imax,imax),b(imax,imax),c(imax,imax))
c.... first touch for a,b,c
c$omp parallel do
         do j=1,imax
         do i=1,imax
            a(i,j) = 0.0d0
            b(i,j) = 0.0d0
            c(i,j) = 0.0d0
         end do
         end do

         call generate_randomno(a,imax,55555)
         call generate_randomno(b,imax,77777)
      end if

      imax_l = imax/dsqrt(dble(np))

      allocate(a_l(imax_l,imax_l),b_l(imax_l,imax_l),
     &     c_l(imax_l,imax_l))
      allocate(a_l_tmp(imax_l,imax_l),b_l_tmp(imax_l,imax_l),
     &     c_l_tmp(imax_l,imax_l))

c.... first touch for a_l,b_l,c_l,a_l_tmp,b_l_tmp,c_l_tmp
c$omp parallel do
      do j=1,imax_l
      do i=1,imax_l
         a_l(i,j) = 0.0d0
         b_l(i,j) = 0.0d0
         c_l(i,j) = 0.0d0
         a_l_tmp(i,j) = 0.0d0
         b_l_tmp(i,j) = 0.0d0
         c_l_tmp(i,j) = 0.0d0
      end do
      end do

c.... distribute glocal matrix into local matrices
      call mpi_type_vector(imax_l,imax_l,imax,mpi_real8,newtype,ierr)
      call mpi_type_commit(newtype,ierr)
      if (iam.eq.0) then
c$omp parallel do
         do j = 1,imax_l
         do i = 1,imax_l
            a_l(i,j) = a(i,j)
            b_l(i,j) = b(i,j)
         end do
         end do
         j=1
         k=1
         do i = 1, np-1
c            index_i=mod(imax_l,i)*imax_l+1
            index_i=k*imax_l+1
            if (i>=imax/imax_l) then
               index_i=(k-1)*imax_l+1
            end if
            if (mod(i,imax/imax_l).eq.0) then
               index_i=1
               j=j+1
               k=1
            end if

            index_j=(j-1)*imax_l+1
            k=k+1
            call mpi_send(a(index_i,index_j),1,newtype,i,
     &           0,mpi_comm_world,ierr)
            call mpi_send(b(index_i,index_j),1,newtype,i,
     &           1,mpi_comm_world,ierr)
         end do
      else
         call mpi_recv(a_l(1,1),imax_l*imax_l,mpi_real8,0,0,
     &        mpi_comm_world,istat,ierr)
         call mpi_recv(b_l(1,1),imax_l*imax_l,mpi_real8,0,1,
     &        mpi_comm_world,istat,ierr)
      end if

c.... define cartesian coordinte
      prds(1) = .true.
      prds(2) = prds(1)
      dims(1) = dsqrt(dble(np))
      dims(2) = dims(1)

      call mpi_cart_create(mpi_comm_world,2,dims,prds,.true.,
     &     cart_comm,ierr)
      call mpi_comm_rank(cart_comm,cart_rank,ierr)
      call mpi_cart_coords(cart_comm,cart_rank,2,coords,ierr)

      call mpi_cart_shift(cart_comm,0,1,left,right,ierr)
      call mpi_cart_shift(cart_comm,1,1,up,down,ierr)

c.... initialize the order
c a11 a12 a13    a11 a12 a13
c a21 a22 a23 => a22 a23 a21
c a31 a32 a33    a33 a31 a32
c
c b11 b12 b13    b11 b22 b33
c b21 b22 b23 => b21 b32 b13
c b31 b32 b33    b31 b12 b23
      call mpi_cart_shift(cart_comm,0,coords(2),src,dest,ierr)
      call mpi_sendrecv_replace(a_l(1,1),imax_l*imax_l,mpi_real8,
     &     src,0,dest,0,cart_comm,istat,ierr)

      call mpi_cart_shift(cart_comm,1,coords(1),src,dest,ierr)
      call mpi_sendrecv_replace(b_l(1,1),imax_l*imax_l,mpi_real8,
     &     src,1,dest,1,cart_comm,istat,ierr)

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
c ex. process 0
c c11 = a11*b11 + a12*b21 + a13*a31 => OK
c process 1
c c21  = a22*b21 + a23*b31 + a21*b11 
c      = a21*b11 + a22*b21 + a23*b31 => OK


c$omp parallel private(i,j,k,l)
      do l=1,dims(1)
c.... overlap
c$omp single
         call mpi_isend(a_l,imax_l*imax_l,mpi_real8,
     &        left,0,cart_comm,ireqs(1),ierr)
         call mpi_isend(b_l,imax_l*imax_l,mpi_real8,
     &        up,0,cart_comm,ireqs(2),ierr)
         call mpi_irecv(a_l_tmp,imax_l*imax_l,mpi_real8,
     &        right,0,cart_comm,ireqs(3),ierr)
         call mpi_irecv(b_l_tmp,imax_l*imax_l,mpi_real8,
     &        down,0,cart_comm,ireqs(4),ierr)
         call mpi_waitall(4,ireqs,istats,ierr)
c$omp end single nowait
c$omp do schedule(dynamic)
      do j=1,imax_l
      do k=1,imax_l
      do i=1,imax_l
         c_l(i,j) = c_l(i,j) + a_l(i,k)*b_l(k,j)
      end do
      end do
      end do
c$omp end do
c..... implicit barrier
c$omp do
      do j = 1,imax_l
      do i = 1,imax_l
         a_l(i,j) = a_l_tmp(i,j)
         b_l(i,j) = b_l_tmp(i,j)
      end do
      end do
c$omp end do
      end do
c$omp end parallel

      call mpi_barrier(mpi_comm_world,ierr)
      time = mpi_wtime() - t0

      if (iam.eq.0) then
         write(6,*) "hybrid time:",
     &        time,flop/time/1000000000, "Gflops" 
      end if

c.... gather local matrices into global one
      if (iam.eq.0) then
c$omp parallel do
         do j=1,imax_l
         do i=1,imax_l
            c(i,j) = c_l(i,j)
         end do
         end do

         j=1
         k=1
         do i = 1, np-1
c            index_i=mod(imax_l,i)*imax_l+1
            index_i=k*imax_l+1
            if (i>=imax/imax_l) then
               index_i=(k-1)*imax_l+1
            end if
            if (mod(i,imax/imax_l).eq.0) then
               index_i=1
               j=j+1
               k=1
            end if
            index_j=(j-1)*imax_l+1
            k=k+1
            call mpi_recv(c(index_i,index_j),1,newtype,i,0,
     &           mpi_comm_world,istat,ierr)
c            j=j+1
         end do
      else
         call mpi_send(c_l(1,1),imax_l*imax_l,mpi_real8,0,0,
     &        mpi_comm_world,ierr)
      end if

c.... check result
      if(iam.eq.0) then
         trace = 0.0d0
c$omp parallel do reduction(+:trace)
         do i = 1, imax
            trace = trace + c(i,i)
         end do
         write(6,*) "trace:",trace
         write(1000) c
      end if

c.... finalize
      call mpi_comm_free(cart_comm,ierr)
      call mpi_type_free(newtype,ierr)
      if (iam.eq.0) then
         deallocate(a,b,c)
      end if
      deallocate(a_l,b_l,c_l)
      deallocate(a_l_tmp,b_l_tmp,c_l_tmp)

      call mpi_finalize(ierr)

      stop
      end program main
