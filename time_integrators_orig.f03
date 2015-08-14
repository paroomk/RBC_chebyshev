module time_integrators

use types,  only: dp
use global, only: set_imex_params, fft_utils, &
                  alloc_err, Bml, Bmx, Tyx, kx, &
                  eye, p, ip
use fftw
use write_pack

implicit none

private
public imex_rk, backward_euler

real(dp) :: c1, c2, c3
real(dp) :: d11, d12, d13
real(dp) :: d21, d22, d23
real(dp) :: d31, d32, d33

contains

subroutine imex_rk(NC, NF, dt, t_final, kappa, PhysCheb, TM, PTM)

integer                 , intent(in)     :: NC, NF
real(dp),                 intent(in)     :: dt, t_final, kappa
real(dp), dimension(:,:), intent(in)     :: PhysCheb, TM, PTM
complex(dp), allocatable, dimension(:,:) :: K1, K2, K3
complex(dp), allocatable, dimension(:,:) :: PBml
real(dp),    allocatable, dimension(:)   :: out_r, out_c
real(dp)                                 :: dt_final, time
real(dp)                                 :: maxTyx
integer                                  :: i, j

! LAPACK and BLAS parameters
real(dp), parameter                      :: scale1=1.0_dp, scale2=0.0_dp
integer , parameter                      :: incx=1, incy=1

call set_imex_params(c1,c2,c3, d11,d21,d31, d12,d22,d32, d13,d23,d33)

allocate(K1(NC,NF/2+1)       , stat=alloc_err)
allocate(K2(NC,NF/2+1)       , stat=alloc_err)
allocate(K3(NC,NF/2+1)       , stat=alloc_err)
allocate(PBml(NC,NF/2+1)     , stat=alloc_err)
allocate(out_r(NC), out_c(NC), stat=alloc_err)

K1    = cmplx(0.0_dp, 0.0_dp)
K2    = cmplx(0.0_dp, 0.0_dp)
K3    = cmplx(0.0_dp, 0.0_dp)
PBml  = cmplx(0.0_dp, 0.0_dp)
out_r = 0.0_dp
out_c = 0.0_dp

time = 0.0_dp

open(unit=8000, file="maxval.txt", action="write", status="unknown", position="append")

do ! while time < t_final

   ! Do some computations using data from previous time step

   ! PBml = PTM*Bml
   do i = 1,NF/2+1
      call dgemv('n', NC, NC, scale1, PTM, NC, real (Bml(:,i)), incx, scale2, out_r, incy)
      call dgemv('n', NC, NC, scale1, PTM, NC, aimag(Bml(:,i)), incx, scale2, out_c, incy)
      PBml(:,i) = cmplx(out_r, out_c)
   end do

   ! Bring to physical space to track decay rate
   call decay(maxTyx, TM, PhysCheb)

   write(8000, fmt=2000) time, maxTyx

   ! Now move on to next time step
   dt_final = t_final - time

   if (dt_final <= dt) then
      time = t_final
   else
      time = time + dt
   end if

   if (time == t_final) then
      exit
   end if

   ! STAGE 1
   call stage1(K1,         dt, kappa, PTM, PBml)

   ! STAGE 2
   call stage2(K2, K1,     dt, kappa, PTM, PBml)

   ! STAGE 3
   call stage3(K3, K1, K2, dt, kappa, PTM, PBml)

   ! Update solution
   Bml = Bml + dt*(c1*K1 + c2*K2 + c3*K3)

end do

close(unit=8000)

2000 format(E25.16E3, E25.16E3)

end subroutine imex_rk

subroutine stage1(K1, dt, kappa, PTM, PBml)

  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  !
  ! Performs first stage of IMEX-RK 3-4-3 method
  !
  ! INPUT:
  !    dt    : time-step size
  !    kappa : thermal diffusivity
  !    PTM   : Projected temperature modes
  !    PBml  : PTM*Bml
  !
  ! OUTPUT:
  !    K1: First stage
  !
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

real(dp)   ,                              intent(in)  :: dt, kappa
real(dp)   , dimension(:,:),              intent(in)  :: PTM
complex(dp), dimension(:,:),              intent(in)  :: PBml
complex(dp), allocatable, dimension(:,:), intent(out) :: K1
real(dp)   , allocatable, dimension(:,:)              :: lhs
real(dp)   , allocatable, dimension(:,:)              :: rhs
real(dp)   , allocatable, dimension(:,:)              :: Kmat
complex(dp), allocatable, dimension(:)                :: temp
real(dp)   , allocatable, dimension(:)                :: K1r, K1c
integer                                               :: i
integer                                               :: NC, NF2

! Some LAPACK and BLAS variables
integer,     allocatable, dimension(:)                :: ipiv
integer                                               :: info
real(dp), parameter                                   :: scale1=1.0_dp, scale2=0.0_dp
integer , parameter                                   :: incx=1, incy=1

NC  = size(Bml,1)
NF2 = size(Bml,2)

allocate(K1(NC,NF2)      , stat=alloc_err)
allocate(lhs(NC,NC)      , stat=alloc_err)
allocate(rhs(NC,2)       , stat=alloc_err)
allocate(Kmat(NC,NC)     , stat=alloc_err)
allocate(K1r(NC), K1c(NC), stat=alloc_err)
allocate(ipiv(NC)        , stat=alloc_err)
allocate(temp(NC)    , stat=alloc_err)

K1   = cmplx(0.0_dp, 0.0_dp)
lhs  = 0.0_dp
rhs  = 0.0_dp
Kmat = 0.0_dp
K1r  = 0.0_dp
K1c  = 0.0_dp
ipiv = 0
temp = cmplx(0.0_dp, 0.0_dp)

do i = 1,NF2
   lhs = (1.0_dp + dt*d11*kappa*kx(i)**2.0_dp)*PTM - dt*d11*kappa*eye
   temp = PBml(:,i)
   rhs(:,1) = real (temp)
   rhs(:,2) = aimag(temp)
   ! Get Bml_1
   call dgesv(NC, 2, lhs, NC, ipiv, rhs, NC, info)
   ! Get K1
   ! Need to do multipliation of the form A(x+iy).  Do this in two parts:
   ! A(x+iy) = Ax + iAy => Two real multiplications.
   Kmat = kappa*(-kx(i)**2.0_dp*PTM + eye)
   call dgemv('n', NC, NC, scale1, Kmat, NC, rhs(:,1), incx, scale2, K1r, incy)
   call dgemv('n', NC, NC, scale1, Kmat, NC, rhs(:,2), incx, scale2, K1c, incy)
   K1(:,i) = cmplx(K1r, K1c)
end do

end subroutine stage1

subroutine stage2(K2, K1, dt, kappa, PTM, PBml)

  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  !
  ! Performs second stage of IMEX-RK 3-4-3 method
  !
  ! INPUT:
  !    K1    : Stage 1 result
  !    dt    : time-step size
  !    kappa : thermal diffusivity
  !    PTM   : Projected temperature modes
  !    PBml  : PTM*Bml
  !
  ! OUTPUT:
  !    K2: Second stage
  !
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

real(dp)   ,                              intent(in)  :: dt, kappa
real(dp)   , dimension(:,:),              intent(in)  :: PTM
complex(dp), dimension(:,:),              intent(in)  :: PBml
complex(dp), dimension(:,:),              intent(in)  :: K1
complex(dp), allocatable, dimension(:,:), intent(out) :: K2
real(dp)   , allocatable, dimension(:,:)              :: lhs
real(dp)   , allocatable, dimension(:,:)              :: rhs
real(dp)   , allocatable, dimension(:,:)              :: Kmat
complex(dp), allocatable, dimension(:)                :: temp
real(dp)   , allocatable, dimension(:)                :: K2r, K2c
integer                                               :: i
integer                                               :: NC, NF2

! Some LAPACK and BLAS variables
integer,     allocatable, dimension(:)                :: ipiv
integer                                               :: info
real(dp), parameter                                   :: scale1=1.0_dp, scale2=0.0_dp
integer , parameter                                   :: incx=1, incy=1

NC  = size(Bml,1)
NF2 = size(Bml,2)

allocate(K2(NC,NF2)      , stat=alloc_err)
allocate(lhs(NC,NC)      , stat=alloc_err)
allocate(rhs(NC,2)       , stat=alloc_err)
allocate(Kmat(NC,NC)     , stat=alloc_err)
allocate(K2r(NC), K2c(NC), stat=alloc_err)
allocate(ipiv(NC)        , stat=alloc_err)
allocate(temp(NC)        , stat=alloc_err)

K2   = cmplx(0.0_dp, 0.0_dp)
lhs  = 0.0_dp
rhs  = 0.0_dp
Kmat = 0.0_dp
K2r  = 0.0_dp
K2c  = 0.0_dp
ipiv = 0
temp = cmplx(0.0_dp, 0.0_dp)

do i = 1,NF2
   lhs = (1.0_dp + dt*d22*kappa*kx(i)**2.0_dp)*PTM - dt*d22*kappa*eye
   temp = PBml(:,i) + dt*d21*K1(:,i)
   rhs(:,1) = real (temp)
   rhs(:,2) = aimag(temp)
   ! Get Bml_2
   call dgesv(NC, 2, lhs, NC, ipiv, rhs, NC, info)
   ! Get K2
   ! Need to do multipliation of the form A(x+iy).  Do this in two parts:
   ! A(x+iy) = Ax + iAy => Two real multiplications.
   Kmat = kappa*(-kx(i)**2.0_dp*PTM + eye)
   call dgemv('n', NC, NC, scale1, Kmat, NC, rhs(:,1), incx, scale2, K2r, incy)
   call dgemv('n', NC, NC, scale1, Kmat, NC, rhs(:,2), incx, scale2, K2c, incy)
   K2(:,i) = cmplx(K2r, K2c)
end do

end subroutine stage2

subroutine stage3(K3, K1, K2, dt, kappa, PTM, PBml)

  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  !
  ! Performs second stage of IMEX-RK 3-4-3 method
  !
  ! INPUT:
  !    K1    : Stage 1 result
  !    K2    : Stage 2 result
  !    dt    : time-step size
  !    kappa : thermal diffusivity
  !    PTM   : Projected temperature modes
  !    PBml  : PTM*Bml
  !
  ! OUTPUT:
  !    K3: Third stage
  !
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

implicit none

real(dp)   ,                              intent(in)  :: dt, kappa
real(dp)   , dimension(:,:),              intent(in)  :: PTM
complex(dp), dimension(:,:),              intent(in)  :: PBml
complex(dp), dimension(:,:),              intent(in)  :: K1, K2
complex(dp), allocatable, dimension(:,:), intent(out) :: K3
real(dp)   , allocatable, dimension(:,:)              :: lhs
real(dp)   , allocatable, dimension(:,:)              :: rhs
real(dp)   , allocatable, dimension(:,:)              :: Kmat
complex(dp), allocatable, dimension(:)                :: temp
real(dp)   , allocatable, dimension(:)                :: K3r, K3c
integer                                               :: i
integer                                               :: NC, NF2

! Some LAPACK and BLAS variables
integer,     allocatable, dimension(:)                :: ipiv
integer                                               :: info
real(dp), parameter                                   :: scale1=1.0_dp, scale2=0.0_dp
integer , parameter                                   :: incx=1, incy=1

NC  = size(Bml,1)
NF2 = size(Bml,2)

allocate(K3(NC,NF2)      , stat=alloc_err)
allocate(lhs(NC,NC)      , stat=alloc_err)
allocate(rhs(NC,2)       , stat=alloc_err)
allocate(Kmat(NC,NC)     , stat=alloc_err)
allocate(K3r(NC), K3c(NC), stat=alloc_err)
allocate(ipiv(NC)        , stat=alloc_err)
allocate(temp(NC)        , stat=alloc_err)

K3   = cmplx(0.0_dp, 0.0_dp)
lhs  = 0.0_dp
rhs  = 0.0_dp
Kmat = 0.0_dp
K3r  = 0.0_dp
K3c  = 0.0_dp
ipiv = 0
temp = cmplx(0.0_dp, 0.0_dp)

do i = 1,NF2
   lhs = (1.0_dp + dt*d33*kappa*kx(i)**2.0_dp)*PTM - dt*d33*kappa*eye
   temp = PBml(:,i) + dt*(d31*K1(:,i) + d32*K2(:,i))
   rhs(:,1) = real (temp)
   rhs(:,2) = aimag(temp)
   ! Get Bml_3
   call dgesv(NC, 2, lhs, NC, ipiv, rhs, NC, info)
   ! Get K3
   ! Need to do multipliation of the form A(x+iy).  Do this in two parts:
   ! A(x+iy) = Ax + iAy => Two real multiplications.
   Kmat = kappa*(-kx(i)**2.0_dp*PTM + eye)
   call dgemv('n', NC, NC, scale1, Kmat, NC, rhs(:,1), incx, scale2, K3r, incy)
   call dgemv('n', NC, NC, scale1, Kmat, NC, rhs(:,2), incx, scale2, K3c, incy)
   K3(:,i) = cmplx(K3r, K3c)
end do

end subroutine stage3

subroutine decay(maxTyx, TM, PhysCheb)

implicit none

integer               :: NP, NF, NC

real(dp), dimension(:,:), intent(in)  :: TM, PhysCheb
real(dp),                 intent(out) :: maxTyx

! Some LAPACK and BLAS parameters
real(dp), parameter :: scale1=1.0_dp, scale2=0.0_dp

NC = size(Bml, 1)
NF = size(Bmx, 2)
NP = size(Tyx, 1)

Tyx = 0.0_dp
call fftw_execute_dft_c2r(ip, Bml, Bmx) ! Fourier to physical
call dgemm('n', 'n', NP, NF, NC, scale1, TM, NP, Bmx, NC, scale2, Tyx, NP) ! Cheb to physical

! Get max of Tyx and return it
maxTyx = maxval(Tyx)

! Now bring things back to Fourier-Chebyshev space
Bmx = 0.0_dp
call dgemm('n', 'n', NC, NF, NP, scale1, PhysCheb, NC, Tyx, NP, scale2, Bmx, NC) ! Phys to Cheb
call fftw_execute_dft_r2c(p, Bmx, Bml) ! Phys to Fourier
Bml = Bml / real(NF, dp)

end subroutine decay

subroutine backward_euler(NC,NF,dt,t_final,kappa, PhysCheb,TM,PTM)

integer                 , intent(in)     :: NC, NF
real(dp),                 intent(in)     :: dt, t_final, kappa
real(dp), dimension(:,:), intent(in)     :: PhysCheb, TM, PTM
real(dp)                                 :: dt_final, time
real(dp)                                 :: maxTyx
real(dp), allocatable, dimension(:,:)    :: lhs, rhs
real(dp), parameter                      :: scale1=1.0_dp, scale2=0.0_dp
integer,  parameter                      :: incx=1, incy=1
integer,  allocatable, dimension(:)      :: ipiv
integer                                  :: info
integer                                  :: i

allocate(lhs(NC,NC), stat=alloc_err)
allocate(rhs(NC,2) , stat=alloc_err)
allocate(ipiv(NC)  , stat=alloc_err)
lhs  = 0.0_dp
rhs  = 0.0_dp
ipiv = 0

time = 0.0_dp

call decay(maxTyx, TM, PhysCheb)

open(unit=8000, file="maxval.txt", action="write", status="unknown", position="append")
write(8000, fmt=2000) time, maxTyx

do ! while time < t_final

   dt_final = t_final - time

   if (dt_final <= dt) then
      time = t_final
   else
      time = time + dt
   end if

   do i = 1,NF/2+1
      lhs = (1.0_dp + kappa*dt*kx(i)**2.0_dp)*PTM - kappa*dt*eye
      ! Compute the RHS
      call dgemv('n', NC, NC, scale1, PTM, NC, real (Bml(:,i)), incx, scale2, rhs(:,1), incy)
      call dgemv('n', NC, NC, scale1, PTM, NC, aimag(Bml(:,i)), incx, scale2, rhs(:,2), incy)
      ! Solve for solution at next time
      call dgesv(NC, 2, lhs, NC, ipiv, rhs, NC, info)
      ! Update solution
      Bml(:,i) = cmplx(rhs(:,1), rhs(:,2))
   end do

   call decay(maxTyx, TM, PhysCheb)
   write(8000, fmt=2000) time, maxTyx

   if (time == t_final) then
      exit
   end if

end do

close(unit=8000)

2000 format(E25.16E3, E25.16E3)

end subroutine backward_euler

end module time_integrators
