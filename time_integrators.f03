module time_integrators

use types,  only: dp
use global, only: set_imex_params, fft_utils,    &
                  alloc_err, Aml, Bml, Amx, Bmx, &
                  Vyx, Tyx, kx, eye,             &
                  pV, ipV, pT, ipT
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

subroutine imex_rk(NC, NF, dt, t_final, nu, kappa,   &
                          PhysChebV, PhysChebT,VM,TM,&
                          GPVM,PTM,GPD2VM,GPD4VM)

integer                 , intent(in)     :: NC, NF
real(dp),                 intent(in)     :: dt, t_final, kappa, nu
real(dp), dimension(:,:), intent(in)     :: PhysChebV, PhysChebT
real(dp), dimension(:,:), intent(in)     :: VM, TM, PTM
real(dp), dimension(:,:), intent(in)     :: GPVM, GPD2VM, GPD4VM
complex(dp), allocatable, dimension(:,:) :: K1V, K2V, K3V
complex(dp), allocatable, dimension(:,:) :: K1T, K2T, K3T
real(dp)                                 :: dt_final, time
real(dp)                                 :: maxVyx, maxTyx
integer                                  :: i, j

! LAPACK and BLAS parameters
real(dp), parameter                      :: scale1=1.0_dp, scale2=0.0_dp
integer , parameter                      :: incx=1, incy=1

call set_imex_params(c1,c2,c3, d11,d21,d31, d12,d22,d32, d13,d23,d33)

allocate(K1T(NC,NF/2+1), K1V(NC,NF/2+1), stat=alloc_err)
allocate(K2T(NC,NF/2+1), K2V(NC,NF/2+1), stat=alloc_err)
allocate(K3T(NC,NF/2+1), K3V(NC,NF/2+1), stat=alloc_err)

K1V = cmplx(0.0_dp, 0.0_dp)
K2V = cmplx(0.0_dp, 0.0_dp)
K3V = cmplx(0.0_dp, 0.0_dp)

K1T = cmplx(0.0_dp, 0.0_dp)
K2T = cmplx(0.0_dp, 0.0_dp)
K3T = cmplx(0.0_dp, 0.0_dp)

time = 0.0_dp

open(unit=8000, file="maxval.txt", action="write", status="unknown", position="append")

do ! while time < t_final

   ! Do some computations using data from previous time step

   ! Bring to physical space to track decay rate
   call decay(maxVyx, maxTyx, VM, TM, PhysChebV, PhysChebT)

   write(8000, fmt=3000) time, maxVyx, maxTyx

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
   call stage1(K1V,K1T,                  dt,nu,kappa,GPVM,GPD2VM,GPD4VM,PTM)

   ! STAGE 2
   call stage2(K2V,K2T, K1V,K1T,         dt,nu,kappa,GPVM,GPD2VM,GPD4VM,PTM)

   ! STAGE 3
   call stage3(K3V,K3T, K1V,K1T,K2V,K2T, dt,nu,kappa,GPVM,GPD2VM,GPD4VM,PTM)

   ! Update solution
   Aml = Aml + dt*(c1*K1V + c2*K2V + c3*K3V)
   Bml = Bml + dt*(c1*K1T + c2*K2T + c3*K3T)

end do

close(unit=8000)

2000 format(E25.16E3, E25.16E3          )
3000 format(E25.16E3, E25.16E3, E25.16E3)

end subroutine imex_rk

subroutine stage1(K1V, K1T, dt,nu,kappa,GPVM,GPD2VM,GPD4VM,PTM)

  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  !
  ! Performs first stage of IMEX-RK 3-4-3 method
  !
  ! INPUT:
  !    dt     : time-step size
  !    nu     : kinematic viscosity
  !    kappa  : thermal diffusivity
  !    GPVM   : Projected temperature modes
  !    GPD2VM : Projected 2nd derivative temperature modes
  !    GPD4VM : Projected 4th derivative temperature modes
  !    PTM    : Projected temperature modes
  !
  ! OUTPUT:
  !    K1V: First stage for vertical velocity equation
  !    K1T: First stage for vertical temperature equation
  !
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

real(dp)   ,                              intent(in)  :: dt, nu, kappa
real(dp)   , dimension(:,:),              intent(in)  :: GPVM, GPD2VM, GPD4VM, PTM
complex(dp), allocatable, dimension(:,:), intent(out) :: K1V, K1T
real(dp)   , allocatable, dimension(:,:)              :: lhs
real(dp)   , allocatable, dimension(:,:)              :: rhs
real(dp)   , allocatable, dimension(:,:)              :: tempmat, ups
real(dp)   , allocatable, dimension(:,:)              :: Kmat
real(dp)   , allocatable, dimension(:)                :: Kr, Kc
real(dp)                                              :: wave, wave2, wave4
integer                                               :: i
integer                                               :: NC, NF2

! Some LAPACK and BLAS variables
integer,     allocatable, dimension(:)                :: ipiv
integer                                               :: info
real(dp), parameter                                   :: scale1=1.0_dp, scale2=0.0_dp
integer , parameter                                   :: incx=1, incy=1

NC  = size(Bml,1)
NF2 = size(Bml,2)

allocate(K1V(NC,NF2), K1T(NC,NF2)  , stat=alloc_err)
allocate(lhs(NC,NC)                , stat=alloc_err)
allocate(rhs(NC,2)                 , stat=alloc_err)
allocate(tempmat(NC,NC), ups(NC,NC), stat=alloc_err)
allocate(Kmat(NC,NC)               , stat=alloc_err)
allocate(Kr(NC), Kc(NC)            , stat=alloc_err)
allocate(ipiv(NC)                  , stat=alloc_err)

K1T     = cmplx(0.0_dp, 0.0_dp)
K1V     = cmplx(0.0_dp, 0.0_dp)
lhs     = 0.0_dp
rhs     = 0.0_dp
tempmat = 0.0_dp
ups     = 0.0_dp
Kmat    = 0.0_dp
Kr      = 0.0_dp
Kc      = 0.0_dp
ipiv    = 0

do i = 1,NF2
   ! Vertical velocity equation
   wave  = kx(i)
   wave2 = wave**2.0_dp
   wave4 = wave2**2.0_dp

   ! Will use this a few times so just pay the memory price
   tempmat = wave4*GPVM - 2.0_dp*wave2*GPD2VM + GPD4VM

   ! Upsilon matrix pre-multiplying solution in vertical velocity equation
   ! time derivative term
   ups = -wave2*GPVM + GPD2VM

   ! RHS for first stage
   rhs(:,1) = real (Aml(:,i))
   rhs(:,2) = aimag(Aml(:,i))

   ! Compute upsilon*rhs
   call dgemv('n', NC, NC, scale1, ups, NC, rhs(:,1), incx, scale2, Kr, incy)
   call dgemv('n', NC, NC, scale1, ups, NC, rhs(:,2), incx, scale2, Kc, incy)
   rhs(:,1) = Kr
   rhs(:,2) = Kc

   ! Now form lhs matrix (lhs*a1 = ups*a_{n-1})
   lhs = ups - dt*d11*nu*tempmat
   
   ! Get Aml_1
   call dgesv(NC, 2, lhs, NC, ipiv, rhs, NC, info)

   ! Get K1 for vertical velocity equation (avoids ever forming inv(upsilon))
   call dgemv('n', NC, NC, scale1, tempmat, NC, rhs(:,1), incx, scale2, Kr, incy)
   call dgemv('n', NC, NC, scale1, tempmat, NC, rhs(:,2), incx, scale2, Kc, incy)
   rhs(:,1) = Kr
   rhs(:,2) = Kc
   call dgesv(NC, 2, ups, NC, ipiv, rhs, NC, info)

   K1V(:,i) = nu*cmplx(rhs(:,1), rhs(:,2))

   ! Temperature equation
   lhs = (1.0_dp + dt*d11*kappa*wave2)*PTM - dt*d11*kappa*eye
   rhs(:,1) = real (Bml(:,i))
   rhs(:,2) = aimag(Bml(:,i))
   ! Get Bml_1
   call dgesv(NC, 2, lhs, NC, ipiv, rhs, NC, info)
   ! Get K1 for temperature equation
   ! Need to do multipliation of the form A(x+iy).  Do this in two parts:
   ! A(x+iy) = Ax + iAy => Two real multiplications.
   Kmat = kappa*(-wave2*PTM + eye)
   call dgemv('n', NC, NC, scale1, Kmat  , NC, rhs(:,1), incx, scale2, Kr, incy)
   call dgemv('n', NC, NC, scale1, Kmat  , NC, rhs(:,2), incx, scale2, Kc, incy)
   K1T(:,i) = cmplx(Kr, Kc)
end do

end subroutine stage1

subroutine stage2(K2V,K2T, K1V,K1T,dt,nu,kappa,GPVM,GPD2VM,GPD4VM,PTM)

  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  ! Performs second stage of IMEX-RK 3-4-3 method
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  !
  ! INPUT:
  !    K1V    : First stage for vertical velocity equation
  !    K1T    : First stage for vertical temperature equation
  !    dt     : time-step size
  !    nu     : kinematic viscosity
  !    kappa  : thermal diffusivity
  !    GPVM   : Projected temperature modes
  !    GPD2VM : Projected 2nd derivative temperature modes
  !    GPD4VM : Projected 4th derivative temperature modes
  !    PTM    : Projected temperature modes
  !
  ! OUTPUT:
  !    K2V: Second stage for vertical velocity equation
  !    K2T: Second stage for vertical temperature equation
  !
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

real(dp)   ,                              intent(in)  :: dt, nu, kappa
real(dp)   , dimension(:,:),              intent(in)  :: GPVM, GPD2VM, GPD4VM, PTM
complex(dp), allocatable, dimension(:,:), intent(in)  :: K1V, K1T
complex(dp), allocatable, dimension(:,:), intent(out) :: K2V, K2T
real(dp)   , allocatable, dimension(:,:)              :: lhs
real(dp)   , allocatable, dimension(:,:)              :: rhs
real(dp)   , allocatable, dimension(:,:)              :: tempmat, ups
real(dp)   , allocatable, dimension(:,:)              :: Kmat
real(dp)   , allocatable, dimension(:)                :: Kr, Kc
complex(dp), allocatable, dimension(:)                :: temp
real(dp)                                              :: wave, wave2, wave4
integer                                               :: i
integer                                               :: NC, NF2

! Some LAPACK and BLAS variables
integer,     allocatable, dimension(:)                :: ipiv
integer                                               :: info
real(dp), parameter                                   :: scale1=1.0_dp, scale2=0.0_dp
integer , parameter                                   :: incx=1, incy=1

NC  = size(Bml,1)
NF2 = size(Bml,2)

allocate(K2V(NC,NF2), K2T(NC,NF2)  , stat=alloc_err)
allocate(lhs(NC,NC)                , stat=alloc_err)
allocate(rhs(NC,2)                 , stat=alloc_err)
allocate(tempmat(NC,NC), ups(NC,NC), stat=alloc_err)
allocate(Kmat(NC,NC)               , stat=alloc_err)
allocate(Kr(NC), Kc(NC)            , stat=alloc_err)
allocate(ipiv(NC)                  , stat=alloc_err)
allocate(temp(NC)                  , stat=alloc_err)

K2T     = cmplx(0.0_dp, 0.0_dp)
K2V     = cmplx(0.0_dp, 0.0_dp)
lhs     = 0.0_dp
rhs     = 0.0_dp
tempmat = 0.0_dp
ups     = 0.0_dp
Kmat    = 0.0_dp
Kr      = 0.0_dp
Kc      = 0.0_dp
ipiv    = 0
temp    = cmplx(0.0_dp, 0.0_dp)

do i = 1,NF2
   ! Vertical velocity equation
   wave  = kx(i)
   wave2 = wave**2.0_dp
   wave4 = wave2**2.0_dp

   ! Will use this a few times so just pay the memory price
   tempmat = wave4*GPVM - 2.0_dp*wave2*GPD2VM + GPD4VM

   ! Upsilon matrix pre-multiplying solution in vertical velocity equation
   ! time derivative term
   ups = -wave2*GPVM + GPD2VM

   ! RHS for first stage
   temp = Aml(:,i) + dt*d21*K1V(:,i)
   rhs(:,1) = real (temp)
   rhs(:,2) = aimag(temp)

   ! Compute upsilon*rhs
   call dgemv('n', NC, NC, scale1, ups, NC, rhs(:,1), incx, scale2, Kr, incy)
   call dgemv('n', NC, NC, scale1, ups, NC, rhs(:,2), incx, scale2, Kc, incy)
   rhs(:,1) = Kr
   rhs(:,2) = Kc

   ! Now form lhs matrix (lhs*a1 = ups*a_{n-1})
   lhs = ups - dt*d22*nu*tempmat
   
   ! Get Aml_2
   call dgesv(NC, 2, lhs, NC, ipiv, rhs, NC, info)

   ! Get K2 for vertical velocity equation (avoids every forming inv(upsilon)
   ! and probably faster too)
   call dgemv('n', NC, NC, scale1, tempmat, NC, rhs(:,1), incx, scale2, Kr, incy)
   call dgemv('n', NC, NC, scale1, tempmat, NC, rhs(:,2), incx, scale2, Kc, incy)
   rhs(:,1) = Kr
   rhs(:,2) = Kc
   call dgesv(NC, 2, ups, NC, ipiv, rhs, NC, info)

   K2V(:,i) = nu*cmplx(rhs(:,1), rhs(:,2))

   ! Temperature equation
   lhs = (1.0_dp + dt*d22*kappa*wave2)*PTM - dt*d22*kappa*eye
   temp = Bml(:,i) + dt*d21*K1T(:,i)
   rhs(:,1) = real (temp)
   rhs(:,2) = aimag(temp)
   ! Get Bml_2
   call dgesv(NC, 2, lhs, NC, ipiv, rhs, NC, info)
   ! Get K2 for temperature equation
   ! Need to do multipliation of the form A(x+iy).  Do this in two parts:
   ! A(x+iy) = Ax + iAy => Two real multiplications.
   Kmat = kappa*(-wave2*PTM + eye)
   call dgemv('n', NC, NC, scale1, Kmat, NC, rhs(:,1), incx, scale2, Kr, incy)
   call dgemv('n', NC, NC, scale1, Kmat, NC, rhs(:,2), incx, scale2, Kc, incy)
   K2T(:,i) = cmplx(Kr, Kc)
end do

end subroutine stage2

subroutine stage3(K3V,K3T, K1V,K1T,K2V,K2T, dt,nu,kappa,GPVM,GPD2VM,GPD4VM,PTM)

  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  ! Performs third stage of IMEX-RK 3-4-3 method
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  !
  ! INPUT:
  !    K1V    : First stage for vertical velocity equation
  !    K1T    : First stage for vertical temperature equation
  !    K2V    : Second stage for vertical velocity equation
  !    K2T    : Second stage for vertical temperature equation
  !    dt     : time-step size
  !    nu     : kinematic viscosity
  !    kappa  : thermal diffusivity
  !    GPVM   : Projected temperature modes
  !    GPD2VM : Projected 2nd derivative temperature modes
  !    GPD4VM : Projected 4th derivative temperature modes
  !    PTM    : Projected temperature modes
  !
  ! OUTPUT:
  !    K3V: Third stage for vertical velocity equation
  !    K3T: Third stage for vertical temperature equation
  !
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

implicit none

real(dp)   ,                              intent(in)  :: dt, nu, kappa
real(dp)   , dimension(:,:),              intent(in)  :: GPVM, GPD2VM, GPD4VM, PTM
complex(dp), allocatable, dimension(:,:), intent(in)  :: K1V, K1T, K2V, K2T
complex(dp), allocatable, dimension(:,:), intent(out) :: K3V, K3T
real(dp)   , allocatable, dimension(:,:)              :: lhs
real(dp)   , allocatable, dimension(:,:)              :: rhs
real(dp)   , allocatable, dimension(:,:)              :: tempmat, ups
real(dp)   , allocatable, dimension(:,:)              :: Kmat
real(dp)   , allocatable, dimension(:)                :: Kr, Kc
complex(dp), allocatable, dimension(:)                :: temp
real(dp)                                              :: wave, wave2, wave4
integer                                               :: i
integer                                               :: NC, NF2

! Some LAPACK and BLAS variables
integer,     allocatable, dimension(:)                :: ipiv
integer                                               :: info
real(dp), parameter                                   :: scale1=1.0_dp, scale2=0.0_dp
integer , parameter                                   :: incx=1, incy=1

NC  = size(Bml,1)
NF2 = size(Bml,2)

allocate(K3V(NC,NF2), K3T(NC,NF2)  , stat=alloc_err)
allocate(lhs(NC,NC)                , stat=alloc_err)
allocate(rhs(NC,2)                 , stat=alloc_err)
allocate(tempmat(NC,NC), ups(NC,NC), stat=alloc_err)
allocate(Kmat(NC,NC)               , stat=alloc_err)
allocate(Kr(NC), Kc(NC)            , stat=alloc_err)
allocate(ipiv(NC)                  , stat=alloc_err)
allocate(temp(NC)                  , stat=alloc_err)

K3T     = cmplx(0.0_dp, 0.0_dp)
K3V     = cmplx(0.0_dp, 0.0_dp)
lhs     = 0.0_dp
rhs     = 0.0_dp
tempmat = 0.0_dp
ups     = 0.0_dp
Kmat    = 0.0_dp
Kr      = 0.0_dp
Kc      = 0.0_dp
ipiv    = 0
temp    = cmplx(0.0_dp, 0.0_dp)

do i = 1,NF2
   ! Vertical velocity equation
   wave  = kx(i)
   wave2 = wave**2.0_dp
   wave4 = wave2**2.0_dp

   ! Will use this a few times so just pay the memory price
   tempmat = wave4*GPVM - 2.0_dp*wave2*GPD2VM + GPD4VM

   ! Upsilon matrix pre-multiplying solution in vertical velocity equation
   ! time derivative term
   ups = -wave2*GPVM + GPD2VM

   ! RHS for first stage
   temp = Aml(:,i) + dt*(d31*K1V(:,i) + d32*K2V(:,i))
   rhs(:,1) = real (temp)
   rhs(:,2) = aimag(temp)

   ! Compute upsilon*rhs
   call dgemv('n', NC, NC, scale1, ups, NC, rhs(:,1), incx, scale2, Kr, incy)
   call dgemv('n', NC, NC, scale1, ups, NC, rhs(:,2), incx, scale2, Kc, incy)
   rhs(:,1) = Kr
   rhs(:,2) = Kc

   ! Now form lhs matrix (lhs*a1 = ups*a_{n-1})
   lhs = ups - dt*d33*nu*tempmat
   
   ! Get Aml_3
   call dgesv(NC, 2, lhs, NC, ipiv, rhs, NC, info)

   ! Get K3 for vertical velocity equation (avoids every forming inv(upsilon)
   ! and probably faster too)
   call dgemv('n', NC, NC, scale1, tempmat, NC, rhs(:,1), incx, scale2, Kr, incy)
   call dgemv('n', NC, NC, scale1, tempmat, NC, rhs(:,2), incx, scale2, Kc, incy)
   rhs(:,1) = Kr
   rhs(:,2) = Kc
   call dgesv(NC, 2, ups, NC, ipiv, rhs, NC, info)

   K3V(:,i) = nu*cmplx(rhs(:,1), rhs(:,2))

   ! Temperature equation
   lhs = (1.0_dp + dt*d33*kappa*wave2)*PTM - dt*d33*kappa*eye
   temp = Bml(:,i) + dt*(d31*K1T(:,i) + d32*K2T(:,i))
   rhs(:,1) = real (temp)
   rhs(:,2) = aimag(temp)
   ! Get Bml_3
   call dgesv(NC, 2, lhs, NC, ipiv, rhs, NC, info)
   ! Get K3 for temperature equation
   ! Need to do multipliation of the form A(x+iy).  Do this in two parts:
   ! A(x+iy) = Ax + iAy => Two real multiplications.
   Kmat = kappa*(-wave2*PTM + eye)
   call dgemv('n', NC, NC, scale1, Kmat, NC, rhs(:,1), incx, scale2, Kr, incy)
   call dgemv('n', NC, NC, scale1, Kmat, NC, rhs(:,2), incx, scale2, Kc, incy)
   K3T(:,i) = cmplx(Kr, Kc)
end do

end subroutine stage3

subroutine decay(maxVyx, maxTyx, VM, TM, PhysChebV, PhysChebT)

implicit none

integer               :: NP, NF, NC

real(dp), dimension(:,:), intent(in)  :: VM, TM
real(dp), dimension(:,:), intent(in)  :: PhysChebV, PhysChebT
real(dp),                 intent(out) :: maxVyx, maxTyx

! Some LAPACK and BLAS parameters
real(dp), parameter :: scale1=1.0_dp, scale2=0.0_dp

! Use temperature arrays to get these sizes.  Could use V arrays too.
NC = size(Bml, 1)
NF = size(Bmx, 2)
NP = size(Tyx, 1)

! Bring fields to physical space
Vyx = 0.0_dp
Tyx = 0.0_dp
call fftw_execute_dft_c2r(ipV, Aml, Amx) ! Vertical velocity Fourier to physical
call fftw_execute_dft_c2r(ipT, Bml, Bmx) ! Temperature Fourier to physical
call dgemm('n', 'n', NP, NF, NC, scale1, VM, NP, Amx, NC, scale2, Vyx, NP) ! Cheb to physical (VV)
call dgemm('n', 'n', NP, NF, NC, scale1, TM, NP, Bmx, NC, scale2, Tyx, NP) ! Cheb to physical (T)

! Get max of Vyx and Tyx and return them
maxVyx = maxval(Vyx)
maxTyx = maxval(Tyx)

! Now bring things back to Fourier-Chebyshev space
Amx = 0.0_dp
Bmx = 0.0_dp
call dgemm('n', 'n', NC, NF, NP, scale1, PhysChebV, NC, Vyx, NP, scale2, Amx, NC) ! Phys to Cheb (VV)
call dgemm('n', 'n', NC, NF, NP, scale1, PhysChebT, NC, Tyx, NP, scale2, Bmx, NC) ! Phys to Cheb (T)
call fftw_execute_dft_r2c(pV, Amx, Aml) ! Phys to Fourier (VV)
call fftw_execute_dft_r2c(pT, Bmx, Bml) ! Phys to Fourier (T)
Aml = Aml / real(NF, dp)
Bml = Bml / real(NF, dp)

end subroutine decay

subroutine backward_euler(NC,NF,dt,t_final,nu,kappa, &
                          PhysChebV, PhysChebT,VM,TM,&
                          GPVM,PTM,GPD2VM,GPD4VM)

integer                 , intent(in)     :: NC, NF
real(dp),                 intent(in)     :: dt, t_final, kappa, nu
real(dp), dimension(:,:), intent(in)     :: PhysChebV, PhysChebT
real(dp), dimension(:,:), intent(in)     :: VM, TM, PTM
real(dp), dimension(:,:), intent(in)     :: GPVM, GPD2VM, GPD4VM
real(dp)                                 :: dt_final, time
real(dp)                                 :: wave, wave2, wave4
real(dp)                                 :: maxVyx, maxTyx
real(dp), allocatable, dimension(:,:)    :: lhs, rhs
real(dp), allocatable, dimension(:,:)    :: upsilon
real(dp), parameter                      :: scale1=1.0_dp, scale2=0.0_dp
integer,  parameter                      :: incx=1, incy=1
integer,  allocatable, dimension(:)      :: ipiv
integer                                  :: info
integer                                  :: i, j, k

allocate(upsilon(NC,NC), stat=alloc_err)
allocate(lhs(NC,NC)    , stat=alloc_err)
allocate(rhs(NC,2)     , stat=alloc_err)
allocate(ipiv(NC)      , stat=alloc_err)
upsilon = 0.0_dp
lhs     = 0.0_dp
rhs     = 0.0_dp
ipiv    = 0

! Initial time
time    = 0.0_dp

call decay(maxVyx, maxTyx, VM, TM, PhysChebV, PhysChebT)

open(unit=8000, file="maxval.txt", action="write", status="unknown", position="append")
write(8000, fmt=3000) time, maxVyx, maxTyx

do ! while time < t_final

   dt_final = t_final - time

   if (dt_final <= dt) then
      time = t_final
   else
      time = time + dt
   end if

   do i = 1,NF/2+1
      ! Vertical velocity equation
      wave  = kx(i)
      wave2 = wave**2.0_dp
      wave4 = wave2**2.0_dp
      upsilon = -wave2*GPVM + GPD2VM
      lhs     = upsilon - nu*dt*(wave4*GPVM - 2.0_dp*wave2*GPD2VM + GPD4VM)
      ! Compute the RHS
      call dgemv('n', NC, NC, scale1, upsilon, NC, real (Aml(:,i)), incx, scale2, rhs(:,1), incy)
      call dgemv('n', NC, NC, scale1, upsilon, NC, aimag(Aml(:,i)), incx, scale2, rhs(:,2), incy)
      ! Solve for the solution at next time step
      call dgesv(NC, 2, lhs, NC, ipiv, rhs, NC, info)
      ! Update solution
      Aml(:,i) = cmplx(rhs(:,1), rhs(:,2))

      ! Temperature Equation
      lhs = (1.0_dp + kappa*dt*wave2)*PTM - kappa*dt*eye
      ! Compute the RHS
      call dgemv('n', NC, NC, scale1, PTM, NC, real (Bml(:,i)), incx, scale2, rhs(:,1), incy)
      call dgemv('n', NC, NC, scale1, PTM, NC, aimag(Bml(:,i)), incx, scale2, rhs(:,2), incy)
      ! Solve for solution at next time
      call dgesv(NC, 2, lhs, NC, ipiv, rhs, NC, info)
      ! Update solution
      Bml(:,i) = cmplx(rhs(:,1), rhs(:,2))
   end do

   call decay(maxVyx, maxTyx, VM, TM, PhysChebV, PhysChebT)
   write(8000, fmt=3000) time, maxVyx, maxTyx

   if (time == t_final) then
      exit
   end if

end do

close(unit=8000)

2000 format(E25.16E3, E25.16E3          )
3000 format(E25.16E3, E25.16E3, E25.16E3)

end subroutine backward_euler

end module time_integrators
