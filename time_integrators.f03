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
real(dp) :: dh11, dh12, dh13, dh14
real(dp) :: dh21, dh22, dh23, dh24
real(dp) :: dh31, dh32, dh33, dh34
real(dp) :: dh41, dh42, dh43, dh44

contains

subroutine imex_rk(NC, NF, dt, t_final, nu, kappa, &
                   PhysChebV, PhysChebT,VM,TM,     &
                   GPVM,GPTM,PVM,PTM,GPD2VM,GPD4VM)

integer                 , intent(in)     :: NC, NF
real(dp),                 intent(in)     :: dt, t_final, kappa, nu
real(dp), dimension(:,:), intent(in)     :: PhysChebV, PhysChebT
real(dp), dimension(:,:), intent(in)     :: VM, TM, PVM, PTM, GPTM
real(dp), dimension(:,:), intent(in)     :: GPVM, GPD2VM, GPD4VM
complex(dp), allocatable, dimension(:,:) :: K1V, K2V, K3V
complex(dp), allocatable, dimension(:,:) :: K1T, K2T, K3T
complex(dp), allocatable, dimension(:,:) :: K1hV, K2hV, K3hV, K4hV
complex(dp), allocatable, dimension(:,:) :: K1hT, K2hT, K3hT, K4hT
real(dp)                                 :: dt_final, time
real(dp)                                 :: maxVyx, maxTyx
integer                                  :: i, j, tstep

! LAPACK and BLAS parameters
real(dp), parameter                      :: scale1=1.0_dp, scale2=0.0_dp
integer , parameter                      :: incx=1, incy=1

call set_imex_params(c1,c2,c3, d11,d21,d31, d12,d22,d32, d13,d23,d33, &
                     dh11,dh12,dh13,dh14, dh21,dh22,dh23,dh24,        &
                     dh31,dh32,dh33,dh34, dh41,dh42,dh43,dh44)

allocate(K1T(NC,NF/2+1), K1V(NC,NF/2+1), stat=alloc_err)
allocate(K2T(NC,NF/2+1), K2V(NC,NF/2+1), stat=alloc_err)
allocate(K3T(NC,NF/2+1), K3V(NC,NF/2+1), stat=alloc_err)

allocate(K1hT(NC,NF/2+1), K1hV(NC,NF/2+1), stat=alloc_err)
allocate(K2hT(NC,NF/2+1), K2hV(NC,NF/2+1), stat=alloc_err)
allocate(K3hT(NC,NF/2+1), K3hV(NC,NF/2+1), stat=alloc_err)
allocate(K4hT(NC,NF/2+1), K4hV(NC,NF/2+1), stat=alloc_err)

K1V = cmplx(0.0_dp, 0.0_dp)
K2V = cmplx(0.0_dp, 0.0_dp)
K3V = cmplx(0.0_dp, 0.0_dp)

K1T = cmplx(0.0_dp, 0.0_dp)
K2T = cmplx(0.0_dp, 0.0_dp)
K3T = cmplx(0.0_dp, 0.0_dp)

K1hV = cmplx(0.0_dp, 0.0_dp)
K2hV = cmplx(0.0_dp, 0.0_dp)
K3hV = cmplx(0.0_dp, 0.0_dp)
K4hV = cmplx(0.0_dp, 0.0_dp)

K1hT = cmplx(0.0_dp, 0.0_dp)
K2hT = cmplx(0.0_dp, 0.0_dp)
K3hT = cmplx(0.0_dp, 0.0_dp)
K4hT = cmplx(0.0_dp, 0.0_dp)

time  = 0.0_dp
tstep = 0

open(unit=8000, file="maxval.txt", action="write", status="unknown", position="append")

do ! while time < t_final

   ! Do some computations using data from previous time step

   ! Bring to physical space to track decay rate
   call decay(maxVyx, maxTyx, VM, TM, PhysChebV, PhysChebT, tstep)
   write(8000, fmt=3000) time, maxVyx, maxTyx

   ! Now move on to next time step
   if (time == t_final) then
      exit
   end if

   dt_final = t_final - time

   if (dt_final <= dt) then
      time = t_final
      tstep = tstep + 1
   else
      time = time + dt
      tstep = tstep + 1
   end if

   call initrk(K1hV,K1hT, GPVM,GPD2VM,GPTM,PVM,PTM)
   call stage1(K1V,K1T,K2hV,K2hT, &
               dt,nu,kappa,GPVM,GPD2VM,GPD4VM,GPTM,PVM,PTM,K1hV,K1hT)
   call stage2(K2V,K2T,K3hV,K3hT, & 
               K1V,K1T,dt,nu,kappa,GPVM,GPD2VM,GPD4VM,GPTM,PVM,PTM,K1hV,K2hV,K1hT,K2hT)
   call stage3(K3V,K3T,K4hV,K4hT, &
               K1V,K1T,K2V,K2T,dt,nu,kappa,GPVM,GPD2VM,GPD4VM,GPTM,PVM,PTM,K1hV,K2hV,K3hV,K1hT,K2hT,K3hT)

   ! Update solution
   Aml = Aml + dt*(c1*(K1V+K2hV) + c2*(K2V+K3hV) + c3*(K3V+K4hV))
   Bml = Bml + dt*(c1*(K1T+K2hT) + c2*(K2T+K3hT) + c3*(K3T+K4hT))

end do

close(unit=8000)

2000 format(E25.16E3, E25.16E3          )
3000 format(E25.16E3, E25.16E3, E25.16E3)

end subroutine imex_rk

subroutine initrk(K1hV,K1hT, GPVM,GPD2VM,GPTM,PVM,PTM)

  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  !
  ! Set up IMEX-RK 3-4-3 method by computing K1hat
  !
  ! INPUT:
  !    GPVM   : Projected vertical velocity modes
  !    GPD2VM : Projected 2nd derivative vertical velocity modes
  !    GPTM   : Projected temperature modes for vertical velocity equation
  !    PVM    : Projected vertical velocity modes for temperature equation
  !    PTM    : Projected temperature modes for temperature equation
  !
  ! OUTPUT:
  !    K1hV: First explicit stage for vertical velocity equation
  !    K1hT: First explicit stage for temperature equation
  !
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

real(dp)   , dimension(:,:),              intent(in)  :: GPVM, GPD2VM, GPTM, PVM, PTM
complex(dp), allocatable, dimension(:,:), intent(out) :: K1hV, K1hT
real(dp)   , allocatable, dimension(:,:)              :: lhs
real(dp)   , allocatable, dimension(:,:)              :: rhs
real(dp)   , allocatable, dimension(:,:)              :: ups
real(dp)   , allocatable, dimension(:)                :: Kre, Kim
real(dp)                                              :: wave, wave2
integer                                               :: i
integer                                               :: NC, NF2

! Some LAPACK and BLAS variables
integer,     allocatable, dimension(:)                :: ipiv
integer                                               :: info
real(dp), parameter                                   :: scale1=1.0_dp, scale2=0.0_dp
integer , parameter                                   :: incx=1, incy=1

NC  = size(Bml,1)
NF2 = size(Bml,2)

allocate(K1hV(NC,NF2), K1hT(NC,NF2), stat=alloc_err)
allocate(lhs(NC,NC)                , stat=alloc_err)
allocate(rhs(NC,2)                 , stat=alloc_err)
allocate(ups(NC,NC)                , stat=alloc_err)
allocate(Kre(NC), Kim(NC)          , stat=alloc_err)
allocate(ipiv(NC)                  , stat=alloc_err)

K1hV    = cmplx(0.0_dp, 0.0_dp)
K1hT    = cmplx(0.0_dp, 0.0_dp)
lhs     = 0.0_dp
rhs     = 0.0_dp
ups     = 0.0_dp
Kre     = 0.0_dp
Kim     = 0.0_dp
ipiv    = 0

do i = 1,NF2
   ! Vertical velocity equation
   wave  = kx(i)
   wave2 = wave**2.0_dp

   ! Upsilon matrix pre-multiplying solution in vertical velocity equation
   ! time derivative term
   ups = -wave2*GPVM + GPD2VM

   ! Compute upsilon*rhs
   call dgemv('n', NC, NC, scale1, GPTM, NC, real (Bml(:,i)), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, GPTM, NC, aimag(Bml(:,i)), incx, scale2, Kim, incy)
   rhs(:,1) = Kre
   rhs(:,2) = Kim
   
   ! Get K1hV
   call dgesv(NC, 2, ups, NC, ipiv, rhs, NC, info)
   K1hV(:,i) = -wave2*cmplx(rhs(:,1), rhs(:,2))

   ! Compute PVM*rhs
   call dgemv('n', NC, NC, scale1, PVM, NC, real (Aml(:,i)), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, PVM, NC, aimag(Aml(:,i)), incx, scale2, Kim, incy)
   rhs(:,1) = Kre
   rhs(:,2) = Kim

   ! Get K1hT
   ups = PTM ! Use ups as a temp array for now
   call dgesv(NC, 2, ups, NC, ipiv, rhs, NC, info)
   K1hT(:,i) = cmplx(rhs(:,1), rhs(:,2))
end do

end subroutine initrk

subroutine stage1(K1V,K1T,K2hV,K2hT, &
                 dt,nu,kappa,GPVM,GPD2VM,GPD4VM,GPTM,PVM,PTM,K1hV,K1hT)

  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  ! Performs first stage of IMEX-RK 3-4-3 method
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  !
  ! INPUT:
  !    dt     : time-step size
  !    nu     : kinematic viscosity
  !    kappa  : thermal diffusivity
  !    GPVM   : Projected vertical velocity modes
  !    GPD2VM : Projected 2nd derivative vertical velocity modes
  !    GPD4VM : Projected 4th derivative vertical velocity modes
  !    GPTM   : Projected temperature modes for vertical velocity equation
  !    PVM    : Projected vertical velocity modes for temperature equation
  !    PTM    : Projected temperature modes
  !    K1hV   : First explicit stage for vertical velocity equation
  !    K1hT   : First explicit stage for temperature equation
  !
  ! OUTPUT:
  !    K1V : First stage for vertical velocity equation
  !    K1T : First stage for vertical temperature equation
  !    K2hV: Second explicit stage for vertical velocity equation
  !    K2hT: Second explicit stage for temperature equation
  !
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

implicit none

real(dp)   ,                              intent(in)  :: dt, nu, kappa
real(dp)   , dimension(:,:),              intent(in)  :: GPVM, GPD2VM, GPD4VM
real(dp)   , dimension(:,:),              intent(in)  :: GPTM, PVM, PTM
complex(dp), allocatable, dimension(:,:), intent(in)  :: K1hV, K1hT
complex(dp), allocatable, dimension(:,:), intent(out) :: K1V, K1T, K2hV, K2hT
real(dp)   , allocatable, dimension(:,:)              :: lhs
real(dp)   , allocatable, dimension(:,:)              :: rhs
real(dp)   , allocatable, dimension(:,:)              :: tempmat, ups
real(dp)   , allocatable, dimension(:,:)              :: Kmat
real(dp)   , allocatable, dimension(:)                :: Kre, Kim
complex(dp), allocatable, dimension(:)                :: temp, Aml1
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
allocate(K2hV(NC,NF2),K2hT(NC,NF2) , stat=alloc_err)
allocate(lhs(NC,NC)                , stat=alloc_err)
allocate(rhs(NC,2)                 , stat=alloc_err)
allocate(tempmat(NC,NC), ups(NC,NC), stat=alloc_err)
allocate(Kmat(NC,NC)               , stat=alloc_err)
allocate(Kre(NC), Kim(NC)          , stat=alloc_err)
allocate(ipiv(NC)                  , stat=alloc_err)
allocate(temp(NC), Aml1(NC)        , stat=alloc_err)

K1T     = cmplx(0.0_dp, 0.0_dp)
K1V     = cmplx(0.0_dp, 0.0_dp)
K2hV    = cmplx(0.0_dp, 0.0_dp)
K2hT    = cmplx(0.0_dp, 0.0_dp)
lhs     = 0.0_dp
rhs     = 0.0_dp
tempmat = 0.0_dp
ups     = 0.0_dp
Kmat    = 0.0_dp
Kre     = 0.0_dp
Kim     = 0.0_dp
ipiv    = 0
temp    = cmplx(0.0_dp, 0.0_dp)
Aml1    = cmplx(0.0_dp, 0.0_dp)

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
   temp     = Aml(:,i) + dt*dh21*K1hV(:,i)
   rhs(:,1) = real (temp)
   rhs(:,2) = aimag(temp)

   ! Compute upsilon*rhs
   call dgemv('n', NC, NC, scale1, ups, NC, rhs(:,1), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, ups, NC, rhs(:,2), incx, scale2, Kim, incy)
   rhs(:,1) = Kre
   rhs(:,2) = Kim

   ! Now form lhs matrix (lhs*a1 = ups*r)
   lhs = ups - dt*d11*nu*tempmat
   
   ! Get Aml_1
   call dgesv(NC, 2, lhs, NC, ipiv, rhs, NC, info)
   Aml1 = cmplx(rhs(:,1), rhs(:,2))

   ! Get K1 for vertical velocity equation (avoids ever forming inv(upsilon))
   call dgemv('n', NC, NC, scale1, tempmat, NC, rhs(:,1), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, tempmat, NC, rhs(:,2), incx, scale2, Kim, incy)
   rhs(:,1) = Kre
   rhs(:,2) = Kim
   tempmat  = ups
   call dgesv(NC, 2, tempmat, NC, ipiv, rhs, NC, info)

   K1V(:,i) = nu*cmplx(rhs(:,1), rhs(:,2))

   ! Temperature equation
   lhs  = (1.0_dp + dt*d11*kappa*wave2)*PTM - dt*d11*kappa*eye
   temp = Bml(:,i) + dt*dh21*K1hT(:,i)
   rhs(:,1) = real (temp)
   rhs(:,2) = aimag(temp)
   ! Get Bml_1 (partially)
   ! Here we just form inv(M)*rhs but Bml_1 = PTM*inv(M)*rhs
   call dgesv(NC, 2, lhs, NC, ipiv, rhs, NC, info)
   ! Get K1 for temperature equation
   Kmat = kappa*(-wave2*PTM + eye)
   call dgemv('n', NC, NC, scale1, Kmat, NC, rhs(:,1), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, Kmat, NC, rhs(:,2), incx, scale2, Kim, incy)
   K1T(:,i) = cmplx(Kre, Kim)

   ! Form explicit terms

   ! First, complete formation of Bml_1
   call dgemv('n', NC, NC, scale1, PTM, NC, rhs(:,1), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, PTM, NC, rhs(:,2), incx, scale2, Kim, incy)
    
   ! Take GPTM * Bml_1
   rhs(:,1) = Kre
   rhs(:,2) = Kim
   call dgemv('n', NC, NC, scale1, GPTM, NC, rhs(:,1), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, GPTM, NC, rhs(:,2), incx, scale2, Kim, incy)

   rhs(:,1) = -wave2*Kre
   rhs(:,2) = -wave2*Kim
   ! Get K2hV
   tempmat = ups
   call dgesv(NC, 2, tempmat, NC, ipiv, rhs, NC, info)
   K2hV(:,i) = cmplx(rhs(:,1), rhs(:,2))

   ! K2hT
   call dgemv('n', NC, NC, scale1, PVM, NC, real (Aml1), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, PVM, NC, aimag(Aml1), incx, scale2, Kim, incy)
   rhs(:,1) = Kre
   rhs(:,2) = Kim
   tempmat = PTM
   call dgesv(NC, 2, tempmat, NC, ipiv, rhs, NC, info)
   K2hT(:,i) = cmplx(rhs(:,1), rhs(:,2))
end do

!call write_out_mat(real (K1V), "ReK1V")
!call write_out_mat(aimag(K1V), "ImK1V")
!call write_out_mat(real (K1T), "ReK1T")
!call write_out_mat(aimag(K1T), "ImK1T")
!call write_out_mat(real (K2hV), "ReK2hV")
!call write_out_mat(aimag(K2hV), "ImK2hV")
!call write_out_mat(real (K2hT), "ReK2hT")
!call write_out_mat(aimag(K2hT), "ImK2hT")
!call write_out_mat(real (Bml), "ReBml")
!call write_out_mat(aimag(Bml), "ImBml")
!
!stop

end subroutine stage1

subroutine stage2(K2V,K2T,K3hV,K3hT, &
               K1V,K1T,dt,nu,kappa,GPVM,GPD2VM,GPD4VM,GPTM,PVM,PTM,K1hV,K2hV,K1hT,K2hT)

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
  !    GPVM   : Projected vertical velocity modes
  !    GPD2VM : Projected 2nd derivative vertical velocity modes
  !    GPD4VM : Projected 4th derivative vertical velocity modes
  !    GPTM   : Projected temperature modes for vertical velocity equation
  !    PVM    : Projected vertical velocity modes for temperature equation
  !    PTM    : Projected temperature modes
  !    K1hV   : First explicit stage for vertical velocity equation
  !    K2hV   : Second explicit stage for vertical velocity equation
  !    K1hT   : First explicit stage for temperature equation
  !    K2hT   : Second explicit stage for temperature equation
  !
  ! OUTPUT:
  !    K2V : Second stage for vertical velocity equation
  !    K2T : Second stage for vertical temperature equation
  !    K3hV: Third explicit stage for vertical velocity equation
  !    K3hT: Third explicit stage for temperature equation
  !
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

implicit none

real(dp)   ,                              intent(in)  :: dt, nu, kappa
real(dp)   , dimension(:,:),              intent(in)  :: GPVM, GPD2VM, GPD4VM
real(dp)   , dimension(:,:),              intent(in)  :: GPTM, PVM, PTM
complex(dp), allocatable, dimension(:,:), intent(in)  :: K1V, K1T
complex(dp), allocatable, dimension(:,:), intent(in)  :: K1hV, K2hV
complex(dp), allocatable, dimension(:,:), intent(in)  :: K1hT, K2hT
complex(dp), allocatable, dimension(:,:), intent(out) :: K2V, K2T, K3hV, K3hT
real(dp)   , allocatable, dimension(:,:)              :: lhs
real(dp)   , allocatable, dimension(:,:)              :: rhs
real(dp)   , allocatable, dimension(:,:)              :: tempmat, ups
real(dp)   , allocatable, dimension(:,:)              :: Kmat
real(dp)   , allocatable, dimension(:)                :: Kre, Kim
complex(dp), allocatable, dimension(:)                :: temp, Aml2
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
allocate(K3hV(NC,NF2),K3hT(NC,NF2) , stat=alloc_err)
allocate(lhs(NC,NC)                , stat=alloc_err)
allocate(rhs(NC,2)                 , stat=alloc_err)
allocate(tempmat(NC,NC), ups(NC,NC), stat=alloc_err)
allocate(Kmat(NC,NC)               , stat=alloc_err)
allocate(Kre(NC), Kim(NC)          , stat=alloc_err)
allocate(ipiv(NC)                  , stat=alloc_err)
allocate(temp(NC), Aml2(NC)        , stat=alloc_err)

K2T     = cmplx(0.0_dp, 0.0_dp)
K2V     = cmplx(0.0_dp, 0.0_dp)
K3hV    = cmplx(0.0_dp, 0.0_dp)
K3hT    = cmplx(0.0_dp, 0.0_dp)
lhs     = 0.0_dp
rhs     = 0.0_dp
tempmat = 0.0_dp
ups     = 0.0_dp
Kmat    = 0.0_dp
Kre     = 0.0_dp
Kim     = 0.0_dp
ipiv    = 0
temp    = cmplx(0.0_dp, 0.0_dp)
Aml2    = cmplx(0.0_dp, 0.0_dp)

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

   ! RHS for second stage
   temp     = Aml(:,i) + dt*(d21*K1V(:,i) + dh31*K1hV(:,i) + dh32*K2hV(:,i))
   rhs(:,1) = real (temp)
   rhs(:,2) = aimag(temp)

   ! Compute upsilon*rhs
   call dgemv('n', NC, NC, scale1, ups, NC, rhs(:,1), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, ups, NC, rhs(:,2), incx, scale2, Kim, incy)
   rhs(:,1) = Kre
   rhs(:,2) = Kim

   ! Now form lhs matrix (lhs*a1 = ups*a_{n-1})
   lhs = ups - dt*d22*nu*tempmat
   
   ! Get Aml_2
   call dgesv(NC, 2, lhs, NC, ipiv, rhs, NC, info)
   Aml2 = cmplx(rhs(:,1), rhs(:,2))

   ! Get K2 for vertical velocity equation (avoids ever forming inv(upsilon))
   call dgemv('n', NC, NC, scale1, tempmat, NC, rhs(:,1), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, tempmat, NC, rhs(:,2), incx, scale2, Kim, incy)
   rhs(:,1) = Kre
   rhs(:,2) = Kim
   tempmat  = ups
   call dgesv(NC, 2, tempmat, NC, ipiv, rhs, NC, info)

   K2V(:,i) = nu*cmplx(rhs(:,1), rhs(:,2))

   ! Temperature equation
   lhs  = (1.0_dp + dt*d22*kappa*wave2)*PTM - dt*d22*kappa*eye
   temp = Bml(:,i) + dt*(d21*K1T(:,i) + dh31*K1hT(:,i) + dh32*K2hT(:,i))
   rhs(:,1) = real (temp)
   rhs(:,2) = aimag(temp)
   ! Get Bml_2 (partially)
   ! Here we just form inv(M)*rhs but Bml_2 = PTM*inv(M)*rhs
   call dgesv(NC, 2, lhs, NC, ipiv, rhs, NC, info)
   ! Get K2 for temperature equation
   Kmat = kappa*(-wave2*PTM + eye)
   call dgemv('n', NC, NC, scale1, Kmat, NC, rhs(:,1), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, Kmat, NC, rhs(:,2), incx, scale2, Kim, incy)
   K2T(:,i) = cmplx(Kre, Kim)

   ! Form explicit terms

   ! First, complete formation of Bml_2
   call dgemv('n', NC, NC, scale1, PTM, NC, rhs(:,1), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, PTM, NC, rhs(:,2), incx, scale2, Kim, incy)
    
   ! Take GPTM * Bml_2
   rhs(:,1) = Kre
   rhs(:,2) = Kim
   call dgemv('n', NC, NC, scale1, GPTM, NC, rhs(:,1), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, GPTM, NC, rhs(:,2), incx, scale2, Kim, incy)

   rhs(:,1) = -wave2*Kre
   rhs(:,2) = -wave2*Kim
   ! Get K3hV
   tempmat = ups
   call dgesv(NC, 2, tempmat, NC, ipiv, rhs, NC, info)
   K3hV(:,i) = cmplx(rhs(:,1), rhs(:,2))

   ! K3hT
   call dgemv('n', NC, NC, scale1, PVM, NC, real (Aml2), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, PVM, NC, aimag(Aml2), incx, scale2, Kim, incy)
   rhs(:,1) = Kre
   rhs(:,2) = Kim
   tempmat = PTM
   call dgesv(NC, 2, tempmat, NC, ipiv, rhs, NC, info)
   K3hT(:,i) = cmplx(rhs(:,1), rhs(:,2))
end do

end subroutine stage2

subroutine stage3(K3V,K3T,K4hV,K4hT, &
               K1V,K1T,K2V,K2T,dt,nu,kappa,GPVM,GPD2VM,GPD4VM,GPTM,PVM,PTM,K1hV,K2hV,K3hV,K1hT,K2hT,K3hT)

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
  !    GPVM   : Projected vertical velocity modes
  !    GPD2VM : Projected 2nd derivative vertical velocity modes
  !    GPD4VM : Projected 4th derivative vertical velocity modes
  !    GPTM   : Projected temperature modes for vertical velocity equation
  !    PVM    : Projected vertical velocity modes for temperature equation
  !    PTM    : Projected temperature modes
  !    K1hV   : First explicit stage for vertical velocity equation
  !    K2hV   : Second explicit stage for vertical velocity equation
  !    K3hV   : Third explicit stage for vertical velocity equation
  !    K1hT   : First explicit stage for temperature equation
  !    K2hT   : Second explicit stage for temperature equation
  !    K3hT   : Third explicit stage for temperature equation
  !
  ! OUTPUT:
  !    K3V : Third stage for vertical velocity equation
  !    K3T : Third stage for vertical temperature equation
  !    K4hV: Fourth explicit stage for vertical velocity equation
  !    K4hT: Fourth explicit stage for temperature equation
  !
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

implicit none

real(dp)   ,                              intent(in)  :: dt, nu, kappa
real(dp)   , dimension(:,:),              intent(in)  :: GPVM, GPD2VM, GPD4VM
real(dp)   , dimension(:,:),              intent(in)  :: GPTM, PVM, PTM
complex(dp), allocatable, dimension(:,:), intent(in)  :: K1V, K1T, K2V, K2T
complex(dp), allocatable, dimension(:,:), intent(in)  :: K1hV, K2hV, K3hV
complex(dp), allocatable, dimension(:,:), intent(in)  :: K1hT, K2hT, K3hT
complex(dp), allocatable, dimension(:,:), intent(out) :: K3V, K3T, K4hV, K4hT
real(dp)   , allocatable, dimension(:,:)              :: lhs
real(dp)   , allocatable, dimension(:,:)              :: rhs
real(dp)   , allocatable, dimension(:,:)              :: tempmat, ups
real(dp)   , allocatable, dimension(:,:)              :: Kmat
real(dp)   , allocatable, dimension(:)                :: Kre, Kim
complex(dp), allocatable, dimension(:)                :: temp, Aml3
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
allocate(K4hV(NC,NF2),K4hT(NC,NF2) , stat=alloc_err)
allocate(lhs(NC,NC)                , stat=alloc_err)
allocate(rhs(NC,2)                 , stat=alloc_err)
allocate(tempmat(NC,NC), ups(NC,NC), stat=alloc_err)
allocate(Kmat(NC,NC)               , stat=alloc_err)
allocate(Kre(NC), Kim(NC)          , stat=alloc_err)
allocate(ipiv(NC)                  , stat=alloc_err)
allocate(temp(NC), Aml3(NC)        , stat=alloc_err)

K3T     = cmplx(0.0_dp, 0.0_dp)
K3V     = cmplx(0.0_dp, 0.0_dp)
K4hV    = cmplx(0.0_dp, 0.0_dp)
K4hT    = cmplx(0.0_dp, 0.0_dp)
lhs     = 0.0_dp
rhs     = 0.0_dp
tempmat = 0.0_dp
ups     = 0.0_dp
Kmat    = 0.0_dp
Kre     = 0.0_dp
Kim     = 0.0_dp
ipiv    = 0
temp    = cmplx(0.0_dp, 0.0_dp)
Aml3    = cmplx(0.0_dp, 0.0_dp)

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
   temp     = Aml(:,i) + dt*(d31*K1V(:,i) + d32*K2V(:,i) + dh41*K1hV(:,i) + dh42*K2hV(:,i) + dh43*K3hV(:,i))
   rhs(:,1) = real (temp)
   rhs(:,2) = aimag(temp)

   ! Compute upsilon*rhs
   call dgemv('n', NC, NC, scale1, ups, NC, rhs(:,1), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, ups, NC, rhs(:,2), incx, scale2, Kim, incy)
   rhs(:,1) = Kre
   rhs(:,2) = Kim

   ! Now form lhs matrix (lhs*a1 = ups*r)
   lhs = ups - dt*d33*nu*tempmat
   
   ! Get Aml_3
   call dgesv(NC, 2, lhs, NC, ipiv, rhs, NC, info)
   Aml3 = cmplx(rhs(:,1), rhs(:,2))

   ! Get K3 for vertical velocity equation (avoids ever forming inv(upsilon))
   call dgemv('n', NC, NC, scale1, tempmat, NC, rhs(:,1), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, tempmat, NC, rhs(:,2), incx, scale2, Kim, incy)
   rhs(:,1) = Kre
   rhs(:,2) = Kim
   tempmat  = ups
   call dgesv(NC, 2, tempmat, NC, ipiv, rhs, NC, info)

   K3V(:,i) = nu*cmplx(rhs(:,1), rhs(:,2))

   ! Temperature equation
   lhs  = (1.0_dp + dt*d33*kappa*wave2)*PTM - dt*d33*kappa*eye
   temp = Bml(:,i) + dt*(d31*K1T(:,i) + d32*K2T(:,i) + dh41*K1hT(:,i) + dh42*K2hT(:,i) + dh43*K3hT(:,i))
   rhs(:,1) = real (temp)
   rhs(:,2) = aimag(temp)
   ! Get Bml_3 (partially)
   ! Here we just form inv(M)*rhs but Bml_3 = PTM*inv(M)*rhs
   call dgesv(NC, 2, lhs, NC, ipiv, rhs, NC, info)
   ! Get K3 for temperature equation
   Kmat = kappa*(-wave2*PTM + eye)
   call dgemv('n', NC, NC, scale1, Kmat, NC, rhs(:,1), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, Kmat, NC, rhs(:,2), incx, scale2, Kim, incy)
   K3T(:,i) = cmplx(Kre, Kim)

   ! Form explicit terms

   ! First, complete formation of Bml_3
   call dgemv('n', NC, NC, scale1, PTM, NC, rhs(:,1), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, PTM, NC, rhs(:,2), incx, scale2, Kim, incy)
    
   ! Take GPTM * Bml_3
   rhs(:,1) = Kre
   rhs(:,2) = Kim
   call dgemv('n', NC, NC, scale1, GPTM, NC, rhs(:,1), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, GPTM, NC, rhs(:,2), incx, scale2, Kim, incy)

   rhs(:,1) = -wave2*Kre
   rhs(:,2) = -wave2*Kim
   ! Get K4hV
   tempmat = ups
   call dgesv(NC, 2, tempmat, NC, ipiv, rhs, NC, info)
   K4hV(:,i) = cmplx(rhs(:,1), rhs(:,2))

   ! K4hT
   call dgemv('n', NC, NC, scale1, PVM, NC, real (Aml3), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, PVM, NC, aimag(Aml3), incx, scale2, Kim, incy)
   rhs(:,1) = Kre
   rhs(:,2) = Kim
   tempmat = PTM
   call dgesv(NC, 2, tempmat, NC, ipiv, rhs, NC, info)
   K4hT(:,i) = cmplx(rhs(:,1), rhs(:,2))
end do

end subroutine stage3

subroutine decay(maxVyx, maxTyx, VM, TM, PhysChebV, PhysChebT, iter)

implicit none

integer,                  intent(in)  :: iter
real(dp), dimension(:,:), intent(in)  :: VM, TM
real(dp), dimension(:,:), intent(in)  :: PhysChebV, PhysChebT
real(dp),                 intent(out) :: maxVyx, maxTyx
integer                               :: NP, NF, NC
character(10)                         :: citer
! Some LAPACK and BLAS parameters
real(dp), parameter :: scale1=1.0_dp, scale2=0.0_dp

character(9)                          :: fiter

write(citer, "(I10)") 1000000000 + iter
fiter = citer(2:10)

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

! Write out the fields
!call write_out_mat(Vyx, "Vyx"//fiter)
!call write_out_mat(Vyx, "Tyx"//fiter)

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
                          GPVM,GPTM,PVM,PTM,GPD2VM,GPD4VM)

integer                 , intent(in)     :: NC, NF
real(dp),                 intent(in)     :: dt, t_final, kappa, nu
real(dp), dimension(:,:), intent(in)     :: PhysChebV, PhysChebT
real(dp), dimension(:,:), intent(in)     :: VM, TM, PVM, PTM, GPTM
real(dp), dimension(:,:), intent(in)     :: GPVM, GPD2VM, GPD4VM
real(dp)                                 :: dt_final, time
real(dp)                                 :: wave, wave2, wave4
real(dp)                                 :: maxVyx, maxTyx
real(dp), allocatable, dimension(:,:)    :: lhs, rhs
real(dp), allocatable, dimension(:,:)    :: upsilon
complex(dp), allocatable, dimension(:)   :: temprhs
real(dp), parameter                      :: scale1=1.0_dp, scale2=0.0_dp
integer,  parameter                      :: incx=1, incy=1
integer,  allocatable, dimension(:)      :: ipiv
integer                                  :: info
integer                                  :: i, j, k, tstep

allocate(upsilon(NC,NC), stat=alloc_err)
allocate(lhs(NC,NC)    , stat=alloc_err)
allocate(rhs(NC,2)     , stat=alloc_err)
allocate(ipiv(NC)      , stat=alloc_err)
allocate(temprhs(NC)   , stat=alloc_err)
upsilon = 0.0_dp
lhs     = 0.0_dp
rhs     = 0.0_dp
ipiv    = 0
temprhs = cmplx(0.0_dp, 0.0_dp)

! Initial time
time   = 0.0_dp
tstep = 0

call decay(maxVyx, maxTyx, VM, TM, PhysChebV, PhysChebT, tstep)

open(unit=8000, file="maxval.txt", action="write", status="unknown", position="append")
write(8000, fmt=3000) time, maxVyx, maxTyx

do ! while time < t_final

   dt_final = t_final - time

   if (dt_final <= dt) then
      time = t_final
      tstep = tstep + 1
   else
      time = time + dt
      tstep = tstep + 1
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
      temprhs = cmplx(rhs(:,1), rhs(:,2))
      call dgemv('n', NC, NC, scale1, GPTM, NC, real (Bml(:,i)), incx, scale2, rhs(:,1), incy)
      call dgemv('n', NC, NC, scale1, GPTM, NC, aimag(Bml(:,i)), incx, scale2, rhs(:,2), incy)
      rhs(:,1) = real(temprhs)  - wave2*dt*rhs(:,1)
      rhs(:,2) = aimag(temprhs) - wave2*dt*rhs(:,2)
      ! Solve for the solution at next time step
      call dgesv(NC, 2, lhs, NC, ipiv, rhs, NC, info)
      ! Update solution
      Aml(:,i) = cmplx(rhs(:,1), rhs(:,2))

      ! Temperature Equation
      lhs = (1.0_dp + kappa*dt*wave2)*PTM - kappa*dt*eye
      ! Compute the RHS
      call dgemv('n', NC, NC, scale1, PTM, NC, real (Bml(:,i)), incx, scale2, rhs(:,1), incy)
      call dgemv('n', NC, NC, scale1, PTM, NC, aimag(Bml(:,i)), incx, scale2, rhs(:,2), incy)
      temprhs = cmplx(rhs(:,1), rhs(:,2))
      call dgemv('n', NC, NC, scale1, PVM, NC, real (Aml(:,i)), incx, scale2, rhs(:,1), incy)
      call dgemv('n', NC, NC, scale1, PVM, NC, aimag(Aml(:,i)), incx, scale2, rhs(:,2), incy)
      rhs(:,1) = real(temprhs)  + dt*rhs(:,1)
      rhs(:,2) = aimag(temprhs) + dt*rhs(:,2)
      ! Solve for solution at next time
      call dgesv(NC, 2, lhs, NC, ipiv, rhs, NC, info)
      ! Update solution
      Bml(:,i) = cmplx(rhs(:,1), rhs(:,2))
   end do

   call decay(maxVyx, maxTyx, VM, TM, PhysChebV, PhysChebT, tstep)
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
