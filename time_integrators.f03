module time_integrators

use types,  only: dp, CI, pi
use global, only: set_imex_params, fft_utils,      &
                  alloc_err, eye, kx,              &
                  Aml, Bml, Uml, Amx, Bmx, Umx,    &
                  Vyx, Tyx, Uyx, fft1_ml, fft1_mx, &
                  fft1_yl, fft1_yx, pU, ipU, pV,   &
                  ipV, pT, ipT, pf1, ipf1, pfy,    &
                  ipfy, y

use fftw
use write_pack
implicit none

private
public imex_rk, backward_euler, nonlinear_terms, update_u, probes

real(dp) :: c1, c2, c3
real(dp) :: d11, d12, d13
real(dp) :: d21, d22, d23
real(dp) :: d31, d32, d33
real(dp) :: dh11, dh12, dh13, dh14
real(dp) :: dh21, dh22, dh23, dh24
real(dp) :: dh31, dh32, dh33, dh34
real(dp) :: dh41, dh42, dh43, dh44

contains

subroutine imex_rk(NC, NF, dt, t_final, nu, kappa,      &
                   PVEL, Pmj,VM,TM,DVM,DTM,D2VM,D3VM,   &
                   GPVM,GPTM,PVM,PDVM,PTM,GPD2VM,GPD4VM)

integer                 , intent(in)     :: NC, NF
real(dp),                 intent(in)     :: dt, t_final, kappa, nu
real(dp), dimension(:,:), intent(in)     :: PVEL, Pmj
real(dp), dimension(:,:), intent(in)     :: VM, TM, DVM, DTM, D2VM, D3VM
real(dp), dimension(:,:), intent(in)     :: PVM, PDVM, PTM, GPTM
real(dp), dimension(:,:), intent(in)     :: GPVM, GPD2VM, GPD4VM
complex(dp), allocatable, dimension(:,:) :: K1V, K2V, K3V
complex(dp), allocatable, dimension(:,:) :: K1T, K2T, K3T
complex(dp), allocatable, dimension(:,:) :: K1hV, K2hV, K3hV, K4hV
complex(dp), allocatable, dimension(:,:) :: K1hT, K2hT, K3hT, K4hT
real(dp)                                 :: dt_final, time
real(dp)                                 :: Vprobe, Tprobe
integer                                  :: tstep
integer                                  :: NF2
integer, parameter                       :: interval=500 ! How often to write fields to file
logical                                  :: iprobe=.true. ! Get time trace of fields

! LAPACK and BLAS parameters
real(dp), parameter                      :: scale1=1.0_dp, scale2=0.0_dp
integer , parameter                      :: incx=1, incy=1

call set_imex_params(c1,c2,c3, d11,d21,d31, d12,d22,d32, d13,d23,d33, &
                     dh11,dh12,dh13,dh14, dh21,dh22,dh23,dh24,        &
                     dh31,dh32,dh33,dh34, dh41,dh42,dh43,dh44)

NF2 = NF/2 + 1

allocate(K1T(NC,NF2), K1V(NC,NF2), stat=alloc_err)
allocate(K2T(NC,NF2), K2V(NC,NF2), stat=alloc_err)
allocate(K3T(NC,NF2), K3V(NC,NF2), stat=alloc_err)
allocate(K1hT(NC,NF2), K1hV(NC,NF2), stat=alloc_err)
allocate(K2hT(NC,NF2), K2hV(NC,NF2), stat=alloc_err)
allocate(K3hT(NC,NF2), K3hV(NC,NF2), stat=alloc_err)
allocate(K4hT(NC,NF2), K4hV(NC,NF2), stat=alloc_err)

K1V = cmplx(0.0_dp, 0.0_dp, kind=dp)
K2V = cmplx(0.0_dp, 0.0_dp, kind=dp)
K3V = cmplx(0.0_dp, 0.0_dp, kind=dp)

K1T = cmplx(0.0_dp, 0.0_dp, kind=dp)
K2T = cmplx(0.0_dp, 0.0_dp, kind=dp)
K3T = cmplx(0.0_dp, 0.0_dp, kind=dp)

K1hV = cmplx(0.0_dp, 0.0_dp, kind=dp)
K2hV = cmplx(0.0_dp, 0.0_dp, kind=dp)
K3hV = cmplx(0.0_dp, 0.0_dp, kind=dp)
K4hV = cmplx(0.0_dp, 0.0_dp, kind=dp)

K1hT = cmplx(0.0_dp, 0.0_dp, kind=dp)
K2hT = cmplx(0.0_dp, 0.0_dp, kind=dp)
K3hT = cmplx(0.0_dp, 0.0_dp, kind=dp)
K4hT = cmplx(0.0_dp, 0.0_dp, kind=dp)

time  = 0.0_dp
tstep = 0

open(unit=1500, file="Tprobe.txt", action="write", status="unknown", position="append")

do ! while time < t_final

   ! Do some computations using data from previous time step
   call probes(Vprobe, Tprobe, VM, TM, tstep, interval, iprobe)
   write(1500, fmt=3000) time, Vprobe, Tprobe

   ! Now move on to next time step
   if (time == t_final) then
      exit
   end if

   dt_final = t_final - time

   if (dt_final <= dt) then
      time = t_final
      tstep = tstep + 1
      write(*,*) "t = ", time
   else
      time = time + dt
      tstep = tstep + 1
      write(*,*) "t = ", time, "dt = ", dt
   end if

   call initrk(K1hV,K1hT, PVEL,Pmj,GPVM,GPD2VM,GPTM,PVM,VM,TM,DVM,DTM,D2VM,D3VM,PTM)
   call stage1(K1V,K1T,K2hV,K2hT,                        &
               dt,nu,kappa,PVEL,Pmj,VM,TM,DVM,DTM,D2VM,D3VM, &
               GPVM,GPD2VM,GPD4VM,GPTM,PVM,PDVM,PTM,     &
               K1hV,K1hT)
   call stage2(K2V,K2T,K3hV,K3hT,                                & 
               K1V,K1T,dt,nu,kappa,PVEL,Pmj,VM,TM,DVM,DTM,D2VM,D3VM, &
               GPVM,GPD2VM,GPD4VM,GPTM,PVM,PDVM,PTM,             &
               K1hV,K2hV,K1hT,K2hT)
   call stage3(K3V,K3T,K4hV,K4hT,                                        &
               K1V,K1T,K2V,K2T,dt,nu,kappa,PVEL,Pmj,VM,TM,DVM,DTM,D2VM,D3VM, &
               GPVM,GPD2VM,GPD4VM,GPTM,PVM,PDVM,PTM,                     &
               K1hV,K2hV,K3hV,K1hT,K2hT,K3hT)

   ! Update solution
   Aml = Aml + dt*(c1*(K1V+K2hV) + c2*(K2V+K3hV) + c3*(K3V+K4hV))
   Bml = Bml + dt*(c1*(K1T+K2hT) + c2*(K2T+K3hT) + c3*(K3T+K4hT))

   ! Update horizontal velocity
   call update_u(Pmj,DVM)

   !Update dt
   call update_dt(VM,TM,dt)
end do

close(unit=1500)

!2000 format(E25.16E3, E25.16E3          )
3000 format(E25.16E3, E25.16E3, E25.16E3)

end subroutine imex_rk

subroutine initrk(K1hV,K1hT, PVEL,Pmj,GPVM,GPD2VM,GPTM,PVM,VM,TM,DVM,DTM,D2VM,D3VM,PTM)

  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  !
  ! Set up IMEX-RK 3-4-3 method by computing K1hat
  !
  ! INPUT:
  !    Pmj    : Chebyshev projection matrix
  !    PVEL   : Cheb-projection matrix for vertical velocity equation
  !    GPVM   : Projected vertical velocity modes
  !    GPD2VM : Projected 2nd derivative vertical velocity modes
  !    GPTM   : Projected temperature modes for vertical velocity equation
  !    PVM    : Projected vertical velocity modes for temperature equation
  !    PTM    : Projected temperature modes for temperature equation
  !    VM     : Vertical velocity modes
  !    TM     : Temperature modes 
  !    DVM    : Modes for 1st derivative of vertical velocity
  !    DTM    : Modes for 1st derivative of temperature
  !    D2VM   : Modes for 2nd derivative of vertical velocity
  !    D3VM   : Modes for 3rd derivative of vertical velocity
  !
  ! OUTPUT:
  !    K1hV: First explicit stage for vertical velocity equation
  !    K1hT: First explicit stage for temperature equation
  !
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

real(dp)   , dimension(:,:),              intent(in)  :: GPVM, GPD2VM, GPTM
real(dp)   , dimension(:,:),              intent(in)  :: PVEL, Pmj, PVM
real(dp)   , dimension(:,:),              intent(in)  :: VM, TM, DVM, DTM, D2VM, D3VM, PTM
complex(dp), allocatable, dimension(:,:), intent(out) :: K1hV, K1hT
real(dp)   , allocatable, dimension(:,:)              :: lhs
real(dp)   , allocatable, dimension(:,:)              :: rhs
real(dp)   , allocatable, dimension(:,:)              :: ups
real(dp)   , allocatable, dimension(:)                :: Kre, Kim
complex(dp), allocatable, dimension(:)                :: temp
complex(dp), allocatable, dimension(:,:)              :: NLTml, NLVml
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

allocate(K1hV(NC,NF2), K1hT(NC,NF2)  , stat=alloc_err)
allocate(lhs(NC,NC)                  , stat=alloc_err)
allocate(rhs(NC,2)                   , stat=alloc_err)
allocate(ups(NC,NC)                  , stat=alloc_err)
allocate(Kre(NC), Kim(NC)            , stat=alloc_err)
allocate(ipiv(NC)                    , stat=alloc_err)
allocate(temp(NC)                    , stat=alloc_err)
allocate(NLTml(NC,NF2), NLVml(NC,NF2), stat=alloc_err)

K1hV  = cmplx(0.0_dp, 0.0_dp, kind=dp)
K1hT  = cmplx(0.0_dp, 0.0_dp, kind=dp)
lhs   = 0.0_dp
rhs   = 0.0_dp
ups   = 0.0_dp
Kre   = 0.0_dp
Kim   = 0.0_dp
temp  = cmplx(0.0_dp, 0.0_dp, kind=dp)
NLTml = cmplx(0.0_dp, 0.0_dp, kind=dp)
NLVml = cmplx(0.0_dp, 0.0_dp, kind=dp)
ipiv  = 0

call nonlinear_terms(NLVml,NLTml, PVEL,Pmj,VM,TM,DVM,DTM,D2VM,D3VM,Aml,Bml)

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

   rhs(:,1) = -(real (NLVml(:,i)) + wave2*Kre)
   rhs(:,2) = -(aimag(NLVml(:,i)) + wave2*Kim)
   
   ! Get K1hV
   call dgesv(NC, 2, ups, NC, ipiv, rhs, NC, info)
   K1hV(:,i) = cmplx(rhs(:,1), rhs(:,2), kind=dp)

   ! Compute PVM*rhs
   call dgemv('n', NC, NC, scale1, PVM, NC, real (Aml(:,i)), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, PVM, NC, aimag(Aml(:,i)), incx, scale2, Kim, incy)

   rhs(:,1) = -real (NLTml(:,i)) + Kre
   rhs(:,2) = -aimag(NLTml(:,i)) + Kim

   ! Get K1hT
   ups = PTM ! Use ups as a temp array for now
   call dgesv(NC, 2, ups, NC, ipiv, rhs, NC, info)
   K1hT(:,i) = cmplx(rhs(:,1), rhs(:,2), kind=dp)
end do

end subroutine initrk

subroutine nonlinear_terms(NLVml,NLTml, PVEL,Pmj,VM,TM,DVM,DTM,D2VM,D3VM,Aiml,Biml)

  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  ! Computation of nonlinear terms
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  !
  ! INPUT:
  !        PVEL: Chebyshev projection for VV equation
  !        Pmj : Chebyshev projection for temperature equation
  !        VM  : VV modes
  !        TM  : Temperature modes
  !        DVM : 1st derivative of VV modes
  !        DTM : 1st derivative of temperature modes
  !        D2VM: 2nd derivative of VV modes
  !        D3VM: 3rd derivative of VV modes
  !        AMli: Horizontal velocity modes at stage i
  !        BMli: Temperature modes at stage i
  !        UMli: Vertical velocity modes at stage i
  !
  ! OUTPUT:
  !        NLVml: Nonlinear term for VV equation
  !        NLTml: Nonlinear term for temperature equation
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

real(dp)   ,              dimension(:,:), intent(in)  :: PVEL, Pmj, TM
real(dp)   ,              dimension(:,:), intent(in)  :: VM, DTM, DVM, D2VM, D3VM 
complex(dp),              dimension(:,:), intent(in)  :: Aiml, Biml
complex(dp), allocatable, dimension(:,:), intent(out) :: NLVml, NLTml
complex(dp), allocatable, dimension(:,:)              :: DTml
complex(dp), allocatable, dimension(:,:)              :: lapV, lapU
real(dp)   , allocatable, dimension(:,:)              :: Uiyx, DTyx
real(dp)   , allocatable, dimension(:,:)              :: NLVyx, NLVmx, NLTyx, NLTmx
real(dp)   , allocatable, dimension(:,:)              :: lap
real(dp)                                              :: alpha, wave, wave2
integer                                               :: NC, NF2, NP, NF
integer                                               :: NF_cut, NC_cut
integer                                               :: i,j
real(dp)                                              :: l

! LAPACK and BLAS variables
real(dp), parameter :: scale1=1.0_dp, scale2=0.0_dp
integer,  parameter :: incx = 1, incy = 1

NC  = size(Bml,1)
NF2 = size(Bml,2)
NP  = size(Tyx,1)
NF  = size(Bmx,2)

allocate(NLVml(NC,NF2)     , stat=alloc_err)
allocate(NLTml(NC,NF2)     , stat=alloc_err)
allocate(DTml(NC,NF2)      , stat=alloc_err)
allocate(Uiyx(NP,NF)       , stat=alloc_err)
allocate(DTyx(NP,NF)       , stat=alloc_err)
allocate(NLVyx(NP,NF)      , stat=alloc_err)
allocate(NLVmx(NC,NF)      , stat=alloc_err)
allocate(NLTyx(NP,NF)      , stat=alloc_err)
allocate(NLTmx(NC,NF)      , stat=alloc_err)
allocate(lap(NP,NC)        , stat=alloc_err)
allocate(lapV(NP,NF2)      , stat=alloc_err)
allocate(lapU(NP,NF2)      , stat=alloc_err)

NLVml = cmplx(0.0_dp, 0.0_dp, kind=dp)
NLTml = cmplx(0.0_dp, 0.0_dp, kind=dp)
DTml  = cmplx(0.0_dp, 0.0_dp, kind=dp)
lapV  = cmplx(0.0_dp, 0.0_dp, kind=dp)
lapU  = cmplx(0.0_dp, 0.0_dp, kind=dp)
Uiyx  = 0.0_dp
DTyx  = 0.0_dp
NLVyx = 0.0_dp
NLVmx = 0.0_dp
NLTyx = 0.0_dp
NLTmx = 0.0_dp
lap   = 0.0_dp

do i = 1,NF2
   wave  = kx(i)
   wave2 = wave**2.0_dp
   DTml(:,i) = CI*wave*Biml(:,i)
   lap = -wave2*VM + D2VM
   lapV(:,i) = matmul(lap, Aiml(:,i))
   lap = -wave2*DVM + D3VM
   if (i /= 1) then
     lapU(:,i)    = CI*matmul(lap, Aiml(:,i))/wave
     fft1_ml(:,i) = CI*Aiml(:,i)/wave
   else
     lapU(:,i)    = cmplx(0.0_dp, 0.0_dp, kind=dp)
     fft1_ml(:,i) = cmplx(0.0_dp, 0.0_dp, kind=dp)
   end if
end do

! Bring fields to physical space
call fftw_execute_dft_c2r(ipf1, fft1_ml, fft1_mx) ! Horizontal velocity Fourier to physical
Uiyx = matmul(DVM, fft1_mx)

fft1_ml = DTml
call fftw_execute_dft_c2r(ipf1, fft1_ml, fft1_mx) ! Temperature deriv Fourier to physical
DTyx = matmul(TM, fft1_mx)

fft1_ml = Biml
call fftw_execute_dft_c2r(ipf1, fft1_ml, fft1_mx) ! Temperature Fourier to physical
Tyx = matmul(DTM, fft1_mx)

fft1_ml = Aiml
call fftw_execute_dft_c2r(ipf1, fft1_ml, fft1_mx) ! Vertical velocity Fourier to physical
Vyx = matmul(VM, fft1_mx)

fft1_yl = lapV 
call fftw_execute_dft_c2r(ipfy, fft1_yl, fft1_yx) ! Laplacian of VV from Fourier to physical

! Calculate nonlinear term
NLTyx = Uiyx*DTyx + Vyx*Tyx ! Temperature equation
NLVyx =  Uiyx*fft1_yx ! Vertical velocity equation (first part)

fft1_yl = lapU
call fftw_execute_dft_c2r(ipfy, fft1_yl, fft1_yx) ! Laplacian of HV from Fourier to physical

NLVyx = NLVyx - Vyx*fft1_yx ! Vertical velocity equation (second part)

! Bring nonlinear term to Chebyshev-Fourier space
NLTmx = matmul(Pmj, NLTyx)
NLVmx = matmul(PVEL, NLVyx)
fft1_mx = NLTmx
call fftw_execute_dft_r2c(pf1, fft1_mx, fft1_ml) ! Phys to Fourier
NLTml = fft1_ml / real(NF, dp)
fft1_mx = NLVmx
call fftw_execute_dft_r2c(pf1, fft1_mx, fft1_ml) ! Phys to Fourier
NLVml = fft1_ml / real(NF, dp)

! Dealias
NF_cut = floor(2.0_dp*NF/3.0_dp)/2
NC_cut = floor(2.0_dp*NC/3.0_dp)
alpha  = kx(2) ! Get alpha rather than passing it through everywhere
do i = 1,NF2
   wave = kx(i)
   l = abs(wave)/alpha
   do j = 1,NC
      if ((l >= NF_cut) .or. (j >= NC_cut)) then
         NLTml(j,i) = cmplx(0.0_dp, 0.0_dp, kind=dp)
         NLVml(j,i) = cmplx(0.0_dp, 0.0_dp, kind=dp)
      else
         NLVml(j,i) = CI*wave*NLVml(j,i)
      end if
   end do
end do

end subroutine nonlinear_terms

subroutine stage1(K1V,K1T,K2hV,K2hT,                           &
                 dt,nu,kappa,PVEL,Pmj,VM,TM,DVM,DTM,D2VM,D3VM, &
                 GPVM,GPD2VM,GPD4VM,GPTM,PVM,PDVM,PTM,         &
                 K1hV,K1hT)

  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  ! Performs first stage of IMEX-RK 3-4-3 method
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  !
  ! INPUT:
  !    dt     : time-step size
  !    nu     : kinematic viscosity
  !    kappa  : thermal diffusivity
  !    PVEL   : Chebyshev projector for VV equation
  !    Pmj    : Chebyshev projector
  !    VM     : Vertical velocity modes
  !    TM     : Temperature modes (and hori-vel modes)
  !    DVM    : 1st derivative of vertical velocity modes
  !    DTM    : 1st derivative of temperature modes
  !    D2VM   : 2nd derivative of vertical velocity modes
  !    D3VM   : 3rd derivative of vertical velocity modes
  !    GPVM   : Projected vertical velocity modes
  !    GPD2VM : Projected 2nd derivative vertical velocity modes
  !    GPD4VM : Projected 4th derivative vertical velocity modes
  !    GPTM   : Projected temperature modes for vertical velocity equation
  !    PVM    : Projected vertical velocity modes for temperature equation
  !    PDVM   : Projected first derivative of vertical velocity modes
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
real(dp)   , dimension(:,:),              intent(in)  :: PVEL, Pmj, VM, TM
real(dp)   , dimension(:,:),              intent(in)  :: DVM, DTM, D2VM, D3VM, GPTM
real(dp)   , dimension(:,:),              intent(in)  :: PVM, PDVM, PTM
complex(dp), allocatable, dimension(:,:), intent(in)  :: K1hV, K1hT
complex(dp), allocatable, dimension(:,:), intent(out) :: K1V, K1T, K2hV, K2hT
real(dp)   , allocatable, dimension(:,:)              :: lhs
real(dp)   , allocatable, dimension(:,:)              :: rhs
real(dp)   , allocatable, dimension(:,:)              :: tempmat, ups
real(dp)   , allocatable, dimension(:,:)              :: Kmat
real(dp)   , allocatable, dimension(:)                :: ore, oim
real(dp)   , allocatable, dimension(:)                :: Kre, Kim
complex(dp), allocatable, dimension(:)                :: temp
complex(dp), allocatable, dimension(:,:)              :: Aml1, Bml1
complex(dp), allocatable, dimension(:,:)              :: NLVml, NLTml
real(dp)                                              :: wave, wave2, wave4
integer                                               :: i
integer                                               :: NC, NP, NF2

! Some LAPACK and BLAS variables
integer,     allocatable, dimension(:)                :: ipiv
integer                                               :: info
real(dp), parameter                                   :: scale1=1.0_dp, scale2=0.0_dp
integer , parameter                                   :: incx=1, incy=1

NC  = size(Bml,1)
NF2 = size(Bml,2)
NP  = size(Tyx,1)

allocate(K1V(NC,NF2), K1T(NC,NF2)    , stat=alloc_err)
allocate(K2hV(NC,NF2),K2hT(NC,NF2)   , stat=alloc_err)
allocate(lhs(NC,NC)                  , stat=alloc_err)
allocate(rhs(NC,2)                   , stat=alloc_err)
allocate(tempmat(NC,NC), ups(NC,NC)  , stat=alloc_err)
allocate(Kmat(NC,NC)                 , stat=alloc_err)
allocate(Kre(NC), Kim(NC)            , stat=alloc_err)
allocate(ore(NP), oim(NP)            , stat=alloc_err)
allocate(ipiv(NC)                    , stat=alloc_err)
allocate(temp(NC)                    , stat=alloc_err)
allocate(Aml1(NC,NF2), Bml1(NC,NF2)  , stat=alloc_err)
allocate(NLVml(NC,NF2), NLTml(NC,NF2), stat=alloc_err)

K1T     = cmplx(0.0_dp, 0.0_dp, kind=dp)
K1V     = cmplx(0.0_dp, 0.0_dp, kind=dp)
K2hV    = cmplx(0.0_dp, 0.0_dp, kind=dp)
K2hT    = cmplx(0.0_dp, 0.0_dp, kind=dp)
lhs     = 0.0_dp
rhs     = 0.0_dp
tempmat = 0.0_dp
ups     = 0.0_dp
Kmat    = 0.0_dp
Kre     = 0.0_dp
Kim     = 0.0_dp
ore     = 0.0_dp
oim     = 0.0_dp
ipiv    = 0
temp    = cmplx(0.0_dp, 0.0_dp, kind=dp)
Aml1    = cmplx(0.0_dp, 0.0_dp, kind=dp)
Bml1    = cmplx(0.0_dp, 0.0_dp, kind=dp)
NLVml   = cmplx(0.0_dp, 0.0_dp, kind=dp)
NLTml   = cmplx(0.0_dp, 0.0_dp, kind=dp)

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
   Aml1(:,i) = cmplx(rhs(:,1), rhs(:,2), kind=dp)

   ! Get K1 for vertical velocity equation (avoids ever forming inv(upsilon))
   call dgemv('n', NC, NC, scale1, tempmat, NC, rhs(:,1), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, tempmat, NC, rhs(:,2), incx, scale2, Kim, incy)
   rhs(:,1) = Kre
   rhs(:,2) = Kim
   tempmat  = ups
   call dgesv(NC, 2, tempmat, NC, ipiv, rhs, NC, info)

   K1V(:,i) = nu*cmplx(rhs(:,1), rhs(:,2), kind=dp)

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
   K1T(:,i) = cmplx(Kre, Kim, kind=dp)

   ! Complete formation of Bml_1
   call dgemv('n', NC, NC, scale1, PTM, NC, rhs(:,1), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, PTM, NC, rhs(:,2), incx, scale2, Kim, incy)
   Bml1(:,i) = cmplx(Kre, Kim, kind=dp)
end do

call nonlinear_terms(NLVml,NLTml, PVEL,Pmj,VM,TM,DVM,DTM,D2VM,D3VM,Aml1,Bml1)

do i = 1,NF2
   wave  = kx(i)
   wave2 = wave**2.0_dp

   ! Tempreature equation
   call dgemv('n', NC, NC, scale1, PVM, NC, real (Aml1(:,i)), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, PVM, NC, aimag(Aml1(:,i)), incx, scale2, Kim, incy)

   rhs(:,1) = -real (NLTml(:,i)) + Kre
   rhs(:,2) = -aimag(NLTml(:,i)) + Kim

   ups = PTM ! Use ups as temp array
   call dgesv(NC, 2, ups, NC, ipiv, rhs, NC, info)
   K2hT(:,i) = cmplx(rhs(:,1), rhs(:,2), kind=dp)

   ! Vertical velocity equation
   call dgemv('n', NC, NC, scale1, GPTM, NC, real (Bml1(:,i)), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, GPTM, NC, aimag(Bml1(:,i)), incx, scale2, Kim, incy)

   rhs(:,1) = -(real (NLVml(:,i)) + wave2*Kre)
   rhs(:,2) = -(aimag(NLVml(:,i)) + wave2*Kim)

   ups = -wave2*GPVM + GPD2VM
   call dgesv(NC, 2, ups, NC, ipiv, rhs, NC, info)
   K2hV(:,i) = cmplx(rhs(:,1), rhs(:,2), kind=dp)
end do

end subroutine stage1

subroutine stage2(K2V,K2T,K3hV,K3hT,                                 &
               K1V,K1T,dt,nu,kappa,PVEL,Pmj,VM,TM,DVM,DTM,D2VM,D3VM, &
               GPVM,GPD2VM,GPD4VM,GPTM,PVM,PDVM,PTM,                 &
               K1hV,K2hV,K1hT,K2hT)

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
  !    PVEL   : Chebyshev projector for VV equation
  !    Pmj    : Chebyshev projector
  !    VM     : Vertical velocity modes
  !    TM     : Temperature modes (and hori-vel modes)
  !    DVM    : 1st derivative of vertical velocity modes
  !    DTM    : 1st derivative of temperature modes
  !    D2VM   : 2nd derivative of vertical velocity modes
  !    D3VM   : 3rd derivative of vertical velocity modes
  !    GPVM   : Projected vertical velocity modes
  !    GPD2VM : Projected 2nd derivative vertical velocity modes
  !    GPD4VM : Projected 4th derivative vertical velocity modes
  !    GPTM   : Projected temperature modes for vertical velocity equation
  !    PVM    : Projected vertical velocity modes for temperature equation
  !    PDVM   : Projected first derivative of vertical velocity modes
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
real(dp)   , dimension(:,:),              intent(in)  :: PVEL, Pmj, TM, VM
real(dp)   , dimension(:,:),              intent(in)  :: DVM,DTM,D2VM, D3VM
real(dp)   , dimension(:,:),              intent(in)  :: GPTM, PVM, PDVM, PTM
complex(dp), allocatable, dimension(:,:), intent(in)  :: K1V, K1T
complex(dp), allocatable, dimension(:,:), intent(in)  :: K1hV, K2hV
complex(dp), allocatable, dimension(:,:), intent(in)  :: K1hT, K2hT
complex(dp), allocatable, dimension(:,:), intent(out) :: K2V, K2T, K3hV, K3hT
real(dp)   , allocatable, dimension(:,:)              :: lhs
real(dp)   , allocatable, dimension(:,:)              :: rhs
real(dp)   , allocatable, dimension(:,:)              :: tempmat, ups
real(dp)   , allocatable, dimension(:,:)              :: Kmat
real(dp)   , allocatable, dimension(:)                :: ore, oim
real(dp)   , allocatable, dimension(:)                :: Kre, Kim
complex(dp), allocatable, dimension(:)                :: temp
complex(dp), allocatable, dimension(:,:)              :: Aml2, Bml2
complex(dp), allocatable, dimension(:,:)              :: NLVml, NLTml
real(dp)                                              :: wave, wave2, wave4
integer                                               :: i
integer                                               :: NC, NP, NF2

! Some LAPACK and BLAS variables
integer,     allocatable, dimension(:)                :: ipiv
integer                                               :: info
real(dp), parameter                                   :: scale1=1.0_dp, scale2=0.0_dp
integer , parameter                                   :: incx=1, incy=1

NC  = size(Bml,1)
NF2 = size(Bml,2)
NP  = size(Tyx,1)

allocate(K2V(NC,NF2), K2T(NC,NF2)    , stat=alloc_err)
allocate(K3hV(NC,NF2),K3hT(NC,NF2)   , stat=alloc_err)
allocate(lhs(NC,NC)                  , stat=alloc_err)
allocate(rhs(NC,2)                   , stat=alloc_err)
allocate(tempmat(NC,NC), ups(NC,NC)  , stat=alloc_err)
allocate(Kmat(NC,NC)                 , stat=alloc_err)
allocate(Kre(NC), Kim(NC)            , stat=alloc_err)
allocate(ore(NP), oim(NP)            , stat=alloc_err)
allocate(ipiv(NC)                    , stat=alloc_err)
allocate(temp(NC)                    , stat=alloc_err)
allocate(Aml2(NC,NF2), Bml2(NC,NF2)  , stat=alloc_err)
allocate(NLVml(NC,NF2), NLTml(NC,NF2), stat=alloc_err)

K2T     = cmplx(0.0_dp, 0.0_dp, kind=dp)
K2V     = cmplx(0.0_dp, 0.0_dp, kind=dp)
K3hV    = cmplx(0.0_dp, 0.0_dp, kind=dp)
K3hT    = cmplx(0.0_dp, 0.0_dp, kind=dp)
lhs     = 0.0_dp
rhs     = 0.0_dp
tempmat = 0.0_dp
ups     = 0.0_dp
Kmat    = 0.0_dp
Kre     = 0.0_dp
Kim     = 0.0_dp
ore     = 0.0_dp
oim     = 0.0_dp
ipiv    = 0
temp    = cmplx(0.0_dp, 0.0_dp, kind=dp)
Aml2    = cmplx(0.0_dp, 0.0_dp, kind=dp)
Bml2    = cmplx(0.0_dp, 0.0_dp, kind=dp)
NLVml   = cmplx(0.0_dp, 0.0_dp, kind=dp)
NLTml   = cmplx(0.0_dp, 0.0_dp, kind=dp)

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
   Aml2(:,i) = cmplx(rhs(:,1), rhs(:,2), kind=dp)

   ! Get K2 for vertical velocity equation (avoids ever forming inv(upsilon))
   call dgemv('n', NC, NC, scale1, tempmat, NC, rhs(:,1), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, tempmat, NC, rhs(:,2), incx, scale2, Kim, incy)
   rhs(:,1) = Kre
   rhs(:,2) = Kim
   tempmat  = ups
   call dgesv(NC, 2, tempmat, NC, ipiv, rhs, NC, info)

   K2V(:,i) = nu*cmplx(rhs(:,1), rhs(:,2), kind=dp)

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
   K2T(:,i) = cmplx(Kre, Kim, kind=dp)

   ! Complete formation of Bml_2
   call dgemv('n', NC, NC, scale1, PTM, NC, rhs(:,1), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, PTM, NC, rhs(:,2), incx, scale2, Kim, incy)
   Bml2(:,i) = cmplx(Kre, Kim, kind=dp)
end do

call nonlinear_terms(NLVml,NLTml, PVEL,Pmj,VM,TM,DVM,DTM,D2VM,D3VM,Aml2,Bml2)

do i = 1,NF2
   wave  = kx(i)
   wave2 = wave**2.0_dp

   ! Tempreature equation
   call dgemv('n', NC, NC, scale1, PVM, NC, real (Aml2(:,i)), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, PVM, NC, aimag(Aml2(:,i)), incx, scale2, Kim, incy)

   rhs(:,1) = -real (NLTml(:,i)) + Kre
   rhs(:,2) = -aimag(NLTml(:,i)) + Kim

   ups = PTM ! Use ups as temp array
   call dgesv(NC, 2, ups, NC, ipiv, rhs, NC, info)
   K3hT(:,i) = cmplx(rhs(:,1), rhs(:,2), kind=dp)

   ! Vertical velocity equation
   call dgemv('n', NC, NC, scale1, GPTM, NC, real (Bml2(:,i)), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, GPTM, NC, aimag(Bml2(:,i)), incx, scale2, Kim, incy)

   rhs(:,1) = -(real (NLVml(:,i)) + wave2*Kre)
   rhs(:,2) = -(aimag(NLVml(:,i)) + wave2*Kim)

   ups = -wave2*GPVM + GPD2VM
   call dgesv(NC, 2, ups, NC, ipiv, rhs, NC, info)
   K3hV(:,i) = cmplx(rhs(:,1), rhs(:,2), kind=dp)
end do

end subroutine stage2

subroutine stage3(K3V,K3T,K4hV,K4hT,                                          &
               K1V,K1T,K2V,K2T,dt,nu,kappa,PVEL,Pmj,VM,TM,DVM,DTM,D2VM,D3VM, &
               GPVM,GPD2VM,GPD4VM,GPTM,PVM,PDVM,PTM,                          &
               K1hV,K2hV,K3hV,K1hT,K2hT,K3hT)

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
  !    PVEL   : Chebyshev projector for VV equation
  !    Pmj    : Chebyshev projector
  !    VM     : Vertical velocity modes
  !    TM     : Temperature modes (and hori-vel modes)
  !    DVM    : 1st derivative of vertical velocity modes
  !    DTM    : 1st derivative of temperature modes
  !    D2VM   : 2nd derivative of vertical velocity modes
  !    D3VM   : 3rd derivative of vertical velocity modes
  !    GPVM   : Projected vertical velocity modes
  !    GPD2VM : Projected 2nd derivative vertical velocity modes
  !    GPD4VM : Projected 4th derivative vertical velocity modes
  !    GPTM   : Projected temperature modes for vertical velocity equation
  !    PVM    : Projected vertical velocity modes for temperature equation
  !    PDVM   : Projected first derivative of vertical velocity modes
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
real(dp)   , dimension(:,:),              intent(in)  :: PVEL, Pmj, VM, TM
real(dp)   , dimension(:,:),              intent(in)  :: DVM, DTM, D2VM, D3VM, GPTM
real(dp)   , dimension(:,:),              intent(in)  :: PVM, PDVM, PTM
complex(dp), allocatable, dimension(:,:), intent(in)  :: K1V, K1T, K2V, K2T
complex(dp), allocatable, dimension(:,:), intent(in)  :: K1hV, K2hV, K3hV
complex(dp), allocatable, dimension(:,:), intent(in)  :: K1hT, K2hT, K3hT
complex(dp), allocatable, dimension(:,:), intent(out) :: K3V, K3T, K4hV, K4hT
real(dp)   , allocatable, dimension(:,:)              :: lhs
real(dp)   , allocatable, dimension(:,:)              :: rhs
real(dp)   , allocatable, dimension(:,:)              :: tempmat, ups
real(dp)   , allocatable, dimension(:,:)              :: Kmat
real(dp)   , allocatable, dimension(:)                :: Kre, Kim
real(dp), allocatable, dimension(:)                   :: ore, oim
complex(dp), allocatable, dimension(:)                :: temp
complex(dp), allocatable, dimension(:,:)              :: Aml3, Bml3
complex(dp), allocatable, dimension(:,:)              :: NLVml, NLTml
real(dp)                                              :: wave, wave2, wave4
integer                                               :: i
integer                                               :: NC, NP, NF2

! Some LAPACK and BLAS variables
integer,     allocatable, dimension(:)                :: ipiv
integer                                               :: info
real(dp), parameter                                   :: scale1=1.0_dp, scale2=0.0_dp
integer , parameter                                   :: incx=1, incy=1

NC  = size(Bml,1)
NF2 = size(Bml,2)
NP  = size(Tyx,1)

allocate(K3V(NC,NF2), K3T(NC,NF2)    , stat=alloc_err)
allocate(K4hV(NC,NF2),K4hT(NC,NF2)   , stat=alloc_err)
allocate(lhs(NC,NC)                  , stat=alloc_err)
allocate(rhs(NC,2)                   , stat=alloc_err)
allocate(tempmat(NC,NC), ups(NC,NC)  , stat=alloc_err)
allocate(Kmat(NC,NC)                 , stat=alloc_err)
allocate(Kre(NC), Kim(NC)            , stat=alloc_err)
allocate(ore(NP), oim(NP)            , stat=alloc_err)
allocate(ipiv(NC)                    , stat=alloc_err)
allocate(temp(NC)                    , stat=alloc_err)
allocate(Aml3(NC,NF2), Bml3(NC,NF2)  , stat=alloc_err)
allocate(NLVml(NC,NF2), NLTml(NC,NF2), stat=alloc_err)

K3T     = cmplx(0.0_dp, 0.0_dp, kind=dp)
K3V     = cmplx(0.0_dp, 0.0_dp, kind=dp)
K4hV    = cmplx(0.0_dp, 0.0_dp, kind=dp)
K4hT    = cmplx(0.0_dp, 0.0_dp, kind=dp)
lhs     = 0.0_dp
rhs     = 0.0_dp
tempmat = 0.0_dp
ups     = 0.0_dp
Kmat    = 0.0_dp
Kre     = 0.0_dp
Kim     = 0.0_dp
ore     = 0.0_dp
oim     = 0.0_dp
ipiv    = 0
temp    = cmplx(0.0_dp, 0.0_dp, kind=dp)
Aml3    = cmplx(0.0_dp, 0.0_dp, kind=dp)
Bml3    = cmplx(0.0_dp, 0.0_dp, kind=dp)
NLVml   = cmplx(0.0_dp, 0.0_dp, kind=dp)
NLTml   = cmplx(0.0_dp, 0.0_dp, kind=dp)

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
   Aml3(:,i) = cmplx(rhs(:,1), rhs(:,2), kind=dp)

   ! Get K3 for vertical velocity equation (avoids ever forming inv(upsilon))
   call dgemv('n', NC, NC, scale1, tempmat, NC, rhs(:,1), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, tempmat, NC, rhs(:,2), incx, scale2, Kim, incy)
   rhs(:,1) = Kre
   rhs(:,2) = Kim
   tempmat  = ups
   call dgesv(NC, 2, tempmat, NC, ipiv, rhs, NC, info)

   K3V(:,i) = nu*cmplx(rhs(:,1), rhs(:,2), kind=dp)

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
   K3T(:,i) = cmplx(Kre, Kim, kind=dp)

   ! Complete formation of Bml_3
   call dgemv('n', NC, NC, scale1, PTM, NC, rhs(:,1), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, PTM, NC, rhs(:,2), incx, scale2, Kim, incy)
   Bml3(:,i) = cmplx(Kre, Kim, kind=dp)
end do

call nonlinear_terms(NLVml,NLTml, PVEL,Pmj,VM,TM,DVM,DTM,D2VM,D3VM,Aml3,Bml3)

do i = 1,NF2
   wave  = kx(i)
   wave2 = wave**2.0_dp

   ! Tempreature equation
   call dgemv('n', NC, NC, scale1, PVM, NC, real (Aml3(:,i)), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, PVM, NC, aimag(Aml3(:,i)), incx, scale2, Kim, incy)

   rhs(:,1) = -real (NLTml(:,i)) + Kre
   rhs(:,2) = -aimag(NLTml(:,i)) + Kim

   ups = PTM ! Use ups as temp array
   call dgesv(NC, 2, ups, NC, ipiv, rhs, NC, info)
   K4hT(:,i) = cmplx(rhs(:,1), rhs(:,2), kind=dp)

   ! Vertical velocity equation
   call dgemv('n', NC, NC, scale1, GPTM, NC, real (Bml3(:,i)), incx, scale2, Kre, incy)
   call dgemv('n', NC, NC, scale1, GPTM, NC, aimag(Bml3(:,i)), incx, scale2, Kim, incy)

   rhs(:,1) = -(real (NLVml(:,i)) + wave2*Kre)
   rhs(:,2) = -(aimag(NLVml(:,i)) + wave2*Kim)

   ups = -wave2*GPVM + GPD2VM
   call dgesv(NC, 2, ups, NC, ipiv, rhs, NC, info)
   K4hV(:,i) = cmplx(rhs(:,1), rhs(:,2), kind=dp)
end do

end subroutine stage3

subroutine update_u(Pmj,DVM)

  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  ! Updates horizontal velocity 
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  !
  ! INPUT:
  !    PVEL: VV Projector
  !    DVM : First derivative of vertical velocity modes
  !
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

implicit none

real(dp),              dimension(:,:), intent(in)  :: Pmj, DVM
real(dp), allocatable, dimension(:)                :: ore, oim
real(dp), allocatable, dimension(:)                :: Kre, Kim
real(dp)                                           :: wave
integer                                            :: i
integer                                            :: NC, NP, NF2

! Some LAPACK and BLAS variables
integer,  allocatable, dimension(:)                :: ipiv
real(dp), parameter                                :: scale1=1.0_dp, scale2=0.0_dp
integer , parameter                                :: incx=1, incy=1

NC  = size(Bml,1)
NF2 = size(Bml,2)
NP  = size(Tyx,1)

allocate(Kre(NC), Kim(NC), stat=alloc_err)
allocate(ore(NP), oim(NP), stat=alloc_err)
allocate(ipiv(NC)        , stat=alloc_err)

Kre     = 0.0_dp
Kim     = 0.0_dp
ore     = 0.0_dp
oim     = 0.0_dp
ipiv    = 0

do i = 1,NF2
   ! Continuity equation
   wave  = kx(i)

   ! Get Uml
   if (i /= 1) then
      call dgemv('n', NP, NC, scale1, DVM,  NP, real (Aml(:,i)), incx, scale2, ore, incy)
      call dgemv('n', NP, NC, scale1, DVM,  NP, aimag(Aml(:,i)), incx, scale2, oim, incy)
      call dgemv('n', NC, NP, scale1, Pmj,  NC, ore, incx, scale2, Kre, incy)
      call dgemv('n', NC, NP, scale1, Pmj,  NC, oim, incx, scale2, Kim, incy)
      Uml(:,i) = CI*cmplx(Kre, Kim, kind=dp)/wave
   else
      Uml(:,i) = cmplx(0.0_dp, 0.0_dp, kind=dp) ! No mean flow for now
   end if

end do

end subroutine update_u

subroutine probes(Vprobe, Tprobe, VM, TM, iter, interval, iprobe)

implicit none

integer,                  intent(in)     :: iter, interval
logical,                  intent(in)     :: iprobe
real(dp), dimension(:,:), intent(in)     :: VM, TM
real(dp),                 intent(out)    :: Vprobe, Tprobe
complex(dp), allocatable, dimension(:,:) :: temp1, temp2
integer                                  :: NP, NF, NC, NF2
integer                                  :: fprint
character(10)                            :: citer
character(9)                             :: fiter

! Some LAPACK and BLAS parameters
real(dp), parameter :: scale1=1.0_dp, scale2=0.0_dp

fprint = mod(iter, interval)
if ( (fprint == 0).or.(iprobe)) then
    ! Use temperature arrays to get these sizes.  Could use V arrays too.
    NC  = size(Bml, 1)
    NF  = size(Bmx, 2)
    NP  = size(Tyx, 1)
    NF2 = NF/2 + 1
    
    allocate(temp1(NC,NF2), temp2(NC,NF2), stat=alloc_err)
    temp1 = cmplx(0.0_dp, 0.0_dp, kind=dp)
    temp2 = cmplx(0.0_dp, 0.0_dp, kind=dp)
    
    ! Bring fields to physical space
    ! temporary storage because ifft overwrites input arrays!
    temp1 = Aml
    temp2 = Bml
    call fftw_execute_dft_c2r(ipV, Aml, Amx) ! Vertical velocity Fourier to physical
    call fftw_execute_dft_c2r(ipT, Bml, Bmx) ! Temperature Fourier to physical
    call dgemm('n', 'n', NP, NF, NC, scale1, VM, NP, Amx, NC, scale2, Vyx, NP) ! Cheb to physical (VV)
    call dgemm('n', 'n', NP, NF, NC, scale1, TM, NP, Bmx, NC, scale2, Tyx, NP) ! Cheb to physical (T)
    Aml = temp1
    Bml = temp2
    
    if (fprint == 0) then
       ! Write out the fields
       write(citer, "(I10)") 1000000000 + iter
       fiter = citer(2:10)
       call write_out_mat(Vyx, "VV/Vyx"//fiter)
       call write_out_mat(Tyx, "Th/Tyx"//fiter)
    end if

    if (iprobe) then
       ! Get trace of Vyx and Tyx
       Vprobe = Vyx(ceiling(NP/2.0_dp), NF/2)
       Tprobe = Tyx(ceiling(NP/2.0_dp), NF/2)
    end if

end if

end subroutine probes

subroutine update_dt(VM,TM,dt)

implicit none
real(dp),parameter                       :: cfl = 1.0_dp
real(dp),parameter                       :: dt_ramp = 5.0_dp
real(dp), dimension(:,:), intent(in)     :: VM, TM
real(dp)                                 :: tmpmax,tmp
real(dp)                                 :: dt,dt_old
real(dp)                                 :: dxmin,dymin,alpha
complex(dp), allocatable, dimension(:,:) :: temp1, temp2
integer                                  :: NP, NF, NC, NF2
integer                                  :: ii,jj
! Some LAPACK and BLAS parameters
real(dp), parameter :: scale1=1.0_dp, scale2=0.0_dp

! Use temperature arrays to get these sizes.  Could use V arrays too.
NC  = size(Bml, 1)
NF  = size(Bmx, 2)
NP  = size(Tyx, 1)
NF2 = NF/2 + 1

allocate(temp1(NC,NF2), temp2(NC,NF2), stat=alloc_err)
temp1 = cmplx(0.0_dp, 0.0_dp, kind=dp)
temp2 = cmplx(0.0_dp, 0.0_dp, kind=dp)

! Bring fields to physical space
! temporary storage because ifft overwrites input arrays!
temp1 = Aml
temp2 = Uml


call fftw_execute_dft_c2r(ipV, Aml, Amx) ! Vertical velocity Fourier to physical
call fftw_execute_dft_c2r(ipU, Uml, Umx) ! Horizontal velocity Fourier to physical

call dgemm('n', 'n', NP, NF, NC, scale1, VM, NP, Amx, NC, scale2, Vyx, NP) ! Cheb to physical (VV)
call dgemm('n', 'n', NP, NF, NC, scale1, TM, NP, Umx, NC, scale2, Uyx, NP) ! Cheb to physical (U)

Aml = temp1
Uml = temp2

alpha = kx(2)

dxmin = 2.0*pi/(alpha*NF)
dymin = abs(y(1,1)-y(2,1))

tmp    = 0.0_dp
tmpmax = 0.0_dp
do jj = 1,NP
   do ii = 1,NF
      tmp = 2.0_dp*pi*abs(Uyx(jj,ii)) / dxmin + abs(Vyx(jj,ii)) / dymin
      if (tmp > tmpmax) then
         tmpmax = tmp
      end if
   end do
end do

dt_old = dt

dt = cfl/tmpmax

if(dt > dt_old*dt_ramp) then
  dt = dt_old*dt_ramp
endif

end subroutine update_dt

subroutine backward_euler(NC,NF,dt,t_final,nu,kappa,        &
                          PVEL,Pmj,VM,TM,DVM,DTM,D2VM,D3VM, &
                          GPVM,GPTM,PVM,PTM,GPD2VM,GPD4VM)

integer                 , intent(in)     :: NC, NF
real(dp),                 intent(in)     :: dt, t_final, kappa, nu
real(dp), dimension(:,:), intent(in)     :: VM, TM, PVM, PTM, GPTM
real(dp), dimension(:,:), intent(in)     :: PVEL, Pmj
real(dp), dimension(:,:), intent(in)     :: DVM, DTM, D2VM, D3VM
real(dp), dimension(:,:), intent(in)     :: GPVM, GPD2VM, GPD4VM
real(dp)                                 :: dt_final, time
real(dp)                                 :: wave, wave2, wave4
real(dp)                                 :: Vprobe, Tprobe
real(dp), allocatable, dimension(:,:)    :: lhs, rhs
real(dp), allocatable, dimension(:,:)    :: upsilon
complex(dp), allocatable, dimension(:)   :: temprhs
complex(dp), allocatable, dimension(:,:) :: rhsT_hold, rhsV_hold
complex(dp), allocatable, dimension(:,:) :: NLVml, NLTml
real(dp), parameter                      :: scale1=1.0_dp, scale2=0.0_dp
integer,  parameter                      :: incx=1, incy=1
integer,  allocatable, dimension(:)      :: ipiv
integer                                  :: info
integer                                  :: i, tstep
integer                                  :: NF2

NF2 = NF/2 + 1

allocate(upsilon(NC,NC)   , stat=alloc_err)
allocate(lhs(NC,NC)       , stat=alloc_err)
allocate(rhs(NC,2)        , stat=alloc_err)
allocate(ipiv(NC)         , stat=alloc_err)
allocate(temprhs(NC)      , stat=alloc_err)
allocate(rhsT_hold(NC,NF2), stat=alloc_err)
allocate(rhsV_hold(NC,NF2), stat=alloc_err)
allocate(NLTml(NC,NF2)    , stat=alloc_err)
allocate(NLVml(NC,NF2)    , stat=alloc_err)
upsilon   = 0.0_dp
lhs       = 0.0_dp
rhs       = 0.0_dp
ipiv      = 0
temprhs   = cmplx(0.0_dp, 0.0_dp, kind=dp)
rhsT_hold = cmplx(0.0_dp, 0.0_dp, kind=dp)
rhsT_hold = cmplx(0.0_dp, 0.0_dp, kind=dp)
NLTml     = cmplx(0.0_dp, 0.0_dp, kind=dp)
NLVml     = cmplx(0.0_dp, 0.0_dp, kind=dp)

! Initial time
time  = 0.0_dp
tstep = 0

!open(unit=1501, file="maxval.txt", action="write", status="unknown", position="append")

do ! while time < t_final

   ! Do some computations using data from previous time step

   ! Bring to physical space to track decay rate
!   call decay(maxVyx, maxTyx, VM, TM, tstep)
!   write(1501, fmt=3000) time, maxVyx, maxTyx

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

   do i = 1,NF/2+1
      wave  = kx(i)
      wave2 = wave**2.0_dp
      ! Vertical velocity equation
      upsilon = -wave2*GPVM + GPD2VM
      ! Compute the RHS
      call dgemv('n', NC, NC, scale1, upsilon, NC, real (Aml(:,i)), incx, scale2, rhs(:,1), incy)
      call dgemv('n', NC, NC, scale1, upsilon, NC, aimag(Aml(:,i)), incx, scale2, rhs(:,2), incy)
      temprhs = cmplx(rhs(:,1), rhs(:,2), kind=dp)
      call dgemv('n', NC, NC, scale1, GPTM, NC, real (Bml(:,i)), incx, scale2, rhs(:,1), incy)
      call dgemv('n', NC, NC, scale1, GPTM, NC, aimag(Bml(:,i)), incx, scale2, rhs(:,2), incy)
      rhsV_hold(:,i) = temprhs - wave2*dt*cmplx(rhs(:,1), rhs(:,2), kind=dp)

      ! Temperature Equation
      ! Compute the RHS
      call dgemv('n', NC, NC, scale1, PTM, NC, real (Bml(:,i)), incx, scale2, rhs(:,1), incy)
      call dgemv('n', NC, NC, scale1, PTM, NC, aimag(Bml(:,i)), incx, scale2, rhs(:,2), incy)
      temprhs = cmplx(rhs(:,1), rhs(:,2), kind=dp)
      call dgemv('n', NC, NC, scale1, PVM, NC, real (Aml(:,i)), incx, scale2, rhs(:,1), incy)
      call dgemv('n', NC, NC, scale1, PVM, NC, aimag(Aml(:,i)), incx, scale2, rhs(:,2), incy)
      rhsT_hold(:,i) = temprhs + dt*cmplx(rhs(:,1), rhs(:,2), kind=dp)
   end do

   call nonlinear_terms(NLVml,NLTml, PVEL,Pmj,VM,TM,DVM,DTM,D2VM,D3VM,Aml,Bml)

   do i = 1,NF2
      wave  = kx(i)
      wave2 = wave**2.0_dp
      wave4 = wave2**2.0_dp

      ! VV equation
      upsilon = -wave2*GPVM + GPD2VM
      lhs = upsilon - nu*dt*(wave4*GPVM - 2.0_dp*wave2*GPD2VM + GPD4VM)
      rhs(:,1) = real (rhsV_hold(:,i)) + dt*real (NLVml(:,i))
      rhs(:,2) = aimag(rhsV_hold(:,i)) + dt*aimag(NLVml(:,i))
      ! Solve for the solution at next time step
      call dgesv(NC, 2, lhs, NC, ipiv, rhs, NC, info)
      ! Update solution
      Aml(:,i) = cmplx(rhs(:,1), rhs(:,2), kind=dp)

      ! Temperature equation
      lhs = (1.0_dp + kappa*dt*wave2)*PTM - kappa*dt*eye
      rhs(:,1) = real (rhsT_hold(:,i)) + dt*real (NLTml(:,i))
      rhs(:,2) = aimag(rhsT_hold(:,i)) + dt*aimag(NLTml(:,i))
      ! Solve for solution at next time
      call dgesv(NC, 2, lhs, NC, ipiv, rhs, NC, info)
      ! Update solution
      Bml(:,i) = cmplx(rhs(:,1), rhs(:,2), kind=dp)
   end do

end do

!close(unit=1501)

!2000 format(E25.16E3, E25.16E3          )
!3000 format(E25.16E3, E25.16E3, E25.16E3)

end subroutine backward_euler

end module time_integrators
