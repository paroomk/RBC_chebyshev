program heat

use types,  only: dp
use global, only: set_imex_params, fft_utils,    &
                  alloc_err, Aml, Bml, Amx, Bmx, &
                  Vyx, Tyx, eye, kx,             &
                  pV, ipV, pT, ipT
use makeICs
use chebutils
use fftw
use write_pack
use time_integrators, only: imex_rk, backward_euler

! Set up basic problem parameters and fields
integer                                  :: NC, NP, NF, LT  ! Problem size
real(dp), allocatable, dimension(:,:)    :: Cjm, Pmj        ! Cheb matrices
real(dp), allocatable, dimension(:,:)    :: TM, DTM         ! Temperature modes
real(dp), allocatable, dimension(:,:)    :: VM, DVM         ! Velocity
real(dp), allocatable, dimension(:,:)    :: D2VM, D3VM      ! modes
real(dp), allocatable, dimension(:,:)    :: PVM, PTM, PD2TM ! Projected Galerkin for T-eq
real(dp), allocatable, dimension(:,:)    :: GPTM, GPVM      ! Projected Galerkin for v-eq
real(dp), allocatable, dimension(:,:)    :: GPD2VM          ! Projected Galerkin for v-eq
real(dp), allocatable, dimension(:,:)    :: GPD4VM          ! Projected Galerkin for v-eq
real(dp), allocatable, dimension(:,:)    :: PVEL            ! Projector for v-eq
real(dp), allocatable, dimension(:,:)    :: PTMinv          ! Inverse of PTM
real(dp), allocatable, dimension(:,:)    :: PVMinv          ! Inverse of GPVM
real(dp), allocatable, dimension(:,:)    :: PhysChebT       ! Transformation from
                                                            ! physical space to
                                                            ! Chebyshev space of Tyx
real(dp), allocatable, dimension(:,:)    :: PhysChebV       ! Transformation from
                                                            ! physical space to
                                                            ! Chebyshev space of Vyx
real(dp), allocatable, dimension(:,:)    :: y               ! y-coordinate
real(dp)                                 :: amp             ! Initial Temperature
                                                            ! amplitude
!real(dp), parameter                      :: alpha   = 1.5585_dp
real(dp)                                 :: alpha
real(dp), parameter                      :: nu      = 0.03_dp
real(dp), parameter                      :: kappa   = 0.1_dp
real(dp), parameter                      :: t_final = 5.0_dp
real(dp), parameter                      :: dt      = 0.1_dp

! Some Lapack and Blas parameters
real(dp), parameter                      :: scale1=1.0_dp, scale2=0.0_dp ! Scalars
integer                                  :: info
integer , allocatable, dimension(:)      :: ipiv
real(dp), allocatable, dimension(:)      :: work

alpha = pi/(2.0_dp*sqrt(2.0))

NC = 10
NP = NC + 4
NF = 8

!NC = 4
!NP = 6
!NF = 4

LT = NF-1

if (NP < NC+4) then
   write(*,*) "NP < NC+4 so setting NP = NC+4"
   NP = NC+4
   write(*,*)
end if

! Allocate temperature and velocity fields
allocate(Tyx(NP,NF), Vyx(NP,NF), stat=alloc_err)
Tyx = 0.0_dp
Vyx = 0.0_dp

! Set up FFT plans
call fft_utils(NC,NF)

Bmx = 0.0_dp
Bml = cmplx(0.0_dp, 0.0_dp)

Amx = 0.0_dp
Aml = cmplx(0.0_dp, 0.0_dp)

! Allocate expansion-projection matrices
allocate(Cjm(NP,NC) , stat=alloc_err)
allocate(Pmj(NC,NP) , stat=alloc_err)

Cjm = 0.0_dp
Pmj = 0.0_dp

! Allocate Galerkin trial functions
allocate(TM(NP,NC)    , DTM(NP,NC)    , stat=alloc_err)
allocate(VM(NP,NC)    , DVM(NP,NC)    , stat=alloc_err)
allocate(D2VM(NP,NC)  , D3VM(NP,NC)   , stat=alloc_err)
allocate(PTM(NC,NC)   , PD2TM(NC,NC)  , stat=alloc_err)
allocate(GPTM(NC,NC)  , GPVM(NC,NC)   , stat=alloc_err)
allocate(GPD2VM(NC,NC), GPD4VM(NC,NC) , stat=alloc_err)
allocate(PVEL(NC,NP)                  , stat=alloc_err)
allocate(PTMinv(NC,NC)                , stat=alloc_err)
allocate(PhysChebV(NC,NP)             , stat=alloc_err)
allocate(PhysChebT(NC,NP)             , stat=alloc_err)

TM        = 0.0_dp
DTM       = 0.0_dp
VM        = 0.0_dp
DVM       = 0.0_dp
D2VM      = 0.0_dp
D3VM      = 0.0_dp
PTM       = 0.0_dp
PD2TM     = 0.0_dp
GPVM      = 0.0_dp
GPD2VM    = 0.0_dp
GPD4VM    = 0.0_dp
PVEL      = 0.0_dp
PTMinv    = 0.0_dp
PhysChebV = 0.0_dp
PhysChebT = 0.0_dp

! Allocate misc arrays
allocate(y(NP,1)   , stat=alloc_err)
allocate(eye(NC,NC), stat=alloc_err)
allocate(kx(NF)    , stat=alloc_err)

y   = 0.0_dp
eye = 0.0_dp
kx  = 0.0_dp

forall(j=1:NC) eye(j,j) = 1.0_dp ! identity matrix

! Wavenumbers
do i = 1,NF/2
   kx(i) = real(i,kind=dp) - 1.0_dp
end do
do i = NF/2+1, NF
   kx(i) = real(i-NF,kind=dp) - 1.0_dp
end do
kx = alpha*kx

! Allocate some LAPACK and BLAS arrays
allocate(ipiv(NC), stat=alloc_err)
allocate(work(NC), stat=alloc_err)
ipiv = 0
work = 0

! Get the Chebyshev Galerkin functions
call makeVTM(VM, TM, DVM, DTM, D2VM, D3VM, Cjm, Pmj, NC,NP)

y(:,1) = Cjm(:,2)

call write_out_mat(y, "y")

! Get the projected Galerkin functions
call projectVT(GPVM,GPD2VM,GPD4VM,GPTM,PVM,PTM,PD2TM,PVEL, Pmj,Cjm,VM,D2VM,TM)

! Compute inv(PTM) for use in Phys-Cheb transformation of Tyx 
PTMinv = PTM
call dgetrf(NC, NC, PTMinv, NC, ipiv, info)       ! LU factorization
call dgetri(NC, PTMinv, NC, ipiv, work, NC, info) ! Inverse

! Compute inv(PTM)*Pmj for Phys-Cheb transformation of Tyx
call dgemm('n', 'n', NC, NP, NC, scale1, PTMinv, NC, Pmj, NC, scale2, PhysChebT, NC)

! Compute inv(PVEL*VM) for use in Phys-Cheb transformation of Vyx 
PVMinv = GPVM
call dgetrf(NC, NC, PVMinv, NC, ipiv, info)       ! LU factorization
call dgetri(NC, PVMinv, NC, ipiv, work, NC, info) ! Inverse

! Compute inv(PVM)*PVEL for Phys-Cheb transformation of Vyx
call dgemm('n', 'n', NC, NP, NC, scale1, PVMinv, NC, PVEL, NC, scale2, PhysChebV, NC)

! Create some initial conditions
amp = 0.03_dp
call initial_conditions(Pmj,y,amp,NC)

! Call time-integrator
call imex_rk(NC, NF, dt, t_final, nu, kappa, &
             PhysChebV, PhysChebT,VM,TM,     &
             GPVM,GPTM,PVM,PTM,GPD2VM,GPD4VM)
!call backward_euler(NC,NF,dt,t_final,nu,kappa, &
!                    PhysChebV, PhysChebT,VM,TM,&
!                    GPVM,GPTM,PVM,PTM,GPD2VM,GPD4VM)

1000 format(E25.16E3                           )
!2000 format(E25.16E3,E25.16E3                  )
!3000 format(E25.16E3,E25.16E3,E25.16E3         )
!4000 format(E25.16E3,E25.16E3,E25.16E3,E25.16E3)

end program heat

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
