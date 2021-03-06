program TC 

use types,  only: dp
use global, only: set_imex_params, fft_utils,      &
                  alloc_err, eye, kx,              &
                  Aml, Bml, Uml, Amx, Bmx, Umx,    &
                  Vyx, Tyx, Uyx, fft1_ml, fft1_mx, &
                  fft1_yl, fft1_yx, pf1,y

use makeICs
use chebutils
use fftw
use write_pack
use time_integrators, only: imex_rk, backward_euler, update_u, probes
use statutils, only: Nusselt 

! Set up basic problem parameters and fields
integer                                  :: NC, NP, NF      ! Problem size
real(dp), allocatable, dimension(:,:)    :: Cjm, Pmj        ! Cheb matrices
real(dp), allocatable, dimension(:,:)    :: TM,DTM,D2TM,DTMb! Temperature modes
real(dp), allocatable, dimension(:,:)    :: VM, DVM         ! Velocity
real(dp), allocatable, dimension(:,:)    :: D2VM, D3VM      ! modes
real(dp), allocatable, dimension(:,:)    :: PVM, PTM, PD2TM ! Projected Galerkin for T-eq
real(dp), allocatable, dimension(:,:)    :: PDVM, PDTM      ! Projected Galerkin for first-derivs 
real(dp), allocatable, dimension(:,:)    :: GPTM, GPVM      ! Projected Galerkin for v-eq
real(dp), allocatable, dimension(:,:)    :: GPD2VM          ! Projected Galerkin for v-eq
real(dp), allocatable, dimension(:,:)    :: GPD4VM          ! Projected Galerkin for v-eq
real(dp), allocatable, dimension(:,:)    :: PVEL            ! Projector for v-eq
!real(dp), allocatable, dimension(:,:)    :: y               ! y-coordinate
real(dp)                                 :: amp             ! Initial Temperature amplitude
real(dp)                                 :: Nuss            ! Nusselt number
real(dp), parameter                      :: alpha   = 1.5581
real(dp)                                 :: nu, kappa
real(dp), parameter                      :: Ra = 5.277302244797398E+06_dp, Pr = 7.0_dp
real(dp), parameter                      :: t_final = 500.0_dp
real(dp)                                 :: dt      = 0.005_dp

logical                                  :: read_ICs = .false.
logical                                  :: mirror_symmetry = .true. ! To apply mirror symmetry [u,v,T] (x,y) = [-u,v,T] (-x,y)
real(dp)                                 :: x, dx

!alpha = pi/(2.0_dp*sqrt(2.0))

nu    = dsqrt(16.0_dp*Pr/Ra)
!nu    = dsqrt(Pr/Ra)
kappa = dsqrt(16.0_dp/(Pr*Ra))
!kappa = dsqrt((Pr*Ra))

! Ra = 4.0e+06
!NC = 100
!NP = NC + 4
!NF = 128

! Ra = 4.0e+06
NC = 100
NP = NC + 4
NF = 128

! Original params
!NC = 45
!NP = NC + 4
!NF = 32

if (NP < NC+4) then
   write(*,*) "NP < NC+4 so setting NP = NC+4"
   NP = NC+4
   write(*,*)
end if

! Allocate temperature and velocity fields
allocate(Tyx(NP,NF), Vyx(NP,NF), Uyx(NP,NF), stat=alloc_err)
Tyx   = 0.0_dp
Vyx   = 0.0_dp
Uyx   = 0.0_dp

! Set up FFT plans
call fft_utils(NP,NC,NF)

Bmx = 0.0_dp
Tyx = 0.0_dp
Bml = cmplx(0.0_dp, 0.0_dp, kind=dp)

Amx = 0.0_dp
Vyx = 0.0_dp
Aml = cmplx(0.0_dp, 0.0_dp, kind=dp)

Umx = 0.0_dp
Uyx = 0.0_dp
Uml = cmplx(0.0_dp, 0.0_dp, kind=dp)

fft1_mx = 0.0_dp
fft1_ml = cmplx(0.0_dp, 0.0_dp, kind=dp)

fft1_yx = 0.0_dp
fft1_yl = cmplx(0.0_dp, 0.0_dp, kind=dp)

! Allocate expansion-projection matrices
allocate(Cjm(NP,NC) , stat=alloc_err)
allocate(Pmj(NC,NP) , stat=alloc_err)

Cjm = 0.0_dp
Pmj = 0.0_dp

! Allocate Galerkin trial functions
allocate(TM(NP,NC)    , DTM(NP,NC)    , stat=alloc_err)
allocate(D2TM(NP,NC)                  , stat=alloc_err)
allocate(VM(NP,NC)    , DVM(NP,NC)    , stat=alloc_err)
allocate(D2VM(NP,NC)  , D3VM(NP,NC)   , stat=alloc_err)
allocate(PTM(NC,NC)   , PD2TM(NC,NC)  , stat=alloc_err)
allocate(PDVM(NC,NC)  , PDTM(NC,NC)   , stat=alloc_err)
allocate(GPTM(NC,NC)  , GPVM(NC,NC)   , stat=alloc_err)
allocate(GPD2VM(NC,NC), GPD4VM(NC,NC) , stat=alloc_err)
allocate(PVEL(NC,NP)                  , stat=alloc_err)
allocate(DTMb(1,NC)                   , stat=alloc_err)

TM        = 0.0_dp
DTM       = 0.0_dp
D2TM      = 0.0_dp
DTMb      = 0.0_dp
VM        = 0.0_dp
DVM       = 0.0_dp
D2VM      = 0.0_dp
D3VM      = 0.0_dp
PTM       = 0.0_dp
PDTM      = 0.0_dp
PDVM      = 0.0_dp
PD2TM     = 0.0_dp
GPVM      = 0.0_dp
GPD2VM    = 0.0_dp
GPD4VM    = 0.0_dp
PVEL      = 0.0_dp

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

! Get the Chebyshev Galerkin functions
call makeVTM(VM,TM,DVM,DTM,DTMb,D2VM,D3VM,Cjm,Pmj, NC,NP)
D2TM = Cjm
y(:,1) = Cjm(:,2)
write(*,*) y(40,1)-y(41,1)
call write_out_mat(y, "y")

! Get the projected Galerkin functions
call projectVT(GPVM,GPD2VM,GPD4VM,GPTM,PVM,PDVM,PTM,PDTM,PD2TM,PVEL, Pmj,Cjm,VM,DVM,D2VM,TM,DTM)

! Create some initial conditions or read in initial values
if (read_ICs) then
   ! Read the fields in (in physical space)
   open(unit=2, file="TICs.txt")
   open(unit=3, file="VICs.txt")
   read(2,*) ((Tyx(jj,ii), ii=1,NF), jj=1,NP)
   read(3,*) ((Vyx(jj,ii), ii=1,NF), jj=1,NP)
   close(2)
   close(3)
   ! Now start bringing them to Chebyshev-Fourier space
   Bmx = matmul(Pmj, Tyx) ! y to Cheb
   Amx = matmul(PVEL, Vyx) ! y to Cheb
   call fftw_execute_dft_r2c(pf1, Bmx, fft1_ml) ! x to Fourier
   Bml = fft1_ml / real(NF, dp)
   call fftw_execute_dft_r2c(pf1, Amx, fft1_ml) ! x to Fourier
   Aml = fft1_ml / real(NF, dp)
   ! Need to compute u (horizontal velocity)
   call update_u(Pmj, DVM)
else 
   amp = -0.03_dp * (pi / 2.0_dp)**2.0_dp
   dx = 2.0_dp * pi / (alpha * NF)
   ! Set ICs in physical space
   do jj = 1,NP
      do ii = 1,NF
         x = real(ii-1, kind=dp)*dx - pi/alpha
         Tyx(jj,ii) = amp*cos(pi*y(jj,1)/2.0_dp)*cos(alpha*x)
         Uyx(jj,ii) = 0 !amp*cos(pi*y(jj,1)/2.0_dp)
      end do
   end do
   ! Bring to Cheb.-Fourier space
   Bmx = matmul(Pmj,Tyx)
   call fftw_execute_dft_r2c(pf1, Bmx, fft1_ml) ! x to Fourier
   Bml = fft1_ml / real(NF, dp)

   Umx = matmul(Pmj,Uyx)
   call fftw_execute_dft_r2c(pf1, Umx, fft1_ml) ! x to Fourier
   Uml = fft1_ml / real(NF, dp)
   ! Write solution to file
   !call initial_conditions(Pmj,y,amp,NC)
end if

! Call time-integrator
call imex_rk(NC, NF, dt, t_final, nu, kappa,    &
             PVEL, Pmj, VM,TM, DVM, DTM,DTMb,D2TM,D2VM, D3VM, &
             GPVM,GPTM,PVM,PDVM,PTM,GPD2VM,GPD4VM,mirror_symmetry)

!call backward_euler(NC,NF,dt,t_final,nu,kappa,        &
!                    PVEL,Pmj,VM,TM,DVM,DTM,DTMb,D2TM,D2VM,D3VM, &
!                    GPVM,GPTM,PVM,PTM,GPD2VM,GPD4VM)

call Nusselt(Nuss, DTMb(1,:), real(Bml(:,1)), NC)

!1000 format(E25.16E3                           )
!2000 format(E25.16E3,E25.16E3                  )
!3000 format(E25.16E3,E25.16E3,E25.16E3         )
!4000 format(E25.16E3,E25.16E3,E25.16E3,E25.16E3)

end program TC

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
