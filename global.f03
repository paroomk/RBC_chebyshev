module global

use types, only: dp
use fftw

implicit none
save
private
public set_imex_params, fft_utils, &
       alloc_err, Aml, Bml, Amx, Bmx, &
       Vyx, Tyx, eye, kx, pV, ipV, pT, ipT

! Utilities
integer                                  :: alloc_err

! Some special global variables
complex(dp), pointer    , dimension(:,:) :: Aml, Bml ! Cheb-Fourier coeffs
                                                     ! Aml: vertical velocity
                                                     ! Bml: Temperature
real(dp)   , pointer    , dimension(:,:) :: Amx, Bmx ! Cheb-Phys coeffs
                                                     ! Amx: vertical velocity
                                                     ! Bmx: Temperature
real(dp)   , allocatable, dimension(:,:) :: Vyx, Tyx ! Physical space fields
real(dp)   , allocatable, dimension(:,:) :: eye      ! Identity matrix
real(dp)   , allocatable, dimension(:)   :: kx       ! Wavenumbers

! FFT plans
type(C_PTR)                              :: pV, ipV   ! plan and inverse plan
type(C_PTR)                              :: pT, ipT   ! plan and inverse plan

contains

subroutine set_imex_params(c1,c2,c3, d11,d21,d31, d12,d22,d32, d13,d23,d33)

   implicit none

   real(dp) :: gmma

   ! Time-integration parameters
   real(dp) :: c1, c2, c3
   real(dp) :: d11, d12, d13
   real(dp) :: d21, d22, d23
   real(dp) :: d31, d32, d33


   gmma = 0.4358665215_dp

   c1 = -3.0_dp*gmma**2.0_dp/2.0_dp + 4.0_dp*gmma - 1.0_dp/4.0_dp
   c2 =  3.0_dp*gmma**2.0_dp/2.0_dp - 5.0_dp*gmma + 5.0_dp/4.0_dp
   c3 =  gmma

   ! d_ij
   d11 = gmma
   d21 = (1.0_dp-gmma)/2.0_dp
   d31 = c1

   d12 = 0.0_dp
   d22 = gmma
   d32 = c2

   d13 = 0.0_dp
   d23 = 0.0_dp
   d33 = c3

end subroutine set_imex_params

subroutine fft_utils(NC,NF)

   use fftw

   implicit none

   integer, intent(in) :: NC, NF

   ! Prepare FFT variables and parameters
   integer                            :: howmany ! how many FFTs to do per array
   integer, allocatable, dimension(:) :: Nxarr_r, Nxarr_c ! Embed parameters

   integer                            :: howmanyT0 ! how many FFTs to do per array
   integer, allocatable, dimension(:) :: Nxarr_rT0, Nxarr_cT0 ! Embed parameters

   ! Allocate for ffts

   ! Create plan for forward FFT
   pV = fftw_alloc_real(int(NC*NF, C_SIZE_T))
   pT = fftw_alloc_real(int(NC*NF, C_SIZE_T))

   ! Allocate variables and associate with plan
   call c_f_pointer(pV, Amx, [NC,NF])
   call c_f_pointer(pT, Bmx, [NC,NF])

   ! Create plan for backward FFT
   ipV = fftw_alloc_complex(int(NC*(NF/2+1), C_SIZE_T))
   ipT = fftw_alloc_complex(int(NC*(NF/2+1), C_SIZE_T))

   ! Allocate variables and associate with plan
   call c_f_pointer(ipV, Aml, [NC, NF/2+1])
   call c_f_pointer(ipT, Bml, [NC, NF/2+1])

   allocate(Nxarr_r(NC), Nxarr_c(NC), stat=alloc_err)
   Nxarr_r = NF
   Nxarr_c = NF/2+1

   howmany=NC

   pV  = fftw_plan_many_dft_r2c(1, Nxarr_r, howmany, Amx, Nxarr_r, howmany, 1,&
                              &                      Aml, Nxarr_c, howmany, 1,&
                              & FFTW_ESTIMATE)

   ipV = fftw_plan_many_dft_c2r(1, Nxarr_r, howmany, Aml, Nxarr_c, howmany, 1,&
                             &                       Amx, Nxarr_r, howmany, 1,&
                             & FFTW_ESTIMATE)

   pT  = fftw_plan_many_dft_r2c(1, Nxarr_r, howmany, Bmx, Nxarr_r, howmany, 1,&
                              &                      Bml, Nxarr_c, howmany, 1,&
                              & FFTW_ESTIMATE)

   ipT = fftw_plan_many_dft_c2r(1, Nxarr_r, howmany, Bml, Nxarr_c, howmany, 1,&
                             &                       Bmx, Nxarr_r, howmany, 1,&
                             & FFTW_ESTIMATE)

end subroutine fft_utils

end module global
