module global

use types, only: dp
use fftw

implicit none
save
private
public set_imex_params, fft_utils,       &
       alloc_err, Aml, Bml, Uml,         &
       Vyx, Amx, Bmx, Umx, Tyx, Uyx,     &
       eye, kx, pU, ipU, pV, ipV, pT,    &
       ipT, pf1, ipf1, fft1_ml, fft1_mx, &
       pfy, ipfy, fft1_yl, fft1_yx

! Utilities
integer                                  :: alloc_err

!! Some special global variables
!  The fields (temperature and velocity)
complex(dp), pointer    , dimension(:,:) :: Aml, Bml, Uml ! Cheb-Fourier coeffs
                                                          ! Aml: vertical velocity
                                                          ! Bml: Temperature
real(dp)   , pointer    , dimension(:,:) :: Amx, Bmx, Umx ! Cheb-Phys coeffs
                                                          ! Amx: vertical velocity
                                                          ! Bmx: Temperature
real(dp)   , allocatable, dimension(:,:) :: Vyx, Tyx, Uyx ! Physical space fields

! Nonlinear fields (temperature)
complex(dp), pointer    , dimension(:,:) :: fft1_ml ! Cheb-Fourier coeffs
real(dp)   , pointer    , dimension(:,:) :: fft1_mx ! Cheb-phys coeffs

! Nonlinear fields (velocity)
complex(dp), pointer    , dimension(:,:) :: fft1_yl ! Cheb-Fourier coeffs
real(dp)   , pointer    , dimension(:,:) :: fft1_yx ! Cheb-phys coeffs

! Basic arrays used by everything
real(dp)   , allocatable, dimension(:,:) :: eye           ! Identity matrix
real(dp)   , allocatable, dimension(:)   :: kx            ! Wavenumbers

! FFT plans
type(C_PTR)                              :: pU, ipU   ! plan and inverse plan
type(C_PTR)                              :: pV, ipV   ! plan and inverse plan
type(C_PTR)                              :: pT, ipT   ! plan and inverse plan
type(C_PTR)                              :: pf1, ipf1 ! plan and inverse plan
type(C_PTR)                              :: pfy, ipfy ! plan and inverse plan

contains

subroutine set_imex_params(c1,c2,c3, d11,d21,d31, d12,d22,d32, d13,d23,d33, &
                           dh11,dh12,dh13,dh14, dh21,dh22,dh23,dh24,        &
                           dh31,dh32,dh33,dh34, dh41,dh42,dh43,dh44)

   implicit none

   real(dp) :: gmma

   ! Time-integration parameters
   real(dp) :: c1, c2, c3
   real(dp) :: d11, d12, d13
   real(dp) :: d21, d22, d23
   real(dp) :: d31, d32, d33
   real(dp) :: dh11, dh12, dh13, dh14
   real(dp) :: dh21, dh22, dh23, dh24
   real(dp) :: dh31, dh32, dh33, dh34
   real(dp) :: dh41, dh42, dh43, dh44


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

   ! dhat_ij
   dh11 = 0.0_dp
   dh21 = 0.4358665215_dp
   dh31 = 0.3212788860_dp
   dh41 = -0.105858296_dp

   dh12 = 0.0_dp
   dh22 = 0.0_dp
   dh32 = 0.3966543747_dp
   dh42 = 0.5529291479_dp

   dh13 = 0.0_dp
   dh23 = 0.0_dp
   dh33 = 0.0_dp
   dh43 = 0.5529291479_dp

   dh14 = 0.0_dp
   dh24 = 0.0_dp
   dh34 = 0.0_dp
   dh44 = 0.0_dp

end subroutine set_imex_params

subroutine fft_utils(NP,NC,NF)

   use fftw

   implicit none

   integer, intent(in) :: NP, NC, NF

   ! Prepare FFT variables and parameters
   integer                            :: howmany ! how many FFTs to do per array
   integer, allocatable, dimension(:) :: Nxarr_r, Nxarr_c ! Embed parameters
   integer, allocatable, dimension(:) :: Nxarr_ry, Nxarr_cy ! Embed parameters

   ! Allocate for ffts

   ! Create plan for forward FFT
   pU   = fftw_alloc_real(int(NC*NF, C_SIZE_T))
   pV   = fftw_alloc_real(int(NC*NF, C_SIZE_T))
   pT   = fftw_alloc_real(int(NC*NF, C_SIZE_T))
   pf1  = fftw_alloc_real(int(NC*NF, C_SIZE_T))
   pfy  = fftw_alloc_real(int(NP*NF, C_SIZE_T))

   ! Allocate variables and associate with plan
   call c_f_pointer(pU,   Umx    , [NC,NF])
   call c_f_pointer(pV,   Amx    , [NC,NF])
   call c_f_pointer(pT,   Bmx    , [NC,NF])
   call c_f_pointer(pf1,  fft1_mx, [NC,NF])
   call c_f_pointer(pfy,  fft1_yx, [NP,NF])

   ! Create plan for backward FFT
   ipU  = fftw_alloc_complex(int(NC*(NF/2+1), C_SIZE_T))
   ipV  = fftw_alloc_complex(int(NC*(NF/2+1), C_SIZE_T))
   ipT  = fftw_alloc_complex(int(NC*(NF/2+1), C_SIZE_T))
   ipf1 = fftw_alloc_complex(int(NC*(NF/2+1), C_SIZE_T))
   ipfy = fftw_alloc_complex(int(NP*(NF/2+1), C_SIZE_T))

   ! Allocate variables and associate with plan
   call c_f_pointer(ipU,  Uml    , [NC, NF/2+1])
   call c_f_pointer(ipV,  Aml    , [NC, NF/2+1])
   call c_f_pointer(ipT,  Bml    , [NC, NF/2+1])
   call c_f_pointer(ipf1, fft1_ml, [NC, NF/2+1])
   call c_f_pointer(ipfy, fft1_yl, [NP, NF/2+1])

   allocate(Nxarr_r(NC),  Nxarr_c(NC),  stat=alloc_err)
   allocate(Nxarr_ry(NP), Nxarr_cy(NP), stat=alloc_err)
   Nxarr_r  = NF
   Nxarr_c  = NF/2+1
   Nxarr_ry = NF
   Nxarr_cy = NF/2+1

   howmany=NC

   ! Plans for standard fields (velocity and temperature)

   pU  = fftw_plan_many_dft_r2c(1, Nxarr_r, howmany, Umx, Nxarr_r, howmany, 1,&
                              &                      Uml, Nxarr_c, howmany, 1,&
                              & FFTW_MEASURE)

   ipU = fftw_plan_many_dft_c2r(1, Nxarr_r, howmany, Uml, Nxarr_c, howmany, 1,&
                             &                       Umx, Nxarr_r, howmany, 1,&
                             & FFTW_MEASURE)

   pV  = fftw_plan_many_dft_r2c(1, Nxarr_r, howmany, Amx, Nxarr_r, howmany, 1,&
                              &                      Aml, Nxarr_c, howmany, 1,&
                              & FFTW_MEASURE)

   ipV = fftw_plan_many_dft_c2r(1, Nxarr_r, howmany, Aml, Nxarr_c, howmany, 1,&
                             &                       Amx, Nxarr_r, howmany, 1,&
                             & FFTW_MEASURE)

   pT  = fftw_plan_many_dft_r2c(1, Nxarr_r, howmany, Bmx, Nxarr_r, howmany, 1,&
                              &                      Bml, Nxarr_c, howmany, 1,&
                              & FFTW_MEASURE)

   ipT = fftw_plan_many_dft_c2r(1, Nxarr_r, howmany, Bml, Nxarr_c, howmany, 1,&
                             &                       Bmx, Nxarr_r, howmany, 1,&
                             & FFTW_MEASURE)


!  Generic plans
   pf1  = fftw_plan_many_dft_r2c(1, Nxarr_r, howmany, fft1_mx, Nxarr_r, howmany, 1,&
                               &                      fft1_ml, Nxarr_c, howmany, 1,&
                               & FFTW_MEASURE)

   ipf1 = fftw_plan_many_dft_c2r(1, Nxarr_r, howmany, fft1_ml, Nxarr_c, howmany, 1,&
                               &                      fft1_mx, Nxarr_r, howmany, 1,&
                               & FFTW_MEASURE)

!  Generic plans
   howmany = NP
   pfy  = fftw_plan_many_dft_r2c(1, Nxarr_ry, howmany, fft1_yx, Nxarr_ry, howmany, 1,&
                               &                       fft1_yl, Nxarr_cy, howmany, 1,&
                               & FFTW_MEASURE)

   ipfy = fftw_plan_many_dft_c2r(1, Nxarr_ry, howmany, fft1_yl, Nxarr_cy, howmany, 1,&
                               &                       fft1_yx, Nxarr_ry, howmany, 1,&
                               & FFTW_MEASURE)

end subroutine fft_utils

end module global
