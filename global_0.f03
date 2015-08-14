module global

implicit none
save

! Get double-precision type
integer, parameter     :: dp = kind(1.0d0)

! Create fundamental constants
real(dp), parameter    :: pi = 3.141592653589793238462643383279502884197_dp
complex(dp), parameter :: CI = (0.0_dp, 1.0_dp)

! Utilities
integer                :: alloc_err

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

! Time-integration parameters
real(dp) :: c1, c2, c3
real(dp) :: d11, d12, d13
real(dp) :: d21, d22, d23
real(dp) :: d31, d32, d33

contains

subroutine set_imex_params

   implicit none

   real(dp) :: gmma

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

end module global
