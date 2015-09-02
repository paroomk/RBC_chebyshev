module statutils

use types, only: dp
use global, only:

implicit none

private
public Nusselt

contains

subroutine Nusselt(Nu, DTMb,Tbar,NC)

  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  !
  ! Compute the Nusselt number
  !
  ! INPUT:
  !    DTMb:  Temperature modes for first derivative at top wall
  !    Tbar:  Temperature averaged in x-direction
  !
  ! OUTPUT:
  !    Nu:  Nusselt number
  !
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

  integer               , intent(in) :: NC
  real(dp), dimension(:), intent(in) :: DTMb, Tbar
  real(dp), intent(out)              :: Nu
  real(dp)                           :: ddot ! BLAS dot product

  ! Some BLAS and LAPACK variables
  real(dp), parameter :: scale1 = 1.0_dp
  integer,  parameter :: incx = 1, incy = 1

  Nu = -ddot(NC, DTMb, incx, Tbar, incy) + 1.0_dp

  write(*,*) "Nu = ", Nu

end subroutine Nusselt

end module statutils
