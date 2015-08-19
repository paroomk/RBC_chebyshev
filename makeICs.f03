module makeICs

use types,  only: dp, pi
use global, only: alloc_err, Aml, Bml

implicit none

contains

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine initial_conditions(Pmj,y,amp,NC)

use write_pack

implicit none

real(dp),              dimension(:,:), intent(in)  :: Pmj
real(dp),              dimension(:,:), intent(in)  :: y
real(dp)                             , intent(in)  :: amp
integer,                               intent(in)  :: NC
integer                                            :: NP
real(dp), allocatable, dimension(:,:)              :: Ty1, D2Ty1

NP = size(y)

allocate(D2Ty1(NP,1), Ty1(Np,1) , stat=alloc_err)
Ty1   = 0.0_dp
D2Ty1 = 0.0_dp

Ty1   = amp*cos(pi*y/2.0_dp)
D2Ty1 = -(pi**2.0_dp/4.0_dp)*Ty1

Bml(:,2) = cmplx(reshape(matmul(Pmj, D2Ty1),[NC]), 0.0_dp)

end subroutine initial_conditions

end module makeICs
