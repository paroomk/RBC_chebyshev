module types

implicit none
private
public dp, pi, CI

! Get double-precision type
integer, parameter     :: dp = kind(1.0d0)

! Create fundamental constants
real(dp), parameter    :: pi = 3.141592653589793238462643383279502884197_dp
complex(dp), parameter :: CI = (0.0_dp, 1.0_dp)

end module types
