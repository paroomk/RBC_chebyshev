module chebutils

use types,  only: dp, pi
use global, only: alloc_err
use write_pack

implicit none

contains

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine makeTP(Tjm, Pmj, NCin, NPin)

  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  !
  ! Returns Chebyshev expansion and projection matrices.
  !    u = Tjm*a, a = Pmj*u
  ! INPUT:
  !    NC:           Number of Chebyshev polynomials
  !    NP (optional) Number of Guass points
  ! OUTPUT:
  !    Tjm: Cheb to Phys transformation matrix
  !    Pmj: Phys to Cheb projection matrix
  !
  ! Adapted from Fabian Waleffe's Matlab code (1999-2005)
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

  integer, intent(in)                                :: NCin
  integer, optional, intent(in)                      :: NPin
  real(dp), allocatable, dimension(:,:), intent(out) :: Tjm, Pmj
  real(dp), allocatable, dimension(:,:)              :: mtj
  integer                                            :: j, m
  integer                                            :: NC, NP

  NC = NCin ! Chebyshev polys

  ! If only one argument is present then set NP = NC
  if (present(NPin)) then
     NP = NPin
  else
     NP = NC
  end if

  ! Guard against user error.  Must have NP >= NC
  if (NP < NC) then
     write(*,*) "Warning:  NP < NC, Setting NP = NC"
     NP = NC
  end if

  ! Allocate transformation and projection arrays
  allocate(Tjm(NP,NC), stat=alloc_err)
  allocate(mtj(NP,NC), stat=alloc_err)
  allocate(Pmj(NC,NP), stat=alloc_err)

  Tjm = 0.0_dp
  mtj = 0.0_dp
  Pmj = 0.0_dp

  ! Form transformation matrix (the following uses all of Wally's ``tricks'')
  do m = 1,NC
     do j = 1,NP
        mtj(j,m) = mod((2*j-1)*(m-1), 4*NP)     ! stay in [0,2pi], 0<mtj<4NP
        mtj(j,m) = 2*NP - abs(2*NP-mtj(j,m))    ! stay in [0, pi], 0<mtj<2NP
        Tjm(j,m) = sin(pi*(NP-mtj(j,m))/(2*NP)) ! transformation matrix (cheb to phys)
     end do
  end do

  ! Now form projection matrix (phys to Cheb)
  Pmj      = 2.0_dp*transpose(Tjm)/NP
  Pmj(1,:) = Pmj(1,:)/2.0_dp ! Fix the first row

end subroutine makeTP

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine chebint(I1T, I0T)

  !:::::::::::::::::::::::::::::::::::::::::::::::::::
  !
  ! CHEB to PHYS integration operator, column-wise
  !
  ! INPUT:
  !   I0T:  Ny X NC, Columns are Cheb polys or kth integral of Cheb polys
  ! OUTPUT:
  !   I1T:  Ny X NC-1, The next integral, column by column
  !         Has 1 fewer column because IT_M cannot be computing w/o T_{M+1}
  !
  ! Note that integration constants will be fixec by the BCs
  !
  ! This code was adapted from Fabian Waleffe's Matlab code (1999, 2003, 2014) 
  !
  !:::::::::::::::::::::::::::::::::::::::::::::::::::

  real(dp), dimension(:,:)              :: I0T
  real(dp), allocatable, dimension(:,:) :: I1T
  integer                               :: ny, nc
  integer                               :: m, mm
  integer, dimension(2)                 :: sz

  ! Get the size of the array to integrate
  sz = shape(I0T)
  ny = sz(1)
  nc = sz(2)

  ! Allocate integration matrix
  allocate(I1T(ny,nc-1), stat=alloc_err)
  I1T = 0.0_dp

  ! Integrals of first two polys
  I1T(:,1) = I0T(:,2)
  I1T(:,2) = I0T(:,3)/4.0_dp ! Must match C2Cint.m for some apps

  ! Integrals of the rest of the polys
  do m = 2,nc-2
     mm = m+1
     I1T(:,mm) = I0T(:,mm+1)/(2.0*(m+1)) - I0T(:,mm-1)/(2.0*(m-1))
  end do

end subroutine chebint

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine makeVTM(VM, TM, DVM, DTM, DTMb, D2VM, D3VM, Cjm, Pmj, NCin,NPin)

  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  !
  ! Returns Chebyshev expansion and projection matrices.
  !         Also provides matrix of temperature modes and
  !         their derivatives.
  !         u = Tjm*a, a = Pmj*u
  !         T(y) = TM(y,m) b(m)
  ! INPUT:
  !    NC:           Number of Chebyshev polynomials
  !    NP (optional) Number of Guass points
  ! OUTPUT:
  !    Cjm : Cheb to Phys transformation matrix
  !    Pmj : Phys to Cheb projection matrix
  !    VM  : Trial functions for vertical velocity (NP x NC)
  !    TM  : Trial functions for temperature (NP x NC)
  !    DVM : Trial functions for derivative of vertical-velocity (NP x NC)
  !    DTM : Trial functions for derivative of temperature (NP x NC)
  !    D2VM: Trial functions for 2nd derivative of vertical-velocity (NP x NC)
  !    D3VM: Trial functions for 3rd derivative of vertical-velocity (NP x NC)
  !
  ! Adapted from Fabian Waleffe's Matlab code (1999-2005)
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

  integer,                               intent(in)  :: NCin
  integer, optional,                     intent(in)  :: NPin
  real(dp), allocatable, dimension(:,:), intent(out) :: TM, VM 
  real(dp), allocatable, dimension(:,:), intent(out) :: DVM, DTM, DTMb
  real(dp), allocatable, dimension(:,:), intent(out) :: D2VM, D3VM
  real(dp), allocatable, dimension(:,:), intent(out) :: Cjm, Pmj
  real(dp), allocatable, dimension(:,:)              :: I1T, I2T, I3T, I4T
  real(dp), allocatable, dimension(:,:)              :: Cjmi, Pmji
  real(dp), allocatable, dimension(:,:)              :: Tpm1, I1Tpm1
  real(dp), allocatable, dimension(:,:)              :: I2Tpm1, I3Tpm1, I4Tpm1
  real(dp), allocatable, dimension(:,:)              :: IbcV, IbcT
  real(dp), allocatable, dimension(:,:)              :: cT0, cT1
  real(dp), allocatable, dimension(:,:)              :: cV0, cV1, cV2, cV3
  real(dp), allocatable, dimension(:,:)              :: y, y2, y3, ones
  real(dp),              dimension(4,4)              :: CM1
  real(dp),              dimension(2,2)              :: CM2
  real(dp)                                           :: at, ab, bt, bb
  integer                                            :: i, j
  integer                                            :: NC, NP
  integer                                            :: info
  integer, allocatable, dimension(:)                 :: ipiv

  NC = NCin ! Chebyshev polys

  ! If only one argument is present then set NP = NC
  if (present(NPin)) then
     NP = NPin
  else
     NP = NC+4
  end if

  ! Guard against user error.  Must have NP >= NC+2 for integration purposes.
  if (NP < NC+4) then
     write(*,*) "Warning:  NP < NC+4, Setting NP = NC+4"
     NP = NC+4
  end if

  ! Allocate y-vectors
  allocate(y(NP,1)   , stat=alloc_err)
  allocate(y2(NP,1)  , stat=alloc_err)
  allocate(y3(NP,1)  , stat=alloc_err)
  allocate(ones(NP,1), stat=alloc_err)
  y = 0.0_dp
  y2 = 0.0_dp
  y3 = 0.0_dp
  ones(:,1) = 1.0_dp

  ! Allocate transformation and projection arrays
  allocate(Cjm(NP,NC), stat=alloc_err)
  allocate(Pmj(NC,NP), stat=alloc_err)

  allocate(Cjmi(NP,NC+4), stat=alloc_err)
  allocate(Pmji(NC+4,NP), stat=alloc_err)

  ! Allocate Galerkin trial functions
  allocate(TM(NP,NC)  , DTM(NP,NC) , stat=alloc_err)
  allocate(VM(NP,NC)  , DVM(NP,NC) , stat=alloc_err)
  allocate(D2VM(NP,NC), D3VM(NP,NC), stat=alloc_err)
  allocate(DTMb(1,NC)              , stat=alloc_err)

  TM   = 0.0_dp
  DTM  = 0.0_dp
  VM   = 0.0_dp
  DVM  = 0.0_dp
  D2VM = 0.0_dp
  D3VM = 0.0_dp
  DTMb = 0.0_dp

  Cjm  = 0.0_dp
  Pmj  = 0.0_dp

  Cjmi = 0.0_dp
  Pmji = 0.0_dp

  ! Allocate integration matrices
  allocate(I1T(NP,NC+3), stat=alloc_err)
  allocate(I2T(NP,NC+2), stat=alloc_err)
  allocate(I3T(NP,NC+1), stat=alloc_err)
  allocate(I4T(NP,NC)  , stat=alloc_err)

  I1T = 0.0_dp
  I2T = 0.0_dp
  I3T = 0.0_dp
  I4T = 0.0_dp

  ! Allocate boundary matrices
  allocate(Tpm1(2,NC+4)        , stat=alloc_err)
  allocate(I1Tpm1(2,NC+3)      , stat=alloc_err)
  allocate(I2Tpm1(2,NC+2)      , stat=alloc_err)
  allocate(I3Tpm1(2,NC+1)      , stat=alloc_err)
  allocate(I4Tpm1(2,NC)        , stat=alloc_err)
  allocate(IbcV(4,NC)          , stat=alloc_err)
  allocate(IbcT(2,NC)          , stat=alloc_err)
  allocate(cT0(1,NC), cT1(1,NC), stat=alloc_err)
  allocate(cV0(1,NC), cV1(1,NC), stat=alloc_err)
  allocate(cV2(1,NC), cV3(1,NC), stat=alloc_err)

  Tpm1   = 0.0_dp
  I1Tpm1 = 0.0_dp
  I2Tpm1 = 0.0_dp
  I3Tpm1 = 0.0_dp
  I4Tpm1 = 0.0_dp
  CM1    = 0.0_dp
  CM2    = 0.0_dp
  IbcV   = 0.0_dp
  IbcT   = 0.0_dp
  cT0    = 0.0_dp
  cT1    = 0.0_dp
  cV0    = 0.0_dp
  cV1    = 0.0_dp
  cV2    = 0.0_dp
  cV3    = 0.0_dp

  ! Allocate some LAPACK and BLAS arrays
  allocate(ipiv(NC), stat=alloc_err)

  ipiv = 0

  ! Boundary conditions
  ! rigid-rigid: at = ab = 1.0
  ! free-free  : at = ab = 0.0
  ! rigid-free : at = 0.0, ab = 1.0
  at = 1.0_dp
  ab = 1.0_dp
  bt = 1.0_dp - at
  bb = 1.0_dp - ab

  ! Make the expansion (Cjm) and projection (Pmj) matrices
  call makeTP(Cjmi, Pmji, NC+4, NP)

  ! y-coordinate
  y(:,1) = Cjmi(:,2)

  ! Create integration matrices
  call chebint(I1T, Cjmi)
  call chebint(I2T, I1T)
  call chebint(I3T, I2T)
  call chebint(I4T, I3T)

  ! Now take care of boundary conditions
  do i = 1,NC+4
     Tpm1(1,i) =   1.0_dp
     Tpm1(2,i) = (-1.0_dp)**(i-1)
  end do

  ! Integrals at the boundaries
  call chebint(I1Tpm1, Tpm1)
  call chebint(I2Tpm1, I1Tpm1)
  call chebint(I3Tpm1, I2Tpm1)
  call chebint(I4Tpm1, I3Tpm1)

! Stokes equation

  ! Get matrix of boundary terms
  CM1(1,1) = Tpm1(1,1)
  CM1(2,1) = Tpm1(2,1)
  CM1(3,1) = 0.0_dp
  CM1(4,1) = 0.0_dp

  CM1(1,2) = I1Tpm1(1,1)
  CM1(2,2) = I1Tpm1(2,1)
  CM1(3,2) = at*Tpm1(1,1)
  CM1(4,2) = ab*Tpm1(2,1)

  CM1(1,3) = I2Tpm1(1,1)
  CM1(2,3) = I2Tpm1(2,1)
  CM1(3,3) = at*I1Tpm1(1,1) + bt*Tpm1(1,1)
  CM1(4,3) = ab*I1Tpm1(2,1) + bb*Tpm1(2,1)

  CM1(1,4) = I3Tpm1(1,1)
  CM1(2,4) = I3Tpm1(2,1)
  CM1(3,4) = at*I2Tpm1(1,1) + bt*I1Tpm1(1,1)
  CM1(4,4) = ab*I1Tpm1(2,1) + bb*I1Tpm1(2,1)

  ! Form RHS of boundary terms
  IbcV(1,:) = I4Tpm1(1,1:NC) ! v(1)
  IbcV(2,:) = I4Tpm1(2,1:NC) ! v(-1)
  IbcV(3,:) = at*I3Tpm1(1,1:NC) + bt*I2Tpm1(1,1:NC) ! at*Dv(1) + (1-at)*D2v(1)
  IbcV(4,:) = ab*I3Tpm1(2,1:NC) + bb*I2Tpm1(2,1:NC) ! ab*Dv(-1) + (1-ab)*D2v(-1)

  ! Solve for boundary matrix: B = inv(CM1)*Ibc
  call dgesv(4, NC, -CM1, 4, ipiv, IbcV, 4, info)

  cV0(1,:) = IbcV(1,:)
  cV1(1,:) = IbcV(2,:)
  cV2(1,:) = IbcV(3,:)
  cV3(1,:) = IbcV(4,:)

  ! Create Galerkin trial functions for velocity
  y2(:,1) = I2T(:,1)
  y3(:,1) = I3T(:,1)
  D3VM = I1T(:,1:NC) + matmul(ones, cV3)
  D2VM = I2T(:,1:NC) + matmul(ones, cV2) + matmul(y, cV3)
  DVM  = I3T(:,1:NC) + matmul(ones, cV1) + matmul(y, cV2) + &
        &              matmul(y2, cV3)
  VM   = I4T(:,1:NC) + matmul(ones, cV0) + matmul(y, cV1) + &
        &              matmul(y2, cV2) + matmul(y3, cV3)

! Temperature equation

  ! Get matrix of boundary terms
  CM2 = reshape([Tpm1(1,1)  , Tpm1(2,1)   , & ! 1st col
               & I1Tpm1(1,1), I1Tpm1(2,1)]  & ! 2nd col
       &,shape(CM2))

  ! Form RHS of boundary terms
  IbcT(1,:) = I2Tpm1(1,1:NC) ! T(1)
  IbcT(2,:) = I2Tpm1(2,1:NC) ! T(-1)

  ! Solve for boundary matrix: B = inv(CM1)*Ibc
  call dgesv(2, NC, -CM2, 2, ipiv, IbcT, 2, info)

  cT0(1,:) = IbcT(1,:)
  cT1(1,:) = IbcT(2,:)

  DTM  = I1T(:,1:NC) + matmul(ones, cT1)
  TM   = I2T(:,1:NC) + matmul(ones, cT0) + matmul(y, cT1)

  DTMb(1,:) = I1Tpm1(1,1:NC) + cT1(1,:) ! Need this for Nu calculation

  Cjm = Cjmi(:,1:NC)
  Pmj = Pmji(1:NC,:)

end subroutine makeVTM

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine projectVT(GPVM,GPD2VM,GPD4VM,GPTM,PVM,PDVM,PTM,PDTM,PD2TM,PVEL, Pmj,Cjm,VM,DVM,D2VM,TM,DTM)

  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
  !
  ! Returns projected Chebyshev-Galerkin trial functions
  !
  ! INPUT:
  !    Pmj :  Chebyshev projection matrix (NC X NP)
  !    Cjm :  Chebyshev expansion matrix (NP X NC)
  !    VM  :  Chebyshev-Galerkin trial function for velocity    (NP X NC)
  !    TM  :  Chebyshev-Galerkin trial function for temperature (NP X NC)
  !    DVM :  Chebyshev-Galerkin trial function for first derivative of velocity
  !           (NP X NC)
  !    D2VM:  Chebyshev-Galerkin trial function for second derivative of
  !           velocity (NP X NC)
  !    DTM :  Chebyshev-Galerkin trial function for first derivative of
  !           temperature (NP X NC)
  ! OUTPUT:
  !    GPVM  : Projected trial function for velocity (NC X NC)
  !    GPD2VM: Projected trial function for 2nd derivative of velocity (NC X NC)
  !    GPD4VM: Projected trial function for 4th derivative of velocity (NC X NC)
  !    GPTM  : Projected trial function for temperature in v-equation (NC X NC)
  !    PVM   : Projected trial function for velocity in temperature equation (NC X NC)
  !    PDVM  : Projected trial function for first derivative of velocity (NC X NC)
  !    PTM   : Projected trial function for temperature (NC X NC)
  !    PDTM  : Projected trial function for first derivative of temperature (NC X NC)
  !    PD2TM : Projected trial function for 2nd derivative of temperature (NC X NC)
  !    PVEL  : Galerkin projector for vertical-velocity equation (NC X NP)
  !
  ! Adapted from Fabian Waleffe's Matlab code (1999-2005)
  !:::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

  real(dp),              dimension(:,:), intent(in)  :: Pmj, Cjm
  real(dp),              dimension(:,:), intent(in)  :: VM, D2VM, TM
  real(dp),              dimension(:,:), intent(in)  :: DVM, DTM
  real(dp), allocatable, dimension(:,:), intent(out) :: GPVM, GPD2VM
  real(dp), allocatable, dimension(:,:), intent(out) :: GPD4VM, GPTM
  real(dp), allocatable, dimension(:,:), intent(out) :: PVM, PTM, PD2TM, PVEL
  real(dp), allocatable, dimension(:,:), intent(out) :: PDVM, PDTM
  real(dp), allocatable, dimension(:)                :: y
  integer                                            :: NC, NP
  integer                                            :: j
  real(dp), parameter                                :: scale1=1.0_dp, scale2=0.0_dp

  NC = size(Pmj,1)
  NP = size(Pmj,2)

  allocate(PVM(NC,NC), PTM(NC,NC), PD2TM(NC,NC), stat=alloc_err)
  allocate(PDVM(NC,NC), PDTM(NC,NC)            , stat=alloc_err)
  allocate(GPTM(NC,NC)  , GPVM(NC,NC)          , stat=alloc_err)
  allocate(GPD2VM(NC,NC), GPD4VM(NC,NC)        , stat=alloc_err)
  allocate(PVEL(NC,NP)                         , stat=alloc_err)
  allocate(y(NP)                               , stat=alloc_err)
  PVM    = 0.0_dp
  PDVM   = 0.0_dp
  PTM    = 0.0_dp
  PDTM   = 0.0_dp
  PD2TM  = 0.0_dp
  GPTM   = 0.0_dp
  GPVM   = 0.0_dp
  GPD2VM = 0.0_dp
  GPD4VM = 0.0_dp
  PVEL   = 0.0_dp
  y      = 0.0_dp

  y = Cjm(:,2)
  ! Not enough to simply project with Pmj.  Need to create a special projector
  ! to avoid spurious positive eigenvalues.
  y = 1.0_dp - y**2.0_dp
  do j = 1,NC
     PVEL(j,:) = y
  end do
  PVEL = Pmj*PVEL

  ! Galerkin projectors for the vertical velocity equation
  call dgemm('n', 'n', NC, NC, NP, scale1, PVEL, NC, TM  , NP, scale2, GPTM  , NC)
  call dgemm('n', 'n', NC, NC, NP, scale1, PVEL, NC, VM  , NP, scale2, GPVM  , NC)
  call dgemm('n', 'n', NC, NC, NP, scale1, PVEL, NC, D2VM, NP, scale2, GPD2VM, NC)
  call dgemm('n', 'n', NC, NC, NP, scale1, PVEL, NC, Cjm , NP, scale2, GPD4VM, NC)

  ! Projectors for temperature equation and continuity equation
  call dgemm('n', 'n', NC, NC, NP, scale1, Pmj , NC, VM   , NP, scale2, PVM , NC)
  call dgemm('n', 'n', NC, NC, NP, scale1, Pmj , NC, TM   , NP, scale2, PTM , NC)
  call dgemm('n', 'n', NC, NC, NP, scale1, Pmj , NC, DVM  , NP, scale2, PDVM, NC)
  call dgemm('n', 'n', NC, NC, NP, scale1, Pmj , NC, DTM  , NP, scale2, PDTM, NC)

  forall(j = 1:NC) PD2TM(j,j) = 1.0_dp ! Just the identity matrix

end subroutine projectVT

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

end module chebutils
