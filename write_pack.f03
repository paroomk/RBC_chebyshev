module write_pack

use types, only: dp

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine write_vec(vname,place, vec)

character(len=*), intent(in) :: vname,place
real(dp), intent(in) :: vec(:)
integer :: ii, n

n = size(vec)

do ii = 1,n
   if (ii == 1) then
      write(*,*) vname//" = ", vec(ii)
   else
      write(*,*) place//"   ", vec(ii)
   end if
end do

end subroutine write_vec

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine write_mat(vname,place, mat)

character(len=*), intent(in) :: vname,place
real(dp), intent(in)         :: mat(:,:)
integer                      :: jj, mx,my

mx = size(mat(1,:))
my = size(mat(:,1))

do jj = 1,my
   if (jj==1) then
      write(*,*) vname//" = ", mat(jj,:)
   else
      write(*,*) place//"   ", mat(jj,:)
   end if
end do

end subroutine write_mat

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine write_vec_cmplx(vname,place,vec)

character(len=*), intent(in) :: vname, place
complex(dp), intent(in) :: vec(:)
integer :: ii, n

n = size(vec)

do ii = 1,n
   if (ii == 1) then
      write(*,*) vname//" = ", vec(ii)
   else
      write(*,*) place//"   ", vec(ii)
   end if
end do

end subroutine write_vec_cmplx

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine write_mat_cmplx(vname,place, mat)

character(len=*), intent(in) :: vname,place
complex(dp), intent(in) :: mat(:,:)
integer :: jj, mx,my

mx = size(mat(1,:))
my = size(mat(:,1))

do jj = 1,my
   if (jj==1) then
      write(*,*) vname//" = ", mat(jj,:)
   else
      write(*,*) place//"   ", mat(jj,:)
   end if
end do

end subroutine write_mat_cmplx

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

end module write_pack
