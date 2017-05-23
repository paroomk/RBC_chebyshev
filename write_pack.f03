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

subroutine write_mat_cmplx(mat, aname)

character(len=*), intent(in) :: aname
complex(dp), intent(in) :: mat(:,:)
integer :: ii, jj, mx,my

mx = size(mat(1,:))
my = size(mat(:,1))

open(unit=8000, file=aname//".txt", action="write", status="unknown", position="append")
do jj = 1,my
   do ii = 1,mx
      write(unit=8000,fmt=*) mat(jj,ii)
   end do
end do
close(unit=8000)

end subroutine write_mat_cmplx

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine write_out_mat(arr, aname)

real(dp), dimension(:,:), intent(in) :: arr
character(len=*),         intent(in) :: aname
integer                              :: rows, cols, i, j

rows = size(arr,1)
cols = size(arr,2)

open(unit=8000, file=aname//".txt", action="write", status="unknown", position="append")
do i = 1,rows
   write(unit=8000, fmt=*) (arr(i,j), j=1,cols)
end do
close(unit=8000)

end subroutine write_out_mat

end module write_pack
