

subroutine calc_force(pos, force, N,L)
  implicit none
  real(8), intent(in) :: L
  integer, intent(in) :: N 
  real(8), intent(in) :: pos(N, 3)
  real(8), intent(inout) :: force(N, 3)
!f2py intent(in, out) force        
  real(8) :: delta_r(3), dr2, F
  real(8), parameter :: rmax = 3.2_8
  integer :: i, j
  force = 0._8
  !print *, 'piep'
  !PRINT '(3F12.6)', pos(1:3,:)
  do i = 1, N
    do j = 1, i-1
      delta_r = pos(i,:) - pos(j,:)
      delta_r = delta_r - nint(delta_r/L)*L
      dr2 = sum(delta_r**2)
      if (dr2<rmax**2) then
        F=24*(2/dr2**7 - 1/dr2**4)
        if (F>1000._8) print *, i, j, dr2
        Force(i,:) = Force(i,:) + F*delta_r 
        Force(j,:) = Force(j,:) - F*delta_r
      end if
    end do
!    print *, SUM(FORCE)
  end do
!  print *, 'haha', sum(force)
end subroutine
