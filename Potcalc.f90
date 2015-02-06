subroutine calc_pot(pos,pot,N,L)
  implicit none
  real(8), intent(in) :: L
  integer, intent(in) :: N 
  real(8), intent(in) :: pos(N, 3)
  real(8), intent(out) :: pot
  real(8) :: delta_r(3), dr2
  real(8), parameter :: rmax = 3.2_8
  integer :: i, j
  pot = 0._8
  do i = 1, N
    do j = 1, i-1
      delta_r = pos(i,:) - pos(j,:)
      delta_r = delta_r - nint(delta_r/L)*L
      dr2 = sum(delta_r**2)
      if (dr2<rmax**2) then
        pot=pot+4*((1/dr2)**6 - (1/dr2)**3)
        if (pot>1000._8) print *, i, j, dr2
      end if
    end do
  end do
end subroutine
