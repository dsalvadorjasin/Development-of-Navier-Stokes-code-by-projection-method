subroutine dump_fields()
! Debug/validation dump of the raw field arrays at full double precision.
! Only active when the environment variable NS_DUMP_FIELDS is set; writes
! fields_<istep>.txt at steps 1,2,3 and every 'isto' steps.
use variables
implicit none

character(len=32) :: env
character(len=32) :: fname
integer :: ienv

call get_environment_variable('NS_DUMP_FIELDS', env, status=ienv)
if (ienv /= 0) return
if (.not. (istep <= 3 .or. mod(istep,isto) == 0)) return

write(fname,'(A,I6.6,A)') 'fields_', istep, '.txt'
open(12, file=fname, status='replace')
write(12,'(A,I8)') 'istep ', istep
write(12,'(A,ES25.17E3)') 'time ', time
write(12,'(A,ES25.17E3)') 'dt ', dt

write(12,'(A,4I6)') 'u ', 0, nx+1, 0, ny+1
do j=0,ny+1
   do i=0,nx+1
      write(12,'(ES25.17E3)') u(i,j)
   end do
end do

write(12,'(A,4I6)') 'v ', 0, nx+1, 0, ny+1
do j=0,ny+1
   do i=0,nx+1
      write(12,'(ES25.17E3)') v(i,j)
   end do
end do

write(12,'(A,4I6)') 'p ', 0, nx+1, 0, ny+1
do j=0,ny+1
   do i=0,nx+1
      write(12,'(ES25.17E3)') p(i,j)
   end do
end do

write(12,'(A,4I6)') 'u_star ', 0, nx+1, 0, ny+1
do j=0,ny+1
   do i=0,nx+1
      write(12,'(ES25.17E3)') u_star(i,j)
   end do
end do

write(12,'(A,4I6)') 'v_star ', 0, nx+1, 0, ny+1
do j=0,ny+1
   do i=0,nx+1
      write(12,'(ES25.17E3)') v_star(i,j)
   end do
end do

write(12,'(A,4I6)') 'rot ', 1, nx, 1, ny
do j=1,ny
   do i=1,nx
      write(12,'(ES25.17E3)') rot(i,j)
   end do
end do

write(12,'(A,4I6)') 'div ', 1, nx, 1, ny
do j=1,ny
   do i=1,nx
      write(12,'(ES25.17E3)') div(i,j)
   end do
end do

close(12)

end subroutine dump_fields
