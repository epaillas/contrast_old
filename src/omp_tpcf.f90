program openmp
use OMP_LIB
implicit none
 
integer*4 :: i, a, tn, nt

 
a = 0
!$OMP PARALLEL DO
do i = 1, 10
  tn = OMP_get_thread_num()
  nt = OMP_get_num_threads()
  print* , nt
!$OMP ATOMIC
a = a + i
end do
!$OMP END PARALLEL DO
 
 write(*,*) a
 
 end program openmp