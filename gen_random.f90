program gen_random
  	use precision_mod, only: dp
  	implicit none
  	integer :: n, i, j, unit_out
  	real(dp), allocatable :: A(:,:), B(:)
  	real(dp) :: r

  	write(*,*) 'Enter system size n:'
  	read(*,*) n

  	allocate(A(n, n), B(n))

  	!инициализация генератора случайных чисел
  	call random_seed()

  	!заполнение матрицы A случайными числами от -10 до 10
  	do i = 1, n
    		do j = 1, n
      			call random_number(r)
      			A(i, j) = -10.0_dp + 20.0_dp * r
    		end do
  	end do

  	!заполнение вектора B случайными числами от -10 до 10
  	do i = 1, n
    		call random_number(r)
    		B(i) = -10.0_dp + 20.0_dp * r
  	end do

  	!запись в data.dat
  	open(newunit=unit_out, file='data.dat', action='write', status='replace')
  	write(unit_out, '(A, 1X, i0)') '#', n
  	do i = 1, n
    		write(unit_out, '(*(g0, 1X))') A(i, :)
  	end do
  	do i = 1, n
    		write(unit_out, '(g0)') B(i)
  	end do
  	close(unit_out)

  	write(*,*) 'Random matrix generated. Size n =', n
  	deallocate(A, B)
end program gen_random
