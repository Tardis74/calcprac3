program main
	use precision_mod, only: dp
	use linear_solver, only: solve
	implicit none
	integer :: method, n, i, unit_in, unit_out, status
	real(dp), allocatable :: A(:,:), B(:), X(:), res(:)
	real(dp) :: norm_res
	character(len=256) :: line, arg

	!чтение метода из командной строки
	if (command_argument_count() < 1) then
		print *, "Usage: program method (0 - Gauss no pivot, 1 - Gauss full pivot, 2 - Jordan)"
		stop
	end if
	call get_command_argument(1, arg)
	read(arg, *, iostat=status) method
	if (status /= 0 .or. method < 0 .or. method > 2) then
		print *, "Invalid method. Must be 0, 1, or 2."
		stop
	end if

	!чтение файла данных
	open(newunit=unit_in, file='data.dat', status='old', action='read', iostat=status)
	if (status /= 0) then
		print *, "Error opening data.dat"
		stop
	end if

	!чтение первой строки
	read(unit_in, '(A)') line
	i = index(line, ' ')
	if (i == 0) then
		print *, "Invalid first line format"
		stop
	end if
	read(line(i+1:), *) n
	print *, "System size n =", n

	allocate(A(n, n), B(n))

	!чтение матрицы A
	do i = 1, n
		read(unit_in, *) A(i, :)
	end do

	!чтение вектора B
	do i = 1, n
		read(unit_in, *) B(i)
	end do

	close(unit_in)

	!решение системы
	X = solve(A, B, method)

	!вычисление модуля невязки
	allocate(res(n))
	res = matmul(A, X) - B
	norm_res = sqrt(sum(res**2))
	print *, "Residual norm =", norm_res

	!запись результата
	open(newunit=unit_out, file='result.dat', action='write', iostat=status)
	if (status /= 0) then
		print *, "Error opening result.dat for writing"
		stop
  	end if

	write(unit_out, '(A, 1X, i0)') 'µ', n
  	do i = 1, n
    		write(unit_out, *) X(i)
  	end do

  	close(unit_out)

	deallocate(A, B, X, res)
end program main
