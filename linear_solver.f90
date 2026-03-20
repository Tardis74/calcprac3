module linear_solver
  	use precision_mod, only: dp
  	implicit none
contains
  	! method = 0 – Гаусс без выбора главного элемента
  	! method = 1 – Гаусс с полным выбором
  	! method = 2 – Жордан
  	function solve(A, B, method) result(X)
    		real(dp), intent(in) :: A(:,:), B(:)
    		integer, intent(in) :: method
    		real(dp), allocatable :: X(:)

    		integer :: n, i, j, k, maxloc_row, maxloc_col
    		real(dp), allocatable :: M(:,:)          !расширенная матрица
    		integer, allocatable :: col_perm(:)      !перестановки столбцов
    		real(dp) :: temp, maxval, pivot
    		real(dp), allocatable :: temp_arr(:)

    		n = size(A, 1)
    		allocate(M(n, n+1), col_perm(n))
    		M(:, 1:n) = A
    		M(:, n+1) = B
    		col_perm = [(i, i = 1, n)]   !исходный порядок столбцов

    		select case (method)
    			case (0)  !Гаусс без выбора
      				do k = 1, n
					pivot = M(k, k)
					if (abs(pivot) < 1e-12_dp) then
					  	print *, "Warning: near-zero pivot at step", k, " value =", pivot
					end if
					!нормализация строки k
					M(k, k:n+1) = M(k, k:n+1) / pivot
					!исключение ниже
					!$OMP PARALLEL DO PRIVATE(i, j, temp) SHARED(M, k, n)
					do i = k+1, n
						temp = M(i, k)
						do j = k, n+1
							M(i, j) = M(i, j) - temp * M(k, j)
						end do
					end do
					!$OMP END PARALLEL DO
				end do
					      !обратная подстановка
				allocate(X(n))
				do i = n, 1, -1
					X(i) = M(i, n+1)
					do j = i+1, n
						X(i) = X(i) - M(i, j) * X(j)
					end do
				end do

		  	case (1)  !Гаусс с полным выбором
		    		do k = 1, n
					!поиск максимального элемента в подматрице M(k:n, k:n)
					maxval = abs(M(k, k))
					maxloc_row = k; maxloc_col = k
					do i = k, n
				  		do j = k, n
				    			if (abs(M(i, j)) > maxval) then
				      				maxval = abs(M(i, j))
				      				maxloc_row = i
				      				maxloc_col = j
				    			end if
				  		end do
					end do
					if (maxval < 1e-12_dp) then
				  		print *, "Warning: near-zero pivot at step", k, " value =", M(maxloc_row, maxloc_col)
					end if
					!перестановка строк
					if (maxloc_row /= k) then
				  		do j = k, n+1
				    			temp = M(k, j)
				    			M(k, j) = M(maxloc_row, j)
				    			M(maxloc_row, j) = temp
				  		end do
					end if
					!перестановка столбцов
					if (maxloc_col /= k) then
				  		do i = 1, n
					    		temp = M(i, k)
					    		M(i, k) = M(i, maxloc_col)
					    		M(i, maxloc_col) = temp
				  		end do
				  		i = col_perm(k)
				  		col_perm(k) = col_perm(maxloc_col)
				  		col_perm(maxloc_col) = i
					end if
					!нормализация строки k
					M(k, k:n+1) = M(k, k:n+1) / M(k, k)
					!исключение ниже
					!$OMP PARALLEL DO PRIVATE(i, j, temp) SHARED(M, k, n)
						do i = k+1, n
				  			temp = M(i, k)
				  			do j = k, n+1
				    				M(i, j) = M(i, j) - temp * M(k, j)
				  			end do
						end do
					!$OMP END PARALLEL DO
		      		end do
		      		!обратная подстановка
		      		allocate(X(n))
		      		do i = n, 1, -1
					X(i) = M(i, n+1)
					do j = i+1, n
			  			X(i) = X(i) - M(i, j) * X(j)
					end do
		      		end do
		      		!восстановление исходного порядка переменных
		      		allocate(temp_arr(n))
		      		temp_arr = X
		      		do i = 1, n
					X(col_perm(i)) = temp_arr(i)
		      		end do
		      		deallocate(temp_arr)

    			case (2)  !Жордан с полным выбором
      				do k = 1, n
					!поиск максимального элемента в подматрице M(k:n, k:n)
					maxval = abs(M(k, k))
					maxloc_row = k; maxloc_col = k
					do i = k, n
		  				do j = k, n
		    					if (abs(M(i, j)) > maxval) then
		      						maxval = abs(M(i, j))
		      						maxloc_row = i
		      						maxloc_col = j
		    					end if
		  				end do
					end do
					if (maxval < 1e-12_dp) then
		  				print *, "Warning: near-zero pivot at step", k, " value =", M(maxloc_row, maxloc_col)
					end if
					!перестановка строк
					if (maxloc_row /= k) then
		  				do j = k, n+1
		    					temp = M(k, j)
		    					M(k, j) = M(maxloc_row, j)
		    					M(maxloc_row, j) = temp
		  				end do
					end if
        				!перестановка столбцов
        				if (maxloc_col /= k) then
          					do i = 1, n
							temp = M(i, k)
						    	M(i, k) = M(i, maxloc_col)
						    	M(i, maxloc_col) = temp
          					end do
						i = col_perm(k)
						col_perm(k) = col_perm(maxloc_col)
						col_perm(maxloc_col) = i
        				end if
        				!нормализация строки k
        				M(k, k:n+1) = M(k, k:n+1) / M(k, k)
        				!исключение во всех остальных строках
					!$OMP PARALLEL DO PRIVATE(i, j, temp) SHARED(M, k, n)
        				do i = 1, n
          					if (i == k) cycle
          					temp = M(i, k)
          					do j = k, n+1
            						M(i, j) = M(i, j) - temp * M(k, j)
          					end do
        				end do
					!$OMP END PARALLEL DO
      				end do
      				!решение – последний столбец
      				allocate(X(n))
      				X = M(:, n+1)
      				!восстановление исходного порядка переменных
      				allocate(temp_arr(n))
      				temp_arr = X
      				do i = 1, n
        				X(col_perm(i)) = temp_arr(i)
      				end do
      				deallocate(temp_arr)

    			case default
      				print *, "Error: unknown method", method
      				allocate(X(0))
      				return
    		end select

		deallocate(M, col_perm)
	end function solve

end module linear_solver
