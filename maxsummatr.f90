module Homework
    contains
        subroutine FindMaxCoordinates(A, x1, y1, x2, y2)
        implicit none
        include "mpif.h"
        integer(4) :: mpiErr, mpiSize, mpiRank, task_step, task_finish, task_start, i, index, err
        real(8), intent(in), dimension(:,:) :: A
        integer(4), intent(out) :: x1, y1, x2, y2
        integer(4) :: n, L, R, Up, Down, m, tmp, border
        real(8), allocatable :: current_column(:), B(:,:), max_sum(:), yetanothersum(:)
        real(8) :: current_sum, optimal
        integer(4), allocatable, dimension(:) :: X_1, Y_1, X_2, Y_2
        logical :: transpos
        real(8) :: start, finish
        call mpi_comm_size(MPI_COMM_WORLD, mpiSize, mpiErr)
        call mpi_comm_rank(MPI_COMM_WORLD, mpiRank, mpiErr)

        m = size(A, dim=1) 
        n = size(A, dim=2) 
        transpos = .FALSE.
        if (m < n) then 
            transpos = .TRUE.   
            B = transpose(A)
            m = size(B, dim=1) 
            n = size(B, dim=2) 
        else
            B = A     
            endif

        allocate(current_column(m))
        allocate(max_sum(mpiSize))
        allocate(yetanothersum(mpiSize))
        allocate(X_1(mpiSize))
        allocate(X_2(mpiSize))
        allocate(Y_1(mpiSize))
        allocate(Y_2(mpiSize))
        max_sum = -huge(0)

        X_1=1
        Y_1=1
        X_2=1
        Y_2=1

!     call cpu_time(start)

        if (n < mpiSize) then
            task_start = 1
            task_finish = n

        else
                  task_step = n/mpiSize**2
                  task_start = 1 + task_step*mpiRank**2 
                  task_finish = task_step*(mpiRank+1)**2
                  if (mpiRank == mpiSize - 1) then
                    task_finish = n
                  endif
        endif
        L = task_start
        do while (L <= task_finish)
            current_column = B(:, L)            
            do R=L,n
                if (R > L) then 
                    current_column = current_column + B(:, R)
                endif
                call FindMaxInArray(current_column, current_sum, Up, Down) 
                if (current_sum > max_sum(mpiRank+1)) then
                    max_sum(mpiRank+1) = current_sum
                    X_1(mpiRank+1) = Up
                    X_2(mpiRank+1) = Down
                    Y_1(mpiRank+1) = L
                    Y_2(mpiRank+1) = R
                endif
            end do
            L = L + 1
        end do

!         call cpu_time(finish)
!         print '("Time = ",f6.3, " seconds.")', (finish-start)

        call mpi_reduce(max_sum, yetanothersum, mpiSize, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, err)   
        call mpi_bcast(yetanothersum, mpiSize, MPI_REAL8, 0, MPI_COMM_WORLD, err)

        index = maxloc(yetanothersum, dim=1) - 1  
        call mpi_bcast(X_1, mpiSize, MPI_INTEGER4, index, MPI_COMM_WORLD, err)
        call mpi_bcast(X_2, mpiSize, MPI_INTEGER4, index, MPI_COMM_WORLD, err)
        call mpi_bcast(Y_1, mpiSize, MPI_INTEGER4, index, MPI_COMM_WORLD, err)
        call mpi_bcast(Y_2, mpiSize, MPI_INTEGER4, index, MPI_COMM_WORLD, err)

        x1 = X_1(maxloc(yetanothersum, dim=1))
        x2 = X_2(maxloc(yetanothersum, dim=1))
        y1 = Y_1(maxloc(yetanothersum, dim=1))
        y2 = Y_2(maxloc(yetanothersum, dim=1))
        if (transpos) then  
            tmp = x1
            x1 = y1
            y1 = tmp
            tmp = y2
            y2 = x2
            x2 = tmp
            endif


        deallocate(B)
        deallocate(current_column)
        deallocate(yetanothersum)
        deallocate(max_sum)        
        deallocate(X_1)
        deallocate(X_2)
        deallocate(Y_1)
        deallocate(Y_2)
    end subroutine

        subroutine FindMaxInArray(a, Sum, Up, Down)

            real(8), intent(in), dimension(:) :: a
            integer(4), intent(out) :: Up, Down
            real(8), intent(out) :: Sum
            real(8) :: cur_sum
            integer(4) :: minus_pos, i

            Sum = a(1)
            Up = 1
            Down = 1
            cur_sum = 0
            minus_pos = 0

            do i=1, size(a)
                cur_sum = cur_sum + a(i)
            if (cur_sum > Sum) then
                Sum = cur_sum
                Up = minus_pos + 1
                Down = i
                endif
            if (cur_sum < 0) then
                cur_sum = 0
                minus_pos = i
                endif
            enddo
        end subroutine FindMaxInArray

end module Homework