
module Homework
    use omp_lib
    implicit none
    contains
        subroutine FindMaxCoordinates(A, x1, y1, x2, y2)
        implicit none
        real(8), intent(in), dimension(:,:) :: A
        integer(4), intent(out) :: x1, y1, x2, y2
        integer(4) :: n, L, R, Up, Down, m, tmp, amount_in_L, L_finish, L_start
        integer(4), allocatable, dimension(:) :: X_1, X_2, Y_1, Y_2
        real(8), allocatable, dimension(:) :: max_sum
        real(8), allocatable :: current_column(:), B(:,:)
        real(8) :: current_sum
        logical :: transposition

        m = size(A, dim=1) 
        n = size(A, dim=2) 
        transposition = .FALSE.



        if (m < n) then 
            transposition = .TRUE.   
            B = transpose(A)
            m = size(B, dim=1) 
            n = size(B, dim=2) 
        else
            B = A     
            endif

        allocate(current_column(m))


        !$omp parallel shared(amount_in_L, X_1, X_2, Y_1, Y_2, max_sum, B), private(current_column, current_sum, Up, Down)


        !$omp single 
        allocate(max_sum(omp_get_num_threads()))
        allocate(X_1(omp_get_num_threads()))
        allocate(X_2(omp_get_num_threads()))
        allocate(Y_1(omp_get_num_threads()))
        allocate(Y_2(omp_get_num_threads()))
        max_sum = B(1,1)
        X_1 = 1
        X_2 = 1
        Y_1 = 1
        Y_2 = 1
        !$omp end single

        !$omp do schedule(guided)
        do L=1, n

            current_column = B(:, L)            
            do R=L,n
 
                if (R > L) then 
                    current_column = current_column + B(:, R)
                endif
                
                call FindMaxInArray(current_column, current_sum, Up, Down) 


                      
                if (current_sum > max_sum(omp_get_thread_num()+1)) then
                    max_sum(omp_get_thread_num()+1) = current_sum
                    X_1(omp_get_thread_num()+1) = Up
                    X_2(omp_get_thread_num()+1) = Down
                    Y_1(omp_get_thread_num()+1) = L
                    Y_2(omp_get_thread_num()+1) = R
                endif
            end do
        end do

       !$omp end parallel 


        x1 = X_1(maxloc(max_sum, dim=1))
        x2 = X_2(maxloc(max_sum, dim=1))
        y1 = Y_1(maxloc(max_sum, dim=1))
        y2 = Y_2(maxloc(max_sum, dim=1))

        deallocate(max_sum)
        deallocate(X_1)
        deallocate(X_2)
        deallocate(Y_1)
        deallocate(Y_2)


        deallocate(current_column)


        if (transposition) then  
            tmp = x1
            x1 = y1
            y1 = tmp
    
            tmp = y2
            y2 = x2
            x2 = tmp
            endif

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



