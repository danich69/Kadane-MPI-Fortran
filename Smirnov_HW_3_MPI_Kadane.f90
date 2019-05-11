module smirnov_HW_3_MPI_Kadane
    use :: mpi
    implicit none
    contains

    subroutine Kande(a, x1, x2, summary)                                                                ! 1D Kadane algorythm
    
        real(8), intent(in), dimension(:) :: a
        integer(4), intent(out) :: x1, x2    
        real(8), intent(out) :: summary
        integer(4) :: i, leftIndex, n, u
        real(8) :: Max_End, possible_1, possible_2

        n = size(a)

        summary = a(1)
        x1 = 1
        x2 = 1
        Max_End = a(1)
        leftIndex = 1
        
        do i = 2, n
            possible_1 = a(i)
            possible_2 = Max_End + a(i)

            if (possible_1 > possible_2) then
                Max_End = possible_1
                leftIndex = i
            else
                Max_End = possible_2
            endif

            if (Max_End > summary) then
                summary = Max_End
                x1 = leftIndex
                x2 = i
            endif
        enddo
        
    end subroutine


    subroutine GetMaxCoordinates(a, x1, y1, x2, y2)                                                       ! a - input matrix; x1, y1 - lower coordinates of submatrix; x2, y2 - upper coordinates of submatrix

        real(8), intent(in), dimension(:,:) :: a
        integer(4), intent(out) :: x1, y1, x2, y2
        
        real(8), dimension(:,:), allocatable :: b                                                         ! auxiliary matrix
        real(8), dimension(:), allocatable :: maximum_S                                                   ! array of maximum sums started at i-th row 
        real(8), dimension(:), allocatable :: Global_maximum_S                                            ! array of maximum sums calculated by i-th thread
        real(8), dimension(:), allocatable :: p                                                           ! array of column inner sums
        real(8) :: maximum, CurrentSum
        
        integer(4), dimension(2) :: coords
        integer(4), dimension(:), allocatable :: maximum_L, maximum_R, maximum_B                          ! array of left, right, bottom coordinates of submatrixes with maximum sums started at i-th row
        integer(4) :: left, right                                                                         ! coordinates of maximum subarray of p 
        integer(4) :: i, j, m, n, k                                                                       ! indexes; m, n - matrix sizes
        integer(4) :: mpiErr, mpiSize, mpiRank                                                            ! mpiErr - error at MPI functions, mpiRank - this thread number, mpiSize - number of threads
        
        call mpi_comm_size(MPI_COMM_WORLD, mpiSize, mpiErr)                                               ! get number of threads
        call mpi_comm_rank(MPI_COMM_WORLD, mpiRank, mpiErr)                                               ! get this thread rank
        
        if(maxval(a) < 0) then
            coords = maxloc(a)
            x2 = coords(1)
            x1 = coords(1)
            y2 = coords(2)
            y1 = coords(2)
            
            return
        endif

        m = size(a, dim = 1)
        n = size(a, dim = 2)
        
        if (m < n) then                                                                                    ! transpose a if needed
        write(*,*) 'here'
            allocate(b(n,m))
            b = transpose(a)
            m = k
            m = n
            n = k   
        else
            allocate(b(m,n))
            b = a            
        endif
        
        allocate(maximum_S(0:m/mpiSize))                                                                  ! array allocations
        allocate(maximum_L(0:m/mpiSize))
        allocate(maximum_R(0:m/mpiSize))
        allocate(maximum_B(0:m/mpiSize))
                
        if (mpiRank == 0) then
            allocate(Global_maximum_S(0:mpiSize-1))
        endif
        
        maximum_S = -1
                
        allocate(p(n))
        
        do i = mpiRank, m-1, mpiSize                                                                      ! main body of algorythm
            p = 0
            do j = i + 1, m
                do k = 1,n
                    p(k) = p(k) + b(j, k)
                enddo
                    
                call Kande(p, left, right, CurrentSum)
    
                if (CurrentSum  >  maximum_S( (i)/mpiSize) .or. i + 1 == j) then
                    maximum_S( (i)/mpiSize ) = CurrentSum;
                    maximum_L( (i)/mpiSize ) = left
                    maximum_R( (i)/mpiSize ) = right
                    maximum_B( (i)/mpiSize ) = j
                endif
            enddo
        enddo
        
        x1 = maxloc(maximum_S, dim = 1) - 1                                                               ! find coords of maximum sum in this thread
        x2 = maximum_B(x1)
        y2 = maximum_R(x1)
        y1 = maximum_L(x1)
        x1 = mpiSize*x1 + mpiRank + 1
        maximum = maxval(maximum_S(:))                                                                    ! find maximum sum in this thread
        
        call mpi_gather(maximum, 1, MPI_REAL8, Global_maximum_S, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpiErr) ! create array of all maximums
        
        deallocate(p) 
        
        if (mpiRank == 0) then
            i = maxloc(Global_maximum_S(0:mpiSize), dim = 1) - 1                                          ! find thread with maximum sum
        endif
        
        if(mpiSize > 1) then                                                                              ! no need to broadcast if there is only 1 thread
            call mpi_bcast(i, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpiErr)                                 ! inform other threads about it                
            call mpi_bcast(x1, 1, MPI_INTEGER4, i, MPI_COMM_WORLD, mpiErr)                                ! send answer to all threads
            call mpi_bcast(x2, 1, MPI_INTEGER4, i, MPI_COMM_WORLD, mpiErr)
            call mpi_bcast(y1, 1, MPI_INTEGER4, i, MPI_COMM_WORLD, mpiErr)
            call mpi_bcast(y2, 1, MPI_INTEGER4, i, MPI_COMM_WORLD, mpiErr)
        endif  
            
        deallocate(maximum_S)                                                                             ! deallocation of all arrays  
        deallocate(maximum_L)
        deallocate(maximum_R)
        deallocate(maximum_B)
        
        if (mpiRank == 0) then
            deallocate(Global_maximum_S)
        endif

    end subroutine

end module
