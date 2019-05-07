module Smirnov_Hw_1_OpenMP
    implicit none
    contains

	subroutine Kande(a, x1, x2, summary)
	
		real(8), intent(in), dimension(:) :: a
		integer(4), intent(out) :: x1, x2	
		real(8), intent(out) :: summary
		integer(4) :: i, leftIndex, n, u
		real(8) :: Max_End, possible_1, possible_2

		n = size(a)

		summary = a(1)
		x1 = 1
		x2 = 1

		Max_End=a(1)
		leftIndex=1
		
		do i=2,n

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


	subroutine GetMaxCoordinates(a, x1, y1, x2, y2)

		real(8), intent(in), dimension(:,:) :: A
		integer(4), intent(out) :: x1, y1, x2, y2
		
		real(8), dimension(:), allocatable :: maximum_S
		real(8), dimension(:), allocatable :: Global_maximum_S
		real(8), dimension(:), allocatable :: p
		real(8) :: maximum, CurrentSum
		
		integer(4), dimension(:), allocatable :: Global_maximum_L, Global_maximum_R, Global_maximum_B
		integer(4), dimension(:), allocatable :: maximum_L, maximum_R, maximum_B
		integer(4) :: left, right, i, j, m, n, k
        integer(4) :: mpiErr, mpiSize, mpiRank
        
        call mpi_comm_size(MPI_COMM_WORLD, mpiSize, mpiErr)
        call mpi_comm_rank(MPI_COMM_WORLD, mpiRank, mpiErr)

		m = size(a, dim = 1)
		n = size(a, dim = 2)
        
        if (m < n .and. mpiRank == 0) then
        
            A = transpose(A)
            m = k
            m = n
            n = k
            
            call mpi_bcast(A, m*n, MPI_REAL8, 0, MPI_COMM_WORLD, mpiErr)
            call mpi_bcast(m, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpiErr)
            call mpi_bcast(n, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpiErr)
        
        endif
		
		allocate(maximum_S(0:m/mpiSize))
		allocate(maximum_L(0:m/mpiSize))
		allocate(maximum_R(0:m/mpiSize))
		allocate(maximum_B(0:m/mpiSize))
        
        if (mpiRank == 0) then
        
            allocate(Global_maximum_S(m/mpiSize + 1))
            allocate(Global_maximum_L(m/mpiSize + 1))
            allocate(Global_maximum_R(m/mpiSize + 1))
            allocate(Global_maximum_B(m/mpiSize + 1))
        
        endif
        
        maximum_S = -1e-38
				
		allocate(p(n))

		
		do i = 1, m
            
            if (mod(i, mpiSize) == mpiSize) then

                p = 0
            
                do j = i, m
            
                    do k = 1,n
                        p(k) = p(k) + a(j, k)
                    enddo
				
                    call Kande(p, left, right, CurrentSum)
    
                    if (CurrentSum  >  maximum_S( (i - 1)/mpiSize) .or. i == j) then
                
                        maximum_S( (i - 1)/mpiSize ) = CurrentSum;
                        maximum_L( (i - 1)/mpiSize ) = left
                        maximum_R( (i - 1)/mpiSize ) = right
                        maximum_B( (i - 1)/mpiSize ) = j
                    
                    endif
                
                enddo
                
            endif
            
		enddo
        
        if (mpiRank == 0) then
            
            do i = 1, m/mpiSize + 1
            
                call mpi_gather(maximum_S(i - 1), 1, MPI_REAL8, Global_maximum_S((i-1)*mpiSize+1:i*mpiSize), mpiSize, MPI_REAL8, 0, MPI_COMM_WORLD, mpiErr)
                call mpi_gather(maximum_L(i - 1), 1, MPI_INTEGER4, Global_maximum_L((i-1)*mpiSize+1:i*mpiSize), mpiSize, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpiErr)
                call mpi_gather(maximum_R(i - 1), 1, MPI_INTEGER4, Global_maximum_L((i-1)*mpiSize+1:i*mpiSize), mpiSize, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpiErr)
                call mpi_gather(maximum_B(i - 1), 1, MPI_INTEGER4, Global_maximum_L((i-1)*mpiSize+1:i*mpiSize), mpiSize, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpiErr)
            
            enddo
            
        call 		

		deallocate(p)
        
        if (mpiRank == 0) then

            x1 = maxloc(Global_maximum_S, dim = 1)
            y1 = Global_maximum_L(x1)
            x2 = Global_maximum_B(x1)
            y2 = Global_maximum_R(x1)
            
            deallocate(Global_maximum_S)
            deallocate(Global_maximum_L)
            deallocate(Global_maximum_R)
            deallocate(Global_maximum_B)
            
        enddo
		
		deallocate(maximum_S)
		deallocate(maximum_L)
		deallocate(maximum_R)
		deallocate(maximum_B)


	end subroutine

end module
