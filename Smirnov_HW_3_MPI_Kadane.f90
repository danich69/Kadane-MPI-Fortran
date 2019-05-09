module smirnov_HW_3_MPI_Kadane
    use :: mpi
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

		real(8), intent(in), dimension(:,:) :: a
		integer(4), intent(out) :: x1, y1, x2, y2
		
		real(8), dimension(:,:), allocatable :: b
		real(8), dimension(:), allocatable :: maximum_S
		real(8), dimension(:), allocatable :: Global_maximum_S
		real(8), dimension(:), allocatable :: p
		real(8) :: maximum, CurrentSum, t2, t1
		
		integer(4), dimension(:), allocatable :: maximum_L, maximum_R, maximum_B
		integer(4) :: left, right, i, j, m, n, k
        integer(4) :: mpiErr, mpiSize, mpiRank
        
        call mpi_init(mpiErr)
        
        call mpi_comm_size(MPI_COMM_WORLD, mpiSize, mpiErr)
        call mpi_comm_rank(MPI_COMM_WORLD, mpiRank, mpiErr)

		m = size(a, dim = 1)
		n = size(a, dim = 2)
        
!        write(*,*) 'Rank is', mpiRank
        b = a
        
        if (m < n .and. mpiRank == 0) then
        
            allocate(b(n,m))
            b = transpose(a)
            m = k
            m = n
            n = k
            
!            write(*,*) 'I am', mpiRank, 'transposing'
            
            call mpi_bcast(b, m*n, MPI_REAL8, 0, MPI_COMM_WORLD, mpiErr)
            call mpi_bcast(m, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpiErr)
            call mpi_bcast(n, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpiErr)
        
        endif
        
		allocate(maximum_S(0:m/mpiSize))
		allocate(maximum_L(0:m/mpiSize))
		allocate(maximum_R(0:m/mpiSize))
		allocate(maximum_B(0:m/mpiSize))
        
!        write(*,*) 'Allocated, size', m/mpiSize
        
        if (mpiRank == 0) then
        
            allocate(Global_maximum_S(0:mpiSize-1))
!            write(*,*) 'Allocated, with ',mpiRank, 'size is', (m/mpiSize + 1)*mpiSize
        
        endif
        
        maximum_S = -1e-38
				
		allocate(p(n))
		
        call cpu_time(t1)
        
		do i = mpiRank, m-1, mpiSize
            
!            if (mod(i, mpiSize) == mpiRank) then

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
                
!            endif
            
		enddo
        
        x1 = mpiSize*(maxloc(maximum_S, dim = 1) - 1) + mpiRank
        x1 = maxloc(maximum_S, dim = 1) - 1
!        write(*,*) maximum_S, mpiRank
!        write(*,*) maximum_B, mpiRank
!        write(*,*) maximum_L, mpiRank
!        write(*,*) maximum_R, mpiRank
        x2 = maximum_B(x1)
        y2 = maximum_R(x1)
        y1 = maximum_L(x1)
        x1 = mpiSize*x1 + mpiRank + 1
        maximum = maxval(maximum_S(:))
        
!        call cpu_time(t2)
        
!        write(*,*) 'I am', mpiRank, 'my time is', t2-t1
!        write(*,*)
!        write(*,*) x1!, x2, y1, y2
        
!        do i = 1, m/mpiSize + 1
!        call cpu_time(t1)    
        call mpi_gather(maximum, 1, MPI_REAL8, Global_maximum_S, &
                    &1, MPI_REAL8, 0, MPI_COMM_WORLD, mpiErr)
!        call cpu_time(t2)            
!        write(*,*) 'I am', mpiRank, 'gather time is', t2-t1
!                call mpi_gather(maximum_L(i - 1), 1, MPI_INTEGER4, Global_maximum_L((i-1)*mpiSize+1:i*mpiSize), &
!                    &1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpiErr)
!                call mpi_gather(maximum_R(i - 1), 1, MPI_INTEGER4, Global_maximum_R((i-1)*mpiSize+1:i*mpiSize), &
!                    &1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpiErr)
!                call mpi_gather(maximum_B(i - 1), 1, MPI_INTEGER4, Global_maximum_B((i-1)*mpiSize+1:i*mpiSize), &
!                    &1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpiErr)
            
!        enddo
         
!         endif  

!        if (mpiRank == 0) then
!        write(*,*) 'i am the 0th thread my max_s', maximum_S
!        endif
!        if (mpiRank == 1) then
!!        write(*,*) 'i am the 1st thread my max_s', maximum_S
!        endif
		deallocate(p) 
        
        if (mpiRank == 0) then
            
!            write(*,*) Global_maximum_S
            i = maxloc(Global_maximum_S(0:mpiSize), dim = 1) - 1
!            write(*,*) i, Global_maximum_S(i)
            
        endif
        
        call mpi_bcast(i, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpiErr)
        
        
        call mpi_bcast(x1, 1, MPI_INTEGER4, i, MPI_COMM_WORLD, mpiErr)
        call mpi_bcast(x2, 1, MPI_INTEGER4, i, MPI_COMM_WORLD, mpiErr)
        call mpi_bcast(y1, 1, MPI_INTEGER4, i, MPI_COMM_WORLD, mpiErr)
        call mpi_bcast(y2, 1, MPI_INTEGER4, i, MPI_COMM_WORLD, mpiErr)
!        write(*,*) 'here', mpiRank
!        write(*,*) maximum_S
		deallocate(maximum_S)
!        write(*,*) 'S', mpiRank
		deallocate(maximum_L)
!        write(*,*) 'L', mpiRank
		deallocate(maximum_R)
!        write(*,*) 'R', mpiRank
		deallocate(maximum_B)
!        write(*,*) 'B', mpiRank
        if (mpiRank == 0) then
            
            deallocate(Global_maximum_S)
            
        endif

        call mpi_finalize(mpiErr)

	end subroutine

end module
