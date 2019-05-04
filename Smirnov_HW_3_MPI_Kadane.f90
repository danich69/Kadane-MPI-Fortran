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
		real(8), dimension(:), allocatable :: p
		real(8) :: maximum, CurrentSum
		
		integer(4), dimension(:), allocatable :: maximum_L, maximum_R, maximum_B
		integer(4) :: left, right, top, bottom, m, n, j

		m = size(a, dim = 1)
		n = size(a, dim = 2)
		
		allocate(maximum_S(m))
		allocate(maximum_L(m))
		allocate(maximum_R(m))
		allocate(maximum_B(m))
				
		allocate(p(n))

		
		do top = 1, m
        
			p = 0
            
			do bottom = top, m
            
				do j = 1,n
					p(j) = p(j) + a(bottom, j)
				enddo
				
				call Kande(p, left, right, CurrentSum)
    
				if (CurrentSum  >  maximum_S(top) .or. top == bottom) then
					maximum_S(top) = CurrentSum;
					maximum_L(top) = left
					maximum_R(top) = right
					maximum_B(top) = bottom
				endif
                
			enddo
            
		enddo
		

		deallocate(p)


		maximum = maximum_S(1)
		x1 = 1
		y1 = maximum_L(1)
		x2 = maximum_B(1)
		y2 = maximum_R(1)
		
		do j = 2, m
        
			if (maximum_S(j)  >  maximum) then
		
        		maximum = maximum_S(j)
				x1 = j
				y1 = maximum_L(j)
				x2 = maximum_B(j)
				y2 = maximum_R(j)
		
            endif
		
        enddo

		deallocate(maximum_S)
		deallocate(maximum_L)
		deallocate(maximum_R)
		deallocate(maximum_B)

	end subroutine

end module
