! --------------------------------------------------------------------
! 
!    This program can sort a set of numbers.  The method used is 
! usually referred to as "selection" method.
! --------------------------------------------------------------------

module sorting_module

  USE parameters

  
  !INTEGER, PARAMETER           :: dp=8

contains 

   subroutine  Sorting(ActualSize, InputData)
   IMPLICIT  NONE
     
   INTEGER, INTENT(in)                        :: ActualSize
   REAL(kind=dp), DIMENSION(ActualSize), INTENT(inout) :: InputData
   INTEGER                        :: i
  
   CALL  sort(ActualSize, InputData)
   
   open(7,file='StateEnergies.dat')
   !WRITE(*,7)
  ! WRITE(*,7) "Sorted Array:"
  do i=1,ActualSize
   WRITE(7,*) InputData(1:ActualSize),'index:',ActualSize,'/n'
 end do
   close(7)
   end subroutine Sorting

! --------------------------------------------------------------------
! INTEGER FUNCTION  FindMinimum():
!    This function returns the location of the minimum in the section
! between Start and End.
! --------------------------------------------------------------------

INTEGER FUNCTION  FindMinimum(ss, x, StaIdx)
        IMPLICIT  NONE
        Integer, INTENT(in) :: ss
        REAL(kind=dp), DIMENSION(ss), INTENT(IN) :: x
         INTEGER, INTENT(IN)                :: StaIdx
        REAL(kind=dp)                            :: Minimum
!       INTEGER                              :: FindMinimum
       INTEGER                            :: Location
       INTEGER                            :: i
      
      Minimum  = x(StaIdx)		! assume the first is the min
      Location = StaIdx			! record its position
      DO i = StaIdx+1, ss		! start with next elements
         IF (x(i) < Minimum) THEN	!   if x(i) less than the min?
            Minimum  = x(i)		!      Yes, a new minimum found
            Location = i                !      record its position
         END IF
      END DO
      FindMinimum = Location        	! return the position
   END FUNCTION  FindMinimum

! --------------------------------------------------------------------
! SUBROUTINE  Swap():
!    This subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------

   SUBROUTINE  Swap(a, b)
      IMPLICIT  NONE
      REAL(kind=dp), INTENT(INOUT) :: a, b
      REAL(kind=dp)                :: Temp

      Temp = a
      a    = b
      b    = Temp
   END SUBROUTINE  Swap

! --------------------------------------------------------------------
! SUBROUTINE  Sort():
!    This subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------

   SUBROUTINE  Sort(ss, x)
      IMPLICIT  NONE
      INTEGER, INTENT(IN)                   :: ss
      REAL(kind=dp), DIMENSION(ss), INTENT(INOUT) :: x
      INTEGER                               :: i
      INTEGER                               :: Location

      DO i = 1, ss-1			! except for the last
         Location = FindMinimum(ss,x, i)	! find min from this to last
         CALL  Swap(x(i), x(Location))	! swap this and the minimum
      END DO
   END SUBROUTINE  Sort


  end module

