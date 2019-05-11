!-------------------------------------------------------------------------- 
    !                        Generalized Eigenvalue Solvers 
    !-------------------------------------------------------------------------- 
module eigreal
use parameters
implicit none


contains 

    subroutine realgeneig(M, doVecs, & 
            A, B, E, V, ierror) 
     
 
        !Declarations---------------------------------------------------------- 
        INTEGER,                                INTENT(IN)    :: M 
        LOGICAL,                                INTENT(IN)    :: doVecs 
        REAL(kind = dp), DIMENSION(M,M),        INTENT(IN)  :: A 
        REAL(kind = dp), DIMENSION(M,M),        INTENT(IN)  :: B 
        REAL(kind = dp), DIMENSION(M),            INTENT(OUT) :: E 
        REAL(kind = dp), DIMENSION(M,M),        INTENT(OUT) :: V 
        INTEGER,                                INTENT(OUT)    :: ierror 
        
        REAL(kind = dp), DIMENSION(M,M)                        :: workMat 
        REAL(kind = dp), DIMENSION(M,M)                        :: workOverlap 
        INTEGER, DIMENSION(M)                                :: IPIV 
        INTEGER                                             :: info 
        INTEGER                                             :: LDA 
        INTEGER                                                :: LDB 
        CHARACTER(LEN = 1)                                    :: JOBVR 
 
        CHARACTER(LEN = string_length)                      :: string1 
        CHARACTER(LEN = string_length)                      :: string2 
        CHARACTER(LEN = string_length)                      :: string3 
 
        !Implementations------------------------------------------------------- 
        ierror = 0 
         
        !Determining if eigenvectors should be computed (argument inspection) 
        if(doVecs.eqv..TRUE.) then 
            JOBVR = 'V' 
        else 
            JOBVR = 'N' 
        end if 
         
        !Making copies of matrices in order not to ruin input system 
        workMat(:,:)     = A(:,:) 
        workOverlap(:,:) = B(:,:) 
         
        !Setting leading rank 
        LDA = M 
        LDB = M 
 
        !Performing Cholesky factorization on Overlap matrix by LAPACK routine 
        call dpotrf('L', M, workOverlap, LDB, info) 
         
        !overlap is now actually the cholesky factor of the previous overlap 
        !Checking that Cholesky factorization went well 
        if(info.ne.0) then 
            ierror = info 
           !! if(verbose) then 
             !!   string1 = 'lapack.f95' 
               !! string2 =  'lapack_diagonalizeGeneralizedDoubleSymetric()' 
                !!WRITE(string3,'("dpotrf is unhappy and returned info = ", 1i9)')& 
                  !!  info 
                !!call print_PrintWarning(string1, string2, & 
                  !!  string3) 
           !! end if 
        end if 
         
        !Performing reduction to ordinary eigenvalue problem by LAPACK routine 
        call dsygst(1, 'L', M, workMat, LDA, workOverlap,  LDB, info) 
     
        !matrix is now the reduced problem L A L^* 
 
        !Checking that reduction to orinary problem went well 
        !!if(info.ne.0) then 
          !!  ierror = info 
            !!if(verbose) then 
              !!  string1 = 'lapack.f95' 
                !!string2 = 'lapack_diagonalizeGeneralizedDoubleSymetric()' 
                !!WRITE(string3,'("dsygst is unhappy and returned info = ",1i3)')& 
                  !!  info 
                !!call print_PrintWarning(string1, string2, & 
                  !!  string3) 
         !!   end if 
        !! end if  
 
        !Diagonalizing the ordinary eigenvalue problem by LAPACK wrapper 
        call lapack_diagSymRea(M, doVecs, workMat, E, V, ierror) 
        !!if(ierror.ne.0 ) then 
          !!  if(VERBOSE) then 
            !!    string1 = 'lapack.f95' 
              !!  string2 = 'lapack_diagGenSymRea()' 
                !!WRITE(string3,& 
                  !!  '("diagSymRea is unhappy and returned info = ", i9)')& 
                    !!    info 
               !! call print_PrintWarning(string1, string2, string3) 
            !!end if 
       !! end if 
         
        !Now matrix contains the eigenvalues of the reduced problem 
        !One must must solve the eigenvectors with the inverse  
        !cholesky factor to obtain eigenvectors to original problem 
             
        !This is only done if eigenvectors are to be computed 
        if(doVecs.eqv..TRUE.) then 
            !Clearing the Upper parts of choleksy factor 
            call lapack_ClearDiagsRea(M, workOverlap, 'U', ierror) 
         !!    if(ierror.ne..0) then 
           !!     if(verbose) then 
             !!       string1 = 'lapack.f95' 
               !!     string2 = 'lapack_diagGenSymRea()' 
                 !!   WRITE(string3,& 
                   !!     '("ClearDiagsRea is unhappy and returned ierror = ", i9)')& 
                     !!   ierror 
                    !!call print_PrintWarning(string1, string2, string3) 
                !!end if 
            !!end if            
 
            !Solving overlap system TODO: Try to make this more efficient, for 
            !instance by back substiturtion 
            call dgesv(M, M, TRANSPOSE(workOverlap),& 
                LDA, IPIV, V, LDB, info) 
 
            !!if(info.ne..0) then 
              !!  if(verbose) then 
                !!    string1 = 'lapack.f95' 
                  !!  string2 = 'lapack_diagonalizeGeneralizedSymtericReal()' 
                 !! 1  WRITE(string3,& 
                    !!    '("dgesv is unhappy and returned info = ", i9)')& 
                      !!  info 
                    !!call print_PrintWarning(string1, string2, string3) 
                !!end if 
            !!end if 
        !!end if 
         
        !!Deallocate work matrices 
         end if
    end subroutine 
 
    subroutine lapack_diagSymRea(M, doVecs, A, E, V, ierror) 
        !Includes-------------------------------------------------------------- 
        USE parameters 
        implicit none 
        !Declarations---------------------------------------------------------- 
        INTEGER,                                INTENT(IN)    :: M 
        LOGICAL,                                INTENT(IN)    :: doVecs 
        REAL(kind = dp),    DIMENSION(M,M),        INTENT(IN)  :: A 
        REAL(kind = dp),    DIMENSION(M),        INTENT(OUT) :: E 
        REAL(kind = dp),    DIMENSION(M,M),        INTENT(OUT) :: V 
        INTEGER,                                INTENT(OUT)    :: ierror 
 
        CHARACTER(LEN = 1)                                    :: jobVec 
 
        REAL(kind = dp),    DIMENSION(M,M)                    :: workMat 
        REAL(kind=dp),        DIMENSION(:), ALLOCATABLE       :: WORK 
        INTEGER                                             :: info 
        INTEGER                                             :: LWORK, LDA 
         
        CHARACTER(LEN = string_length)                      :: string1 
        CHARACTER(LEN = string_length)                      :: string2 
        CHARACTER(LEN = string_length)                      :: string3 
 
        !Implementations------------------------------------------------------- 
        ierror = 0 
 
        !Only compute eigenvectors if present in function call 
        if(doVecs.eqv..TRUE.) then 
            jobVec = 'V' 
        else 
            jobVec = 'E' 
        end if 
         
        !Allocating eigenvector matrix if eigenvectors are to be computed 
        if(jobVec.eq.'V') then 
            V(:,:) = A(:,:) 
        !Elsewise, allocate a work matrix to be discarded  
        else 
            V(:,:) = A(:,:) 
            !workMat(:,:) = A(:,:) 
        end if 
 
        !Setting leading rank 
        LDA = M 
 
        !Calculating efficient worksize by passing -1 to LAPACK algorithm 
        LWORK =  -1 
        ALLOCATE(WORK(1)) 
         
        !Calling LAPACK diagonalization routine 
        !If eigenvectors are to be computed, perform proper call  
        if(jobVec.eq.'V') then 
            call dsyev(jobVec, 'L', M, V, LDA,& 
                E, WORK, LWORK, info) 
            !Reallocate workspace to optimal size 
            LWORK = WORK(1) 
            DEALLOCATE(WORK) 
            ALLOCATE(WORK(LWORK)) 
            call dsyev(jobVec, 'L', M, V, LDA, E,& 
                WORK, LWORK, info) 
        !Elsewise, use workMat 
        else 
            call dsyev(jobVec, 'L', M, V, LDA,& 
                E, WORK, LWORK, info) 
            !Reallocate workspace to optimal size 
            LWORK = WORK(1) 
            DEALLOCATE(WORK) 
            ALLOCATE(WORK(LWORK)) 
            call dsyev(jobVec, 'L', M, V, LDA, E,& 
                WORK, LWORK, info) 
        end if 
 
        !Performing check that everything went well 
        if(info.ne.0) then 
            ierror = info 
            if(verbose) then 
                string1 = 'lapack.f95' 
                string2 = 'lapack_diagonalizeDoubleSymetric()' 
                WRITE(string3, '("dsyev returned info = ",1i4)') info 
                call print_PrintWarning(string1, string2, & 
                    string3) 
            end if 
        end if 
 
        !Performing check that workspace was satisfactory 
        if(WORK(1).ne.LWORK) then 
            if(verbose) then 
                string1 = 'lapack.f95' 
                string2 = 'lapack_diagonalizeDoubleSymetric()' 
                WRITE(string3, & 
        '("Unoptimal workspace. Lapack used: ",i3,", optimal is: ",f8.1)') & 
                LWORK, WORK(1) 
                call print_PrintWarning(string1, string2, & 
                    string3) 
            end if 
        end if 
 
        !Deallocating work arrays 
        DEALLOCATE(WORK) 
     
    end subroutine

    subroutine print_Print(string)
  use parameters
 implicit none    
    CHARACTER (LEN = string_length)                            :: string
   
  write(*,*) string
   
    end subroutine
   
  subroutine print_PrintWarning(string1, string2, string3)
  use parameters
  implicit none  
  CHARACTER(LEN = string_length)                      :: string1 
        CHARACTER(LEN = string_length)                      :: string2 
        CHARACTER(LEN = string_length)                      :: string3 
   write(*,*) string1, '/n',string2,'/n',string3,'/n'
  end subroutine

subroutine lapack_ClearDiagsRea(M, A, which, ierror)
!Includes--------------------------------------------------------------
        use parameters
        implicit none
!Declarations----------------------------------------------------------
        INTEGER,                         INTENT(IN)        :: M
        REAL(kind = dp), DIMENSION(M,M), INTENT(INOUT)  :: A
        CHARACTER(LEN = 1),              INTENT(IN)     :: which
        INTEGER,                         INTENT(OUT)    :: ierror
        INTEGER                                         :: idx

        CHARACTER(LEN = string_length)                  :: string1
        CHARACTER(LEN = string_length)                  :: string2
        CHARACTER(LEN = string_length)                  :: string3
!Implementations-------------------------------------------------------
        ierror = 0
        !Looping over columns
        do idx = 1, M
            !If which is L, clear lower part
            if(which.eq.'L') then
                A(idx+1:,idx) = 0.d0
            !elsewise if which is U, clear upper part
            else if(which.eq.'U') then
                A(:idx-1,idx) = 0.d0
            !elsewise if which is D, clear diagonal
            else if(which.eq.'D') then
                A(idx, idx) = 0.d0
            !Else, provide warning that nothing is done
            else
                ierror = -2
                !!if(verbose) then
                  !!  string1 = 'matrixtools.f95'
                    !!string2 = 'matrixtools_ClearDiagonalsReal()'
                    !!string3 = 'Cannot determine what to remove (must be L,U,D)'
                    !!call print_PrintWarning(string1, &
                      !!  string2, string3)
                !!end if
            end if
        end do
    end subroutine

    subroutine lapack_ClearDiagsCplx(M, A, which, ierror)
!Includes--------------------------------------------------------------
        use parameters
        implicit none
!Declarations----------------------------------------------------------
        INTEGER,                                INTENT(IN)      :: M
        COMPLEX(kind = dp), DIMENSION(M,M),     INTENT(INOUT) :: A
        CHARACTER(LEN = 1),                     INTENT(IN)    :: which
        INTEGER,                                INTENT(OUT)      :: ierror
        INTEGER                                               :: N
        INTEGER                                               :: idx
        CHARACTER(LEN = string_length)                        :: string1
        CHARACTER(LEN = string_length)                        :: string2
        CHARACTER(LEN = string_length)                        :: string3
!Implementations-------------------------------------------------------
        ierror = 0

        !Looping over columns
        do idx = 1, M
            !If which is L, clear lower part (sub diags)
            if(which.eq.'L') then
                A(idx+1:,idx) = 0.d0
            !elsewise if which i U, clea upper part (super diags)
            else if(which.eq.'U') then
                A(:idx-1,idx) = 0.d0
            !Elsewise if which i D clear the main diagonal
            else if(which.eq.'D') then
                A(idx, idx) = 0.d0
            !Else, provide warning that nothing is done
            else
                !Something is wrong with wich
                ierror = -2
                !!if(verbose) then
                 !!   string1 = 'matrixtools.f95'
                   !! string2 = 'matrixtools_ClearDiagonalsComplex()'
                    !!string3 = 'Cannot determine what to remove (must be L,U,D)'
                    !!call print_PrintWarning(string1, string2, &
                      !!  string3)
                !!end if
            end if
        end do
    end subroutine

end module
