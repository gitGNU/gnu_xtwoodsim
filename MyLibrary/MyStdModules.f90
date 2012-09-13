! *********************************************************************************
! *  XTwoodSim, A Solver for the Finite Element Software Elmer
! *  
! *  Copyright 11th September 2012 - , Uwe Jaschke
! *  
! *  This file is part of XTwoodSim.
! *  
! *      XTwoodSim is free software: you can redistribute it and/or modify
! *      it under the terms of the GNU General Public License as published by
! *      the Free Software Foundation, either version 3 of the License, or
! *      (at your option) any later version.
! *  
! *      XTwoodSim is distributed in the hope that it will be useful,
! *      but WITHOUT ANY WARRANTY; without even the implied warranty of
! *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *      GNU General Public License for more details.
! *  
! *      You should have received a copy of the GNU General Public License
! *      along with XTwoodSim.  If not, see <http://www.gnu.org/licenses/>.
! *  
! *      Diese Datei ist Teil von XTwoodSim.
! *  
! *      XTwoodSim ist Freie Software: Sie können es unter den Bedingungen
! *      der GNU General Public License, wie von der Free Software Foundation,
! *      Version 3 der Lizenz oder (nach Ihrer Option) jeder späteren
! *      veröffentlichten Version, weiterverbreiten und/oder modifizieren.
! *  
! *      XTwoodSim wird in der Hoffnung, dass es nützlich sein wird, aber
! *      OHNE JEDE GEWÄHELEISTUNG, bereitgestellt; sogar ohne die implizite
! *      Gewährleistung der MARKTFÄHIGKEIT oder EIGNUNG FÜR EINEN BESTIMMTEN ZWECK.
! *      Siehe die GNU General Public License für weitere Details.
! *  
! *      Sie sollten eine Kopie der GNU General Public License zusammen mit diesem
! *      Programm erhalten haben. Wenn nicht, siehe <http://www.gnu.org/licenses/>.
! *********************************************************************************
! My Standard Routines are collected 
!Tue Nov 29 12:28:01 CET 2011

MODULE MyStdModules
  USE DefUtils, ONLY:&
                     & INFO, &
                     & GetSolverParams, &
                     & GetInteger, &
                     & GetConstReal, &
                     & FATAL
  USE Types, Only: MESSAGE, & !For MESSAGE Type to write
                  & dp, &
                  & ValueList_t, &
                  & Solver_t
  IMPLICIT NONE

 Type StdSolverParams_t
  ! collects standard values used in every simulation
  ! in a single structure
   REAL(kind=dp) :: NonlinearTol
   INTEGER:: NonlinearIter
   INTEGER:: STDOFs
   Type(ValueList_t), POINTER :: SolverParams
  END Type StdSolverParams_t

  Type LocalSystemMatrices_t
  ! collects the tree system matrices Mass, Stiff and Force
   Real(kind=dp), Allocatable :: Mass(:,:), Stiff(:,:), Force(:)
  END Type LocalSystemMatrices_t

  !Overload subroutine print Array to print Vectors and Arrays
  INTERFACE printArray
   MODULE PROCEDURE print1dArray, print2DArray
  END INTERFACE


CONTAINS
 Function ReturnStdSolverParams(Solver) RESULT(StdSolverParams)
  Type(Solver_t) :: Solver
  Type(StdSolverParams_t):: StdSolverParams
  Type(ValueList_t), POINTER:: SolverParams
  LOGICAL :: FOUND 

  StdSolverParams % STDOFs = Solver % Variable % DOFs
  !---------------------------------------------------------
  !Read in solver parameters
  !---------------------------------------------------
  StdSolverParams % SolverParams => GetSolverParams(Solver)
  If ( .NOT. ASSOCIATED(StdSolverParams % SolverParams))&
   CALL FATAL('MoistureSolver','No Solver section FOUND')
  StdSolverParams % NonlinearIter = GetInteger(StdSolverParams % SolverParams, &
   'Nonlinear System Max Iterations', FOUND)
  IF (.NOT.FOUND) StdSolverParams % NonlinearIter = 1
  StdSolverParams % NonlinearTol = GetConstReal(StdSolverParams %  SolverParams, &
   'Nonlinear System Convergence Tolerance', FOUND)
  IF ( .NOT.FOUND) StdSolverParams % NonlinearTol = 1.0D-03
 END FUNCTION ReturnStdSolverParams

 SUBROUTINE printStdInfo(str,iter,NonlinearIter)
  INTEGER :: iter, NonlinearIter
  CHARACTER(len=*) :: str

        CALL Info( str, ' ', Level=4 ) 
        CALL Info( str, ' ', Level=4 ) 
        CALL Info( str, '-------------------------------------',Level=4 )
        WRITE( MESSAGE,* ) str, iter 
        CALL Info( str, MESSAGE, Level=4 ) 
        CALL Info( str, '-------------------------------------',Level=4 ) 
        CALL Info( str, ' ', Level=4 ) 
        CALL Info( str, 'Starting Assembly...', Level=4 ) 
        
        WRITE(MESSAGE,'(A,I5,A,I5)') 'Nonlinear iteration no.',&
              iter, ' of max. ', NonlinearIter
        CALL INFO(str, MESSAGE, level=1)
 END SUBROUTINE printStdInfo

 SUBROUTINE print2DArray(array, str)
  Real(kind=dp) :: array(:,:)
  CHARACTER(LEN=*) :: str
  INTEGER :: w, p, i, j

  w = size(array(1,:))
  p = size(array(:,1))

  print *, 'Array size of ', str, ' is ', p, 'x', w

    do i = 1, p
        write (*,10) (array(i,j), j = 1, w)
    end do  

10  format (24e12.3)

 END SUBROUTINE print2DArray

 SUBROUTINE print1dArray(array, str)
  Real(kind=dp) :: array(:)
  CHARACTER(LEN=*) :: str
  INTEGER :: w, p, i

  w = size(array(:))

  print *, 'Array size of ', str, ' is ',  w, ' x 1'
  write (*,20) (array(i), i = 1, w)

20  format (24e12.3)

 END SUBROUTINE print1dArray
 ! Routines taken from Elmer, thx
  SUBROUTINE RotateElasticityMatrix(C,T,dim)
!------------------------------------------------------------------------------
    INTEGER :: dim
    REAL(KIND=dp) :: T(:,:), C(:,:)
!------------------------------------------------------------------------------
    CALL Info("StressSolver", "Insode RotateElasticityMatrix", level=5)
    SELECT CASE(dim)
    CASE(2)
      CALL RotateElasticityMatrix2D(C,T)
    CASE(3)
      CALL RotateElasticityMatrix3D(C,T)
    END SELECT
!------------------------------------------------------------------------------
  END SUBROUTINE RotateElasticityMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE RotateElasticityMatrix2D(C,T)
!------------------------------------------------------------------------------
    IMPLICIT NONE

    REAL(KIND=dp) :: T(:,:), C(:,:), CT(2,2,2,2)
    INTEGER :: i,j,p,q,r,s
    INTEGER :: I1(3) = (/ 1,2,1 /), I2(3) = (/ 1,2,2 /)

    !
    ! Convert C-matrix to 4 index elasticity tensor:
    ! ----------------------------------------------
    CT = 0.0d0
    DO i=1,2
      p = I1(i)
      q = I2(i)
      DO j=1,2
        r = I1(j)
        s = I2(j)
        CT(p,q,r,s) = C(i,j)
        CT(p,q,s,r) = C(i,j)
        CT(q,p,r,s) = C(i,j)
        CT(q,p,s,r) = C(i,j)
      END DO
    END DO

    !
    ! Rotate the tensor:
    ! ------------------
    CALL Rotate4IndexTensor( CT, T, 2 )

    !
    ! Convert back to matrix form:
    ! ----------------------------
    DO i=1,2
      p = I1(i)
      q = I2(i)
      DO j=1,2
        r = I1(j)
        s = I2(j)
        C(i,j) = CT(p,q,r,s)
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE RotateElasticityMatrix2D
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE RotateElasticityMatrix3D(C,T)
!------------------------------------------------------------------------------
    IMPLICIT NONE

    REAL(KIND=dp) :: T(:,:), C(:,:), CT(3,3,3,3)
    INTEGER :: i,j,p,q,r,s
    INTEGER :: I1(6) = (/ 1,2,3,1,2,1 /), I2(6) = (/ 1,2,3,2,3,3 /)
    !
    ! Convert C-matrix to 4 index elasticity tensor:
    ! ----------------------------------------------
    CALL Info("StressSolver", "Inside RotateElasticityMatrix3D", level=5)
    CT = 0.0d0
    DO i=1,6
      p = I1(i)
      q = I2(i)
      DO j=1,6
        r = I1(j)
        s = I2(j)
        CT(p,q,r,s) = C(i,j)
        CT(p,q,s,r) = C(i,j)
        CT(q,p,r,s) = C(i,j)
        CT(q,p,s,r) = C(i,j)
      END DO
    END DO

    !
    ! Rotate the tensor:
    ! ------------------
    CALL Rotate4IndexTensor( CT, T, 3 )

    !
    ! Convert back to matrix form:
    ! ----------------------------
    DO i=1,6
      p = I1(i)
      q = I2(i)
      DO j=1,6
        r = I1(j)
        s = I2(j)
        C(i,j) = CT(p,q,r,s)
      END DO
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE RotateElasticityMatrix3D
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE Rotate2IndexTensor( C, T, dim )
!------------------------------------------------------------------------------
     INTEGER :: dim
     REAL(KIND=dp) :: C(:,:),T(:,:)
!------------------------------------------------------------------------------
     INTEGER :: i,j
     REAL(KIND=dp) :: C1(dim,dim)
!------------------------------------------------------------------------------
     C1 = 0
     DO i=1,dim
       DO j=1,dim
         C1(:,i) = C1(:,i) + T(i,j)*C(:,j)
       END DO
     END DO

     C = 0
     DO i=1,dim
       DO j=1,dim
         C(i,:) = C(i,:) + T(i,j)*C1(j,:)
       END DO
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE Rotate2IndexTensor
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE Rotate4IndexTensor( C, T, dim )
!------------------------------------------------------------------------------
     INTEGER :: dim
     REAL(KIND=dp) :: C(:,:,:,:),T(:,:)
!------------------------------------------------------------------------------
     INTEGER :: i,j
     REAL(KIND=dp) :: C1(dim,dim,dim,dim)
!------------------------------------------------------------------------------
     CALL INFO("StressSolver", "inside Rotate4IndexTensor", level=5)
     C1 = 0
     DO i=1,dim
       DO j=1,dim
         C1(:,:,:,i) = C1(:,:,:,i) + T(i,j)*C(:,:,:,j)
       END DO
     END DO

     C = 0
     DO i=1,dim
       DO j=1,dim
         C(:,:,i,:) = C(:,:,i,:) + T(i,j)*C1(:,:,j,:)
       END DO
     END DO

     C1 = 0
     DO i=1,dim
       DO j=1,dim
         C1(:,i,:,:) = C1(:,i,:,:) + T(i,j)*C(:,j,:,:)
       END DO
     END DO

     C = 0
     DO i=1,dim
       DO j=1,dim
         C(i,:,:,:) = C(i,:,:,:) + T(i,j)*C1(j,:,:,:)
       END DO
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE Rotate4IndexTensor
!------------------------------------------------------------------------------
! Ende Routines from Elmer
END MODULE MyStdModules
