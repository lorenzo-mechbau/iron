
!> \file
!> \author Benjamin Maier
!> \brief This module contains routines that print the composite values of each type defined in types.f90
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand, the University of Oxford, Oxford, United
!> Kingdom and King's College, London, United Kingdom. Portions created
!> by the University of Auckland, the University of Oxford and King's
!> College, London are Copyright (C) 2007-2010 by the University of
!> Auckland, the University of Oxford and King's College, London.
!> All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>
!> This module contains a print routine for each type in types.f90
!> The routines were created automatically at 2018-01-11 07:29:39.

MODULE PRINT_TYPES_ROUTINES_SELECTION
  USE Types
  
  IMPLICIT NONE
  
CONTAINS

FUNCTION GetPrintIndent(Depth)
  INTEGER(INTG), INTENT(IN) :: Depth
  CHARACTER(LEN=400) :: GetPrintIndent
  INTEGER(INTG) :: I
  
  WRITE(GetPrintIndent, "(I3)") Depth
  GetPrintIndent = TRIM(GetPrintIndent) // "|"
  DO I = 1, Depth
    GetPrintIndent = TRIM(GetPrintIndent) // " ."
  END DO

END FUNCTION GetPrintIndent


  !
  !================================================================================================================================
  !
  SUBROUTINE Print_DOMAIN_MAPPING(Variable, MaxDepth, MaxArrayLength)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER, INTENT(IN) :: Variable   !< the variable to be printed
    INTEGER(INTG), INTENT(IN) :: MaxDepth           !< the maximum recursion depth down to which data is printed
    INTEGER(INTG), INTENT(IN) :: MaxArrayLength   !< the maximum array length that is printed  


    PRINT*, "Print TYPE(DOMAIN_MAPPING_TYPE), POINTER :: Variable with maximum depth ",MaxDepth,", maximum array length:", &
      & MaxArrayLength
    CALL Print_DOMAIN_MAPPING_TYPE(Variable, 1, MaxDepth, MaxArrayLength)
  END SUBROUTINE Print_DOMAIN_MAPPING

!
!================================================================================================================================
!
RECURSIVE SUBROUTINE Print_DOMAIN_MAPPING_TYPE(Variable, Depth, MaxDepth, MaxArrayLength)
  TYPE(DOMAIN_MAPPING_TYPE), POINTER, INTENT(IN) :: Variable
  INTEGER(INTG), INTENT(IN) :: Depth   !< the recursion depth
  INTEGER(INTG), INTENT(IN) :: MaxDepth   !< the maximum recursion depth for which data is printed
  INTEGER(INTG), INTENT(IN) :: MaxArrayLength   !< the maximum array length that is printed  
  CHARACTER(LEN=400) :: PrintIndent   !< a string containing 2*Depth space characters
  LOGICAL :: IsAllocated
  LOGICAL :: IsAssociated
  INTEGER(INTG) :: Dimension

  ! iterator variables
  INTEGER(INTG) :: I, I0
  CHARACTER(LEN=30) :: I0_STR

  ! null pointers


  TYPE(DOMAIN_MAPPING_TYPE), POINTER :: Variable2



  ! return if maximum allowed depth is reached
  IF (DEPTH > MaxDepth) THEN
    RETURN
  END IF

  PrintIndent = TRIM(GetPrintIndent(Depth))
  Variable2 => Variable


  PRINT*, TRIM(PrintIndent),"INTEGER(INTG) :: " // &
    & "NUMBER_OF_LOCAL:                 ", &
    & Variable%NUMBER_OF_LOCAL
  PRINT*, TRIM(PrintIndent),"INTEGER(INTG) :: " // &
    & "TOTAL_NUMBER_OF_LOCAL:           ", &
    & Variable%TOTAL_NUMBER_OF_LOCAL
  
  ! Variable%NUMBER_OF_DOMAIN_LOCAL(:)
  ! case allocatable
  IsAllocated = ALLOCATED(Variable%NUMBER_OF_DOMAIN_LOCAL)
  IF (IsAllocated) THEN 
    Dimension = SIZE(SHAPE(Variable%NUMBER_OF_DOMAIN_LOCAL),1)
    WRITE(*,'(3A,I2,A)',advance='no') &
      & " ",TRIM(PrintIndent), "INTEGER(INTG), ALLOCATABLE :: " // &
      & "NUMBER_OF_DOMAIN_LOCAL(:) " // &
      & "(allocated, rank=", Dimension, ", size="

    ! loop over dimensions
    DO I = 1,Dimension
      WRITE(*,'(I5)',advance='no') &
        & SIZE(Variable%NUMBER_OF_DOMAIN_LOCAL, I)
      IF (I /= Dimension) WRITE(*,'(A)',advance='no') " x "
    ENDDO
    WRITE(*,'(A)',advance='no') "): "
    PRINT*," "
    DO I0 = LBOUND(Variable%NUMBER_OF_DOMAIN_LOCAL,1), MIN(LBOUND(Variable% &
      & NUMBER_OF_DOMAIN_LOCAL,1)+MaxArrayLength-1, UBOUND(Variable%NUMBER_OF_DOMAIN_LOCAL,1))
      WRITE(I0_STR,"(I4)") I0 

      PRINT*, TRIM(PrintIndent), TRIM(I0_STR), ":", Variable%NUMBER_OF_DOMAIN_LOCAL(I0)

    ENDDO  ! I0
  ELSE
    PRINT*, TRIM(PrintIndent),"INTEGER(INTG), ALLOCATABLE :: " // &
      & "NUMBER_OF_DOMAIN_LOCAL(:) (not allocated)"
  ENDIF ! IF (IsAllocated)
  
  ! Variable%NUMBER_OF_DOMAIN_GHOST(:)
  ! case allocatable
  IsAllocated = ALLOCATED(Variable%NUMBER_OF_DOMAIN_GHOST)
  IF (IsAllocated) THEN 
    Dimension = SIZE(SHAPE(Variable%NUMBER_OF_DOMAIN_GHOST),1)
    WRITE(*,'(3A,I2,A)',advance='no') &
      & " ",TRIM(PrintIndent), "INTEGER(INTG), ALLOCATABLE :: " // &
      & "NUMBER_OF_DOMAIN_GHOST(:) " // &
      & "(allocated, rank=", Dimension, ", size="

    ! loop over dimensions
    DO I = 1,Dimension
      WRITE(*,'(I5)',advance='no') &
        & SIZE(Variable%NUMBER_OF_DOMAIN_GHOST, I)
      IF (I /= Dimension) WRITE(*,'(A)',advance='no') " x "
    ENDDO
    WRITE(*,'(A)',advance='no') "): "
    PRINT*," "
    DO I0 = LBOUND(Variable%NUMBER_OF_DOMAIN_GHOST,1), MIN(LBOUND(Variable% &
      & NUMBER_OF_DOMAIN_GHOST,1)+MaxArrayLength-1, UBOUND(Variable%NUMBER_OF_DOMAIN_GHOST,1))
      WRITE(I0_STR,"(I4)") I0 

      PRINT*, TRIM(PrintIndent), TRIM(I0_STR), ":", Variable%NUMBER_OF_DOMAIN_GHOST(I0)

    ENDDO  ! I0
  ELSE
    PRINT*, TRIM(PrintIndent),"INTEGER(INTG), ALLOCATABLE :: " // &
      & "NUMBER_OF_DOMAIN_GHOST(:) (not allocated)"
  ENDIF ! IF (IsAllocated)
  PRINT*, TRIM(PrintIndent),"INTEGER(INTG) :: " // &
    & "NUMBER_OF_GLOBAL:                ", &
    & Variable%NUMBER_OF_GLOBAL
  PRINT*, TRIM(PrintIndent),"INTEGER(INTG) :: " // &
    & "NUMBER_OF_DOMAINS:               ", &
    & Variable%NUMBER_OF_DOMAINS
  PRINT*, TRIM(PrintIndent),"INTEGER(INTG) :: " // &
    & "NUMBER_OF_INTERNAL:              ", &
    & Variable%NUMBER_OF_INTERNAL
  PRINT*, TRIM(PrintIndent),"INTEGER(INTG) :: " // &
    & "NUMBER_OF_BOUNDARY:              ", &
    & Variable%NUMBER_OF_BOUNDARY
  PRINT*, TRIM(PrintIndent),"INTEGER(INTG) :: " // &
    & "NUMBER_OF_GHOST:                 ", &
    & Variable%NUMBER_OF_GHOST
  PRINT*, TRIM(PrintIndent),"INTEGER(INTG) :: " // &
    & "INTERNAL_START:                  ", &
    & Variable%INTERNAL_START
  PRINT*, TRIM(PrintIndent),"INTEGER(INTG) :: " // &
    & "INTERNAL_FINISH:                 ", &
    & Variable%INTERNAL_FINISH
  PRINT*, TRIM(PrintIndent),"INTEGER(INTG) :: " // &
    & "BOUNDARY_START:                  ", &
    & Variable%BOUNDARY_START
  PRINT*, TRIM(PrintIndent),"INTEGER(INTG) :: " // &
    & "BOUNDARY_FINISH:                 ", &
    & Variable%BOUNDARY_FINISH
  PRINT*, TRIM(PrintIndent),"INTEGER(INTG) :: " // &
    & "GHOST_START:                     ", &
    & Variable%GHOST_START
  PRINT*, TRIM(PrintIndent),"INTEGER(INTG) :: " // &
    & "GHOST_FINISH:                    ", &
    & Variable%GHOST_FINISH
  
  ! Variable%DOMAIN_LIST(:)
  ! case allocatable
  IsAllocated = ALLOCATED(Variable%DOMAIN_LIST)
  IF (IsAllocated) THEN 
    Dimension = SIZE(SHAPE(Variable%DOMAIN_LIST),1)
    WRITE(*,'(3A,I2,A)',advance='no') &
      & " ",TRIM(PrintIndent), "INTEGER(INTG), ALLOCATABLE :: " // &
      & "DOMAIN_LIST(:) " // &
      & "(allocated, rank=", Dimension, ", size="

    ! loop over dimensions
    DO I = 1,Dimension
      WRITE(*,'(I5)',advance='no') &
        & SIZE(Variable%DOMAIN_LIST, I)
      IF (I /= Dimension) WRITE(*,'(A)',advance='no') " x "
    ENDDO
    WRITE(*,'(A)',advance='no') "): "
    PRINT*," "
    DO I0 = LBOUND(Variable%DOMAIN_LIST,1), MIN(LBOUND(Variable%DOMAIN_LIST,1)+MaxArrayLength-1, UBOUND(Variable%DOMAIN_LIST,1))
      WRITE(I0_STR,"(I4)") I0 

      PRINT*, TRIM(PrintIndent), TRIM(I0_STR), ":", Variable%DOMAIN_LIST(I0)

    ENDDO  ! I0
  ELSE
    PRINT*, TRIM(PrintIndent),"INTEGER(INTG), ALLOCATABLE :: " // &
      & "DOMAIN_LIST(:) (not allocated)"
  ENDIF ! IF (IsAllocated)
  
  ! Variable%LOCAL_TO_GLOBAL_MAP(:)
  ! case allocatable
  IsAllocated = ALLOCATED(Variable%LOCAL_TO_GLOBAL_MAP)
  IF (IsAllocated) THEN 
    Dimension = SIZE(SHAPE(Variable%LOCAL_TO_GLOBAL_MAP),1)
    WRITE(*,'(3A,I2,A)',advance='no') &
      & " ",TRIM(PrintIndent), "INTEGER(INTG), ALLOCATABLE :: " // &
      & "LOCAL_TO_GLOBAL_MAP(:) " // &
      & "(allocated, rank=", Dimension, ", size="

    ! loop over dimensions
    DO I = 1,Dimension
      WRITE(*,'(I5)',advance='no') &
        & SIZE(Variable%LOCAL_TO_GLOBAL_MAP, I)
      IF (I /= Dimension) WRITE(*,'(A)',advance='no') " x "
    ENDDO
    WRITE(*,'(A)',advance='no') "): "
    PRINT*," "
    DO I0 = LBOUND(Variable%LOCAL_TO_GLOBAL_MAP,1), MIN(LBOUND(Variable% &
      & LOCAL_TO_GLOBAL_MAP,1)+MaxArrayLength-1, UBOUND(Variable%LOCAL_TO_GLOBAL_MAP,1))
      WRITE(I0_STR,"(I4)") I0 

      PRINT*, TRIM(PrintIndent), TRIM(I0_STR), ":", Variable%LOCAL_TO_GLOBAL_MAP(I0)

    ENDDO  ! I0
  ELSE
    PRINT*, TRIM(PrintIndent),"INTEGER(INTG), ALLOCATABLE :: " // &
      & "LOCAL_TO_GLOBAL_MAP(:) (not allocated)"
  ENDIF ! IF (IsAllocated)
  PRINT*, TRIM(PrintIndent),"INTEGER(INTG) :: " // &
    & "NUMBER_OF_ADJACENT_DOMAINS:      ", &
    & Variable%NUMBER_OF_ADJACENT_DOMAINS
  
  ! Variable%ADJACENT_DOMAINS(:)
  ! case allocatable
  IsAllocated = ALLOCATED(Variable%ADJACENT_DOMAINS)
  IF (IsAllocated) THEN 
    Dimension = SIZE(SHAPE(Variable%ADJACENT_DOMAINS),1)
    WRITE(*,'(3A,I2,A)',advance='no') &
      & " ",TRIM(PrintIndent), "TYPE(DOMAIN_ADJACENT_DOMAIN_TYPE), ALLOCATABLE :: " // &
      & "ADJACENT_DOMAINS(:) " // &
      & "(allocated, rank=", Dimension, ", size="

    ! loop over dimensions
    DO I = 1,Dimension
      WRITE(*,'(I5)',advance='no') &
        & SIZE(Variable%ADJACENT_DOMAINS, I)
      IF (I /= Dimension) WRITE(*,'(A)',advance='no') " x "
    ENDDO
    WRITE(*,'(A)',advance='no') "): "
  
    DO I0 = LBOUND(Variable%ADJACENT_DOMAINS,1), MIN(LBOUND(Variable%ADJACENT_DOMAINS,1)+MaxArrayLength-1, UBOUND(Variable% &
      & ADJACENT_DOMAINS,1))
      WRITE(I0_STR,"(I4)") I0 

! Signature of Print_DOMAIN_ADJACENT_DOMAIN_TYPE has the following CheckVariable types: 
! The type to handle in this routine is: DOMAIN_MAPPING_TYPE (Variable)
! Available CheckVariables are: 
        PRINT*, NEW_LINE('A') // " " // TRIM(PrintIndent),"ADJACENT_DOMAINS("//TRIM(ADJUSTL(I0_STR))//")"
        CALL Print_DOMAIN_ADJACENT_DOMAIN_TYPE(Variable%ADJACENT_DOMAINS(I0), Depth+1, MaxDepth, MaxArrayLength)

    ENDDO  ! I0
  ELSE
    PRINT*, TRIM(PrintIndent),"TYPE(DOMAIN_ADJACENT_DOMAIN_TYPE), ALLOCATABLE :: " // &
      & "ADJACENT_DOMAINS(:) (not allocated)"
  ENDIF ! IF (IsAllocated)
  
  ! Variable%ADJACENT_DOMAINS_PTR(:)
  ! case allocatable
  IsAllocated = ALLOCATED(Variable%ADJACENT_DOMAINS_PTR)
  IF (IsAllocated) THEN 
    Dimension = SIZE(SHAPE(Variable%ADJACENT_DOMAINS_PTR),1)
    WRITE(*,'(3A,I2,A)',advance='no') &
      & " ",TRIM(PrintIndent), "INTEGER(INTG), ALLOCATABLE :: " // &
      & "ADJACENT_DOMAINS_PTR(:) " // &
      & "(allocated, rank=", Dimension, ", size="

    ! loop over dimensions
    DO I = 1,Dimension
      WRITE(*,'(I5)',advance='no') &
        & SIZE(Variable%ADJACENT_DOMAINS_PTR, I)
      IF (I /= Dimension) WRITE(*,'(A)',advance='no') " x "
    ENDDO
    WRITE(*,'(A)',advance='no') "): "
    PRINT*," "
    DO I0 = LBOUND(Variable%ADJACENT_DOMAINS_PTR,1), MIN(LBOUND(Variable% &
      & ADJACENT_DOMAINS_PTR,1)+MaxArrayLength-1, UBOUND(Variable%ADJACENT_DOMAINS_PTR,1))
      WRITE(I0_STR,"(I4)") I0 

      PRINT*, TRIM(PrintIndent), TRIM(I0_STR), ":", Variable%ADJACENT_DOMAINS_PTR(I0)

    ENDDO  ! I0
  ELSE
    PRINT*, TRIM(PrintIndent),"INTEGER(INTG), ALLOCATABLE :: " // &
      & "ADJACENT_DOMAINS_PTR(:) (not allocated)"
  ENDIF ! IF (IsAllocated)
  
  ! Variable%ADJACENT_DOMAINS_LIST(:)
  ! case allocatable
  IsAllocated = ALLOCATED(Variable%ADJACENT_DOMAINS_LIST)
  IF (IsAllocated) THEN 
    Dimension = SIZE(SHAPE(Variable%ADJACENT_DOMAINS_LIST),1)
    WRITE(*,'(3A,I2,A)',advance='no') &
      & " ",TRIM(PrintIndent), "INTEGER(INTG), ALLOCATABLE :: " // &
      & "ADJACENT_DOMAINS_LIST(:) " // &
      & "(allocated, rank=", Dimension, ", size="

    ! loop over dimensions
    DO I = 1,Dimension
      WRITE(*,'(I5)',advance='no') &
        & SIZE(Variable%ADJACENT_DOMAINS_LIST, I)
      IF (I /= Dimension) WRITE(*,'(A)',advance='no') " x "
    ENDDO
    WRITE(*,'(A)',advance='no') "): "
    PRINT*," "
    DO I0 = LBOUND(Variable%ADJACENT_DOMAINS_LIST,1), MIN(LBOUND(Variable% &
      & ADJACENT_DOMAINS_LIST,1)+MaxArrayLength-1, UBOUND(Variable%ADJACENT_DOMAINS_LIST,1))
      WRITE(I0_STR,"(I4)") I0 

      PRINT*, TRIM(PrintIndent), TRIM(I0_STR), ":", Variable%ADJACENT_DOMAINS_LIST(I0)

    ENDDO  ! I0
  ELSE
    PRINT*, TRIM(PrintIndent),"INTEGER(INTG), ALLOCATABLE :: " // &
      & "ADJACENT_DOMAINS_LIST(:) (not allocated)"
  ENDIF ! IF (IsAllocated)
  
  ! Variable%GLOBAL_TO_LOCAL_MAP(:)
  ! case allocatable
  IsAllocated = ALLOCATED(Variable%GLOBAL_TO_LOCAL_MAP)
  IF (IsAllocated) THEN 
    Dimension = SIZE(SHAPE(Variable%GLOBAL_TO_LOCAL_MAP),1)
    WRITE(*,'(3A,I2,A)',advance='no') &
      & " ",TRIM(PrintIndent), "TYPE(DOMAIN_GLOBAL_MAPPING_TYPE), ALLOCATABLE :: " // &
      & "GLOBAL_TO_LOCAL_MAP(:) " // &
      & "(allocated, rank=", Dimension, ", size="

    ! loop over dimensions
    DO I = 1,Dimension
      WRITE(*,'(I5)',advance='no') &
        & SIZE(Variable%GLOBAL_TO_LOCAL_MAP, I)
      IF (I /= Dimension) WRITE(*,'(A)',advance='no') " x "
    ENDDO
    WRITE(*,'(A)',advance='no') "): "
  
    DO I0 = LBOUND(Variable%GLOBAL_TO_LOCAL_MAP,1), MIN(LBOUND(Variable% &
      & GLOBAL_TO_LOCAL_MAP,1)+MaxArrayLength-1, UBOUND(Variable%GLOBAL_TO_LOCAL_MAP,1))
      WRITE(I0_STR,"(I4)") I0 

! Signature of Print_DOMAIN_GLOBAL_MAPPING_TYPE has the following CheckVariable types: 
! The type to handle in this routine is: DOMAIN_MAPPING_TYPE (Variable)
! Available CheckVariables are: 
        PRINT*, NEW_LINE('A') // " " // TRIM(PrintIndent),"GLOBAL_TO_LOCAL_MAP("//TRIM(ADJUSTL(I0_STR))//")"
        CALL Print_DOMAIN_GLOBAL_MAPPING_TYPE(Variable%GLOBAL_TO_LOCAL_MAP(I0), Depth+1, MaxDepth, MaxArrayLength)

    ENDDO  ! I0
  ELSE
    PRINT*, TRIM(PrintIndent),"TYPE(DOMAIN_GLOBAL_MAPPING_TYPE), ALLOCATABLE :: " // &
      & "GLOBAL_TO_LOCAL_MAP(:) (not allocated)"
  ENDIF ! IF (IsAllocated)

  
END SUBROUTINE Print_DOMAIN_MAPPING_TYPE

!
!================================================================================================================================
!
RECURSIVE SUBROUTINE Print_DOMAIN_ADJACENT_DOMAIN_TYPE(Variable, Depth, MaxDepth, MaxArrayLength)
  TYPE(DOMAIN_ADJACENT_DOMAIN_TYPE), POINTER, INTENT(IN) :: Variable
  INTEGER(INTG), INTENT(IN) :: Depth   !< the recursion depth
  INTEGER(INTG), INTENT(IN) :: MaxDepth   !< the maximum recursion depth for which data is printed
  INTEGER(INTG), INTENT(IN) :: MaxArrayLength   !< the maximum array length that is printed  
  CHARACTER(LEN=400) :: PrintIndent   !< a string containing 2*Depth space characters
  LOGICAL :: IsAllocated
  LOGICAL :: IsAssociated
  INTEGER(INTG) :: Dimension

  ! iterator variables
  INTEGER(INTG) :: I, I0
  CHARACTER(LEN=30) :: I0_STR

  ! null pointers


  TYPE(DOMAIN_ADJACENT_DOMAIN_TYPE), POINTER :: Variable2



  ! return if maximum allowed depth is reached
  IF (DEPTH > MaxDepth) THEN
    RETURN
  END IF

  PrintIndent = TRIM(GetPrintIndent(Depth))
  Variable2 => Variable


  PRINT*, TRIM(PrintIndent),"INTEGER(INTG) :: " // &
    & "DOMAIN_NUMBER:                   ", &
    & Variable%DOMAIN_NUMBER
  PRINT*, TRIM(PrintIndent),"INTEGER(INTG) :: " // &
    & "NUMBER_OF_SEND_GHOSTS:           ", &
    & Variable%NUMBER_OF_SEND_GHOSTS
  PRINT*, TRIM(PrintIndent),"INTEGER(INTG) :: " // &
    & "NUMBER_OF_RECEIVE_GHOSTS:        ", &
    & Variable%NUMBER_OF_RECEIVE_GHOSTS
  
  ! Variable%LOCAL_GHOST_SEND_INDICES(:)
  ! case allocatable
  IsAllocated = ALLOCATED(Variable%LOCAL_GHOST_SEND_INDICES)
  IF (IsAllocated) THEN 
    Dimension = SIZE(SHAPE(Variable%LOCAL_GHOST_SEND_INDICES),1)
    WRITE(*,'(3A,I2,A)',advance='no') &
      & " ",TRIM(PrintIndent), "INTEGER(INTG), ALLOCATABLE :: " // &
      & "LOCAL_GHOST_SEND_INDICES(:) " // &
      & "(allocated, rank=", Dimension, ", size="

    ! loop over dimensions
    DO I = 1,Dimension
      WRITE(*,'(I5)',advance='no') &
        & SIZE(Variable%LOCAL_GHOST_SEND_INDICES, I)
      IF (I /= Dimension) WRITE(*,'(A)',advance='no') " x "
    ENDDO
    WRITE(*,'(A)',advance='no') "): "
    PRINT*," "
    DO I0 = LBOUND(Variable%LOCAL_GHOST_SEND_INDICES,1), MIN(LBOUND(Variable% &
      & LOCAL_GHOST_SEND_INDICES,1)+MaxArrayLength-1, UBOUND(Variable%LOCAL_GHOST_SEND_INDICES,1))
      WRITE(I0_STR,"(I4)") I0 

      PRINT*, TRIM(PrintIndent), TRIM(I0_STR), ":", Variable%LOCAL_GHOST_SEND_INDICES(I0)

    ENDDO  ! I0
  ELSE
    PRINT*, TRIM(PrintIndent),"INTEGER(INTG), ALLOCATABLE :: " // &
      & "LOCAL_GHOST_SEND_INDICES(:) (not allocated)"
  ENDIF ! IF (IsAllocated)
  
  ! Variable%LOCAL_GHOST_RECEIVE_INDICES(:)
  ! case allocatable
  IsAllocated = ALLOCATED(Variable%LOCAL_GHOST_RECEIVE_INDICES)
  IF (IsAllocated) THEN 
    Dimension = SIZE(SHAPE(Variable%LOCAL_GHOST_RECEIVE_INDICES),1)
    WRITE(*,'(3A,I2,A)',advance='no') &
      & " ",TRIM(PrintIndent), "INTEGER(INTG), ALLOCATABLE :: " // &
      & "LOCAL_GHOST_RECEIVE_INDICES(:) " // &
      & "(allocated, rank=", Dimension, ", size="

    ! loop over dimensions
    DO I = 1,Dimension
      WRITE(*,'(I5)',advance='no') &
        & SIZE(Variable%LOCAL_GHOST_RECEIVE_INDICES, I)
      IF (I /= Dimension) WRITE(*,'(A)',advance='no') " x "
    ENDDO
    WRITE(*,'(A)',advance='no') "): "
    PRINT*," "
    DO I0 = LBOUND(Variable%LOCAL_GHOST_RECEIVE_INDICES,1), MIN(LBOUND(Variable% &
      & LOCAL_GHOST_RECEIVE_INDICES,1)+MaxArrayLength-1, UBOUND(Variable%LOCAL_GHOST_RECEIVE_INDICES,1))
      WRITE(I0_STR,"(I4)") I0 

      PRINT*, TRIM(PrintIndent), TRIM(I0_STR), ":", Variable%LOCAL_GHOST_RECEIVE_INDICES(I0)

    ENDDO  ! I0
  ELSE
    PRINT*, TRIM(PrintIndent),"INTEGER(INTG), ALLOCATABLE :: " // &
      & "LOCAL_GHOST_RECEIVE_INDICES(:) (not allocated)"
  ENDIF ! IF (IsAllocated)
  PRINT*, TRIM(PrintIndent),"INTEGER(INTG) :: " // &
    & "NUMBER_OF_FURTHER_LINKED_GHOSTS: ", &
    & Variable%NUMBER_OF_FURTHER_LINKED_GHOSTS
  
  ! Variable%LOCAL_GHOST_FURTHER_INDICES(:)
  ! case allocatable
  IsAllocated = ALLOCATED(Variable%LOCAL_GHOST_FURTHER_INDICES)
  IF (IsAllocated) THEN 
    Dimension = SIZE(SHAPE(Variable%LOCAL_GHOST_FURTHER_INDICES),1)
    WRITE(*,'(3A,I2,A)',advance='no') &
      & " ",TRIM(PrintIndent), "INTEGER(INTG), ALLOCATABLE :: " // &
      & "LOCAL_GHOST_FURTHER_INDICES(:) " // &
      & "(allocated, rank=", Dimension, ", size="

    ! loop over dimensions
    DO I = 1,Dimension
      WRITE(*,'(I5)',advance='no') &
        & SIZE(Variable%LOCAL_GHOST_FURTHER_INDICES, I)
      IF (I /= Dimension) WRITE(*,'(A)',advance='no') " x "
    ENDDO
    WRITE(*,'(A)',advance='no') "): "
    PRINT*," "
    DO I0 = LBOUND(Variable%LOCAL_GHOST_FURTHER_INDICES,1), MIN(LBOUND(Variable% &
      & LOCAL_GHOST_FURTHER_INDICES,1)+MaxArrayLength-1, UBOUND(Variable%LOCAL_GHOST_FURTHER_INDICES,1))
      WRITE(I0_STR,"(I4)") I0 

      PRINT*, TRIM(PrintIndent), TRIM(I0_STR), ":", Variable%LOCAL_GHOST_FURTHER_INDICES(I0)

    ENDDO  ! I0
  ELSE
    PRINT*, TRIM(PrintIndent),"INTEGER(INTG), ALLOCATABLE :: " // &
      & "LOCAL_GHOST_FURTHER_INDICES(:) (not allocated)"
  ENDIF ! IF (IsAllocated)

  
END SUBROUTINE Print_DOMAIN_ADJACENT_DOMAIN_TYPE

!
!================================================================================================================================
!
RECURSIVE SUBROUTINE Print_DOMAIN_GLOBAL_MAPPING_TYPE(Variable, Depth, MaxDepth, MaxArrayLength)
  TYPE(DOMAIN_GLOBAL_MAPPING_TYPE), POINTER, INTENT(IN) :: Variable
  INTEGER(INTG), INTENT(IN) :: Depth   !< the recursion depth
  INTEGER(INTG), INTENT(IN) :: MaxDepth   !< the maximum recursion depth for which data is printed
  INTEGER(INTG), INTENT(IN) :: MaxArrayLength   !< the maximum array length that is printed  
  CHARACTER(LEN=400) :: PrintIndent   !< a string containing 2*Depth space characters
  LOGICAL :: IsAllocated
  LOGICAL :: IsAssociated
  INTEGER(INTG) :: Dimension

  ! iterator variables
  INTEGER(INTG) :: I, I0
  CHARACTER(LEN=30) :: I0_STR

  ! null pointers


  TYPE(DOMAIN_GLOBAL_MAPPING_TYPE), POINTER :: Variable2



  ! return if maximum allowed depth is reached
  IF (DEPTH > MaxDepth) THEN
    RETURN
  END IF

  PrintIndent = TRIM(GetPrintIndent(Depth))
  Variable2 => Variable


  PRINT*, TRIM(PrintIndent),"INTEGER(INTG) :: " // &
    & "NUMBER_OF_DOMAINS:               ", &
    & Variable%NUMBER_OF_DOMAINS
  
  ! Variable%LOCAL_NUMBER(:)
  ! case allocatable
  IsAllocated = ALLOCATED(Variable%LOCAL_NUMBER)
  IF (IsAllocated) THEN 
    Dimension = SIZE(SHAPE(Variable%LOCAL_NUMBER),1)
    WRITE(*,'(3A,I2,A)',advance='no') &
      & " ",TRIM(PrintIndent), "INTEGER(INTG), ALLOCATABLE :: " // &
      & "LOCAL_NUMBER(:) " // &
      & "(allocated, rank=", Dimension, ", size="

    ! loop over dimensions
    DO I = 1,Dimension
      WRITE(*,'(I5)',advance='no') &
        & SIZE(Variable%LOCAL_NUMBER, I)
      IF (I /= Dimension) WRITE(*,'(A)',advance='no') " x "
    ENDDO
    WRITE(*,'(A)',advance='no') "): "
    PRINT*," "
    DO I0 = LBOUND(Variable%LOCAL_NUMBER,1), MIN(LBOUND(Variable%LOCAL_NUMBER,1)+MaxArrayLength-1, UBOUND(Variable%LOCAL_NUMBER,1))
      WRITE(I0_STR,"(I4)") I0 

      PRINT*, TRIM(PrintIndent), TRIM(I0_STR), ":", Variable%LOCAL_NUMBER(I0)

    ENDDO  ! I0
  ELSE
    PRINT*, TRIM(PrintIndent),"INTEGER(INTG), ALLOCATABLE :: " // &
      & "LOCAL_NUMBER(:) (not allocated)"
  ENDIF ! IF (IsAllocated)
  
  ! Variable%DOMAIN_NUMBER(:)
  ! case allocatable
  IsAllocated = ALLOCATED(Variable%DOMAIN_NUMBER)
  IF (IsAllocated) THEN 
    Dimension = SIZE(SHAPE(Variable%DOMAIN_NUMBER),1)
    WRITE(*,'(3A,I2,A)',advance='no') &
      & " ",TRIM(PrintIndent), "INTEGER(INTG), ALLOCATABLE :: " // &
      & "DOMAIN_NUMBER(:) " // &
      & "(allocated, rank=", Dimension, ", size="

    ! loop over dimensions
    DO I = 1,Dimension
      WRITE(*,'(I5)',advance='no') &
        & SIZE(Variable%DOMAIN_NUMBER, I)
      IF (I /= Dimension) WRITE(*,'(A)',advance='no') " x "
    ENDDO
    WRITE(*,'(A)',advance='no') "): "
    PRINT*," "
    DO I0 = LBOUND(Variable%DOMAIN_NUMBER,1), MIN(LBOUND(Variable%DOMAIN_NUMBER,1)+MaxArrayLength-1, UBOUND(Variable% &
      & DOMAIN_NUMBER,1))
      WRITE(I0_STR,"(I4)") I0 

      PRINT*, TRIM(PrintIndent), TRIM(I0_STR), ":", Variable%DOMAIN_NUMBER(I0)

    ENDDO  ! I0
  ELSE
    PRINT*, TRIM(PrintIndent),"INTEGER(INTG), ALLOCATABLE :: " // &
      & "DOMAIN_NUMBER(:) (not allocated)"
  ENDIF ! IF (IsAllocated)
  
  ! Variable%LOCAL_TYPE(:)
  ! case allocatable
  IsAllocated = ALLOCATED(Variable%LOCAL_TYPE)
  IF (IsAllocated) THEN 
    Dimension = SIZE(SHAPE(Variable%LOCAL_TYPE),1)
    WRITE(*,'(3A,I2,A)',advance='no') &
      & " ",TRIM(PrintIndent), "INTEGER(INTG), ALLOCATABLE :: " // &
      & "LOCAL_TYPE(:) " // &
      & "(allocated, rank=", Dimension, ", size="

    ! loop over dimensions
    DO I = 1,Dimension
      WRITE(*,'(I5)',advance='no') &
        & SIZE(Variable%LOCAL_TYPE, I)
      IF (I /= Dimension) WRITE(*,'(A)',advance='no') " x "
    ENDDO
    WRITE(*,'(A)',advance='no') "): "
    PRINT*," "
    DO I0 = LBOUND(Variable%LOCAL_TYPE,1), MIN(LBOUND(Variable%LOCAL_TYPE,1)+MaxArrayLength-1, UBOUND(Variable%LOCAL_TYPE,1))
      WRITE(I0_STR,"(I4)") I0 

      PRINT*, TRIM(PrintIndent), TRIM(I0_STR), ":", Variable%LOCAL_TYPE(I0)

    ENDDO  ! I0
  ELSE
    PRINT*, TRIM(PrintIndent),"INTEGER(INTG), ALLOCATABLE :: " // &
      & "LOCAL_TYPE(:) (not allocated)"
  ENDIF ! IF (IsAllocated)

  
END SUBROUTINE Print_DOMAIN_GLOBAL_MAPPING_TYPE

  !
  !================================================================================================================================
  !
  SUBROUTINE Print_FIELD_DOF_TO_PARAM_MAP(Variable, MaxDepth, MaxArrayLength)
    TYPE(FIELD_DOF_TO_PARAM_MAP_TYPE), POINTER, INTENT(IN) :: Variable   !< the variable to be printed
    INTEGER(INTG), INTENT(IN) :: MaxDepth           !< the maximum recursion depth down to which data is printed
    INTEGER(INTG), INTENT(IN) :: MaxArrayLength   !< the maximum array length that is printed  


    PRINT*, "Print TYPE(FIELD_DOF_TO_PARAM_MAP_TYPE), POINTER :: Variable with maximum depth ",MaxDepth,"," // &
      & "maximum array length:",MaxArrayLength
    CALL Print_FIELD_DOF_TO_PARAM_MAP_TYPE(Variable, 1, MaxDepth, MaxArrayLength)
    
  END SUBROUTINE Print_FIELD_DOF_TO_PARAM_MAP

!
!================================================================================================================================
!
RECURSIVE SUBROUTINE Print_FIELD_DOF_TO_PARAM_MAP_TYPE(Variable, Depth, MaxDepth, MaxArrayLength)
  TYPE(FIELD_DOF_TO_PARAM_MAP_TYPE), POINTER, INTENT(IN) :: Variable
  INTEGER(INTG), INTENT(IN) :: Depth   !< the recursion depth
  INTEGER(INTG), INTENT(IN) :: MaxDepth   !< the maximum recursion depth for which data is printed
  INTEGER(INTG), INTENT(IN) :: MaxArrayLength   !< the maximum array length that is printed  
  CHARACTER(LEN=400) :: PrintIndent   !< a string containing 2*Depth space characters
  LOGICAL :: IsAllocated
  LOGICAL :: IsAssociated
  INTEGER(INTG) :: Dimension

  ! iterator variables
  INTEGER(INTG) :: I, I0, I1
  CHARACTER(LEN=30) :: I0_STR, I1_STR

  ! null pointers


  TYPE(FIELD_DOF_TO_PARAM_MAP_TYPE), POINTER :: Variable2



  ! return if maximum allowed depth is reached
  IF (DEPTH > MaxDepth) THEN
    RETURN
  END IF

  PrintIndent = TRIM(GetPrintIndent(Depth))
  Variable2 => Variable


  PRINT*, TRIM(PrintIndent),"INTEGER(INTG) :: " // &
    & "NUMBER_OF_DOFS:                  ", &
    & Variable%NUMBER_OF_DOFS
  
  ! Variable%DOF_TYPE(:,:)
  ! case allocatable
  IsAllocated = ALLOCATED(Variable%DOF_TYPE)
  IF (IsAllocated) THEN 
    Dimension = SIZE(SHAPE(Variable%DOF_TYPE),1)
    WRITE(*,'(3A,I2,A)',advance='no') &
      & " ",TRIM(PrintIndent), "INTEGER(INTG), ALLOCATABLE :: " // &
      & "DOF_TYPE(:,:) " // &
      & "(allocated, rank=", Dimension, ", size="

    ! loop over dimensions
    DO I = 1,Dimension
      WRITE(*,'(I5)',advance='no') &
        & SIZE(Variable%DOF_TYPE, I)
      IF (I /= Dimension) WRITE(*,'(A)',advance='no') " x "
    ENDDO
    WRITE(*,'(A)',advance='no') "): "
    PRINT*," "
    DO I0 = LBOUND(Variable%DOF_TYPE,1), MIN(LBOUND(Variable%DOF_TYPE,1)+MaxArrayLength-1, UBOUND(Variable%DOF_TYPE,1))
      WRITE(I0_STR,"(I4)") I0 

      PRINT*, TRIM(PrintIndent),"("//TRIM(ADJUSTL(I0_STR))//",*): "
      DO I1 = LBOUND(Variable%DOF_TYPE,2), MIN(LBOUND(Variable%DOF_TYPE,2)+MaxArrayLength-1, UBOUND(Variable%DOF_TYPE,2))
        WRITE(I1_STR,"(I4)") I1 

        PRINT*, TRIM(PrintIndent), TRIM(I1_STR), ":", Variable%DOF_TYPE(I0,I1)

      ENDDO  ! I1
    ENDDO  ! I0
  ELSE
    PRINT*, TRIM(PrintIndent),"INTEGER(INTG), ALLOCATABLE :: " // &
      & "DOF_TYPE(:,:) (not allocated)"
  ENDIF ! IF (IsAllocated)
  PRINT*, TRIM(PrintIndent),"INTEGER(INTG) :: " // &
    & "NUMBER_OF_CONSTANT_DOFS:         ", &
    & Variable%NUMBER_OF_CONSTANT_DOFS
  PRINT*, TRIM(PrintIndent),"INTEGER(INTG) :: " // &
    & "NUMBER_OF_ELEMENT_DOFS:          ", &
    & Variable%NUMBER_OF_ELEMENT_DOFS
  PRINT*, TRIM(PrintIndent),"INTEGER(INTG) :: " // &
    & "NUMBER_OF_NODE_DOFS:             ", &
    & Variable%NUMBER_OF_NODE_DOFS
  PRINT*, TRIM(PrintIndent),"INTEGER(INTG) :: " // &
    & "NUMBER_OF_GRID_POINT_DOFS:       ", &
    & Variable%NUMBER_OF_GRID_POINT_DOFS
  PRINT*, TRIM(PrintIndent),"INTEGER(INTG) :: " // &
    & "NUMBER_OF_GAUSS_POINT_DOFS:      ", &
    & Variable%NUMBER_OF_GAUSS_POINT_DOFS
  PRINT*, TRIM(PrintIndent),"INTEGER(INTG) :: " // &
    & "NUMBER_OF_DATA_POINT_DOFS:       ", &
    & Variable%NUMBER_OF_DATA_POINT_DOFS
  
  ! Variable%CONSTANT_DOF2PARAM_MAP(:)
  ! case allocatable
  IsAllocated = ALLOCATED(Variable%CONSTANT_DOF2PARAM_MAP)
  IF (IsAllocated) THEN 
    Dimension = SIZE(SHAPE(Variable%CONSTANT_DOF2PARAM_MAP),1)
    WRITE(*,'(3A,I2,A)',advance='no') &
      & " ",TRIM(PrintIndent), "INTEGER(INTG), ALLOCATABLE :: " // &
      & "CONSTANT_DOF2PARAM_MAP(:) " // &
      & "(allocated, rank=", Dimension, ", size="

    ! loop over dimensions
    DO I = 1,Dimension
      WRITE(*,'(I5)',advance='no') &
        & SIZE(Variable%CONSTANT_DOF2PARAM_MAP, I)
      IF (I /= Dimension) WRITE(*,'(A)',advance='no') " x "
    ENDDO
    WRITE(*,'(A)',advance='no') "): "
    PRINT*," "
    DO I0 = LBOUND(Variable%CONSTANT_DOF2PARAM_MAP,1), MIN(LBOUND(Variable% &
      & CONSTANT_DOF2PARAM_MAP,1)+MaxArrayLength-1, UBOUND(Variable%CONSTANT_DOF2PARAM_MAP,1))
      WRITE(I0_STR,"(I4)") I0 

      PRINT*, TRIM(PrintIndent), TRIM(I0_STR), ":", Variable%CONSTANT_DOF2PARAM_MAP(I0)

    ENDDO  ! I0
  ELSE
    PRINT*, TRIM(PrintIndent),"INTEGER(INTG), ALLOCATABLE :: " // &
      & "CONSTANT_DOF2PARAM_MAP(:) (not allocated)"
  ENDIF ! IF (IsAllocated)
  
  ! Variable%ELEMENT_DOF2PARAM_MAP(:,:)
  ! case allocatable
  IsAllocated = ALLOCATED(Variable%ELEMENT_DOF2PARAM_MAP)
  IF (IsAllocated) THEN 
    Dimension = SIZE(SHAPE(Variable%ELEMENT_DOF2PARAM_MAP),1)
    WRITE(*,'(3A,I2,A)',advance='no') &
      & " ",TRIM(PrintIndent), "INTEGER(INTG), ALLOCATABLE :: " // &
      & "ELEMENT_DOF2PARAM_MAP(:,:) " // &
      & "(allocated, rank=", Dimension, ", size="

    ! loop over dimensions
    DO I = 1,Dimension
      WRITE(*,'(I5)',advance='no') &
        & SIZE(Variable%ELEMENT_DOF2PARAM_MAP, I)
      IF (I /= Dimension) WRITE(*,'(A)',advance='no') " x "
    ENDDO
    WRITE(*,'(A)',advance='no') "): "
    PRINT*," "
    DO I0 = LBOUND(Variable%ELEMENT_DOF2PARAM_MAP,1), MIN(LBOUND(Variable% &
      & ELEMENT_DOF2PARAM_MAP,1)+MaxArrayLength-1, UBOUND(Variable%ELEMENT_DOF2PARAM_MAP,1))
      WRITE(I0_STR,"(I4)") I0 

      PRINT*, TRIM(PrintIndent),"("//TRIM(ADJUSTL(I0_STR))//",*): "
      DO I1 = LBOUND(Variable%ELEMENT_DOF2PARAM_MAP,2), MIN(LBOUND(Variable% &
        & ELEMENT_DOF2PARAM_MAP,2)+MaxArrayLength-1, UBOUND(Variable%ELEMENT_DOF2PARAM_MAP,2))
        WRITE(I1_STR,"(I4)") I1 

        PRINT*, TRIM(PrintIndent), TRIM(I1_STR), ":", Variable%ELEMENT_DOF2PARAM_MAP(I0,I1)

      ENDDO  ! I1
    ENDDO  ! I0
  ELSE
    PRINT*, TRIM(PrintIndent),"INTEGER(INTG), ALLOCATABLE :: " // &
      & "ELEMENT_DOF2PARAM_MAP(:,:) (not allocated)"
  ENDIF ! IF (IsAllocated)
  
  ! Variable%NODE_DOF2PARAM_MAP(:,:)
  ! case allocatable
  IsAllocated = ALLOCATED(Variable%NODE_DOF2PARAM_MAP)
  IF (IsAllocated) THEN 
    Dimension = SIZE(SHAPE(Variable%NODE_DOF2PARAM_MAP),1)
    WRITE(*,'(3A,I2,A)',advance='no') &
      & " ",TRIM(PrintIndent), "INTEGER(INTG), ALLOCATABLE :: " // &
      & "NODE_DOF2PARAM_MAP(:,:) " // &
      & "(allocated, rank=", Dimension, ", size="

    ! loop over dimensions
    DO I = 1,Dimension
      WRITE(*,'(I5)',advance='no') &
        & SIZE(Variable%NODE_DOF2PARAM_MAP, I)
      IF (I /= Dimension) WRITE(*,'(A)',advance='no') " x "
    ENDDO
    WRITE(*,'(A)',advance='no') "): "
    PRINT*," "
    DO I0 = LBOUND(Variable%NODE_DOF2PARAM_MAP,1), MIN(LBOUND(Variable% &
      & NODE_DOF2PARAM_MAP,1)+MaxArrayLength-1, UBOUND(Variable%NODE_DOF2PARAM_MAP,1))
      WRITE(I0_STR,"(I4)") I0 

      PRINT*, TRIM(PrintIndent),"("//TRIM(ADJUSTL(I0_STR))//",*): "
      DO I1 = LBOUND(Variable%NODE_DOF2PARAM_MAP,2), MIN(LBOUND(Variable% &
        & NODE_DOF2PARAM_MAP,2)+MaxArrayLength-1, UBOUND(Variable%NODE_DOF2PARAM_MAP,2))
        WRITE(I1_STR,"(I4)") I1 

        PRINT*, TRIM(PrintIndent), TRIM(I1_STR), ":", Variable%NODE_DOF2PARAM_MAP(I0,I1)

      ENDDO  ! I1
    ENDDO  ! I0
  ELSE
    PRINT*, TRIM(PrintIndent),"INTEGER(INTG), ALLOCATABLE :: " // &
      & "NODE_DOF2PARAM_MAP(:,:) (not allocated)"
  ENDIF ! IF (IsAllocated)
  
  ! Variable%GRID_POINT_DOF2PARAM_MAP(:,:)
  ! case allocatable
  IsAllocated = ALLOCATED(Variable%GRID_POINT_DOF2PARAM_MAP)
  IF (IsAllocated) THEN 
    Dimension = SIZE(SHAPE(Variable%GRID_POINT_DOF2PARAM_MAP),1)
    WRITE(*,'(3A,I2,A)',advance='no') &
      & " ",TRIM(PrintIndent), "INTEGER(INTG), ALLOCATABLE :: " // &
      & "GRID_POINT_DOF2PARAM_MAP(:,:) " // &
      & "(allocated, rank=", Dimension, ", size="

    ! loop over dimensions
    DO I = 1,Dimension
      WRITE(*,'(I5)',advance='no') &
        & SIZE(Variable%GRID_POINT_DOF2PARAM_MAP, I)
      IF (I /= Dimension) WRITE(*,'(A)',advance='no') " x "
    ENDDO
    WRITE(*,'(A)',advance='no') "): "
    PRINT*," "
    DO I0 = LBOUND(Variable%GRID_POINT_DOF2PARAM_MAP,1), MIN(LBOUND(Variable% &
      & GRID_POINT_DOF2PARAM_MAP,1)+MaxArrayLength-1, UBOUND(Variable%GRID_POINT_DOF2PARAM_MAP,1))
      WRITE(I0_STR,"(I4)") I0 

      PRINT*, TRIM(PrintIndent),"("//TRIM(ADJUSTL(I0_STR))//",*): "
      DO I1 = LBOUND(Variable%GRID_POINT_DOF2PARAM_MAP,2), MIN(LBOUND(Variable% &
        & GRID_POINT_DOF2PARAM_MAP,2)+MaxArrayLength-1, UBOUND(Variable%GRID_POINT_DOF2PARAM_MAP,2))
        WRITE(I1_STR,"(I4)") I1 

        PRINT*, TRIM(PrintIndent), TRIM(I1_STR), ":", Variable%GRID_POINT_DOF2PARAM_MAP(I0,I1)

      ENDDO  ! I1
    ENDDO  ! I0
  ELSE
    PRINT*, TRIM(PrintIndent),"INTEGER(INTG), ALLOCATABLE :: " // &
      & "GRID_POINT_DOF2PARAM_MAP(:,:) (not allocated)"
  ENDIF ! IF (IsAllocated)
  
  ! Variable%GAUSS_POINT_DOF2PARAM_MAP(:,:)
  ! case allocatable
  IsAllocated = ALLOCATED(Variable%GAUSS_POINT_DOF2PARAM_MAP)
  IF (IsAllocated) THEN 
    Dimension = SIZE(SHAPE(Variable%GAUSS_POINT_DOF2PARAM_MAP),1)
    WRITE(*,'(3A,I2,A)',advance='no') &
      & " ",TRIM(PrintIndent), "INTEGER(INTG), ALLOCATABLE :: " // &
      & "GAUSS_POINT_DOF2PARAM_MAP(:,:) " // &
      & "(allocated, rank=", Dimension, ", size="

    ! loop over dimensions
    DO I = 1,Dimension
      WRITE(*,'(I5)',advance='no') &
        & SIZE(Variable%GAUSS_POINT_DOF2PARAM_MAP, I)
      IF (I /= Dimension) WRITE(*,'(A)',advance='no') " x "
    ENDDO
    WRITE(*,'(A)',advance='no') "): "
    PRINT*," "
    DO I0 = LBOUND(Variable%GAUSS_POINT_DOF2PARAM_MAP,1), MIN(LBOUND(Variable% &
      & GAUSS_POINT_DOF2PARAM_MAP,1)+MaxArrayLength-1, UBOUND(Variable%GAUSS_POINT_DOF2PARAM_MAP,1))
      WRITE(I0_STR,"(I4)") I0 

      PRINT*, TRIM(PrintIndent),"("//TRIM(ADJUSTL(I0_STR))//",*): "
      DO I1 = LBOUND(Variable%GAUSS_POINT_DOF2PARAM_MAP,2), MIN(LBOUND(Variable% &
        & GAUSS_POINT_DOF2PARAM_MAP,2)+MaxArrayLength-1, UBOUND(Variable%GAUSS_POINT_DOF2PARAM_MAP,2))
        WRITE(I1_STR,"(I4)") I1 

        PRINT*, TRIM(PrintIndent), TRIM(I1_STR), ":", Variable%GAUSS_POINT_DOF2PARAM_MAP(I0,I1)

      ENDDO  ! I1
    ENDDO  ! I0
  ELSE
    PRINT*, TRIM(PrintIndent),"INTEGER(INTG), ALLOCATABLE :: " // &
      & "GAUSS_POINT_DOF2PARAM_MAP(:,:) (not allocated)"
  ENDIF ! IF (IsAllocated)
  
  ! Variable%DATA_POINT_DOF2PARAM_MAP(:,:)
  ! case allocatable
  IsAllocated = ALLOCATED(Variable%DATA_POINT_DOF2PARAM_MAP)
  IF (IsAllocated) THEN 
    Dimension = SIZE(SHAPE(Variable%DATA_POINT_DOF2PARAM_MAP),1)
    WRITE(*,'(3A,I2,A)',advance='no') &
      & " ",TRIM(PrintIndent), "INTEGER(INTG), ALLOCATABLE :: " // &
      & "DATA_POINT_DOF2PARAM_MAP(:,:) " // &
      & "(allocated, rank=", Dimension, ", size="

    ! loop over dimensions
    DO I = 1,Dimension
      WRITE(*,'(I5)',advance='no') &
        & SIZE(Variable%DATA_POINT_DOF2PARAM_MAP, I)
      IF (I /= Dimension) WRITE(*,'(A)',advance='no') " x "
    ENDDO
    WRITE(*,'(A)',advance='no') "): "
    PRINT*," "
    DO I0 = LBOUND(Variable%DATA_POINT_DOF2PARAM_MAP,1), MIN(LBOUND(Variable% &
      & DATA_POINT_DOF2PARAM_MAP,1)+MaxArrayLength-1, UBOUND(Variable%DATA_POINT_DOF2PARAM_MAP,1))
      WRITE(I0_STR,"(I4)") I0 

      PRINT*, TRIM(PrintIndent),"("//TRIM(ADJUSTL(I0_STR))//",*): "
      DO I1 = LBOUND(Variable%DATA_POINT_DOF2PARAM_MAP,2), MIN(LBOUND(Variable% &
        & DATA_POINT_DOF2PARAM_MAP,2)+MaxArrayLength-1, UBOUND(Variable%DATA_POINT_DOF2PARAM_MAP,2))
        WRITE(I1_STR,"(I4)") I1 

        PRINT*, TRIM(PrintIndent), TRIM(I1_STR), ":", Variable%DATA_POINT_DOF2PARAM_MAP(I0,I1)

      ENDDO  ! I1
    ENDDO  ! I0
  ELSE
    PRINT*, TRIM(PrintIndent),"INTEGER(INTG), ALLOCATABLE :: " // &
      & "DATA_POINT_DOF2PARAM_MAP(:,:) (not allocated)"
  ENDIF ! IF (IsAllocated)

  
END SUBROUTINE Print_FIELD_DOF_TO_PARAM_MAP_TYPE


END MODULE PRINT_TYPES_ROUTINES_SELECTION
