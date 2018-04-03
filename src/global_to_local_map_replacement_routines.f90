!> \file
!> \author Benjamin Maier
!> \brief This module in a hack to temporary replace GLOBAL_TO_LOCAL_MAP.
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

!> This module handles all domain mappings routines.
MODULE GLOBAL_TO_LOCAL_MAP_REPLACEMENT_ROUTINES

  USE BASE_ROUTINES
  USE COMP_ENVIRONMENT
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE LISTS
  USE STRINGS
  USE TYPES

#include "macros.h"  
  
  IMPLICIT NONE
  
  LOGICAL :: GLOBAL_MAPPING_ACCESS_IN_PROGRESS = .FALSE.
  
  PUBLIC GET_GLOBAL_MAPPING, END_GLOBAL_MAPPING_ACCESS
  
  PUBLIC GET_GLOBAL_MAPPING_LOCAL_NUMBER, GET_GLOBAL_MAPPING_DOMAIN_NUMBER, GET_GLOBAL_MAPPING_LOCAL_TYPE, &
    & GET_GLOBAL_MAPPING_NUMBER_OF_DOMAINS
  
CONTAINS
  
  !>Returns the information on which subdomains the node/element with given global number is stored and their local number and type 
  FUNCTION GET_GLOBAL_MAPPING(DomainMapping,GlobalNumber,ERR,ERROR)

    !Argument variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DomainMapping !<A pointer to the domain mapping to get the global mapping from
    INTEGER(INTG), INTENT(IN) :: GlobalNumber !<The global number to get the global mapping from
    TYPE(DOMAIN_GLOBAL_MAPPING_TYPE) :: GET_GLOBAL_MAPPING   !< the returned global mapping
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    !Local Variables
    TYPE(DOMAIN_GLOBAL_MAPPING_TYPE) :: OldStoredGlobalMapping
    TYPE(DOMAIN_GLOBAL_MAPPING_TYPE) :: RetrievedGlobalMapping
    INTEGER(INTG) :: I
    INTEGER(INTG) :: ComputationalNodeNumber, NumberOfComputationalNodes
    LOGICAL :: Mismatch = .FALSE.
    
    ENTERS("GET_GLOBAL_MAPPING",ERR,ERROR,*999)
    
    ComputationalNodeNumber = COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
    NumberOfComputationalNodes = COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
    
    
    ! - check if global number is already in INCOMPLETE_GLOBAL_TO_LOCAL_MAP, if yes, done
    ! - check if global number resides on own process (all local numbers should be in INCOMPLETE_GLOBAL_TO_LOCAL_MAP) 
    
    ! - if not yet done, MPI_Irecv for all other processes (GLOBAL_MAPPING_ACCESS_IN_PROGRESS)
    ! - loop
    !   - MPI_TEST for Irecv operation from all other processes
    !     - if received data available, handle request
    !   - MPI_Ibcast to all other processes to request global number
    !   - if all required information was gathered, exit loop, set GLOBAL_MAPPING_ACCESS_IN_PROGRESS to false
    
    ! types of requests:
    ! 1) request entry (LOCAL_NUMBER, LOCAL_TYPE) if present
    ! 2) end global mapping access
    
    ! the following would be the correct access to GLOBAL_TO_LOCAL_MAP
    OldStoredGlobalMapping = DomainMapping%GLOBAL_TO_LOCAL_MAP(GlobalNumber)
    RetrievedGlobalMapping = DomainMapping%GLOBAL_TO_LOCAL_MAP(GlobalNumber)
    
    ! check if retrieved value and old stored value are the same
    IF (OldStoredGlobalMapping%NUMBER_OF_DOMAINS /= RetrievedGlobalMapping%NUMBER_OF_DOMAINS) THEN
      PRINT *, ComputationalNodeNumber," GlobalNumber: ",GlobalNumber,", NUMBER_OF_DOMAINS does not match. ",&
        & ", retrieved: ",RetrievedGlobalMapping%NUMBER_OF_DOMAINS,&
        & ", correct: ",OldStoredGlobalMapping%NUMBER_OF_DOMAINS
      Mismatch = .TRUE.
    ELSE 
      DO I=1, OldStoredGlobalMapping%NUMBER_OF_DOMAINS
        IF (OldStoredGlobalMapping%LOCAL_NUMBER(I) /= RetrievedGlobalMapping%LOCAL_NUMBER(I)) THEN
          PRINT *, ComputationalNodeNumber," GlobalNumber: ",GlobalNumber,", LOCAL_NUMBER(",I,") does not match. ",&
            & ", retrieved: ",RetrievedGlobalMapping%LOCAL_NUMBER(I),&
            & ", correct: ",OldStoredGlobalMapping%LOCAL_NUMBER(I)
          Mismatch = .TRUE.
        ENDIF
        IF (OldStoredGlobalMapping%DOMAIN_NUMBER(I) /= RetrievedGlobalMapping%DOMAIN_NUMBER(I)) THEN
          PRINT *, ComputationalNodeNumber," GlobalNumber: ",GlobalNumber,", DOMAIN_NUMBER(",I,") does not match. ",&
            & ", retrieved: ",RetrievedGlobalMapping%DOMAIN_NUMBER(I),&
            & ", correct: ",OldStoredGlobalMapping%DOMAIN_NUMBER(I)
          Mismatch = .TRUE.
        ENDIF
        IF (OldStoredGlobalMapping%LOCAL_TYPE(I) /= RetrievedGlobalMapping%LOCAL_TYPE(I)) THEN
          PRINT *, ComputationalNodeNumber," GlobalNumber: ",GlobalNumber,", LOCAL_TYPE(",I,") does not match. ",&
            & ", retrieved: ",RetrievedGlobalMapping%LOCAL_TYPE(I),&
            & ", correct: ",OldStoredGlobalMapping%LOCAL_TYPE(I)
          Mismatch = .TRUE.
        ENDIF
      ENDDO
    ENDIF
    
    IF(Mismatch) CALL FlagError("Mismatch between retrieved and stored Global mapping.",ERR,ERROR,*999)
    
    GET_GLOBAL_MAPPING = OldStoredGlobalMapping
    
    EXITS("GET_GLOBAL_MAPPING")
    RETURN
999 ERRORSEXITS("GET_GLOBAL_MAPPING_LOCAL_TYPE",ERR,ERROR)
  END FUNCTION GET_GLOBAL_MAPPING
  
  !
  !================================================================================================================================
  !

  !> Terminates access to global mapping variable, this needs to be called collectively after GET_GLOBAL_MAPPING was called for different processes, acts as a barrier
  SUBROUTINE END_GLOBAL_MAPPING_ACCESS(DOMAIN_MAPPING,ERR,ERROR,*)
    !Argument variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING !<A pointer to the domain mapping to get the global mapping from
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    ENTERS("END_GLOBAL_MAPPING_ACCESS",ERR,ERROR,*999)
    
    ! - loop
    !   - MPI_TEST for Irecv operation from all other processes
    !     - if received data available, handle request
    !     - if request is "end global mapping access", exit loop
    
    
    EXITS("END_GLOBAL_MAPPING_ACCESS")
    RETURN
999 ERRORSEXITS("GET_GLOBAL_MAPPING_LOCAL_TYPE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE END_GLOBAL_MAPPING_ACCESS
  
  !
  !================================================================================================================================
  !
  
  !>Returns the information on which subdomains the node/element with given global number is stored and their local number and type 
  FUNCTION GET_GLOBAL_MAPPING_LOCAL_NUMBER(DomainMapping,GlobalNumber,DomainIdx,ERR,ERROR)
    !Argument variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DomainMapping !<A pointer to the domain mapping to get the global mapping from
    INTEGER(INTG), INTENT(IN) :: GlobalNumber !<The global number to get the global mapping from
    INTEGER(INTG), INTENT(IN) :: DomainIdx !<The domain index of the local number array
    INTEGER(INTG) :: GET_GLOBAL_MAPPING_LOCAL_NUMBER   !< the returned local number
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(DOMAIN_GLOBAL_MAPPING_TYPE) :: GlobalMapping
    
    ENTERS("GET_GLOBAL_MAPPING_LOCAL_NUMBER",ERR,ERROR,*999)
    
    ! get the global mapping
    GlobalMapping = GET_GLOBAL_MAPPING(DomainMapping,GlobalNumber,ERR,ERROR)
    
    ! extract the local number
    GET_GLOBAL_MAPPING_LOCAL_NUMBER = GlobalMapping%LOCAL_NUMBER(DomainIdx)
    
    EXITS("GET_GLOBAL_MAPPING_LOCAL_NUMBER")
    RETURN
999 ERRORSEXITS("GET_GLOBAL_MAPPING_LOCAL_TYPE",ERR,ERROR)
  END FUNCTION
  
  !
  !================================================================================================================================
  !
  
  !>Returns the information on which subdomains the node/element with given global number is stored and their local number and type 
  FUNCTION GET_GLOBAL_MAPPING_DOMAIN_NUMBER(DomainMapping,GlobalNumber,DomainIdx,ERR,ERROR)
    !Argument variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DomainMapping !<A pointer to the domain mapping to get the global mapping from
    INTEGER(INTG), INTENT(IN) :: GlobalNumber !<The global number to get the global mapping from
    INTEGER(INTG), INTENT(IN) :: DomainIdx !<The domain index of the local number array
    INTEGER(INTG) :: GET_GLOBAL_MAPPING_DOMAIN_NUMBER   !< the returned local number
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(DOMAIN_GLOBAL_MAPPING_TYPE) :: GlobalMapping
    
    ENTERS("GET_GLOBAL_MAPPING_DOMAIN_NUMBER",ERR,ERROR,*999)
    
    ! get the global mapping
    GlobalMapping = GET_GLOBAL_MAPPING(DomainMapping,GlobalNumber,ERR,ERROR)
    
    ! extract the local number
    GET_GLOBAL_MAPPING_DOMAIN_NUMBER = GlobalMapping%DOMAIN_NUMBER(DomainIdx)
    
    EXITS("GET_GLOBAL_MAPPING_DOMAIN_NUMBER")
    RETURN
999 ERRORSEXITS("GET_GLOBAL_MAPPING_LOCAL_TYPE",ERR,ERROR)
  END FUNCTION
  
  !
  !================================================================================================================================
  !
  
  !>Returns the information on which subdomains the node/element with given global number is stored and their local number and type 
  FUNCTION GET_GLOBAL_MAPPING_LOCAL_TYPE(DomainMapping,GlobalNumber,DomainIdx,ERR,ERROR)
    !Argument variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DomainMapping !<A pointer to the domain mapping to get the global mapping from
    INTEGER(INTG), INTENT(IN) :: GlobalNumber !<The global number to get the global mapping from
    INTEGER(INTG), INTENT(IN) :: DomainIdx !<The domain index of the local type array
    INTEGER(INTG) :: GET_GLOBAL_MAPPING_LOCAL_TYPE   !< the returned local type
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(DOMAIN_GLOBAL_MAPPING_TYPE) :: GlobalMapping
    
    ENTERS("GET_GLOBAL_MAPPING_LOCAL_TYPE",ERR,ERROR,*999)
    
    ! get the global mapping
    GlobalMapping = GET_GLOBAL_MAPPING(DomainMapping,GlobalNumber,ERR,ERROR)
    
    ! extract the local type
    GET_GLOBAL_MAPPING_LOCAL_TYPE = GlobalMapping%LOCAL_TYPE(DomainIdx)
    
    EXITS("GET_GLOBAL_MAPPING_LOCAL_TYPE")
    RETURN
999 ERRORSEXITS("GET_GLOBAL_MAPPING_LOCAL_TYPE",ERR,ERROR)
  END FUNCTION
  
  !
  !================================================================================================================================
  !
  
  !>Returns the information on which subdomains the node/element with given global number is stored and their local number and type 
  FUNCTION GET_GLOBAL_MAPPING_NUMBER_OF_DOMAINS(DomainMapping,GlobalNumber,ERR,ERROR)
    !Argument variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DomainMapping !<A pointer to the domain mapping to get the global mapping from
    INTEGER(INTG), INTENT(IN) :: GlobalNumber !<The global number to get the global mapping from
    INTEGER(INTG) :: GET_GLOBAL_MAPPING_NUMBER_OF_DOMAINS   !< the returned number of domains
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(DOMAIN_GLOBAL_MAPPING_TYPE) :: GlobalMapping
    
    ENTERS("GET_GLOBAL_MAPPING_NUMBER_OF_DOMAINS",ERR,ERROR,*999)
    
    ! get the global mapping
    GlobalMapping = GET_GLOBAL_MAPPING(DomainMapping,GlobalNumber,ERR,ERROR)
    
    ! extract number of domains
    GET_GLOBAL_MAPPING_NUMBER_OF_DOMAINS = GlobalMapping%NUMBER_OF_DOMAINS
    
    EXITS("GET_GLOBAL_MAPPING_NUMBER_OF_DOMAINS")
    RETURN
999 ERRORSEXITS("GET_GLOBAL_MAPPING_LOCAL_TYPE",ERR,ERROR)
  END FUNCTION
  
END MODULE GLOBAL_TO_LOCAL_MAP_REPLACEMENT_ROUTINES