!> \file
!> \author Chris Bradley
!> \brief This module handles all mesh (node and element) routines.
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

!> This module handles all mesh (node and element) routines.
MODULE MESH_ROUTINES

  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE CmissMPI
  USE CMISS_PARMETIS
  USE ComputationRoutines
  USE ComputationAccessRoutines
  USE Constants
  USE ContextAccessRoutines
  USE COORDINATE_ROUTINES
  USE DataProjectionAccessRoutines
  USE DOMAIN_MAPPINGS
  USE Kinds
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE Lists
  USE MeshAccessRoutines
#ifndef NOMPIMOD
  USE MPI
#endif
  USE NODE_ROUTINES
  USE PRINT_TYPES_ROUTINES_SELECTION
  USE RegionAccessRoutines
  USE Strings
  USE Trees
  USE Types

#include "macros.h"  

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  PRIVATE

  !Module parameters

  !> \addtogroup MESH_ROUTINES_MeshBoundaryTypes MESH_ROUTINES::MeshBoundaryTypes
  !> \brief The types of whether or not a node/element is on a mesh domain boundary.
  !> \see MESH_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: MESH_OFF_DOMAIN_BOUNDARY=0 !<The node/element is not on the mesh domain boundary. \see MESH_ROUTINES_MeshBoundaryTypes,MESH_ROUTINES
  INTEGER(INTG), PARAMETER :: MESH_ON_DOMAIN_BOUNDARY=1 !<The node/element is on the mesh domain boundary. \see MESH_ROUTINES_MeshBoundaryTypes,MESH_ROUTINES
  !>@}
   
  !> \addtogroup MESH_ROUTINES_DecompositionTypes MESH_ROUTINES::DecompositionTypes
  !> \brief The Decomposition types parameters
  !> \see MESH_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: DECOMPOSITION_ALL_TYPE=1 !<The decomposition contains all elements. \see MESH_ROUTINES_DecompositionTypes,MESH_ROUTINES
  INTEGER(INTG), PARAMETER :: DECOMPOSITION_CALCULATED_TYPE=2 !<The element decomposition is calculated by graph partitioning. \see MESH_ROUTINES_DecompositionTypes,MESH_ROUTINES
  INTEGER(INTG), PARAMETER :: DECOMPOSITION_USER_DEFINED_TYPE=3 !<The user will set the element decomposition. \see MESH_ROUTINES_DecompositionTypes,MESH_ROUTINES
  !>@}
  
  !Module types

  !Module variables

  !Interfaces

  INTERFACE DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS
    MODULE PROCEDURE DecompositionTopology_ElementCheckExists
  END INTERFACE DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS
  
  !>Starts the process of creating a mesh
  INTERFACE MESH_CREATE_START
    MODULE PROCEDURE MESH_CREATE_START_INTERFACE
    MODULE PROCEDURE MESH_CREATE_START_REGION
  END INTERFACE !MESH_CREATE_START

  !>Initialises the meshes for a region or interface.
  INTERFACE MESHES_INITIALISE
    MODULE PROCEDURE MESHES_INITIALISE_INTERFACE
    MODULE PROCEDURE MESHES_INITIALISE_REGION
  END INTERFACE !MESHES_INITIALISE

  INTERFACE MeshTopology_ElementCheckExists
    MODULE PROCEDURE MeshTopology_ElementCheckExistsMesh
    MODULE PROCEDURE MeshTopology_ElementCheckExistsMeshElements
  END INTERFACE MeshTopology_ElementCheckExists
  
  INTERFACE MeshTopology_ElementGet
    MODULE PROCEDURE MeshTopology_ElementGetMeshElements
  END INTERFACE MeshTopology_ElementGet
  
  INTERFACE MeshTopology_NodeCheckExists
    MODULE PROCEDURE MeshTopology_NodeCheckExistsMesh
    MODULE PROCEDURE MeshTopology_NodeCheckExistsMeshNodes    
  END INTERFACE MeshTopology_NodeCheckExists  
  
  INTERFACE MeshTopology_NodeGet
    MODULE PROCEDURE MeshTopology_NodeGetMeshNodes    
  END INTERFACE MeshTopology_NodeGet
  
  PUBLIC DECOMPOSITION_ALL_TYPE,DECOMPOSITION_CALCULATED_TYPE,DECOMPOSITION_USER_DEFINED_TYPE

  PUBLIC MESH_ON_DOMAIN_BOUNDARY,MESH_OFF_DOMAIN_BOUNDARY

  PUBLIC DECOMPOSITIONS_INITIALISE,DECOMPOSITIONS_FINALISE

  PUBLIC DECOMPOSITION_CREATE_START,DECOMPOSITION_CREATE_FINISH

  PUBLIC DECOMPOSITION_DESTROY

  PUBLIC DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE

  PUBLIC DECOMPOSITION_ELEMENT_DOMAIN_GET,DECOMPOSITION_ELEMENT_DOMAIN_SET

  PUBLIC DECOMPOSITION_MESH_COMPONENT_NUMBER_GET,DECOMPOSITION_MESH_COMPONENT_NUMBER_SET
  
  PUBLIC DECOMPOSITION_NUMBER_OF_DOMAINS_GET,DECOMPOSITION_NUMBER_OF_DOMAINS_SET

  PUBLIC Decomposition_WorkGroupSet

  PUBLIC DecompositionTopology_DataPointCheckExists
  
  PUBLIC DecompositionTopology_DataProjectionCalculate

  PUBLIC DecompositionTopology_ElementCheckExists,DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS

  PUBLIC DecompositionTopology_ElementDataPointLocalNumberGet

  PUBLIC DecompositionTopology_ElementDataPointUserNumberGet

  PUBLIC DecompositionTopology_ElementGet

  PUBLIC DecompositionTopology_NumberOfElementDataPointsGet
  
  PUBLIC DECOMPOSITION_TYPE_GET,DECOMPOSITION_TYPE_SET
  
  PUBLIC DECOMPOSITION_NODE_DOMAIN_GET

  PUBLIC DECOMPOSITION_CALCULATE_LINES_SET,DECOMPOSITION_CALCULATE_FACES_SET

  PUBLIC DomainTopology_ElementBasisGet, DOMAIN_MAPPINGS_COMPARE, DOMAIN_MAPPINGS_PRINT_NUMB

  PUBLIC DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS

  PUBLIC MeshTopology_ElementCheckExists,MeshTopology_NodeCheckExists
  
  PUBLIC MESH_CREATE_START,MESH_CREATE_FINISH

  PUBLIC MESH_DESTROY
  
  PUBLIC MESH_NUMBER_OF_COMPONENTS_GET,MESH_NUMBER_OF_COMPONENTS_SET

  PUBLIC MESH_NUMBER_OF_ELEMENTS_GET,MESH_NUMBER_OF_ELEMENTS_SET

  PUBLIC MESH_TOPOLOGY_ELEMENTS_CREATE_START,MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH

  PUBLIC MESH_TOPOLOGY_ELEMENTS_DESTROY

  PUBLIC MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET,MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET

  PUBLIC MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENT_GET

  PUBLIC MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET

  PUBLIC MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET,MeshElements_ElementNodeVersionSet

  PUBLIC MeshTopology_DataPointsCalculateProjection

  PUBLIC MeshTopology_ElementOnBoundaryGet

  PUBLIC MeshElements_ElementUserNumberGet,MeshElements_ElementUserNumberSet
  
  PUBLIC MeshTopology_ElementsUserNumbersAllSet

  PUBLIC MeshTopology_NodeOnBoundaryGet

  PUBLIC MeshTopology_NodeDerivativesGet

  PUBLIC MeshTopology_NodeNumberOfDerivativesGet

  PUBLIC MeshTopology_NodeNumberOfVersionsGet

  PUBLIC MeshTopology_NodesNumberOfNodesGet

  PUBLIC MeshTopology_NodesDestroy
  
  PUBLIC MESH_SURROUNDING_ELEMENTS_CALCULATE_SET

  PUBLIC MESH_EMBEDDING_CREATE,MESH_EMBEDDING_SET_CHILD_NODE_POSITION

  PUBLIC MESH_EMBEDDING_SET_GAUSS_POINT_DATA

  PUBLIC MESHES_INITIALISE,MESHES_FINALISE

  LOGICAL, PARAMETER :: USE_OLD_GLOBAL_IMPLEMENTATION = .TRUE.     ! code that doesn't get executed when this is set to false should be removed when removal of GLOBAL_TO_LOCAL_MAP
  LOGICAL, PARAMETER :: USE_NEW_LOCAL_IMPLEMENTATION = .TRUE.     ! new code 

CONTAINS

  !
  !================================================================================================================================
  !

  !>Finialises a decomposition adjacent element information and deallocates all memory
  SUBROUTINE DECOMPOSITION_ADJACENT_ELEMENT_FINALISE(DECOMPOSITION_ADJACENT_ELEMENT,ERR,ERROR,*)
    
    !Argument variables
    TYPE(DECOMPOSITION_ADJACENT_ELEMENT_TYPE) :: DECOMPOSITION_ADJACENT_ELEMENT !<The decomposition adjacent element to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_ADJACENT_ELEMENT_FINALISE",ERR,ERROR,*999)

    DECOMPOSITION_ADJACENT_ELEMENT%NUMBER_OF_ADJACENT_ELEMENTS=0
    IF(ALLOCATED(DECOMPOSITION_ADJACENT_ELEMENT%ADJACENT_ELEMENTS)) DEALLOCATE(DECOMPOSITION_ADJACENT_ELEMENT%ADJACENT_ELEMENTS)
       
    EXITS("DECOMPOSITION_ADJACENT_ELEMENT_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_ADJACENT_ELEMENT_FINALISE",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE DECOMPOSITION_ADJACENT_ELEMENT_FINALISE

  !
  !================================================================================================================================
  !
  !>Initalises a decomposition adjacent element information.
  SUBROUTINE DECOMPOSITION_ADJACENT_ELEMENT_INITIALISE(DECOMPOSITION_ADJACENT_ELEMENT,ERR,ERROR,*)
    
    !Argument variables
    TYPE(DECOMPOSITION_ADJACENT_ELEMENT_TYPE) :: DECOMPOSITION_ADJACENT_ELEMENT !<The decomposition adjacent element to initialise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_ADJACENT_ELEMENT_INITIALISE",ERR,ERROR,*999)

    DECOMPOSITION_ADJACENT_ELEMENT%NUMBER_OF_ADJACENT_ELEMENTS=0
       
    EXITS("DECOMPOSITION_ADJACENT_ELEMENT_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_ADJACENT_ELEMENT_INITIALISE",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE DECOMPOSITION_ADJACENT_ELEMENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a domain decomposition on a given mesh. \see OPENCMISS::Iron::cmfe_DecompositionCreateFinish
  SUBROUTINE DECOMPOSITION_CREATE_FINISH(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to finish creating
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: decomposition_no
    TYPE(MESH_TYPE), POINTER :: MESH

    ENTERS("DECOMPOSITION_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      !Calculate which elements belong to which domain
      CALL DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE(DECOMPOSITION,ERR,ERROR,*999)
      !Initialise the topology information for this decomposition
      CALL DECOMPOSITION_TOPOLOGY_INITIALISE(DECOMPOSITION,ERR,ERROR,*999)
      !Initialise the domain for this computation node
      CALL DOMAIN_INITIALISE(DECOMPOSITION,ERR,ERROR,*999)
      !Calculate the decomposition topology
      CALL DECOMPOSITION_TOPOLOGY_CALCULATE(DECOMPOSITION,ERR,ERROR,*999)
      DECOMPOSITION%DECOMPOSITION_FINISHED=.TRUE.
      ! 
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      MESH=>DECOMPOSITION%MESH
      IF(ASSOCIATED(MESH)) THEN
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Mesh = ",MESH%USER_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of decompositions = ", &
          & MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS,ERR,ERROR,*999)
        DO decomposition_no=1,MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Decomposition number = ",decomposition_no,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global number        = ", &
            & MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_no)%PTR%GLOBAL_NUMBER,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    User number          = ", &
            & MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_no)%PTR%USER_NUMBER,ERR,ERROR,*999)
        ENDDO !decomposition_no
      ELSE
        CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
      ENDIF
    ENDIF
    
    EXITS("DECOMPOSITION_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_CREATE_FINISH",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Starts the creation of a domain decomposition for a given mesh. \see OPENCMISS::Iron::cmfe_DecompositionCreateStart
  SUBROUTINE DECOMPOSITION_CREATE_START(USER_NUMBER,MESH,DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the decomposition 
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to decompose
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<On return a pointer to the created decomposition. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: decompositionIdx,dummyErr
    TYPE(ComputationEnvironmentType), POINTER :: computationEnvironment
    TYPE(ContextType), POINTER :: context
    TYPE(DECOMPOSITION_TYPE), POINTER :: newDecomposition
    TYPE(DECOMPOSITION_PTR_TYPE), ALLOCATABLE :: newDecompositions(:)
    TYPE(REGION_TYPE), POINTER :: region    
    TYPE(WorkGroupType), POINTER :: worldWorkGroup
    TYPE(LIST_TYPE), POINTER :: GlobalElementDomain 

    NULLIFY(newDecomposition)

    ENTERS("DECOMPOSITION_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(mesh)) THEN
      IF(mesh%MESH_FINISHED) THEN
        IF(ASSOCIATED(mesh%topology)) THEN
          IF(ASSOCIATED(mesh%decompositions)) THEN
            IF(ASSOCIATED(decomposition)) THEN
              CALL FlagError("Decomposition is already associated.",err,error,*999)
            ELSE
              NULLIFY(decomposition)
              CALL Decomposition_UserNumberFind(USER_NUMBER,MESH,decomposition,ERR,ERROR,*999)
              IF(ASSOCIATED(DECOMPOSITION)) THEN
                LOCAL_ERROR="Decomposition number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
                  & " has already been created on mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))//"."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              ELSE
                NULLIFY(newDecomposition)
                CALL Decomposition_Initialise(newDecomposition,err,error,*999)
                !Set default decomposition properties              
                newDecomposition%GLOBAL_NUMBER=MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS+1
                newDecomposition%USER_NUMBER=USER_NUMBER
                newDecomposition%DECOMPOSITIONS=>MESH%DECOMPOSITIONS
                newDecomposition%MESH=>MESH
                newDecomposition%region=>mesh%region
                newDecomposition%INTERFACE=>mesh%INTERFACE
                newDecomposition%numberOfDimensions=mesh%NUMBER_OF_DIMENSIONS
                newDecomposition%numberOfComponents=mesh%NUMBER_OF_COMPONENTS
                !By default, the process of decompostion was done on the first mesh components. But the decomposition is the
                !same for all mesh components, since the decomposition is element-based.
                newDecomposition%MESH_COMPONENT_NUMBER=1
                !Default decomposition is all the mesh with one domain.
                newDecomposition%DECOMPOSITION_TYPE=DECOMPOSITION_ALL_TYPE
                NULLIFY(region)
                CALL Mesh_RegionGet(mesh,region,err,error,*999)
                NULLIFY(context)
                CALL Region_ContextGet(region,context,err,error,*999)
                NULLIFY(computationEnvironment)
                CALL Context_ComputationEnvironmentGet(context,computationEnvironment,err,error,*999)
                NULLIFY(worldWorkGroup)
                CALL ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err,error,*999)
                newDecomposition%workGroup=>worldWorkGroup
                newDecomposition%NUMBER_OF_DOMAINS=1
                newDecomposition%numberOfElements=mesh%NUMBER_OF_ELEMENTS
                ALLOCATE(newDecomposition%ELEMENT_DOMAIN(MESH%NUMBER_OF_ELEMENTS),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate new decomposition element domain.",ERR,ERROR,*999)
                newDecomposition%ELEMENT_DOMAIN=0  

               ! create DECOMPOSITION%GlobalElementDomain, with standard size (of 10 entries)
                NULLIFY(GlobalElementDomain)
                CALL LIST_CREATE_START(GlobalElementDomain,ERR,ERROR,*999)
                CALL LIST_DATA_TYPE_SET(GlobalElementDomain,LIST_INTG_TYPE,ERR,ERROR,*999)
                CALL LIST_DATA_DIMENSION_SET(GlobalElementDomain,2,ERR,ERROR,*999)
                CALL LIST_CREATE_FINISH(GlobalElementDomain,ERR,ERROR,*999)
                newDecomposition%GlobalElementDomain=>GlobalElementDomain
               
                NULLIFY(newDecomposition%ELEMENTS_MAPPING) 
                ! Nullify the domain
                NULLIFY(newDecomposition%DOMAIN)
                ! Nullify the topology
                NULLIFY(newDecomposition%TOPOLOGY)

                !Add new decomposition into list of decompositions on the mesh
                ALLOCATE(newDecompositions(mesh%decompositions%NUMBER_OF_DECOMPOSITIONS+1),STAT=err)
                IF(ERR/=0) CALL FlagError("Could not allocate new decompositions.",ERR,ERROR,*999)
                DO decompositionIdx=1,MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS
                  newDecompositions(decompositionIdx)%ptr=>mesh%decompositions%decompositions(decompositionIdx)%ptr
                ENDDO !decompositionIdx
                newDecompositions(mesh%decompositions%NUMBER_OF_DECOMPOSITIONS+1)%ptr=>newDecomposition
                CALL MOVE_ALLOC(newDecompositions,mesh%decompositions%decompositions)
                mesh%decompositions%NUMBER_OF_DECOMPOSITIONS=mesh%decompositions%NUMBER_OF_DECOMPOSITIONS+1        
                decomposition=>newDecomposition

              ENDIF
            ENDIF
          ELSE
            LOCAL_ERROR="The decompositions on mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))// &
              & " are not associated."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Mesh topology is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Mesh has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_CREATE_START")
    RETURN
999 IF(ASSOCIATED(newDecomposition)) CALL Decomposition_Finalise(newDecomposition,dummyErr,dummyError,*998)
998 IF(ALLOCATED(newDecompositions)) DEALLOCATE(newDecompositions)
    ERRORSEXITS("DECOMPOSITION_CREATE_START",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE DECOMPOSITION_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroys a domain decomposition identified by a user number and deallocates all memory. \see OPENCMISS::Iron::cmfe_DecompositionDestroy
  SUBROUTINE DECOMPOSITION_DESTROY_NUMBER(USER_NUMBER,MESH,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the decomposition to destroy.
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh containing the decomposition.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: decomposition_idx,decomposition_position
    LOGICAL :: FOUND    
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(DECOMPOSITION_PTR_TYPE), ALLOCATABLE :: newDecompositions(:)

    ENTERS("DECOMPOSITION_DESTROY_NUMBER",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(ASSOCIATED(MESH%DECOMPOSITIONS)) THEN

        !Find the decomposition identified by the user number
        FOUND=.FALSE.
        decomposition_position=0
        DO WHILE(decomposition_position<MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS.AND..NOT.FOUND)
          decomposition_position=decomposition_position+1
          IF(MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_position)%PTR%USER_NUMBER==USER_NUMBER) FOUND=.TRUE.
        ENDDO
        
        IF(FOUND) THEN
          
          decomposition=>MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_position)%PTR

          !Finalise the decomposition
          CALL Decomposition_Finalise(decomposition,err,error,*999)

          !Remove the decomposition from the list of decompositions
          IF(MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS>1) THEN
            ALLOCATE(newDecompositions(MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS-1),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate new decompositions.",ERR,ERROR,*999)
            DO decomposition_idx=1,MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS
              IF(decomposition_idx<decomposition_position) THEN
                newDecompositions(decomposition_idx)%ptr=>mesh%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR
              ELSE IF(decomposition_idx>decomposition_position) THEN
                MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR%GLOBAL_NUMBER= &
                  & MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR%GLOBAL_NUMBER-1
                newDecompositions(decomposition_idx-1)%PTR=>MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR
              ENDIF
            ENDDO !decomposition_idx
            CALL MOVE_ALLOC(newDecompositions,mesh%decompositions%decompositions)
            MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS=MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS-1
          ELSE
            DEALLOCATE(MESH%DECOMPOSITIONS%DECOMPOSITIONS)
            MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS=0
          ENDIF
          
        ELSE
          LOCAL_ERROR="Decomposition number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
            & " has not been created on mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The decompositions on mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))// &
          & " are not associated."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_DESTROY_NUMBER")
    RETURN
999 IF(ALLOCATED(newDecompositions)) DEALLOCATE(newDecompositions)
    ERRORSEXITS("DECOMPOSITION_DESTROY_NUMBER",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE DECOMPOSITION_DESTROY_NUMBER

  !
  !================================================================================================================================
  !

  !>Destroys a domain decomposition identified by an object and deallocates all memory. \see OPENCMISS::Iron::cmfe_DecompositionDestroy
  SUBROUTINE DECOMPOSITION_DESTROY(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITON !<The decomposition to destroy.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: decomposition_idx,decomposition_position
    LOGICAL :: FOUND    
    TYPE(MESH_TYPE), POINTER :: MESH !
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(DECOMPOSITION_PTR_TYPE), ALLOCATABLE :: newDecompositions(:)
    TYPE(DECOMPOSITIONS_TYPE), POINTER :: decompositions

    ENTERS("DECOMPOSITION_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      DECOMPOSITIONS=>DECOMPOSITION%DECOMPOSITIONS
      IF(ASSOCIATED(DECOMPOSITIONS)) THEN
        mesh=>decomposition%mesh
        IF(ASSOCIATED(mesh)) THEN
          !Find the decomposition identified by the user number
          FOUND=.FALSE.
          decomposition_position=0
          DO WHILE(decomposition_position<MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS.AND..NOT.FOUND)
            decomposition_position=decomposition_position+1
            IF(MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_position)%PTR%USER_NUMBER==DECOMPOSITION%USER_NUMBER) FOUND=.TRUE.
          ENDDO

          !Finalise the decomposition
          CALL Decomposition_Finalise(decomposition,err,error,*999)
          
          !Remove the decomposition from the list of decompositions
          IF(DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS>1) THEN
            ALLOCATE(newDecompositions(DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS-1),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate new decompositions.",ERR,ERROR,*999)
            DO decomposition_idx=1,DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS
              IF(decomposition_idx<decomposition_position) THEN
                newDecompositions(decomposition_idx)%PTR=>DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR
              ELSE IF(decomposition_idx>decomposition_position) THEN
                DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR%GLOBAL_NUMBER= &
                  & DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR%GLOBAL_NUMBER-1
                newDecompositions(decomposition_idx-1)%PTR=>DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR
              ENDIF
            ENDDO !decomposition_idx
            CALL MOVE_ALLOC(newDecompositions,DECOMPOSITIONS%DECOMPOSITIONS)
            DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS=DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS-1
          ELSE
            DEALLOCATE(DECOMPOSITIONS%DECOMPOSITIONS)
            DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS=0
          ENDIF
        ELSE
          CALL FlagError("Decomposition mesh is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition decompositions is not associated.",ERR,ERROR,*999)
       ENDIF
    ELSE
      CALL FlagError("Decompositions is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_DESTROY")
    RETURN
999 IF(ALLOCATED(newDecompositions)) DEALLOCATE(newDecompositions)
    ERRORSEXITS("DECOMPOSITION_DESTROY",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_DESTROY

  !
  !================================================================================================================================
  !

  !>Calculates the element domains for a decomposition of a mesh. \see OPENCMISS::Iron::cmfe_DecompositionElementDomainCalculate
  SUBROUTINE DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to calculate the element domains for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    !INTEGER(INTG) :: number_elem_indicies,elem_index,elem_count,ne,nn,my_computational_node_number,number_computational_nodes, &
    !  & no_computational_node,ELEMENT_START,ELEMENT_STOP,MY_ELEMENT_START,MY_ELEMENT_STOP,NUMBER_OF_ELEMENTS, &
    !  & MY_NUMBER_OF_ELEMENTS,MPI_IERROR,MAX_NUMBER_ELEMENTS_PER_NODE,component_idx,minNumberXi
    !INTEGER(INTG), ALLOCATABLE :: ELEMENT_COUNT(:),ELEMENT_PTR(:),ELEMENT_INDICIES(:),ELEMENT_DISTANCE(:),DISPLACEMENTS(:), &
    !  & RECEIVE_COUNTS(:)
    INTEGER(INTG) :: NumberElementIndices,elem_index,elem_count,ne,nn,ComputationalNodeNo,MyComputationalNodeNumber,&
      & NumberComputationalNodes,NumberOfElements,NumberElementsOnOwnComputationalNode, &
      & ElementStart,ElementStop,MyElementStart,MyElementStop, I,J,ElementToSendNo,GlobalElementNo, &
      & MyNumberOfElements,MPI_IERROR,MaxNumberElementsPerNode,component_idx,MinNumberXi,RequestHandle,&
      & MaximumNumberElementsToSendToOneComputationalNode, GlobalElementNoDomainPair(2)
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_COUNT(:),ADJACENT_NODE_PTR(:),ADJACENT_NODE_INDICES(:),ELEMENT_DISTRIBUTION(:), &
      ! & DISPLACEMENTS(:),RECEIVE_COUNTS(:),ElementDomain(:),NumberElementsDomainTemporary(:),RequestHandles(:),&
      ! & NumberElementsToReceive(:),ReceiveBuffer(:)
      & RECEIVE_COUNTS(:),ElementDomain(:) 

    INTEGER(INTG) :: ELEMENT_WEIGHT(1),WEIGHT_FLAG,NUMBER_FLAG,NUMBER_OF_CONSTRAINTS, &
      & NUMBER_OF_COMMON_NODES,PARMETIS_OPTIONS(0:2)
    !ParMETIS now has double for these
    !REAL(SP) :: UBVEC(1)
    !REAL(SP), ALLOCATABLE :: TPWGTS(:)

    ! LOGICAL, ALLOCATABLE :: IRecvCompleted(:)
    TYPE(LIST_TYPE), POINTER :: GlobalElementDomain 

    INTEGER(INTG) :: groupCommunicator

    REAL(DP) :: UBVEC(1)
    REAL(DP), ALLOCATABLE :: TPWGTS(:)
    !REAL(DP) :: NUMBER_ELEMENTS_PER_NODE
    REAL(DP) :: NumberElementsPerNode 
    TYPE(ContextType), POINTER :: context    
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(ASSOCIATED(DECOMPOSITION%MESH)) THEN
        MESH=>DECOMPOSITION%MESH
        IF(ASSOCIATED(MESH%TOPOLOGY)) THEN
          IF(ASSOCIATED(DECOMPOSITION%GlobalElementDomain)) THEN 
            !component_idx=DECOMPOSITION%MESH_COMPONENT_NUMBER
          
            ! get rank count and own rank 
            !NumberComputationalNodes=ComputationalEnvironment_NumberOfNodesGet(ERR,ERROR)
            !IF(ERR/=0) GOTO 999
            !MyComputationalNodeNumber=ComputationalEnvironment_NodeNumberGet(ERR,ERROR)
            !IF(ERR/=0) GOTO 999
          
            !SELECT CASE(DECOMPOSITION%DECOMPOSITION_TYPE)
            !CASE(DECOMPOSITION_ALL_TYPE)
              !Do nothing. Decomposition checked below.
            !CASE(DECOMPOSITION_CALCULATED_TYPE)
            !Calculate the general decomposition

            !IF(DECOMPOSITION%NUMBER_OF_DOMAINS==1) THEN
            component_idx=DECOMPOSITION%MESH_COMPONENT_NUMBER

            CALL WorkGroup_GroupCommunicatorGet(decomposition%workGroup,groupCommunicator,err,error,*999)
            CALL WorkGroup_NumberOfGroupNodesGet(decomposition%workGroup,numberOfComputationNodes,err,error,*999)
            CALL WorkGroup_GroupNodeNumberGet(decomposition%workGroup,myComputationNodeNumber,err,error,*999)

            ! Replaced by work groups above!!   
            ! get rank count and own rank
            !NumberComputationalNodes=ComputationalEnvironment_NumberOfNodesGet(ERR,ERROR)
            !IF(ERR/=0) GOTO 999
            !MyComputationalNodeNumber=ComputationalEnvironment_NodeNumberGet(ERR,ERROR)
            !IF(ERR/=0) GOTO 999

            SELECT CASE(DECOMPOSITION%DECOMPOSITION_TYPE)
            CASE(DECOMPOSITION_ALL_TYPE)
              !Do nothing. Decomposition checked below.
            CASE(DECOMPOSITION_CALCULATED_TYPE)
              !Calculate the general decomposition

              IF(DECOMPOSITION%NUMBER_OF_DOMAINS==1) THEN 
                IF (USE_OLD_GLOBAL_IMPLEMENTATION) THEN 
                  ! set domain of every element to 0
                  DECOMPOSITION%ELEMENT_DOMAIN=0
                ENDIF


                IF (USE_NEW_LOCAL_IMPLEMENTATION) THEN
                  ! get total number of elements from mesh
                  MyNumberOfElements = MESH%NUMBER_OF_ELEMENTS 

                  NULLIFY(GlobalElementDomain)
                  CALL LIST_CREATE_START(GlobalElementDomain,ERR,ERROR,*999)
                  CALL LIST_DATA_TYPE_SET(GlobalElementDomain,LIST_INTG_TYPE,ERR,ERROR,*999)
                  CALL LIST_DATA_DIMENSION_SET(GlobalElementDomain,2,ERR,ERROR,*999)
                  CALL LIST_INITIAL_SIZE_SET(GlobalElementDomain,MyNumberOfElements,ERR,ERROR,*999)
                  CALL LIST_CREATE_FINISH(GlobalElementDomain,ERR,ERROR,*999)
                  
                  ! fill DECOMPOSITION%GlobalElementDomain list with data from ElementDomain
                  ! the list contains pairs of (global el. no., domain no)
                  DO GlobalElementNo=0,MyNumberOfElements
                    GlobalElementNoDomainPair(1) = GlobalElementNo
                    GlobalElementNoDomainPair(2) = MyComputationalNodeNumber    ! this is 0
                    
                    CALL LIST_ITEM_ADD(GlobalElementDomain,GlobalElementNoDomainPair,ERR,ERROR,*999)
                  ENDDO
                  
                  ! create new list with known number of elements
                  IF (ASSOCIATED(DECOMPOSITION%GlobalElementDomain)) THEN
                    CALL LIST_DESTROY(DECOMPOSITION%GlobalElementDomain,ERR,ERROR,*999)
                  ENDIF
                  DECOMPOSITION%GlobalElementDomain=>GlobalElementDomain                   

                ENDIF
                
              ELSE
                ! Compute average number of elements per node
                NumberElementsPerNode=REAL(MESH%NUMBER_OF_ELEMENTS,DP)/REAL(NumberComputationalNodes,DP)
                ElementStart=1    ! first index of element range for a particular rank
                ElementStop=0     ! last index of element range for a particular rank
                MaxNumberElementsPerNode=-1 
                ! Allocate arrays
                ! note for merge: DISPLACEMENTS(i) == ELEMENT_DISTRIBUTION(i), therefore variable DISPLACEMENTS removed
                !ALLOCATE(DISPLACEMENTS(0:NumberComputationalNodes-1),STAT=ERR)
                !IF(ERR/=0) CALL FlagError("Could not allocate displacements.",ERR,ERROR,*999) 


                ! ELEMENT_DISTRIBUTION(np): first global element number for process np
                ALLOCATE(ELEMENT_DISTRIBUTION(0:NumberComputationalNodes),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate element distance.",ERR,ERROR,*999) 

                ! RECEIVE_COUNTS(np): number of global elements for process np
                ALLOCATE(RECEIVE_COUNTS(0:NumberComputationalNodes-1),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate recieve counts.",ERR,ERROR,*999) 

                ELEMENT_DISTRIBUTION(0)=0 

                ! Loop over computational nodes
                DO ComputationalNodeNo=0,NumberComputationalNodes-1 
                  ! Compute range [ElementStart, ElementStop] of elements that will be on ComputationalNodeNo
                  ElementStart=ElementStop+1
                  IF(ComputationalNodeNo==NumberComputationalNodes-1) THEN
                    ElementStop=MESH%NUMBER_OF_ELEMENTS
                  ELSE
                    ElementStop=ElementStart+NINT(NumberElementsPerNode,INTG)-1
                  ENDIF
                  
                  ! If after current range there are more computational nodes left than elements, reduce current range
                  IF((NumberComputationalNodes-1-ComputationalNodeNo)>(MESH%NUMBER_OF_ELEMENTS-ElementStop)) &
                    & ElementStop=MESH%NUMBER_OF_ELEMENTS-(NumberComputationalNodes-1-ComputationalNodeNo)
                  IF(ElementStart>MESH%NUMBER_OF_ELEMENTS) ElementStart=MESH%NUMBER_OF_ELEMENTS
                  IF(ElementStop>MESH%NUMBER_OF_ELEMENTS) ElementStop=MESH%NUMBER_OF_ELEMENTS
                  
                  ! assign beginning of element range for next process (C numbering)
                  ELEMENT_DISTRIBUTION(ComputationalNodeNo+1)=ElementStop !C numbering
                  
                  ! assign number of elements in range for current rank to RECEIVE_COUNTS
                  NumberOfElements=ElementStop-ElementStart+1
                  RECEIVE_COUNTS(ComputationalNodeNo)=NumberOfElements
                  
                  ! store maximum value of NumberOfElements on a rank
                  IF(NumberOfElements>MaxNumberElementsPerNode) MaxNumberElementsPerNode=NumberOfElements
                  
                  ! on local rank
                  IF(ComputationalNodeNo==MyComputationalNodeNumber) THEN
                    MyElementStart=ElementStart
                    MyElementStop=ElementStop
                    MyNumberOfElements=NumberOfElements
                    NumberElementIndices=0
                    
                    ! count number of nodes that are assigned to this rank
                    ! loop over local elements
                    DO ne=MyElementStart,MyElementStop
                      BASIS=>MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%ELEMENTS(ne)%BASIS
                      NumberElementIndices=NumberElementIndices+BASIS%NUMBER_OF_NODES
                    ENDDO !ne
                  ENDIF
                ENDDO !ComputationalNodeNo
                
                ALLOCATE(ADJACENT_NODE_PTR(0:MyNumberOfElements),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate element pointer list.",ERR,ERROR,*999)
                ALLOCATE(ADJACENT_NODE_INDICES(0:NumberElementIndices-1),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate element indicies list.",ERR,ERROR,*999)
                ALLOCATE(TPWGTS(1:DECOMPOSITION%NUMBER_OF_DOMAINS),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate tpwgts.",ERR,ERROR,*999)
                elem_index=0
                elem_count=0
                ADJACENT_NODE_PTR(0)=0    ! ADJACENT_NODE_PTR(localElementNo) is the first index in ADJACENT_NODE_INDICES of the element (C numbering)
                MinNumberXi=99999
                
                ! loop over local elements
                DO ne=MyElementStart,MyElementStop
                  elem_count=elem_count+1
                  BASIS=>MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%ELEMENTS(ne)%BASIS
                  
                  ! Store minimum number of xi-directions of the elements basis
                  IF(BASIS%NUMBER_OF_XI<MinNumberXi) MinNumberXi=BASIS%NUMBER_OF_XI
                  
                  ! loop over nodes of element
                  DO nn=1,BASIS%NUMBER_OF_NODES
                    ADJACENT_NODE_INDICES(elem_index)=MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%ELEMENTS(ne)% &
                      & MESH_ELEMENT_NODES(nn)-1 !C numbering
                      
                    ! ADJACENT_NODE_INDICES(ADJACENT_NODE_PTR(elem_count-1)+1:ADJACENT_NODE_PTR(elem_count)+1) contains
                    ! the indices of adjacent nodes that define the graph to be partioned
                
                    elem_index=elem_index+1
                  ENDDO !nn
                  ADJACENT_NODE_PTR(elem_count)=elem_index !C numbering
                ENDDO !ne
                
                !Set up ParMETIS variables
                NULLIFY(context)
                CALL WorkGroup_ContextGet(decomposition%workGroup,context,err,error,*999)
                CALL Context_RandomSeedsSizeGet(context,randomSeedsSize,err,error,*999)
                ALLOCATE(randomSeeds(randomSeedsSize),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate random seeds.",err,error,*999)
                CALL Context_RandomSeedsGet(context,randomSeeds,err,error,*999)              
                WEIGHT_FLAG=0 !No weights
                ELEMENT_WEIGHT(1)=1 !Isn't used due to weight flag
                NUMBER_FLAG=0 !C Numbering as there is a bug with Fortran numbering
                NUMBER_OF_CONSTRAINTS=1   ! number of weight per vertex in the graph that will be partitioned
                
                ! NUMBER_OF_COMMON_NODES is the number of nodes that two elements have to share to be considered adjacent
                IF(MinNumberXi==1) THEN     ! for one-dimensional meshes it is 1 node
                  NUMBER_OF_COMMON_NODES=1    
                ELSE
                  NUMBER_OF_COMMON_NODES=2  ! for 2D or 3D nodes it is 2 nodes, i.e. in a 3D mesh also elements sharing an edge are considered adjacent 
                ENDIF

                TPWGTS=1.0_DP/REAL(DECOMPOSITION%NUMBER_OF_DOMAINS,DP)  ! target value for weights fraction per partition
                UBVEC=1.05_DP ! imbalance tolerance: 5% 

                PARMETIS_OPTIONS(0)=1 !If zero, defaults are used, otherwise next two values are used
                PARMETIS_OPTIONS(1)=7 !Level of information to output
                !PARMETIS_OPTIONS(2)=cmissRandomSeeds(1) !Seed for random number generator
                PARMETIS_OPTIONS(2)=randomSeeds(1) !Seed for random number generator

                IF(ALLOCATED(randomSeeds)) DEALLOCATE(randomSeeds)

                IF (USE_OLD_GLOBAL_IMPLEMENTATION) THEN 
                ! use variable DECOMPOSITION%ELEMENT_DOMAIN, distribute information such that every rank has all information
                
                !Call ParMETIS to calculate the partitioning of the mesh graph. The graph is the dual graph to the node graph.
                !see parmetis reference manual: http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/manual.pdf
                  CALL PARMETIS_PARTMESHKWAY(&
                    & ELEMENT_DISTRIBUTION,&    ! starting indices of the elements per process
                    & ADJACENT_NODE_PTR,ADJACENT_NODE_INDICES,&  ! local adjacent nodes in CSR format: ADJACENT_NODE_PTR is the starting indices to sequences in ADJACENT_NODE_INDICES that contain global node numbers that are adjacent to each of the local element
                    & ELEMENT_WEIGHT,WEIGHT_FLAG, &  ! WEIGHT_FLAG=0, i.e. no weighting of the graph
                    & NUMBER_FLAG, &      ! NUMBER_FLAG=0 means c-style numbering
                    & NUMBER_OF_CONSTRAINTS,NUMBER_OF_COMMON_NODES,DECOMPOSITION%NUMBER_OF_DOMAINS, & ! number of contraints=1, NUMBER_OF_COMMON_NODES: see above, number of generated partitions
                    & TPWGTS,UBVEC,PARMETIS_OPTIONS, &    ! UBVEC: imbalance tolerance
                    & DECOMPOSITION%NUMBER_OF_EDGES_CUT, &  ! output: number of edge cuts
                    & DECOMPOSITION%ELEMENT_DOMAIN(ELEMENT_DISTRIBUTION(MyComputationalNodeNumber)+1:), &    ! output: ELEMENT_DOMAIN(ne) contains newly determined domain number of local element ne
                    & groupCommunicator,ERR,ERROR,*999)
                
                !Transfer all the element domain information to the other computational nodes so that each rank has all the info
                  IF(NumberComputationalNodes>1) THEN
                  !This should work on a single processor but doesn't for mpich2 under windows. Maybe a bug? Avoid for now.
                  CALL MPI_ALLGATHERV(MPI_IN_PLACE,MaxNumberElementsPerNode,MPI_INTEGER,DECOMPOSITION%ELEMENT_DOMAIN, &
                    & RECEIVE_COUNTS,ELEMENT_DISTRIBUTION,MPI_INTEGER,groupCommunicator,MPI_IERROR)

                  CALL MPI_ERROR_CHECK("MPI_ALLGATHERV",MPI_IERROR,ERR,ERROR,*999)
                ! note: allgatherV: complete fractioned vector on each rank, where each rank has only its own part
                  ENDIF

                ENDIF   ! if USE_OLD_GLOBAL_IMPLEMENTATION
              
                IF (USE_NEW_LOCAL_IMPLEMENTATION) THEN               

                  ! use variable DECOMPOSITION%GLOBAL_ELEMENT_NUMBER which is the same as DOMAIN%LOCAL_TO_GLOBAL_MAP for elements
                  ALLOCATE(ElementDomain(0:MyNumberOfElements-1),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate ElementDomain variable.",ERR,ERROR,*999)
                  
                  !Call ParMETIS to calculate the partitioning of the mesh graph. The graph is the dual graph to the node graph.
                  !see parmetis reference manual: http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/manual.pdf
                  ! Store result into local variable ElementDomain
                  CALL PARMETIS_PARTMESHKWAY(&
                    & ELEMENT_DISTRIBUTION,&    ! starting indices of the elements per process
                    & ADJACENT_NODE_PTR,ADJACENT_NODE_INDICES,&  ! local adjacent nodes in CSR format: ADJACENT_NODE_PTR is the starting indices to sequences in ADJACENT_NODE_INDICES that contain global node numbers that are adjacent to each of the local element
                    & ELEMENT_WEIGHT,WEIGHT_FLAG, &  ! WEIGHT_FLAG=0, i.e. no weighting of the graph
                    & NUMBER_FLAG, &      ! NUMBER_FLAG=0 means c-style numbering
                    & NUMBER_OF_CONSTRAINTS,NUMBER_OF_COMMON_NODES,DECOMPOSITION%NUMBER_OF_DOMAINS, & ! number of contraints=1, NUMBER_OF_COMMON_NODES: see above, number of generated partitions
                    & TPWGTS,UBVEC,PARMETIS_OPTIONS, &    ! UBVEC: imbalance tolerance
                    & DECOMPOSITION%NUMBER_OF_EDGES_CUT, &  ! output: number of edge cuts
                    & ElementDomain, &    ! output: ELEMENT_DOMAIN(ne) contains newly determined domain number of local element ne, c style numbering
                    & groupCommunicator,ERR,ERROR,*999)
                  
                  ! fill DECOMPOSITION%GlobalElementDomain list with data from ElementDomain
                  ! the list contains pairs of (global el. no., domain no)
                  
                  NULLIFY(GlobalElementDomain)
                  CALL LIST_CREATE_START(GlobalElementDomain,ERR,ERROR,*999)
                  CALL LIST_DATA_TYPE_SET(GlobalElementDomain,LIST_INTG_TYPE,ERR,ERROR,*999)
                  CALL LIST_DATA_DIMENSION_SET(GlobalElementDomain,2,ERR,ERROR,*999)
                  CALL LIST_INITIAL_SIZE_SET(GlobalElementDomain,MyNumberOfElements,ERR,ERROR,*999)
                  CALL LIST_CREATE_FINISH(GlobalElementDomain,ERR,ERROR,*999)                 
              
                  DO GlobalElementNo=MyElementStart,MyElementStop
                    GlobalElementNoDomainPair(1) = GlobalElementNo
                    GlobalElementNoDomainPair(2) = MyComputationalNodeNumber
                    
                    CALL LIST_ITEM_ADD(GlobalElementDomain,GlobalElementNoDomainPair,ERR,ERROR,*999)
                  ENDDO
                  
                  ! create new list with known number of elements
                  IF (ASSOCIATED(DECOMPOSITION%GlobalElementDomain)) THEN
                    CALL LIST_DESTROY(DECOMPOSITION%GlobalElementDomain,ERR,ERROR,*999)
                  ENDIF
                  DECOMPOSITION%GlobalElementDomain=>GlobalElementDomain
                  
                  DEALLOCATE(ElementDomain) 
                ENDIF

              ! deallocate temporary arrays
                DEALLOCATE(RECEIVE_COUNTS)
!               DEALLOCATE(ELEMENT_DISTANCE)
!               DEALLOCATE(ELEMENT_PTR)
!               DEALLOCATE(ELEMENT_INDICIES)
                DEALLOCATE(ELEMENT_DISTRIBUTION)
                DEALLOCATE(ADJACENT_NODE_PTR)
                DEALLOCATE(ADJACENT_NODE_INDICES)
                DEALLOCATE(TPWGTS)

              ENDIF
              
            CASE(DECOMPOSITION_USER_DEFINED_TYPE)
            !Do nothing. Decomposition checked below.          
            CASE DEFAULT
              CALL FlagError("Invalid domain decomposition type.",ERR,ERROR,*999)            
            END SELECT

          !Check decomposition and check that each domain has an element in it.
!          ALLOCATE(ELEMENT_COUNT(0:number_computational_nodes-1),STAT=ERR)
!          IF(ERR/=0) CALL FlagError("Could not allocate element count.",ERR,ERROR,*999)
!          ELEMENT_COUNT=0
!          DO elem_index=1,MESH%NUMBER_OF_ELEMENTS
!            no_computational_node=DECOMPOSITION%ELEMENT_DOMAIN(elem_index)
!            IF(no_computational_node>=0.AND.no_computational_node<number_computational_nodes) THEN
!              ELEMENT_COUNT(no_computational_node)=ELEMENT_COUNT(no_computational_node)+1
!            ELSE
!              LOCAL_ERROR="The computational node number of "//TRIM(NUMBER_TO_VSTRING(no_computational_node,"*",ERR,ERROR))// &
!                & " for element number "//TRIM(NUMBER_TO_VSTRING(elem_index,"*",ERR,ERROR))// &
!                & " is invalid. The computational node number must be between 0 and "// &
!                & TRIM(NUMBER_TO_VSTRING(number_computational_nodes-1,"*",ERR,ERROR))//"."
!              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
!            ENDIF
!          ENDDO !elem_index
!          DO no_computational_node=0,number_computational_nodes-1
!            IF(ELEMENT_COUNT(no_computational_node)==0) THEN
!              LOCAL_ERROR="Invalid decomposition. There are no elements in computational node "// &
!                & TRIM(NUMBER_TO_VSTRING(no_computational_node,"*",ERR,ERROR))//"."
!              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
!            ENDIF
!          ENDDO !no_computational_node
!          DEALLOCATE(ELEMENT_COUNT)
            IF (USE_NEW_LOCAL_IMPLEMENTATION) THEN
          ! distribute information in GlobalElementDomain such that afterwards each rank has only its local information in DECOMPOSITION%ELEMENTS_MAPPING
              CALL DECOMPOSITION_ELEMENTS_MAPPING_MAP_CALCULATE(DECOMPOSITION,ERR,ERROR,*999) 
              CALL DECOMPOSITION_ELEMENTS_MAPPING_CALCULATE(DECOMPOSITION,ERR,ERROR,*999)

              CALL DECOMPOSITION_ADJACENT_DOMAINS_CALCULATE(DECOMPOSITION,ERR,ERROR,*999) 
            END IF
 
            IF (USE_OLD_GLOBAL_IMPLEMENTATION) THEN 
              !Check decomposition and check that each domain has an element in it.
              ALLOCATE(ELEMENT_COUNT(0:NumberComputationalNodes-1),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate element count.",ERR,ERROR,*999)
              ELEMENT_COUNT=0
              ! count number of elements per rank
              ! loop over global elements
              DO elem_index=1,MESH%NUMBER_OF_ELEMENTS
                ComputationalNodeNo=DECOMPOSITION%ELEMENT_DOMAIN(elem_index)
                IF(ComputationalNodeNo>=0.AND.ComputationalNodeNo<NumberComputationalNodes) THEN
                  ELEMENT_COUNT(ComputationalNodeNo)=ELEMENT_COUNT(ComputationalNodeNo)+1
                ELSE
                  LOCAL_ERROR="The computational node number of "//TRIM(NUMBER_TO_VSTRING(ComputationalNodeNo,"*",ERR,ERROR))// &
                    & " for element number "//TRIM(NUMBER_TO_VSTRING(elem_index,"*",ERR,ERROR))// &
                    & " is invalid. The computational node number must be between 0 and "// &
                    & TRIM(NUMBER_TO_VSTRING(NumberComputationalNodes-1,"*",ERR,ERROR))//"."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDDO !elem_index
            
              ! Check that each rank has at least 1 element
              DO ComputationalNodeNo=0,NumberComputationalNodes-1
                IF(ELEMENT_COUNT(ComputationalNodeNo)==0) THEN
                  LOCAL_ERROR="Invalid decomposition. There are no elements in computational node "// &
                    & TRIM(NUMBER_TO_VSTRING(ComputationalNodeNo,"*",ERR,ERROR))//"."
                  CALL FlagWarning(LOCAL_ERROR,ERR,ERROR,*999) 
                ENDIF
              END DO !ComputationalNodeNo
              DEALLOCATE(ELEMENT_COUNT)
          
            ENDIF  ! IF USE_OLD_GLOBAL_IMPLEMENTATION          
          ELSE
            CALL FlagError("Decomposition GlobalElementDomain is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Decomposition mesh topology is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
      ENDIF

    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Decomposition for mesh number ",DECOMPOSITION%MESH%USER_NUMBER, &
        & ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of domains = ", DECOMPOSITION%NUMBER_OF_DOMAINS, &
        & ERR,ERROR,*999)
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Element domains:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Decomposition type = ", DECOMPOSITION%DECOMPOSITION_TYPE, &
        & ERR,ERROR,*999)
      IF(DECOMPOSITION%DECOMPOSITION_TYPE==DECOMPOSITION_CALCULATED_TYPE) THEN
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of edges cut = ",DECOMPOSITION%NUMBER_OF_EDGES_CUT, &
          & ERR,ERROR,*999)
      ENDIF
    
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of elements = ",DECOMPOSITION%numberOfElements, &
        & ERR,ERROR,*999)

      DO ne=1,DECOMPOSITION%numberOfElements
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Element = ",ne,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Domain = ",DECOMPOSITION%ELEMENT_DOMAIN(ne), &
          & ERR,ERROR,*999)
      ENDDO !ne
    
    ENDIF
    
    EXITS("DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE")
    RETURN
999 IF(ALLOCATED(RECEIVE_COUNTS)) DEALLOCATE(RECEIVE_COUNTS)
!   IF(ALLOCATED(DISPLACEMENTS)) DEALLOCATE(DISPLACEMENTS)
!   IF(ALLOCATED(ELEMENT_DISTANCE)) DEALLOCATE(ELEMENT_DISTANCE)
!   IF(ALLOCATED(ELEMENT_PTR)) DEALLOCATE(ELEMENT_PTR)
!   IF(ALLOCATED(ELEMENT_INDICIES)) DEALLOCATE(ELEMENT_INDICIES)
    IF(ALLOCATED(ELEMENT_DISTRIBUTION)) DEALLOCATE(ELEMENT_DISTRIBUTION)
    IF(ALLOCATED(ADJACENT_NODE_PTR)) DEALLOCATE(ADJACENT_NODE_PTR)
    IF(ALLOCATED(ADJACENT_NODE_INDICES)) DEALLOCATE(ADJACENT_NODE_INDICES) 
    IF(ALLOCATED(TPWGTS)) DEALLOCATE(TPWGTS)
    IF(ALLOCATED(randomSeeds)) DEALLOCATE(randomSeeds)
    ERRORSEXITS("DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE
  
  !
  !================================================================================================================================
  !

!>Calculates for decomposition the ELEMENTS_MAPPING%NUMBER_OF_LOCAL, ELEMENTS_MAPPING%NUMBER_OF_DOMAIN_LOCAL and ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP from GlobalElementDomain
  SUBROUTINE DECOMPOSITION_ELEMENTS_MAPPING_MAP_CALCULATE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to calculate the element domains for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: MyComputationalNodeNumber,NumberComputationalNodes,MPI_IERROR,MpiMemoryWindow,ComputationalNodeNo,&
      & ListSize, I, J, GlobalElementNo, Domain,GlobalElementNoDomainPair(2), NumberComputationalNodesToSendTo, &
      & MaximumNumberElementsToSendToOneComputationalNode, NumberComputationalNodesToReceiveFrom, &
      & MaximumNumberElementsToReceiveFromOneComputationalNode, ComputationalNodesToSendToNo, ComputationalNodesToReceiveFromNo, &
      & NumberElementsToComputationalNode, NumberElementsFromComputationalNode, NumberElementsInReceivedElementsList, &
      & ElementToSendNo, NumberOwnLocalElements
    INTEGER(INTG) :: groupCommunicator
    INTEGER(INTG), ALLOCATABLE :: NumberLocalElementsOnRank(:), NumberStoredElementsForRank(:), ElementSendBuffers(:,:), &
      & SendRequestHandle(:), ReceiveRequestHandle(:)
    INTEGER(INTG), ALLOCATABLE, DIMENSION(:,:) :: ElementReceiveBuffers
    INTEGER(KIND=MPI_ADDRESS_KIND) LB, EXTENT, MpiWindowSizeBytes, MpiTargetDisplacement, MpiDisplacementsUnitBytes
    TYPE(LIST_TYPE), POINTER :: ReceivedElementsList

    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(MeshComponentTopologyType), POINTER :: TOPOLOGY
   ! LOGICAL, PARAMETER :: DEBUGGING = .FALSE.

    ENTERS("DECOMPOSITION_ELEMENTS_MAPPING_MAP",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION%ELEMENTS_MAPPING)) THEN
      CALL FlagError("Decomposition Elements mapping is already allocated.",ERR,ERROR,*999)
    ENDIF

    CALL WorkGroup_GroupCommunicatorGet(decomposition%workGroup,groupCommunicator,err,error,*999)
    CALL WorkGroup_NumberOfGroupNodesGet(decomposition%workGroup,numberOfComputationNodes,err,error,*999)
    CALL WorkGroup_GroupNodeNumberGet(decomposition%workGroup,myComputationNodeNumber,err,error,*999)

    ! Replaced by work groups above!!   
    ! get rank count and own rank
    ! NumberComputationalNodes=ComputationalEnvironment_NumberOfNodesGet(ERR,ERROR)
    ! IF(ERR/=0) GOTO 999
    ! MyComputationalNodeNumber=ComputationalEnvironment_NodeNumberGet(ERR,ERROR)
    ! IF(ERR/=0) GOTO 999

    ! Allocate arrays
    ALLOCATE(DECOMPOSITION%ELEMENTS_MAPPING,STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate ELEMENTS_MAPPING.",ERR,ERROR,*999)
     
    ALLOCATE(DECOMPOSITION%ELEMENTS_MAPPING%NUMBER_OF_DOMAIN_LOCAL(0:NumberComputationalNodes-1),STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate NUMBER_OF_DOMAIN_LOCAL.",ERR,ERROR,*999)
     
    CALL DOMAIN_MAPPINGS_MAPPING_INITIALISE(DECOMPOSITION%ELEMENTS_MAPPING, DECOMPOSITION%workGroup,ERR,ERROR,*999)
    ELEMENTS_MAPPING=>DECOMPOSITION%ELEMENTS_MAPPING

    
    IF (NumberComputationalNodes == 1) THEN    
      ! count number of elements that we have from each rank
      CALL List_NumberOfItemsGet(DECOMPOSITION%GlobalElementDomain,NumberOwnLocalElements,ERR,ERROR,*999)
      ELEMENTS_MAPPING%NUMBER_OF_LOCAL = NumberOwnLocalElements
      ELEMENTS_MAPPING%NUMBER_OF_DOMAIN_LOCAL(0) = ELEMENTS_MAPPING%NUMBER_OF_LOCAL 
 
      IF (NumberOwnLocalElements == 0) CALL FlagError("No elements set for decomposition.",ERR,ERROR,*999) 

     ! create list for elements
      NULLIFY(ReceivedElementsList)
      CALL LIST_CREATE_START(ReceivedElementsList,ERR,ERROR,*999)
      CALL LIST_DATA_TYPE_SET(ReceivedElementsList,LIST_INTG_TYPE,ERR,ERROR,*999)
      CALL LIST_INITIAL_SIZE_SET(ReceivedElementsList,ELEMENTS_MAPPING%NUMBER_OF_LOCAL,ERR,ERROR,*999) 
      CALL LIST_CREATE_FINISH(ReceivedElementsList,ERR,ERROR,*999)
    
      ! add elements that are already stored on the own rank to list
      DO I=1,NumberOwnLocalElements
        CALL List_ItemGet(DECOMPOSITION%GlobalElementDomain,I,GlobalElementNoDomainPair,ERR,ERROR,*999)
        GlobalElementNo = GlobalElementNoDomainPair(1)
        Domain = GlobalElementNoDomainPair(2)
        CALL LIST_ITEM_ADD(ReceivedElementsList,GlobalElementNo,ERR,ERROR,*999)
      ENDDO
      
      ! Sort list by global element number and store in local elements storage
      CALL List_Sort(ReceivedElementsList,ERR,ERROR,*999)
      CALL LIST_DETACH_AND_DESTROY(ReceivedElementsList,NumberElementsInReceivedElementsList,&
        & ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP,ERR,ERROR,*999) 

    ELSE    ! multiple computational nodes
      ! figure out which rank has how many of my local elements
      ALLOCATE(NumberLocalElementsOnRank(0:NumberComputationalNodes-1),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate NumberLocalElementsOnRank variable.",ERR,ERROR,*999)
      NumberLocalElementsOnRank = 0
            
      ALLOCATE(NumberStoredElementsForRank(0:NumberComputationalNodes-1),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate NumberStoredElementsForRank variable.",ERR,ERROR,*999)
      NumberStoredElementsForRank = 0
      ! count number of elements that we have from each rank
      CALL List_NumberOfItemsGet(DECOMPOSITION%GlobalElementDomain,ListSize,ERR,ERROR,*999)
      
      DO I=1,ListSize
        CALL List_ItemGet(DECOMPOSITION%GlobalElementDomain,I,GlobalElementNoDomainPair,ERR,ERROR,*999)
        GlobalElementNo = GlobalElementNoDomainPair(1)
        Domain = GlobalElementNoDomainPair(2) 

        IF (DIAGNOSTICS2) THEN
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"GlobalElementDomain("//TRIM(NumberToVString(I,"*",ERR,ERROR))//")=("//&
            & TRIM(NumberToVString(GlobalElementNo,"*",ERR,ERROR))//","//TRIM(NumberToVString(Domain,"*",ERR,ERROR))//")",&
            & ERR,ERROR,*999)
          !WRITE(*,"(I2,A,I3,A,I3,A,I3,A)") MyComputationalNodeNumber,": GlobalElementDomain(",I, &
          !  & ")=(", GlobalElementNo,",",Domain,")" 
        ENDIF

        NumberStoredElementsForRank(Domain) = NumberStoredElementsForRank(Domain) + 1
      END DO
      
      IF (DIAGNOSTICS2) THEN 
        DO ComputationalNodeNo=0,NumberComputationalNodes-1

          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"ComputationalNodeNo: "// &
            & TRIM(NumberToVString(ComputationalNodeNo,"*",ERR,ERROR))//", NumberStoredElementsForRank: "// &
            & TRIM(NumberToVString(NumberStoredElementsForRank(ComputationalNodeNo),"*",ERR,ERROR)),ERR,ERROR,*999)
          !WRITE(*,"(I2,A,I2,A,I3)") MyComputationalNodeNumber,": ComputationalNodeNo: ",ComputationalNodeNo, &
          !  & ", NumberStoredElementsForRank: ", NumberStoredElementsForRank(ComputationalNodeNo) 
        ENDDO 

      ENDIF

      ! get size in bytes of MPI_INT
      CALL MPI_TYPE_GET_EXTENT(MPI_INT, LB, EXTENT, MPI_IERROR)
      CALL MPI_ERROR_CHECK("MPI_TYPE_GET_EXTENT",MPI_IERROR,ERR,ERROR,*999)
      
      ! start mpi RMA
      MpiWindowSizeBytes = NumberComputationalNodes*EXTENT
      MpiDisplacementsUnitBytes = EXTENT
      CALL MPI_WIN_CREATE(NumberLocalElementsOnRank(0:NumberComputationalNodes-1), MpiWindowSizeBytes, &
        & MpiDisplacementsUnitBytes, MPI_INFO_NULL, groupCommunicator, MpiMemoryWindow, MPI_IERROR)
      CALL MPI_ERROR_CHECK("MPI_WIN_CREATE",MPI_IERROR,ERR,ERROR,*999)
      
      CALL MPI_WIN_FENCE(0, MpiMemoryWindow, MPI_IERROR)
      CALL MPI_ERROR_CHECK("MPI_WIN_FENCE",MPI_IERROR,ERR,ERROR,*999)
            
      ! write information of elements to remote ranks
      DO ComputationalNodeNo=0,NumberComputationalNodes-1
        IF (NumberStoredElementsForRank(ComputationalNodeNo) /= 0) THEN
          MpiTargetDisplacement = MyComputationalNodeNumber
          CALL MPI_PUT(NumberStoredElementsForRank(ComputationalNodeNo), 1, MPI_INT, ComputationalNodeNo,&
            & MpiTargetDisplacement, 1, MPI_INT, MpiMemoryWindow, MPI_IERROR)
          CALL MPI_ERROR_CHECK("MPI_WIN_FENCE",MPI_IERROR,ERR,ERROR,*999)
        ENDIF
      ENDDO
            
      ! synchronize RMA
      CALL MPI_WIN_FENCE(0, MpiMemoryWindow, MPI_IERROR)
      CALL MPI_ERROR_CHECK("MPI_WIN_FENCE",MPI_IERROR,ERR,ERROR,*999)
            
      ! deallocate window 
      CALL MPI_WIN_FREE(MpiMemoryWindow, MPI_IERROR)
      CALL MPI_ERROR_CHECK("MPI_WIN_FREE",MPI_IERROR,ERR,ERROR,*999)
            
      IF (DIAGNOSTICS2) THEN 
        DO ComputationalNodeNo=0,NumberComputationalNodes-1

          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"ComputationalNodeNo: "// &
            & TRIM(NumberToVString(ComputationalNodeNo,"*",ERR,ERROR))//", NumberLocalElementsOnRank: "// &
            & TRIM(NumberToVString(NumberLocalElementsOnRank(ComputationalNodeNo),"*",ERR,ERROR)),ERR,ERROR,*999)
          !WRITE(*,"(I2,A,I2,A,I3)") MyComputationalNodeNumber,": ComputationalNodeNo: ",ComputationalNodeNo, &
          !  & ", NumberLocalElementsOnRank: ", NumberLocalElementsOnRank(ComputationalNodeNo) 
        ENDDO 

      ENDIF

      ! count number of computational nodes to send elements to
      NumberComputationalNodesToSendTo = 0
      MaximumNumberElementsToSendToOneComputationalNode = 0
      DO ComputationalNodeNo=0,NumberComputationalNodes-1
        IF (ComputationalNodeNo == MyComputationalNodeNumber) CYCLE 
   
        IF (NumberStoredElementsForRank(ComputationalNodeNo) /= 0) THEN
          NumberComputationalNodesToSendTo = NumberComputationalNodesToSendTo + 1
          MaximumNumberElementsToSendToOneComputationalNode = MAX(MaximumNumberElementsToSendToOneComputationalNode, &
            & NumberStoredElementsForRank(ComputationalNodeNo))
        ENDIF
      ENDDO 
      ! count number of computational nodes to receive elements from
      NumberComputationalNodesToReceiveFrom = 0
      MaximumNumberElementsToReceiveFromOneComputationalNode = 0
      DO ComputationalNodeNo=0,NumberComputationalNodes-1
        IF (ComputationalNodeNo == MyComputationalNodeNumber) CYCLE
      
        IF (NumberLocalElementsOnRank(ComputationalNodeNo) /= 0) THEN
          NumberComputationalNodesToReceiveFrom = NumberComputationalNodesToReceiveFrom + 1
          MaximumNumberElementsToReceiveFromOneComputationalNode = MAX(MaximumNumberElementsToReceiveFromOneComputationalNode, &
            & NumberLocalElementsOnRank(ComputationalNodeNo))
        ENDIF
      ENDDO
      
      ! get actual elements by message passing
      !transfer information on actual elements to the ranks that own them.
      
      ! allocate send buffers
      ALLOCATE(ElementSendBuffers(1:MaximumNumberElementsToSendToOneComputationalNode,1:NumberComputationalNodesToSendTo), &
        & STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate ElementSendBuffers array with size "//&
        & TRIM(NUMBER_TO_VSTRING(MaximumNumberElementsToSendToOneComputationalNode,"*",ERR,ERROR))//" x "//&
        & TRIM(NUMBER_TO_VSTRING(NumberComputationalNodesToSendTo,"*",ERR,ERROR))//".",ERR,ERROR,*999)
      
      ! allocate receive buffers
      ALLOCATE(ElementReceiveBuffers(1:MaximumNumberElementsToReceiveFromOneComputationalNode, &
        & 1:NumberComputationalNodesToReceiveFrom), STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate ElementReceiveBuffers array with size "//&
        & TRIM(NUMBER_TO_VSTRING(MaximumNumberElementsToReceiveFromOneComputationalNode,"*",ERR,ERROR))//" x "//&
        & TRIM(NUMBER_TO_VSTRING(NumberComputationalNodesToReceiveFrom,"*",ERR,ERROR))//".",ERR,ERROR,*999)
      
      ! allocate request handles
      ALLOCATE(SendRequestHandle(1:NumberComputationalNodesToSendTo), STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate SendRequestHandle array with size "//&
        & TRIM(NUMBER_TO_VSTRING(NumberComputationalNodesToSendTo,"*",ERR,ERROR))//".",ERR,ERROR,*999)
      
      ALLOCATE(ReceiveRequestHandle(1:NumberComputationalNodesToReceiveFrom), STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate ReceiveRequestHandle array with size "//&
        & TRIM(NUMBER_TO_VSTRING(NumberComputationalNodesToReceiveFrom,"*",ERR,ERROR))//".",ERR,ERROR,*999)
      
      ! Commit all send commands
      ComputationalNodesToSendToNo = 1
      
      ! Loop over all computational nodes, except the own, starting with the next to the own computational node
      ComputationalNodeNo = MyComputationalNodeNumber
      DO I=1,NumberComputationalNodes-1
        ComputationalNodeNo=ComputationalNodeNo+1
        IF (ComputationalNodeNo == NumberComputationalNodes) ComputationalNodeNo=0         

        
        ! find out how many elements this computational node will get from MyComputationalNodeNumber
        NumberElementsToComputationalNode = NumberStoredElementsForRank(ComputationalNodeNo) 

        ! send elements to computational node ComputationalNodeNo
        IF (NumberElementsToComputationalNode /= 0) THEN 

          ! fill send buffer 
          ElementToSendNo = 1
          DO J=1,ListSize
            CALL List_ItemGet(DECOMPOSITION%GlobalElementDomain,J,GlobalElementNoDomainPair,ERR,ERROR,*999)
            GlobalElementNo = GlobalElementNoDomainPair(1)
            Domain = GlobalElementNoDomainPair(2)
          
            IF (Domain == ComputationalNodeNo) THEN
              ElementSendBuffers(ElementToSendNo,ComputationalNodesToSendToNo) = GlobalElementNo
              ElementToSendNo = ElementToSendNo + 1
            ENDIF
          ENDDO
          
          CALL MPI_ISEND(ElementSendBuffers(1:NumberElementsToComputationalNode,ComputationalNodesToSendToNo), &
            & NumberElementsToComputationalNode, MPI_INT, ComputationalNodeNo, 0, &
            & groupCommunicator, SendRequestHandle(ComputationalNodesToSendToNo), MPI_IERROR)
          CALL MPI_ERROR_CHECK("MPI_ISEND",MPI_IERROR,ERR,ERROR,*999)
          
          ComputationalNodesToSendToNo = ComputationalNodesToSendToNo+1
        ENDIF
      ENDDO        

      ! commit all receive commands
    
      ! Loop over all computational nodes backwards, except the own, starting with the previous to the own computational node
      ComputationalNodeNo = MyComputationalNodeNumber
      ComputationalNodesToReceiveFromNo = 1
      DO I=1,NumberComputationalNodes-1
        ComputationalNodeNo=ComputationalNodeNo-1
        IF (ComputationalNodeNo == -1) ComputationalNodeNo=NumberComputationalNodes-1
        
        ! find out how many elements this computational node will get from ComputationalNodeNo
        NumberElementsFromComputationalNode = NumberLocalElementsOnRank(ComputationalNodeNo)
        
        IF (NumberElementsFromComputationalNode /= 0) THEN
          
          CALL MPI_IRECV(ElementReceiveBuffers(1:NumberElementsFromComputationalNode,ComputationalNodesToReceiveFromNo), &
            & NumberElementsFromComputationalNode, MPI_INT, ComputationalNodeNo, 0, &
            & groupCommunicator, ReceiveRequestHandle(ComputationalNodesToReceiveFromNo), MPI_IERROR)
          CALL MPI_ERROR_CHECK("MPI_IRECV",MPI_IERROR,ERR,ERROR,*999)
          
          ComputationalNodesToReceiveFromNo = ComputationalNodesToReceiveFromNo + 1
        ENDIF
      ENDDO 
      

      ! create list for received elements
      NULLIFY(ReceivedElementsList)
      CALL LIST_CREATE_START(ReceivedElementsList,ERR,ERROR,*999)
      CALL LIST_DATA_TYPE_SET(ReceivedElementsList,LIST_INTG_TYPE,ERR,ERROR,*999)
      CALL LIST_INITIAL_SIZE_SET(ReceivedElementsList,ListSize,ERR,ERROR,*999)
      CALL LIST_CREATE_FINISH(ReceivedElementsList,ERR,ERROR,*999)
    
      ! add elements that are already stored on the own rank to list
      DO I=1,ListSize
        CALL List_ItemGet(DECOMPOSITION%GlobalElementDomain,I,GlobalElementNoDomainPair,ERR,ERROR,*999)
        GlobalElementNo = GlobalElementNoDomainPair(1)
        Domain = GlobalElementNoDomainPair(2)      
        IF (Domain == MyComputationalNodeNumber) THEN
          CALL LIST_ITEM_ADD(ReceivedElementsList,GlobalElementNo,ERR,ERROR,*999)
        ENDIF
      ENDDO
    
      ! wait for all receive requests to terminate and store received values in receivebuffer list
      ComputationalNodeNo = MyComputationalNodeNumber
      ComputationalNodesToReceiveFromNo = 1
      DO I=1,NumberComputationalNodes-1
        ComputationalNodeNo=ComputationalNodeNo-1
        IF (ComputationalNodeNo == -1) ComputationalNodeNo=NumberComputationalNodes-1

        CALL MPI_WAIT(ReceiveRequestHandle(ComputationalNodesToReceiveFromNo), MPI_STATUS_IGNORE, MPI_IERROR)
        CALL MPI_ERROR_CHECK("MPI_WAIT",MPI_IERROR,ERR,ERROR,*999)
        
        NumberElementsFromComputationalNode = NumberLocalElementsOnRank(ComputationalNodeNo)
        IF (NumberElementsFromComputationalNode /= 0) THEN
          !DO I=1,NumberElementsFromComputationalNode
          !  CALL LIST_ITEM_ADD(ReceivedElementsList,ElementReceiveBuffers(I,ComputationalNodesToReceiveFromNo),ERR,ERROR,*999)
          !ENDDO
          CALL LIST_ITEMS_ADD(ReceivedElementsList,ElementReceiveBuffers(1:NumberElementsFromComputationalNode,&
            & ComputationalNodesToReceiveFromNo),ERR,ERROR,*999)
            
          IF (DIAGNOSTICS2) THEN
            
            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"index "// &
              & TRIM(NumberToVString(ComputationalNodesToReceiveFromNo,"*",ERR,ERROR))//", "// &
              & TRIM(NumberToVString(NumberElementsFromComputationalNode,"*",ERR,ERROR))//" received element(s) from "//&
              & TRIM(NumberToVString(ComputationalNodeNo,"*",ERR,ERROR))//": ",ERR,ERROR,*999)
            !WRITE(*,"(I2,3(A,I2),A)") MyComputationalNodeNumber,": index ",ComputationalNodesToReceiveFromNo,&
            !  & ", ",NumberElementsFromComputationalNode," received element(s) from ", ComputationalNodeNo, ": " 
            DO J=1,NumberElementsFromComputationalNode

              CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,TRIM(NumberToVString(ElementReceiveBuffers(J, &
                & ComputationalNodesToReceiveFromNo),"(I10)",ERR,ERROR)),ERR,ERROR,*999)
              !WRITE(*,"(I10)") ElementReceiveBuffers(J,ComputationalNodesToReceiveFromNo) 
            ENDDO
          ENDIF
          
          ComputationalNodesToReceiveFromNo = ComputationalNodesToReceiveFromNo + 1
        ENDIF
      ENDDO         
      ! Sort list by global element number and store in local elements storage
      CALL List_Sort(ReceivedElementsList,ERR,ERROR,*999)
      CALL LIST_REMOVE_DUPLICATES(ReceivedElementsList,ERR,ERROR,*999)
      CALL LIST_DETACH_AND_DESTROY(ReceivedElementsList,NumberElementsInReceivedElementsList,&
        & ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP,ERR,ERROR,*999) 
        
      ! Eventually shrink array
      ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP = RESHAPE(ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP, &
        & [NumberElementsInReceivedElementsList], [0]) 

      ELEMENTS_MAPPING%NUMBER_OF_LOCAL = NumberElementsInReceivedElementsList
      ELEMENTS_MAPPING%NUMBER_OF_DOMAIN_LOCAL(MyComputationalNodeNumber) = NumberElementsInReceivedElementsList 
        
      ! Terminate all send commands
      DO ComputationalNodesToSendToNo=1,NumberComputationalNodesToSendTo
        CALL MPI_WAIT(SendRequestHandle(ComputationalNodesToSendToNo), MPI_STATUS_IGNORE, MPI_IERROR)
        CALL MPI_ERROR_CHECK("MPI_WAIT",MPI_IERROR,ERR,ERROR,*999)
      ENDDO
      ! fill NUMBER_OF_DOMAIN_LOCAL
      CALL MPI_ALLGATHER(MPI_IN_PLACE,1,MPI_INTEGER,ELEMENTS_MAPPING%& 
        & NUMBER_OF_DOMAIN_LOCAL(0:NumberComputationalNodes-1), &
        & 1,MPI_INTEGER,groupCommunicator,MPI_IERROR)
      CALL MPI_ERROR_CHECK("MPI_ALLGATHER",MPI_IERROR,ERR,ERROR,*999)
         

      !IF (DIAGNOSTICS2) THEN
      !  WRITE(*,"(I2,A)") MyComputationalNodeNumber,": NUMBER_OF_DOMAIN_LOCAL after ALLGATHER: "
      !  DO I=0,NumberComputationalNodes-1
      !    PRINT *, ELEMENTS_MAPPING%NUMBER_OF_DOMAIN_LOCAL(I)
      !  ENDDO
      !ENDIF  
    ENDIF
    

    ! compute total number of elements
    ELEMENTS_MAPPING%NUMBER_OF_GLOBAL = 0
    DO ComputationalNodeNo=0,NumberComputationalNodes-1
      ELEMENTS_MAPPING%NUMBER_OF_GLOBAL = ELEMENTS_MAPPING%NUMBER_OF_GLOBAL &
        & + ELEMENTS_MAPPING%NUMBER_OF_DOMAIN_LOCAL(ComputationalNodeNo) 
    ENDDO
    
    IF (ELEMENTS_MAPPING%NUMBER_OF_GLOBAL /= DECOMPOSITION%MESH%NUMBER_OF_ELEMENTS) THEN
      CALL FlagError("The mesh contains "//&
        & TRIM(NUMBER_TO_VSTRING(DECOMPOSITION%MESH%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//" elements, but you only specified "//&
        & TRIM(NUMBER_TO_VSTRING(ELEMENTS_MAPPING%NUMBER_OF_GLOBAL,"*",ERR,ERROR))//&
        & " elements in the decomposition.",ERR,ERROR,*999)
    ENDIF
             
    EXITS("DECOMPOSITION_ELEMENTS_MAPPING_MAP")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_ELEMENTS_MAPPING_MAP", ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_ELEMENTS_MAPPING_MAP_CALCULATE

  !
  !================================================================================================================================
  !



  !>Calculates for decomposition the internal, boundary and ghost numberings in ELEMENTS_MAPPING
  SUBROUTINE DECOMPOSITION_ELEMENTS_MAPPING_CALCULATE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to calculate the element domains for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ElementLocalNo, ElementGlobalNo, NumberComputationalNodes, MyComputationalNodeNumber, nn, NodeGlobalNo, &
      & AdjacentElementIdx, AdjacentElementGlobalNo, I, MPI_IERROR
    INTEGER(INTG) :: groupCommunicator
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(MeshComponentTopologyType), POINTER :: TOPOLOGY
    TYPE(LIST_TYPE), POINTER :: GhostElementsList, InternalElementsList, BoundaryElementsList
    LOGICAL :: ElementIsInterior
    INTEGER(INTG), ALLOCATABLE :: InternalElements(:), GhostElements(:), BoundaryElements(:)

    ENTERS("DECOMPOSITION_ELEMENTS_MAPPING_CALCULATE",ERR,ERROR,*999)

    CALL WorkGroup_GroupCommunicatorGet(decomposition%workGroup,groupCommunicator,err,error,*999)
    CALL WorkGroup_NumberOfGroupNodesGet(decomposition%workGroup,numberOfComputationNodes,err,error,*999)
    CALL WorkGroup_GroupNodeNumberGet(decomposition%workGroup,myComputationNodeNumber,err,error,*999)

    ! Replaced by work groups above!!   
    ! get rank count and own rank
    !NumberComputationalNodes=ComputationalEnvironment_NumberOfNodesGet(ERR,ERROR)
    !IF(ERR/=0) GOTO 999
    !MyComputationalNodeNumber=ComputationalEnvironment_NodeNumberGet(ERR,ERROR)
    !IF(ERR/=0) GOTO 999
    
    ELEMENTS_MAPPING=>DECOMPOSITION%ELEMENTS_MAPPING
    TOPOLOGY=>DECOMPOSITION%MESH%TOPOLOGY(DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
      
    ! set boundary types in domain mapping
    
    ALLOCATE(GhostElementsList,STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate ghost elements list.",ERR,ERROR,*999)
         
    ! create list to hold ghost elements (global numbers)
    NULLIFY(GhostElementsList)
    CALL LIST_CREATE_START(GhostElementsList,ERR,ERROR,*999)
    CALL LIST_DATA_DIMENSION_SET(GhostElementsList,1,ERR,ERROR,*999)
    CALL LIST_DATA_TYPE_SET(GhostElementsList,LIST_INTG_TYPE,ERR,ERROR,*999)
    CALL LIST_INITIAL_SIZE_SET(GhostElementsList,CEILING(ELEMENTS_MAPPING%NUMBER_OF_LOCAL*0.2),ERR,ERROR,*999)
    CALL LIST_CREATE_FINISH(GhostElementsList,ERR,ERROR,*999)
    
    ! create list to hold internal elements (local numbers)
    NULLIFY(InternalElementsList)
    CALL LIST_CREATE_START(InternalElementsList,ERR,ERROR,*999)
    CALL LIST_DATA_DIMENSION_SET(InternalElementsList,1,ERR,ERROR,*999)
    CALL LIST_DATA_TYPE_SET(InternalElementsList,LIST_INTG_TYPE,ERR,ERROR,*999)
    CALL LIST_INITIAL_SIZE_SET(InternalElementsList,ELEMENTS_MAPPING%NUMBER_OF_LOCAL,ERR,ERROR,*999)
    CALL LIST_CREATE_FINISH(InternalElementsList,ERR,ERROR,*999)
    
    ! create list to hold boundary elements (local numbers)
    NULLIFY(BoundaryElementsList)
    CALL LIST_CREATE_START(BoundaryElementsList,ERR,ERROR,*999)
    CALL LIST_DATA_DIMENSION_SET(BoundaryElementsList,1,ERR,ERROR,*999)
    CALL LIST_DATA_TYPE_SET(BoundaryElementsList,LIST_INTG_TYPE,ERR,ERROR,*999)
    CALL LIST_INITIAL_SIZE_SET(BoundaryElementsList,CEILING(ELEMENTS_MAPPING%NUMBER_OF_LOCAL*0.2),ERR,ERROR,*999)
    CALL LIST_CREATE_FINISH(BoundaryElementsList,ERR,ERROR,*999)
    
    ! loop over elements on local rank
    DO ElementLocalNo = 1,ELEMENTS_MAPPING%NUMBER_OF_LOCAL
      ElementGlobalNo = ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(ElementLocalNo)
      
      ElementIsInterior = .TRUE.
      ! loop over nodes of element
      DO nn=1,TOPOLOGY%ELEMENTS%ELEMENTS(ElementGlobalNo)%BASIS%NUMBER_OF_NODES
        NodeGlobalNo=TOPOLOGY%ELEMENTS%ELEMENTS(ElementGlobalNo)%MESH_ELEMENT_NODES(nn)
        
        ! loop over adjacent elements of node
        DO AdjacentElementIdx=1,TOPOLOGY%NODES%NODES(NodeGlobalNo)%numberOfSurroundingElements
          AdjacentElementGlobalNo=TOPOLOGY%NODES%NODES(NodeGlobalNo)%surroundingElements(AdjacentElementIdx)
          
          ! if adjacent element is not local
          IF (.NOT. DECOMPOSITION_ELEMENT_DOMAIN_IS_LOCAL(DECOMPOSITION,AdjacentElementGlobalNo,ERR,ERROR)) THEN
          
            CALL LIST_ITEM_ADD(GhostElementsList,AdjacentElementGlobalNo,ERR,ERROR,*999)
            ElementIsInterior = .FALSE.
          ENDIF
        ENDDO  ! AdjacentElementIdx
      ENDDO   ! nn
      
      IF (ElementIsInterior) THEN
        CALL LIST_ITEM_ADD(InternalElementsList,ElementLocalNo,ERR,ERROR,*999)
      ELSE
        CALL LIST_ITEM_ADD(BoundaryElementsList,ElementLocalNo,ERR,ERROR,*999)
      ENDIF
    ENDDO
    
    CALL LIST_SORT(InternalElementsList,ERR,ERROR,*999)
    CALL LIST_DETACH_AND_DESTROY(InternalElementsList,ELEMENTS_MAPPING%NUMBER_OF_INTERNAL,InternalElements,ERR,ERROR,*999)
    CALL LIST_SORT(BoundaryElementsList,ERR,ERROR,*999)
    CALL LIST_DETACH_AND_DESTROY(BoundaryElementsList,ELEMENTS_MAPPING%NUMBER_OF_BOUNDARY,BoundaryElements,ERR,ERROR,*999)
    
    CALL List_Sort(GhostElementsList,ERR,ERROR,*999)
    CALL LIST_REMOVE_DUPLICATES(GhostElementsList,ERR,ERROR,*999)
    CALL LIST_DETACH_AND_DESTROY(GhostElementsList,ELEMENTS_MAPPING%NUMBER_OF_GHOST,GhostElements,ERR,ERROR,*999)
    
    ! compute total number of local elements and reshape LOCAL_TO_GLOBAL_MAP
    ELEMENTS_MAPPING%TOTAL_NUMBER_OF_LOCAL = ELEMENTS_MAPPING%NUMBER_OF_LOCAL + ELEMENTS_MAPPING%NUMBER_OF_GHOST
        
    ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP = RESHAPE(ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP, &
      & [ELEMENTS_MAPPING%TOTAL_NUMBER_OF_LOCAL], [0])
    
    ! add global numbers of GhostElements at the end of LOCAL_TO_GLOBAL_MAP
    IF (ELEMENTS_MAPPING%NUMBER_OF_GHOST > 0) THEN
      ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(ELEMENTS_MAPPING%NUMBER_OF_LOCAL+1:ELEMENTS_MAPPING%TOTAL_NUMBER_OF_LOCAL)=&
        & GhostElements(1:ELEMENTS_MAPPING%NUMBER_OF_GHOST)
    ENDIF
    
    ! assign numbers in domain mapping
    ELEMENTS_MAPPING%INTERNAL_START=1
    ELEMENTS_MAPPING%INTERNAL_FINISH=ELEMENTS_MAPPING%INTERNAL_START+ELEMENTS_MAPPING%NUMBER_OF_INTERNAL-1
    ELEMENTS_MAPPING%BOUNDARY_START=ELEMENTS_MAPPING%INTERNAL_FINISH+1
    ELEMENTS_MAPPING%BOUNDARY_FINISH=ELEMENTS_MAPPING%BOUNDARY_START+ELEMENTS_MAPPING%NUMBER_OF_BOUNDARY-1
    ELEMENTS_MAPPING%GHOST_START=ELEMENTS_MAPPING%BOUNDARY_FINISH+1
    ELEMENTS_MAPPING%GHOST_FINISH=ELEMENTS_MAPPING%GHOST_START+ELEMENTS_MAPPING%NUMBER_OF_GHOST-1
    
    ! assign local element numbers in domain_list
    ALLOCATE(ELEMENTS_MAPPING%DOMAIN_LIST(ELEMENTS_MAPPING%INTERNAL_START:ELEMENTS_MAPPING%GHOST_FINISH),STAT=ERR) 
    IF(ERR/=0) CALL FlagError("Could not allocate DOMAIN_LIST.",ERR,ERROR,*999)
    
    ! internal elements
    IF (ELEMENTS_MAPPING%NUMBER_OF_INTERNAL > 0) THEN
      ELEMENTS_MAPPING%DOMAIN_LIST(ELEMENTS_MAPPING%INTERNAL_START:ELEMENTS_MAPPING%INTERNAL_FINISH) = &
        & InternalElements(1:ELEMENTS_MAPPING%NUMBER_OF_INTERNAL)
    ENDIF
    
    ! boundary elements
    IF (ELEMENTS_MAPPING%NUMBER_OF_BOUNDARY > 0) THEN
      ELEMENTS_MAPPING%DOMAIN_LIST(ELEMENTS_MAPPING%BOUNDARY_START:ELEMENTS_MAPPING%BOUNDARY_FINISH) = &
        & BoundaryElements(1:ELEMENTS_MAPPING%NUMBER_OF_BOUNDARY)
    ENDIF
    
    ! ghost elements
    IF (ELEMENTS_MAPPING%NUMBER_OF_GHOST > 0) THEN
      ELEMENTS_MAPPING%DOMAIN_LIST(ELEMENTS_MAPPING%GHOST_START:ELEMENTS_MAPPING%GHOST_FINISH) = &
        & [(I, I=ELEMENTS_MAPPING%NUMBER_OF_LOCAL+1,ELEMENTS_MAPPING%TOTAL_NUMBER_OF_LOCAL,1)]
    ENDIF
           
    ! allocate number_of_domain_ghost
    ALLOCATE(ELEMENTS_MAPPING%NUMBER_OF_DOMAIN_GHOST(0:NumberComputationalNodes-1),STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate NUMBER_OF_DOMAIN_GHOST.",ERR,ERROR,*999)

    ELEMENTS_MAPPING%NUMBER_OF_DOMAIN_GHOST(MyComputationalNodeNumber) = ELEMENTS_MAPPING%NUMBER_OF_GHOST
          
    CALL MPI_ALLGATHER(MPI_IN_PLACE,1,MPI_INTEGER,ELEMENTS_MAPPING%&
      & NUMBER_OF_DOMAIN_GHOST(0:NumberComputationalNodes-1), &
      & 1,MPI_INTEGER,groupCommunicator,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_ALLGATHER",MPI_IERROR,ERR,ERROR,*999)
                   
    EXITS("DECOMPOSITION_ELEMENTS_MAPPING_CALCULATE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_ELEMENTS_MAPPING_CALCULATE",ERR,ERROR) 
    RETURN 1
  END SUBROUTINE DECOMPOSITION_ELEMENTS_MAPPING_CALCULATE

   !
   !================================================================================================================================
   !


  !>Calculates for decomposition ELEMENTS_MAPPING the variable NUMBER_OF_ADJACENT_DOMAINS and ADJACENT_DOMAINS(:) array
  SUBROUTINE DECOMPOSITION_ADJACENT_DOMAINS_CALCULATE(DECOMPOSITION,ERR,ERROR,*)
  
    !Argument variables
     TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to calculate the element domains for.  
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: MPI_IERROR, I, NumberComputationalNodes, MyComputationalNodeNumber, FirstBoundaryElementGlobalNo, &
      & LastBoundaryElementGlobalNo, ElementLocalNo, GhostElementLocalNo, ElementGlobalNo, GhostElementGlobalNo, &
      & NumberRemoteBoundaryElements, RemoteDomainNo, ComputationalNodeNo, MpiMemoryWindow, Sign, &
      & AdjacentDomainIdx, DomainIdx, ElementIdx, MpiRequest, AssertValue, MpiGroup, GhostSendIndexGlobalNo, LocalElementNo, &
      & MaximumReceiveBufferSize, MaximumSendBufferSize, FoundElementsCount, GhostIdx, &
      & FirstGhostElementGlobalNo, LastGhostElementGlobalNo, NumberRemoteGhostElements 
    INTEGER(INTG) :: groupCommunicator
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    INTEGER(INTG), ALLOCATABLE :: NumberBoundaryElements(:), FirstLastBoundaryElement(:), BoundaryElementGlobalNumber(:), &
      & RemoteBoundaryElements(:), SendRequestHandle(:), ReceiveRequestHandle(:), SendBuffers(:,:), ReceiveBuffers(:,:), &
      & NumberGhostElements(:), FirstLastGhostElement(:), GhostElementGlobalNumber(:), &
      & RemoteGhostElements(:) 
    INTEGER(KIND=MPI_ADDRESS_KIND) LB, EXTENT, MpiLocalWindowSizeBytes, MpiTargetDisplacement, MpiDisplacementsUnitBytes
    TYPE(LIST_TYPE), POINTER :: AdjacentDomainsList
    TYPE(DOMAIN_ADJACENT_DOMAIN_TYPE), POINTER :: AdjacentDomain
    TYPE(DOMAIN_ADJACENT_DOMAIN_TYPE), POINTER :: MappingAdjacentDomain
    TYPE(DOMAIN_ADJACENT_DOMAIN_TYPE), ALLOCATABLE :: AdjacentDomains(:)

    LOGICAL :: ElementFound, AdjacentDomainAlreadyStored, NotYetCheckedForLastRemoteBoundaryElements, &
      & NotYetCheckedForLastRemoteGhostElements, GhostIsAlreadyStored 

    ENTERS("DECOMPOSITION_ADJACENT_DOMAINS_CALCULATE",ERR,ERROR,*999) 

    CALL WorkGroup_GroupCommunicatorGet(decomposition%workGroup,groupCommunicator,err,error,*999)
    CALL WorkGroup_NumberOfGroupNodesGet(decomposition%workGroup,numberOfComputationNodes,err,error,*999)
    CALL WorkGroup_GroupNodeNumberGet(decomposition%workGroup,myComputationNodeNumber,err,error,*999)

    ! Replaced by work groups above!!  
    ! get rank count and own rank
    !NumberComputationalNodes=ComputationalEnvironment_NumberOfNodesGet(ERR,ERROR)
    !IF(ERR/=0) GOTO 999
    !MyComputationalNodeNumber=ComputationalEnvironment_NodeNumberGet(ERR,ERROR)
    !IF(ERR/=0) GOTO 999
    
    ELEMENTS_MAPPING=>DECOMPOSITION%ELEMENTS_MAPPING
    
    ! communicate number of boundary elements on each rank
    ALLOCATE(NumberBoundaryElements(0:NumberComputationalNodes-1),STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate NumberBoundaryElements.",ERR,ERROR,*999)
    NumberBoundaryElements(MyComputationalNodeNumber) = ELEMENTS_MAPPING%NUMBER_OF_BOUNDARY
    
    CALL MPI_ALLGATHER(MPI_IN_PLACE,1,MPI_INTEGER, &
      & NumberBoundaryElements(0:NumberComputationalNodes-1),1,MPI_INTEGER, &
      & groupCommunicator,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_ALLGATHER",MPI_IERROR,ERR,ERROR,*999)
    
    ! communicate bounds of boundary nodes, i.e. the first and last global element number of the contained boundary nodes
    ALLOCATE(FirstLastBoundaryElement(0:NumberComputationalNodes*2-1),STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate FirstLastBoundaryElement.",ERR,ERROR,*999)
    
    IF (ELEMENTS_MAPPING%NUMBER_OF_BOUNDARY == 0) THEN
      FirstBoundaryElementGlobalNo = 0
      LastBoundaryElementGlobalNo = 0
    ELSE
      FirstBoundaryElementGlobalNo = ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(ELEMENTS_MAPPING%DOMAIN_LIST(&
        &ELEMENTS_MAPPING%BOUNDARY_START))
      LastBoundaryElementGlobalNo = ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(ELEMENTS_MAPPING%DOMAIN_LIST(&
        & ELEMENTS_MAPPING%BOUNDARY_FINISH))
    ENDIF
    
    FirstLastBoundaryElement(2*MyComputationalNodeNumber) = FirstBoundaryElementGlobalNo
    FirstLastBoundaryElement(2*MyComputationalNodeNumber+1) = LastBoundaryElementGlobalNo
    
    CALL MPI_ALLGATHER(MPI_IN_PLACE,2,MPI_INTEGER, &
      & FirstLastBoundaryElement(0:NumberComputationalNodes*2-1),2,MPI_INTEGER, &
      & groupCommunicator,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_ALLGATHER",MPI_IERROR,ERR,ERROR,*999)
    
    !PRINT *, MyComputationalNodeNumber, ": NumberBoundaryElements: ", NumberBoundaryElements
    !PRINT *, MyComputationalNodeNumber, ": FirstLastBoundaryElement: ", FirstLastBoundaryElement
    
    ! provide own boundary element global numbers to all processes via RMA
    ALLOCATE(BoundaryElementGlobalNumber(ELEMENTS_MAPPING%NUMBER_OF_BOUNDARY),STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate BoundaryElementGlobalNumber.",ERR,ERROR,*999)
    
    I = 1
    DO ElementIdx=ELEMENTS_MAPPING%BOUNDARY_START,ELEMENTS_MAPPING%BOUNDARY_FINISH 
      ElementLocalNo = ELEMENTS_MAPPING%DOMAIN_LIST(ElementIdx)
      ElementGlobalNo = ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(ElementLocalNo)
      BoundaryElementGlobalNumber(I) = ElementGlobalNo
      I = I+1
    ENDDO
    
    !PRINT *, MyComputationalNodeNumber, ": BoundaryElementGlobalNumber: ", BoundaryElementGlobalNumber
    
    ! get size in bytes of MPI_INT
    CALL MPI_TYPE_GET_EXTENT(MPI_INT, LB, EXTENT, MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_TYPE_GET_EXTENT",MPI_IERROR,ERR,ERROR,*999)
    
    ! start mpi RMA
    MpiLocalWindowSizeBytes   = ELEMENTS_MAPPING%NUMBER_OF_BOUNDARY*EXTENT
    MpiDisplacementsUnitBytes = EXTENT
    CALL MPI_WIN_CREATE(BoundaryElementGlobalNumber(1:ELEMENTS_MAPPING%NUMBER_OF_BOUNDARY), MpiLocalWindowSizeBytes, &
      & MpiDisplacementsUnitBytes, MPI_INFO_NULL, groupCommunicator, MpiMemoryWindow, MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_WIN_CREATE",MPI_IERROR,ERR,ERROR,*999)
    
    !CALL MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, MpiMemoryWindow, MPI_IERROR)
    !CALL MPI_ERROR_CHECK("MPI_WIN_FENCE",MPI_IERROR,ERR,ERROR,*999)
        
    AssertValue = 0
    !CALL MPI_COMM_GROUP(groupCommunicator, MpiGroup, MPI_IERROR)
    !CALL MPI_ERROR_CHECK("MPI_WIN_START",MPI_IERROR,ERR,ERROR,*999)
    !CALL MPI_WIN_START(MpiGroup, AssertValue, MpiMemoryWindow, MPI_IERROR)
    !CALL MPI_ERROR_CHECK("MPI_WIN_START",MPI_IERROR,ERR,ERROR,*999)
        
    ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS = 0
    
    ! search for own ghost elements at remote ranks
    ! loop over ghost elements
    DO ElementIdx=ELEMENTS_MAPPING%GHOST_START,ELEMENTS_MAPPING%GHOST_FINISH
      GhostElementLocalNo=ELEMENTS_MAPPING%DOMAIN_LIST(ElementIdx)
      GhostElementGlobalNo=ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(GhostElementLocalNo)
      

      IF (DIAGNOSTICS2) THEN

        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"GhostElement local "// &
            & TRIM(NumberToVString(GhostElementLocalNo,"*",ERR,ERROR))//", global "// &
            & TRIM(NumberToVString(GhostElementGlobalNo,"*",ERR,ERROR)),ERR,ERROR,*999)
        !PRINT *, MyComputationalNodeNumber, ": GhostElement local ",GhostElementLocalNo, ", global ",GhostElementGlobalNo 
      ENDIF
      
      NotYetCheckedForLastRemoteBoundaryElements = .TRUE. 
      
      ! loop over other ranks to see which one could have the element as boundary element
      ! iteration order is r-1,r+1,r-2,r+2,... with r being own computational node
      DO I=1,NumberComputationalNodes
        DO Sign=-1,1,2
          ComputationalNodeNo=MyComputationalNodeNumber+I*Sign 
          ComputationalNodeNo=MODULO(ComputationalNodeNo,NumberComputationalNodes)
        
          ElementFound = .FALSE.
          IF (ALLOCATED(RemoteBoundaryElements) .AND. NotYetCheckedForLastRemoteBoundaryElements) THEN
          
            ! test if the previously received remote boundary nodes list satisfies the current ghost element request
            IF (SORTED_ARRAY_CONTAINS_ELEMENT(RemoteBoundaryElements,GhostElementGlobalNo,ERR,ERROR)) THEN
              ElementFound = .TRUE.
              NotYetCheckedForLastRemoteBoundaryElements = .FALSE.


              IF (DIAGNOSTICS2) THEN              
                CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"element is in list from last rank ("// &
                    & TRIM(NumberToVString(RemoteDomainNo,"*",ERR,ERROR))//") -> element found",ERR,ERROR,*999)
                !PRINT *, MyComputationalNodeNumber, ": element is in list from last rank (",RemoteDomainNo,") -> element found" 
              ENDIF 
           ENDIF
          ENDIF
         
          IF (.NOT. ElementFound) THEN
            IF (DIAGNOSTICS2) THEN

              CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"query rank "// &
                  & TRIM(NumberToVString(ComputationalNodeNo,"*",ERR,ERROR)),ERR,ERROR,*999)
              !PRINT *, MyComputationalNodeNumber, ": query rank ", ComputationalNodeNo 
            ENDIF

            FirstBoundaryElementGlobalNo = FirstLastBoundaryElement(2*ComputationalNodeNo)
            LastBoundaryElementGlobalNo = FirstLastBoundaryElement(2*ComputationalNodeNo+1)
            IF (FirstBoundaryElementGlobalNo <= GhostElementGlobalNo .AND. GhostElementGlobalNo <= LastBoundaryElementGlobalNo) THEN
                
              ! communicate to see if rank contains node

              IF (DIAGNOSTICS2) THEN

                CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"element lies in interval ["// &
                    & TRIM(NumberToVString(FirstBoundaryElementGlobalNo,"*",ERR,ERROR))//","// &
                    & TRIM(NumberToVString(LastBoundaryElementGlobalNo,"*",ERR,ERROR))//"]",ERR,ERROR,*999)
                !PRINT *, MyComputationalNodeNumber, ": element lies in interval [",FirstBoundaryElementGlobalNo,",",&
                !  & LastBoundaryElementGlobalNo,"]" 
              ENDIF 
                
              IF (ALLOCATED(RemoteBoundaryElements)) DEALLOCATE(RemoteBoundaryElements)
                
              ! read list with all boundary elements of rank ComputationalNodeNo into RemoteBoundaryElements
              NumberRemoteBoundaryElements = NumberBoundaryElements(ComputationalNodeNo)
              ALLOCATE(RemoteBoundaryElements(NumberRemoteBoundaryElements),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate RemoteBoundaryElements.",ERR,ERROR,*999)
              RemoteBoundaryElements = 0

              !CALL MPI_WIN_POST(MpiGroup,0,MpiMemoryWindow, MPI_IERROR)
              !CALL MPI_ERROR_CHECK("MPI_WIN_POST",MPI_IERROR,ERR,ERROR,*999)               

              ! start passive target communication (see http://www.mcs.anl.gov/research/projects/mpi/mpi-standard/mpi-report-2.0/node126.htm for introduction)
              CALL MPI_WIN_LOCK(MPI_LOCK_SHARED, ComputationalNodeNo, 0, MpiMemoryWindow, MPI_IERROR)
              CALL MPI_ERROR_CHECK("MPI_WIN_LOCK",MPI_IERROR,ERR,ERROR,*999)               


              MpiTargetDisplacement = 0
              CALL MPI_GET(RemoteBoundaryElements, NumberRemoteBoundaryElements, MPI_INT, &
                & ComputationalNodeNo, MpiTargetDisplacement, NumberRemoteBoundaryElements, MPI_INT, MpiMemoryWindow, &
                & MPI_IERROR)
              CALL MPI_ERROR_CHECK("MPI_GET",MPI_IERROR,ERR,ERROR,*999) 

              !CALL MPI_WIN_FLUSH(ComputationalNodeNo, MpiMemoryWindow, MPI_IERROR)
              !CALL MPI_ERROR_CHECK("MPI_WIN_FLUSH",MPI_IERROR,ERR,ERROR,*999) 

              CALL MPI_WIN_UNLOCK(ComputationalNodeNo, MpiMemoryWindow, MPI_IERROR)
              CALL MPI_ERROR_CHECK("MPI_WIN_UNLOCK",MPI_IERROR,ERR,ERROR,*999)
              
              !CALL MPI_WIN_WAIT(MpiMemoryWindow, MPI_IERROR)
              !CALL MPI_ERROR_CHECK("MPI_WIN_WAIT",MPI_IERROR,ERR,ERROR,*999)
              
              ! synchronize RMA
              !CALL MPI_WIN_FENCE(IOR(MPI_MODE_NOSTORE, MPI_MODE_NOPUT), MpiMemoryWindow, MPI_IERROR)
              !CALL MPI_ERROR_CHECK("MPI_WIN_FENCE",MPI_IERROR,ERR,ERROR,*999) 
           

              IF (DIAGNOSTICS2) THEN

                CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"find ghost "// &
                    & TRIM(NumberToVString(GhostElementGlobalNo,"*",ERR,ERROR))//" at rank "// &
                    & TRIM(NumberToVString(ComputationalNodeNo,"*",ERR,ERROR))//" which has "// &
                    & TRIM(NumberToVString(NumberRemoteBoundaryElements,"*",ERR,ERROR))//" boundary nodes",ERR,ERROR,*999)
                !PRINT *, MyComputationalNodeNumber, ": find ghost ",GhostElementGlobalNo," at rank ", ComputationalNodeNo, &
                !  & " which has ", NumberRemoteBoundaryElements, " boundary nodes :", RemoteBoundaryElements 
              ENDIF
              
              RemoteDomainNo = ComputationalNodeNo
              NotYetCheckedForLastRemoteBoundaryElements = .TRUE. 

              ! if element was found 
              IF (SORTED_ARRAY_CONTAINS_ELEMENT(RemoteBoundaryElements,GhostElementGlobalNo,ERR,ERROR)) THEN

                IF (DIAGNOSTICS2) THEN

                  CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"RemoteBoundaryElements contains ghost element -> element found",&
                    & ERR,ERROR,*999)
                  !PRINT *, MyComputationalNodeNumber, ": RemoteBoundaryElements contains ghost element -> element found" 
                ENDIF 
                ElementFound = .TRUE.
              ELSE
                IF (DIAGNOSTICS2) THEN

                  CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"RemoteBoundaryElements does not contain ghost element",&
                    & ERR,ERROR,*999)
                  !PRINT *, MyComputationalNodeNumber, ": RemoteBoundaryElements does not contain ghost element" 
                ENDIF
              ENDIF
            ELSE
              IF (DIAGNOSTICS2) THEN

                CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"boundary element does not lie in range of computational node",&
                  & ERR,ERROR,*999)
                !PRINT *, MyComputationalNodeNumber, ": boundary element does not lie in range of computational node" 
              ENDIF
              
            ENDIF ! boundary element no. lies in range of computational node
          ENDIF   ! .NOT. ElementFound
              
              
          IF (ElementFound) THEN
            ! loop over already collected adjacent domains to see if this one was already stored
            AdjacentDomainAlreadyStored = .FALSE.
            DO AdjacentDomainIdx=1,ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
              MappingAdjacentDomain=>ELEMENTS_MAPPING%ADJACENT_DOMAINS(AdjacentDomainIdx)
              IF (MappingAdjacentDomain%DOMAIN_NUMBER == RemoteDomainNo) THEN

                IF (DIAGNOSTICS2) THEN

                  CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"add ghost to AdjacentDomainIdx "// &
                      & TRIM(NumberToVString(AdjacentDomainIdx,"*",ERR,ERROR)),ERR,ERROR,*999)
                  !PRINT *, MyComputationalNodeNumber, ": add ghost to AdjacentDomainIdx ",AdjacentDomainIdx 
                ENDIF
            
                ! add ghost 
                MappingAdjacentDomain%NUMBER_OF_RECEIVE_GHOSTS = MappingAdjacentDomain%NUMBER_OF_RECEIVE_GHOSTS+1
                
                MappingAdjacentDomain%LOCAL_GHOST_RECEIVE_INDICES = RESHAPE(MappingAdjacentDomain%LOCAL_GHOST_RECEIVE_INDICES, &
                  & [MappingAdjacentDomain%NUMBER_OF_RECEIVE_GHOSTS], [GhostElementLocalNo])
                
                AdjacentDomainAlreadyStored = .TRUE.
                EXIT               

              ENDIF
            ENDDO
            
            IF (.NOT.AdjacentDomainAlreadyStored) THEN

              IF (DIAGNOSTICS2) THEN

                CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"create new AdjacentDomain",ERR,ERROR,*999)
                !PRINT *, MyComputationalNodeNumber, ": create new AdjacentDomain" 
              ENDIF

              ! create new AdjacentDomain entry in ELEMENTS_MAPPING%ADJACENT_DOMAINS  
            
              NULLIFY(AdjacentDomain)
              ALLOCATE(AdjacentDomain,STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate AdjacentDomain.",ERR,ERROR,*999) 

              AdjacentDomain%DOMAIN_NUMBER = RemoteDomainNo
              AdjacentDomain%NUMBER_OF_SEND_GHOSTS = 0
              AdjacentDomain%NUMBER_OF_RECEIVE_GHOSTS = 0
              AdjacentDomain%NUMBER_OF_FURTHER_LINKED_GHOSTS = 0 
              ALLOCATE(AdjacentDomain%LOCAL_GHOST_RECEIVE_INDICES(1),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate LOCAL_GHOST_RECEIVE_INDICES.",ERR,ERROR,*999)
              AdjacentDomain%NUMBER_OF_RECEIVE_GHOSTS = 1
              AdjacentDomain%LOCAL_GHOST_RECEIVE_INDICES(1) = GhostElementLocalNo              
              
              ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS = ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS+1
              

              IF (DIAGNOSTICS2) THEN               
                CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"reshape to size "// &
                  & TRIM(NumberToVString(ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS,"*",ERR,ERROR)),ERR,ERROR,*999)
                !PRINT *, MyComputationalNodeNumber, ": reshape to size ",ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS 
              ENDIF

              ALLOCATE(AdjacentDomains(ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate AdjacentDomains.",ERR,ERROR,*999)
              
              ! if there is a previous entry in ELEMENTS_MAPPING%ADJACENT_DOMAINS, save them
              IF (ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS > 1) THEN
                AdjacentDomains(1:ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS-1)=ELEMENTS_MAPPING%ADJACENT_DOMAINS
              ENDIF
              AdjacentDomains(ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS)=AdjacentDomain 

              CALL MOVE_ALLOC(AdjacentDomains,ELEMENTS_MAPPING%ADJACENT_DOMAINS)
          
            ENDIF ! .NOT.AdjacentDomainAlreadyStored 

          ELSE
            IF (DIAGNOSTICS2) THEN

              CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"element not found on rank "// &
                & TRIM(NumberToVString(ComputationalNodeNo,"*",ERR,ERROR)),ERR,ERROR,*999)
              !PRINT *, MyComputationalNodeNumber, ": element not found on rank ", ComputationalNodeNo 
            ENDIF 
          ENDIF ! element was found
            
          IF (ElementFound) EXIT
        ENDDO ! sign
        IF (ElementFound) EXIT
      ENDDO ! I
    ENDDO ! GhostElementLocalNo
        
          
    ! synchronize RMA
    CALL MPI_WIN_FENCE(0, MpiMemoryWindow, MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_WIN_FENCE",MPI_IERROR,ERR,ERROR,*999)
          
    ! deallocate window 
    CALL MPI_WIN_FREE(MpiMemoryWindow, MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_WIN_FREE",MPI_IERROR,ERR,ERROR,*999)
            
    ! exchange LOCAL_GHOST_SEND_INDICES
    ALLOCATE(SendRequestHandle(ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS),STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate SendRequestHandle.",ERR,ERROR,*999)
    
    ALLOCATE(ReceiveRequestHandle(ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS),STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate ReceiveRequestHandle.",ERR,ERROR,*999)
    
    MaximumSendBufferSize = 0
    
    ! commit send commands for number of elements to transfer
    ! loop over adjacent domains
    DO AdjacentDomainIdx=1,ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
      MappingAdjacentDomain=>ELEMENTS_MAPPING%ADJACENT_DOMAINS(AdjacentDomainIdx)
      
      CALL MPI_ISEND(MappingAdjacentDomain%NUMBER_OF_RECEIVE_GHOSTS, &
        & 1, MPI_INT, MappingAdjacentDomain%DOMAIN_NUMBER, 0, &
        & groupCommunicator, SendRequestHandle(AdjacentDomainIdx), MPI_IERROR)
      CALL MPI_ERROR_CHECK("MPI_ISEND",MPI_IERROR,ERR,ERROR,*999)
      MaximumSendBufferSize = MAX(MaximumSendBufferSize, MappingAdjacentDomain%NUMBER_OF_RECEIVE_GHOSTS)
    ENDDO
    
    ! commit receive commands for number of elements to transfer
    ! loop over adjacent domains
    DO AdjacentDomainIdx=ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS,1,-1
      MappingAdjacentDomain=>ELEMENTS_MAPPING%ADJACENT_DOMAINS(AdjacentDomainIdx)
      
      CALL MPI_IRECV(MappingAdjacentDomain%NUMBER_OF_SEND_GHOSTS, &
        & 1, MPI_INT, MappingAdjacentDomain%DOMAIN_NUMBER, 0, &
        & groupCommunicator, ReceiveRequestHandle(AdjacentDomainIdx), MPI_IERROR)
      CALL MPI_ERROR_CHECK("MPI_IRECV",MPI_IERROR,ERR,ERROR,*999)
    ENDDO
    
    CALL MPI_WAITALL(ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS, SendRequestHandle, MPI_STATUSES_IGNORE, MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_WAITALL",MPI_IERROR,ERR,ERROR,*999)
    CALL MPI_WAITALL(ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS, ReceiveRequestHandle, MPI_STATUSES_IGNORE, MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_WAITALL",MPI_IERROR,ERR,ERROR,*999)
    
    ! determine maximum receive buffer size
    MaximumReceiveBufferSize = 0
    DO AdjacentDomainIdx=ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS,1,-1
      MappingAdjacentDomain=>ELEMENTS_MAPPING%ADJACENT_DOMAINS(AdjacentDomainIdx)
      
      MaximumReceiveBufferSize = MAX(MaximumReceiveBufferSize, MappingAdjacentDomain%NUMBER_OF_SEND_GHOSTS)
      
      ! allocate space for LOCAL_GHOST_SEND_INDICES
      ALLOCATE(MappingAdjacentDomain%LOCAL_GHOST_SEND_INDICES(1:MappingAdjacentDomain%NUMBER_OF_SEND_GHOSTS),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate LOCAL_GHOST_SEND_INDICES of size "//&
        & TRIM(NUMBER_TO_VSTRING(MappingAdjacentDomain%NUMBER_OF_SEND_GHOSTS,"*",ERR,ERROR))//".",ERR,ERROR,*999)
    ENDDO
    
    ! allocate send and receive buffers
    ALLOCATE(SendBuffers(MaximumSendBufferSize,ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS),STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate SendBuffers.",ERR,ERROR,*999)
    ALLOCATE(ReceiveBuffers(MaximumReceiveBufferSize,ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS),STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate ReceiveBuffers.",ERR,ERROR,*999)
    
    
    ! commit send commands for actual ghost element global numbers
    ! loop over adjacent domains
    DO AdjacentDomainIdx=1,ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
      MappingAdjacentDomain=>ELEMENTS_MAPPING%ADJACENT_DOMAINS(AdjacentDomainIdx)
              
      ! fill send buffer
      DO I=1,MappingAdjacentDomain%NUMBER_OF_RECEIVE_GHOSTS
        SendBuffers(I,AdjacentDomainIdx)=&
          & ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(MappingAdjacentDomain%LOCAL_GHOST_RECEIVE_INDICES(I))     
      ENDDO
      
      CALL MPI_ISEND(SendBuffers(1:MappingAdjacentDomain%NUMBER_OF_RECEIVE_GHOSTS,AdjacentDomainIdx), &
        & MappingAdjacentDomain%NUMBER_OF_RECEIVE_GHOSTS, MPI_INT, MappingAdjacentDomain%DOMAIN_NUMBER, 0, &
        & groupCommunicator, SendRequestHandle(AdjacentDomainIdx), MPI_IERROR)
      CALL MPI_ERROR_CHECK("MPI_ISEND",MPI_IERROR,ERR,ERROR,*999)
    ENDDO
    
    ! commit receive commands for number of elements to transfer
    ! loop over adjacent domains
    DO AdjacentDomainIdx=ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS,1,-1
      MappingAdjacentDomain=>ELEMENTS_MAPPING%ADJACENT_DOMAINS(AdjacentDomainIdx)
      
      CALL MPI_IRECV(ReceiveBuffers(1:MappingAdjacentDomain%NUMBER_OF_SEND_GHOSTS,AdjacentDomainIdx), &
        & MappingAdjacentDomain%NUMBER_OF_SEND_GHOSTS, MPI_INT, MappingAdjacentDomain%DOMAIN_NUMBER, 0, &
        & groupCommunicator, ReceiveRequestHandle(AdjacentDomainIdx), MPI_IERROR)
      CALL MPI_ERROR_CHECK("MPI_IRECV",MPI_IERROR,ERR,ERROR,*999)
    ENDDO 
    

    CALL MPI_WAITALL(ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS, SendRequestHandle, MPI_STATUSES_IGNORE, MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_WAITALL",MPI_IERROR,ERR,ERROR,*999)
    CALL MPI_WAITALL(ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS, ReceiveRequestHandle, MPI_STATUSES_IGNORE, MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_WAITALL",MPI_IERROR,ERR,ERROR,*999)
    
    ! convert received global element numbers back to local element numbers
    ! loop over adjacent domains
    DO AdjacentDomainIdx=1,ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
      MappingAdjacentDomain=>ELEMENTS_MAPPING%ADJACENT_DOMAINS(AdjacentDomainIdx)
      
      LocalElementNo = 1
      DO I=1,MappingAdjacentDomain%NUMBER_OF_SEND_GHOSTS
        GhostSendIndexGlobalNo = ReceiveBuffers(I,AdjacentDomainIdx)
        
        ! convert global element number GhostSendIndexGlobalNo to local element number LocalElementNo
        DO WHILE (ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(LocalELementNo) < GhostSendIndexGlobalNo)
          LocalElementNo = LocalElementNo+1
        ENDDO
      
        MappingAdjacentDomain%LOCAL_GHOST_SEND_INDICES(I) = LocalElementNo
      ENDDO
    ENDDO
    
    ! --------------------------------------------------------
    ! find all ghost elements on other ranks that are ghost elements on this rank as well and store them to LOCAL_GHOST_FURTHER_INDICES
    
    ! communicate number of ghost elements on each rank
    ALLOCATE(NumberGhostElements(0:NumberComputationalNodes-1),STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate NumberGhostElements.",ERR,ERROR,*999)
    NumberGhostElements(MyComputationalNodeNumber) = ELEMENTS_MAPPING%NUMBER_OF_GHOST
    
    CALL MPI_ALLGATHER(MPI_IN_PLACE,1,MPI_INTEGER, &
      & NumberGhostElements(0:NumberComputationalNodes-1),1,MPI_INTEGER, &
      & groupCommunicator,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_ALLGATHER",MPI_IERROR,ERR,ERROR,*999)
    
    ! communicate bounds of ghost nodes, i.e. the first and last global element number of the contained ghost nodes
    ALLOCATE(FirstLastGhostElement(0:NumberComputationalNodes*2-1),STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate FirstLastGhostElement.",ERR,ERROR,*999)
    
    IF (ELEMENTS_MAPPING%NUMBER_OF_GHOST == 0) THEN
      FirstGhostElementGlobalNo = 0
      LastGhostElementGlobalNo = 0
    ELSE
      FirstGhostElementGlobalNo = ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(ELEMENTS_MAPPING%DOMAIN_LIST(&
        &ELEMENTS_MAPPING%GHOST_START))
      LastGhostElementGlobalNo = ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(ELEMENTS_MAPPING%DOMAIN_LIST(&
        & ELEMENTS_MAPPING%GHOST_FINISH))
    ENDIF
    
    FirstLastGhostElement(2*MyComputationalNodeNumber) = FirstGhostElementGlobalNo
    FirstLastGhostElement(2*MyComputationalNodeNumber+1) = LastGhostElementGlobalNo
    
    CALL MPI_ALLGATHER(MPI_IN_PLACE,2,MPI_INTEGER, &
      & FirstLastGhostElement(0:NumberComputationalNodes*2-1),2,MPI_INTEGER, &
      & groupCommunicator,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_ALLGATHER",MPI_IERROR,ERR,ERROR,*999)
    
    ! provide own ghost element global numbers to all processes via RMA
    ALLOCATE(GhostElementGlobalNumber(ELEMENTS_MAPPING%NUMBER_OF_GHOST),STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate GhostElementGlobalNumber.",ERR,ERROR,*999)
    
    I = 1
    DO ElementIdx=ELEMENTS_MAPPING%GHOST_START,ELEMENTS_MAPPING%GHOST_FINISH 
      ElementLocalNo = ELEMENTS_MAPPING%DOMAIN_LIST(ElementIdx)
      ElementGlobalNo = ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(ElementLocalNo)
      GhostElementGlobalNumber(I) = ElementGlobalNo
      I = I+1
    ENDDO
    
      IF (DIAGNOSTICS2) THEN

        !CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"*my GhostElements global: "// &
        !  & TRIM(NumberToVString(GhostElementGlobalNumber,"*",ERR,ERROR)),ERR,ERROR,*999)
        !PRINT *, MyComputationalNodeNumber, ": *my GhostElements global: ",GhostElementGlobalNumber 
      ENDIF
    
    ! start mpi RMA
    MpiLocalWindowSizeBytes = ELEMENTS_MAPPING%NUMBER_OF_GHOST*EXTENT
    MpiDisplacementsUnitBytes = EXTENT
    CALL MPI_WIN_CREATE(GhostElementGlobalNumber(1:ELEMENTS_MAPPING%NUMBER_OF_GHOST), MpiLocalWindowSizeBytes, &
      & MpiDisplacementsUnitBytes, MPI_INFO_NULL, groupCommunicator, MpiMemoryWindow, MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_WIN_CREATE",MPI_IERROR,ERR,ERROR,*999)
            
    AssertValue = 0
    
    ! search for own ghost elements at remote ranks
    ! loop over ghost elements
    DO ElementIdx=ELEMENTS_MAPPING%GHOST_START,ELEMENTS_MAPPING%GHOST_FINISH
      GhostElementLocalNo=ELEMENTS_MAPPING%DOMAIN_LIST(ElementIdx)
      GhostElementGlobalNo=ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(GhostElementLocalNo)
      
      IF (DIAGNOSTICS2) THEN

        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"*GhostElement local "// &
          & TRIM(NumberToVString(GhostElementLocalNo,"*",ERR,ERROR))//", global "// &
          & TRIM(NumberToVString(GhostElementGlobalNo,"*",ERR,ERROR)),ERR,ERROR,*999)
        !PRINT *, MyComputationalNodeNumber, ": *GhostElement local ",GhostElementLocalNo, ", global ",GhostElementGlobalNo 
      ENDIF
      
      NotYetCheckedForLastRemoteGhostElements = .TRUE.
      FoundElementsCount = 0
      
      ! loop over other ranks to see which one has the element as ghost element
      ! iteration order is r-1,r+1,r-2,r+2,... with r being own computational node
      DO I=1,NumberComputationalNodes
        DO Sign=-1,1,2
          ComputationalNodeNo=MyComputationalNodeNumber+I*Sign 
          ComputationalNodeNo=MODULO(ComputationalNodeNo,NumberComputationalNodes)
        
          IF (ComputationalNodeNo == MyComputationalNodeNumber) CYCLE
        
          ElementFound = .FALSE.
          IF (ALLOCATED(RemoteGhostElements) .AND. NotYetCheckedForLastRemoteGhostElements) THEN
          
            ! test if the previously received remote ghost elements list satisfies the current ghost element request
            IF (SORTED_ARRAY_CONTAINS_ELEMENT(RemoteGhostElements,GhostElementGlobalNo,ERR,ERROR)) THEN
              ElementFound = .TRUE.
              NotYetCheckedForLastRemoteGhostElements = .FALSE.
              IF (DIAGNOSTICS2) THEN

                CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"*element is in list from last rank ("// &
                  & TRIM(NumberToVString(RemoteDomainNo,"*",ERR,ERROR))//") -> element found",ERR,ERROR,*999)
                !PRINT *, MyComputationalNodeNumber, ": *element is in list from last rank (",RemoteDomainNo,") -> element found" 
              ENDIF
            ENDIF
          ENDIF
          
          
          IF (.NOT. ElementFound) THEN
            IF (DIAGNOSTICS2) THEN

              CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"*query rank "// &
                & TRIM(NumberToVString(ComputationalNodeNo,"*",ERR,ERROR)),ERR,ERROR,*999)
              !PRINT *, MyComputationalNodeNumber, ": *query rank ", ComputationalNodeNo 
            ENDIF
          
            FirstGhostElementGlobalNo = FirstLastGhostElement(2*ComputationalNodeNo)
            LastGhostElementGlobalNo = FirstLastGhostElement(2*ComputationalNodeNo+1)
            IF (FirstGhostElementGlobalNo <= GhostElementGlobalNo .AND. GhostElementGlobalNo <= LastGhostElementGlobalNo) THEN
                
              ! communicate to see if rank contains element
              
              IF (DIAGNOSTICS2) THEN

                CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"*element lies in interval ["// &
                  & TRIM(NumberToVString(FirstGhostElementGlobalNo,"*",ERR,ERROR))//","// &
                  & TRIM(NumberToVString(LastGhostElementGlobalNo,"*",ERR,ERROR))//"]" &
                  & ,ERR,ERROR,*999)
                !PRINT *, MyComputationalNodeNumber, ": *element lies in interval [",FirstGhostElementGlobalNo,",",&
                !  & LastGhostElementGlobalNo,"]" 
              ENDIF
                
              IF (ALLOCATED(RemoteGhostElements)) DEALLOCATE(RemoteGhostElements)
                
              ! read list with all ghost elements of rank ComputationalNodeNo into RemoteGhostElements
              NumberRemoteGhostElements = NumberGhostElements(ComputationalNodeNo)
              ALLOCATE(RemoteGhostElements(NumberRemoteGhostElements),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate RemoteGhostElements.",ERR,ERROR,*999)
              RemoteGhostElements = 0
              
              ! start passive target communication (see http://www.mcs.anl.gov/research/projects/mpi/mpi-standard/mpi-report-2.0/node126.htm for introduction)
              CALL MPI_WIN_LOCK(MPI_LOCK_SHARED, ComputationalNodeNo, 0, MpiMemoryWindow, MPI_IERROR)
              CALL MPI_ERROR_CHECK("MPI_WIN_LOCK",MPI_IERROR,ERR,ERROR,*999)
              
              MpiTargetDisplacement = 0
              CALL MPI_GET(RemoteGhostElements, NumberRemoteGhostElements, MPI_INT, &
                & ComputationalNodeNo, MpiTargetDisplacement, NumberRemoteGhostElements, MPI_INT, MpiMemoryWindow, &
                & MPI_IERROR)
              CALL MPI_ERROR_CHECK("MPI_GET",MPI_IERROR,ERR,ERROR,*999)
                  
              CALL MPI_WIN_UNLOCK(ComputationalNodeNo, MpiMemoryWindow, MPI_IERROR)
              CALL MPI_ERROR_CHECK("MPI_WIN_UNLOCK",MPI_IERROR,ERR,ERROR,*999)

    
              IF (DIAGNOSTICS2) THEN

                CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"*find ghost "// &
                  & TRIM(NumberToVString(GhostElementGlobalNo,"*",ERR,ERROR))//" at rank "// &
                  & TRIM(NumberToVString(ComputationalNodeNo,"*",ERR,ERROR))//" which has "// &
                  & TRIM(NumberToVString(NumberRemoteGhostElements,"*",ERR,ERROR))//" *ghost nodes" &
                  & ,ERR,ERROR,*999)
                !PRINT *, MyComputationalNodeNumber, ": *find ghost ",GhostElementGlobalNo," at rank ", ComputationalNodeNo, &
                !  & " which has ", NumberRemoteGhostElements, " *ghost nodes :", RemoteGhostElements 
              ENDIF
              
              RemoteDomainNo = ComputationalNodeNo
              NotYetCheckedForLastRemoteGhostElements = .TRUE.
                
              ! if element was found 
              IF (SORTED_ARRAY_CONTAINS_ELEMENT(RemoteGhostElements,GhostElementGlobalNo,ERR,ERROR)) THEN
                IF (DIAGNOSTICS2) THEN

                  CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"*RemoteGhostElements contains ghost element -> element found" &
                    & ,ERR,ERROR,*999)
                  !PRINT *, MyComputationalNodeNumber, ": *RemoteGhostElements contains ghost element -> element found" 
                ENDIF
                ElementFound = .TRUE.
              ELSE
                IF (DIAGNOSTICS2) THEN

                  CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"*RemoteGhostElements does not contain ghost element" &
                    & ,ERR,ERROR,*999)
                  !PRINT *, MyComputationalNodeNumber, ": *RemoteGhostElements does not contain ghost element" 
                ENDIF
              ENDIF
            ELSE
              IF (DIAGNOSTICS2) THEN

                CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"*ghost element does not lie in range of computational node" &
                  & ,ERR,ERROR,*999)
                !PRINT *, MyComputationalNodeNumber, ": *ghost element does not lie in range of computational node" 
              ENDIF
              
            ENDIF ! ghost element no. lies in range of computational node
          ENDIF   ! .NOT. ElementFound
              
              
          IF (ElementFound) THEN
            FoundElementsCount = FoundElementsCount + 1
          
            ! loop over already collected adjacent domains to see if this one was already stored
            AdjacentDomainAlreadyStored = .FALSE.
            DO AdjacentDomainIdx=1,ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
              MappingAdjacentDomain=>ELEMENTS_MAPPING%ADJACENT_DOMAINS(AdjacentDomainIdx)
              IF (MappingAdjacentDomain%DOMAIN_NUMBER == RemoteDomainNo) THEN
              
                IF (DIAGNOSTICS2) THEN

                  CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"*add ghost to AdjacentDomainIdx "// &
                    & TRIM(NumberToVString(AdjacentDomainIdx,"*",ERR,ERROR)) &
                    & ,ERR,ERROR,*999)
                  !PRINT *, MyComputationalNodeNumber, ": *add ghost to AdjacentDomainIdx ",AdjacentDomainIdx 
                ENDIF
            
                ! check if ghost is already stored
                GhostIsAlreadyStored = .FALSE.
                DO GhostIdx=1,MappingAdjacentDomain%NUMBER_OF_FURTHER_LINKED_GHOSTS
                  IF (MappingAdjacentDomain%LOCAL_GHOST_FURTHER_INDICES(GhostIdx) == GhostElementLocalNo) THEN
                    GhostIsAlreadyStored = .TRUE.
                    EXIT
                  ENDIF
                ENDDO
            
                IF (.NOT. GhostIsAlreadyStored) THEN
                  ! add ghost 
                  MappingAdjacentDomain%NUMBER_OF_FURTHER_LINKED_GHOSTS = MappingAdjacentDomain%NUMBER_OF_FURTHER_LINKED_GHOSTS+1
                  
                  IF (DIAGNOSTICS2) THEN

                    CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"*NUMBER_OF_FURTHER_LINKED_GHOSTS "// &
                      & TRIM(NumberToVString(MappingAdjacentDomain%NUMBER_OF_FURTHER_LINKED_GHOSTS,"*",ERR,ERROR)) &
                      & ,ERR,ERROR,*999)
                    !PRINT *, MyComputationalNodeNumber, ": *NUMBER_OF_FURTHER_LINKED_GHOSTS ", &
                    !  & MappingAdjacentDomain%NUMBER_OF_FURTHER_LINKED_GHOSTS 
                  ENDIF
                  
                  IF (.NOT.ALLOCATED(MappingAdjacentDomain%LOCAL_GHOST_FURTHER_INDICES) .AND. &
                    & MappingAdjacentDomain%NUMBER_OF_FURTHER_LINKED_GHOSTS == 1) THEN
                                  
                    ALLOCATE(MappingAdjacentDomain%LOCAL_GHOST_FURTHER_INDICES(1),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate LOCAL_GHOST_FURTHER_INDICES.",ERR,ERROR,*999)
                    
                    MappingAdjacentDomain%LOCAL_GHOST_FURTHER_INDICES(1) = GhostElementLocalNo
                  ELSE
                    MappingAdjacentDomain%LOCAL_GHOST_FURTHER_INDICES = RESHAPE(MappingAdjacentDomain%LOCAL_GHOST_FURTHER_INDICES, &
                      & [MappingAdjacentDomain%NUMBER_OF_FURTHER_LINKED_GHOSTS], [GhostElementLocalNo])
                      
                  ENDIF
                ENDIF
                
                AdjacentDomainAlreadyStored = .TRUE.
                EXIT
              ENDIF
            ENDDO
            
            IF (.NOT.AdjacentDomainAlreadyStored) THEN
            
              IF (DIAGNOSTICS2) THEN

                CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"*create new AdjacentDomain" &
                  & ,ERR,ERROR,*999)
                !PRINT *, MyComputationalNodeNumber, ": *create new AdjacentDomain" 
              ENDIF
            
              ! create new AdjacentDomain entry in ELEMENTS_MAPPING%ADJACENT_DOMAINS
              NULLIFY(AdjacentDomain)
              ALLOCATE(AdjacentDomain,STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate AdjacentDomain.",ERR,ERROR,*999)
              
              AdjacentDomain%DOMAIN_NUMBER = RemoteDomainNo
              AdjacentDomain%NUMBER_OF_SEND_GHOSTS = 0
              AdjacentDomain%NUMBER_OF_RECEIVE_GHOSTS = 0
              AdjacentDomain%NUMBER_OF_FURTHER_LINKED_GHOSTS = 0
              ALLOCATE(AdjacentDomain%LOCAL_GHOST_FURTHER_INDICES(1),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate LOCAL_GHOST_FURTHER_INDICES.",ERR,ERROR,*999)
              AdjacentDomain%NUMBER_OF_FURTHER_LINKED_GHOSTS = 1
              AdjacentDomain%LOCAL_GHOST_FURTHER_INDICES(1) = GhostElementLocalNo
              
              ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS = ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS+1
              
              IF (DIAGNOSTICS2) THEN

                CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"*reshape to size " // &
                  & TRIM(NumberToVString(ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS,"*",ERR,ERROR)) &
                  & ,ERR,ERROR,*999)
                !PRINT *, MyComputationalNodeNumber, ": *reshape to size ",ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS 
              ENDIF
              
              ALLOCATE(AdjacentDomains(ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate AdjacentDomains.",ERR,ERROR,*999)
              
              ! if there is a previous entry in ELEMENTS_MAPPING%ADJACENT_DOMAINS, save them
              IF (ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS > 1) THEN
                AdjacentDomains(1:ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS-1)=ELEMENTS_MAPPING%ADJACENT_DOMAINS
              ENDIF
              AdjacentDomains(ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS)=AdjacentDomain
              
              CALL MOVE_ALLOC(AdjacentDomains,ELEMENTS_MAPPING%ADJACENT_DOMAINS)
          
            ENDIF ! .NOT.AdjacentDomainAlreadyStored
          ELSE
            IF (DIAGNOSTICS2) THEN

              CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"*element not found on rank " // &
                & TRIM(NumberToVString(ComputationalNodeNo,"*",ERR,ERROR)) &
                & ,ERR,ERROR,*999)
              !PRINT *, MyComputationalNodeNumber, ": *element not found on rank ", ComputationalNodeNo 
            ENDIF
          ENDIF ! element was found
            
          IF (DIAGNOSTICS2) THEN

            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"*FoundElementsCount: " // &
              & TRIM(NumberToVString(FoundElementsCount,"*",ERR,ERROR)) &
              & ,ERR,ERROR,*999)
            !PRINT *, MyComputationalNodeNumber, ": *FoundElementsCount: ", FoundElementsCount 
          ENDIF
          IF (FoundElementsCount >= 6) THEN
            EXIT
          ENDIF
        ENDDO ! sign
        IF (FoundElementsCount >= 6) EXIT
      ENDDO ! I
    ENDDO ! GhostElementLocalNo
        
    ! synchronize RMA
    CALL MPI_WIN_FENCE(MPI_MODE_NOSUCCEED, MpiMemoryWindow, MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_WIN_FENCE",MPI_IERROR,ERR,ERROR,*999)
          
    ! deallocate window 
    CALL MPI_WIN_FREE(MpiMemoryWindow, MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_WIN_FREE",MPI_IERROR,ERR,ERROR,*999)
             
             
           
    EXITS("DECOMPOSITION_ADJACENT_DOMAINS_CALCULATE") 
    RETURN
999 ERRORSEXITS("DECOMPOSITION_ADJACENT_DOMAINS_CALCULATE",ERR,ERROR) 
    RETURN 1
  END SUBROUTINE DECOMPOSITION_ADJACENT_DOMAINS_CALCULATE 
  
  !
  !================================================================================================================================
  !
  !>Checks if a global element is contained in the local domain 
  FUNCTION DECOMPOSITION_ELEMENT_DOMAIN_IS_LOCAL(DECOMPOSITION,GLOBAL_ELEMENT_NO,ERR,ERROR) 
  
    !Argument variables

    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition
    INTEGER(INTG), INTENT(IN) :: GLOBAL_ELEMENT_NO !<The global element number to determine if it is on the local domain 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Function variable
    LOGICAL :: DECOMPOSITION_ELEMENT_DOMAIN_IS_LOCAL 

    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: ElementLocalNo, SearchElementGlobalNo, LowerBoundElementLocalNo, UpperBoundElementLocalNo, &
      & PreviousElementLocalNo 

    ENTERS("DECOMPOSITION_ELEMENT_DOMAIN_IS_LOCAL",ERR,ERROR,*999) 


    DECOMPOSITION_ELEMENT_DOMAIN_IS_LOCAL = .FALSE.
    
    IF(ASSOCIATED(DECOMPOSITION)) THEN 
      IF(ASSOCIATED(DECOMPOSITION%ELEMENTS_MAPPING)) THEN 
        DECOMPOSITION_ELEMENT_DOMAIN_IS_LOCAL = SORTED_ARRAY_CONTAINS_ELEMENT(DECOMPOSITION%ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP, &
          & GLOBAL_ELEMENT_NO,ERR,ERROR) 
      ELSE
        CALL FlagError("Decomposition elements mapping is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_ELEMENT_DOMAIN_IS_LOCAL") 
    RETURN

999 ERRORSEXITS("DECOMPOSITION_ELEMENT_DOMAIN_IS_LOCAL",ERR,ERROR)
    RETURN
  END FUNCTION DECOMPOSITION_ELEMENT_DOMAIN_IS_LOCAL 
     
  
  !
  !================================================================================================================================
  !

  !>Determines if a sorted integer array contains a particular element, if yes, returns the index where that element was found
  FUNCTION SORTED_ARRAY_CONTAINS_ELEMENT_INDEX(ARRAY,ELEMENT,FOUND_INDEX,ERR,ERROR) 
    !Argument variables
    INTEGER(INTG), DIMENSION(:) :: ARRAY !<The array in which to search for the element, needs to be sorted
    INTEGER(INTG), INTENT(IN) :: ELEMENT !<The element to search for
    INTEGER(INTG), INTENT(OUT) :: FOUND_INDEX !< the search index of the array, where the element was found 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    LOGICAL :: SORTED_ARRAY_CONTAINS_ELEMENT_INDEX 
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: SearchIndex, ValueAtSearchIndex, LowerBoundSearchIndex, UpperBoundSearchIndex, &
      & PreviousSearchIndex

    ENTERS("SORTED_ARRAY_CONTAINS_ELEMENT_INDEX",ERR,ERROR,*999) 

    SORTED_ARRAY_CONTAINS_ELEMENT_INDEX = .FALSE. 

    IF (SIZE(ARRAY) == 0) RETURN
    
    SearchIndex = (LBOUND(ARRAY,1)+UBOUND(ARRAY,1))/2
    ! bound values are first indices where the elements can't be
    LowerBoundSearchIndex = LBOUND(ARRAY,1)-1
    UpperBoundSearchIndex = UBOUND(ARRAY,1)+1
    PreviousSearchIndex = -1
    
    !PRINT *, "find ",ELEMENT," in array: ",ARRAY,LBOUND(ARRAY,1),UBOUND(ARRAY,1)
    
    DO WHILE(LowerBoundSearchIndex < UpperBoundSearchIndex)
      ValueAtSearchIndex = ARRAY(SearchIndex)
      
      !PRINT *, "  [",LowerBoundSearchIndex,SearchIndex,UpperBoundSearchIndex,"]"
      !PRINT *, "  current pointer ValueAtSearchIndex=",ValueAtSearchIndex,", to search ELEMENT=", &
      !  & ELEMENT
    
      IF (ValueAtSearchIndex == ELEMENT) THEN
        ! global element number was found among local elements

        SORTED_ARRAY_CONTAINS_ELEMENT_INDEX = .TRUE.
        FOUND_INDEX = SearchIndex 
        !PRINT *, "found"
        EXITS("SORTED_ARRAY_CONTAINS_ELEMENT_INDEX") 
        RETURN
      
      ELSEIF (ValueAtSearchIndex > ELEMENT) THEN
        ! search position was too high
        UpperBoundSearchIndex = SearchIndex
        SearchIndex = NINT((SearchIndex+LowerBoundSearchIndex)/2.)
        
      ELSEIF (ValueAtSearchIndex < ELEMENT) THEN
        ! search position was too low
        LowerBoundSearchIndex = SearchIndex
        SearchIndex = NINT((SearchIndex+UpperBoundSearchIndex)/2.)
      ENDIF
      
      ! If search position did not change or new position is invalid, terminate search
      IF (SearchIndex == PreviousSearchIndex .OR. SearchIndex <= LBOUND(ARRAY,1)-1 .OR. SearchIndex >= UBOUND(ARRAY,1)+1) EXIT
      
      PreviousSearchIndex = SearchIndex
    ENDDO
    
    EXITS("SORTED_ARRAY_CONTAINS_ELEMENT_INDEX")
    RETURN
999 ERRORSEXITS("SORTED_ARRAY_CONTAINS_ELEMENT_INDEX",ERR,ERROR)
    RETURN
  END FUNCTION SORTED_ARRAY_CONTAINS_ELEMENT_INDEX
  
  !
  !================================================================================================================================
  !

  !>Determines if a sorted integer array contains a particular element
  FUNCTION SORTED_ARRAY_CONTAINS_ELEMENT(ARRAY,ELEMENT,ERR,ERROR)
    !Argument variables
    INTEGER(INTG), DIMENSION(:) :: ARRAY !<The array in which to search for the element, needs to be sorted
    INTEGER(INTG), INTENT(IN) :: ELEMENT !<The element to search for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    LOGICAL :: SORTED_ARRAY_CONTAINS_ELEMENT
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: FoundIndex

    ENTERS("SORTED_ARRAY_CONTAINS_ELEMENT",ERR,ERROR,*999)

    SORTED_ARRAY_CONTAINS_ELEMENT = SORTED_ARRAY_CONTAINS_ELEMENT_INDEX(ARRAY,ELEMENT,FoundIndex,ERR,ERROR)
    
    EXITS("SORTED_ARRAY_CONTAINS_ELEMENT")
    RETURN
999 ERRORSEXITS("SORTED_ARRAY_CONTAINS_ELEMENT",ERR,ERROR)
    RETURN
  END FUNCTION SORTED_ARRAY_CONTAINS_ELEMENT

  !
  !================================================================================================================================
  !

  !>Gets the domain for a given element in a decomposition of a mesh. \todo should be able to specify lists of elements. \see OPENCMISS::Iron::cmfe_DecompositionElementDomainGet
  SUBROUTINE DECOMPOSITION_ELEMENT_DOMAIN_GET(DECOMPOSITION,USER_ELEMENT_NUMBER,DOMAIN_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to set the element domain for
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NUMBER !<The user element number to set the domain for.
    INTEGER(INTG), INTENT(OUT) :: DOMAIN_NUMBER !<On return, the domain of the global element.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(MeshComponentTopologyType), POINTER :: MESH_TOPOLOGY
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: GLOBAL_ELEMENT_NUMBER
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE
    TYPE(MeshElementsType), POINTER :: MESH_ELEMENTS


    ENTERS("DECOMPOSITION_ELEMENT_DOMAIN_GET",ERR,ERROR,*999)

    GLOBAL_ELEMENT_NUMBER=0
    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        MESH=>DECOMPOSITION%MESH
        IF(ASSOCIATED(MESH)) THEN
          MESH_TOPOLOGY=>MESH%TOPOLOGY(DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
          IF(ASSOCIATED(MESH_TOPOLOGY)) THEN
            MESH_ELEMENTS=>MESH_TOPOLOGY%ELEMENTS
            IF(ASSOCIATED(MESH_ELEMENTS)) THEN
              NULLIFY(TREE_NODE)
              CALL TREE_SEARCH(MESH_ELEMENTS%ELEMENTS_TREE,USER_ELEMENT_NUMBER,TREE_NODE,ERR,ERROR,*999)
              IF(ASSOCIATED(TREE_NODE)) THEN
                CALL TREE_NODE_VALUE_GET(MESH_ELEMENTS%ELEMENTS_TREE,TREE_NODE,GLOBAL_ELEMENT_NUMBER,ERR,ERROR,*999)
                IF(GLOBAL_ELEMENT_NUMBER>0.AND.GLOBAL_ELEMENT_NUMBER<=MESH_TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS) THEN
                  DOMAIN_NUMBER=DECOMPOSITION%ELEMENT_DOMAIN(GLOBAL_ELEMENT_NUMBER)
                ELSE
                  LOCAL_ERROR="Global element number found "//TRIM(NUMBER_TO_VSTRING(GLOBAL_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                    & " is invalid. The limits are 1 to "// &
                    & TRIM(NUMBER_TO_VSTRING(MESH_TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="Decomposition mesh element corresponding to user element number "// &
                  & TRIM(NumberToVString(USER_ELEMENT_NUMBER,"*",err,error))//" is not found."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Decomposition mesh elements are not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Decomposition mesh topology is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_ELEMENT_DOMAIN_GET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_ELEMENT_DOMAIN_GET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_ELEMENT_DOMAIN_GET
  
  !
  !================================================================================================================================
  !

  !>Sets the domain for a given element in a decomposition of a mesh. \todo move to user number, should be able to specify lists of elements. \see OPENCMISS::Iron::cmfe_DecompositionElementDomainSet 
  SUBROUTINE DECOMPOSITION_ELEMENT_DOMAIN_SET(DECOMPOSITION,GLOBAL_ELEMENT_NUMBER,DOMAIN_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to set the element domain for
    INTEGER(INTG), INTENT(IN) :: GLOBAL_ELEMENT_NUMBER !<The global element number to set the domain for.
    INTEGER(INTG), INTENT(IN) :: DOMAIN_NUMBER !<The domain of the global element.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: number_computational_nodes
    INTEGER(INTG) :: GlobalElementNoDomainPair(2)
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(MeshComponentTopologyType), POINTER :: MESH_TOPOLOGY
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("DECOMPOSITION_ELEMENT_DOMAIN_SET",ERR,ERROR,*999)

!!TODO: interface should specify user element number ???
    
    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        CALL FlagError("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        MESH=>DECOMPOSITION%MESH
        IF(ASSOCIATED(MESH)) THEN
          MESH_TOPOLOGY=>MESH%TOPOLOGY(DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
          IF(ASSOCIATED(MESH_TOPOLOGY)) THEN
            IF(GLOBAL_ELEMENT_NUMBER>0.AND.GLOBAL_ELEMENT_NUMBER<=MESH_TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS) THEN
              CALL WorkGroup_NumberOfGroupNodesGet(decomposition%workGroup,number_computational_nodes,err,error,*999)
              !number_computational_nodes=ComputationalEnvironment_NumberOfNodesGet(ERR,ERROR)
              IF(ERR/=0) GOTO 999
              IF(DOMAIN_NUMBER>=0.AND.DOMAIN_NUMBER<number_computational_nodes) THEN
!              DECOMPOSITION%ELEMENT_DOMAIN(GLOBAL_ELEMENT_NUMBER)=DOMAIN_NUMBER
                IF (USE_OLD_GLOBAL_IMPLEMENTATION) THEN
                  DECOMPOSITION%ELEMENT_DOMAIN(GLOBAL_ELEMENT_NUMBER)=DOMAIN_NUMBER
                ENDIF
                
                IF (USE_NEW_LOCAL_IMPLEMENTATION) THEN
                  GlobalElementNoDomainPair(1) = GLOBAL_ELEMENT_NUMBER
                  GlobalElementNoDomainPair(2) = DOMAIN_NUMBER
                  CALL LIST_ITEM_ADD(DECOMPOSITION%GlobalElementDomain,GlobalElementNoDomainPair,ERR,ERROR,*999)
                ENDIF 
              ELSE
                LOCAL_ERROR="Domain number "//TRIM(NUMBER_TO_VSTRING(DOMAIN_NUMBER,"*",ERR,ERROR))// &
                  & " is invalid. The limits are 0 to "//TRIM(NUMBER_TO_VSTRING(number_computational_nodes,"*",ERR,ERROR))//"."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                & " is invalid. The limits are 1 to "// &
                & TRIM(NUMBER_TO_VSTRING(MESH_TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Decomposition mesh topology is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_ELEMENT_DOMAIN_SET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_ELEMENT_DOMAIN_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_ELEMENT_DOMAIN_SET
  
  !
  !================================================================================================================================
  !

  !>Finalises a decomposition and deallocates all memory.
  SUBROUTINE Decomposition_Finalise(decomposition,err,error,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition !<A pointer to the decomposition to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Decomposition_Finalise",err,error,*999)

    
    IF(ASSOCIATED(decomposition)) THEN
      !Destroy all the decomposition components          
      IF(ALLOCATED(decomposition%ELEMENT_DOMAIN)) DEALLOCATE(decomposition%ELEMENT_DOMAIN)
      CALL DECOMPOSITION_TOPOLOGY_FINALISE(decomposition,err,error,*999)
      CALL DOMAIN_FINALISE(decomposition,err,error,*999)          
      DEALLOCATE(decomposition)
    ENDIF
    
    EXITS("Decomposition_Finalise")
    RETURN
999 ERRORSEXITS("Decomposition_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_Finalise
  
  !
  !================================================================================================================================
  !

  !>Inintialises a decomposition and allocates all memory.
  SUBROUTINE Decomposition_Initialise(decomposition,err,error,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition !<A pointer to the decomposition to initialise. Must not be allocated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Decomposition_Initialise",err,error,*998)

    IF(ASSOCIATED(decomposition)) THEN
      CALL FlagError("Decomposition is already associated.",err,error,*998)
    ELSE
      ALLOCATE(decomposition,STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate decomposition.",err,error,*999)
      decomposition%GLOBAL_NUMBER=0
      decomposition%USER_NUMBER=0
      decomposition%DECOMPOSITION_FINISHED=.FALSE.
      NULLIFY(decomposition%decompositions)
      NULLIFY(decomposition%mesh)
      NULLIFY(decomposition%region)
      NULLIFY(decomposition%INTERFACE)
      decomposition%numberOfDimensions=0
      decomposition%numberOfComponents=0
      decomposition%DECOMPOSITION_TYPE=DECOMPOSITION_ALL_TYPE
      NULLIFY(decomposition%workGroup)
      decomposition%NUMBER_OF_DOMAINS=0
      decomposition%NUMBER_OF_EDGES_CUT=0
      decomposition%numberOfElements=0
      NULLIFY(decomposition%topology)
      NULLIFY(decomposition%domain)
      decomposition%CALCULATE_LINES=.TRUE. 
      decomposition%CALCULATE_FACES=.TRUE.
    ENDIF
    
    EXITS("Decomposition_Initialise")
    RETURN
999 CALL Decomposition_Finalise(decomposition,dummyErr,dummyError,*998)
998 ERRORSEXITS("Decomposition_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_Initialise
  
  !
  !================================================================================================================================
  !

  !!MERGE: ditto
  
  !>Gets the mesh component number which will be used for the decomposition of a mesh. \see OPENCMISS::Iron::cmfe_DecompositionMeshComponentGet
  SUBROUTINE DECOMPOSITION_MESH_COMPONENT_NUMBER_GET(DECOMPOSITION,MESH_COMPONENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to get the mesh component for
    INTEGER(INTG), INTENT(OUT) :: MESH_COMPONENT_NUMBER !<On return, the mesh component number to get
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_MESH_COMPONENT_NUMBER_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN     
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        IF(ASSOCIATED(DECOMPOSITION%MESH)) THEN
          MESH_COMPONENT_NUMBER=DECOMPOSITION%MESH_COMPONENT_NUMBER
        ELSE
          CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition has been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_MESH_COMPONENT_NUMBER_GET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_MESH_COMPONENT_NUMBER_GET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_MESH_COMPONENT_NUMBER_GET
  
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the mesh component number which will be used for the decomposition of a mesh. \see OPENCMISS::Iron::cmfe_DecompositionMeshComponentSet
  SUBROUTINE DECOMPOSITION_MESH_COMPONENT_NUMBER_SET(DECOMPOSITION,MESH_COMPONENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to set the mesh component for
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT_NUMBER !<The mesh component number to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("DECOMPOSITION_MESH_COMPONENT_NUMBER_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN     
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        CALL FlagError("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(DECOMPOSITION%MESH)) THEN
          IF(MESH_COMPONENT_NUMBER>0.AND.MESH_COMPONENT_NUMBER<=DECOMPOSITION%numberOfComponents) THEN
            DECOMPOSITION%MESH_COMPONENT_NUMBER=MESH_COMPONENT_NUMBER
          ELSE
            LOCAL_ERROR="The specified mesh component number of "//TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT_NUMBER,"*",ERR,ERROR))// &
              & "is invalid. The component number must be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(DECOMPOSITION%numberOfComponents,"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_MESH_COMPONENT_NUMBER_SET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_MESH_COMPONENT_NUMBER_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_MESH_COMPONENT_NUMBER_SET
  
  !
  !================================================================================================================================
  !

  !!MERGE: ditto
  
  !>Gets the number of domains for a decomposition. \see OPENCMISS::Iron::cmfe_DecompositionNumberOfDomainsGet
  SUBROUTINE DECOMPOSITION_NUMBER_OF_DOMAINS_GET(DECOMPOSITION,NUMBER_OF_DOMAINS,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to get the number of domains for.
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_DOMAINS !<On return, the number of domains to get.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_NUMBER_OF_DOMAINS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        CALL FlagError("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        NUMBER_OF_DOMAINS=DECOMPOSITION%NUMBER_OF_DOMAINS
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_NUMBER_OF_DOMAINS_GET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_NUMBER_OF_DOMAINS_GET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_NUMBER_OF_DOMAINS_GET
  

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of domains for a decomposition. \see OPENCMISS::Iron::cmfe_DecompositionNumberOfDomainsSet
  SUBROUTINE DECOMPOSITION_NUMBER_OF_DOMAINS_SET(DECOMPOSITION,NUMBER_OF_DOMAINS,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to set the number of domains for.
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DOMAINS !<The number of domains to set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfComputationalNodes
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("DECOMPOSITION_NUMBER_OF_DOMAINS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        CALL FlagError("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(DECOMPOSITION%DECOMPOSITION_TYPE)
        CASE(DECOMPOSITION_ALL_TYPE)
          IF(NUMBER_OF_DOMAINS==1) THEN
            DECOMPOSITION%NUMBER_OF_DOMAINS=1
          ELSE
            CALL FlagError("Can only have one domain for all decomposition type.",ERR,ERROR,*999)
          ENDIF
        CASE(DECOMPOSITION_CALCULATED_TYPE,DECOMPOSITION_USER_DEFINED_TYPE)
          IF(NUMBER_OF_DOMAINS>=1) THEN
            !wolfye???<=?
            IF(NUMBER_OF_DOMAINS<=DECOMPOSITION%numberOfElements) THEN
              !Get the number of computational nodes
              !numberOfComputationalNodes=ComputationalEnvironment_NumberOfNodesGet(ERR,ERROR)
              CALL WorkGroup_NumberOfGroupNodesGet(decomposition%workGroup,numberOfComputationalNodes, &
                & err,error,*999)
              IF(ERR/=0) GOTO 999
              !!TODO: relax this later
              !IF(NUMBER_OF_DOMAINS==numberOfComputationalNodes) THEN
                DECOMPOSITION%NUMBER_OF_DOMAINS=NUMBER_OF_DOMAINS             
              !ELSE
              !  LOCAL_ERROR="The number of domains ("//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DOMAINSS,"*",ERR,ERROR))// &
              !    & ") is not equal to the number of computational nodes ("// &
              !    & TRIM(NUMBER_TO_VSTRING(numberOfComputationalNodes,"*",ERR,ERROR))//")"
              !  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              !ENDIF
            ELSE
              LOCAL_ERROR="The number of domains ("//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DOMAINS,"*",ERR,ERROR))// &
                & ") must be <= the number of global elements ("// &
                & TRIM(NUMBER_TO_VSTRING(DECOMPOSITION%numberOfElements,"*",ERR,ERROR))//") in the mesh."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
             CALL FlagError("Number of domains must be >= 1.",ERR,ERROR,*999)
           ENDIF
         CASE DEFAULT
          LOCAL_ERROR="Decomposition type "//TRIM(NUMBER_TO_VSTRING(DECOMPOSITION%DECOMPOSITION_TYPE,"*",ERR,ERROR))// &
            & " is not valid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_NUMBER_OF_DOMAINS_SET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_NUMBER_OF_DOMAINS_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_NUMBER_OF_DOMAINS_SET

  !
  !================================================================================================================================
  !



  !>Sets the workgroup to use for a decomposition on a given mesh. \see OpenCMISS::Iron::cmfe_Decomposition_WorkGroupSet
  SUBROUTINE Decomposition_WorkGroupSet(decomposition,workGroup,err,error,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition !<A pointer to the decomposition to set the work group for
    TYPE(WorkGroupType), POINTER :: workGroup !<A pointer to the work group to use for the decomposition
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: decomposition_no
    TYPE(MESH_TYPE), POINTER :: MESH

    ENTERS("Decomposition_WorkGroupSet",err,error,*999)

    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
    IF(decomposition%DECOMPOSITION_FINISHED) CALL FlagError("Decomposition has already been finished.",err,error,*999)
    IF(.NOT.ASSOCIATED(workGroup)) CALL FlagError("Work group is not associated.",err,error,*999)
    IF(.NOT.workGroup%workGroupFinished) CALL FlagError("Work group has not been finished.",err,error,*999)

    decomposition%workGroup=>workGroup
       
    EXITS("Decomposition_WorkGroupSet")
    RETURN
999 ERRORSEXITS("Decomposition_WorkGroupSet",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_WorkGroupSet

  !
  !================================================================================================================================
  !

  !>Calculates the topology for a decomposition.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_CALCULATE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to finish creating.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: meshComponentNumber

    ENTERS("DECOMPOSITION_TOPOLOGY_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION%TOPOLOGY)) THEN
      !Calculate the elements topology
      CALL DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      !Calculate the line topology
      IF(DECOMPOSITION%CALCULATE_LINES)THEN
        CALL DECOMPOSITION_TOPOLOGY_LINES_CALCULATE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      ENDIF
      !Calculate the face topology
      IF(DECOMPOSITION%CALCULATE_FACES) THEN
        CALL DECOMPOSITION_TOPOLOGY_FACES_CALCULATE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      ENDIF
      meshComponentNumber=DECOMPOSITION%MESH_COMPONENT_NUMBER
      IF(ALLOCATED(DECOMPOSITION%MESH%TOPOLOGY(meshComponentNumber)%PTR%dataPoints%dataPoints)) THEN
          CALL DecompositionTopology_DataPointsCalculate(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
        ENDIF   
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TOPOLOGY_CALCULATE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_CALCULATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_CALCULATE
  

  !
  !================================================================================================================================
  !

  !>Calculates the decomposition element topology.
  SUBROUTINE DecompositionTopology_DataPointsCalculate(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to calculate the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: localElement,globalElement,dataPointIdx,localData,meshComponentNumber
    INTEGER(INTG) :: INSERT_STATUS,MPI_IERROR,NUMBER_OF_COMPUTATIONAL_NODES,myComputationalNodeNumber,NUMBER_OF_GHOST_DATA, &
      & NUMBER_OF_LOCAL_DATA, groupCommunicator
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: decompositionElements
    TYPE(DecompositionDataPointsType), POINTER :: decompositionData
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: elementsMapping
    TYPE(MeshDataPointsType), POINTER :: meshData

    ENTERS("DecompositionTopology_DataPointsCalculate",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      decompositionData=>TOPOLOGY%dataPoints
      IF(ASSOCIATED(decompositionData)) THEN
        decomposition=>decompositionData%DECOMPOSITION
        IF(ASSOCIATED(decomposition)) THEN
         CALL WorkGroup_GroupCommunicatorGet(decomposition%workGroup,groupCommunicator,err,error,*999)
         decompositionElements=>TOPOLOGY%ELEMENTS
         IF(ASSOCIATED(decompositionElements)) THEN
           elementsMapping=>decomposition%DOMAIN(decomposition%MESH_COMPONENT_NUMBER)%PTR%MAPPINGS%ELEMENTS
           IF(ASSOCIATED(elementsMapping)) THEN
              meshComponentNumber=decomposition%MESH_COMPONENT_NUMBER
              meshData=>decomposition%MESH%TOPOLOGY(meshComponentNumber)%PTR%dataPoints
              IF(ASSOCIATED(meshData)) THEN
                CALL WorkGroup_NumberOfGroupNodesGet(decomposition%workGroup,NUMBER_OF_COMPUTATIONAL_NODES,err,error,*999)
                CALL WorkGroup_GroupNodeNumberGet(decomposition%workGroup,myComputationalNodeNumber,err,error,*999)
                !NUMBER_OF_COMPUTATIONAL_NODES=ComputationalEnvironment_NumberOfNodesGet(ERR,ERROR)
                !IF(ERR/=0) GOTO 999
                !myComputationalNodeNumber=ComputationalEnvironment_NodeNumberGet(ERR,ERROR)
                !IF(ERR/=0) GOTO 999
                ALLOCATE(decompositionData%numberOfDomainLocal(0:NUMBER_OF_COMPUTATIONAL_NODES-1),STAT=ERR)
                ALLOCATE(decompositionData%numberOfDomainGhost(0:NUMBER_OF_COMPUTATIONAL_NODES-1),STAT=ERR)
                ALLOCATE(decompositionData%numberOfElementDataPoints(decompositionElements%NUMBER_OF_GLOBAL_ELEMENTS),STAT=ERR)
                ALLOCATE(decompositionData%elementDataPoint(decompositionElements%TOTAL_NUMBER_OF_ELEMENTS),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate decomposition element data points.",ERR,ERROR,*999)
                CALL TREE_CREATE_START(decompositionData%dataPointsTree,ERR,ERROR,*999)
                CALL TREE_INSERT_TYPE_SET(decompositionData%dataPointsTree,TREE_NO_DUPLICATES_ALLOWED,ERR,ERROR,*999)
                CALL TREE_CREATE_FINISH(decompositionData%dataPointsTree,ERR,ERROR,*999)
                decompositionData%numberOfGlobalDataPoints=meshData%totalNumberOfProjectedData 
                DO globalElement=1,decompositionElements%NUMBER_OF_GLOBAL_ELEMENTS
                  decompositionData%numberOfElementDataPoints(globalElement)= &
                    & meshData%elementDataPoint(globalElement)%numberOfProjectedData
                ENDDO !globalElement
                localData=0;
                DO localElement=1,decompositionElements%TOTAL_NUMBER_OF_ELEMENTS
                  globalElement=decompositionElements%ELEMENTS(localElement)%GLOBAL_NUMBER
                  decompositionData%elementDataPoint(localElement)%numberOfProjectedData= &
                    & meshData%elementDataPoint(globalElement)%numberOfProjectedData
                  decompositionData%elementDataPoint(localElement)%globalElementNumber=globalElement
                  IF(localElement<elementsMapping%GHOST_START) THEN
                    decompositionData%numberOfDataPoints=decompositionData%numberOfDataPoints+ &
                      & decompositionData%elementDataPoint(localElement)%numberOfProjectedData
                  ENDIF               
                  decompositionData%totalNumberOfDataPoints=decompositionData%totalNumberOfDataPoints+ &
                    & decompositionData%elementDataPoint(localElement)%numberOfProjectedData
                  ALLOCATE(decompositionData%elementDataPoint(localElement)%dataIndices(decompositionData% &
                    & elementDataPoint(localElement)%numberOfProjectedData),STAT=ERR)
                  DO dataPointIdx=1,decompositionData%elementDataPoint(localElement)%numberOfProjectedData
                    decompositionData%elementDataPoint(localElement)%dataIndices(dataPointIdx)%userNumber= &
                      & meshData%elementDataPoint(globalElement)%dataIndices(dataPointIdx)%userNumber
                    decompositionData%elementDataPoint(localElement)%dataIndices(dataPointIdx)%globalNumber= &
                      & meshData%elementDataPoint(globalElement)%dataIndices(dataPointIdx)%globalNumber
                    localData=localData+1
                    decompositionData%elementDataPoint(localElement)%dataIndices(dataPointIdx)%localNumber=localData
                    CALL TREE_ITEM_INSERT(decompositionData%dataPointsTree,decompositionData% &
                      & elementDataPoint(localElement)%dataIndices(dataPointIdx)%userNumber,localData, &
                      & INSERT_STATUS,ERR,ERROR,*999)
                  ENDDO !dataPointIdx
                ENDDO !localElement   
                !Calculate number of ghost data points on the current computational domain
                NUMBER_OF_LOCAL_DATA=decompositionData%numberOfDataPoints
                NUMBER_OF_GHOST_DATA=decompositionData%totalNumberOfDataPoints-decompositionData%numberOfDataPoints
                !Gather number of local data points on all computational nodes
                CALL MPI_ALLGATHER(NUMBER_OF_LOCAL_DATA,1,MPI_INTEGER,decompositionData% &
                  & numberOfDomainLocal,1,MPI_INTEGER,groupCommunicator,MPI_IERROR)
                CALL MPI_ERROR_CHECK("MPI_ALLGATHER",MPI_IERROR,ERR,ERROR,*999)
                !Gather number of ghost data points on all computational nodes
                CALL MPI_ALLGATHER(NUMBER_OF_GHOST_DATA,1,MPI_INTEGER,decompositionData% &
                  & numberOfDomainGhost,1,MPI_INTEGER,groupCommunicator,MPI_IERROR)
                CALL MPI_ERROR_CHECK("MPI_ALLGATHER",MPI_IERROR,ERR,ERROR,*999)
              ELSE
                CALL FlagError("Mesh data points topology is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Element mapping  is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Decomposition elements topology is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition data points topology is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DecompositionTopology_DataPointsCalculate")
    RETURN
999 ERRORSEXITS("DecompositionTopology_DataPointsCalculate",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DecompositionTopology_DataPointsCalculate

  !
  !================================================================================================================================
  !

  !>Calculates the decomposition element topology for a data projection (for data projections on fields).
  SUBROUTINE DecompositionTopology_DataProjectionCalculate(decompositionTopology,err,error,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology !<A pointer to the decomposition topology to calculate the elements for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("DecompositionTopology_DataProjectionCalculate",err,error,*999)

    IF(ASSOCIATED(decompositionTopology)) THEN
      CALL DecompositionTopology_DataPointsInitialise(decompositionTopology,err,error,*999)
      CALL DecompositionTopology_DataPointsCalculate(decompositionTopology,err,error,*999)
    ELSE
      CALL FlagError("Decomposition topology is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DecompositionTopology_DataProjectionCalculate")
    RETURN
999 ERRORS("DecompositionTopology_DataProjectionCalculate",err,error)
    EXITS("DecompositionTopology_DataProjectionCalculate")
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_DataProjectionCalculate

  !
  !================================================================================================================================
  !

  !>Gets the local data point number for data points projected on an element
  SUBROUTINE DecompositionTopology_ElementDataPointLocalNumberGet(decompositionTopology,elementNumber,dataPointIndex, &
       & dataPointLocalNumber,err,error,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology !<A pointer to the decomposition topology to calculate the elements for
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to get the data point for
    INTEGER(INTG), INTENT(IN) :: dataPointIndex !<The data point index to get the number for
    INTEGER(INTG), INTENT(OUT) :: dataPointLocalNumber !<The data point local number to return
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(DecompositionDataPointsType), POINTER :: decompositionData
    INTEGER(INTG) :: numberOfDataPoints
    TYPE(VARYING_STRING) :: localError

    ENTERS("DecompositionTopology_ElementDataPointLocalNumberGet",err,error,*999)

    IF(ASSOCIATED(decompositionTopology)) THEN
      decompositionData=>decompositionTopology%dataPoints
      IF(ASSOCIATED(decompositionData)) THEN
        numberOfDataPoints = decompositionData%elementDataPoint(elementNumber)%numberOfProjectedData
        IF(dataPointIndex > 0 .AND. dataPointIndex <= numberOfDataPoints) THEN
          dataPointLocalNumber = decompositionData%elementDataPoint(elementNumber)%dataIndices(dataPointIndex)%localNumber
        ELSE
          localError="Element data point index "//TRIM(NUMBER_TO_VSTRING(dataPointIndex,"*",ERR,ERROR))// &
           & " out of range for element "//TRIM(NUMBER_TO_VSTRING(elementNumber,"*",ERR,ERROR))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition topology data points are not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition topology is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DecompositionTopology_ElementDataPointLocalNumberGet")
    RETURN
999 ERRORS("DecompositionTopology_ElementDataPointLocalNumberGet",err,error)
    EXITS("DecompositionTopology_ElementDataPointLocalNumberGet")
    RETURN 1
  END SUBROUTINE DecompositionTopology_ElementDataPointLocalNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the user (global) data point number for data points projected on an element
  SUBROUTINE DecompositionTopology_ElementDataPointUserNumberGet(decompositionTopology,userElementNumber,dataPointIndex, &
       & dataPointUserNumber,err,error,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology !<A pointer to the decomposition topology to calculate the elements for
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The element number to get the data point for
    INTEGER(INTG), INTENT(IN) :: dataPointIndex !<The data point index to get the number for
    INTEGER(INTG), INTENT(OUT) :: dataPointUserNumber !<The data point user number to return
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(DecompositionDataPointsType), POINTER :: decompositionData
    INTEGER(INTG) :: numberOfDataPoints,decompositionLocalElementNumber
    LOGICAL :: ghostElement,userElementExists
    TYPE(VARYING_STRING) :: localError

    ENTERS("DecompositionTopology_ElementDataPointUserNumberGet",err,error,*999)

    IF(ASSOCIATED(decompositionTopology)) THEN
      decompositionData=>decompositionTopology%dataPoints
      IF(ASSOCIATED(decompositionData)) THEN
        CALL DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS(decompositionTopology,userElementNumber, &
          & userElementExists,decompositionLocalElementNumber,ghostElement,err,error,*999)
        IF(userElementExists) THEN
          IF(ghostElement) THEN
            localError="Cannot update by data point for user element "// &
              & TRIM(NUMBER_TO_VSTRING(userElementNumber,"*",err,error))//" as it is a ghost element."
            CALL FlagError(localError,err,error,*999)
          ELSE
            numberOfDataPoints = decompositionData%elementDataPoint(decompositionLocalElementNumber)%numberOfProjectedData
            IF(dataPointIndex > 0 .AND. dataPointIndex <= numberOfDataPoints) THEN
              dataPointUserNumber = decompositionData%elementDataPoint(decompositionLocalElementNumber)% &
                & dataIndices(dataPointIndex)%userNumber
            ELSE
              localError="Element data point index "//TRIM(NUMBER_TO_VSTRING(dataPointIndex,"*",ERR,ERROR))// &
               & " out of range for element "//TRIM(NUMBER_TO_VSTRING(userElementNumber,"*",ERR,ERROR))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ENDIF
        ELSE
          localError="The specified user element number of "// &
            & TRIM(NUMBER_TO_VSTRING(userElementNumber,"*",err,error))// &
            & " does not exist."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition topology data points are not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition topology is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DecompositionTopology_ElementDataPointUserNumberGet")
    RETURN
999 ERRORS("DecompositionTopology_ElementDataPointUserNumberGet",err,error)
    EXITS("DecompositionTopology_ElementDataPointUserNumberGet")
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_ElementDataPointUserNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the number of data points projected on an element
  SUBROUTINE DecompositionTopology_NumberOfElementDataPointsGet(decompositionTopology,userElementNumber, &
       & numberOfDataPoints,err,error,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology !<A pointer to the decomposition topology to calculate the elements for
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The element number to get the data point for
    INTEGER(INTG), INTENT(OUT) :: numberOfDataPoints !<The data point local number to return
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(DecompositionDataPointsType), POINTER :: decompositionData
    INTEGER(INTG) :: decompositionLocalElementNumber
    LOGICAL :: ghostElement,userElementExists
    TYPE(VARYING_STRING) :: localError

    ENTERS("DecompositionTopology_NumberOfElementDataPointsGet",err,error,*999)

    IF(ASSOCIATED(decompositionTopology)) THEN
      decompositionData=>decompositionTopology%dataPoints
      IF(ASSOCIATED(decompositionData)) THEN
        CALL DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS(decompositionTopology,userElementNumber, &
          & userElementExists,decompositionLocalElementNumber,ghostElement,err,error,*999)
        IF(userElementExists) THEN
          IF(ghostElement) THEN
            localError="Cannot update by data point for user element "// &
              & TRIM(NUMBER_TO_VSTRING(userElementNumber,"*",err,error))//" as it is a ghost element."
            CALL FlagError(localError,err,error,*999)
          ELSE
            numberOfDataPoints = decompositionData%elementDataPoint(decompositionLocalElementNumber)%numberOfProjectedData
          ENDIF
        ELSE
          localError="The specified user element number of "// &
            & TRIM(NUMBER_TO_VSTRING(userElementNumber,"*",err,error))// &
            & " does not exist."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition topology data points are not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition topology is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DecompositionTopology_NumberOfElementDataPointsGet")
    RETURN
999 ERRORS("DecompositionTopology_NumberOfElementDataPointsGet",err,error)
    EXITS("DecompositionTopology_NumberOfElementDataPointsGet")
    RETURN 1
  END SUBROUTINE DecompositionTopology_NumberOfElementDataPointsGet
  
  !
  !================================================================================================================================
  !

  !>Checks that a user element number exists in a decomposition. 
  SUBROUTINE DecompositionTopology_DataPointCheckExists(decompositionTopology,userDataPointNumber,userDataPointExists, &
        & decompositionLocalDataPointNumber,ghostDataPoint,err,error,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology !<A pointer to the decomposition topology to check the data point exists on
    INTEGER(INTG), INTENT(IN) :: userDataPointNumber !<The user data point number to check if it exists
    LOGICAL, INTENT(OUT) :: userDataPointExists !<On exit, is .TRUE. if the data point user number exists in the decomposition topology, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: decompositionLocalDataPointNumber !<On exit, if the data point exists the local number corresponding to the user data point number. If the data point does not exist then local number will be 0.
    LOGICAL, INTENT(OUT) :: ghostDataPoint !<On exit, is .TRUE. if the local data point (if it exists) is a ghost data point, .FALSE. if not.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DecompositionDataPointsType), POINTER :: decompositionData
    TYPE(TREE_NODE_TYPE), POINTER :: treeNode
    
    ENTERS("DecompositionTopology_DataPointCheckExists",ERR,error,*999)

    userDataPointExists=.FALSE.
    decompositionLocalDataPointNumber=0
    ghostDataPoint=.FALSE.
    IF(ASSOCIATED(decompositionTopology)) THEN
      decompositionData=>decompositionTopology%dataPoints
      IF(ASSOCIATED(decompositionData)) THEN
        NULLIFY(treeNode)
        CALL TREE_SEARCH(decompositionData%dataPointsTree,userDataPointNumber,treeNode,err,error,*999)
        IF(ASSOCIATED(treeNode)) THEN
          CALL TREE_NODE_VALUE_GET(decompositionData%dataPointsTree,treeNode,decompositionLocalDataPointNumber,err,error,*999)
          userDataPointExists=.TRUE.
          ghostDataPoint=decompositionLocalDataPointNumber>decompositionData%numberOfDataPoints
        ENDIF
      ELSE
        CALL FlagError("Decomposition data point topology is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition topology is not associated.",err,error,*999)
    ENDIF
    
    EXITS("DecompositionTopology_DataPointCheckExists")
    RETURN
999 ERRORS("DecompositionTopology_DataPointCheckExists",err,error)
    EXITS("DecompositionTopology_DataPointCheckExists")
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_DataPointCheckExists
  
  !
  !================================================================================================================================
  !

  !>Checks that a user element number exists in a decomposition. 
  SUBROUTINE DecompositionTopology_ElementCheckExists(decompositionTopology,userElementNumber,elementExists,localElementNumber, &
    & ghostElement,err,error,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology !<A pointer to the decomposition topology to check the element exists on
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number to check if it exists
    LOGICAL, INTENT(OUT) :: elementExists!<On exit, is .TRUE. if the element user number exists in the decomposition topology, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: localElementNumber !<On exit, if the element exists the local number corresponding to the user element number. If the element does not exist then local number will be 0.
    LOGICAL, INTENT(OUT) :: ghostElement !<On exit, is .TRUE. if the local element (if it exists) is a ghost element, .FALSE. if not.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: decompositionElements
    TYPE(TREE_NODE_TYPE), POINTER :: treeNode
    
    ENTERS("DecompositionTopology_ElementCheckExists",err,error,*999)

    elementExists=.FALSE.
    localElementNumber=0
    ghostElement=.FALSE.
    IF(.NOT.ASSOCIATED(decompositionTopology)) CALL FlagError("Decomposition topology is not associated.",err,error,*999)
    
    NULLIFY(decompositionElements)
    CALL DecompositionTopology_ElementsGet(decompositionTopology,decompositionElements,err,error,*999)
    NULLIFY(treeNode)
    CALL Tree_Search(decompositionElements%ELEMENTS_TREE,userElementNumber,treeNode,err,error,*999)
    IF(ASSOCIATED(treeNode)) THEN
      CALL Tree_NodeValueGet(decompositionElements%ELEMENTS_TREE,treeNode,localElementNumber,err,error,*999)
      elementExists=.TRUE.
      ghostElement=localElementNumber>decompositionElements%NUMBER_OF_ELEMENTS
    ENDIF
    
    EXITS("DecompositionTopology_ElementCheckExists")
    RETURN
999 ERRORSEXITS("DecompositionTopology_ElementCheckExists",err,error)
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_ElementCheckExists

  !
  !================================================================================================================================
  !

!!TODO: Should this be decompositionElements???
  !>Gets a local element number that corresponds to a user element number from a decomposition. An error will be raised if the user element number does not exist.
  SUBROUTINE DecompositionTopology_ElementGet(decompositionTopology,userElementNumber,localElementNumber,ghostElement,err,error,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology !<A pointer to the decomposition topology to get the element on
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number to get
    INTEGER(INTG), INTENT(OUT) :: localElementNumber !<On exit, the local number corresponding to the user element number.
    LOGICAL, INTENT(OUT) :: ghostElement !<On exit, is .TRUE. if the local element is a ghost element, .FALSE. if not.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    LOGICAL :: elementExists
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DecompositionTopology_ElementGet",err,error,*999)

    CALL DecompositionTopology_ElementCheckExists(decompositionTopology,userElementNumber,elementExists,localElementNumber, &
      & ghostElement,err,error,*999)
    IF(.NOT.elementExists) THEN
      decomposition=>decompositionTopology%decomposition
      IF(ASSOCIATED(decomposition)) THEN
        localError="The user element number "//TRIM(NumberToVString(userElementNumber,"*",err,error))// &
          & " does not exist in decomposition number "//TRIM(NumberToVString(decomposition%USER_NUMBER,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ELSE
        localError="The user element number "//TRIM(NumberToVString(userElementNumber,"*",err,error))//" does not exist."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDIF
   
    EXITS("DecompositionTopology_ElementGet")
    RETURN
999 ERRORSEXITS("DecompositionTopology_ElementGet",err,error)
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_ElementGet

  !
  !================================================================================================================================
  !

  !>Finalises the given decomposition topology element.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE(ELEMENT,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_ELEMENT_TYPE) :: ELEMENT !<The decomposition element to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nic !\todo add comment

    ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(ELEMENT%ADJACENT_ELEMENTS)) THEN
      DO nic=LBOUND(ELEMENT%ADJACENT_ELEMENTS,1),UBOUND(ELEMENT%ADJACENT_ELEMENTS,1)
        CALL DECOMPOSITION_ADJACENT_ELEMENT_FINALISE(ELEMENT%ADJACENT_ELEMENTS(nic),ERR,ERROR,*999)
      ENDDO !nic
      DEALLOCATE(ELEMENT%ADJACENT_ELEMENTS)
    ENDIF
    IF(ALLOCATED(ELEMENT%ELEMENT_LINES)) DEALLOCATE(ELEMENT%ELEMENT_LINES)
    IF(ALLOCATED(ELEMENT%ELEMENT_FACES)) DEALLOCATE(ELEMENT%ELEMENT_FACES)
 
    EXITS("DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the given decomposition topology element.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENT_INITIALISE(ELEMENT,ERR,ERROR,*)
 
    !Argument variables
    TYPE(DECOMPOSITION_ELEMENT_TYPE) :: ELEMENT !<The decomposition element to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENT_INITIALISE",ERR,ERROR,*999)

    ELEMENT%USER_NUMBER=0
    ELEMENT%LOCAL_NUMBER=0
    ELEMENT%GLOBAL_NUMBER=0
    ELEMENT%BOUNDARY_ELEMENT=.FALSE.
  
    EXITS("DECOMPOSITION_TOPOLOGY_ELEMENT_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_ELEMENT_INITIALISE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Calculates the element numbers adjacent to an element in a decomposition topology.
  SUBROUTINE DecompositionTopology_ElementAdjacentElementCalculate(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to calculate the adjacent elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: j,ne,ne1,nep1,ni,nic,nn,nn1,nn2,nn3,np,np1,DUMMY_ERR,FACE_XI(2),FACE_XIC(3),NODE_POSITION_INDEX(4)
    INTEGER(INTG) :: xi_direction,direction_index,xi_dir_check,xi_dir_search,NUMBER_NODE_MATCHES
    INTEGER(INTG) :: candidate_idx,face_node_idx,node_idx,surrounding_el_idx,candidate_el,idx
    INTEGER(INTG) :: NUMBER_SURROUNDING,NUMBER_OF_NODES_XIC(4),numberSurroundingElements
    INTEGER(INTG), ALLOCATABLE :: NODE_MATCHES(:),ADJACENT_ELEMENTS(:), surroundingElements(:)
    LOGICAL :: XI_COLLAPSED,FACE_COLLAPSED(-3:3),SUBSET
    TYPE(LIST_TYPE), POINTER :: NODE_MATCH_LIST, surroundingElementsList
    TYPE(LIST_PTR_TYPE) :: ADJACENT_ELEMENTS_LIST(-4:4)
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: DECOMPOSITION_ELEMENTS
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: DOMAIN_TOPOLOGY
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    NULLIFY(NODE_MATCH_LIST)
    DO nic=-4,4
      NULLIFY(ADJACENT_ELEMENTS_LIST(nic)%PTR)
    ENDDO !nic
    
    ENTERS("DecompositionTopology_ElementAdjacentElementCalculate",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
      IF(ASSOCIATED(DECOMPOSITION)) THEN
        DECOMPOSITION_ELEMENTS=>TOPOLOGY%ELEMENTS
        IF(ASSOCIATED(DECOMPOSITION_ELEMENTS)) THEN
          DOMAIN=>DECOMPOSITION%DOMAIN(DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
          IF(ASSOCIATED(DOMAIN)) THEN
            DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
            IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
              DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
              IF(ASSOCIATED(DOMAIN_NODES)) THEN
                DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS
                IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                  !Loop over the elements in the decomposition
                  DO ne=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                    BASIS=>DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS
                    !Create a list for every xi direction (plus and minus)
                    DO nic=-BASIS%NUMBER_OF_XI_COORDINATES,BASIS%NUMBER_OF_XI_COORDINATES
                      NULLIFY(ADJACENT_ELEMENTS_LIST(nic)%PTR)
                      CALL LIST_CREATE_START(ADJACENT_ELEMENTS_LIST(nic)%PTR,ERR,ERROR,*999)
                      CALL LIST_DATA_TYPE_SET(ADJACENT_ELEMENTS_LIST(nic)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                      CALL LIST_INITIAL_SIZE_SET(ADJACENT_ELEMENTS_LIST(nic)%PTR,5,ERR,ERROR,*999)
                      CALL LIST_CREATE_FINISH(ADJACENT_ELEMENTS_LIST(nic)%PTR,ERR,ERROR,*999)
                    ENDDO !nic
                    NUMBER_OF_NODES_XIC=1
                    NUMBER_OF_NODES_XIC(1:BASIS%NUMBER_OF_XI_COORDINATES)= &
                      & BASIS%NUMBER_OF_NODES_XIC(1:BASIS%NUMBER_OF_XI_COORDINATES)
                    !Place the current element in the surrounding list
                    CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(0)%PTR,DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%LOCAL_NUMBER, &
                      & ERR,ERROR,*999)
                    SELECT CASE(BASIS%TYPE)
                    CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
!!TODO: Calculate this and set it as part of the basis type
                      !Determine the collapsed "faces" if any
                      NODE_POSITION_INDEX=1
                      !Loop over the face normals of the element
                      DO ni=1,BASIS%NUMBER_OF_XI
                        !Determine the face xi directions that lie in this xi direction
                        FACE_XI(1)=OTHER_XI_DIRECTIONS3(ni,2,1)
                        FACE_XI(2)=OTHER_XI_DIRECTIONS3(ni,3,1)
                        !Reset the node_position_index in this xi direction
                        NODE_POSITION_INDEX(ni)=1
                        !Loop over the two faces with this normal
                        DO direction_index=-1,1,2
                          xi_direction=direction_index*ni
                          FACE_COLLAPSED(xi_direction)=.FALSE.
                          DO j=1,2
                            xi_dir_check=FACE_XI(j)
                            IF(xi_dir_check<=BASIS%NUMBER_OF_XI) THEN
                              xi_dir_search=FACE_XI(3-j)
                              NODE_POSITION_INDEX(xi_dir_search)=1
                              XI_COLLAPSED=.TRUE.
                              DO WHILE(NODE_POSITION_INDEX(xi_dir_search)<=NUMBER_OF_NODES_XIC(xi_dir_search).AND.XI_COLLAPSED)
                                !Get the first local node along the xi check direction
                                NODE_POSITION_INDEX(xi_dir_check)=1
                                nn1=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2), &
                                  & NODE_POSITION_INDEX(3),1)
                                !Get the second local node along the xi check direction
                                NODE_POSITION_INDEX(xi_dir_check)=2
                                nn2=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2), &
                                  & NODE_POSITION_INDEX(3),1)
                                IF(nn1/=0.AND.nn2/=0) THEN
                                  IF(DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn1)/= &
                                    & DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn2)) XI_COLLAPSED=.FALSE.
                                ENDIF
                                NODE_POSITION_INDEX(xi_dir_search)=NODE_POSITION_INDEX(xi_dir_search)+1
                              ENDDO !xi_dir_search
                              IF(XI_COLLAPSED) FACE_COLLAPSED(xi_direction)=.TRUE.
                            ENDIF
                          ENDDO !j
                          NODE_POSITION_INDEX(ni)=NUMBER_OF_NODES_XIC(ni)
                        ENDDO !direction_index
                      ENDDO !ni
                      !Loop over the xi directions and calculate the surrounding elements
                      DO ni=1,BASIS%NUMBER_OF_XI
                        !Determine the xi directions that lie in this xi direction
                        FACE_XI(1)=OTHER_XI_DIRECTIONS3(ni,2,1)
                        FACE_XI(2)=OTHER_XI_DIRECTIONS3(ni,3,1)
                        !Loop over the two faces
                        DO direction_index=-1,1,2
                          xi_direction=direction_index*ni                  
                          !Find nodes in the element on the appropriate face/line/point
                          NULLIFY(NODE_MATCH_LIST)
                          CALL LIST_CREATE_START(NODE_MATCH_LIST,ERR,ERROR,*999)
                          CALL LIST_DATA_TYPE_SET(NODE_MATCH_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                          CALL LIST_INITIAL_SIZE_SET(NODE_MATCH_LIST,16,ERR,ERROR,*999)
                          CALL LIST_CREATE_FINISH(NODE_MATCH_LIST,ERR,ERROR,*999)
                          IF(direction_index==-1) THEN
                            NODE_POSITION_INDEX(ni)=1
                          ELSE
                            NODE_POSITION_INDEX(ni)=NUMBER_OF_NODES_XIC(ni)
                          ENDIF
                          !If the face is collapsed then don't look in this xi direction. The exception is if the opposite face is
                          !also collapsed. This may indicate that we have a funny element in non-rc coordinates that goes around the
                          !central axis back to itself
                          IF(FACE_COLLAPSED(xi_direction).AND..NOT.FACE_COLLAPSED(-xi_direction)) THEN
                            !Do nothing - the match lists are already empty
                          ELSE
                            !Find the nodes to match and add them to the node match list
                            SELECT CASE(BASIS%NUMBER_OF_XI)
                            CASE(1)
                              nn=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),1,1,1)
                              IF(nn/=0) THEN
                                np=DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn)
                                CALL LIST_ITEM_ADD(NODE_MATCH_LIST,np,ERR,ERROR,*999)
                              ENDIF
                            CASE(2)
                              DO nn1=1,NUMBER_OF_NODES_XIC(FACE_XI(1)),NUMBER_OF_NODES_XIC(FACE_XI(1))-1
                                NODE_POSITION_INDEX(FACE_XI(1))=nn1
                                nn=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2),1,1)
                                IF(nn/=0) THEN
                                  np=DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn)
                                  CALL LIST_ITEM_ADD(NODE_MATCH_LIST,np,ERR,ERROR,*999)
                                ENDIF
                              ENDDO !nn1
                            CASE(3)
                              DO nn1=1,NUMBER_OF_NODES_XIC(FACE_XI(1)),NUMBER_OF_NODES_XIC(FACE_XI(1))-1
                                NODE_POSITION_INDEX(FACE_XI(1))=nn1
                                DO nn2=1,NUMBER_OF_NODES_XIC(FACE_XI(2)),NUMBER_OF_NODES_XIC(FACE_XI(2))-1
                                  NODE_POSITION_INDEX(FACE_XI(2))=nn2
                                  nn=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2), &
                                    & NODE_POSITION_INDEX(3),1)
                                  IF(nn/=0) THEN
                                    np=DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn)
                                    CALL LIST_ITEM_ADD(NODE_MATCH_LIST,np,ERR,ERROR,*999)
                                  ENDIF
                                ENDDO !nn2
                              ENDDO !nn1
                            CASE DEFAULT
                              LOCAL_ERROR="The number of xi directions in the basis of "// &
                                & TRIM(NUMBER_TO_VSTRING(BASIS%NUMBER_OF_XI,"*",ERR,ERROR))//" is invalid."
                              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                            END SELECT
                          ENDIF
                          CALL LIST_REMOVE_DUPLICATES(NODE_MATCH_LIST,ERR,ERROR,*999)
                          CALL LIST_DETACH_AND_DESTROY(NODE_MATCH_LIST,NUMBER_NODE_MATCHES,NODE_MATCHES,ERR,ERROR,*999)
                          NUMBER_SURROUNDING=0
                          IF(NUMBER_NODE_MATCHES>0) THEN
                            !NODE_MATCHES now contain the list of corner nodes in the current face with normal_xi of ni.
                            !Look at the surrounding elements of each of these nodes, if there is a repeated element that
                            !is not the current element ne, it's an adjacent element.
                            candidate_idx=0
                            NULLIFY(surroundingElementsList)
                            CALL LIST_CREATE_START(surroundingElementsList,ERR,ERROR,*999)
                            CALL LIST_DATA_TYPE_SET(surroundingElementsList,LIST_INTG_TYPE,ERR,ERROR,*999)
                            CALL LIST_INITIAL_SIZE_SET(surroundingElementsList,2,ERR,ERROR,*999)
                            CALL LIST_CREATE_FINISH(surroundingElementsList,ERR,ERROR,*999)
                            DO face_node_idx=1,NUMBER_NODE_MATCHES
                              !Dump all the surrounding elements into an array, see if any are repeated
                              node_idx=NODE_MATCHES(face_node_idx)
                              DO surrounding_el_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_SURROUNDING_ELEMENTS
                                candidate_el=DOMAIN_NODES%NODES(node_idx)%SURROUNDING_ELEMENTS(surrounding_el_idx)
                                IF(candidate_el/=ne) THEN
                                  candidate_idx=candidate_idx+1
                                  CALL LIST_ITEM_ADD(surroundingElementsList,candidate_el,ERR,ERROR,*999)
                                ENDIF
                              ENDDO
                            ENDDO !face_node_idx
                            CALL LIST_DETACH_AND_DESTROY(surroundingElementsList,numberSurroundingElements,surroundingElements, &
                              & ERR,ERROR,*999)
                            DO idx=1,candidate_idx
                              ne1=surroundingElements(idx)
                              IF(COUNT(surroundingElements(1:numberSurroundingElements)==ne1)>=BASIS%NUMBER_OF_XI) THEN
                                !Found it, just exit
                                CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(xi_direction)%PTR,ne1,ERR,ERROR,*999)
                                NUMBER_SURROUNDING=NUMBER_SURROUNDING+1
                                EXIT
                              ENDIF
                            ENDDO
                          ENDIF
                          IF(ALLOCATED(NODE_MATCHES)) DEALLOCATE(NODE_MATCHES)
                          IF(ALLOCATED(surroundingElements)) DEALLOCATE(surroundingElements)
                        ENDDO !direction_index
                      ENDDO !ni
                    CASE(BASIS_SIMPLEX_TYPE)
                      !Loop over the xi coordinates and calculate the surrounding elements
                      DO nic=1,BASIS%NUMBER_OF_XI_COORDINATES
                        !Find the other coordinates of the face/line/point
                        FACE_XIC(1)=OTHER_XI_DIRECTIONS4(nic,1)
                        FACE_XIC(2)=OTHER_XI_DIRECTIONS4(nic,2)
                        FACE_XIC(3)=OTHER_XI_DIRECTIONS4(nic,3)
                        !Find nodes in the element on the appropriate face/line/point
                        NULLIFY(NODE_MATCH_LIST)
                        CALL LIST_CREATE_START(NODE_MATCH_LIST,ERR,ERROR,*999)
                        CALL LIST_DATA_TYPE_SET(NODE_MATCH_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                        CALL LIST_INITIAL_SIZE_SET(NODE_MATCH_LIST,16,ERR,ERROR,*999)
                        CALL LIST_CREATE_FINISH(NODE_MATCH_LIST,ERR,ERROR,*999)
                        NODE_POSITION_INDEX(nic)=1 !Furtherest away from node with the nic'th coordinate
                        !Find the nodes to match and add them to the node match list
                        DO nn1=1,NUMBER_OF_NODES_XIC(FACE_XIC(1))
                          NODE_POSITION_INDEX(FACE_XIC(1))=nn1
                          DO nn2=1,NUMBER_OF_NODES_XIC(FACE_XIC(2))
                            NODE_POSITION_INDEX(FACE_XIC(2))=nn2
                            DO nn3=1,NUMBER_OF_NODES_XIC(FACE_XIC(3))
                              NODE_POSITION_INDEX(FACE_XIC(3))=nn3
                              nn=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2), &
                                & NODE_POSITION_INDEX(3),NODE_POSITION_INDEX(4))
                              IF(nn/=0) THEN
                                np=DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn)
                                CALL LIST_ITEM_ADD(NODE_MATCH_LIST,np,ERR,ERROR,*999)
                              ENDIF
                            ENDDO !nn3
                          ENDDO !nn2
                        ENDDO !nn1
                        CALL LIST_REMOVE_DUPLICATES(NODE_MATCH_LIST,ERR,ERROR,*999)
                        CALL LIST_DETACH_AND_DESTROY(NODE_MATCH_LIST,NUMBER_NODE_MATCHES,NODE_MATCHES,ERR,ERROR,*999)
                        IF(NUMBER_NODE_MATCHES>0) THEN
                          !Find list of elements surrounding those nodes
                          DO node_idx=1,NUMBER_NODE_MATCHES
                            np1=NODE_MATCHES(node_idx)
                            DO nep1=1,DOMAIN_NODES%NODES(np1)%NUMBER_OF_SURROUNDING_ELEMENTS
                              ne1=DOMAIN_NODES%NODES(np1)%SURROUNDING_ELEMENTS(nep1)
                              IF(ne1/=ne) THEN !Don't want the current element
                                ! grab the nodes list for current and this surrouding elements
                                ! current face : NODE_MATCHES
                                ! candidate elem : TOPOLOGY%ELEMENTS%ELEMENTS(ne1)%MESH_ELEMENT_NODES 
                                ! if all of current face belongs to the candidate element, we will have found the neighbour
                                CALL LIST_SUBSET_OF(NODE_MATCHES(1:NUMBER_NODE_MATCHES),DOMAIN_ELEMENTS%ELEMENTS(ne1)% &
                                  & ELEMENT_NODES,SUBSET,ERR,ERROR,*999)
                                IF(SUBSET) THEN
                                  CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(nic)%PTR,ne1,ERR,ERROR,*999)
                                ENDIF
                              ENDIF
                            ENDDO !nep1
                          ENDDO !node_idx
                        ENDIF
                        IF(ALLOCATED(NODE_MATCHES)) DEALLOCATE(NODE_MATCHES)
                      ENDDO !nic
                    CASE(BASIS_SERENDIPITY_TYPE)
                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                    CASE(BASIS_AUXILLIARY_TYPE)
                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                    CASE(BASIS_B_SPLINE_TP_TYPE)
                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                    CASE(BASIS_FOURIER_LAGRANGE_HERMITE_TP_TYPE)
                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                    CASE(BASIS_EXTENDED_LAGRANGE_TP_TYPE)
                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The basis type of "//TRIM(NUMBER_TO_VSTRING(BASIS%TYPE,"*",ERR,ERROR))// &
                        & " is invalid."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)                     
                    END SELECT
                    !Set the surrounding elements for this element
                    ALLOCATE(DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(-BASIS%NUMBER_OF_XI_COORDINATES: &
                      BASIS%NUMBER_OF_XI_COORDINATES),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate adjacent elements.",ERR,ERROR,*999)
                    DO nic=-BASIS%NUMBER_OF_XI_COORDINATES,BASIS%NUMBER_OF_XI_COORDINATES
                      CALL DECOMPOSITION_ADJACENT_ELEMENT_INITIALISE(DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic), &
                        & ERR,ERROR,*999)
                      CALL LIST_DETACH_AND_DESTROY(ADJACENT_ELEMENTS_LIST(nic)%PTR,DECOMPOSITION_ELEMENTS%ELEMENTS(ne)% &
                        & ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS,ADJACENT_ELEMENTS,ERR,ERROR,*999)
                      ALLOCATE(DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%ADJACENT_ELEMENTS( &
                        DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS),STAT=ERR)
                      IF(ERR/=0) CALL FlagError("Could not allocate element adjacent elements.",ERR,ERROR,*999)
                      DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%ADJACENT_ELEMENTS(1: &
                        & DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS)= &
                        ADJACENT_ELEMENTS(1:DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS)
                      IF(ALLOCATED(ADJACENT_ELEMENTS)) DEALLOCATE(ADJACENT_ELEMENTS)
                    ENDDO !nic
                  ENDDO !ne           
                ELSE
                  CALL FlagError("Domain topology elements is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Domain topology nodes is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Topology decomposition domain topology is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Topology decomposition domain is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Topology elements is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Topology decomposition is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Topology is not allocated.",ERR,ERROR,*999)
    ENDIF

    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Total number of elements = ",DECOMPOSITION_ELEMENTS% &
        & TOTAL_NUMBER_OF_ELEMENTS,ERR,ERROR,*999)
      DO ne=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
        BASIS=>DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Local element number : ",ne,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of xi coordinates = ",BASIS%NUMBER_OF_XI_COORDINATES, &
          & ERR,ERROR,*999)
        DO nic=-BASIS%NUMBER_OF_XI_COORDINATES,BASIS%NUMBER_OF_XI_COORDINATES
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Xi coordinate : ",nic,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of adjacent elements = ", &
            & DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS,ERR,ERROR,*999)
          IF(DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS>0) THEN
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DECOMPOSITION_ELEMENTS%ELEMENTS(ne)% &
              & ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS,8,8,DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)% &
              & ADJACENT_ELEMENTS,'("        Adjacent elements :",8(X,I6))','(30x,8(X,I6))',ERR,ERROR,*999)
          ENDIF
        ENDDO !nic
      ENDDO !ne
    ENDIF

    EXITS("DecompositionTopology_ElementAdjacentElementCalculate")
    RETURN
999 IF(ALLOCATED(NODE_MATCHES)) DEALLOCATE(NODE_MATCHES)
    IF(ALLOCATED(surroundingElements)) DEALLOCATE(surroundingElements)
    IF(ALLOCATED(ADJACENT_ELEMENTS)) DEALLOCATE(ADJACENT_ELEMENTS)
    IF(ASSOCIATED(NODE_MATCH_LIST)) CALL LIST_DESTROY(NODE_MATCH_LIST,DUMMY_ERR,DUMMY_ERROR,*998)
998 IF(ASSOCIATED(surroundingElementsList)) CALL LIST_DESTROY(surroundingElementsList,DUMMY_ERR,DUMMY_ERROR,*997)
997 DO nic=-4,4
      IF(ASSOCIATED(ADJACENT_ELEMENTS_LIST(nic)%PTR)) CALL LIST_DESTROY(ADJACENT_ELEMENTS_LIST(nic)%PTR,DUMMY_ERR,DUMMY_ERROR,*996)
    ENDDO !ni
996 ERRORS("DecompositionTopology_ElementAdjacentElementCalculate",ERR,ERROR)
    EXITS("DecompositionTopology_ElementAdjacentElementCalculate")
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_ElementAdjacentElementCalculate
  
  !
  !================================================================================================================================
  !

  !>Calculates the decomposition element topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to calculate the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: global_element,INSERT_STATUS,local_element
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: DECOMPOSITION_ELEMENTS
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_ELEMENTS_MAPPING
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: DOMAIN_TOPOLOGY
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(MeshElementsType), POINTER :: MESH_ELEMENTS
    TYPE(MeshComponentTopologyType), POINTER :: MESH_TOPOLOGY

    ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      DECOMPOSITION_ELEMENTS=>TOPOLOGY%ELEMENTS
      IF(ASSOCIATED(DECOMPOSITION_ELEMENTS)) THEN
        DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
        IF(ASSOCIATED(DECOMPOSITION)) THEN
          DOMAIN=>DECOMPOSITION%DOMAIN(DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
          IF(ASSOCIATED(DOMAIN)) THEN
            DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
            IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
              DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS
              IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                DOMAIN_MAPPINGS=>DOMAIN%MAPPINGS
                IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
                  DOMAIN_ELEMENTS_MAPPING=>DOMAIN_MAPPINGS%ELEMENTS
                  IF(ASSOCIATED(DOMAIN_ELEMENTS_MAPPING)) THEN
                    MESH=>DECOMPOSITION%MESH
                    IF(ASSOCIATED(MESH)) THEN
                      MESH_TOPOLOGY=>MESH%TOPOLOGY(DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
                      IF(ASSOCIATED(MESH_TOPOLOGY)) THEN
                        MESH_ELEMENTS=>MESH_TOPOLOGY%ELEMENTS
                        IF(ASSOCIATED(MESH_ELEMENTS)) THEN
                          !Allocate the element topology arrays
                          ALLOCATE(DECOMPOSITION_ELEMENTS%ELEMENTS(DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS),STAT=ERR)
                          IF(ERR/=0) CALL FlagError("Could not allocate decomposition elements elements.",ERR,ERROR,*999)
                          DECOMPOSITION_ELEMENTS%NUMBER_OF_ELEMENTS=DOMAIN_ELEMENTS%NUMBER_OF_ELEMENTS
                          DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS=DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                          DECOMPOSITION_ELEMENTS%NUMBER_OF_GLOBAL_ELEMENTS=DOMAIN_ELEMENTS%NUMBER_OF_GLOBAL_ELEMENTS
                          CALL TREE_CREATE_START(DECOMPOSITION_ELEMENTS%ELEMENTS_TREE,ERR,ERROR,*999)
                          CALL TREE_INSERT_TYPE_SET(DECOMPOSITION_ELEMENTS%ELEMENTS_TREE,TREE_NO_DUPLICATES_ALLOWED,ERR,ERROR,*999)
                          CALL TREE_CREATE_FINISH(DECOMPOSITION_ELEMENTS%ELEMENTS_TREE,ERR,ERROR,*999)
                          DO local_element=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                            CALL DECOMPOSITION_TOPOLOGY_ELEMENT_INITIALISE(DECOMPOSITION_ELEMENTS%ELEMENTS(local_element), &
                              & ERR,ERROR,*999)
                            global_element=DOMAIN_ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(local_element)
                            DECOMPOSITION_ELEMENTS%ELEMENTS(local_element)%USER_NUMBER=MESH_ELEMENTS%ELEMENTS(global_element)% &
                              & USER_NUMBER
                            DECOMPOSITION_ELEMENTS%ELEMENTS(local_element)%LOCAL_NUMBER=local_element
                            CALL TREE_ITEM_INSERT(DECOMPOSITION_ELEMENTS%ELEMENTS_TREE,DECOMPOSITION_ELEMENTS% &
                              & ELEMENTS(local_element)%USER_NUMBER,local_element,INSERT_STATUS,ERR,ERROR,*999)
                            DECOMPOSITION_ELEMENTS%ELEMENTS(local_element)%GLOBAL_NUMBER=global_element
                            DECOMPOSITION_ELEMENTS%ELEMENTS(local_element)%BOUNDARY_ELEMENT=MESH_ELEMENTS% &
                              & ELEMENTS(global_element)%BOUNDARY_ELEMENT
                          ENDDO !local_element
                          !Calculate the elements surrounding the elements in the decomposition topology
                          CALL DecompositionTopology_ElementAdjacentElementCalculate(TOPOLOGY,ERR,ERROR,*999)
                        ELSE
                          CALL FlagError("Mesh elements is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Mesh topology is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Domain mappings elements is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Domain mappings is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Domain topology elements is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Topology decomposition domain topology is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Topology decomposition domain is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Topology decomposition is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Topology elements is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Finalises the elements in the given decomposition topology. \todo Pass in the decomposition elements pointer.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to finalise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ne

    ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        DO ne=1,TOPOLOGY%ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
          CALL DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE(TOPOLOGY%ELEMENTS%ELEMENTS(ne),ERR,ERROR,*999)
        ENDDO !ne
        IF(ASSOCIATED(TOPOLOGY%ELEMENTS%ELEMENTS)) DEALLOCATE(TOPOLOGY%ELEMENTS%ELEMENTS)
        IF(ASSOCIATED(TOPOLOGY%ELEMENTS%ELEMENTS_TREE)) CALL TREE_DESTROY(TOPOLOGY%ELEMENTS%ELEMENTS_TREE,ERR,ERROR,*999)
        DEALLOCATE(TOPOLOGY%ELEMENTS)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the element data structures for a decomposition topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to initialise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        CALL FlagError("Decomposition already has topology elements associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%ELEMENTS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology elements.",ERR,ERROR,*999)
        TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS=0
        TOPOLOGY%ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS=0
        TOPOLOGY%ELEMENTS%NUMBER_OF_GLOBAL_ELEMENTS=0
        TOPOLOGY%ELEMENTS%DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
        NULLIFY(TOPOLOGY%ELEMENTS%ELEMENTS)
        NULLIFY(TOPOLOGY%ELEMENTS%ELEMENTS_TREE)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the topology in the given decomposition. \todo Pass in a pointer to the decomposition topology
  SUBROUTINE DECOMPOSITION_TOPOLOGY_FINALISE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to finalise the topology for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      CALL DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      IF(DECOMPOSITION%CALCULATE_LINES) THEN
        CALL DECOMPOSITION_TOPOLOGY_LINES_FINALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      ENDIF
      IF(DECOMPOSITION%CALCULATE_FACES) THEN
        CALL DECOMPOSITION_TOPOLOGY_FACES_FINALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      ENDIF
      DEALLOCATE(DECOMPOSITION%TOPOLOGY)
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DECOMPOSITION_TOPOLOGY_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the topology for a given decomposition.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_INITIALISE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to initialise the topology for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: meshComponentNumber

    ENTERS("DECOMPOSITION_TOPOLOGY_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(ASSOCIATED(DECOMPOSITION%TOPOLOGY)) THEN
        CALL FlagError("Decomposition already has topology associated.",ERR,ERROR,*999)
      ELSE
        !Allocate decomposition topology
        ALLOCATE(DECOMPOSITION%TOPOLOGY,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Decomposition topology could not be allocated.",ERR,ERROR,*999)
        DECOMPOSITION%TOPOLOGY%DECOMPOSITION=>DECOMPOSITION
        NULLIFY(DECOMPOSITION%TOPOLOGY%ELEMENTS)
        NULLIFY(DECOMPOSITION%TOPOLOGY%LINES)
        NULLIFY(DECOMPOSITION%TOPOLOGY%FACES)
        NULLIFY(DECOMPOSITION%TOPOLOGY%dataPoints)
        !Initialise the topology components
        CALL DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
        IF(DECOMPOSITION%CALCULATE_LINES) THEN !Default is currently true
          CALL DECOMPOSITION_TOPOLOGY_LINES_INITIALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
        ENDIF
        IF(DECOMPOSITION%CALCULATE_FACES) THEN !Default is currently false
          CALL DECOMPOSITION_TOPOLOGY_FACES_INITIALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
        ENDIF      
        meshComponentNumber=DECOMPOSITION%MESH_COMPONENT_NUMBER
        IF(ALLOCATED(DECOMPOSITION%MESH%TOPOLOGY(meshComponentNumber)%PTR%dataPoints%dataPoints)) THEN
          CALL DecompositionTopology_DataPointsInitialise(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
        ENDIF   
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TOPOLOGY_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises a line in the given decomposition topology and deallocates all memory.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_LINE_FINALISE(LINE,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_LINE_TYPE) :: LINE !<The decomposition line to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_LINE_FINALISE",ERR,ERROR,*999)

    LINE%NUMBER=0
    LINE%XI_DIRECTION=0
    LINE%NUMBER_OF_SURROUNDING_ELEMENTS=0
    IF(ALLOCATED(LINE%SURROUNDING_ELEMENTS)) DEALLOCATE(LINE%SURROUNDING_ELEMENTS)
    IF(ALLOCATED(LINE%ELEMENT_LINES)) DEALLOCATE(LINE%ELEMENT_LINES)
    LINE%ADJACENT_LINES=0
 
    EXITS("DECOMPOSITION_TOPOLOGY_LINE_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_LINE_FINALISE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_LINE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the line data structure for a decomposition topology line.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_LINE_INITIALISE(LINE,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_LINE_TYPE) :: LINE !<The decomposition line to initialise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_LINE_INITIALISE",ERR,ERROR,*999)

    LINE%NUMBER=0
    LINE%XI_DIRECTION=0
    LINE%NUMBER_OF_SURROUNDING_ELEMENTS=0
    LINE%ADJACENT_LINES=0
    LINE%BOUNDARY_LINE=.FALSE.

    EXITS("DECOMPOSITION_TOPOLOGY_LINE_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_LINE_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_LINE_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Calculates the lines in the given decomposition topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to calculate the lines for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,element_idx,surrounding_element_idx,basis_local_line_idx, &
      & surrounding_element_basis_local_line_idx,element_local_node_idx,basis_local_line_node_idx,derivative_idx,version_idx, &
      & local_line_idx,surrounding_element_local_line_idx,node_idx,local_node_idx,elem_idx,line_end_node_idx,basis_node_idx, &
      & NODES_IN_LINE(4),NUMBER_OF_LINES,MAX_NUMBER_OF_LINES,NEW_MAX_NUMBER_OF_LINES,LINE_NUMBER,COUNT
    INTEGER(INTG), ALLOCATABLE :: NODES_NUMBER_OF_LINES(:)
    INTEGER(INTG), POINTER :: TEMP_LINES(:,:),NEW_TEMP_LINES(:,:)
    REAL(DP) :: APPROX_DIMENSION
    LOGICAL :: FOUND
    TYPE(BASIS_TYPE), POINTER :: BASIS,BASIS2
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DECOMPOSITION_ELEMENT_TYPE), POINTER :: DECOMPOSITION_ELEMENT
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: DECOMPOSITION_ELEMENTS
    TYPE(DECOMPOSITION_LINE_TYPE), POINTER :: DECOMPOSITION_LINE,DECOMPOSITION_LINE2
    TYPE(DECOMPOSITION_LINES_TYPE), POINTER :: DECOMPOSITION_LINES
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_ELEMENT_TYPE), POINTER :: DOMAIN_ELEMENT
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS
    TYPE(DOMAIN_LINE_TYPE), POINTER :: DOMAIN_LINE,DOMAIN_LINE2
    TYPE(DOMAIN_LINES_TYPE), POINTER :: DOMAIN_LINES
    TYPE(DOMAIN_NODE_TYPE), POINTER :: DOMAIN_NODE
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: DOMAIN_TOPOLOGY
    TYPE(MESH_TYPE), POINTER :: MESH

    NULLIFY(TEMP_LINES)
    NULLIFY(NEW_TEMP_LINES)
    
    ENTERS("DECOMPOSITION_TOPOLOGY_LINES_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      DECOMPOSITION_LINES=>TOPOLOGY%LINES
      IF(ASSOCIATED(DECOMPOSITION_LINES)) THEN
        DECOMPOSITION_ELEMENTS=>TOPOLOGY%ELEMENTS
        IF(ASSOCIATED(DECOMPOSITION_ELEMENTS)) THEN
          DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
          IF(ASSOCIATED(DECOMPOSITION)) THEN
            !Process the mesh component number (component number the decomposition was calculated from) first to establish line
            !topology then process the other mesh components.
            DOMAIN=>DECOMPOSITION%DOMAIN(DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
            IF(ASSOCIATED(DOMAIN)) THEN
              DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
              IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
                DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
                IF(ASSOCIATED(DOMAIN_NODES)) THEN
                  DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS
                  IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                    !Guestimate the number of lines
                    SELECT CASE(DOMAIN%NUMBER_OF_DIMENSIONS)
                    CASE(1)
                      MAX_NUMBER_OF_LINES=DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                    CASE(2)
                      APPROX_DIMENSION=SQRT(REAL(DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS,DP))
                      !This should give the maximum and will over estimate the number of lines for a "square mesh" by approx 33%
                      MAX_NUMBER_OF_LINES=NINT(3.0_DP*APPROX_DIMENSION*(APPROX_DIMENSION+1),INTG)
                    CASE(3)
                      !This should give the maximum and will over estimate the number of lines for a "cube mesh" by approx 73%
                      APPROX_DIMENSION=REAL(DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS,DP)**(1.0_DP/3.0_DP)
                      MAX_NUMBER_OF_LINES=NINT(11.0_DP*APPROX_DIMENSION*APPROX_DIMENSION*(APPROX_DIMENSION+1),INTG)
                    CASE DEFAULT
                      CALL FlagError("Invalid number of dimensions for a topology domain.",ERR,ERROR,*999)
                    END SELECT
                    DOMAIN_LINES=>DOMAIN_TOPOLOGY%LINES
                    IF(ASSOCIATED(DOMAIN_LINES)) THEN
                      ALLOCATE(TEMP_LINES(4,MAX_NUMBER_OF_LINES),STAT=ERR)
                      IF(ERR/=0) CALL FlagError("Could not allocate temporary lines array.",ERR,ERROR,*999)
                      ALLOCATE(NODES_NUMBER_OF_LINES(DOMAIN_NODES%TOTAL_NUMBER_OF_NODES),STAT=ERR)
                      IF(ERR/=0) CALL FlagError("Could not allocate nodes number of lines array.",ERR,ERROR,*999)
                      NODES_NUMBER_OF_LINES=0
                      NUMBER_OF_LINES=0
                      TEMP_LINES=0
                      !Loop over the elements in the topology
                      DO element_idx=1,DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                        DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(element_idx)
                        DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(element_idx)
                        BASIS=>DOMAIN_ELEMENT%BASIS
                        ALLOCATE(DECOMPOSITION_ELEMENT%ELEMENT_LINES(BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate element element lines.",ERR,ERROR,*999)
                        !Loop over the local lines of the element
                        DO basis_local_line_idx=1,BASIS%NUMBER_OF_LOCAL_LINES
                          !Calculate the topology node numbers that make up the line
                          NODES_IN_LINE=0
                          DO basis_local_line_node_idx=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx)
                            NODES_IN_LINE(basis_local_line_node_idx)=DOMAIN_ELEMENT%ELEMENT_NODES( &
                              & BASIS%NODE_NUMBERS_IN_LOCAL_LINE(basis_local_line_node_idx,basis_local_line_idx))
                          ENDDO !basis_local_line_node_idx
                          !Try and find a previously created line that matches in the adjacent elements
                          FOUND=.FALSE.
                          node_idx=NODES_IN_LINE(1)
                          DO elem_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_SURROUNDING_ELEMENTS
                            surrounding_element_idx=DOMAIN_NODES%NODES(node_idx)%SURROUNDING_ELEMENTS(elem_idx)
                            IF(surrounding_element_idx/=element_idx) THEN
                              IF(ALLOCATED(DECOMPOSITION_ELEMENTS%ELEMENTS(surrounding_element_idx)%ELEMENT_LINES)) THEN
                                BASIS2=>DOMAIN_ELEMENTS%ELEMENTS(surrounding_element_idx)%BASIS
                                DO surrounding_element_basis_local_line_idx=1,BASIS2%NUMBER_OF_LOCAL_LINES
                                  local_line_idx=DECOMPOSITION_ELEMENTS%ELEMENTS(surrounding_element_idx)% &
                                    & ELEMENT_LINES(surrounding_element_basis_local_line_idx)
                                  IF(ALL(NODES_IN_LINE(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx))== &
                                    & TEMP_LINES(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx),local_line_idx))) THEN
                                    FOUND=.TRUE.
                                    EXIT
                                  ENDIF
                                ENDDO !surrounding_element_basis_local_line_idx
                                IF(FOUND) EXIT
                              ENDIF
                            ENDIF
                          ENDDO !elem_idx
                          IF(FOUND) THEN
                            !Line has already been created
                            DECOMPOSITION_ELEMENT%ELEMENT_LINES(basis_local_line_idx)=local_line_idx
                          ELSE
                            !Line has not been created
                            IF(NUMBER_OF_LINES==MAX_NUMBER_OF_LINES) THEN
                              !We are at maximum. Reallocate the LINES array to be 20% bigger and try again.
                              NEW_MAX_NUMBER_OF_LINES=NINT(1.20_DP*REAL(MAX_NUMBER_OF_LINES,DP),INTG)
                              ALLOCATE(NEW_TEMP_LINES(4,NEW_MAX_NUMBER_OF_LINES),STAT=ERR)
                              IF(ERR/=0) CALL FlagError("Could not allocate new number of lines.",ERR,ERROR,*999)
                              NEW_TEMP_LINES(:,1:NUMBER_OF_LINES)=TEMP_LINES(:,1:NUMBER_OF_LINES)
                              NEW_TEMP_LINES(:,NUMBER_OF_LINES+1:NEW_MAX_NUMBER_OF_LINES)=0
                              DEALLOCATE(TEMP_LINES)
                              TEMP_LINES=>NEW_TEMP_LINES
                              NULLIFY(NEW_TEMP_LINES)
                              MAX_NUMBER_OF_LINES=NEW_MAX_NUMBER_OF_LINES
                            ENDIF
                            NUMBER_OF_LINES=NUMBER_OF_LINES+1
                            TEMP_LINES(:,NUMBER_OF_LINES)=NODES_IN_LINE
                            DECOMPOSITION_ELEMENT%ELEMENT_LINES(basis_local_line_idx)=NUMBER_OF_LINES
                            DO basis_local_line_node_idx=1,SIZE(NODES_IN_LINE,1)
                              IF(NODES_IN_LINE(basis_local_line_node_idx)/=0) &
                                & NODES_NUMBER_OF_LINES(NODES_IN_LINE(basis_local_line_node_idx))= &
                                & NODES_NUMBER_OF_LINES(NODES_IN_LINE(basis_local_line_node_idx))+1
                            ENDDO !basis_local_line_node_idx
                          ENDIF
                        ENDDO !basis_local_line_idx
                      ENDDO !element_idx
                      !Allocate the line arrays and set them from the LINES and NODE_LINES arrays
                      DO node_idx=1,DOMAIN_NODES%TOTAL_NUMBER_OF_NODES
                        ALLOCATE(DOMAIN_NODES%NODES(node_idx)%NODE_LINES(NODES_NUMBER_OF_LINES(node_idx)),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate node lines array.",ERR,ERROR,*999)
                        DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_LINES=0
                      ENDDO !node_idx
                      DEALLOCATE(NODES_NUMBER_OF_LINES)
                      ALLOCATE(DECOMPOSITION_LINES%LINES(NUMBER_OF_LINES),STAT=ERR)
                      IF(ERR/=0) CALL FlagError("Could not allocate decomposition topology lines.",ERR,ERROR,*999)
                      DECOMPOSITION_LINES%NUMBER_OF_LINES=NUMBER_OF_LINES
                      ALLOCATE(DOMAIN_LINES%LINES(NUMBER_OF_LINES),STAT=ERR)
                      IF(ERR/=0) CALL FlagError("Could not allocate domain topology lines.",ERR,ERROR,*999)
                      DOMAIN_LINES%NUMBER_OF_LINES=NUMBER_OF_LINES
                      DO local_line_idx=1,DOMAIN_LINES%NUMBER_OF_LINES
                        CALL DECOMPOSITION_TOPOLOGY_LINE_INITIALISE(DECOMPOSITION_LINES%LINES(local_line_idx),ERR,ERROR,*999)
                        CALL DOMAIN_TOPOLOGY_LINE_INITIALISE(DOMAIN_LINES%LINES(local_line_idx),ERR,ERROR,*999)
                        DO basis_local_line_node_idx=1,SIZE(TEMP_LINES,1)
                          IF(TEMP_LINES(basis_local_line_node_idx,local_line_idx)/=0) THEN
                            node_idx=TEMP_LINES(basis_local_line_node_idx,local_line_idx)
                            DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_LINES=DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_LINES+1
                            DOMAIN_NODES%NODES(node_idx)%NODE_LINES(DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_LINES)= &
                              & local_line_idx
                          ENDIF
                        ENDDO !basis_local_line_node_idx  
                      ENDDO !local_line_idx
                      DO element_idx=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                        DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(element_idx)
                        DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(element_idx)
                        BASIS=>DOMAIN_ELEMENT%BASIS
                        DO basis_local_line_idx=1,BASIS%NUMBER_OF_LOCAL_LINES
                          LINE_NUMBER=DECOMPOSITION_ELEMENT%ELEMENT_LINES(basis_local_line_idx)
                          DECOMPOSITION_LINE=>DECOMPOSITION_LINES%LINES(LINE_NUMBER)
                          DOMAIN_LINE=>DOMAIN_LINES%LINES(LINE_NUMBER)
                          DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS=DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS+1
                          IF(.NOT.ASSOCIATED(DOMAIN_LINE%BASIS)) THEN
                            DECOMPOSITION_LINE%NUMBER=LINE_NUMBER
                            DOMAIN_LINE%NUMBER=LINE_NUMBER
                            DOMAIN_LINE%ELEMENT_NUMBER=element_idx !Needs checking
                            DECOMPOSITION_LINE%XI_DIRECTION=BASIS%localLineXiDirection(basis_local_line_idx)
                            IF(ALLOCATED(BASIS%localLineBasis)) THEN
                              DOMAIN_LINE%BASIS=>BASIS%localLineBasis(basis_local_line_idx)%PTR
                            ELSE
                              !Basis is only 1D
                              DOMAIN_LINE%BASIS=>BASIS
                            ENDIF
                            ALLOCATE(DOMAIN_LINE%NODES_IN_LINE(BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx)),STAT=ERR)
                            IF(ERR/=0) CALL FlagError("Could not allocate line nodes in line.",ERR,ERROR,*999)
                            ALLOCATE(DOMAIN_LINE%DERIVATIVES_IN_LINE(2,DOMAIN_LINE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                              & BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx)),STAT=ERR)
                            IF(ERR/=0) CALL FlagError("Could not allocate line derivatives in line.",ERR,ERROR,*999)
                            DOMAIN_LINE%NODES_IN_LINE(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx))= &
                              & TEMP_LINES(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx),LINE_NUMBER)
                            DO basis_local_line_node_idx=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx)
                              !Set derivative number of u (NO_GLOBAL_DERIV) for the domain line
                              DOMAIN_LINE%DERIVATIVES_IN_LINE(1,1,basis_local_line_node_idx)=NO_GLOBAL_DERIV
                              !Set version number of u (NO_GLOBAL_DERIV) for the domain line
                              version_idx=DOMAIN_ELEMENT%elementVersions(1,BASIS%NODE_NUMBERS_IN_LOCAL_LINE( &
                                & basis_local_line_node_idx,basis_local_line_idx))
                              DOMAIN_LINE%DERIVATIVES_IN_LINE(2,1,basis_local_line_node_idx)=version_idx
                              IF(DOMAIN_LINE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES>1) THEN
                                derivative_idx=DOMAIN_ELEMENT%ELEMENT_DERIVATIVES( &
                                  & BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(basis_local_line_node_idx,basis_local_line_idx), &
                                  & BASIS%NODE_NUMBERS_IN_LOCAL_LINE(basis_local_line_node_idx,basis_local_line_idx))
                                DOMAIN_LINE%DERIVATIVES_IN_LINE(1,2,basis_local_line_node_idx)=derivative_idx
                                version_idx=DOMAIN_ELEMENT%elementVersions(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE( &
                                  & basis_local_line_node_idx,basis_local_line_idx),BASIS%NODE_NUMBERS_IN_LOCAL_LINE( &
                                  & basis_local_line_node_idx,basis_local_line_idx))
                                DOMAIN_LINE%DERIVATIVES_IN_LINE(2,2,basis_local_line_node_idx)=version_idx
                              ENDIF
                            ENDDO !basis_local_line_node_idx
                          ENDIF
                        ENDDO !basis_local_line_idx
                      ENDDO !element_idx
                      DEALLOCATE(TEMP_LINES)
                      !Calculate adjacent lines and the surrounding elements for each line
                      DO local_line_idx=1,DECOMPOSITION_LINES%NUMBER_OF_LINES
                        DECOMPOSITION_LINE=>DECOMPOSITION_LINES%LINES(local_line_idx)
                        DOMAIN_LINE=>DOMAIN_LINES%LINES(local_line_idx)
                        BASIS=>DOMAIN_LINE%BASIS
                        IF(DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS==1) THEN
                          DECOMPOSITION_LINE%BOUNDARY_LINE=.TRUE.
                          DOMAIN_LINE%BOUNDARY_LINE=.TRUE.
                        ENDIF
                        !Allocate the elements surrounding the line
                        ALLOCATE(DECOMPOSITION_LINE%SURROUNDING_ELEMENTS(DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS), &
                          & STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate line surrounding elements.",ERR,ERROR,*999)
                        ALLOCATE(DECOMPOSITION_LINE%ELEMENT_LINES(DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS), &
                          & STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate line element lines.",ERR,ERROR,*999)
                        DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS=0
                        DECOMPOSITION_LINE%ADJACENT_LINES=0
                        !Loop over the nodes at each end of the line
                        DO line_end_node_idx=0,1
                          FOUND=.FALSE.
                          node_idx=DOMAIN_LINE%NODES_IN_LINE(line_end_node_idx*(BASIS%NUMBER_OF_NODES-1)+1)
                          !Loop over the elements surrounding the node.
                          DO elem_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_SURROUNDING_ELEMENTS
                            element_idx=DOMAIN_NODES%NODES(node_idx)%SURROUNDING_ELEMENTS(elem_idx)
                            DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(element_idx)
                            DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(element_idx)
                            !Loop over the local lines of the element
                            DO basis_local_line_idx=1,DOMAIN_ELEMENT%BASIS%NUMBER_OF_LOCAL_LINES
                              surrounding_element_local_line_idx=DECOMPOSITION_ELEMENT%ELEMENT_LINES(basis_local_line_idx)
                              IF(surrounding_element_local_line_idx/=local_line_idx) THEN
                                DECOMPOSITION_LINE2=>DECOMPOSITION_LINES%LINES(surrounding_element_local_line_idx)
                                DOMAIN_LINE2=>DOMAIN_LINES%LINES(surrounding_element_local_line_idx)
                                IF(DECOMPOSITION_LINE2%XI_DIRECTION==DECOMPOSITION_LINE%XI_DIRECTION) THEN
                                  !Lines run in the same direction.
                                  BASIS2=>DOMAIN_LINE2%BASIS
                                  IF(line_end_node_idx==0) THEN
                                    local_node_idx=DOMAIN_LINE2%NODES_IN_LINE(BASIS2%NUMBER_OF_NODES)
                                  ELSE
                                    local_node_idx=DOMAIN_LINE2%NODES_IN_LINE(1)
                                  ENDIF
                                  IF(local_node_idx==node_idx) THEN
                                    !The node at the 'other' end of this line matches the node at the current end of the line.
                                    !Check it is not a coexistant line running the other way
                                    IF(BASIS2%INTERPOLATION_ORDER(1)==BASIS%INTERPOLATION_ORDER(1)) THEN
                                      COUNT=0
                                      DO basis_node_idx=1,BASIS%NUMBER_OF_NODES
                                        IF(DOMAIN_LINE2%NODES_IN_LINE(basis_node_idx)== &
                                          & DOMAIN_LINE%NODES_IN_LINE(BASIS2%NUMBER_OF_NODES-basis_node_idx+1)) &
                                          & COUNT=COUNT+1
                                      ENDDO !basis_node_idx
                                      IF(COUNT<BASIS%NUMBER_OF_NODES) THEN
                                        FOUND=.TRUE.
                                        EXIT
                                      ENDIF
                                    ELSE
                                      FOUND=.TRUE.
                                      EXIT
                                    ENDIF
                                  ENDIF
                                ENDIF
                              ENDIF
                            ENDDO !basis_local_line_idx
                            IF(FOUND) EXIT
                          ENDDO !element_idx
                          IF(FOUND) DECOMPOSITION_LINE%ADJACENT_LINES(line_end_node_idx)=surrounding_element_local_line_idx
                        ENDDO !line_end_node_idx
                      ENDDO !local_line_idx
                      !Set the surrounding elements
                      DO element_idx=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                        DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(element_idx)
                        DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(element_idx)
                        BASIS=>DOMAIN_ELEMENT%BASIS
                        DO basis_local_line_idx=1,BASIS%NUMBER_OF_LOCAL_LINES
                          LINE_NUMBER=DECOMPOSITION_ELEMENT%ELEMENT_LINES(basis_local_line_idx)
                          DECOMPOSITION_LINE=>DECOMPOSITION_LINES%LINES(LINE_NUMBER)
                          DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS=DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS+1
                          DECOMPOSITION_LINE%SURROUNDING_ELEMENTS(DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS)=element_idx
                          DECOMPOSITION_LINE%ELEMENT_LINES(DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS)=basis_local_line_idx
                        ENDDO !basis_local_line_idx
                      ENDDO !element_idx
                    ELSE
                      CALL FlagError("Domain topology lines is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Domain topology elements is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Domain topology nodes is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Topology decomposition domain topology is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Topology decomposition domain is not associated.",ERR,ERROR,*999)
            ENDIF
            !Now loop over the other mesh components in the decomposition and calculate the domain lines
            MESH=>DECOMPOSITION%MESH
            IF(ASSOCIATED(MESH)) THEN
              DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
                IF(component_idx/=DECOMPOSITION%MESH_COMPONENT_NUMBER) THEN
                  DOMAIN=>DECOMPOSITION%DOMAIN(component_idx)%PTR
                  IF(ASSOCIATED(DOMAIN)) THEN
                    DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
                    IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
                      DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
                      IF(ASSOCIATED(DOMAIN_NODES)) THEN
                        DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS
                        IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                          DOMAIN_LINES=>DOMAIN_TOPOLOGY%LINES                      
                          IF(ASSOCIATED(DOMAIN_LINES)) THEN
                            ALLOCATE(DOMAIN_LINES%LINES(DECOMPOSITION_LINES%NUMBER_OF_LINES),STAT=ERR)
                            IF(ERR/=0) CALL FlagError("Could not allocate domain lines lines.",ERR,ERROR,*999)
                            DOMAIN_LINES%NUMBER_OF_LINES=DECOMPOSITION_LINES%NUMBER_OF_LINES
                            ALLOCATE(NODES_NUMBER_OF_LINES(DOMAIN_NODES%TOTAL_NUMBER_OF_NODES),STAT=ERR)
                            IF(ERR/=0) CALL FlagError("Could not allocate nodes number of lines array.",ERR,ERROR,*999)
                            NODES_NUMBER_OF_LINES=0
                            !Loop over the lines in the topology
                            DO local_line_idx=1,DECOMPOSITION_LINES%NUMBER_OF_LINES
                              DECOMPOSITION_LINE=>DECOMPOSITION_LINES%LINES(local_line_idx)
                              DOMAIN_LINE=>DOMAIN_LINES%LINES(local_line_idx)
                              IF(DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS>0) THEN
                                element_idx=DECOMPOSITION_LINE%SURROUNDING_ELEMENTS(1)
                                basis_local_line_idx=DECOMPOSITION_LINE%ELEMENT_LINES(1)
                                CALL DOMAIN_TOPOLOGY_LINE_INITIALISE(DOMAIN_LINES%LINES(local_line_idx),ERR,ERROR,*999)
                                DOMAIN_LINE%NUMBER=local_line_idx
                                DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(element_idx)
                                BASIS=>DOMAIN_ELEMENT%BASIS
                                DOMAIN_LINE%ELEMENT_NUMBER=DOMAIN_ELEMENT%NUMBER
                                IF(ALLOCATED(BASIS%localLineBasis)) THEN
                                  DOMAIN_LINE%BASIS=>BASIS%localLineBasis(basis_local_line_idx)%PTR
                                ELSE
                                  !Basis is only 1D
                                  DOMAIN_LINE%BASIS=>BASIS
                                ENDIF
                                ALLOCATE(DOMAIN_LINE%NODES_IN_LINE(BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx)), &
                                  & STAT=ERR)
                                IF(ERR/=0) CALL FlagError("Could not allocate nodes in line.",ERR,ERROR,*999)
                                ALLOCATE(DOMAIN_LINE%DERIVATIVES_IN_LINE(2,DOMAIN_LINE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                                  & BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx)),STAT=ERR)
                                IF(ERR/=0) CALL FlagError("Could not allocate derivatives in line.",ERR,ERROR,*999)
                                DO basis_local_line_node_idx=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx)
                                  element_local_node_idx=BASIS%NODE_NUMBERS_IN_LOCAL_LINE(basis_local_line_node_idx, &
                                    & basis_local_line_idx)
                                  node_idx=DOMAIN_ELEMENT%ELEMENT_NODES(element_local_node_idx)
                                  DOMAIN_LINE%NODES_IN_LINE(basis_local_line_node_idx)=node_idx
                                  !Set derivative number of u (NO_GLOBAL_DERIV) for the domain line
                                  DOMAIN_LINE%DERIVATIVES_IN_LINE(1,1,basis_local_line_node_idx)=NO_GLOBAL_DERIV
                                  !Set version number of u (NO_GLOBAL_DERIV) for the domain line
                                  version_idx=DOMAIN_ELEMENT%elementVersions(1,BASIS%NODE_NUMBERS_IN_LOCAL_LINE( &
                                    & basis_local_line_node_idx,basis_local_line_idx))
                                  DOMAIN_LINE%DERIVATIVES_IN_LINE(2,1,basis_local_line_node_idx)=version_idx
                                  IF(DOMAIN_LINE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES>1) THEN
                                    derivative_idx=DOMAIN_ELEMENT%ELEMENT_DERIVATIVES(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE( &
                                      & basis_local_line_node_idx,basis_local_line_idx),element_local_node_idx)
                                    DOMAIN_LINE%DERIVATIVES_IN_LINE(1,2,basis_local_line_node_idx)=derivative_idx
                                    version_idx=DOMAIN_ELEMENT%elementVersions(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE( &
                                      & basis_local_line_node_idx,basis_local_line_idx),BASIS%NODE_NUMBERS_IN_LOCAL_LINE( &
                                      & basis_local_line_node_idx,basis_local_line_idx))
                                    DOMAIN_LINE%DERIVATIVES_IN_LINE(2,2,basis_local_line_node_idx)=version_idx
                                  ENDIF
                                  NODES_NUMBER_OF_LINES(node_idx)=NODES_NUMBER_OF_LINES(node_idx)+1
                                ENDDO !basis_local_line_node_idx
                              ELSE
                                CALL FlagError("Line is not surrounded by any elements?",ERR,ERROR,*999)
                              ENDIF
                            ENDDO !local_line_idx
                            DO node_idx=1,DOMAIN_NODES%TOTAL_NUMBER_OF_NODES
                              ALLOCATE(DOMAIN_NODES%NODES(node_idx)%NODE_LINES(NODES_NUMBER_OF_LINES(node_idx)),STAT=ERR)
                              IF(ERR/=0) CALL FlagError("Could not allocate node lines.",ERR,ERROR,*999)
                              DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_LINES=0
                            ENDDO !node_idx
                            DEALLOCATE(NODES_NUMBER_OF_LINES)
                            DO local_line_idx=1,DOMAIN_LINES%NUMBER_OF_LINES
                              DOMAIN_LINE=>DOMAIN_LINES%LINES(local_line_idx)
                              BASIS=>DOMAIN_LINE%BASIS
                              DO basis_local_line_node_idx=1,BASIS%NUMBER_OF_NODES
                                node_idx=DOMAIN_LINE%NODES_IN_LINE(basis_local_line_node_idx)
                                DOMAIN_NODE=>DOMAIN_NODES%NODES(node_idx)
                                DOMAIN_NODE%NUMBER_OF_NODE_LINES=DOMAIN_NODE%NUMBER_OF_NODE_LINES+1
                                DOMAIN_NODE%NODE_LINES(DOMAIN_NODE%NUMBER_OF_NODE_LINES)=local_line_idx
                              ENDDO !basis_local_line_node_idx
                            ENDDO !local_line_idx
                          ELSE
                            CALL FlagError("Domain lines is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FlagError("Domain elements is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Domain nodes is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Domain topology is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Decomposition mesh is not associated",ERR,ERROR,*999)
                  ENDIF
                ENDIF
              ENDDO !component_idx
            ELSE
              CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
            ENDIF                        
          ELSE
            CALL FlagError("Topology decomposition is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Topology decomposition elements is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Topology lines is not associated.",ERR,ERROR,*999)

      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Decomposition topology lines:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of mesh components = ",MESH%NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of lines = ",DECOMPOSITION_LINES%NUMBER_OF_LINES,ERR,ERROR,*999)
      DO local_line_idx=1,DECOMPOSITION_LINES%NUMBER_OF_LINES
        DECOMPOSITION_LINE=>DECOMPOSITION_LINES%LINES(local_line_idx)
        DOMAIN_LINE=>DOMAIN_LINES%LINES(local_line_idx)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Line number = ",DECOMPOSITION_LINE%NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Xi direction = ",DECOMPOSITION_LINE%XI_DIRECTION,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of surrounding elements = ", &
          & DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS,4,4, &
          & DECOMPOSITION_LINE%SURROUNDING_ELEMENTS,'("      Surrounding elements :",4(X,I8))','(28X,4(X,I8))',ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS,4,4, &
          & DECOMPOSITION_LINE%ELEMENT_LINES,'("      Element lines        :",4(X,I8))','(28X,4(X,I8))',ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,2,2,2,DECOMPOSITION_LINE%ADJACENT_LINES, &
          & '("      Adjacent lines       :",2(X,I8))','(28X,2(X,I8))',ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Boundary line = ",DECOMPOSITION_LINE%BOUNDARY_LINE,ERR,ERROR,*999)
        DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Mesh component : ",component_idx,ERR,ERROR,*999)
          DOMAIN=>DECOMPOSITION%DOMAIN(component_idx)%PTR
          DOMAIN_LINE=>DOMAIN%TOPOLOGY%LINES%LINES(local_line_idx)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis user number = ",DOMAIN_LINE%BASIS%USER_NUMBER, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis family number = ",DOMAIN_LINE%BASIS%FAMILY_NUMBER, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis interpolation type = ",DOMAIN_LINE%BASIS% &
            & INTERPOLATION_TYPE(1),ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis interpolation order = ",DOMAIN_LINE%BASIS% &
            & INTERPOLATION_ORDER(1),ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of nodes in lines = ",DOMAIN_LINE%BASIS%NUMBER_OF_NODES, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_LINE%BASIS%NUMBER_OF_NODES,4,4,DOMAIN_LINE%NODES_IN_LINE, &
            & '("        Nodes in line        :",4(X,I8))','(30X,4(X,I8))',ERR,ERROR,*999)
          DO basis_local_line_node_idx=1,DOMAIN_LINE%BASIS%NUMBER_OF_NODES
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Node : ",basis_local_line_node_idx,ERR,ERROR,*999)
            !/TODO::Loop over local_derivative index so this output makes more sense !<DERIVATIVES_IN_LINE(i,local_derivative_idx,local_node_idx)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
              & DOMAIN_LINE%BASIS%NUMBER_OF_DERIVATIVES(basis_local_line_node_idx),4,4, &
              & DOMAIN_LINE%DERIVATIVES_IN_LINE(1,:,basis_local_line_node_idx),'("            Derivatives in line  :",4(X,I8))', &
              & '(34X,4(X,I8))',ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
              & DOMAIN_LINE%BASIS%NUMBER_OF_DERIVATIVES(basis_local_line_node_idx),4,4, &
              & DOMAIN_LINE%DERIVATIVES_IN_LINE(2,:,basis_local_line_node_idx), &
              & '("            Derivatives Versions in line  :",4(X,I8))','(34X,4(X,I8))',ERR,ERROR,*999)
          ENDDO !basis_local_line_node_idx
        ENDDO !component_idx
      ENDDO !local_line_idx
    ENDIF
    
    EXITS("DECOMPOSITION_TOPOLOGY_LINES_CALCULATE")
    RETURN
999 IF(ASSOCIATED(TEMP_LINES)) DEALLOCATE(TEMP_LINES)
    IF(ASSOCIATED(NEW_TEMP_LINES)) DEALLOCATE(NEW_TEMP_LINES)
    IF(ALLOCATED(NODES_NUMBER_OF_LINES)) DEALLOCATE(NODES_NUMBER_OF_LINES)
    ERRORSEXITS("DECOMPOSITION_TOPOLOGY_LINES_CALCULATE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finalises the lines in the given decomposition topology. \todo Pass in the topology lines
  SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nl
    
    ENTERS("DECOMPOSITION_TOPOLOGY_LINES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%LINES)) THEN
        DO nl=1,TOPOLOGY%LINES%NUMBER_OF_LINES
          CALL DECOMPOSITION_TOPOLOGY_LINE_FINALISE(TOPOLOGY%LINES%LINES(nl),ERR,ERROR,*999)
        ENDDO !nl
        IF(ALLOCATED(TOPOLOGY%LINES%LINES)) DEALLOCATE(TOPOLOGY%LINES%LINES)
        DEALLOCATE(TOPOLOGY%LINES)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DECOMPOSITION_TOPOLOGY_LINES_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_LINES_FINALISE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the line data structures for a decomposition topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to initialise the lines for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_LINES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%LINES)) THEN
        CALL FlagError("Decomposition already has topology lines associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%LINES,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology lines.",ERR,ERROR,*999)
        TOPOLOGY%LINES%NUMBER_OF_LINES=0
        TOPOLOGY%LINES%DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TOPOLOGY_LINES_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_LINES_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the line data structures for a decomposition topology.
  SUBROUTINE DecompositionTopology_DataPointsInitialise(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to initialise the lines for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DecompositionTopology_DataPointsInitialise",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%dataPoints)) THEN
        CALL FlagError("Decomposition already has topology data points associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%dataPoints,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology data points.",ERR,ERROR,*999)
        TOPOLOGY%dataPoints%numberOfDataPoints=0
        TOPOLOGY%dataPoints%totalNumberOfDataPoints=0
        TOPOLOGY%dataPoints%numberOfGlobalDataPoints=0
        NULLIFY(TOPOLOGY%dataPoints%dataPointsTree)
        TOPOLOGY%dataPoints%DECOMPOSITION=>TOPOLOGY%DECOMPOSITION      
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DecompositionTopology_DataPointsInitialise")
    RETURN
999 ERRORSEXITS("DecompositionTopology_DataPointsInitialise",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DecompositionTopology_DataPointsInitialise
  
  !
  !================================================================================================================================
  !
 
  !>Finalises a face in the given decomposition topology and deallocates all memory.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_FACE_FINALISE(FACE,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_FACE_TYPE) :: FACE !<The decomposition face to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_FACE_FINALISE",ERR,ERROR,*999)

    FACE%NUMBER=0
    FACE%XI_DIRECTION=0
    FACE%NUMBER_OF_SURROUNDING_ELEMENTS=0
    IF(ALLOCATED(FACE%SURROUNDING_ELEMENTS)) DEALLOCATE(FACE%SURROUNDING_ELEMENTS)
    IF(ALLOCATED(FACE%ELEMENT_FACES)) DEALLOCATE(FACE%ELEMENT_FACES)
!    FACE%ADJACENT_FACES=0
 
    EXITS("DECOMPOSITION_TOPOLOGY_FACE_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_FACE_FINALISE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_FACE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the face data structure for a decomposition topology face.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_FACE_INITIALISE(FACE,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_FACE_TYPE) :: FACE !<The decomposition face to initialise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_FACE_INITIALISE",ERR,ERROR,*999)

    FACE%NUMBER=0
    FACE%XI_DIRECTION=0
    FACE%NUMBER_OF_SURROUNDING_ELEMENTS=0
!    FACE%ADJACENT_FACES=0
    FACE%BOUNDARY_FACE=.FALSE.
    
    EXITS("DECOMPOSITION_TOPOLOGY_FACE_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_FACE_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_FACE_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Calculates the faces in the given decomposition topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_FACES_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to calculate the faces for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,ne,surrounding_element_idx,basis_local_face_idx,surrounding_element_basis_local_face_idx, &
      & element_local_node_idx,basis_local_face_node_idx,basis_local_face_derivative_idx,derivative_idx,version_idx,face_idx, &
      & node_idx,elem_idx,NODES_IN_FACE(16),NUMBER_OF_FACES,MAX_NUMBER_OF_FACES,NEW_MAX_NUMBER_OF_FACES,FACE_NUMBER
    INTEGER(INTG), ALLOCATABLE :: NODES_NUMBER_OF_FACES(:)
    INTEGER(INTG), POINTER :: TEMP_FACES(:,:),NEW_TEMP_FACES(:,:)
    LOGICAL :: FOUND
    TYPE(BASIS_TYPE), POINTER :: BASIS,BASIS2
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DECOMPOSITION_ELEMENT_TYPE), POINTER :: DECOMPOSITION_ELEMENT
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: DECOMPOSITION_ELEMENTS
    TYPE(DECOMPOSITION_FACE_TYPE), POINTER :: DECOMPOSITION_FACE!,DECOMPOSITION_FACE2
    TYPE(DECOMPOSITION_FACES_TYPE), POINTER :: DECOMPOSITION_FACES
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_ELEMENT_TYPE), POINTER :: DOMAIN_ELEMENT
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS
    TYPE(DOMAIN_FACE_TYPE), POINTER :: DOMAIN_FACE!,DOMAIN_FACE2
    TYPE(DOMAIN_FACES_TYPE), POINTER :: DOMAIN_FACES
    TYPE(DOMAIN_NODE_TYPE), POINTER :: DOMAIN_NODE
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: DOMAIN_TOPOLOGY
    TYPE(MESH_TYPE), POINTER :: MESH

    NULLIFY(TEMP_FACES)
    NULLIFY(NEW_TEMP_FACES)
    
    ENTERS("DECOMPOSITION_TOPOLOGY_FACES_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      DECOMPOSITION_FACES=>TOPOLOGY%FACES
      IF(ASSOCIATED(DECOMPOSITION_FACES)) THEN
        DECOMPOSITION_ELEMENTS=>TOPOLOGY%ELEMENTS
        IF(ASSOCIATED(DECOMPOSITION_ELEMENTS)) THEN
          DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
          IF(ASSOCIATED(DECOMPOSITION)) THEN
            !Process the mesh component number (component number the decomposition was calculated from) first to establish face
            !topology then process the other mesh components.
            DOMAIN=>DECOMPOSITION%DOMAIN(DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
            IF(ASSOCIATED(DOMAIN)) THEN
              DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
              IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
                DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
                IF(ASSOCIATED(DOMAIN_NODES)) THEN
                  DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS
                  IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                    !Estimate the number of faces
                    SELECT CASE(DOMAIN%NUMBER_OF_DIMENSIONS)
                    CASE(1)
                      ! Faces not calculated in 1D 
                    CASE(2)
                      ! Faces not calculated in 2D
                    CASE(3)
                      !This should give the maximum and will over estimate the number of faces for a "cube mesh" by approx 33%
                      MAX_NUMBER_OF_FACES= &
                        & NINT(((REAL(DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS,DP)*5.0_DP)+1.0_DP)*(4.0_DP/3.0_DP),INTG)

                      DOMAIN_FACES=>DOMAIN_TOPOLOGY%FACES
                      IF(ASSOCIATED(DOMAIN_FACES)) THEN
                        ALLOCATE(TEMP_FACES(16,MAX_NUMBER_OF_FACES),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate temporary faces array",ERR,ERROR,*999)
                        ALLOCATE(NODES_NUMBER_OF_FACES(DOMAIN_NODES%TOTAL_NUMBER_OF_NODES),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate nodes number of faces array",ERR,ERROR,*999)
                        NODES_NUMBER_OF_FACES=0
                        NUMBER_OF_FACES=0
                        TEMP_FACES=0
                        !Loop over the elements in the topology and fill temp_faces with node numbers for each element
                        DO ne=1,DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                          DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                          DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(ne)
                          BASIS=>DOMAIN_ELEMENT%BASIS
                          ALLOCATE(DECOMPOSITION_ELEMENT%ELEMENT_FACES(BASIS%NUMBER_OF_LOCAL_FACES),STAT=ERR)
                          IF(ERR/=0) CALL FlagError("Could not allocate element faces of element",ERR,ERROR,*999)
                          !Loop over the local faces of the element
                          DO basis_local_face_idx=1,BASIS%NUMBER_OF_LOCAL_FACES
                            !Calculate the topology node numbers that make up the face
                            NODES_IN_FACE=0
                            !Check whether face has already been read out
                            DO basis_local_face_node_idx=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx)
                              !Read out node numbers of local face from ELEMENT_NODES
                              NODES_IN_FACE(basis_local_face_node_idx)=DOMAIN_ELEMENT%ELEMENT_NODES( &
                                & BASIS%NODE_NUMBERS_IN_LOCAL_FACE(basis_local_face_node_idx,basis_local_face_idx))
                            ENDDO !basis_local_face_node_idx
                            !Try and find a previously created face that matches in the adjacent elements
                            FOUND=.FALSE.
                            node_idx=NODES_IN_FACE(1)
                            DO elem_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_SURROUNDING_ELEMENTS
                              surrounding_element_idx=DOMAIN_NODES%NODES(node_idx)%SURROUNDING_ELEMENTS(elem_idx)
                              IF(surrounding_element_idx/=ne) THEN
                                IF(ALLOCATED(DECOMPOSITION_ELEMENTS%ELEMENTS(surrounding_element_idx)%ELEMENT_FACES)) THEN
                                  BASIS2=>DOMAIN_ELEMENTS%ELEMENTS(surrounding_element_idx)%BASIS
                                  DO surrounding_element_basis_local_face_idx=1,BASIS2%NUMBER_OF_LOCAL_FACES
                                    face_idx=DECOMPOSITION_ELEMENTS%ELEMENTS(surrounding_element_idx)%ELEMENT_FACES( &
                                      & surrounding_element_basis_local_face_idx)
                                    IF(ALL(NODES_IN_FACE(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx))== &
                                      & TEMP_FACES(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx),face_idx))) THEN
                                      FOUND=.TRUE.
                                      EXIT
                                    ENDIF
                                  ENDDO !surrounding_element_basis_local_face_idx
                                  IF(FOUND) EXIT
                                ENDIF
                              ENDIF
                            ENDDO !elem_idx
                            IF(FOUND) THEN
                              !Face has already been created
                              DECOMPOSITION_ELEMENT%ELEMENT_FACES(basis_local_face_idx)=face_idx
                            ELSE
                              !Face has not been created
                              IF(NUMBER_OF_FACES==MAX_NUMBER_OF_FACES) THEN
                                !We are at maximum. Reallocate the FACES array to be 20% bigger and try again.
                                NEW_MAX_NUMBER_OF_FACES=NINT(1.20_DP*REAL(MAX_NUMBER_OF_FACES,DP),INTG)
                                !\todo: Change 16 to a variable and above for NODES_IN_FACE
                                ALLOCATE(NEW_TEMP_FACES(16,NEW_MAX_NUMBER_OF_FACES),STAT=ERR)
                                IF(ERR/=0) CALL FlagError("Could not allocate new number of faces",ERR,ERROR,*999)
                                NEW_TEMP_FACES(:,1:NUMBER_OF_FACES)=TEMP_FACES(:,1:NUMBER_OF_FACES)
                                NEW_TEMP_FACES(:,NUMBER_OF_FACES+1:NEW_MAX_NUMBER_OF_FACES)=0
                                DEALLOCATE(TEMP_FACES)
                                TEMP_FACES=>NEW_TEMP_FACES
                                NULLIFY(NEW_TEMP_FACES)
                                MAX_NUMBER_OF_FACES=NEW_MAX_NUMBER_OF_FACES
                              ENDIF
                              NUMBER_OF_FACES=NUMBER_OF_FACES+1
                              TEMP_FACES(:,NUMBER_OF_FACES)=NODES_IN_FACE(:)
                              DECOMPOSITION_ELEMENT%ELEMENT_FACES(basis_local_face_idx)=NUMBER_OF_FACES
                              DO basis_local_face_node_idx=1,SIZE(NODES_IN_FACE,1)
                                IF(NODES_IN_FACE(basis_local_face_node_idx)/=0) &
                                  & NODES_NUMBER_OF_FACES(NODES_IN_FACE(basis_local_face_node_idx))= &
                                  & NODES_NUMBER_OF_FACES(NODES_IN_FACE(basis_local_face_node_idx))+1
                              ENDDO !basis_local_face_node_idx
                            ENDIF
                          ENDDO !basis_local_face_idx
                        ENDDO !ne

                        !Allocate the face arrays and set them from the FACES and NODE_FACES arrays
                        DO node_idx=1,DOMAIN_NODES%TOTAL_NUMBER_OF_NODES
                          ALLOCATE(DOMAIN_NODES%NODES(node_idx)%NODE_FACES(NODES_NUMBER_OF_FACES(node_idx)),STAT=ERR)
                          IF(ERR/=0) CALL FlagError("Could not allocate node faces array",ERR,ERROR,*999)
                          DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_FACES=0
                        ENDDO !node_idx
                        DEALLOCATE(NODES_NUMBER_OF_FACES)
                        ALLOCATE(DECOMPOSITION_FACES%FACES(NUMBER_OF_FACES),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate decomposition topology faces",ERR,ERROR,*999)
                        DECOMPOSITION_FACES%NUMBER_OF_FACES=NUMBER_OF_FACES
                        ALLOCATE(DOMAIN_FACES%FACES(NUMBER_OF_FACES),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate domain topology faces",ERR,ERROR,*999)
                        DOMAIN_FACES%NUMBER_OF_FACES=NUMBER_OF_FACES
                        DO face_idx=1,DOMAIN_FACES%NUMBER_OF_FACES
                          CALL DECOMPOSITION_TOPOLOGY_FACE_INITIALISE(DECOMPOSITION_FACES%FACES(face_idx),ERR,ERROR,*999)
                          CALL DOMAIN_TOPOLOGY_FACE_INITIALISE(DOMAIN_FACES%FACES(face_idx),ERR,ERROR,*999)
                          DO basis_local_face_node_idx=1,SIZE(TEMP_FACES,1)
                            IF(TEMP_FACES(basis_local_face_node_idx,face_idx)/=0) THEN
                              node_idx=TEMP_FACES(basis_local_face_node_idx,face_idx)
                              DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_FACES=DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_FACES+1
                              DOMAIN_NODES%NODES(node_idx)%NODE_FACES(DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_FACES)=face_idx
                            ENDIF
                          ENDDO !basis_local_face_node_idx  
                        ENDDO !face_idx

                        !Set nodes in face and derivatives of nodes in face for domain faces
                        DO ne=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                          DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(ne)
                          DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                          BASIS=>DOMAIN_ELEMENT%BASIS
                          !Loop over local faces of element
                          DO basis_local_face_idx=1,BASIS%NUMBER_OF_LOCAL_FACES
                            FACE_NUMBER=DECOMPOSITION_ELEMENT%ELEMENT_FACES(basis_local_face_idx)
                            DECOMPOSITION_FACE=>DECOMPOSITION_FACES%FACES(FACE_NUMBER)
                            DOMAIN_FACE=>DOMAIN_FACES%FACES(FACE_NUMBER)
                            DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS=DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS+1
                            IF(.NOT.ASSOCIATED(DOMAIN_FACE%BASIS)) THEN
                              DECOMPOSITION_FACE%NUMBER=FACE_NUMBER
                              DOMAIN_FACE%NUMBER=FACE_NUMBER
                              DOMAIN_FACE%ELEMENT_NUMBER=ne !! Needs checking
!                              DECOMPOSITION_FACE%ELEMENT_NUMBER=DECOMPOSITION_ELEMENT%NUMBER
!                              DOMAIN_FACE%ELEMENT_NUMBER=DOMAIN_ELEMENT%NUMBER
                              DECOMPOSITION_FACE%XI_DIRECTION=BASIS%localFaceXiNormal(basis_local_face_idx)
                              IF(ALLOCATED(BASIS%localFaceBasis)) THEN
                                DOMAIN_FACE%BASIS=>BASIS%localFaceBasis(basis_local_face_idx)%PTR
                              ELSE
                                !Basis is only 2D
                                DOMAIN_FACE%BASIS=>BASIS
                              ENDIF
                              ALLOCATE(DOMAIN_FACE%NODES_IN_FACE(BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx)), &
                                & STAT=ERR)
                              IF(ERR/=0) CALL FlagError("Could not allocate face nodes in face",ERR,ERROR,*999)
                              ALLOCATE(DOMAIN_FACE%DERIVATIVES_IN_FACE(2,DOMAIN_FACE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                                & BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx)),STAT=ERR)
                              IF(ERR/=0) CALL FlagError("Could not allocate face derivatives in face",ERR,ERROR,*999)
                              DOMAIN_FACE%DERIVATIVES_IN_FACE=0
                              !Set nodes in face based upon face number
                              DOMAIN_FACE%NODES_IN_FACE(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx))= &
                                & TEMP_FACES(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx),FACE_NUMBER)
                              !Set derivatives of nodes in domain face from derivatives of nodes in element
                              DO basis_local_face_node_idx=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx)
                                element_local_node_idx=BASIS%NODE_NUMBERS_IN_LOCAL_FACE(basis_local_face_node_idx, &
                                  & basis_local_face_idx)
                                !Set derivative number of u (NO_GLOBAL_DERIV) for the domain face
                                DOMAIN_FACE%DERIVATIVES_IN_FACE(1,1,basis_local_face_node_idx)=NO_GLOBAL_DERIV
                                !Set version number of u (NO_GLOBAL_DERIV) for the domain face
                                version_idx=DOMAIN_ELEMENT%elementVersions(1,BASIS%NODE_NUMBERS_IN_LOCAL_FACE( &
                                  & basis_local_face_node_idx,basis_local_face_idx))
                                DOMAIN_FACE%DERIVATIVES_IN_FACE(2,1,basis_local_face_node_idx)=version_idx
                                IF(DOMAIN_FACE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES>1) THEN
                                  DO basis_local_face_derivative_idx=2,DOMAIN_FACE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES
                                    derivative_idx=DOMAIN_ELEMENT%ELEMENT_DERIVATIVES(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE( &
                                      & basis_local_face_derivative_idx,basis_local_face_node_idx,basis_local_face_idx), &
                                      & element_local_node_idx)
                                    DOMAIN_FACE%DERIVATIVES_IN_FACE(1,basis_local_face_derivative_idx, &
                                      & basis_local_face_node_idx)=derivative_idx
                                    version_idx=DOMAIN_ELEMENT%elementVersions(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE( &
                                      & basis_local_face_derivative_idx,basis_local_face_node_idx,basis_local_face_idx), &
                                      & element_local_node_idx)
                                    DOMAIN_FACE%DERIVATIVES_IN_FACE(2,basis_local_face_derivative_idx, &
                                      & basis_local_face_node_idx)=version_idx
                                  ENDDO !basis_local_face_derivative_idx
                                ENDIF
                              ENDDO !basis_local_face_node_idx
                            ENDIF
                          ENDDO !basis_local_face_idx
                        ENDDO !ne

                        DEALLOCATE(TEMP_FACES)
                        !\todo Note: Adjacency will be left out of faces calculation for the time being
                        !Calculate adjacent faces and the surrounding elements for each face
                        DO face_idx=1,DECOMPOSITION_FACES%NUMBER_OF_FACES
                          DECOMPOSITION_FACE=>DECOMPOSITION_FACES%FACES(face_idx)
                          DOMAIN_FACE=>DOMAIN_FACES%FACES(face_idx)
                          BASIS=>DOMAIN_FACE%BASIS
                          IF(DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS==1) THEN
                            DECOMPOSITION_FACE%BOUNDARY_FACE=.TRUE.
                            DOMAIN_FACE%BOUNDARY_FACE=.TRUE.
                          ENDIF
                          !Allocate the elements surrounding the face
                          ALLOCATE(DECOMPOSITION_FACE%SURROUNDING_ELEMENTS(DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS), &
                            & STAT=ERR)
                          IF(ERR/=0) CALL FlagError("Could not allocate face surrounding elements",ERR,ERROR,*999)

                          ALLOCATE(DECOMPOSITION_FACE%ELEMENT_FACES(DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS), &
                            & STAT=ERR)
                          IF(ERR/=0) CALL FlagError("Could not allocate face element faces",ERR,ERROR,*999)
!                          DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS=0
!                          DECOMPOSITION_FACE%ADJACENT_FACES=0

                           !Loop over the nodes at each end of the face
!                          DO node_idx1=0,1
!                           DO node_idx2=0,1
!                            FOUND=.FALSE.
!                            node_idx=DOMAIN_FACE%NODES_IN_FACE((node_idx2*BASIS%NUMBER_OF_NODES_IN_XI_DIRECTION*(BASIS%NUMBER_OF_FACES-1))&
!                                                                             &+(node_idx1*(BASIS%NUMBER_OF_NODES_IN_XI_DIRECTION-1))+1)
                             !Loop over the elements surrounding the node.
!                            DO elem_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_SURROUNDING_ELEMENTS
!                              ne=DOMAIN_NODES%NODES(node_idx)%SURROUNDING_ELEMENTS(elem_idx)
!                              DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(ne)
!                              DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                               !Loop over the local faces of the element
!                              DO basis_local_face_idx=1,DOMAIN_ELEMENT%BASIS%NUMBER_OF_LOCAL_FACES
!                                nf2=DECOMPOSITION_ELEMENT%ELEMENT_FACES(basis_local_face_idx)
!                                IF(nf2/=face_idx) THEN
!                                  DECOMPOSITION_FACE2=>DECOMPOSITION_FACES%FACES(nf2)
!                                  DOMAIN_FACE2=>DOMAIN_FACES%FACES(nf2)
                                   !Check whether XI of face have same direction
!                                  IF ((OTHER_XI_DIRECTIONS3(BASIS%localFaceXiNormal(basis_local_face_idx),2,1)==&
!                                     &OTHER_XI_DIRECTIONS3(BASIS2%localFaceXiNormal(basis_local_face_idx),2,1)).OR.&
!                                     &(OTHER_XI_DIRECTIONS3(BASIS%localFaceXiNormal(basis_local_face_idx),3,1)==&
!                                     &OTHER_XI_DIRECTIONS3(BASIS2%localFaceXiNormal(basis_local_face_idx),3,1))) THEN
                                     !Loop over nodes in face of surrounding element
!                                    BASIS2=>DOMAIN_FACE2%BASIS
!                                    IF(BASIS2%INTERPOLATION_ORDER(1)==BASIS%INTERPOLATION_ORDER(1)) THEN
!                                      NODE_COUNT=0
!                                      DO node_idx3=1,BASIS%NUMBER_OF_NODES_IN_XI_DIRECTION
!                                        DO node_idx4=1,BASIS%NUMBER_OF_NODES_IN_XI_DIRECTION
!                                          np2=DOMAIN_FACE2%NODES_IN_FACE((node_idx4*(BASIS2%NUMBER_OF_FACES-1))&
!                                                                      &+(node_idx3*(BASIS2%NUMBER_OF_NODES_IN_XI_DIRECTION-1))+1)
!                                          IF(np2==node_idx) NODE_COUNT=NODE_COUNT+1
!                                        ENDDO !node_idx4
!                                      ENDDO !node_idx3
!                                      IF(NODE_COUNT<BASIS%NUMBER_OF_NODES) THEN
!                                        FOUND=.TRUE.
!                                        EXIT
!                                      ENDIF
!                                    ENDIF
!                                  ENDIF
!                                ENDIF
!                              ENDDO !basis_local_face_idx
!                                IF(FOUND) EXIT
!                            ENDDO !elem_idx
!                            IF(FOUND) DECOMPOSITION_FACE%ADJACENT_FACES(node_idx2)=nf2
!                           ENDDO !node_idx2
!                           IF(FOUND) DECOMPOSITION_FACE%ADJACENT_FACES(node_idx1)=nf2
!                          ENDDO !node_idx1
                        ENDDO !face_idx

                        !Set the surrounding elements
                        DO ne=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                          DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(ne)
                          DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                          BASIS=>DOMAIN_ELEMENT%BASIS
                          DO basis_local_face_idx=1,BASIS%NUMBER_OF_LOCAL_FACES
                            FACE_NUMBER=DECOMPOSITION_ELEMENT%ELEMENT_FACES(basis_local_face_idx)
                            DECOMPOSITION_FACE=>DECOMPOSITION_FACES%FACES(FACE_NUMBER)
                            DO face_idx=1,DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS
                              DECOMPOSITION_FACE%SURROUNDING_ELEMENTS(face_idx)=ne
                              DECOMPOSITION_FACE%ELEMENT_FACES(face_idx)=basis_local_face_idx
                            ENDDO
                          ENDDO !basis_local_face_idx
                        ENDDO !ne
                      ELSE
                        CALL FlagError("Domain topology faces is not associated",ERR,ERROR,*999)
                      ENDIF
                    CASE DEFAULT
                      CALL FlagError("Invalid number of dimensions for a topology domain",ERR,ERROR,*999)
                    END SELECT
                 ELSE
                    CALL FlagError("Domain topology elements is not associated",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Domain topology nodes is not associated",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Topology decomposition domain topology is not associated",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Topology decomposition domain is not associated",ERR,ERROR,*999)
            ENDIF
            !Now loop over the other mesh components in the decomposition and calculate the domain faces
            MESH=>DECOMPOSITION%MESH
            IF(ASSOCIATED(MESH)) THEN
              DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
                IF(component_idx/=DECOMPOSITION%MESH_COMPONENT_NUMBER) THEN
                  DOMAIN=>DECOMPOSITION%DOMAIN(component_idx)%PTR
                  IF(ASSOCIATED(DOMAIN)) THEN
                    DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
                    IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
                      DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
                      IF(ASSOCIATED(DOMAIN_NODES)) THEN
                        DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS
                        IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                          DOMAIN_FACES=>DOMAIN_TOPOLOGY%FACES
                          IF(ASSOCIATED(DOMAIN_FACES)) THEN
                            ALLOCATE(DOMAIN_FACES%FACES(DECOMPOSITION_FACES%NUMBER_OF_FACES),STAT=ERR)
                            IF(ERR/=0) CALL FlagError("Could not allocate domain faces faces",ERR,ERROR,*999)
                            DOMAIN_FACES%NUMBER_OF_FACES=DECOMPOSITION_FACES%NUMBER_OF_FACES
                            ALLOCATE(NODES_NUMBER_OF_FACES(DOMAIN_NODES%TOTAL_NUMBER_OF_NODES),STAT=ERR)
                            IF(ERR/=0) CALL FlagError("Could not allocate nodes number of faces array",ERR,ERROR,*999)
                            NODES_NUMBER_OF_FACES=0
                            !Loop over the faces in the topology
                            DO face_idx=1,DECOMPOSITION_FACES%NUMBER_OF_FACES
                              DECOMPOSITION_FACE=>DECOMPOSITION_FACES%FACES(face_idx)
                              DOMAIN_FACE=>DOMAIN_FACES%FACES(face_idx)
                              IF(DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS>0) THEN
                                ne=DECOMPOSITION_FACE%SURROUNDING_ELEMENTS(1)
                                basis_local_face_idx=DECOMPOSITION_FACE%ELEMENT_FACES(1)
                                CALL DOMAIN_TOPOLOGY_FACE_INITIALISE(DOMAIN_FACES%FACES(face_idx),ERR,ERROR,*999)
                                DOMAIN_FACE%NUMBER=face_idx
                                DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                                BASIS=>DOMAIN_ELEMENT%BASIS
                                IF(ALLOCATED(BASIS%localFaceBasis)) THEN
                                  DOMAIN_FACE%BASIS=>BASIS%localFaceBasis(basis_local_face_idx)%PTR
                                ELSE
                                  !Basis is only 2D
                                  DOMAIN_FACE%BASIS=>BASIS
                                ENDIF
                                ALLOCATE(DOMAIN_FACE%NODES_IN_FACE(BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx)), &
                                  & STAT=ERR)
                                IF(ERR/=0) CALL FlagError("Could not allocate nodes in face",ERR,ERROR,*999)
                                ALLOCATE(DOMAIN_FACE%DERIVATIVES_IN_FACE(2,DOMAIN_FACE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                                  & BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx)),STAT=ERR)
                                IF(ERR/=0) CALL FlagError("Could not allocate derivatives in face",ERR,ERROR,*999)
                                !Set derivatives of nodes in domain face from derivatives of nodes in element
                                DO basis_local_face_node_idx=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx)
                                  element_local_node_idx=BASIS%NODE_NUMBERS_IN_LOCAL_FACE(basis_local_face_node_idx, &
                                    & basis_local_face_idx)
                                  node_idx=DOMAIN_ELEMENT%ELEMENT_NODES(element_local_node_idx)
                                  DOMAIN_FACE%NODES_IN_FACE(basis_local_face_node_idx)=node_idx
                                  !Set derivative number of u (NO_GLOBAL_DERIV) for the domain face
                                  DOMAIN_FACE%DERIVATIVES_IN_FACE(1,1,basis_local_face_node_idx)=NO_GLOBAL_DERIV
                                  !Set version number of u (NO_GLOBAL_DERIV) for the domain face
                                  version_idx=DOMAIN_ELEMENT%elementVersions(1,BASIS%NODE_NUMBERS_IN_LOCAL_FACE( &
                                    & basis_local_face_node_idx,basis_local_face_idx))
                                  DOMAIN_FACE%DERIVATIVES_IN_FACE(2,1,basis_local_face_node_idx)=version_idx
                                  IF(DOMAIN_FACE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES>1) THEN
                                    DO basis_local_face_derivative_idx=2,DOMAIN_FACE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES
                                      derivative_idx=DOMAIN_ELEMENT%ELEMENT_DERIVATIVES(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE( &
                                        & basis_local_face_derivative_idx,basis_local_face_node_idx,basis_local_face_idx), &
                                        & element_local_node_idx)
                                      DOMAIN_FACE%DERIVATIVES_IN_FACE(1,basis_local_face_derivative_idx, &
                                        & basis_local_face_node_idx)=derivative_idx
                                      version_idx=DOMAIN_ELEMENT%elementVersions(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE( &
                                        & basis_local_face_derivative_idx,basis_local_face_node_idx,basis_local_face_idx), &
                                        & element_local_node_idx)
                                      DOMAIN_FACE%DERIVATIVES_IN_FACE(2,basis_local_face_derivative_idx, &
                                        & basis_local_face_node_idx)=version_idx
                                    ENDDO !basis_local_face_derivative_idx
                                  ENDIF
                                  NODES_NUMBER_OF_FACES(node_idx)=NODES_NUMBER_OF_FACES(node_idx)+1
                                ENDDO !basis_local_face_node_idx
                              ELSE
                                CALL FlagError("Face is not surrounded by any elements?",ERR,ERROR,*999)
                              ENDIF
                            ENDDO !face_idx
                            DO node_idx=1,DOMAIN_NODES%TOTAL_NUMBER_OF_NODES
                              ALLOCATE(DOMAIN_NODES%NODES(node_idx)%NODE_FACES(NODES_NUMBER_OF_FACES(node_idx)),STAT=ERR)
                              IF(ERR/=0) CALL FlagError("Could not allocate node faces",ERR,ERROR,*999)
                              DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_FACES=0
                            ENDDO !node_idx
                            DEALLOCATE(NODES_NUMBER_OF_FACES)
                            DO face_idx=1,DOMAIN_FACES%NUMBER_OF_FACES
                              DOMAIN_FACE=>DOMAIN_FACES%FACES(face_idx)
                              BASIS=>DOMAIN_FACE%BASIS
                              DO basis_local_face_node_idx=1,BASIS%NUMBER_OF_NODES
                                node_idx=DOMAIN_FACE%NODES_IN_FACE(basis_local_face_node_idx)
                                DOMAIN_NODE=>DOMAIN_NODES%NODES(node_idx)
                                DOMAIN_NODE%NUMBER_OF_NODE_FACES=DOMAIN_NODE%NUMBER_OF_NODE_FACES+1
                                !Set the face numbers a node is on
                                DOMAIN_NODE%NODE_FACES(DOMAIN_NODE%NUMBER_OF_NODE_FACES)=face_idx
                              ENDDO !basis_local_face_node_idx
                            ENDDO !face_idx
                          ELSE
                            CALL FlagError("Domain faces is not associated",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FlagError("Domain elements is not associated",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Domain nodes is not associated",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Domain topology is not associated",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Decomposition mesh is not associated",ERR,ERROR,*999)
                  ENDIF
                ENDIF
              ENDDO !component_idx
            ELSE
              CALL FlagError("Decomposition mesh is not associated",ERR,ERROR,*999)
            ENDIF                        
          ELSE
            CALL FlagError("Topology decomposition is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Topology decomposition elements is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Topology faces is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Decomposition topology faces:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of mesh components = ",MESH%NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of faces = ",DECOMPOSITION_FACES%NUMBER_OF_FACES,ERR,ERROR,*999)
      DO face_idx=1,DECOMPOSITION_FACES%NUMBER_OF_FACES
        DECOMPOSITION_FACE=>DECOMPOSITION_FACES%FACES(face_idx)
        DOMAIN_FACE=>DOMAIN_FACES%FACES(face_idx)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Face number = ",DECOMPOSITION_FACE%NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Xi direction (Normal to Face) = &
                                                                         &",DECOMPOSITION_FACE%XI_DIRECTION,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of surrounding elements = ", &
          & DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS,4,4, &
          & DECOMPOSITION_FACE%SURROUNDING_ELEMENTS,'("      Surrounding elements :",4(X,I8))','(28X,4(X,I8))',ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS,4,4, &
          & DECOMPOSITION_FACE%ELEMENT_FACES,'("      Element faces        :",4(X,I8))','(28X,4(X,I8))',ERR,ERROR,*999)
!        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,2,2,2,DECOMPOSITION_FACE%ADJACENT_FACES, &
!          & '("      Adjacent faces       :",2(X,I8))','(28X,2(X,I8))',ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Boundary face = ",DECOMPOSITION_FACE%BOUNDARY_FACE,ERR,ERROR,*999)
        DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Mesh component : ",component_idx,ERR,ERROR,*999)
          DOMAIN=>DECOMPOSITION%DOMAIN(component_idx)%PTR
          DOMAIN_FACE=>DOMAIN%TOPOLOGY%FACES%FACES(face_idx)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis user number = ",DOMAIN_FACE%BASIS%USER_NUMBER, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis family number = ",DOMAIN_FACE%BASIS%FAMILY_NUMBER, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis interpolation type = ",DOMAIN_FACE%BASIS% &
            & INTERPOLATION_TYPE(1),ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis interpolation order = ",DOMAIN_FACE%BASIS% &
            & INTERPOLATION_ORDER(1),ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of nodes in faces = ",DOMAIN_FACE%BASIS%NUMBER_OF_NODES, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_FACE%BASIS%NUMBER_OF_NODES,4,4,DOMAIN_FACE%NODES_IN_FACE, &
            & '("        Nodes in face        :",4(X,I8))','(30X,4(X,I8))',ERR,ERROR,*999)
          DO basis_local_face_node_idx=1,DOMAIN_FACE%BASIS%NUMBER_OF_NODES
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Node : ",basis_local_face_node_idx,ERR,ERROR,*999)
            !/TODO::Loop over local_derivative index so this output makes more sense !<DERIVATIVES_IN_LINE(i,local_derivative_idx,local_node_idx)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
              & DOMAIN_FACE%BASIS%NUMBER_OF_DERIVATIVES(basis_local_face_node_idx),4,4,DOMAIN_FACE% &
              & DERIVATIVES_IN_FACE(1,:,basis_local_face_node_idx),'("            Derivatives in face  :",4(X,I8))', &
              & '(34X,4(X,I8))',ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
              & DOMAIN_FACE%BASIS%NUMBER_OF_DERIVATIVES(basis_local_face_node_idx),4,4,DOMAIN_FACE% &
              & DERIVATIVES_IN_FACE(2,:,basis_local_face_node_idx),'("            Derivatives Versions in face  :",4(X,I8))', &
              & '(34X,4(X,I8))',ERR,ERROR,*999)
          ENDDO !basis_local_face_node_idx
        ENDDO !component_idx
      ENDDO !face_idx
    ENDIF

    EXITS("DECOMPOSITION_TOPOLOGY_FACES_CALCULATE")
    RETURN
999 IF(ASSOCIATED(TEMP_FACES)) DEALLOCATE(TEMP_FACES)
    IF(ASSOCIATED(NEW_TEMP_FACES)) DEALLOCATE(NEW_TEMP_FACES)
    IF(ALLOCATED(NODES_NUMBER_OF_FACES)) DEALLOCATE(NODES_NUMBER_OF_FACES)
    ERRORSEXITS("DECOMPOSITION_TOPOLOGY_FACES_CALCULATE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_FACES_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finalises the faces in the given decomposition topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_FACES_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nf
    
    ENTERS("DECOMPOSITION_TOPOLOGY_FACES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%FACES)) THEN
        DO nf=1,TOPOLOGY%FACES%NUMBER_OF_FACES
          CALL DECOMPOSITION_TOPOLOGY_FACE_FINALISE(TOPOLOGY%FACES%FACES(nf),ERR,ERROR,*999)
        ENDDO !nf
        IF(ALLOCATED(TOPOLOGY%FACES%FACES)) DEALLOCATE(TOPOLOGY%FACES%FACES)
        DEALLOCATE(TOPOLOGY%FACES)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DECOMPOSITION_TOPOLOGY_FACES_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_FACES_FINALISE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_FACES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the face data structures for a decomposition topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_FACES_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to initialise the faces for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TOPOLOGY_FACES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%FACES)) THEN
        CALL FlagError("Decomposition already has topology faces associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%FACES,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology faces",ERR,ERROR,*999)
        TOPOLOGY%FACES%NUMBER_OF_FACES=0
        TOPOLOGY%FACES%DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TOPOLOGY_FACES_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TOPOLOGY_FACES_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_FACES_INITIALISE

  !
  !================================================================================================================================
  !
  
  !>Gets the decomposition type for a decomposition. \see OPENCMISS::Iron::cmfe_DecompositionTypeGet
  SUBROUTINE DECOMPOSITION_TYPE_GET(DECOMPOSITION,TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to get the type for
    INTEGER(INTG), INTENT(OUT) :: TYPE !<On return, the decomposition type for the specified decomposition \see MESH_ROUTINES_DecompositionTypes,MESH_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITION_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        TYPE=DECOMPOSITION%DECOMPOSITION_TYPE
      ELSE
        CALL FlagError("Decomposition has not finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TYPE_GET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TYPE_GET",ERR,ERROR)
    RETURN
  END SUBROUTINE DECOMPOSITION_TYPE_GET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the decomposition type for a decomposition.  \see OPENCMISS::Iron::cmfe_DecompositionTypeSet
  SUBROUTINE DECOMPOSITION_TYPE_SET(DECOMPOSITION,TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to set the type for
    INTEGER(INTG), INTENT(IN) :: TYPE !<The decomposition type to set \see MESH_ROUTINES_DecompositionTypes,MESH_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("DECOMPOSITION_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        CALL FlagError("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(TYPE)
        CASE(DECOMPOSITION_ALL_TYPE)
          !heye: three types for decomposition--decompostion_all_type means no decomposition
          DECOMPOSITION%DECOMPOSITION_TYPE=DECOMPOSITION_ALL_TYPE
        CASE(DECOMPOSITION_CALCULATED_TYPE)
          DECOMPOSITION%DECOMPOSITION_TYPE=DECOMPOSITION_CALCULATED_TYPE
        CASE(DECOMPOSITION_USER_DEFINED_TYPE)
          DECOMPOSITION%DECOMPOSITION_TYPE=DECOMPOSITION_USER_DEFINED_TYPE
        CASE DEFAULT
          LOCAL_ERROR="Decomposition type "//TRIM(NUMBER_TO_VSTRING(TYPE,"*",ERR,ERROR))//" is not valid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_TYPE_SET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_TYPE_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TYPE_SET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes whether lines should be calculated in the the decomposition. \see OPENCMISS::Iron::cmfe_DecompositionCalculateLinesSet
  SUBROUTINE DECOMPOSITION_CALCULATE_LINES_SET(DECOMPOSITION,CALCULATE_LINES_FLAG,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition
    LOGICAL, INTENT(IN) :: CALCULATE_LINES_FLAG !<The boolean flag to determine whether the lines should be calculated or not
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    ENTERS("DECOMPOSITION_CALCULATE_LINES_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        CALL FlagError("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        DECOMPOSITION%CALCULATE_LINES=CALCULATE_LINES_FLAG
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_CALCULATE_LINES_SET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_CALCULATE_LINES_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_CALCULATE_LINES_SET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes whether faces should be calculated in the the decomposition. \see OPENCMISS::Iron::cmfe_DecompositionCalculateFacesSet
  SUBROUTINE DECOMPOSITION_CALCULATE_FACES_SET(DECOMPOSITION,CALCULATE_FACES_FLAG,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition
    LOGICAL, INTENT(IN) :: CALCULATE_FACES_FLAG !<The boolean flag to determine whether the faces should be calculated or not
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    ENTERS("DECOMPOSITION_CALCULATE_FACES_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        CALL FlagError("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        DECOMPOSITION%CALCULATE_FACES=CALCULATE_FACES_FLAG
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITION_CALCULATE_FACES_SET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_CALCULATE_FACES_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_CALCULATE_FACES_SET
  
  !
  !================================================================================================================================
  !

  !>Finalises the domain decompositions for a given mesh. \todo Pass in a pointer to the decompositions.
  SUBROUTINE DECOMPOSITIONS_FINALISE(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to finalise the decomposition for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITIONS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(ASSOCIATED(MESH%DECOMPOSITIONS)) THEN
        DO WHILE(MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS>0)
          CALL DECOMPOSITION_DESTROY(MESH%DECOMPOSITIONS%DECOMPOSITIONS(1)%PTR,ERR,ERROR,*999)
        ENDDO !no_decomposition
       DEALLOCATE(MESH%DECOMPOSITIONS)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITIONS_FINALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITIONS_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITIONS_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the domain decompositions for a given mesh.
  SUBROUTINE DECOMPOSITIONS_INITIALISE(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to initialise the decompositions for

    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DECOMPOSITIONS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(ASSOCIATED(MESH%DECOMPOSITIONS)) THEN
        CALL FlagError("Mesh already has decompositions associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(MESH%DECOMPOSITIONS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Mesh decompositions could not be allocated.",ERR,ERROR,*999)
        MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS=0
        MESH%DECOMPOSITIONS%MESH=>MESH
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DECOMPOSITIONS_INITIALISE")
    RETURN
999 ERRORSEXITS("DECOMPOSITIONS_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITIONS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the domain for a given decomposition and deallocates all memory. \todo Pass in a pointer to the domain
  SUBROUTINE DOMAIN_FINALISE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to finalise the domain for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx
    
    ENTERS("DOMAIN_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(ASSOCIATED(DECOMPOSITION%MESH)) THEN
        IF(ASSOCIATED(DECOMPOSITION%DOMAIN)) THEN
          DO component_idx=1,SIZE(DECOMPOSITION%DOMAIN,1)
            IF(ALLOCATED(DECOMPOSITION%DOMAIN(component_idx)%PTR%NODE_DOMAIN))  &
              & DEALLOCATE(DECOMPOSITION%DOMAIN(component_idx)%PTR%NODE_DOMAIN)
            CALL DOMAIN_MAPPINGS_FINALISE(DECOMPOSITION%DOMAIN(component_idx)%PTR,ERR,ERROR,*999)        
            CALL DOMAIN_TOPOLOGY_FINALISE(DECOMPOSITION%DOMAIN(component_idx)%PTR,ERR,ERROR,*999)
            DEALLOCATE(DECOMPOSITION%DOMAIN(component_idx)%PTR)
          ENDDO !component_idx
          DEALLOCATE(DECOMPOSITION%DOMAIN)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the domain for a given decomposition.
  SUBROUTINE DOMAIN_INITIALISE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to initialise the domain for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx

    ENTERS("DOMAIN_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(ASSOCIATED(DECOMPOSITION%MESH)) THEN
        IF(ASSOCIATED(DECOMPOSITION%DOMAIN)) THEN
          CALL FlagError("Decomposition already has a domain associated.",ERR,ERROR,*999)
        ELSE
          ALLOCATE(DECOMPOSITION%DOMAIN(DECOMPOSITION%numberOfComponents),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Decomposition domain could not be allocated.",ERR,ERROR,*999)
          DO component_idx=1,DECOMPOSITION%numberOfComponents!Mesh component
            ALLOCATE(DECOMPOSITION%DOMAIN(component_idx)%PTR,STAT=ERR)
            IF(ERR/=0) CALL FlagError("Decomposition domain component could not be allocated.",ERR,ERROR,*999)
            DECOMPOSITION%DOMAIN(component_idx)%PTR%DECOMPOSITION=>DECOMPOSITION
            DECOMPOSITION%DOMAIN(component_idx)%PTR%MESH=>DECOMPOSITION%MESH
            DECOMPOSITION%DOMAIN(component_idx)%PTR%MESH_COMPONENT_NUMBER=component_idx
            DECOMPOSITION%DOMAIN(component_idx)%PTR%REGION=>DECOMPOSITION%region
            DECOMPOSITION%DOMAIN(component_idx)%PTR%NUMBER_OF_DIMENSIONS=DECOMPOSITION%numberOfDimensions
            !DECOMPOSITION%DOMAIN(component_idx)%PTR%NUMBER_OF_ELEMENTS=0
            !DECOMPOSITION%DOMAIN(component_idx)%PTR%NUMBER_OF_FACES=0
            !DECOMPOSITION%DOMAIN(component_idx)%PTR%NUMBER_OF_LINES=0
            !DECOMPOSITION%DOMAIN(component_idx)%PTR%NUMBER_OF_NODES=0
            !DECOMPOSITION%DOMAIN(component_idx)%PTR%NUMBER_OF_MESH_DOFS=0
            NULLIFY(DECOMPOSITION%DOMAIN(component_idx)%PTR%MAPPINGS)
            NULLIFY(DECOMPOSITION%DOMAIN(component_idx)%PTR%TOPOLOGY)
            CALL DOMAIN_MAPPINGS_INITIALISE(DECOMPOSITION%DOMAIN(component_idx)%PTR,ERR,ERROR,*999)
            CALL DOMAIN_TOPOLOGY_INITIALISE(DECOMPOSITION%DOMAIN(component_idx)%PTR,ERR,ERROR,*999)
          ENDDO !component_idx
        ENDIF
      ELSE
        CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_INITIALISE",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE DOMAIN_INITIALISE
  
 !
  !================================================================================================================================
  !

  !>Finalises the dofs mapping in the given domain mappings. \todo Pass in the domain mappings dofs
  SUBROUTINE DOMAIN_MAPPINGS_DOFS_FINALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS !<A pointer to the domain mappings to finalise the dofs for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("DOMAIN_MAPPINGS_DOFS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%DOFS)) THEN
        CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(DOMAIN_MAPPINGS%DOFS,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_MAPPINGS_DOFS_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_DOFS_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_DOFS_FINALISE

  !
  !================================================================================================================================
  !

  !>Intialises the dofs mapping in the given domain mapping.
  SUBROUTINE DOMAIN_MAPPINGS_DOFS_INITIALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS !<A pointer to the domain mappings to initialise the dofs for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("DOMAIN_MAPPINGS_DOFS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%DOFS)) THEN
        CALL FlagError("Domain dofs mappings are already associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(DOMAIN_MAPPINGS%DOFS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate domain mappings dofs.",ERR,ERROR,*999)
        CALL DOMAIN_MAPPINGS_MAPPING_INITIALISE(DOMAIN_MAPPINGS%DOFS,DOMAIN_MAPPINGS%DOMAIN%DECOMPOSITION%workGroup, &
          & ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_MAPPINGS_DOFS_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_DOFS_INITIALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_DOFS_INITIALISE

  !
  !================================================================================================================================
  !

  !>Calculates the local/global element mappings for a domain decomposition.
  SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_CALCULATE(DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<A pointer to the domain to calculate the element mappings for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,no_adjacent_element,adjacent_element,domain_no,domain_idx,ne,nn,np,NUMBER_OF_DOMAINS, &
!      & NUMBER_OF_ADJACENT_ELEMENTS,myComputationalNodeNumber,component_idx
       & NUMBER_OF_ADJACENT_ELEMENTS,component_idx 
    INTEGER(INTG), ALLOCATABLE :: ADJACENT_ELEMENTS(:),DOMAINS(:),LOCAL_ELEMENT_NUMBERS(:)
    TYPE(LIST_TYPE), POINTER :: ADJACENT_DOMAINS_LIST
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: ADJACENT_ELEMENTS_LIST(:)
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    ENTERS("DOMAIN_MAPPINGS_ELEMENTS_CALCULATE",ERR,ERROR,*999)

    ! This routine uses the following information:
    !   Adjacent elements stored at nodes
    !   DECOMPOSITION%ELEMENT_DOMAIN(global element number)
    !
    ! It computes and sets the following variables:
    !   ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP
    !   it then calls DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE
    !   to set everything else in ELEMENTS_MAPPING
    !  

    IF(ASSOCIATED(DOMAIN)) THEN
      IF(ASSOCIATED(DOMAIN%MAPPINGS)) THEN
        IF(ASSOCIATED(DOMAIN%MAPPINGS%ELEMENTS)) THEN
          ELEMENTS_MAPPING=>DOMAIN%MAPPINGS%ELEMENTS
          IF(ASSOCIATED(DOMAIN%DECOMPOSITION)) THEN
            DECOMPOSITION=>DOMAIN%DECOMPOSITION
            IF(ASSOCIATED(DOMAIN%MESH)) THEN
              MESH=>DOMAIN%MESH
              component_idx=DOMAIN%MESH_COMPONENT_NUMBER
             
              !Calculate the local and global numbers and set up the mappings
              ALLOCATE(ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(MESH%NUMBER_OF_ELEMENTS),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate element mapping global to local map.",ERR,ERROR,*999)
              ELEMENTS_MAPPING%NUMBER_OF_GLOBAL=MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%NUMBER_OF_ELEMENTS
              !Loop over the global elements and calculate local numbers
              ! LOCAL_ELEMENT_NUMBERS(nd) contains a counter for the local element number on domain nd, it is a temporary array used to create GLOBAL_TO_LOCAL_MAP 
              ALLOCATE(LOCAL_ELEMENT_NUMBERS(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate local element numbers.",ERR,ERROR,*999)
              LOCAL_ELEMENT_NUMBERS=0
              ALLOCATE(ADJACENT_ELEMENTS_LIST(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate adjacent elements list.",ERR,ERROR,*999)
              DO domain_idx=0,DECOMPOSITION%NUMBER_OF_DOMAINS-1
                NULLIFY(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR)
                CALL LIST_CREATE_START(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,ERR,ERROR,*999)
                CALL LIST_DATA_TYPE_SET(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                CALL LIST_INITIAL_SIZE_SET(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,MAX(INT(MESH%NUMBER_OF_ELEMENTS/2),1), &
                  & ERR,ERROR,*999)
                CALL LIST_CREATE_FINISH(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,ERR,ERROR,*999)
              ENDDO !domain_idx
            
              ! loop over global elements 
              DO ne=1,MESH%NUMBER_OF_ELEMENTS
                !Calculate the local numbers (increment counter of respective domain) 
                domain_no=DECOMPOSITION%ELEMENT_DOMAIN(ne)
                LOCAL_ELEMENT_NUMBERS(domain_no)=LOCAL_ELEMENT_NUMBERS(domain_no)+1
                !Calculate the adjacent elements to the computational domains and the adjacent domain numbers themselves
                BASIS=>MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%ELEMENTS(ne)%BASIS
                NULLIFY(ADJACENT_DOMAINS_LIST)
                CALL LIST_CREATE_START(ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                CALL LIST_DATA_TYPE_SET(ADJACENT_DOMAINS_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                CALL LIST_INITIAL_SIZE_SET(ADJACENT_DOMAINS_LIST,DECOMPOSITION%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
                CALL LIST_CREATE_FINISH(ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                CALL LIST_ITEM_ADD(ADJACENT_DOMAINS_LIST,domain_no,ERR,ERROR,*999)

                ! loop over the number of nodes of the element 
                DO nn=1,BASIS%NUMBER_OF_NODES
                  np=MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%ELEMENTS(ne)%MESH_ELEMENT_NODES(nn)
                  DO no_adjacent_element=1,MESH%TOPOLOGY(component_idx)%PTR%NODES%NODES(np)%numberOfSurroundingElements
                    adjacent_element=MESH%TOPOLOGY(component_idx)%PTR%NODES%NODES(np)%surroundingElements(no_adjacent_element)
                    IF(DECOMPOSITION%ELEMENT_DOMAIN(adjacent_element)/=domain_no) THEN
                      CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(domain_no)%PTR,adjacent_element,ERR,ERROR,*999)
                      CALL LIST_ITEM_ADD(ADJACENT_DOMAINS_LIST,DECOMPOSITION%ELEMENT_DOMAIN(adjacent_element),ERR,ERROR,*999)
                    ENDIF
                  ENDDO !no_adjacent_element
                ENDDO !nn
                CALL LIST_REMOVE_DUPLICATES(ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                CALL LIST_DETACH_AND_DESTROY(ADJACENT_DOMAINS_LIST,NUMBER_OF_DOMAINS,DOMAINS,ERR,ERROR,*999)
                DEALLOCATE(DOMAINS)
                ! The previously filled list is destroyed here and the underlying array is deallocated. What remains is the number of entries
                ! (NUMBER_OF_DOMAINS). Setting up the list was necessary to count non-duplicate occurences of the entries.
                ! The entries were the domains on which adjacent elements were located, i.e. now we know how many adjacent domains to the current element exist. 
                CALL DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE(ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne),ERR,ERROR,*999)
                ALLOCATE(ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_NUMBER(NUMBER_OF_DOMAINS),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate element global to local map local number.",ERR,ERROR,*999)
                ALLOCATE(ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%DOMAIN_NUMBER(NUMBER_OF_DOMAINS),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate element global to local map domain number.",ERR,ERROR,*999)
                ALLOCATE(ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_TYPE(NUMBER_OF_DOMAINS),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate element global to local map local type.",ERR,ERROR,*999)
                ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%NUMBER_OF_DOMAINS=1
                ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_NUMBER(1)=LOCAL_ELEMENT_NUMBERS(domain_no)
                ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%DOMAIN_NUMBER(1)=DECOMPOSITION%ELEMENT_DOMAIN(ne) 
                IF(NUMBER_OF_DOMAINS==1) THEN   ! There are no adjacent element on another domain 
                  !Element is an internal element
                  ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_TYPE(1)=DOMAIN_LOCAL_INTERNAL
                ELSE ! There are adjacent elements on other domains
                  !Element is on the boundary of computational domains
                  ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_TYPE(1)=DOMAIN_LOCAL_BOUNDARY
                ENDIF
              ENDDO !ne
              
              !Compute ghost element mappings
              DO domain_idx=0,DECOMPOSITION%NUMBER_OF_DOMAINS-1
                ! get number of adjacent elements 
                CALL LIST_REMOVE_DUPLICATES(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,ERR,ERROR,*999)
                CALL LIST_DETACH_AND_DESTROY(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,NUMBER_OF_ADJACENT_ELEMENTS, &
                  & ADJACENT_ELEMENTS,ERR,ERROR,*999)
                ! loop over adjacent elements of any element in the current domain 
                DO no_adjacent_element=1,NUMBER_OF_ADJACENT_ELEMENTS
                  adjacent_element=ADJACENT_ELEMENTS(no_adjacent_element)

                  ! compute ghost element number by incrementing the local element counter of this domain 
                  LOCAL_ELEMENT_NUMBERS(domain_idx)=LOCAL_ELEMENT_NUMBERS(domain_idx)+1
                  ! for the adjacent element increase the number of domains where it resides (since it is now a ghost at domain domain_idx) 
                  ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%NUMBER_OF_DOMAINS= &
                    & ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%NUMBER_OF_DOMAINS+1
                  ! set the new element number to the new counter value 
                  ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%LOCAL_NUMBER( &
                    & ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%NUMBER_OF_DOMAINS)=LOCAL_ELEMENT_NUMBERS(domain_idx)
                  ! set the domain to be domain_idx 
                  ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%DOMAIN_NUMBER( &
                    & ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%NUMBER_OF_DOMAINS)=domain_idx
                  ! set the type to be "ghost" 
                  ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%LOCAL_TYPE( &
                    & ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%NUMBER_OF_DOMAINS)= &
                    & DOMAIN_LOCAL_GHOST
                ENDDO !no_adjacent_element
                IF(ALLOCATED(ADJACENT_ELEMENTS)) DEALLOCATE(ADJACENT_ELEMENTS)
              ENDDO !domain_idx
              
              DEALLOCATE(ADJACENT_ELEMENTS_LIST)
              DEALLOCATE(LOCAL_ELEMENT_NUMBERS)
              
              !Calculate element local to global maps from global to local map
              CALL DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(ELEMENTS_MAPPING,ERR,ERROR,*999)
                            
            ELSE
              CALL FlagError("Domain mesh is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Domain decomposition is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Domain mappings elements is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Domain mappings is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain is not associated.",ERR,ERROR,*998)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Element mappings :",ERR,ERROR,*999)
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Global to local map :",ERR,ERROR,*999)
      DO ne=1,MESH%NUMBER_OF_ELEMENTS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global element = ",ne,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of domains  = ", &
          & ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)% &
          & NUMBER_OF_DOMAINS,8,8,ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_NUMBER, &
          & '("      Local number :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)% &
          & NUMBER_OF_DOMAINS,8,8,ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%DOMAIN_NUMBER, &
          & '("      Domain number:",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)% &
          & NUMBER_OF_DOMAINS,8,8,ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_TYPE, &
          & '("      Local type   :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
      ENDDO !ne
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Local to global map :",ERR,ERROR,*999)
      DO ne=1,ELEMENTS_MAPPING%TOTAL_NUMBER_OF_LOCAL
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Local element = ",ne,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Global element = ", &
          & ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(ne),ERR,ERROR,*999)
      ENDDO !ne
      IF(DIAGNOSTICS2) THEN
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Internal elements :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of internal elements = ", &
          & ELEMENTS_MAPPING%NUMBER_OF_INTERNAL,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%NUMBER_OF_INTERNAL,8,8, &
          & ELEMENTS_MAPPING%DOMAIN_LIST(ELEMENTS_MAPPING%INTERNAL_START:ELEMENTS_MAPPING%INTERNAL_FINISH), &
          & '("    Internal elements:",8(X,I7))','(22X,8(X,I7))',ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary elements :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of boundary elements = ", &
          & ELEMENTS_MAPPING%NUMBER_OF_BOUNDARY,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%NUMBER_OF_BOUNDARY,8,8, &
          & ELEMENTS_MAPPING%DOMAIN_LIST(ELEMENTS_MAPPING%BOUNDARY_START:ELEMENTS_MAPPING%BOUNDARY_FINISH), &
          & '("    Boundary elements:",8(X,I7))','(22X,8(X,I7))',ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Ghost elements :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of ghost elements = ", &
          & ELEMENTs_MAPPING%NUMBER_OF_GHOST,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%NUMBER_OF_GHOST,8,8, &
          & ELEMENTS_MAPPING%DOMAIN_LIST(ELEMENTS_MAPPING%GHOST_START:ELEMENTS_MAPPING%GHOST_FINISH), &
          & '("    Ghost elements   :",8(X,I7))','(22X,8(X,I7))',ERR,ERROR,*999)
      ENDIF
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Adjacent domains :",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of adjacent domains = ", &
        & ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS,ERR,ERROR,*999)

      DO domain_idx=1,ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Adjacent domain idx : ",domain_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Domain number = ", &
          & ELEMENTS_MAPPING%ADJACENT_DOMAINS(domain_idx)%DOMAIN_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of send ghosts    = ", &
          & ELEMENTS_MAPPING%ADJACENT_DOMAINS(domain_idx)%NUMBER_OF_SEND_GHOSTS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%ADJACENT_DOMAINS(domain_idx)% &
          & NUMBER_OF_SEND_GHOSTS,6,6,ELEMENTS_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_SEND_INDICES, &
          & '("      Local send ghost indicies       :",6(X,I7))','(39X,6(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of recieve ghosts = ", &
          & ELEMENTS_MAPPING%ADJACENT_DOMAINS(domain_idx)%NUMBER_OF_RECEIVE_GHOSTS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%ADJACENT_DOMAINS(domain_idx)% &
          & NUMBER_OF_RECEIVE_GHOSTS,6,6,ELEMENTS_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_RECEIVE_INDICES, &
          & '("      Local receive ghost indicies    :",6(X,I7))','(39X,6(X,I7))',ERR,ERROR,*999)              
      ENDDO !domain_idx
    ENDIF
    
    EXITS("DOMAIN_MAPPINGS_ELEMENTS_CALCULATE")
    RETURN
999 IF(ALLOCATED(DOMAINS)) DEALLOCATE(DOMAINS)
    IF(ALLOCATED(ADJACENT_ELEMENTS)) DEALLOCATE(ADJACENT_ELEMENTS)    
    IF(ASSOCIATED(DOMAIN%MAPPINGS%ELEMENTS)) CALL DOMAIN_MAPPINGS_ELEMENTS_FINALISE(DOMAIN%MAPPINGS,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("DOMAIN_MAPPINGS_ELEMENTS_CALCULATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Finalises the mappings in the given domain. \todo pass in the domain mappings
  SUBROUTINE DOMAIN_MAPPINGS_FINALISE(DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<A pointer to the domain to finalise the mappings for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_MAPPINGS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      CALL DOMAIN_MAPPINGS_ELEMENTS_FINALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)
      CALL DOMAIN_MAPPINGS_NODES_FINALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)
      CALL DOMAIN_MAPPINGS_DOFS_FINALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)
      DEALLOCATE(DOMAIN%MAPPINGS)
    ELSE
      CALL FlagError("Domain is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_MAPPINGS_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_FINALISE

  !
  !================================================================================================================================
  !

  !>Finalises the element mapping in the given domain mapping. \todo pass in the domain mappings elements
  SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_FINALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS !<A pointer to the domain mappings to finalise the elements for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("DOMAIN_MAPPINGS_ELEMENTS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%ELEMENTS)) THEN
        CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(DOMAIN_MAPPINGS%ELEMENTS,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_MAPPINGS_ELEMENTS_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_ELEMENTS_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the element mapping in the given domain mapping.
  SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_INITIALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS !<A pointer to the domain mappings to initialise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("DOMAIN_MAPPINGS_ELEMENTS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%ELEMENTS)) THEN
        CALL FlagError("Domain elements mappings are already associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(DOMAIN_MAPPINGS%ELEMENTS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate domain mappings elements.",ERR,ERROR,*999)
        CALL DOMAIN_MAPPINGS_MAPPING_INITIALISE(DOMAIN_MAPPINGS%ELEMENTS,DOMAIN_MAPPINGS%DOMAIN%DECOMPOSITION%workGroup, &
          & ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_MAPPINGS_ELEMENTS_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_ELEMENTS_INITIALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the mappings for a domain decomposition. \todo finalise on error.
  SUBROUTINE DOMAIN_MAPPINGS_COMPARE(DOMAIN_MAPPING1,DOMAIN_MAPPING2,NAME,ERR,ERROR,*) 

    !Argument variables

    TYPE(DOMAIN_MAPPING_TYPE), INTENT(IN) :: DOMAIN_MAPPING1 !<A pointer to the first domain mapping
    TYPE(DOMAIN_MAPPING_TYPE), INTENT(IN) :: DOMAIN_MAPPING2 !<A pointer to the second domain mapping to compare to
    CHARACTER(*), INTENT(IN) :: NAME !< A name for the domain mappings to be printed when the mappings are equal 

    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    INTEGER(INTG) :: I,J
    LOGICAL :: MAPPINGS_EQUAL
    MAPPINGS_EQUAL = .TRUE.

    ENTERS("DOMAIN_MAPPINGS_COMPARE",ERR,ERROR,*999)

    IF (DOMAIN_MAPPING1%NUMBER_OF_LOCAL /= DOMAIN_MAPPING2%NUMBER_OF_LOCAL) THEN
      CALL FlagWarning(TRIM(NAME)//" Domain mapping NUMBER_OF_LOCAL is different: "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING1%NUMBER_OF_LOCAL,"*",ERR,ERROR))//" vs. "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%NUMBER_OF_LOCAL,"*",ERR,ERROR)),ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.
    ENDIF
    IF (DOMAIN_MAPPING1%TOTAL_NUMBER_OF_LOCAL /= DOMAIN_MAPPING2%TOTAL_NUMBER_OF_LOCAL) THEN
      CALL FlagWarning(TRIM(NAME)//" Domain mapping TOTAL_NUMBER_OF_LOCAL is different: "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING1%TOTAL_NUMBER_OF_LOCAL,"*",ERR,ERROR))//" vs. "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%TOTAL_NUMBER_OF_LOCAL,"*",ERR,ERROR)),ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.
    ENDIF
    IF (LBOUND(DOMAIN_MAPPING1%NUMBER_OF_DOMAIN_LOCAL,1) /= LBOUND(DOMAIN_MAPPING2%NUMBER_OF_DOMAIN_LOCAL,1) .OR. &
      & UBOUND(DOMAIN_MAPPING1%NUMBER_OF_DOMAIN_LOCAL,1) /= UBOUND(DOMAIN_MAPPING2%NUMBER_OF_DOMAIN_LOCAL,1)) THEN          
      CALL FlagWarning(TRIM(NAME)//" Domain mapping NUMBER_OF_DOMAIN_LOCAL has different bounds: ["//&
              & TRIM(NUMBER_TO_VSTRING(LBOUND(DOMAIN_MAPPING1%NUMBER_OF_DOMAIN_LOCAL,1),"*",ERR,ERROR))//","//&
              & TRIM(NUMBER_TO_VSTRING(UBOUND(DOMAIN_MAPPING1%NUMBER_OF_DOMAIN_LOCAL,1),"*",ERR,ERROR))//"] vs. ["//&
              & TRIM(NUMBER_TO_VSTRING(LBOUND(DOMAIN_MAPPING2%NUMBER_OF_DOMAIN_LOCAL,1),"*",ERR,ERROR))//","//&
              & TRIM(NUMBER_TO_VSTRING(UBOUND(DOMAIN_MAPPING2%NUMBER_OF_DOMAIN_LOCAL,1),"*",ERR,ERROR))//"].",ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.
    ELSE
      DO I=LBOUND(DOMAIN_MAPPING1%NUMBER_OF_DOMAIN_LOCAL,1),UBOUND(DOMAIN_MAPPING1%NUMBER_OF_DOMAIN_LOCAL,1)
        IF (DOMAIN_MAPPING1%NUMBER_OF_DOMAIN_LOCAL(I) /= DOMAIN_MAPPING2%NUMBER_OF_DOMAIN_LOCAL(I)) THEN
          CALL FlagWarning(TRIM(NAME)//" Domain mapping NUMBER_OF_DOMAIN_LOCAL("//&
                & TRIM(NUMBER_TO_VSTRING(I,"*",ERR,ERROR))//") is different: "//&
                  & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING1%NUMBER_OF_DOMAIN_LOCAL(I),"*",ERR,ERROR))//" vs. "//&
                  & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%NUMBER_OF_DOMAIN_LOCAL(I),"*",ERR,ERROR)),ERR,ERROR,*999)
        MAPPINGS_EQUAL = .FALSE.
        ENDIF
      ENDDO
    ENDIF
    IF (LBOUND(DOMAIN_MAPPING1%NUMBER_OF_DOMAIN_GHOST,1) /= LBOUND(DOMAIN_MAPPING2%NUMBER_OF_DOMAIN_GHOST,1) .OR. &
      & UBOUND(DOMAIN_MAPPING1%NUMBER_OF_DOMAIN_GHOST,1) /= UBOUND(DOMAIN_MAPPING2%NUMBER_OF_DOMAIN_GHOST,1)) THEN          
      CALL FlagWarning(TRIM(NAME)//" Domain mapping NUMBER_OF_DOMAIN_GHOST has different bounds: ["//&
              & TRIM(NUMBER_TO_VSTRING(LBOUND(DOMAIN_MAPPING1%NUMBER_OF_DOMAIN_GHOST,1),"*",ERR,ERROR))//","//&
              & TRIM(NUMBER_TO_VSTRING(UBOUND(DOMAIN_MAPPING1%NUMBER_OF_DOMAIN_GHOST,1),"*",ERR,ERROR))//"] vs. ["//&
              & TRIM(NUMBER_TO_VSTRING(LBOUND(DOMAIN_MAPPING2%NUMBER_OF_DOMAIN_GHOST,1),"*",ERR,ERROR))//","//&
              & TRIM(NUMBER_TO_VSTRING(UBOUND(DOMAIN_MAPPING2%NUMBER_OF_DOMAIN_GHOST,1),"*",ERR,ERROR))//"].",ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.
    ELSE
      DO I=LBOUND(DOMAIN_MAPPING1%NUMBER_OF_DOMAIN_GHOST,1),UBOUND(DOMAIN_MAPPING1%NUMBER_OF_DOMAIN_GHOST,1)
        IF (DOMAIN_MAPPING1%NUMBER_OF_DOMAIN_GHOST(I) /= DOMAIN_MAPPING2%NUMBER_OF_DOMAIN_GHOST(I)) THEN
          CALL FlagWarning(TRIM(NAME)//" Domain mapping NUMBER_OF_DOMAIN_GHOST("//&
                & TRIM(NUMBER_TO_VSTRING(I,"*",ERR,ERROR))//") is different: "//&
                  & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING1%NUMBER_OF_DOMAIN_GHOST(I),"*",ERR,ERROR))//" vs. "//&
                  & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%NUMBER_OF_DOMAIN_GHOST(I),"*",ERR,ERROR)),ERR,ERROR,*999)
        MAPPINGS_EQUAL = .FALSE.
        ENDIF
      ENDDO
    ENDIF
    IF (DOMAIN_MAPPING1%NUMBER_OF_GLOBAL /= DOMAIN_MAPPING2%NUMBER_OF_GLOBAL) THEN
      CALL FlagWarning(TRIM(NAME)//" Domain mapping NUMBER_OF_GLOBAL is different: "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING1%NUMBER_OF_GLOBAL,"*",ERR,ERROR))//" vs. "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%NUMBER_OF_GLOBAL,"*",ERR,ERROR)),ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.
    ENDIF
    IF (DOMAIN_MAPPING1%NUMBER_OF_DOMAINS /= DOMAIN_MAPPING2%NUMBER_OF_DOMAINS) THEN
      CALL FlagWarning(TRIM(NAME)//" Domain mapping NUMBER_OF_DOMAINS is different: "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING1%NUMBER_OF_DOMAINS,"*",ERR,ERROR))//" vs. "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%NUMBER_OF_DOMAINS,"*",ERR,ERROR)),ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.
    ENDIF
    IF (DOMAIN_MAPPING1%NUMBER_OF_INTERNAL /= DOMAIN_MAPPING2%NUMBER_OF_INTERNAL) THEN
      CALL FlagWarning(TRIM(NAME)//" Domain mapping NUMBER_OF_INTERNAL is different: "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING1%NUMBER_OF_INTERNAL,"*",ERR,ERROR))//" vs. "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%NUMBER_OF_INTERNAL,"*",ERR,ERROR)),ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.
    ENDIF
    IF (DOMAIN_MAPPING1%NUMBER_OF_BOUNDARY /= DOMAIN_MAPPING2%NUMBER_OF_BOUNDARY) THEN
      CALL FlagWarning(TRIM(NAME)//" Domain mapping NUMBER_OF_BOUNDARY is different: "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING1%NUMBER_OF_BOUNDARY,"*",ERR,ERROR))//" vs. "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%NUMBER_OF_BOUNDARY,"*",ERR,ERROR)),ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.
    ENDIF
    IF (DOMAIN_MAPPING1%NUMBER_OF_GHOST /= DOMAIN_MAPPING2%NUMBER_OF_GHOST) THEN
      CALL FlagWarning(TRIM(NAME)//" Domain mapping NUMBER_OF_GHOST is different: "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING1%NUMBER_OF_GHOST,"*",ERR,ERROR))//" vs. "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%NUMBER_OF_GHOST,"*",ERR,ERROR)),ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.
    ENDIF
    IF (DOMAIN_MAPPING1%INTERNAL_START /= DOMAIN_MAPPING2%INTERNAL_START) THEN
      CALL FlagWarning(TRIM(NAME)//" Domain mapping INTERNAL_START is different: "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING1%INTERNAL_START,"*",ERR,ERROR))//" vs. "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%INTERNAL_START,"*",ERR,ERROR)),ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.
    ENDIF
    IF (DOMAIN_MAPPING1%INTERNAL_FINISH /= DOMAIN_MAPPING2%INTERNAL_FINISH) THEN
      CALL FlagWarning(TRIM(NAME)//" Domain mapping INTERNAL_FINISH is different: "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING1%INTERNAL_FINISH,"*",ERR,ERROR))//" vs. "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%INTERNAL_FINISH,"*",ERR,ERROR)),ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.
    ENDIF
    IF (DOMAIN_MAPPING1%BOUNDARY_START /= DOMAIN_MAPPING2%BOUNDARY_START) THEN
      CALL FlagWarning(TRIM(NAME)//" Domain mapping BOUNDARY_START is different: "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING1%BOUNDARY_START,"*",ERR,ERROR))//" vs. "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%BOUNDARY_START,"*",ERR,ERROR)),ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.
    ENDIF
    IF (DOMAIN_MAPPING1%BOUNDARY_FINISH /= DOMAIN_MAPPING2%BOUNDARY_FINISH) THEN
      CALL FlagWarning(TRIM(NAME)//" Domain mapping BOUNDARY_FINISH is different: "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING1%BOUNDARY_FINISH,"*",ERR,ERROR))//" vs. "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%BOUNDARY_FINISH,"*",ERR,ERROR)),ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.
    ENDIF
    IF (DOMAIN_MAPPING1%GHOST_START /= DOMAIN_MAPPING2%GHOST_START) THEN
      CALL FlagWarning(TRIM(NAME)//" Domain mapping GHOST_START is different: "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING1%GHOST_START,"*",ERR,ERROR))//" vs. "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%GHOST_START,"*",ERR,ERROR)),ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.
    ENDIF
    IF (DOMAIN_MAPPING1%GHOST_FINISH /= DOMAIN_MAPPING2%GHOST_FINISH) THEN
      CALL FlagWarning(TRIM(NAME)//" Domain mapping GHOST_FINISH is different: "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING1%GHOST_FINISH,"*",ERR,ERROR))//" vs. "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%GHOST_FINISH,"*",ERR,ERROR)),ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.
    ENDIF
    IF (LBOUND(DOMAIN_MAPPING1%DOMAIN_LIST,1) /= LBOUND(DOMAIN_MAPPING2%DOMAIN_LIST,1) .OR. &
      & UBOUND(DOMAIN_MAPPING1%DOMAIN_LIST,1) /= UBOUND(DOMAIN_MAPPING2%DOMAIN_LIST,1)) THEN          
      CALL FlagWarning(TRIM(NAME)//" Domain mapping DOMAIN_LIST has different bounds: ["//&
              & TRIM(NUMBER_TO_VSTRING(LBOUND(DOMAIN_MAPPING1%DOMAIN_LIST,1),"*",ERR,ERROR))//","//&
              & TRIM(NUMBER_TO_VSTRING(UBOUND(DOMAIN_MAPPING1%DOMAIN_LIST,1),"*",ERR,ERROR))//"] vs. ["//&
              & TRIM(NUMBER_TO_VSTRING(LBOUND(DOMAIN_MAPPING2%DOMAIN_LIST,1),"*",ERR,ERROR))//","//&
              & TRIM(NUMBER_TO_VSTRING(UBOUND(DOMAIN_MAPPING2%DOMAIN_LIST,1),"*",ERR,ERROR))//"].",ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.
    ELSE
      DO I=LBOUND(DOMAIN_MAPPING1%DOMAIN_LIST,1),UBOUND(DOMAIN_MAPPING1%DOMAIN_LIST,1)
        IF (DOMAIN_MAPPING1%DOMAIN_LIST(I) /= DOMAIN_MAPPING2%DOMAIN_LIST(I)) THEN
          CALL FlagWarning(TRIM(NAME)//" Domain mapping DOMAIN_LIST("//&
                & TRIM(NUMBER_TO_VSTRING(I,"*",ERR,ERROR))//") is different: "//&
                  & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING1%DOMAIN_LIST(I),"*",ERR,ERROR))//" vs. "//&
                  & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%DOMAIN_LIST(I),"*",ERR,ERROR)),ERR,ERROR,*999)
        MAPPINGS_EQUAL = .FALSE.
        ENDIF
      ENDDO
    ENDIF
    IF (LBOUND(DOMAIN_MAPPING1%LOCAL_TO_GLOBAL_MAP,1) /= LBOUND(DOMAIN_MAPPING2%LOCAL_TO_GLOBAL_MAP,1) .OR. &
      & UBOUND(DOMAIN_MAPPING1%LOCAL_TO_GLOBAL_MAP,1) /= UBOUND(DOMAIN_MAPPING2%LOCAL_TO_GLOBAL_MAP,1)) THEN          
      CALL FlagWarning(TRIM(NAME)//" Domain mapping LOCAL_TO_GLOBAL_MAP has different bounds: ["//&
              & TRIM(NUMBER_TO_VSTRING(LBOUND(DOMAIN_MAPPING1%LOCAL_TO_GLOBAL_MAP,1),"*",ERR,ERROR))//","//&
              & TRIM(NUMBER_TO_VSTRING(UBOUND(DOMAIN_MAPPING1%LOCAL_TO_GLOBAL_MAP,1),"*",ERR,ERROR))//"] vs. ["//&
              & TRIM(NUMBER_TO_VSTRING(LBOUND(DOMAIN_MAPPING2%LOCAL_TO_GLOBAL_MAP,1),"*",ERR,ERROR))//","//&
              & TRIM(NUMBER_TO_VSTRING(UBOUND(DOMAIN_MAPPING2%LOCAL_TO_GLOBAL_MAP,1),"*",ERR,ERROR))//"].",ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.
    ELSE
      DO I=LBOUND(DOMAIN_MAPPING1%LOCAL_TO_GLOBAL_MAP,1),UBOUND(DOMAIN_MAPPING1%LOCAL_TO_GLOBAL_MAP,1)
        IF (DOMAIN_MAPPING1%LOCAL_TO_GLOBAL_MAP(I) /= DOMAIN_MAPPING2%LOCAL_TO_GLOBAL_MAP(I)) THEN
          CALL FlagWarning(TRIM(NAME)//" Domain mapping LOCAL_TO_GLOBAL_MAP("//&
                & TRIM(NUMBER_TO_VSTRING(I,"*",ERR,ERROR))//") is different: "//&
                  & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING1%LOCAL_TO_GLOBAL_MAP(I),"*",ERR,ERROR))//" vs. "//&
                  & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%LOCAL_TO_GLOBAL_MAP(I),"*",ERR,ERROR)),ERR,ERROR,*999)
        MAPPINGS_EQUAL = .FALSE.
        ENDIF
      ENDDO
    ENDIF
    IF (DOMAIN_MAPPING1%NUMBER_OF_ADJACENT_DOMAINS /= DOMAIN_MAPPING2%NUMBER_OF_ADJACENT_DOMAINS) THEN
      CALL FlagWarning(TRIM(NAME)//" Domain mapping NUMBER_OF_ADJACENT_DOMAINS is different: "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING1%NUMBER_OF_ADJACENT_DOMAINS,"*",ERR,ERROR))//" vs. "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%NUMBER_OF_ADJACENT_DOMAINS,"*",ERR,ERROR)),ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.
    ENDIF
    IF (LBOUND(DOMAIN_MAPPING1%ADJACENT_DOMAINS,1) /= LBOUND(DOMAIN_MAPPING2%ADJACENT_DOMAINS,1) .OR. &
      & UBOUND(DOMAIN_MAPPING1%ADJACENT_DOMAINS,1) /= UBOUND(DOMAIN_MAPPING2%ADJACENT_DOMAINS,1)) THEN          
      CALL FlagWarning(TRIM(NAME)//" Domain mapping ADJACENT_DOMAINS has different bounds: ["//&
              & TRIM(NUMBER_TO_VSTRING(LBOUND(DOMAIN_MAPPING1%ADJACENT_DOMAINS,1),"*",ERR,ERROR))//","//&
              & TRIM(NUMBER_TO_VSTRING(UBOUND(DOMAIN_MAPPING1%ADJACENT_DOMAINS,1),"*",ERR,ERROR))//"] vs. ["//&
              & TRIM(NUMBER_TO_VSTRING(LBOUND(DOMAIN_MAPPING2%ADJACENT_DOMAINS,1),"*",ERR,ERROR))//","//&
              & TRIM(NUMBER_TO_VSTRING(UBOUND(DOMAIN_MAPPING2%ADJACENT_DOMAINS,1),"*",ERR,ERROR))//"].",ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.
    ELSE
      DO I=LBOUND(DOMAIN_MAPPING1%ADJACENT_DOMAINS,1),UBOUND(DOMAIN_MAPPING1%ADJACENT_DOMAINS,1)
        IF (DOMAIN_MAPPING1%ADJACENT_DOMAINS(I)%DOMAIN_NUMBER /= DOMAIN_MAPPING2%ADJACENT_DOMAINS(I)%DOMAIN_NUMBER) THEN
          CALL FlagWarning(TRIM(NAME)//" Domain mapping ADJACENT_DOMAINS("//&
                & TRIM(NUMBER_TO_VSTRING(I,"*",ERR,ERROR))//")%DOMAIN_NUMBER is different: "//&
                  & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING1%ADJACENT_DOMAINS(I)%DOMAIN_NUMBER,"*",ERR,ERROR))//" vs. "//&
                  & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%ADJACENT_DOMAINS(I)%DOMAIN_NUMBER,"*",ERR,ERROR)),ERR,ERROR,*999)
        MAPPINGS_EQUAL = .FALSE.
        ENDIF
        IF (DOMAIN_MAPPING1%ADJACENT_DOMAINS(I)%NUMBER_OF_SEND_GHOSTS /= DOMAIN_MAPPING2%ADJACENT_DOMAINS(I)%&
          & NUMBER_OF_SEND_GHOSTS) THEN
          CALL FlagWarning(TRIM(NAME)//" Domain mapping ADJACENT_DOMAINS("//&
                & TRIM(NUMBER_TO_VSTRING(I,"*",ERR,ERROR))//")%NUMBER_OF_SEND_GHOSTS is different: "//&
                  & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING1%ADJACENT_DOMAINS(I)%NUMBER_OF_SEND_GHOSTS,"*",ERR,ERROR))&
                  & //" vs. "//&
                  & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%ADJACENT_DOMAINS(I)%NUMBER_OF_SEND_GHOSTS,"*",ERR,ERROR)),&
                  & ERR,ERROR,*999)
        MAPPINGS_EQUAL = .FALSE.
        ENDIF
        IF (DOMAIN_MAPPING1%ADJACENT_DOMAINS(I)%NUMBER_OF_RECEIVE_GHOSTS /= DOMAIN_MAPPING2%ADJACENT_DOMAINS(I)%&
          & NUMBER_OF_RECEIVE_GHOSTS) THEN
          CALL FlagWarning(TRIM(NAME)//" Domain mapping ADJACENT_DOMAINS("//&
                & TRIM(NUMBER_TO_VSTRING(I,"*",ERR,ERROR))//")%NUMBER_OF_RECEIVE_GHOSTS is different: "//&
                & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING1%ADJACENT_DOMAINS(I)%NUMBER_OF_RECEIVE_GHOSTS,"*",ERR,ERROR))&
                & //" vs. "//&
                & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%ADJACENT_DOMAINS(I)%NUMBER_OF_RECEIVE_GHOSTS,"*",ERR,ERROR)),&
                & ERR,ERROR,*999)
        MAPPINGS_EQUAL = .FALSE.
        ENDIF
        IF (LBOUND(DOMAIN_MAPPING1%ADJACENT_DOMAINS(I)%LOCAL_GHOST_SEND_INDICES,1) /= &
          & LBOUND(DOMAIN_MAPPING2%ADJACENT_DOMAINS(I)%LOCAL_GHOST_SEND_INDICES,1) .OR. &
          & UBOUND(DOMAIN_MAPPING1%ADJACENT_DOMAINS(I)%LOCAL_GHOST_SEND_INDICES,1) /= &
          & UBOUND(DOMAIN_MAPPING2%ADJACENT_DOMAINS(I)%LOCAL_GHOST_SEND_INDICES,1)) THEN          
          CALL FlagWarning(TRIM(NAME)//" Domain mapping ADJACENT_DOMAINS("//&
                & TRIM(NUMBER_TO_VSTRING(I,"*",ERR,ERROR))//")%LOCAL_GHOST_SEND_INDICES has different bounds: ["//&
                  & TRIM(NUMBER_TO_VSTRING(LBOUND(DOMAIN_MAPPING1%ADJACENT_DOMAINS(I)%LOCAL_GHOST_SEND_INDICES,1),&
                  & "*",ERR,ERROR))//","//&
                  & TRIM(NUMBER_TO_VSTRING(UBOUND(DOMAIN_MAPPING1%ADJACENT_DOMAINS(I)%LOCAL_GHOST_SEND_INDICES,1),&
                  & "*",ERR,ERROR))//"] vs. ["//&
                  & TRIM(NUMBER_TO_VSTRING(LBOUND(DOMAIN_MAPPING2%ADJACENT_DOMAINS(I)%LOCAL_GHOST_SEND_INDICES,1),&
                  & "*",ERR,ERROR))//","//&
                  & TRIM(NUMBER_TO_VSTRING(UBOUND(DOMAIN_MAPPING2%ADJACENT_DOMAINS(I)%LOCAL_GHOST_SEND_INDICES,1),&
                  & "*",ERR,ERROR))//"].",ERR,ERROR,*999)
        MAPPINGS_EQUAL = .FALSE.
        ELSE
          DO J=LBOUND(DOMAIN_MAPPING1%ADJACENT_DOMAINS(I)%LOCAL_GHOST_SEND_INDICES,1),&
            & UBOUND(DOMAIN_MAPPING1%ADJACENT_DOMAINS(I)%LOCAL_GHOST_SEND_INDICES,1)
            IF (DOMAIN_MAPPING1%ADJACENT_DOMAINS(I)%LOCAL_GHOST_SEND_INDICES(J) /= &
              & DOMAIN_MAPPING2%ADJACENT_DOMAINS(I)%LOCAL_GHOST_SEND_INDICES(J)) THEN
              CALL FlagWarning(TRIM(NAME)//" Domain mapping ADJACENT_DOMAINS("//&
                    & TRIM(NUMBER_TO_VSTRING(I,"*",ERR,ERROR))//")%LOCAL_GHOST_SEND_INDICES("//&
                    & TRIM(NUMBER_TO_VSTRING(J,"*",ERR,ERROR))//") is different: "//&
                    & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING1%ADJACENT_DOMAINS(I)%LOCAL_GHOST_SEND_INDICES(J),&
                    & "*",ERR,ERROR))//" vs. "//&
                    & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%ADJACENT_DOMAINS(I)%LOCAL_GHOST_SEND_INDICES(J),&
                    & "*",ERR,ERROR)),ERR,ERROR,*999)
              MAPPINGS_EQUAL = .FALSE.
            ENDIF
          ENDDO
        ENDIF
        IF (LBOUND(DOMAIN_MAPPING1%ADJACENT_DOMAINS(I)%LOCAL_GHOST_RECEIVE_INDICES,1) /= &
          & LBOUND(DOMAIN_MAPPING2%ADJACENT_DOMAINS(I)%LOCAL_GHOST_RECEIVE_INDICES,1) .OR. &
          & UBOUND(DOMAIN_MAPPING1%ADJACENT_DOMAINS(I)%LOCAL_GHOST_RECEIVE_INDICES,1) /= &
          & UBOUND(DOMAIN_MAPPING2%ADJACENT_DOMAINS(I)%LOCAL_GHOST_RECEIVE_INDICES,1)) THEN          
          CALL FlagWarning(TRIM(NAME)//" Domain mapping ADJACENT_DOMAINS("//&
                & TRIM(NUMBER_TO_VSTRING(I,"*",ERR,ERROR))//")%LOCAL_GHOST_RECEIVE_INDICES has different bounds: ["//&
                  & TRIM(NUMBER_TO_VSTRING(LBOUND(DOMAIN_MAPPING1%ADJACENT_DOMAINS(I)%LOCAL_GHOST_RECEIVE_INDICES,1),&
                  & "*",ERR,ERROR))//","//&
                  & TRIM(NUMBER_TO_VSTRING(UBOUND(DOMAIN_MAPPING1%ADJACENT_DOMAINS(I)%LOCAL_GHOST_RECEIVE_INDICES,1),&
                  & "*",ERR,ERROR))//"] vs. ["//&
                  & TRIM(NUMBER_TO_VSTRING(LBOUND(DOMAIN_MAPPING2%ADJACENT_DOMAINS(I)%LOCAL_GHOST_RECEIVE_INDICES,1),&
                  & "*",ERR,ERROR))//","//&
                  & TRIM(NUMBER_TO_VSTRING(UBOUND(DOMAIN_MAPPING2%ADJACENT_DOMAINS(I)%LOCAL_GHOST_RECEIVE_INDICES,1),&
                  & "*",ERR,ERROR))//"].",ERR,ERROR,*999)
          MAPPINGS_EQUAL = .FALSE.
        ELSE
          DO J=LBOUND(DOMAIN_MAPPING1%ADJACENT_DOMAINS(I)%LOCAL_GHOST_RECEIVE_INDICES,1),&
            & UBOUND(DOMAIN_MAPPING1%ADJACENT_DOMAINS(I)%LOCAL_GHOST_RECEIVE_INDICES,1)
            IF (DOMAIN_MAPPING1%ADJACENT_DOMAINS(I)%LOCAL_GHOST_RECEIVE_INDICES(J) /= &
              & DOMAIN_MAPPING2%ADJACENT_DOMAINS(I)%LOCAL_GHOST_RECEIVE_INDICES(J)) THEN
              CALL FlagWarning(TRIM(NAME)//" Domain mapping ADJACENT_DOMAINS("//&
                    & TRIM(NUMBER_TO_VSTRING(I,"*",ERR,ERROR))//")%LOCAL_GHOST_RECEIVE_INDICES("//&
                    & TRIM(NUMBER_TO_VSTRING(J,"*",ERR,ERROR))//") is different: "//&
                    & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING1%ADJACENT_DOMAINS(I)%LOCAL_GHOST_RECEIVE_INDICES(J),&
                    & "*",ERR,ERROR))//" vs. "//&
                    & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%ADJACENT_DOMAINS(I)%LOCAL_GHOST_RECEIVE_INDICES(J),&
                    & "*",ERR,ERROR)),ERR,ERROR,*999)
              MAPPINGS_EQUAL = .FALSE.
            ENDIF
          ENDDO
        ENDIF
      ENDDO
    ENDIF
    
    IF (MAPPINGS_EQUAL) CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,TRIM(NAME)//" Mappings are equal",ERR,ERROR,*999)
        
    
    EXITS("DOMAIN_MAPPINGS_COMPARE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_COMPARE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_COMPARE

  !
  !================================================================================================================================
  !
  !>Print out one mappings characteristics of a domain decomposition. \todo finalise on error.
  SUBROUTINE DOMAIN_MAPPINGS_PRINT_NUMB(DOMAIN_MAPPING2,NAME,ERR,ERROR,*) 

    !Argument variables

   ! TYPE(DOMAIN_MAPPING_TYPE), INTENT(IN) :: DOMAIN_MAPPING1 !<A pointer to the first domain mapping (new)
    TYPE(DOMAIN_MAPPING_TYPE), INTENT(IN) :: DOMAIN_MAPPING2 !<A pointer to the second domain mapping to compare to (old)
    CHARACTER(*), INTENT(IN) :: NAME !< A name for the domain mappings to be printed when the mappings are equal 

    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    INTEGER(INTG) :: I,J
    LOGICAL :: MAPPINGS_EQUAL
    MAPPINGS_EQUAL = .TRUE.

    ENTERS("DOMAIN_MAPPINGS_PRINT_NUMB",ERR,ERROR,*999)

      CALL FlagWarning(TRIM(NAME)//" Domain mapping NUMBER_OF_LOCAL is : "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%NUMBER_OF_LOCAL,"*",ERR,ERROR)),ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.

      CALL FlagWarning(TRIM(NAME)//" Domain mapping TOTAL_NUMBER_OF_LOCAL is : "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%TOTAL_NUMBER_OF_LOCAL,"*",ERR,ERROR)),ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.

      CALL FlagWarning(TRIM(NAME)//" Domain mapping NUMBER_OF_DOMAIN_LOCAL has  bounds: ["//&
              & TRIM(NUMBER_TO_VSTRING(LBOUND(DOMAIN_MAPPING2%NUMBER_OF_DOMAIN_LOCAL,1),"*",ERR,ERROR))//","//&
              & TRIM(NUMBER_TO_VSTRING(UBOUND(DOMAIN_MAPPING2%NUMBER_OF_DOMAIN_LOCAL,1),"*",ERR,ERROR))//"].",ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.

      DO I=LBOUND(DOMAIN_MAPPING2%NUMBER_OF_DOMAIN_LOCAL,1),UBOUND(DOMAIN_MAPPING2%NUMBER_OF_DOMAIN_LOCAL,1)
          CALL FlagWarning(TRIM(NAME)//" Domain mapping NUMBER_OF_DOMAIN_LOCAL("//&
                & TRIM(NUMBER_TO_VSTRING(I,"*",ERR,ERROR))//") is : "//&
                  & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%NUMBER_OF_DOMAIN_LOCAL(I),"*",ERR,ERROR)),ERR,ERROR,*999)
        MAPPINGS_EQUAL = .FALSE.
      ENDDO

      CALL FlagWarning(TRIM(NAME)//" Domain mapping NUMBER_OF_DOMAIN_GHOST has  bounds: ["//&
              & TRIM(NUMBER_TO_VSTRING(LBOUND(DOMAIN_MAPPING2%NUMBER_OF_DOMAIN_GHOST,1),"*",ERR,ERROR))//","//&
              & TRIM(NUMBER_TO_VSTRING(UBOUND(DOMAIN_MAPPING2%NUMBER_OF_DOMAIN_GHOST,1),"*",ERR,ERROR))//"].",ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.

      DO I=LBOUND(DOMAIN_MAPPING2%NUMBER_OF_DOMAIN_GHOST,1),UBOUND(DOMAIN_MAPPING2%NUMBER_OF_DOMAIN_GHOST,1)
          CALL FlagWarning(TRIM(NAME)//" Domain mapping NUMBER_OF_DOMAIN_GHOST("//&
                & TRIM(NUMBER_TO_VSTRING(I,"*",ERR,ERROR))//") is : "//&
                  & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%NUMBER_OF_DOMAIN_GHOST(I),"*",ERR,ERROR)),ERR,ERROR,*999)
        MAPPINGS_EQUAL = .FALSE.
      ENDDO

      CALL FlagWarning(TRIM(NAME)//" Domain mapping NUMBER_OF_GLOBAL is : "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%NUMBER_OF_GLOBAL,"*",ERR,ERROR)),ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.

      CALL FlagWarning(TRIM(NAME)//" Domain mapping NUMBER_OF_DOMAINS is : "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%NUMBER_OF_DOMAINS,"*",ERR,ERROR)),ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.

      CALL FlagWarning(TRIM(NAME)//" Domain mapping NUMBER_OF_INTERNAL is : "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%NUMBER_OF_INTERNAL,"*",ERR,ERROR)),ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.

      CALL FlagWarning(TRIM(NAME)//" Domain mapping NUMBER_OF_BOUNDARY is : "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%NUMBER_OF_BOUNDARY,"*",ERR,ERROR)),ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.

      CALL FlagWarning(TRIM(NAME)//" Domain mapping NUMBER_OF_GHOST is : "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%NUMBER_OF_GHOST,"*",ERR,ERROR)),ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.

      CALL FlagWarning(TRIM(NAME)//" Domain mapping INTERNAL_START is : "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%INTERNAL_START,"*",ERR,ERROR)),ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.

      CALL FlagWarning(TRIM(NAME)//" Domain mapping INTERNAL_FINISH is : "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%INTERNAL_FINISH,"*",ERR,ERROR)),ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.

      CALL FlagWarning(TRIM(NAME)//" Domain mapping BOUNDARY_START is : "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%BOUNDARY_START,"*",ERR,ERROR)),ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.

      CALL FlagWarning(TRIM(NAME)//" Domain mapping BOUNDARY_FINISH is : "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%BOUNDARY_FINISH,"*",ERR,ERROR)),ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.

      CALL FlagWarning(TRIM(NAME)//" Domain mapping GHOST_START is : "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%GHOST_START,"*",ERR,ERROR)),ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.

      CALL FlagWarning(TRIM(NAME)//" Domain mapping GHOST_FINISH is : "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%GHOST_FINISH,"*",ERR,ERROR)),ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.

      CALL FlagWarning(TRIM(NAME)//" Domain mapping DOMAIN_LIST has  bounds: ["//&
              & TRIM(NUMBER_TO_VSTRING(LBOUND(DOMAIN_MAPPING2%DOMAIN_LIST,1),"*",ERR,ERROR))//","//&
              & TRIM(NUMBER_TO_VSTRING(UBOUND(DOMAIN_MAPPING2%DOMAIN_LIST,1),"*",ERR,ERROR))//"].",ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.

      DO I=LBOUND(DOMAIN_MAPPING2%DOMAIN_LIST,1),UBOUND(DOMAIN_MAPPING2%DOMAIN_LIST,1)
          CALL FlagWarning(TRIM(NAME)//" Domain mapping DOMAIN_LIST("//&
                & TRIM(NUMBER_TO_VSTRING(I,"*",ERR,ERROR))//") is : "//&
                  & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%DOMAIN_LIST(I),"*",ERR,ERROR)),ERR,ERROR,*999)
        MAPPINGS_EQUAL = .FALSE.
      ENDDO

      CALL FlagWarning(TRIM(NAME)//" Domain mapping LOCAL_TO_GLOBAL_MAP has  bounds: ["//&
              & TRIM(NUMBER_TO_VSTRING(LBOUND(DOMAIN_MAPPING2%LOCAL_TO_GLOBAL_MAP,1),"*",ERR,ERROR))//","//&
              & TRIM(NUMBER_TO_VSTRING(UBOUND(DOMAIN_MAPPING2%LOCAL_TO_GLOBAL_MAP,1),"*",ERR,ERROR))//"].",ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.

      DO I=LBOUND(DOMAIN_MAPPING2%LOCAL_TO_GLOBAL_MAP,1),UBOUND(DOMAIN_MAPPING2%LOCAL_TO_GLOBAL_MAP,1)
          CALL FlagWarning(TRIM(NAME)//" Domain mapping LOCAL_TO_GLOBAL_MAP("//&
                & TRIM(NUMBER_TO_VSTRING(I,"*",ERR,ERROR))//") is : "//&
                  & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%LOCAL_TO_GLOBAL_MAP(I),"*",ERR,ERROR)),ERR,ERROR,*999)
        MAPPINGS_EQUAL = .FALSE.
      ENDDO

      CALL FlagWarning(TRIM(NAME)//" Domain mapping NUMBER_OF_ADJACENT_DOMAINS is : "//&
              & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%NUMBER_OF_ADJACENT_DOMAINS,"*",ERR,ERROR)),ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.

      CALL FlagWarning(TRIM(NAME)//" Domain mapping ADJACENT_DOMAINS has  bounds: ["//&
              & TRIM(NUMBER_TO_VSTRING(LBOUND(DOMAIN_MAPPING2%ADJACENT_DOMAINS,1),"*",ERR,ERROR))//","//&
              & TRIM(NUMBER_TO_VSTRING(UBOUND(DOMAIN_MAPPING2%ADJACENT_DOMAINS,1),"*",ERR,ERROR))//"].",ERR,ERROR,*999)
      MAPPINGS_EQUAL = .FALSE.

      DO I=LBOUND(DOMAIN_MAPPING2%ADJACENT_DOMAINS,1),UBOUND(DOMAIN_MAPPING2%ADJACENT_DOMAINS,1)
          CALL FlagWarning(TRIM(NAME)//" Domain mapping ADJACENT_DOMAINS("//&
                & TRIM(NUMBER_TO_VSTRING(I,"*",ERR,ERROR))//")%DOMAIN_NUMBER is : "//&
                  & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%ADJACENT_DOMAINS(I)%DOMAIN_NUMBER,"*",ERR,ERROR)),ERR,ERROR,*999)
        MAPPINGS_EQUAL = .FALSE.

          CALL FlagWarning(TRIM(NAME)//" Domain mapping ADJACENT_DOMAINS("//&
                & TRIM(NUMBER_TO_VSTRING(I,"*",ERR,ERROR))//")%NUMBER_OF_SEND_GHOSTS is : "//&
                  & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%ADJACENT_DOMAINS(I)%NUMBER_OF_SEND_GHOSTS,"*",ERR,ERROR)),&
                  & ERR,ERROR,*999)
        MAPPINGS_EQUAL = .FALSE.

          CALL FlagWarning(TRIM(NAME)//" Domain mapping ADJACENT_DOMAINS("//&
                & TRIM(NUMBER_TO_VSTRING(I,"*",ERR,ERROR))//")%NUMBER_OF_RECEIVE_GHOSTS is : "//&
                & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%ADJACENT_DOMAINS(I)%NUMBER_OF_RECEIVE_GHOSTS,"*",ERR,ERROR)),&
                & ERR,ERROR,*999)
        MAPPINGS_EQUAL = .FALSE.
          CALL FlagWarning(TRIM(NAME)//" Domain mapping ADJACENT_DOMAINS("//&
                & TRIM(NUMBER_TO_VSTRING(I,"*",ERR,ERROR))//")%LOCAL_GHOST_SEND_INDICES has  bounds: ["//&
                  & TRIM(NUMBER_TO_VSTRING(LBOUND(DOMAIN_MAPPING2%ADJACENT_DOMAINS(I)%LOCAL_GHOST_SEND_INDICES,1),&
                  & "*",ERR,ERROR))//","//&
                  & TRIM(NUMBER_TO_VSTRING(UBOUND(DOMAIN_MAPPING2%ADJACENT_DOMAINS(I)%LOCAL_GHOST_SEND_INDICES,1),&
                  & "*",ERR,ERROR))//"].",ERR,ERROR,*999)
        MAPPINGS_EQUAL = .FALSE.

          DO J=LBOUND(DOMAIN_MAPPING2%ADJACENT_DOMAINS(I)%LOCAL_GHOST_SEND_INDICES,1),&
            & UBOUND(DOMAIN_MAPPING2%ADJACENT_DOMAINS(I)%LOCAL_GHOST_SEND_INDICES,1)
              CALL FlagWarning(TRIM(NAME)//" Domain mapping ADJACENT_DOMAINS("//&
                    & TRIM(NUMBER_TO_VSTRING(I,"*",ERR,ERROR))//")%LOCAL_GHOST_SEND_INDICES("//&
                    & TRIM(NUMBER_TO_VSTRING(J,"*",ERR,ERROR))//") is : "//&
                    & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%ADJACENT_DOMAINS(I)%LOCAL_GHOST_SEND_INDICES(J),&
                    & "*",ERR,ERROR)),ERR,ERROR,*999)
              MAPPINGS_EQUAL = .FALSE.
          ENDDO

          CALL FlagWarning(TRIM(NAME)//" Domain mapping ADJACENT_DOMAINS("//&
                & TRIM(NUMBER_TO_VSTRING(I,"*",ERR,ERROR))//")%LOCAL_GHOST_RECEIVE_INDICES has  bounds: ["//&
                  & TRIM(NUMBER_TO_VSTRING(LBOUND(DOMAIN_MAPPING2%ADJACENT_DOMAINS(I)%LOCAL_GHOST_RECEIVE_INDICES,1),&
                  & "*",ERR,ERROR))//","//&
                  & TRIM(NUMBER_TO_VSTRING(UBOUND(DOMAIN_MAPPING2%ADJACENT_DOMAINS(I)%LOCAL_GHOST_RECEIVE_INDICES,1),&
                  & "*",ERR,ERROR))//"].",ERR,ERROR,*999)
          MAPPINGS_EQUAL = .FALSE.

          DO J=LBOUND(DOMAIN_MAPPING2%ADJACENT_DOMAINS(I)%LOCAL_GHOST_RECEIVE_INDICES,1),&
            & UBOUND(DOMAIN_MAPPING2%ADJACENT_DOMAINS(I)%LOCAL_GHOST_RECEIVE_INDICES,1)
              CALL FlagWarning(TRIM(NAME)//" Domain mapping ADJACENT_DOMAINS("//&
                    & TRIM(NUMBER_TO_VSTRING(I,"*",ERR,ERROR))//")%LOCAL_GHOST_RECEIVE_INDICES("//&
                    & TRIM(NUMBER_TO_VSTRING(J,"*",ERR,ERROR))//") is : "//&
                    & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING2%ADJACENT_DOMAINS(I)%LOCAL_GHOST_RECEIVE_INDICES(J),&
                    & "*",ERR,ERROR)),ERR,ERROR,*999)
              MAPPINGS_EQUAL = .FALSE.
          ENDDO
      ENDDO
    
        
    
    EXITS("DOMAIN_MAPPINGS_PRINT_NUMB")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_PRINT_NUMB",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_PRINT_NUMB  
  !
  !================================================================================================================================
  !

  !>Initialises the mappings for a domain decomposition. \todo finalise on error.
  SUBROUTINE DOMAIN_MAPPINGS_INITIALISE(DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<A pointer to the domain to initialise the mappings for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(DOMAIN_MAPPING_TYPE) :: NODES_MAPPING_NEW
    TYPE(DOMAIN_MAPPING_TYPE) :: DOFS_MAPPING_NEW 
    TYPE(DOMAIN_MAPPING_TYPE),POINTER :: NODES_MAPPING_NEW_PTR
    INTEGER(INTG) :: DOFS_NUMBER_OF_LOCAL_NEW, DOFS_TOTAL_NUMBER_OF_LOCAL_NEW, DOFS_NUMBER_OF_GLOBAL_NEW, &
      & MaxDepth,MaxArrayLength,I,MyComputationalNodeNumber,NumberComputationalNodes, groupCommunicator


    ENTERS("DOMAIN_MAPPINGS_INITIALISE",ERR,ERROR,*999)

    CALL WorkGroup_GroupCommunicatorGet(decomposition%workGroup,groupCommunicator,err,error,*999)
    CALL WorkGroup_NumberOfGroupNodesGet(decomposition%workGroup,NumberComputationalNodes,err,error,*999)
    CALL WorkGroup_GroupNodeNumberGet(decomposition%workGroup,MyComputationalNodeNumber,err,error,*999)

    ! Replaced by work groups above!!   
    ! get rank count and own rank
    !NumberComputationalNodes=ComputationalEnvironment_NumberOfNodesGet(ERR,ERROR)
    !IF(ERR/=0) GOTO 999
    !MyComputationalNodeNumber=ComputationalEnvironment_NodeNumberGet(ERR,ERROR)
    !IF(ERR/=0) GOTO 999
        

    IF(ASSOCIATED(DOMAIN)) THEN
      IF(ASSOCIATED(DOMAIN%MAPPINGS)) THEN
        CALL FlagError("Domain already has mappings associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(DOMAIN%MAPPINGS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate domain mappings.",ERR,ERROR,*999)
        DOMAIN%MAPPINGS%DOMAIN=>DOMAIN
        NULLIFY(DOMAIN%MAPPINGS%ELEMENTS)
        NULLIFY(DOMAIN%MAPPINGS%NODES)
        NULLIFY(DOMAIN%MAPPINGS%DOFS)
        !Calculate the node and element mappings
        CALL DOMAIN_MAPPINGS_ELEMENTS_INITIALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)
        CALL DOMAIN_MAPPINGS_NODES_INITIALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)
        CALL DOMAIN_MAPPINGS_DOFS_INITIALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)

        IF (USE_NEW_LOCAL_IMPLEMENTATION) THEN
          CALL DOMAIN_MAPPINGS_NODES_CALCULATE(DOMAIN,ERR,ERROR,*999) ! also sets number of dofs
          
          ! output nodes mapping from new implementation

          IF (DIAGNOSTICS3) THEN
            IF (USE_NEW_LOCAL_IMPLEMENTATION) THEN
              MaxDepth = 2
              MaxArrayLength = 100
              IF (NumberComputationalNodes == 1) THEN   ! serial execution
                PRINT *, ""
                PRINT *, "------------------------------------------"
                WRITE(*,'(A,I3,A,I3,A)') "Rank ",I," of ",NumberComputationalNodes, &
                  & ": DOMAIN%MAPPINGS%NODES (new implementation)"
                CALL Print_DOMAIN_MAPPING(DOMAIN%MAPPINGS%NODES, MaxDepth, MaxArrayLength)
              ELSE       
                ! Loop over computational nodes
                DO I = 0,NumberComputationalNodes-1
                  CALL MPI_BARRIER(MPI_COMM_WORLD, ERR)
                  IF (MyComputationalNodeNumber == I) THEN
                    PRINT *, ""
                    IF (MyComputationalNodeNumber == 0) PRINT *, "------------------------------------------"
                    WRITE(*,'(A,I3,A,I3,A)') "Rank ",I," of ",NumberComputationalNodes, &
                      & ": DOMAIN%MAPPINGS%NODES (new implementation)"
                    CALL Print_DOMAIN_MAPPING(DOMAIN%MAPPINGS%NODES, MaxDepth, MaxArrayLength)
                    CALL FLUSH()   ! flush stdout
                  ENDIF
                ENDDO 
                CALL MPI_BARRIER(MPI_COMM_WORLD, ERR)

              ENDIF
              ! output dofs mapping from new implementation
              IF (NumberComputationalNodes == 1) THEN   ! serial execution
                PRINT *, ""
                PRINT *, "------------------------------------------"
                WRITE(*,'(A,I3,A,I3,A)') "Rank ",I," of ",NumberComputationalNodes, &
                  & ": DOMAIN%MAPPINGS%DOFS (new implementation)"
                CALL Print_DOMAIN_MAPPING(DOMAIN%MAPPINGS%DOFS, MaxDepth, MaxArrayLength)
              ELSE       
                ! Loop over computational nodes
                DO I = 0,NumberComputationalNodes-1
                  CALL MPI_BARRIER(MPI_COMM_WORLD, ERR)
                  IF (MyComputationalNodeNumber == I) THEN
                    PRINT *, ""
                    IF (MyComputationalNodeNumber == 0) PRINT *, "------------------------------------------"
                    WRITE(*,'(A,I3,A,I3,A)') "Rank ",I," of ",NumberComputationalNodes, &
                      & ": DOMAIN%MAPPINGS%DOFS (new implementation)"
                    CALL Print_DOMAIN_MAPPING(DOMAIN%MAPPINGS%DOFS, MaxDepth, MaxArrayLength)
                    CALL FLUSH()   ! flush stdout
                  ENDIF
                ENDDO
                CALL MPI_BARRIER(MPI_COMM_WORLD, ERR)
              ENDIF
            ENDIF  ! USE_NEW_LOCAL_IMPLEMENTATION
          ENDIF           
        ENDIF
        
        IF (USE_NEW_LOCAL_IMPLEMENTATION .AND. USE_OLD_GLOBAL_IMPLEMENTATION) THEN        
          ! backup the DOMAIN%MAPPINGS%NODES and %DOFS to compare to old implementation 
          NODES_MAPPING_NEW = DOMAIN%MAPPINGS%NODES
          DOFS_MAPPING_NEW = DOMAIN%MAPPINGS%DOFS 
          DOFS_NUMBER_OF_LOCAL_NEW = DOMAIN%MAPPINGS%DOFS%NUMBER_OF_LOCAL
          DOFS_TOTAL_NUMBER_OF_LOCAL_NEW = DOMAIN%MAPPINGS%DOFS%TOTAL_NUMBER_OF_LOCAL
          DOFS_NUMBER_OF_GLOBAL_NEW = DOMAIN%MAPPINGS%DOFS%NUMBER_OF_GLOBAL
          
          ! reset nodes mapping
          NULLIFY(DOMAIN%MAPPINGS%NODES)
          CALL DOMAIN_MAPPINGS_NODES_INITIALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)

          ! reset dofs mapping
          NULLIFY(DOMAIN%MAPPINGS%DOFS)
          CALL DOMAIN_MAPPINGS_DOFS_INITIALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999) 
        ENDIF
        
        IF (USE_OLD_GLOBAL_IMPLEMENTATION) THEN
          CALL DOMAIN_MAPPINGS_ELEMENTS_CALCULATE(DOMAIN,ERR,ERROR,*999)    ! sets DECOMPOSITION%DOMAIN(.)%MAPPINGS%ELEMENTS, new implementation sets DECOMPOSITION%ELEMENTS_MAPPING in DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE
          CALL DOMAIN_MAPPINGS_NODES_DOFS_CALCULATE(DOMAIN,ERR,ERROR,*999)
        ENDIF
            
        ! output elements and nodes mappings
        IF (DIAGNOSTICS3) THEN
          ! DO NOT output elements mapping from old implementation
          !IF (USE_OLD_GLOBAL_IMPLEMENTATION) THEN
          IF (.FALSE.) THEN
            MaxDepth = 2
            MaxArrayLength = 100
            IF (NumberComputationalNodes == 1) THEN   ! serial execution

              PRINT *, ""
              PRINT *, "------------------------------------------"
              WRITE(*,'(A,I3,A,I3,A)') "Rank ",I," of ",NumberComputationalNodes, &
                    & ": DOMAIN%MAPPINGS%ELEMENTS (old implementation)" 
              CALL Print_DOMAIN_MAPPING(DOMAIN%MAPPINGS%ELEMENTS, MaxDepth, MaxArrayLength)
            ELSE       
              ! Loop over computational nodes
              DO I = 0,NumberComputationalNodes-1
                CALL MPI_BARRIER(MPI_COMM_WORLD, ERR)
                IF (MyComputationalNodeNumber == I) THEN
                  PRINT *, ""
                  IF (MyComputationalNodeNumber == 0) PRINT *, "------------------------------------------"
                  WRITE(*,'(A,I3,A,I3,A)') "Rank ",I," of ",NumberComputationalNodes, &
                    & ": DOMAIN%MAPPINGS%ELEMENTS (old implementation)"
                  CALL Print_DOMAIN_MAPPING(DOMAIN%MAPPINGS%ELEMENTS, MaxDepth, MaxArrayLength)
                  CALL FLUSH()   ! flush stdout
                ENDIF
              ENDDO
              CALL MPI_BARRIER(MPI_COMM_WORLD, ERR)
            ENDIF
          ENDIF   ! USE_OLD_GLOBAL_IMPLEMENTATION
            
          ! output elements mapping from new implementation
          IF (USE_NEW_LOCAL_IMPLEMENTATION) THEN
            MaxDepth = 2
            MaxArrayLength = 100
            IF (NumberComputationalNodes == 1) THEN   ! serial execution

              PRINT *, ""
              PRINT *, "------------------------------------------"
              WRITE(*,'(A,I3,A,I3,A)') "Rank ",I," of ",NumberComputationalNodes, &
                    & ": DECOMPOSITION%ELEMENTS_MAPPING (new implementation)" 
              CALL Print_DOMAIN_MAPPING(DOMAIN%DECOMPOSITION%ELEMENTS_MAPPING, MaxDepth, MaxArrayLength)
            ELSE       
              ! Loop over computational nodes
              DO I = 0,NumberComputationalNodes-1
                CALL MPI_BARRIER(MPI_COMM_WORLD, ERR)
                IF (MyComputationalNodeNumber == I) THEN
                  PRINT *, ""
                  IF (MyComputationalNodeNumber == 0) PRINT *, "------------------------------------------"
                  WRITE(*,'(A,I3,A,I3,A)') "Rank ",I," of ",NumberComputationalNodes, &
                    & ": DECOMPOSITION%ELEMENTS_MAPPING (new implementation)"
                  CALL Print_DOMAIN_MAPPING(DOMAIN%DECOMPOSITION%ELEMENTS_MAPPING, MaxDepth, MaxArrayLength)
                  CALL FLUSH()   ! flush stdout
                ENDIF
              ENDDO
              CALL MPI_BARRIER(MPI_COMM_WORLD, ERR)
            ENDIF
          ENDIF  ! USE_NEW_LOCAL_IMPLEMENTATION
          
          ! output nodes mapping from old implementation
          !IF (USE_OLD_GLOBAL_IMPLEMENTATION) THEN
          IF (.FALSE.) THEN
            MaxDepth = 2
            MaxArrayLength = 100
            IF (NumberComputationalNodes == 1) THEN   ! serial execution
              PRINT *, ""
              PRINT *, "------------------------------------------"
              WRITE(*,'(A,I3,A,I3,A)') "Rank ",I," of ",NumberComputationalNodes, &
                    & ": DOMAIN%MAPPINGS%NODES (old implementation)" 
              CALL Print_DOMAIN_MAPPING(DOMAIN%MAPPINGS%NODES, MaxDepth, MaxArrayLength)
            ELSE       
              ! Loop over computational nodes
              DO I = 0,NumberComputationalNodes-1
                CALL MPI_BARRIER(MPI_COMM_WORLD, ERR)
                IF (MyComputationalNodeNumber == I) THEN
                  PRINT *, ""
                  IF (MyComputationalNodeNumber == 0) PRINT *, "------------------------------------------"
                  WRITE(*,'(A,I3,A,I3,A)') "Rank ",I," of ",NumberComputationalNodes, &
                    & ": DOMAIN%MAPPINGS%NODES (old implementation)"
                  CALL Print_DOMAIN_MAPPING(DOMAIN%MAPPINGS%NODES, MaxDepth, MaxArrayLength)
                  CALL FLUSH()   ! flush stdout
                ENDIF
              ENDDO
              CALL MPI_BARRIER(MPI_COMM_WORLD, ERR)
            ENDIF
          ENDIF   ! USE_OLD_GLOBAL_IMPLEMENTATION

          ! output dofs mapping from old implementation
          !IF (USE_OLD_GLOBAL_IMPLEMENTATION) THEN
          IF (.FALSE.) THEN
            MaxDepth = 2
            MaxArrayLength = 100
            IF (NumberComputationalNodes == 1) THEN   ! serial execution
              PRINT *, ""
              PRINT *, "------------------------------------------"
              WRITE(*,'(A,I3,A,I3,A)') "Rank ",I," of ",NumberComputationalNodes, &
                    & ": DOMAIN%MAPPINGS%DOFS (old implementation)"
              CALL Print_DOMAIN_MAPPING(DOMAIN%MAPPINGS%DOFS, MaxDepth, MaxArrayLength)
            ELSE       
              ! Loop over computational nodes
              DO I = 0,NumberComputationalNodes-1
                CALL MPI_BARRIER(MPI_COMM_WORLD, ERR)
                IF (MyComputationalNodeNumber == I) THEN
                  PRINT *, ""
                  IF (MyComputationalNodeNumber == 0) PRINT *, "------------------------------------------"
                  WRITE(*,'(A,I3,A,I3,A)') "Rank ",I," of ",NumberComputationalNodes, &
                    & ": DOMAIN%MAPPINGS%DOFS (old implementation)"
                  CALL Print_DOMAIN_MAPPING(DOMAIN%MAPPINGS%DOFS, MaxDepth, MaxArrayLength)
                  CALL FLUSH()   ! flush stdout
                ENDIF
              ENDDO
              CALL MPI_BARRIER(MPI_COMM_WORLD, ERR)
            ENDIF
          ENDIF   ! USE_OLD_GLOBAL_IMPLEMENTATION
            
        ENDIF
        
        ! compare elements and nodes mapping
        IF (USE_NEW_LOCAL_IMPLEMENTATION .AND. USE_OLD_GLOBAL_IMPLEMENTATION) THEN        
          ! compare elements and nodes mapping
          !CALL DOMAIN_MAPPINGS_COMPARE(DOMAIN%DECOMPOSITION%ELEMENTS_MAPPING, DOMAIN%MAPPINGS%ELEMENTS, "ELEMENTS",ERR,ERROR,*999)
          !CALL DOMAIN_MAPPINGS_COMPARE(NODES_MAPPING_NEW, DOMAIN%MAPPINGS%NODES, "NODES",ERR,ERROR,*999)
          
          !CALL DOMAIN_MAPPINGS_PRINT_NUMB(NODES_MAPPING_NEW, DOMAIN%MAPPINGS%NODES, "NODES",ERR,ERROR,*999)
          CALL DOMAIN_MAPPINGS_PRINT_NUMB(NODES_MAPPING_NEW, "NODES",ERR,ERROR,*999)

          ! compare number of dofs
          IF (DOFS_NUMBER_OF_LOCAL_NEW /= DOMAIN%MAPPINGS%DOFS%NUMBER_OF_LOCAL) THEN
          CALL FlagWarning("Domain dofs NUMBER_OF_LOCAL is different: "//&
            & TRIM(NUMBER_TO_VSTRING(DOFS_NUMBER_OF_LOCAL_NEW,"*",ERR,ERROR))//" vs. "//&
            & TRIM(NUMBER_TO_VSTRING(DOMAIN%MAPPINGS%DOFS%NUMBER_OF_LOCAL,"*",ERR,ERROR)),ERR,ERROR,*999)
          ENDIF
          IF (DOFS_TOTAL_NUMBER_OF_LOCAL_NEW /= DOMAIN%MAPPINGS%DOFS%TOTAL_NUMBER_OF_LOCAL) THEN
          CALL FlagWarning("Domain dofs TOTAL_NUMBER_OF_LOCAL is different: "//&
            & TRIM(NUMBER_TO_VSTRING(DOFS_TOTAL_NUMBER_OF_LOCAL_NEW,"*",ERR,ERROR))//" vs. "//&
            & TRIM(NUMBER_TO_VSTRING(DOMAIN%MAPPINGS%DOFS%TOTAL_NUMBER_OF_LOCAL,"*",ERR,ERROR)),ERR,ERROR,*999)
          ENDIF
          IF (DOFS_NUMBER_OF_GLOBAL_NEW /= DOMAIN%MAPPINGS%DOFS%NUMBER_OF_GLOBAL) THEN
          CALL FlagWarning("Domain dofs NUMBER_OF_GLOBAL is different: "//&
            & TRIM(NUMBER_TO_VSTRING(DOFS_NUMBER_OF_GLOBAL_NEW,"*",ERR,ERROR))//" vs. "//&
            & TRIM(NUMBER_TO_VSTRING(DOMAIN%MAPPINGS%DOFS%NUMBER_OF_GLOBAL,"*",ERR,ERROR)),ERR,ERROR,*999)
          ENDIF
        ENDIF 

      ENDIF
    ELSE
      CALL FlagError("Domain is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_MAPPINGS_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_INITIALISE


  !
  !================================================================================================================================
  !

  !>Calculates the local/global node mappings for a domain decomposition.
  SUBROUTINE DOMAIN_MAPPINGS_NODES_CALCULATE(DOMAIN,ERR,ERROR,*)
  
    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<A pointer to the domain to calculate the node dofs for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: MyComputationalNodeNumber,NumberComputationalNodes,component_idx,NumberBoundaryAndGhostNodes, &
      & BoundaryNodesListIndex, AdjacentElementLocalNo, BoundaryAndGhostNodeIdx, AdjacentElementIdx, AdjacentElementGlobalNo, &
      & AdjacentDomainIdx, ElementGlobalNo, ElementLocalNo, I, J, nn, NodeGlobalNo, DomainNo, NumberDomains, DomainIdx, &
      & ElementIdx, ArrayIndex, TotalNumberInternalNodes, TotalNumberBoundaryAndGhostNodes, NumberAdjacentDomains, &
      & AdjacentDomain, MPI_IERROR, SharedNode, SharedNodeIdx, NumberSharedNodes, NumberLocalNodes, SharedNodeIdx2, &
      & NumberInteriorAndLocalNodes(2), NumberSharedNodesOnRank, NumberLocalAndAdjacentDomains, DomainToAssignNodesToIdx,&
      & DomainToAssignNodesTo, BoundaryNodeGlobalNo, BoundaryNodesIdx, DomainListBoundaryIdx, DomainListGhostIdx, DomainListIdx, &
      & DomainListInternalIdx, GhostNodeGlobalNo, GhostNodeIdx, GlobalNodeNo, InternalNodeGlobalNo, InternalNodesIdx, &
      & LocalNodeNo, NumberInCurrentSetToClaim, NumberInteriorNodesOnRank, NumberLocalNodesOnRank, NumberToSend, &
      & NumberNodesWithThatSetOfDomains, BoundaryAndGhostNodeIdx2, MaximumNumberToSend, NumberOfNodesToReceive, NumberGhostNodes, &
      & GhostDomain, GhostElementIdx, DomainNo2, AdjacentDomainIdx2, NodeIdx, DerivativeIdx, NodeLocalNo 
      
    INTEGER(INTG), ALLOCATABLE :: BoundaryAndGhostNodes(:), AdjacentDomains(:), SendRequestHandle(:), ReceiveRequestHandle(:), &
      & NumberInteriorAndLocalNodesOnRank(:,:), LocalAndAdjacentDomains(:), BoundaryAndGhostNodesDomain(:), &
      & NumberNodesInDomain(:), GhostNodes(:), InternalNodes(:), BoundaryNodes(:), IntegerArray(:), ReceiveBuffer(:), &
      & SendRequestHandle0(:), SendRequestHandle1(:), GhostNodesDomains(:,:), SendBuffer(:,:), NumberToSendToDomain(:)

    INTEGER(INTG) :: groupCommunicator

    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(MeshComponentTopologyType), POINTER :: TOPOLOGY
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: NODES_MAPPING
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOFS_MAPPING 
    TYPE(LIST_TYPE), POINTER :: InternalNodesList, BoundaryAndGhostNodesList, AdjacentDomainsList, LocalAndAdjacentDomainsList, &
      & BoundaryNodesList, GhostNodesList, GhostNodesDomainsList
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: DomainsOfNodeList(:), SharedNodesList(:), LocalGhostSendIndices(:), &
      & LocalGhostReceiveIndices(:), SendBufferList(:), NewLocalGhostSendIndices, NewLocalGhostReceiveIndices
    REAL(DP) :: OptimalNumberNodesPerDomain, NumberNodes, TotalNumberNodes, NumberNodesAboveOptimum, PortionToDistribute

    LOGICAL :: OnOtherDomain,NodeIsOnBorderBetweenBoundaryAndGhostElementsOnThatDomain,AdacentDomainEntryFound 
    TYPE(DOMAIN_ADJACENT_DOMAIN_TYPE), POINTER :: NewAdjacentDomain 
    !LOGICAL :: DIAGNOSTICS1 = .TRUE.
    !LOGICAL :: DIAGNOSTICS2 = .TRUE.

    ENTERS("DOMAIN_MAPPINGS_NODES_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      IF(ASSOCIATED(DOMAIN%MAPPINGS)) THEN
        IF(ASSOCIATED(DOMAIN%MAPPINGS%NODES)) THEN
          NODES_MAPPING=>DOMAIN%MAPPINGS%NODES
          IF(ASSOCIATED(DOMAIN%DECOMPOSITION)) THEN
            CALL WorkGroup_GroupCommunicatorGet(decomposition%workGroup,groupCommunicator,err,error,*999)
            DECOMPOSITION=>DOMAIN%DECOMPOSITION
            IF(ASSOCIATED(DECOMPOSITION%ELEMENTS_MAPPING)) THEN
              ELEMENTS_MAPPING=>DECOMPOSITION%ELEMENTS_MAPPING
              IF(ASSOCIATED(DOMAIN%MESH)) THEN
                MESH=>DOMAIN%MESH
                component_idx=DOMAIN%MESH_COMPONENT_NUMBER
                TOPOLOGY=>MESH%TOPOLOGY(component_idx)%PTR

                CALL WorkGroup_NumberOfGroupNodesGet(decomposition%workGroup,NumberComputationalNodes,err,error,*999)
                CALL WorkGroup_GroupNodeNumberGet(decomposition%workGroup,MyComputationalNodeNumber,err,error,*999)

            ! Replaced by work groups above!!                   
                !NumberComputationalNodes=ComputationalEnvironment_NumberOfNodesGet(ERR,ERROR)
                !IF(ERR/=0) GOTO 999
                !MyComputationalNodeNumber=ComputationalEnvironment_NodeNumberGet(ERR,ERROR)
                !IF(ERR/=0) GOTO 999
                
                IF(DIAGNOSTICS1) THEN
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,&
                    & "  NumberComputationalNodes: ",NumberComputationalNodes,ERR,ERROR,*999)
                    
                ENDIF
                
                ! create list for internal nodes
                NULLIFY(InternalNodesList)
                CALL LIST_CREATE_START(InternalNodesList,ERR,ERROR,*999)
                CALL LIST_DATA_TYPE_SET(InternalNodesList,LIST_INTG_TYPE,ERR,ERROR,*999)
                CALL LIST_INITIAL_SIZE_SET(InternalNodesList,ELEMENTS_MAPPING%NUMBER_OF_LOCAL*2,ERR,ERROR,*999)
                CALL LIST_CREATE_FINISH(InternalNodesList,ERR,ERROR,*999)
                
                ! determine interior and boundary nodes
                ! loop over interior elements
                DO ElementIdx = ELEMENTS_MAPPING%INTERNAL_START,ELEMENTS_MAPPING%INTERNAL_FINISH
                  ElementLocalNo = ELEMENTS_MAPPING%DOMAIN_LIST(ElementIdx)
                  ElementGlobalNo = ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(ElementLocalNo)
                  
                  ! loop over nodes of element
                  DO nn=1,TOPOLOGY%ELEMENTS%ELEMENTS(ElementGlobalNo)%BASIS%NUMBER_OF_NODES
                    NodeGlobalNo=TOPOLOGY%ELEMENTS%ELEMENTS(ElementGlobalNo)%MESH_ELEMENT_NODES(nn)
                    
                    ! add node to internal nodes list
                    CALL LIST_ITEM_ADD(InternalNodesList, NodeGlobalNo,ERR,ERROR,*999)
                  ENDDO
                ENDDO
                
                ! loop over boundary and ghost elements
                DO ElementIdx = ELEMENTS_MAPPING%BOUNDARY_START,ELEMENTS_MAPPING%BOUNDARY_FINISH
                  ElementLocalNo = ELEMENTS_MAPPING%DOMAIN_LIST(ElementIdx)
                  ElementGlobalNo = ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(ElementLocalNo)
                  
                  ! loop over nodes of element
                  DO nn=1,TOPOLOGY%ELEMENTS%ELEMENTS(ElementGlobalNo)%BASIS%NUMBER_OF_NODES
                    NodeGlobalNo=TOPOLOGY%ELEMENTS%ELEMENTS(ElementGlobalNo)%MESH_ELEMENT_NODES(nn)
                      

                    IF (DIAGNOSTICS2) &
                      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"node global "//TRIM(NumberToVString(NodeGlobalNo,"*",ERR,ERROR)), &
                      & ERR,ERROR,*999)
                      !PRINT *, MyComputationalNodeNumber, ": node global ",NodeGlobalNo 
      
                    OnOtherDomain = .FALSE.
                      
                    ! loop over adjacent elements of this node
                    DO AdjacentElementIdx=1,TOPOLOGY%NODES%NODES(NodeGlobalNo)%numberOfSurroundingElements
                      AdjacentElementGlobalNo=TOPOLOGY%NODES%NODES(NodeGlobalNo)%surroundingElements(AdjacentElementIdx)
                      

                      IF (DIAGNOSTICS2) &
                        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"     adjacent element global " //&
                          & TRIM(NumberToVString(NodeGlobalNo,"*",ERR,ERROR)),ERR,ERROR,*999)
                        !PRINT *, MyComputationalNodeNumber, ":     adjacent element global ",AdjacentElementGlobalNo 
                    
                      ! determine if adjacent element is on a different domain, that is the case if it is a ghost
                      IF (SORTED_ARRAY_CONTAINS_ELEMENT(ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP( &
                        & ELEMENTS_MAPPING%GHOST_START:ELEMENTS_MAPPING%GHOST_FINISH),AdjacentElementGlobalNo,ERR,ERROR)) THEN
                        

                        IF (DIAGNOSTICS2) &
                          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"     on other domain",ERR,ERROR,*999)
                          !PRINT *, MyComputationalNodeNumber, ":     on other domain" 
                    
                        OnOtherDomain = .TRUE.
                        EXIT
                      ENDIF
                    
                    ENDDO
                    
                     !IF (DIAGNOSTICS2) PRINT *, MyComputationalNodeNumber, ": is node internal? ",(.NOT. OnOtherDomain) 
                    
                    IF (.NOT. OnOtherDomain) THEN
                      ! add node to internal nodes list
                      CALL LIST_ITEM_ADD(InternalNodesList, NodeGlobalNo,ERR,ERROR,*999)
                    ENDIF
                    
                  ENDDO
                  
                ENDDO
                
                ! Sort list by global node number and store in local nodes storage
                CALL LIST_REMOVE_DUPLICATES(InternalNodesList,ERR,ERROR,*999)
                CALL LIST_DETACH_AND_DESTROY(InternalNodesList,NODES_MAPPING%NUMBER_OF_INTERNAL,&
                  & InternalNodes,ERR,ERROR,*999)
                
                IF (DIAGNOSTICS1) THEN
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,&
                    & "number internal nodes: ",NODES_MAPPING%NUMBER_OF_INTERNAL,ERR,ERROR,*999)
                  DO I=1,NODES_MAPPING%NUMBER_OF_INTERNAL
                    CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  global no. ", &
                      & InternalNodes(I),ERR,ERROR,*999)
                  ENDDO
                ENDIF    
                
                ! create list for boundary nodes (global node numbers)
                NULLIFY(BoundaryAndGhostNodesList)
                CALL LIST_CREATE_START(BoundaryAndGhostNodesList,ERR,ERROR,*999)
                CALL LIST_DATA_TYPE_SET(BoundaryAndGhostNodesList,LIST_INTG_TYPE,ERR,ERROR,*999)
                CALL LIST_INITIAL_SIZE_SET(BoundaryAndGhostNodesList,ELEMENTS_MAPPING%NUMBER_OF_LOCAL,ERR,ERROR,*999)
                CALL LIST_CREATE_FINISH(BoundaryAndGhostNodesList,ERR,ERROR,*999)
                
                ! loop over boundary elements
                DO ElementIdx = ELEMENTS_MAPPING%BOUNDARY_START,ELEMENTS_MAPPING%BOUNDARY_FINISH
                  ElementLocalNo = ELEMENTS_MAPPING%DOMAIN_LIST(ElementIdx)
                  ElementGlobalNo = ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(ElementLocalNo)
                  
                  ! loop over nodes of element
                  DO nn=1,TOPOLOGY%ELEMENTS%ELEMENTS(ElementGlobalNo)%BASIS%NUMBER_OF_NODES
                    NodeGlobalNo=TOPOLOGY%ELEMENTS%ELEMENTS(ElementGlobalNo)%MESH_ELEMENT_NODES(nn)
            
                    ! if global node is not contained in internal nodes list
                    IF (.NOT. SORTED_ARRAY_CONTAINS_ELEMENT(InternalNodes,NodeGlobalNo,ERR,ERROR)) THEN
                      
                      ! add node to boundary nodes list
                      CALL LIST_ITEM_ADD(BoundaryAndGhostNodesList,NodeGlobalNo,ERR,ERROR,*999)
                    ENDIF
                  ENDDO
                ENDDO
                
                ! Sort list by global node number and remove duplicates
                CALL LIST_REMOVE_DUPLICATES(BoundaryAndGhostNodesList,ERR,ERROR,*999)
                
                ! count number of nodes that lie on boundary and are shared with other domains
                CALL LIST_DETACH_AND_DESTROY(BoundaryAndGhostNodesList,NumberBoundaryAndGhostNodes,IntegerArray,&
                  & ERR,ERROR,*999)
                BoundaryAndGhostNodes = IntegerArray(1:NumberBoundaryAndGhostNodes)
                DEALLOCATE(IntegerArray)
            
                IF (DIAGNOSTICS1) THEN
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,&
                    & "number boundary and ghost nodes: ",NumberBoundaryAndGhostNodes,ERR,ERROR,*999)
                  DO I=1,NumberBoundaryAndGhostNodes
                    CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,&
                      & "  global no. ",BoundaryAndGhostNodes(I),ERR,ERROR,*999)
                  ENDDO
                ENDIF
                
                ! allocate lists of foreign domains in which boundary nodes are present
                ALLOCATE(DomainsOfNodeList(NumberBoundaryAndGhostNodes),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate domain list.",ERR,ERROR,*999)
                
                DO I=1, NumberBoundaryAndGhostNodes
                  NULLIFY(DomainsOfNodeList(I)%PTR)
                  CALL LIST_CREATE_START(DomainsOfNodeList(I)%PTR,ERR,ERROR,*999)
                  CALL LIST_DATA_TYPE_SET(DomainsOfNodeList(I)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                  CALL LIST_INITIAL_SIZE_SET(DomainsOfNodeList(I)%PTR,7,ERR,ERROR,*999)
                  CALL LIST_CREATE_FINISH(DomainsOfNodeList(I)%PTR,ERR,ERROR,*999)
                ENDDO
                
                ! fill DomainsOfNodeList with all foreign domains for each boundary node
                
                ! loop over boundary and ghost nodes
                DO BoundaryAndGhostNodeIdx=1,NumberBoundaryAndGhostNodes
                  NodeGlobalNo=BoundaryAndGhostNodes(BoundaryAndGhostNodeIdx)
                

                  IF (DIAGNOSTICS2) &
                      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"node global "//TRIM(NumberToVString(NodeGlobalNo,"*",ERR,ERROR)), &
                      & ERR,ERROR,*999)
                      ! PRINT *, MyComputationalNodeNumber, ": node global ",NodeGlobalNo 
      
                  ! determine the other domains which have this node as boundary node
                  
                  ! loop over adjacent elements of this node
                  DO AdjacentElementIdx=1,TOPOLOGY%NODES%NODES(NodeGlobalNo)%numberOfSurroundingElements
                    AdjacentElementGlobalNo=TOPOLOGY%NODES%NODES(NodeGlobalNo)%surroundingElements(AdjacentElementIdx)
                    

                    IF (DIAGNOSTICS2) &
                      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"     adjacent element global "// &
                        & TRIM(NumberToVString(NodeGlobalNo,"*",ERR,ERROR)),ERR,ERROR,*999)
                      ! PRINT *, MyComputationalNodeNumber, ":     adjacent element global ",AdjacentElementGlobalNo 
                    
                    ! determine local (ghost) number
                    ! if adjacent element is among the ghost elements, i.e. on a different domain (and not a boundary element)
                    IF (SORTED_ARRAY_CONTAINS_ELEMENT_INDEX(ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP( &
                      & ELEMENTS_MAPPING%GHOST_START:ELEMENTS_MAPPING%GHOST_FINISH),AdjacentElementGlobalNo,ArrayIndex,&
                      & ERR,ERROR)) THEN
                      
                      AdjacentElementLocalNo = ELEMENTS_MAPPING%GHOST_START + ArrayIndex-1
                      

                      IF (DIAGNOSTICS2) &
                        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"     adjacent element local "// &
                          & TRIM(NumberToVString(AdjacentElementLocalNo,"*",ERR,ERROR)),ERR,ERROR,*999)
                        ! PRINT *, MyComputationalNodeNumber, ":     adjacent element local ",AdjacentElementLocalNo 
                      
                      ! determine domain of AdjacentElementLocalNo
                      DO AdjacentDomainIdx=1,ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
                        IF (SORTED_ARRAY_CONTAINS_ELEMENT(ELEMENTS_MAPPING%ADJACENT_DOMAINS(AdjacentDomainIdx)%&
                          & LOCAL_GHOST_RECEIVE_INDICES,AdjacentElementLocalNo,ERR,ERROR)) THEN
                          
                          DomainNo = ELEMENTS_MAPPING%ADJACENT_DOMAINS(AdjacentDomainIdx)%DOMAIN_NUMBER

                          IF (DIAGNOSTICS2) &
                            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"     on domain "// &
                              & TRIM(NumberToVString(DomainNo,"*",ERR,ERROR)),ERR,ERROR,*999)
                            ! PRINT *, MyComputationalNodeNumber, ":     on domain ",DomainNo 
                          
                          ! now it was found that NodeGlobalNo is also on domain AdjacentDomainIdx, add to list
                          CALL LIST_ITEM_ADD(DomainsOfNodeList(BoundaryAndGhostNodeIdx)%PTR,DomainNo,ERR,ERROR,*999)
                          EXIT
                        ENDIF
                      ENDDO  ! AdjacentDomainIdx
                    ELSE 

                      IF (DIAGNOSTICS2) &
                        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"     is on same domain",ERR,ERROR,*999)
                        ! PRINT *, MyComputationalNodeNumber, ":     is on same domain" 
                    ENDIF
                  ENDDO  ! AdjacentElementIdx
                ENDDO  ! BoundaryAndGhostNodeIdx
              
                
                NULLIFY(AdjacentDomainsList)
                CALL LIST_CREATE_START(AdjacentDomainsList,ERR,ERROR,*999)
                CALL LIST_DATA_TYPE_SET(AdjacentDomainsList,LIST_INTG_TYPE,ERR,ERROR,*999)
                CALL LIST_CREATE_FINISH(AdjacentDomainsList,ERR,ERROR,*999)
                
                ! remove duplicate entries for the domains for the nodes
                DO I=1, NumberBoundaryAndGhostNodes
                  CALL LIST_REMOVE_DUPLICATES(DomainsOfNodeList(I)%PTR,ERR,ERROR,*999)
                  
                  ! add adjacent domains to AdjacentDomainsList
                  CALL LIST_NUMBER_OF_ITEMS_GET(DomainsOfNodeList(I)%PTR,NumberDomains,ERR,ERROR,*999)
                  DO J=1,NumberDomains
                    CALL LIST_ITEM_GET(DomainsOfNodeList(I)%PTR,J,AdjacentDomain,ERR,ERROR,*999)
                    CALL LIST_ITEM_ADD(AdjacentDomainsList,AdjacentDomain,ERR,ERROR,*999)
                  ENDDO
                ENDDO
                
                ! determine distinct list of adjacent domains
                CALL LIST_REMOVE_DUPLICATES(AdjacentDomainsList,ERR,ERROR,*999)
                
                ! create sorted list that contains local domain and all adjacent domains
                NULLIFY(LocalAndAdjacentDomainsList)
                CALL LIST_CREATE_START(LocalAndAdjacentDomainsList,ERR,ERROR,*999)
                CALL LIST_DATA_TYPE_SET(LocalAndAdjacentDomainsList,LIST_INTG_TYPE,ERR,ERROR,*999)
                CALL LIST_CREATE_FINISH(LocalAndAdjacentDomainsList,ERR,ERROR,*999)
                CALL LIST_ITEM_ADD(LocalAndAdjacentDomainsList,MyComputationalNodeNumber,ERR,ERROR,*999)
                CALL LIST_APPENDLIST(LocalAndAdjacentDomainsList,AdjacentDomainsList,ERR,ERROR,*999)
                CALL LIST_SORT(LocalAndAdjacentDomainsList,ERR,ERROR,*999)
                
                CALL LIST_DETACH_AND_DESTROY(AdjacentDomainsList,NumberAdjacentDomains,IntegerArray,ERR,ERROR,*999)
                AdjacentDomains = IntegerArray(1:NumberAdjacentDomains)
                DEALLOCATE(IntegerArray)
                
                CALL LIST_DETACH_AND_DESTROY(LocalAndAdjacentDomainsList,NumberLocalAndAdjacentDomains,IntegerArray,&
                  & ERR,ERROR,*999)
                LocalAndAdjacentDomains = IntegerArray(1:NumberLocalAndAdjacentDomains)
                DEALLOCATE(IntegerArray)
                
                IF (DIAGNOSTICS1) THEN

                  CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE, TRIM(NumberToVString(NumberAdjacentDomains,"*",ERR,ERROR)) // &
                    & " AdjacentDomains",ERR,ERROR,*999)
                  CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE, TRIM(NumberToVString(NumberLocalAndAdjacentDomains,"*",ERR,ERROR)) // &
                    & " LocalAndAdjacentDomains",ERR,ERROR,*999)
                  !PRINT *, MyComputationalNodeNumber, ": ",NumberAdjacentDomains," AdjacentDomains:",AdjacentDomains
                  !PRINT *, MyComputationalNodeNumber, ": ",NumberLocalAndAdjacentDomains," LocalAndAdjacentDomains:", &
                  !  & LocalAndAdjacentDomains 
                ENDIF
              
                ! create list for each adjacent domain with shared nodes
                ALLOCATE(SharedNodesList(NumberAdjacentDomains),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate shared nodes list.",ERR,ERROR,*999)
                  
                DO AdjacentDomainIdx=1,NumberAdjacentDomains
                  NULLIFY(SharedNodesList(AdjacentDomainIdx)%PTR)
                  CALL LIST_CREATE_START(SharedNodesList(AdjacentDomainIdx)%PTR,ERR,ERROR,*999)
                  CALL LIST_DATA_TYPE_SET(SharedNodesList(AdjacentDomainIdx)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                  CALL LIST_CREATE_FINISH(SharedNodesList(AdjacentDomainIdx)%PTR,ERR,ERROR,*999)
                ENDDO
                
                ! fill SharedNodesList
                DO I=1,NumberBoundaryAndGhostNodes
                
                  ! loop over adjacent domains of node I
                  CALL LIST_NUMBER_OF_ITEMS_GET(DomainsOfNodeList(I)%PTR,NumberDomains,ERR,ERROR,*999)
                  DO J=1,NumberDomains
                    CALL LIST_ITEM_GET(DomainsOfNodeList(I)%PTR,J,AdjacentDomain,ERR,ERROR,*999)
                    
                    ! get index of SharedNodesList
                    IF (SORTED_ARRAY_CONTAINS_ELEMENT_INDEX(AdjacentDomains(1:NumberAdjacentDomains),&
                      & AdjacentDomain,AdjacentDomainIdx,ERR,ERROR)) THEN
                    
                      ! add node to sharedNodeslist of domain AdjacentDomainIdx
                      CALL LIST_ITEM_ADD(SharedNodesList(AdjacentDomainIdx)%PTR,I,ERR,ERROR,*999)
                    ENDIF
                  ENDDO
                ENDDO
                
                IF (DIAGNOSTICS1) THEN
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,&
                    & "number boundary and ghost nodes: ",NumberBoundaryAndGhostNodes,ERR,ERROR,*999)
                  DO I=1,NumberBoundaryAndGhostNodes
                    CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,&
                      & "  global no. ",BoundaryAndGhostNodes(I),ERR,ERROR,*999)
                    CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  domains:",ERR,ERROR,*999)
                    
                    CALL LIST_NUMBER_OF_ITEMS_GET(DomainsOfNodeList(I)%PTR,NumberDomains,ERR,ERROR,*999)
                    DO DomainIdx=1,NumberDomains
                      CALL LIST_ITEM_GET(DomainsOfNodeList(I)%PTR,DomainIdx,DomainNo,ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,&
                        & "    ",DomainNo,ERR,ERROR,*999)
                    ENDDO
                  ENDDO
                  
                  CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"adjacent domains with their nodes: ",ERR,ERROR,*999)
                  
                  DO AdjacentDomainIdx=1,NumberAdjacentDomains                
                    AdjacentDomain = AdjacentDomains(AdjacentDomainIdx)
                    CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,&
                        & "    rank ",AdjacentDomain,ERR,ERROR,*999)
                    CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    nodes: ",ERR,ERROR,*999)
                    CALL LIST_NUMBER_OF_ITEMS_GET(SharedNodesList(AdjacentDomainIdx)%PTR,NumberSharedNodes,ERR,ERROR,*999)
                    DO SharedNodeIdx = 1,NumberSharedNodes
                      CALL LIST_ITEM_GET(SharedNodesList(AdjacentDomainIdx)%PTR,SharedNodeIdx,BoundaryAndGhostNodeIdx,&
                        & ERR,ERROR,*999)
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    ",BoundaryAndGhostNodes(BoundaryAndGhostNodeIdx), &
                        & ERR,ERROR,*999)
                  
                    ENDDO
                  ENDDO
                ENDIF
                
                ! initialize NODES_MAPPING%ADJACENT_DOMAINS
                NODES_MAPPING%NUMBER_OF_ADJACENT_DOMAINS = NumberAdjacentDomains
                
                ! allocate ADJACENT_DOMAINS
                ALLOCATE(NODES_MAPPING%ADJACENT_DOMAINS(NumberAdjacentDomains), STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate adjacent domains array with size "//&
                  & TRIM(NUMBER_TO_VSTRING(NumberAdjacentDomains,"*",ERR,ERROR))//".",ERR,ERROR,*999)
                
                DO AdjacentDomainIdx=1,NumberAdjacentDomains                
                  AdjacentDomain = AdjacentDomains(AdjacentDomainIdx)
                  NODES_MAPPING%ADJACENT_DOMAINS(AdjacentDomainIdx)%DOMAIN_NUMBER = AdjacentDomain
                  NODES_MAPPING%ADJACENT_DOMAINS(AdjacentDomainIdx)%NUMBER_OF_SEND_GHOSTS = 0
                  NODES_MAPPING%ADJACENT_DOMAINS(AdjacentDomainIdx)%NUMBER_OF_RECEIVE_GHOSTS = 0
                  NODES_MAPPING%ADJACENT_DOMAINS(AdjacentDomainIdx)%NUMBER_OF_FURTHER_LINKED_GHOSTS = 0
                ENDDO
                
                ! determine total number of nodes on all ranks
                
                ! count number of nodes where boundary/ghost nodes are counted as the resp. fraction such that the sum is 1 for all processes
                NumberNodes = INT(NODES_MAPPING%NUMBER_OF_INTERNAL)
                
                DO I=1,NumberBoundaryAndGhostNodes
                  CALL LIST_NUMBER_OF_ITEMS_GET(DomainsOfNodeList(I)%PTR,NumberDomains,ERR,ERROR,*999)
                  NumberNodes = NumberNodes + 1.0/(NumberDomains+1)
                ENDDO
                
                ! allreduce number of nodes
                CALL MPI_ALLREDUCE(NumberNodes,TotalNumberNodes,1,MPI_DOUBLE,MPI_SUM, &
                  & groupCommunicator,MPI_IERROR)
                CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
                NODES_MAPPING%NUMBER_OF_GLOBAL = NINT(TotalNumberNodes)
                

                IF (DIAGNOSTICS1) THEN
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"NODES_MAPPING%NUMBER_OF_INTERNAL: ",&
                    & NODES_MAPPING%NUMBER_OF_INTERNAL,ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"NumberBoundaryAndGhostNodes: ", &
                    & NumberBoundaryAndGhostNodes,ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"NumberNodes: ",NumberNodes,ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"TotalNumberNodes: ",TotalNumberNodes,ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"NODES_MAPPING%NUMBER_OF_GLOBAL: ", &
                    & NODES_MAPPING%NUMBER_OF_GLOBAL,ERR,ERROR,*999)
                  !PRINT *, MyComputationalNodeNumber, ": NODES_MAPPING%NUMBER_OF_INTERNAL: ",NODES_MAPPING%NUMBER_OF_INTERNAL
                  !PRINT *, MyComputationalNodeNumber, ": NumberBoundaryAndGhostNodes: ",NumberBoundaryAndGhostNodes
                  !PRINT *, MyComputationalNodeNumber, ": NumberNodes: ",NumberNodes
                  !PRINT *, MyComputationalNodeNumber, ": TotalNumberNodes: ",TotalNumberNodes
                  !PRINT *, MyComputationalNodeNumber, ": NODES_MAPPING%NUMBER_OF_GLOBAL: ",NODES_MAPPING%NUMBER_OF_GLOBAL 
                ENDIF
                
                ! compute average number per domain
                NODES_MAPPING%NUMBER_OF_DOMAINS = ELEMENTS_MAPPING%NUMBER_OF_DOMAINS
                OptimalNumberNodesPerDomain = REAL(NODES_MAPPING%NUMBER_OF_GLOBAL) / NODES_MAPPING%NUMBER_OF_DOMAINS
                
                IF (DIAGNOSTICS1) THEN
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,&
                    & "global number of nodes: ",NODES_MAPPING%NUMBER_OF_GLOBAL,ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,&
                    & "OptimalNumberNodesPerDomain: ",OptimalNumberNodesPerDomain,ERR,ERROR,*999)
                ENDIF
                
                ! exchange number of local nodes among adjacent ranks
                NumberLocalNodes = NODES_MAPPING%NUMBER_OF_INTERNAL + NumberBoundaryAndGhostNodes
                NumberInteriorAndLocalNodes(1) = NODES_MAPPING%NUMBER_OF_INTERNAL
                NumberInteriorAndLocalNodes(2) = NumberLocalNodes
                          
                ! allocate request handles
                ALLOCATE(SendRequestHandle(NumberAdjacentDomains+1), STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate SendRequestHandle array with size "//&
                  & TRIM(NUMBER_TO_VSTRING(NumberAdjacentDomains+1,"*",ERR,ERROR))//".",ERR,ERROR,*999)
                
                ALLOCATE(ReceiveRequestHandle(NumberAdjacentDomains+1), STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate ReceiveRequestHandle array with size "//&
                  & TRIM(NUMBER_TO_VSTRING(NumberAdjacentDomains+1,"*",ERR,ERROR))//".",ERR,ERROR,*999)
                
                ! allocate receive buffer (2 entries for interior and totally stored local nodes, adjacent domains + local domain)
                ALLOCATE(NumberInteriorAndLocalNodesOnRank(2,NumberAdjacentDomains+1), STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate NumberInteriorAndLocalNodesOnRank array with size "//&
                  & TRIM(NUMBER_TO_VSTRING(2*(NumberAdjacentDomains+1),"*",ERR,ERROR))//".",ERR,ERROR,*999)
                
                ! commit send commands
                DO AdjacentDomainIdx=1,NumberAdjacentDomains+1 
                  AdjacentDomain = LocalAndAdjacentDomains(AdjacentDomainIdx)
                
                  CALL MPI_ISEND(NumberInteriorAndLocalNodes, 2, MPI_INT, AdjacentDomain, 0, &
                    & groupCommunicator, SendRequestHandle(AdjacentDomainIdx), MPI_IERROR)
                  CALL MPI_ERROR_CHECK("MPI_ISEND",MPI_IERROR,ERR,ERROR,*999)
                ENDDO
                
                ! commit receive commands
                DO AdjacentDomainIdx=1,NumberAdjacentDomains+1 
                  AdjacentDomain = LocalAndAdjacentDomains(AdjacentDomainIdx)
                  CALL MPI_IRECV(NumberInteriorAndLocalNodesOnRank(:,AdjacentDomainIdx), 2, MPI_INT, AdjacentDomain, 0, &
                    & groupCommunicator, ReceiveRequestHandle(AdjacentDomainIdx), MPI_IERROR)
                  CALL MPI_ERROR_CHECK("MPI_IRECV",MPI_IERROR,ERR,ERROR,*999)            
                ENDDO
                
                ! wait for all communication to finish
                CALL MPI_WAITALL(NumberAdjacentDomains+1, SendRequestHandle, MPI_STATUSES_IGNORE, MPI_IERROR)
                CALL MPI_ERROR_CHECK("MPI_WAITALL",MPI_IERROR,ERR,ERROR,*999)
                
                CALL MPI_WAITALL(NumberAdjacentDomains+1, ReceiveRequestHandle, MPI_STATUSES_IGNORE, MPI_IERROR)
                CALL MPI_ERROR_CHECK("MPI_WAITALL",MPI_IERROR,ERR,ERROR,*999)
        
                IF (DIAGNOSTICS1) THEN
                  DO AdjacentDomainIdx=1,NumberAdjacentDomains+1 
                    AdjacentDomain = LocalAndAdjacentDomains(AdjacentDomainIdx)

                    CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE, "rank "//TRIM(NumberToVString(AdjacentDomain,"*",ERR,ERROR)) // &
                      & ", interior:"//TRIM(NumberToVString(AdjacentDomain,"*",ERR,ERROR)) // &
                      & ", total:"//TRIM(NumberToVString(NumberInteriorAndLocalNodesOnRank(2,AdjacentDomainIdx),"*",ERR,ERROR)) //&
                      & ", AdjacentDomainIdx="//TRIM(NumberToVString(AdjacentDomainIdx,"*",ERR,ERROR)),ERR,ERROR,*999)
                    !PRINT *, MyComputationalNodeNumber, ": rank ",AdjacentDomain,&
                    !  & " interior:",NumberInteriorAndLocalNodesOnRank(1,AdjacentDomainIdx), &
                    !  & ", total:",NumberInteriorAndLocalNodesOnRank(2,AdjacentDomainIdx),", AdjacentDomainIdx=",AdjacentDomainIdx 
                  ENDDO
                ENDIF
        
                ! assign boundary nodes to domains
                ! allocate BoundaryAndGhostNodesDomain
                ALLOCATE(BoundaryAndGhostNodesDomain(NumberBoundaryAndGhostNodes), STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate BoundaryAndGhostNodesDomain array with size "//&
                  & TRIM(NUMBER_TO_VSTRING(NumberBoundaryAndGhostNodes,"*",ERR,ERROR))//".",ERR,ERROR,*999)
                BoundaryAndGhostNodesDomain = -1
        
                ALLOCATE(NumberNodesInDomain(NumberAdjacentDomains+1), STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate NumberNodesInDomain array with size "//&
                  & TRIM(NUMBER_TO_VSTRING(NumberAdjacentDomains+1,"*",ERR,ERROR))//".",ERR,ERROR,*999)
                NumberNodesInDomain = 0
      
                ! add local domain to domain list of each node
                DO BoundaryAndGhostNodeIdx=1,NumberBoundaryAndGhostNodes
                  CALL LIST_ITEM_ADD(DomainsOfNodeList(BoundaryAndGhostNodeIdx)%PTR,MyComputationalNodeNumber,ERR,ERROR,*999)
                  CALL LIST_SORT(DomainsOfNodeList(BoundaryAndGhostNodeIdx)%PTR,ERR,ERROR,*999)
                ENDDO
        
                ! loop over adjacent domains
                DO AdjacentDomainIdx=1,NumberAdjacentDomains
                  AdjacentDomain = AdjacentDomains(AdjacentDomainIdx)
                          
                  IF (DIAGNOSTICS2) THEN
                    CALL LIST_NUMBER_OF_ITEMS_GET(SharedNodesList(AdjacentDomainIdx)%PTR,NumberSharedNodes,ERR,ERROR,*999)

                    CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE, ": domain "//TRIM(NumberToVString(AdjacentDomain,"*",ERR,ERROR)) // &
                      & " has "//TRIM(NumberToVString(NumberSharedNodes,"*",ERR,ERROR)) // &
                      & " nodes ",ERR,ERROR,*999)
                    !PRINT *, MyComputationalNodeNumber, ": domain ", AdjacentDomain, " has ",NumberSharedNodes," nodes " 
                    DO SharedNodeIdx = 1,NumberSharedNodes
                      CALL LIST_ITEM_GET(SharedNodesList(AdjacentDomainIdx)%PTR,SharedNodeIdx,BoundaryAndGhostNodeIdx,&
                        & ERR,ERROR,*999)
                      GlobalNodeNo = BoundaryAndGhostNodes(BoundaryAndGhostNodeIdx)

                      !PRINT *, MyComputationalNodeNumber, ":     global no. ", GlobalNodeNo
                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE, "     global no. ",GlobalNodeNo,ERR,ERROR,*999) 
                    ENDDO
                        
                  ENDIF
                          
                  ! loop over adjacent nodes on that adjacent domain
                  CALL LIST_NUMBER_OF_ITEMS_GET(SharedNodesList(AdjacentDomainIdx)%PTR,NumberSharedNodes,ERR,ERROR,*999)
                  DO SharedNodeIdx = 1,NumberSharedNodes
                    CALL LIST_ITEM_GET(SharedNodesList(AdjacentDomainIdx)%PTR,SharedNodeIdx,BoundaryAndGhostNodeIdx,&
                      & ERR,ERROR,*999)
                    GlobalNodeNo = BoundaryAndGhostNodes(BoundaryAndGhostNodeIdx)
                    
                    IF (DIAGNOSTICS2) THEN
                      CALL LIST_NUMBER_OF_ITEMS_GET(DomainsOfNodeList(BoundaryAndGhostNodeIdx)%PTR, NumberDomains, ERR,ERROR,*999)

                      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE, " consider node "//&
                        & TRIM(NumberToVString(GlobalNodeNo,"*",ERR,ERROR)) // &
                        & " which is shared by the "//TRIM(NumberToVString(NumberDomains,"*",ERR,ERROR)) // &
                        & " domains ",ERR,ERROR,*999)
                      !PRINT *, MyComputationalNodeNumber, ": consider node ",GlobalNodeNo," which is shared by the ",&
                      !  NumberDomains," domains " 
                      DO I=1,NumberDomains
                        CALL LIST_ITEM_GET(DomainsOfNodeList(BoundaryAndGhostNodeIdx)%PTR,I,DomainNo,ERR,ERROR,*999)

                        !PRINT *, MyComputationalNodeNumber, ":   ",DomainNo
                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE, "   ",DomainNo,ERR,ERROR,*999) 
                      ENDDO
                    ENDIF
                    
                    IF (DIAGNOSTICS2) THEN

                      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE, " node assigned to domain ",BoundaryAndGhostNodesDomain( &
                        & BoundaryAndGhostNodeIdx),ERR,ERROR,*999)
                      !PRINT *, MyComputationalNodeNumber, ": node assigned to domain ", &
                      !  & BoundaryAndGhostNodesDomain(BoundaryAndGhostNodeIdx) 
                    ENDIF
                      
                    ! if node was not yet assigned to a domain
                    IF (BoundaryAndGhostNodesDomain(BoundaryAndGhostNodeIdx) == -1) THEN
                      ! note: this node is shared by the domains in DomainsOfNodeList(BoundaryAndGhostNodeIdx) (includes local domain)
                      
                      NumberNodesWithThatSetOfDomains = 0
                      
                      ! loop over all further nodes to find other nodes that are shared by the same domains, count total number of them
                      DO SharedNodeIdx2 = SharedNodeIdx,NumberSharedNodes
                        CALL LIST_ITEM_GET(SharedNodesList(AdjacentDomainIdx)%PTR,SharedNodeIdx2,BoundaryAndGhostNodeIdx2,&
                          & ERR,ERROR,*999)
                          
                        IF (DIAGNOSTICS2) THEN
                          GlobalNodeNo = BoundaryAndGhostNodes(BoundaryAndGhostNodeIdx2)
                    

                          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE, " further node ",GlobalNodeNo,ERR,ERROR,*999)
                          !PRINT *, MyComputationalNodeNumber, ": further node ", GlobalNodeNo 
                          CALL LIST_NUMBER_OF_ITEMS_GET(DomainsOfNodeList(BoundaryAndGhostNodeIdx2)%PTR,NumberDomains,&
                            & ERR,ERROR,*999)
                          DO I=1,NumberDomains
                            CALL LIST_ITEM_GET(DomainsOfNodeList(BoundaryAndGhostNodeIdx2)%PTR,I,DomainNo,ERR,ERROR,*999)

                            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE, "     domain set  ",DomainNo,ERR,ERROR,*999)
                            !PRINT *, MyComputationalNodeNumber, ":     domain set  ",DomainNo 
                          ENDDO
                        ENDIF
                        
                        IF (BoundaryAndGhostNodesDomain(BoundaryAndGhostNodeIdx2) == -1) THEN
                          
                          IF (LIST_EQUAL(DomainsOfNodeList(BoundaryAndGhostNodeIdx)%PTR, &
                            & DomainsOfNodeList(BoundaryAndGhostNodeIdx2)%PTR,ERR,ERROR)) THEN

                            IF (DIAGNOSTICS2) CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE, "     domain set is equal",ERR,ERROR,*999)
                            !PRINT *, MyComputationalNodeNumber, ":     domain set is equal" 
                          ENDIF
                        ENDIF
                      ENDDO
                      IF (DIAGNOSTICS2) THEN

                        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE, " there are "//&
                          & TRIM(NumberToVString(NumberNodesWithThatSetOfDomains,"*",ERR,ERROR)) // &
                          & " nodes with that set of domains. Now assign the nodes to domains. ",ERR,ERROR,*999)
                        !PRINT *, MyComputationalNodeNumber, ": there are ", NumberNodesWithThatSetOfDomains, &
                        !  & " nodes with that set of domains. Now assign the nodes to domains. DIAGNOSTICS2:",DIAGNOSTICS2 
                      ENDIF
                      
                      ! get first of the shared domains
                      NumberNodesInDomain = 0
                      DomainToAssignNodesToIdx = 1
                      CALL LIST_ITEM_GET(DomainsOfNodeList(BoundaryAndGhostNodeIdx2)%PTR,DomainToAssignNodesToIdx, &
                        & DomainToAssignNodesTo,ERR,ERROR,*999)
                      CALL LIST_NUMBER_OF_ITEMS_GET(DomainsOfNodeList(BoundaryAndGhostNodeIdx2)%PTR,NumberDomains,ERR,ERROR,*999)
                      
                      ! again loop over all further nodes with the same domains, assign them to domains
                      DO SharedNodeIdx2 = SharedNodeIdx,NumberSharedNodes
                        CALL LIST_ITEM_GET(SharedNodesList(AdjacentDomainIdx)%PTR,SharedNodeIdx2,BoundaryAndGhostNodeIdx2,&
                            & ERR,ERROR,*999)
                        
                        IF (BoundaryAndGhostNodesDomain(BoundaryAndGhostNodeIdx2) /= -1) CYCLE
                        
                        IF (LIST_EQUAL(DomainsOfNodeList(BoundaryAndGhostNodeIdx)%PTR, &
                          & DomainsOfNodeList(BoundaryAndGhostNodeIdx2)%PTR,ERR,ERROR)) THEN
                          
                          DO
                          
                            ! find index for NumberInteriorAndLocalNodesOnRank that corresponds to domain DomainToAssignNodesTo
                            AdjacentDomainIdx2 = 1
                            DO I=1,NumberAdjacentDomains+1 
                              AdjacentDomain = LocalAndAdjacentDomains(I)
                              IF (AdjacentDomain == DomainToAssignNodesTo) THEN
                                AdjacentDomainIdx2 = I
                                EXIT
                              ENDIF
                            ENDDO
                          
                            ! compute the number of nodes that this domain will get from the current set of nodes
                            NumberLocalNodesOnRank = NumberInteriorAndLocalNodesOnRank(2,AdjacentDomainIdx2)
                            NumberInteriorNodesOnRank = NumberInteriorAndLocalNodesOnRank(1,AdjacentDomainIdx2)
                            
                            NumberNodesAboveOptimum = NumberLocalNodesOnRank - OptimalNumberNodesPerDomain
                            NumberSharedNodesOnRank = NumberLocalNodesOnRank - NumberInteriorNodesOnRank
                              
                            ! compute the portion of each node that adjacent domain should get from other domains
                            ! E.g. if PortionToDistribute=0.4 that means that 4 out of 10 shared nodes of the adjacent domain should be assigned to other domains and 6 should remain their own
                            PortionToDistribute = NumberNodesAboveOptimum / NumberSharedNodesOnRank
                            
                            ! compute the number of nodes from the current set of nodes that have the same domains that should belong to adjacent domain
                            NumberInCurrentSetToClaim = NINT((1.0 - PortionToDistribute) * NumberNodesWithThatSetOfDomains)
                            
                            
                            IF (DIAGNOSTICS2) THEN

                              CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE, " Domain "//&
                                & TRIM(NumberToVString(DomainToAssignNodesTo,"*",ERR,ERROR)),ERR,ERROR,*999)
                              CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE, "   has "//&
                                & TRIM(NumberToVString(NumberInteriorNodesOnRank,"*",ERR,ERROR))//" interior, "//&
                                & TRIM(NumberToVString(NumberLocalNodesOnRank,"*",ERR,ERROR))//" total, "//&
                                & TRIM(NumberToVString(NumberSharedNodesOnRank,"*",ERR,ERROR))//" shared nodes, "//&
                                & "AdjacentDomainIdx2="//&
                                & TRIM(NumberToVString(AdjacentDomainIdx2,"*",ERR,ERROR)),ERR,ERROR,*999)
                              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE, "   NumberNodesAboveOptimum=", &
                                & NumberNodesAboveOptimum,ERR,ERROR,*999)
                              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE, "   PortionToDistribute=",PortionToDistribute, &
                                & ERR,ERROR,*999)
                              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE, "   NumberInCurrentSetToClaim=", &
                                & NumberInCurrentSetToClaim,ERR,ERROR,*999)
                              CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE, "   already got "// &
                                & TRIM(NumberToVString(NumberNodesInDomain(DomainToAssignNodesToIdx),"*",ERR,ERROR))// &
                                & " nodes",ERR,ERROR,*999)
                            
                              !PRINT *, MyComputationalNodeNumber, ": Domain ", DomainToAssignNodesTo
                              !PRINT *, MyComputationalNodeNumber, ":   has ", NumberInteriorNodesOnRank, " interior, ",&
                              !  & NumberLocalNodesOnRank, " total, ", NumberSharedNodesOnRank, " shared nodes",&
                              !  & ", AdjacentDomainIdx2=",AdjacentDomainIdx2
                              !PRINT *, MyComputationalNodeNumber, ":   NumberNodesAboveOptimum=",NumberNodesAboveOptimum
                              !PRINT *, MyComputationalNodeNumber, ":   PortionToDistribute=",PortionToDistribute
                              !PRINT *, MyComputationalNodeNumber, ":   NumberInCurrentSetToClaim=",NumberInCurrentSetToClaim
                              !PRINT *, MyComputationalNodeNumber, ":   already got ",&
                              !  & NumberNodesInDomain(DomainToAssignNodesToIdx)," nodes" 
                            ENDIF
                            
                            ! if there weren't enough nodes assigned to adjacent domain DomainToAssignNodesToIdx yet, the domain which will claim this node is found
                            IF (NumberNodesInDomain(DomainToAssignNodesToIdx) < NumberInCurrentSetToClaim) THEN
                              EXIT
                            ELSE
                              ! if there are no further domains, exit trying to find the next domain
                              IF (DomainToAssignNodesToIdx == NumberDomains) EXIT
                            
                              ! now that domain has got enough nodes, get the next domain
                              DomainToAssignNodesToIdx = DomainToAssignNodesToIdx + 1
                              CALL LIST_ITEM_GET(DomainsOfNodeList(BoundaryAndGhostNodeIdx2)%PTR,DomainToAssignNodesToIdx, &
                                & DomainToAssignNodesTo,ERR,ERROR,*999)
                            ENDIF
                          ENDDO
                          
                          ! assign node SharedNodeIdx2 to domain DomainToAssignNodesTo
                          BoundaryAndGhostNodesDomain(BoundaryAndGhostNodeIdx2) = DomainToAssignNodesTo
                            
                          ! increment count of nodes that this domain has claimed in the current set
                          NumberNodesInDomain(DomainToAssignNodesToIdx) = NumberNodesInDomain(DomainToAssignNodesToIdx) + 1
                          
                          IF (DIAGNOSTICS2) THEN
                            GlobalNodeNo = BoundaryAndGhostNodes(BoundaryAndGhostNodeIdx2)

                            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE, " assign node "//&
                              & TRIM(NumberToVString(GlobalNodeNo,"*",ERR,ERROR))//" to domain "//&
                              & TRIM(NumberToVString(DomainToAssignNodesTo,"*",ERR,ERROR))//", now it has "//&
                              & TRIM(NumberToVString(NumberNodesInDomain(DomainToAssignNodesToIdx),"*",ERR,ERROR))//&
                              & " nodes",ERR,ERROR,*999)
                              
                            !PRINT *, MyComputationalNodeNumber, ": assign node ",GlobalNodeNo," to domain ",DomainToAssignNodesTo, &
                            !  & ", now it has ",NumberNodesInDomain(DomainToAssignNodesToIdx)," nodes" 
                          ENDIF
                          
                        ENDIF
                      ENDDO  ! SharedNodeIdx2                                          
                    ELSE
                      IF (DIAGNOSTICS2) THEN

                        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE, " node already assigned to domain "//&
                          & TRIM(NumberToVString(BoundaryAndGhostNodesDomain(BoundaryAndGhostNodeIdx),"*",ERR,ERROR))&
                          & ,ERR,ERROR,*999)
                        !PRINT *, MyComputationalNodeNumber, ": node already assigned to domain ", &
                        !  & BoundaryAndGhostNodesDomain(BoundaryAndGhostNodeIdx) 
                      ENDIF
                    ENDIF  ! node not yet assigned to domain
                    
                  ENDDO  ! SharedNodeIdx
                ENDDO  ! AdjacentDomainIdx
                
                ! fill local numbering arrays
                
                ! create list to hold ghost nodes (BoundaryAndGhostNodeIdx)
                NULLIFY(GhostNodesList)
                CALL LIST_CREATE_START(GhostNodesList,ERR,ERROR,*999)
                CALL LIST_DATA_DIMENSION_SET(GhostNodesList,1,ERR,ERROR,*999)
                CALL LIST_DATA_TYPE_SET(GhostNodesList,LIST_INTG_TYPE,ERR,ERROR,*999)
                CALL LIST_INITIAL_SIZE_SET(GhostNodesList,CEILING(ELEMENTS_MAPPING%NUMBER_OF_LOCAL*0.2),ERR,ERROR,*999)
                CALL LIST_CREATE_FINISH(GhostNodesList,ERR,ERROR,*999)
                
                ! create list to hold boundary nodes (BoundaryAndGhostNodeIdx)
                NULLIFY(BoundaryNodesList)
                CALL LIST_CREATE_START(BoundaryNodesList,ERR,ERROR,*999)
                CALL LIST_DATA_DIMENSION_SET(BoundaryNodesList,1,ERR,ERROR,*999)
                CALL LIST_DATA_TYPE_SET(BoundaryNodesList,LIST_INTG_TYPE,ERR,ERROR,*999)
                CALL LIST_INITIAL_SIZE_SET(BoundaryNodesList,CEILING(ELEMENTS_MAPPING%NUMBER_OF_LOCAL*0.2),ERR,ERROR,*999)
                CALL LIST_CREATE_FINISH(BoundaryNodesList,ERR,ERROR,*999)
                 
                ! count number of boundary and ghost nodes
                DO BoundaryAndGhostNodeIdx=1,NumberBoundaryAndGhostNodes
                  DomainNo = BoundaryAndGhostNodesDomain(BoundaryAndGhostNodeIdx)
                  
                  ! if node is a boundary nodes and stays on own domain
                  IF (DomainNo == MyComputationalNodeNumber) THEN
                    CALL LIST_ITEM_ADD(BoundaryNodesList,BoundaryAndGhostNodeIdx,ERR,ERROR,*999)
                  ELSE
                    CALL LIST_ITEM_ADD(GhostNodesList,BoundaryAndGhostNodeIdx,ERR,ERROR,*999)
                  ENDIF
                ENDDO
                
                ! extract arrays of local boundary and ghost nodes
                CALL LIST_DETACH_AND_DESTROY(BoundaryNodesList, NODES_MAPPING%NUMBER_OF_BOUNDARY, BoundaryNodes, ERR,ERROR,*999)
                CALL LIST_DETACH_AND_DESTROY(GhostNodesList, NumberGhostNodes, GhostNodes, ERR,ERROR,*999)
                
                ! create list of ghost nodes with their home domain
                NULLIFY(GhostNodesDomainsList)
                CALL LIST_CREATE_START(GhostNodesDomainsList,ERR,ERROR,*999)
                CALL LIST_DATA_DIMENSION_SET(GhostNodesDomainsList,2,ERR,ERROR,*999)
                CALL LIST_DATA_TYPE_SET(GhostNodesDomainsList,LIST_INTG_TYPE,ERR,ERROR,*999)
                CALL LIST_INITIAL_SIZE_SET(GhostNodesDomainsList,MAX(1,2*NumberGhostNodes),ERR,ERROR,*999) 
                CALL LIST_CREATE_FINISH(GhostNodesDomainsList,ERR,ERROR,*999)
                
                DO GhostNodeIdx=1,NumberGhostNodes
                  BoundaryAndGhostNodeIdx = GhostNodes(GhostNodeIdx)
                  GhostNodeGlobalNo = BoundaryAndGhostNodes(BoundaryAndGhostNodeIdx)
                  GhostDomain = BoundaryAndGhostNodesDomain(BoundaryAndGhostNodeIdx)
                  CALL LIST_ITEM_ADD(GhostNodesDomainsList,[GhostNodeGlobalNo,GhostDomain],ERR,ERROR,*999)
                ENDDO
                
                ! exchange outer ghost nodes that are not nodes of a boundary element
                ! allocate send buffers
                ALLOCATE(SendBufferList(ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS), STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate SendBufferList array with size "//&
                  & TRIM(NUMBER_TO_VSTRING(ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS,"*",ERR,ERROR))//".",ERR,ERROR,*999)
                
                ALLOCATE(SendRequestHandle0(ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS), STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate SendRequestHandle0 array with size "//&
                  & TRIM(NUMBER_TO_VSTRING(ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS,"*",ERR,ERROR))//".",ERR,ERROR,*999)
                
                ALLOCATE(SendRequestHandle1(ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS), STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate SendRequestHandle1 array with size "//&
                  & TRIM(NUMBER_TO_VSTRING(ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS,"*",ERR,ERROR))//".",ERR,ERROR,*999)
                
                DO I=1,ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
                  NULLIFY(SendBufferList(I)%PTR)
                  CALL LIST_CREATE_START(SendBufferList(I)%PTR,ERR,ERROR,*999)
                  CALL LIST_DATA_DIMENSION_SET(SendBufferList(I)%PTR,1,ERR,ERROR,*999)
                  CALL LIST_DATA_TYPE_SET(SendBufferList(I)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                  CALL LIST_INITIAL_SIZE_SET(SendBufferList(I)%PTR,10,ERR,ERROR,*999)
                  CALL LIST_CREATE_FINISH(SendBufferList(I)%PTR,ERR,ERROR,*999)
                ENDDO
                
                ! add own nodes that are ghost nodes on other domains and do not lie on the border between boundary and ghost elements at the other domain 
                ! to send buffer to the send to that domain
                
                ! loop over boundary elements
                DO ElementIdx = ELEMENTS_MAPPING%BOUNDARY_START,ELEMENTS_MAPPING%GHOST_FINISH
                  ElementLocalNo = ELEMENTS_MAPPING%DOMAIN_LIST(ElementIdx)
                  ElementGlobalNo = ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(ElementLocalNo)
                  
                    
                  ! loop over nodes of element
                  DO nn=1,TOPOLOGY%ELEMENTS%ELEMENTS(ElementGlobalNo)%BASIS%NUMBER_OF_NODES
                    NodeGlobalNo=TOPOLOGY%ELEMENTS%ELEMENTS(ElementGlobalNo)%MESH_ELEMENT_NODES(nn)
                   
                    ! if node is at border between boundary and ghost elements, i.e. it is boundary node at some domain
                    IF (SORTED_ARRAY_CONTAINS_ELEMENT_INDEX(BoundaryAndGhostNodes, NodeGlobalNo,BoundaryAndGhostNodeIdx,&
                      & ERR,ERROR)) THEN
                      
                     
                      ! get the domain that owns the node                      
                      DomainNo = BoundaryAndGhostNodesDomain(BoundaryAndGhostNodeIdx)
                      
                      IF (DIAGNOSTICS1) THEN

                        !PRINT *, MyComputationalNodeNumber, ": boundary element local ", ElementLocalNo, &
                        !  & ", node global ",NodeGlobalNo, " is at border boundary/ghost element, ", &
                        !  & "DomainNo: ",DomainNo 
                      ENDIF
                      
                      ! only add to send buffer if it is owned by own domain
                      IF (DomainNo == MyComputationalNodeNumber) THEN
                    
                        ! do not add to send buffer for domains where it lies on the border between boundary and ghost elements (these are already known at the other domain)
                        ! loop over all domains that have the current element as ghost element and only consider those domains where the node does not have the domain in DomainsOfNodeList, sorry if confusing
                        
                        ! loop over domains where boundary element is ghost element
                        DO DomainIdx=1,ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
                        
                          DO GhostElementIdx=1,ELEMENTS_MAPPING%ADJACENT_DOMAINS(DomainIdx)%NUMBER_OF_SEND_GHOSTS
                           ! if element is a ghost on domain DomainIdx
                            IF (ELEMENTS_MAPPING%ADJACENT_DOMAINS(DomainIdx)%LOCAL_GHOST_SEND_INDICES(GhostElementIdx) &
                              & == ElementLocalNo) THEN
                              DomainNo = ELEMENTS_MAPPING%ADJACENT_DOMAINS(DomainIdx)%DOMAIN_NUMBER
                              
                              IF (DIAGNOSTICS1) THEN
                                 !PRINT *, MyComputationalNodeNumber, ":      boundary element is ghost on domain ",DomainNo 
                              ENDIF
                        
                              NodeIsOnBorderBetweenBoundaryAndGhostElementsOnThatDomain = .FALSE.
                              
                              ! loop over domains that share the node 
                              CALL LIST_NUMBER_OF_ITEMS_GET(DomainsOfNodeList(BoundaryAndGhostNodeIdx)%PTR,NumberDomains, &
                                & ERR,ERROR,*999)
                              DO I=1,NumberDomains
                                CALL LIST_ITEM_GET(DomainsOfNodeList(BoundaryAndGhostNodeIdx)%PTR,I, &
                                  & DomainNo2,ERR,ERROR,*999)
                                  
                                IF (DIAGNOSTICS1) THEN
                                   !PRINT *, MyComputationalNodeNumber, ":        node is on domain ",DomainNo2 
                                ENDIF
                                IF (DomainNo == DomainNo2) THEN
                                  NodeIsOnBorderBetweenBoundaryAndGhostElementsOnThatDomain = .TRUE. 
                                  EXIT
                                ENDIF
                              ENDDO  ! I
                              
                              IF (DIAGNOSTICS1) THEN

                                !PRINT *, MyComputationalNodeNumber, ":        -> NodeIsOnBorderBetweenBoundary"//&
                                !  & "AndGhostElementsOnThatDomain=",NodeIsOnBorderBetweenBoundaryAndGhostElementsOnThatDomain 
                              ENDIF
                              
                              IF (.NOT.NodeIsOnBorderBetweenBoundaryAndGhostElementsOnThatDomain) THEN
                              
                                ! add node to list to be send to domain of boundary element
                                CALL LIST_ITEM_ADD(SendBufferList(DomainIdx)%PTR, NodeGlobalNo, ERR,ERROR,*999)
                              
                              ENDIF
                            ENDIF
                          ENDDO ! GhostElementIdx, Send Ghosts
                          
                          
                          DO GhostElementIdx=1,ELEMENTS_MAPPING%ADJACENT_DOMAINS(DomainIdx)%NUMBER_OF_RECEIVE_GHOSTS
                            ! if element is a ghost on domain DomainIdx
                            IF (ELEMENTS_MAPPING%ADJACENT_DOMAINS(DomainIdx)%LOCAL_GHOST_RECEIVE_INDICES(GhostElementIdx) &
                              & == ElementLocalNo) THEN
                              DomainNo = ELEMENTS_MAPPING%ADJACENT_DOMAINS(DomainIdx)%DOMAIN_NUMBER
                              
                              IF (DIAGNOSTICS1) THEN
                                 !PRINT *, MyComputationalNodeNumber, ":      boundary element is (recv) ghost on domain ",DomainNo 
                              ENDIF
                        
                              NodeIsOnBorderBetweenBoundaryAndGhostElementsOnThatDomain = .FALSE.
                              
                              ! loop over domains that share the node 
                              CALL LIST_NUMBER_OF_ITEMS_GET(DomainsOfNodeList(BoundaryAndGhostNodeIdx)%PTR,NumberDomains, &
                                & ERR,ERROR,*999)
                              DO I=1,NumberDomains
                                CALL LIST_ITEM_GET(DomainsOfNodeList(BoundaryAndGhostNodeIdx)%PTR,I, &
                                  & DomainNo2,ERR,ERROR,*999)
                                  
                                IF (DIAGNOSTICS1) THEN
                                   !PRINT *, MyComputationalNodeNumber, ":        node is on domain ",DomainNo2 
                                ENDIF
                                IF (DomainNo == DomainNo2) THEN
                                  NodeIsOnBorderBetweenBoundaryAndGhostElementsOnThatDomain = .TRUE. 
                                  EXIT
                                ENDIF
                              ENDDO  ! I
                              
                              IF (DIAGNOSTICS1) THEN

                                !PRINT *, MyComputationalNodeNumber, ":        -> NodeIsOnBorderBetweenBoundary"//&
                                !  & "AndGhostElementsOnThatDomain=",NodeIsOnBorderBetweenBoundaryAndGhostElementsOnThatDomain 
                              ENDIF
                              
                              IF (.NOT.NodeIsOnBorderBetweenBoundaryAndGhostElementsOnThatDomain) THEN
                              
                                ! add node to list to be send to domain of boundary element
                                CALL LIST_ITEM_ADD(SendBufferList(DomainIdx)%PTR, NodeGlobalNo, ERR,ERROR,*999)
                              
                              ENDIF
                            ENDIF
                          ENDDO ! GhostElementIdx, Receive Ghosts
                          
                          
                          DO GhostElementIdx=1,ELEMENTS_MAPPING%ADJACENT_DOMAINS(DomainIdx)%NUMBER_OF_FURTHER_LINKED_GHOSTS
                            ! if element is a ghost on domain DomainIdx
                           IF (ELEMENTS_MAPPING%ADJACENT_DOMAINS(DomainIdx)%LOCAL_GHOST_FURTHER_INDICES(GhostElementIdx) &
                              & == ElementLocalNo) THEN
                              DomainNo = ELEMENTS_MAPPING%ADJACENT_DOMAINS(DomainIdx)%DOMAIN_NUMBER
                              
                              IF (DIAGNOSTICS1) THEN
                                 !PRINT *, MyComputationalNodeNumber, ":      boundary element is (further) ghost on domain ",DomainNo 
                              ENDIF
                        
                              NodeIsOnBorderBetweenBoundaryAndGhostElementsOnThatDomain = .FALSE.
                              
                              ! loop over domains that share the node 
                              CALL LIST_NUMBER_OF_ITEMS_GET(DomainsOfNodeList(BoundaryAndGhostNodeIdx)%PTR,NumberDomains, &
                                & ERR,ERROR,*999)
                              DO I=1,NumberDomains
                                CALL LIST_ITEM_GET(DomainsOfNodeList(BoundaryAndGhostNodeIdx)%PTR,I, &
                                  & DomainNo2,ERR,ERROR,*999)
                                 
                                IF (DIAGNOSTICS1) THEN
                                   !PRINT *, MyComputationalNodeNumber, ":        node is on domain ",DomainNo2 
                                ENDIF
                                IF (DomainNo == DomainNo2) THEN
                                  NodeIsOnBorderBetweenBoundaryAndGhostElementsOnThatDomain = .TRUE. 
                                  EXIT
                                ENDIF
                              ENDDO  ! I
                              
                              IF (DIAGNOSTICS1) THEN

                                !PRINT *, MyComputationalNodeNumber, ":        -> NodeIsOnBorderBetweenBoundary"//&
                                !  & "AndGhostElementsOnThatDomain=",NodeIsOnBorderBetweenBoundaryAndGhostElementsOnThatDomain 
                              ENDIF
                              
                              IF (.NOT.NodeIsOnBorderBetweenBoundaryAndGhostElementsOnThatDomain) THEN
                              
                                ! add node to list to be send to domain of boundary element
                                CALL LIST_ITEM_ADD(SendBufferList(DomainIdx)%PTR, NodeGlobalNo, ERR,ERROR,*999)
                              
                              ENDIF
                            ENDIF
                          ENDDO ! GhostElementIdx, further linked ghosts
                          
                        ENDDO ! DomainIdx
                      ENDIF
                    ELSE
                      ! if node is not at border between boundary and ghost elements, i.e. it is not in BoundaryAndGhostNodes
                      
                      IF (DIAGNOSTICS1) THEN

                       ! PRINT *, MyComputationalNodeNumber, ": boundary element local ", ElementLocalNo, &
                       !   & ", node global ",NodeGlobalNo, " is not at border boundary/ghost element, ", &
                       !   & "BoundaryAndGhostNodes: ", BoundaryAndGhostNodes 
                      ENDIF
                      
                      ! get domains where boundary element is ghost element
                      DO DomainIdx=1,ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
                        DO GhostElementIdx=1,ELEMENTS_MAPPING%ADJACENT_DOMAINS(DomainIdx)%NUMBER_OF_SEND_GHOSTS
                          ! if element is a ghost on domain DomainIdx
                          IF (ELEMENTS_MAPPING%ADJACENT_DOMAINS(DomainIdx)%LOCAL_GHOST_SEND_INDICES(GhostElementIdx) &
                            & == ElementLocalNo) THEN
                            DomainNo = ELEMENTS_MAPPING%ADJACENT_DOMAINS(DomainIdx)%DOMAIN_NUMBER
                            
                            IF (DIAGNOSTICS1) THEN
                               !PRINT *, MyComputationalNodeNumber, ":    ghost on domain ", DomainNo, ", DomainIdx=",DomainIdx 
                            ENDIF
                            
                            ! add node to list to be send to domain of boundary element
                            CALL LIST_ITEM_ADD(SendBufferList(DomainIdx)%PTR, NodeGlobalNo, ERR,ERROR,*999)
                            
                            EXIT
                          ENDIF
                        ENDDO  ! GhostElementIdx
                        
                        
                      ENDDO  ! DomainIdx
                    ENDIF
                  ENDDO
                ENDDO
                
                ! remove duplicates and sort list
                MaximumNumberToSend = 0
                DO DomainIdx=1,ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
                  DomainNo = ELEMENTS_MAPPING%ADJACENT_DOMAINS(DomainIdx)%DOMAIN_NUMBER
                  
                  CALL LIST_REMOVE_DUPLICATES(SendBufferList(DomainIdx)%PTR,ERR,ERROR,*999)
                  CALL LIST_NUMBER_OF_ITEMS_GET(SendBufferList(DomainIdx)%PTR,NumberToSend,ERR,ERROR,*999)
                  MaximumNumberToSend = MAX(MaximumNumberToSend,NumberToSend)
                ENDDO
                
                ALLOCATE(SendBuffer(MaximumNumberToSend,ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS), STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate SendBuffer array with size "//&
                  & TRIM(NUMBER_TO_VSTRING(MaximumNumberToSend,"*",ERR,ERROR))//"x"//&
                  & TRIM(NUMBER_TO_VSTRING(ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS,"*",ERR,ERROR))//".",ERR,ERROR,*999)
                
                ALLOCATE(NumberToSendToDomain(ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS), STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate NumberToSendToDomain array with size "//&
                  & TRIM(NUMBER_TO_VSTRING(ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS,"*",ERR,ERROR))//".",ERR,ERROR,*999)
                
                ! commit send calls
                DO DomainIdx=1,ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
                  DomainNo = ELEMENTS_MAPPING%ADJACENT_DOMAINS(DomainIdx)%DOMAIN_NUMBER
                  CALL LIST_DETACH_AND_DESTROY(SendBufferList(DomainIdx)%PTR,NumberToSendToDomain(DomainIdx),&
                    & IntegerArray,ERR,ERROR,*999)
                  SendBuffer(1:NumberToSendToDomain(DomainIdx),DomainIdx) = IntegerArray(1:NumberToSendToDomain(DomainIdx))
                  DEALLOCATE(IntegerArray)
                  
                  IF (DIAGNOSTICS1) THEN

                    !PRINT *, MyComputationalNodeNumber, ": send to ", DomainNo, ": ", &
                    !  & SendBuffer(1:NumberToSendToDomain(DomainIdx),DomainIdx) 
                  ENDIF
                  
                  ! number of nodes to be send to domain
                  CALL MPI_ISEND(NumberToSendToDomain(DomainIdx),1,MPI_INT,DomainNo,0, &
                    & groupCommunicator, SendRequestHandle0(DomainIdx), MPI_IERROR)
                  CALL MPI_ERROR_CHECK("MPI_ISEND",MPI_IERROR,ERR,ERROR,*999)
                  
                  ! actual global node numbers to send
                  IF (NumberToSendToDomain(DomainIdx) > 0) THEN
                    CALL MPI_ISEND(SendBuffer(1:NumberToSendToDomain(DomainIdx),DomainIdx),NumberToSendToDomain(DomainIdx), &
                      & MPI_INT,DomainNo,1,groupCommunicator, SendRequestHandle1(DomainIdx), MPI_IERROR)
                    CALL MPI_ERROR_CHECK("MPI_ISEND",MPI_IERROR,ERR,ERROR,*999)
                  ENDIF
                ENDDO
                
                ! receive number of nodes to receive
                DO DomainIdx=1,ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
                  DomainNo = ELEMENTS_MAPPING%ADJACENT_DOMAINS(DomainIdx)%DOMAIN_NUMBER
                          
                  CALL MPI_RECV(NumberOfNodesToReceive,1,MPI_INT,DomainNo,0, &
                    & groupCommunicator,MPI_STATUS_IGNORE, MPI_IERROR)
                  CALL MPI_ERROR_CHECK("MPI_RECV",MPI_IERROR,ERR,ERROR,*999)
            
                  IF (NumberOfNodesToReceive>0) THEN  
                    ALLOCATE(ReceiveBuffer(NumberOfNodesToReceive), STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate ReceiveBuffer array with size "//&
                      & TRIM(NUMBER_TO_VSTRING(NumberOfNodesToReceive,"*",ERR,ERROR))//".",ERR,ERROR,*999)
                  
                    CALL MPI_RECV(ReceiveBuffer,NumberOfNodesToReceive,MPI_INT,DomainNo,1, &
                      & groupCommunicator,MPI_STATUS_IGNORE,MPI_IERROR)
                    CALL MPI_ERROR_CHECK("MPI_RECV",MPI_IERROR,ERR,ERROR,*999)
              
                    ! store received ghost nodes with domain in GhostNodesDomainsList
                    DO I=1,NumberOfNodesToReceive
                      GhostNodeGlobalNo = ReceiveBuffer(I)
                      GhostDomain = DomainNo
                      CALL LIST_ITEM_ADD(GhostNodesDomainsList,[GhostNodeGlobalNo,GhostDomain],ERR,ERROR,*999)
                    ENDDO
                    
                    IF (DIAGNOSTICS1) THEN

                      !PRINT *, MyComputationalNodeNumber, ": Received ", NumberOfNodesToReceive, " ghosts from rank ",DomainNo, &
                      !  & ": ",ReceiveBuffer 
                    ENDIF
                    
                    DEALLOCATE(ReceiveBuffer)
              
                  ENDIF
                ENDDO
                
                ! wait for all communication to finish
                !CALL MPI_WAITALL(ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS, SendRequestHandle0, MPI_STATUSES_IGNORE, MPI_IERROR)
                !CALL MPI_ERROR_CHECK("MPI_WAITALL",MPI_IERROR,ERR,ERROR,*999)
                
                !CALL MPI_WAITALL(ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS, SendRequestHandle1, MPI_STATUSES_IGNORE, MPI_IERROR)
                !CALL MPI_ERROR_CHECK("MPI_WAITALL",MPI_IERROR,ERR,ERROR,*999)
                
                CALL LIST_SORT(GhostNodesDomainsList,ERR,ERROR,*999)
                CALL LIST_DETACH_AND_DESTROY(GhostNodesDomainsList,NODES_MAPPING%NUMBER_OF_GHOST,GhostNodesDomains,ERR,ERROR,*999)
                
                ! store number of local and total number of local nodes
                NODES_MAPPING%NUMBER_OF_LOCAL = NODES_MAPPING%NUMBER_OF_INTERNAL + NODES_MAPPING%NUMBER_OF_BOUNDARY 
                NODES_MAPPING%TOTAL_NUMBER_OF_LOCAL = NODES_MAPPING%NUMBER_OF_LOCAL + NODES_MAPPING%NUMBER_OF_GHOST
                
                NODES_MAPPING%INTERNAL_START=1
                NODES_MAPPING%INTERNAL_FINISH=NODES_MAPPING%INTERNAL_START+NODES_MAPPING%NUMBER_OF_INTERNAL-1
                NODES_MAPPING%BOUNDARY_START=NODES_MAPPING%INTERNAL_FINISH+1
                NODES_MAPPING%BOUNDARY_FINISH=NODES_MAPPING%BOUNDARY_START+NODES_MAPPING%NUMBER_OF_BOUNDARY-1
                NODES_MAPPING%GHOST_START=NODES_MAPPING%BOUNDARY_FINISH+1
                NODES_MAPPING%GHOST_FINISH=NODES_MAPPING%GHOST_START+NODES_MAPPING%NUMBER_OF_GHOST-1
                
                ! allocate LOCAL_TO_GLOBAL_MAP
                ALLOCATE(NODES_MAPPING%LOCAL_TO_GLOBAL_MAP(NODES_MAPPING%TOTAL_NUMBER_OF_LOCAL), STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate nodes mapping local to global map with size "//&
                  & TRIM(NUMBER_TO_VSTRING(NODES_MAPPING%TOTAL_NUMBER_OF_LOCAL,"*",ERR,ERROR))//".",ERR,ERROR,*999)
                
                ! allocate DOMAIN_LIST
                ALLOCATE(NODES_MAPPING%DOMAIN_LIST(NODES_MAPPING%TOTAL_NUMBER_OF_LOCAL), STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate nodes mapping domain list with size "//&
                  & TRIM(NUMBER_TO_VSTRING(NODES_MAPPING%TOTAL_NUMBER_OF_LOCAL,"*",ERR,ERROR))//".",ERR,ERROR,*999)
                
                ! create lists for local ghost send and receive indices
                ALLOCATE(LocalGhostSendIndices(NumberAdjacentDomains), STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate LocalGhostSendIndices with size "//&
                  & TRIM(NUMBER_TO_VSTRING(NumberAdjacentDomains,"*",ERR,ERROR))//".",ERR,ERROR,*999)
                
                ALLOCATE(LocalGhostReceiveIndices(NumberAdjacentDomains), STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate LocalGhostReceiveIndices with size "//&
                  & TRIM(NUMBER_TO_VSTRING(NumberAdjacentDomains,"*",ERR,ERROR))//".",ERR,ERROR,*999)
                
                DO AdjacentDomainIdx=1,NumberAdjacentDomains

                  ! create list for local ghost send indices
                  NULLIFY(LocalGhostSendIndices(AdjacentDomainIdx)%PTR)
                  CALL LIST_CREATE_START(LocalGhostSendIndices(AdjacentDomainIdx)%PTR,ERR,ERROR,*999)
                  CALL LIST_DATA_DIMENSION_SET(LocalGhostSendIndices(AdjacentDomainIdx)%PTR,1,ERR,ERROR,*999)
                  CALL LIST_DATA_TYPE_SET(LocalGhostSendIndices(AdjacentDomainIdx)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                  CALL LIST_INITIAL_SIZE_SET(LocalGhostSendIndices(AdjacentDomainIdx)%PTR,NumberBoundaryAndGhostNodes,&
                    & ERR,ERROR,*999)
                  CALL LIST_CREATE_FINISH(LocalGhostSendIndices(AdjacentDomainIdx)%PTR,ERR,ERROR,*999)
                  
                  ! create list for local ghost receive indices
                  NULLIFY(LocalGhostReceiveIndices(AdjacentDomainIdx)%PTR)
                  CALL LIST_CREATE_START(LocalGhostReceiveIndices(AdjacentDomainIdx)%PTR,ERR,ERROR,*999)
                  CALL LIST_DATA_DIMENSION_SET(LocalGhostReceiveIndices(AdjacentDomainIdx)%PTR,1,ERR,ERROR,*999)
                  CALL LIST_DATA_TYPE_SET(LocalGhostReceiveIndices(AdjacentDomainIdx)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                  CALL LIST_INITIAL_SIZE_SET(LocalGhostReceiveIndices(AdjacentDomainIdx)%PTR,NumberBoundaryAndGhostNodes,&
                    & ERR,ERROR,*999)
                  CALL LIST_CREATE_FINISH(LocalGhostReceiveIndices(AdjacentDomainIdx)%PTR,ERR,ERROR,*999)
                ENDDO
                
                ! InternalNodes contains the global numbers of the internal nodes
                
                ! loop over internal and boundary nodes and assign local node numbers
                InternalNodesIdx = 1
                BoundaryNodesIdx = 1
                LocalNodeNo = 1
                DomainListInternalIdx = NODES_MAPPING%INTERNAL_START
                DomainListBoundaryIdx = NODES_MAPPING%BOUNDARY_START
                DO WHILE(InternalNodesIdx <= NODES_MAPPING%NUMBER_OF_INTERNAL &
                  & .OR. BoundaryNodesIdx <= NODES_MAPPING%NUMBER_OF_BOUNDARY)
                  
                  ! get next internal or boundary node, ascending with global node number, if there are no more internal or boundary nodes, set to -1
                  IF (BoundaryNodesIdx <= NODES_MAPPING%NUMBER_OF_BOUNDARY) THEN
                    BoundaryAndGhostNodeIdx = BoundaryNodes(BoundaryNodesIdx)
                   BoundaryNodeGlobalNo = BoundaryAndGhostNodes(BoundaryAndGhostNodeIdx)
                  ELSE
                    BoundaryNodeGlobalNo = -1
                  ENDIF
                  
                  IF (InternalNodesIdx <= NODES_MAPPING%NUMBER_OF_INTERNAL) THEN
                    InternalNodeGlobalNo = InternalNodes(InternalNodesIdx)
                  ELSE
                    InternalNodeGlobalNo = -1
                  ENDIF
                  
                  IF ((InternalNodeGlobalNo < BoundaryNodeGlobalNo .AND. InternalNodeGlobalNo /= -1) &
                    & .OR. BoundaryNodeGlobalNo == -1) THEN
                    ! the next node is an internal node
                    NODES_MAPPING%DOMAIN_LIST(DomainListInternalIdx) = LocalNodeNo
                    NODES_MAPPING%LOCAL_TO_GLOBAL_MAP(LocalNodeNo) = InternalNodeGlobalNo
                    DomainListInternalIdx = DomainListInternalIdx + 1
                    InternalNodesIdx = InternalNodesIdx + 1
                    
                    IF (DIAGNOSTICS1) THEN

                      !PRINT *, MyComputationalNodeNumber, ": Node global ",InternalNodeGlobalNo," local ",LocalNodeNo,&
                      !  & ", interior" 
                    ENDIF
                  ELSE
                    ! the next node is a boundary node
                    NODES_MAPPING%DOMAIN_LIST(DomainListBoundaryIdx) = LocalNodeNo
                    NODES_MAPPING%LOCAL_TO_GLOBAL_MAP(LocalNodeNo) = BoundaryNodeGlobalNo
                    DomainListBoundaryIdx = DomainListBoundaryIdx + 1
                    BoundaryNodesIdx = BoundaryNodesIdx + 1
                    
                    IF (DIAGNOSTICS1) THEN

                      !PRINT *, MyComputationalNodeNumber, ": Node global ",BoundaryNodeGlobalNo," local ",LocalNodeNo,&
                      !  & ", boundary" 
                    ENDIF
                    
                    ! update LOCAL_GHOST_SEND_INDICES
                    AdjacentDomainIdx = 1
                    ! loop over the domains where this node is present
                    CALL LIST_NUMBER_OF_ITEMS_GET(DomainsOfNodeList(BoundaryAndGhostNodeIdx)%PTR, NumberDomains, ERR,ERROR,*999)
                    DO I=1,NumberDomains
                      ! get domain no
                      CALL LIST_ITEM_GET(DomainsOfNodeList(BoundaryAndGhostNodeIdx)%PTR,I,DomainNo,ERR,ERROR,*999)
                      IF (DomainNo /= MyComputationalNodeNumber) THEN
                            
                        ! go to entry in ADJACENT_DOMAINS array for DomainNo, at index AdjacentDomainIdx
                        DO
                          IF (AdjacentDomains(AdjacentDomainIdx) == DomainNo) THEN
                            
                            IF (DIAGNOSTICS1) THEN

                              !PRINT *, MyComputationalNodeNumber, ":    add local number to LOCAL_GHOST_SEND_INDICES of domain ", &
                              !  & DomainNo,", AdjacentDomainIdx=",AdjacentDomainIdx 
                            ENDIF
                            
                            CALL LIST_ITEM_ADD(LocalGhostSendIndices(AdjacentDomainIdx)%PTR,LocalNodeNo,ERR,ERROR,*999)
                            EXIT
                          ENDIF
                          AdjacentDomainIdx = AdjacentDomainIdx+1
                        ENDDO
                      
                      ENDIF
                    ENDDO
                    
                    DomainNo = BoundaryAndGhostNodesDomain(BoundaryAndGhostNodeIdx)
                    
                  ENDIF
                  LocalNodeNo = LocalNodeNo + 1
                ENDDO
                
                DomainListGhostIdx = NODES_MAPPING%GHOST_START
                ! loop over ghost nodes
                DO GhostNodeIdx = 1,NODES_MAPPING%NUMBER_OF_GHOST
                  GhostNodeGlobalNo = GhostNodesDomains(1,GhostNodeIdx)
                  GhostDomain = GhostNodesDomains(2,GhostNodeIdx)
                  
                  NODES_MAPPING%DOMAIN_LIST(DomainListGhostIdx) = LocalNodeNo
                  NODES_MAPPING%LOCAL_TO_GLOBAL_MAP(DomainListGhostIdx) = GhostNodeGlobalNo
                  DomainListGhostIdx = DomainListGhostIdx + 1
                  LocalNodeNo = LocalNodeNo + 1
                  
                  ! update LOCAL_GHOST_RECEIVE_INDICES
                  
                  IF (DIAGNOSTICS1) THEN

                    !PRINT *, MyComputationalNodeNumber, ": Node global ",GhostNodeGlobalNo," local ",LocalNodeNo-1,&
                    !  & ", ghost, is at domain ",GhostDomain 
                  ENDIF
                  
                  ! find entry ADJACENT_DOMAINS array for domainNo
                  AdacentDomainEntryFound = .FALSE. 
                  DO AdjacentDomainIdx=1,NODES_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
                    IF (NODES_MAPPING%ADJACENT_DOMAINS(AdjacentDomainIdx)%DOMAIN_NUMBER == GhostDomain) THEN
                      
                      IF (DIAGNOSTICS1) THEN

                        !PRINT *, MyComputationalNodeNumber, ":    add local number to NUMBER_OF_RECEIVE_GHOSTS of domain ",&
                        !  & GhostDomain,", AdjacentDomainIdx=",AdjacentDomainIdx 
                      ENDIF

                      CALL LIST_ITEM_ADD(LocalGhostReceiveIndices(AdjacentDomainIdx)%PTR,LocalNodeNo-1,ERR,ERROR,*999)
                      AdacentDomainEntryFound = .TRUE.                             
                      EXIT
                    ENDIF

                  ENDDO ! AdjacentDomainIdx
                  
                  ! TODO: fix following
                  IF (.NOT. AdacentDomainEntryFound.AND..TRUE.) THEN
                    IF (DIAGNOSTICS1) THEN
                       !PRINT *, MyComputationalNodeNumber, ": create new AdjacentDomainEntry" 
                    ENDIF
                  
                    ! create new AdjacentDomain entry in NODES_MAPPING%ADJACENT_DOMAINS
                    NULLIFY(NewAdjacentDomain)
                    ALLOCATE(NewAdjacentDomain,STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate NewAdjacentDomain.",ERR,ERROR,*999)
                    
                    NewAdjacentDomain%DOMAIN_NUMBER = GhostDomain
                    NewAdjacentDomain%NUMBER_OF_SEND_GHOSTS = 0
                    NewAdjacentDomain%NUMBER_OF_RECEIVE_GHOSTS = 0
                    NewAdjacentDomain%NUMBER_OF_FURTHER_LINKED_GHOSTS = 0
                    
                    NODES_MAPPING%NUMBER_OF_ADJACENT_DOMAINS = NODES_MAPPING%NUMBER_OF_ADJACENT_DOMAINS+1
                    
                    
                    IF (DIAGNOSTICS1) THEN

                      !PRINT *, MyComputationalNodeNumber, ": NODES_MAPPING%NUMBER_OF_ADJACENT_DOMAINS=",&
                      !  &NODES_MAPPING%NUMBER_OF_ADJACENT_DOMAINS 
                    ENDIF                    
                    
                    IF (NODES_MAPPING%NUMBER_OF_ADJACENT_DOMAINS == 1) THEN
                      ALLOCATE(NODES_MAPPING%ADJACENT_DOMAINS(1),STAT=ERR)
                      IF(ERR/=0) CALL FlagError("Could not allocate NODES_MAPPING%ADJACENT_DOMAINS.",ERR,ERROR,*999)
                      
                      NODES_MAPPING%ADJACENT_DOMAINS(1) = NewAdjacentDomain
                    ELSE
                      NODES_MAPPING%ADJACENT_DOMAINS = RESHAPE(NODES_MAPPING%ADJACENT_DOMAINS, &
                        & [NODES_MAPPING%NUMBER_OF_ADJACENT_DOMAINS], [NewAdjacentDomain])
                    ENDIF
                    
                    IF (DIAGNOSTICS1) THEN
                       !PRINT *, MyComputationalNodeNumber, ": create lists" 
                    ENDIF
                    
                    
                    ! create list for local ghost send indices
                    ALLOCATE(NewLocalGhostSendIndices,STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate NewLocalGhostSendIndices.",ERR,ERROR,*999)
                    
                    ALLOCATE(NewLocalGhostReceiveIndices,STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate NewLocalGhostReceiveIndices.",ERR,ERROR,*999)
                    
                    NULLIFY(NewLocalGhostSendIndices%PTR)
                    CALL LIST_CREATE_START(NewLocalGhostSendIndices%PTR,ERR,ERROR,*999)
                    CALL LIST_DATA_DIMENSION_SET(NewLocalGhostSendIndices%PTR,1,ERR,ERROR,*999)
                    CALL LIST_DATA_TYPE_SET(NewLocalGhostSendIndices%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                    CALL LIST_CREATE_FINISH(NewLocalGhostSendIndices%PTR,ERR,ERROR,*999)
                    
                    ! create list for local ghost receive indices
                    NULLIFY(NewLocalGhostReceiveIndices%PTR)
                    CALL LIST_CREATE_START(NewLocalGhostReceiveIndices%PTR,ERR,ERROR,*999)
                    CALL LIST_DATA_DIMENSION_SET(NewLocalGhostReceiveIndices%PTR,1,ERR,ERROR,*999)
                    CALL LIST_DATA_TYPE_SET(NewLocalGhostReceiveIndices%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                    CALL LIST_CREATE_FINISH(NewLocalGhostReceiveIndices%PTR,ERR,ERROR,*999)
                  
                    
                    IF (DIAGNOSTICS1) THEN
                       !PRINT *, MyComputationalNodeNumber, ": reshape" 
                    ENDIF
                    IF (NODES_MAPPING%NUMBER_OF_ADJACENT_DOMAINS == 1) THEN
                      ALLOCATE(LocalGhostReceiveIndices(1),STAT=ERR)
                      IF(ERR/=0) CALL FlagError("Could not allocate LocalGhostReceiveIndices.",ERR,ERROR,*999)
                      
                      LocalGhostReceiveIndices(1) = NewLocalGhostReceiveIndices
                    ELSE
                      LocalGhostReceiveIndices = RESHAPE(LocalGhostReceiveIndices, [NODES_MAPPING%NUMBER_OF_ADJACENT_DOMAINS], &
                        & [NewLocalGhostReceiveIndices])
                    ENDIF
                    
                    IF (NODES_MAPPING%NUMBER_OF_ADJACENT_DOMAINS == 1) THEN
                      ALLOCATE(LocalGhostSendIndices(1),STAT=ERR)
                      IF(ERR/=0) CALL FlagError("Could not allocate LocalGhostSendIndices.",ERR,ERROR,*999)
                      
                      LocalGhostSendIndices(1) = NewLocalGhostSendIndices
                    ELSE
                      LocalGhostSendIndices = RESHAPE(LocalGhostSendIndices, [NODES_MAPPING%NUMBER_OF_ADJACENT_DOMAINS], &
                        & [NewLocalGhostSendIndices])
                    ENDIF
                    IF (NODES_MAPPING%NUMBER_OF_ADJACENT_DOMAINS == 1) THEN
                      ALLOCATE(AdjacentDomains(1),STAT=ERR)
                      IF(ERR/=0) CALL FlagError("Could not allocate AdjacentDomains.",ERR,ERROR,*999)
                      
                      AdjacentDomains(1) = GhostDomain
                    ELSE
                      AdjacentDomains = RESHAPE(AdjacentDomains, [NODES_MAPPING%NUMBER_OF_ADJACENT_DOMAINS], &
                        & [GhostDomain])
                    ENDIF
                    
                    IF (DIAGNOSTICS1) THEN
                       !PRINT *, MyComputationalNodeNumber, ": add to list" 
                    ENDIF
                    AdjacentDomainIdx = NODES_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
                    NumberAdjacentDomains = NODES_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
                    IF (DIAGNOSTICS1) THEN

                      !PRINT *, MyComputationalNodeNumber, ":    add local number to NUMBER_OF_RECEIVE_GHOSTS of domain ",&
                      !  & GhostDomain,", AdjacentDomainIdx=",AdjacentDomainIdx 
                    ENDIF
                          
                    CALL LIST_ITEM_ADD(LocalGhostReceiveIndices(AdjacentDomainIdx)%PTR,LocalNodeNo-1,ERR,ERROR,*999)
                  ENDIF
                  
                ENDDO  ! GhostNodeIdx 
                
                DO AdjacentDomainIdx=1,NumberAdjacentDomains
                  ! create LocalGhostSendIndices for the adjacent domain
                  CALL LIST_DETACH_AND_DESTROY(LocalGhostSendIndices(AdjacentDomainIdx)%PTR, &
                    & NODES_MAPPING%ADJACENT_DOMAINS(AdjacentDomainIdx)%NUMBER_OF_SEND_GHOSTS,IntegerArray,ERR,ERROR,*999)
                  NODES_MAPPING%ADJACENT_DOMAINS(AdjacentDomainIdx)%LOCAL_GHOST_SEND_INDICES = &
                    & IntegerArray(1:NODES_MAPPING%ADJACENT_DOMAINS(AdjacentDomainIdx)%NUMBER_OF_SEND_GHOSTS)
                  DEALLOCATE(IntegerArray)
                
                  ! create LocalGhostReceiveIndices for the adjacent domain
                  CALL LIST_DETACH_AND_DESTROY(LocalGhostReceiveIndices(AdjacentDomainIdx)%PTR, &
                    & NODES_MAPPING%ADJACENT_DOMAINS(AdjacentDomainIdx)%NUMBER_OF_RECEIVE_GHOSTS,IntegerArray,ERR,ERROR,*999)
                  NODES_MAPPING%ADJACENT_DOMAINS(AdjacentDomainIdx)%LOCAL_GHOST_RECEIVE_INDICES = &
                    & IntegerArray(1:NODES_MAPPING%ADJACENT_DOMAINS(AdjacentDomainIdx)%NUMBER_OF_RECEIVE_GHOSTS)
                  DEALLOCATE(IntegerArray)
                      
                  IF (DIAGNOSTICS1) THEN
                    AdjacentDomain = AdjacentDomains(AdjacentDomainIdx)

                    !PRINT *, MyComputationalNodeNumber, ": adjacent domain ",AdjacentDomain,", LOCAL_GHOST_SEND_INDICES:", &
                    !  & NODES_MAPPING%ADJACENT_DOMAINS(AdjacentDomainIdx)%LOCAL_GHOST_SEND_INDICES(1: &
                    !  & NODES_MAPPING%ADJACENT_DOMAINS(AdjacentDomainIdx)%NUMBER_OF_SEND_GHOSTS), &
                    !  & ", LOCAL_GHOST_RECEIVE_INDICES: ",&
                    !  & NODES_MAPPING%ADJACENT_DOMAINS(AdjacentDomainIdx)%LOCAL_GHOST_RECEIVE_INDICES(1: &
                    !  & NODES_MAPPING%ADJACENT_DOMAINS(AdjacentDomainIdx)%NUMBER_OF_RECEIVE_GHOSTS) 
                  ENDIF
                ENDDO
                
                
                ! exchange NUMBER_OF_DOMAIN_LOCAL and NUMBER_OF_DOMAIN_GHOST
                ! allocate number_of_domain_local
                ALLOCATE(NODES_MAPPING%NUMBER_OF_DOMAIN_LOCAL(0:NumberComputationalNodes-1),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate NUMBER_OF_DOMAIN_LOCAL.",ERR,ERROR,*999)

                NODES_MAPPING%NUMBER_OF_DOMAIN_LOCAL(MyComputationalNodeNumber) = NODES_MAPPING%NUMBER_OF_LOCAL
                      
                CALL MPI_ALLGATHER(MPI_IN_PLACE,1,MPI_INTEGER,NODES_MAPPING%&
                  & NUMBER_OF_DOMAIN_LOCAL(0:NumberComputationalNodes-1), &
                  & 1,MPI_INTEGER,groupCommunicator,MPI_IERROR)
                CALL MPI_ERROR_CHECK("MPI_ALLGATHER",MPI_IERROR,ERR,ERROR,*999)
                      
                ! allocate number_of_domain_ghost
                ALLOCATE(NODES_MAPPING%NUMBER_OF_DOMAIN_GHOST(0:NumberComputationalNodes-1),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate NUMBER_OF_DOMAIN_GHOST.",ERR,ERROR,*999)

                NODES_MAPPING%NUMBER_OF_DOMAIN_GHOST(MyComputationalNodeNumber) = NODES_MAPPING%NUMBER_OF_GHOST
                      
                CALL MPI_ALLGATHER(MPI_IN_PLACE,1,MPI_INTEGER,NODES_MAPPING%&
                  & NUMBER_OF_DOMAIN_GHOST(0:NumberComputationalNodes-1), &
                  & 1,MPI_INTEGER,groupCommunicator,MPI_IERROR)
                CALL MPI_ERROR_CHECK("MPI_ALLGATHER",MPI_IERROR,ERR,ERROR,*999)
                      
              
                IF (DIAGNOSTICS1) THEN

                  !PRINT *, MyComputationalNodeNumber, ": NUMBER_OF_DOMAIN_LOCAL:"
                  !DO I=0,NumberComputationalNodes-1
                  !  PRINT *, MyComputationalNodeNumber, ": rank ",I,": ",NODES_MAPPING%NUMBER_OF_DOMAIN_LOCAL(I)
                  !ENDDO
                  !PRINT *, MyComputationalNodeNumber, ": NUMBER_OF_DOMAIN_GHOST:"
                  !DO I=0,NumberComputationalNodes-1
                  !  PRINT *, MyComputationalNodeNumber, ": rank ",I,": ",NODES_MAPPING%NUMBER_OF_DOMAIN_GHOST(I)
                  !ENDDO 
                ENDIF
                
                ! ----------------------------------------------
                ! set some values in DOFS_MAPPING
                
                IF(ASSOCIATED(DOMAIN%MAPPINGS%DOFS)) THEN
                  DOFS_MAPPING=>DOMAIN%MAPPINGS%DOFS
          
                  ! count number of local dofs
                  ! loop over internal nodes
                  DO NodeIdx = NODES_MAPPING%INTERNAL_START,NODES_MAPPING%INTERNAL_FINISH
                    NodeLocalNo = NODES_MAPPING%DOMAIN_LIST(NodeIdx)
                    NodeGlobalNo = NODES_MAPPING%LOCAL_TO_GLOBAL_MAP(NodeLocalNo)
                    
                    DO DerivativeIdx=1,TOPOLOGY%NODES%NODES(NodeGlobalNo)%numberOfDerivatives
                      DOFS_MAPPING%NUMBER_OF_INTERNAL = DOFS_MAPPING%NUMBER_OF_INTERNAL + &
                        & TOPOLOGY%NODES%NODES(NodeGlobalNo)%DERIVATIVES(DerivativeIdx)%numberOfVersions
                    ENDDO
                  ENDDO
                  
                  ! count number of boundary dofs
                  ! loop over boundary nodes
                  DO NodeIdx = NODES_MAPPING%BOUNDARY_START,NODES_MAPPING%BOUNDARY_FINISH
                    NodeLocalNo = NODES_MAPPING%DOMAIN_LIST(NodeIdx)
                    NodeGlobalNo = NODES_MAPPING%LOCAL_TO_GLOBAL_MAP(NodeLocalNo)
                    
                    DO DerivativeIdx=1,TOPOLOGY%NODES%NODES(NodeGlobalNo)%numberOfDerivatives
                      DOFS_MAPPING%NUMBER_OF_BOUNDARY = DOFS_MAPPING%NUMBER_OF_BOUNDARY + &
                        & TOPOLOGY%NODES%NODES(NodeGlobalNo)%DERIVATIVES(DerivativeIdx)%numberOfVersions
                    ENDDO
                  ENDDO
                  
                  ! count number of ghost dofs
                  ! loop over ghost nodes
                  DO NodeIdx = NODES_MAPPING%GHOST_START,NODES_MAPPING%GHOST_FINISH
                    NodeLocalNo = NODES_MAPPING%DOMAIN_LIST(NodeIdx)
                    NodeGlobalNo = NODES_MAPPING%LOCAL_TO_GLOBAL_MAP(NodeLocalNo)
                    
                    DO DerivativeIdx=1,TOPOLOGY%NODES%NODES(NodeGlobalNo)%numberOfDerivatives
                      DOFS_MAPPING%NUMBER_OF_GHOST = DOFS_MAPPING%NUMBER_OF_GHOST + &
                        & TOPOLOGY%NODES%NODES(NodeGlobalNo)%DERIVATIVES(DerivativeIdx)%numberOfVersions
                    ENDDO
                  ENDDO
      
                  
                  ! sum up number of local and total number of local
                  DOFS_MAPPING%NUMBER_OF_LOCAL = DOFS_MAPPING%NUMBER_OF_INTERNAL + DOFS_MAPPING%NUMBER_OF_BOUNDARY
                  DOFS_MAPPING%TOTAL_NUMBER_OF_LOCAL = DOFS_MAPPING%NUMBER_OF_LOCAL + DOFS_MAPPING%NUMBER_OF_GHOST
               
                  ! allreduce number of dofs
                  CALL MPI_ALLREDUCE(DOFS_MAPPING%NUMBER_OF_LOCAL,DOFS_MAPPING%NUMBER_OF_GLOBAL,1,MPI_INT,MPI_SUM, &
                    & groupCommunicator,MPI_IERROR)
                  CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
                  
                  ! compute index positions for DOMAIN_LIST
                  DOFS_MAPPING%INTERNAL_START=1
                  DOFS_MAPPING%INTERNAL_FINISH=DOFS_MAPPING%INTERNAL_START+DOFS_MAPPING%NUMBER_OF_INTERNAL-1
                  DOFS_MAPPING%BOUNDARY_START=DOFS_MAPPING%INTERNAL_FINISH+1
                  DOFS_MAPPING%BOUNDARY_FINISH=DOFS_MAPPING%BOUNDARY_START+DOFS_MAPPING%NUMBER_OF_BOUNDARY-1
                  DOFS_MAPPING%GHOST_START=DOFS_MAPPING%BOUNDARY_FINISH+1
                  DOFS_MAPPING%GHOST_FINISH=DOFS_MAPPING%GHOST_START+DOFS_MAPPING%NUMBER_OF_GHOST-1
            
                  ! exchange NUMBER_OF_DOMAIN_LOCAL and NUMBER_OF_DOMAIN_GHOST
                  ! allocate number_of_domain_local
                  ALLOCATE(DOFS_MAPPING%NUMBER_OF_DOMAIN_LOCAL(0:NumberComputationalNodes-1),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate NUMBER_OF_DOMAIN_LOCAL.",ERR,ERROR,*999)

                  DOFS_MAPPING%NUMBER_OF_DOMAIN_LOCAL(MyComputationalNodeNumber) = DOFS_MAPPING%NUMBER_OF_LOCAL
                 
                  CALL MPI_ALLGATHER(MPI_IN_PLACE,1,MPI_INTEGER,DOFS_MAPPING%&
                    & NUMBER_OF_DOMAIN_LOCAL(0:NumberComputationalNodes-1), &
                    & 1,MPI_INTEGER,groupCommunicator,MPI_IERROR)
                  CALL MPI_ERROR_CHECK("MPI_ALLGATHER",MPI_IERROR,ERR,ERROR,*999)
                         
                  ! allocate number_of_domain_ghost
                  ALLOCATE(DOFS_MAPPING%NUMBER_OF_DOMAIN_GHOST(0:NumberComputationalNodes-1),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate NUMBER_OF_DOMAIN_GHOST.",ERR,ERROR,*999)

                  DOFS_MAPPING%NUMBER_OF_DOMAIN_GHOST(MyComputationalNodeNumber) = DOFS_MAPPING%NUMBER_OF_GHOST
                        
                  CALL MPI_ALLGATHER(MPI_IN_PLACE,1,MPI_INTEGER,DOFS_MAPPING%&
                    & NUMBER_OF_DOMAIN_GHOST(0:NumberComputationalNodes-1), &
                    & 1,MPI_INTEGER,groupCommunicator,MPI_IERROR)
                  CALL MPI_ERROR_CHECK("MPI_ALLGATHER",MPI_IERROR,ERR,ERROR,*999)
                        
                  ! adjacent domains from nodes mapping
                  DOFS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS = NODES_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
                  
                  ALLOCATE(DOFS_MAPPING%ADJACENT_DOMAINS(1:DOFS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate ADJACENT_DOMAINS.",ERR,ERROR,*999)

                  DOFS_MAPPING%ADJACENT_DOMAINS = NODES_MAPPING%ADJACENT_DOMAINS
                  DOFS_MAPPING%NUMBER_OF_DOMAINS = NODES_MAPPING%NUMBER_OF_DOMAINS
       
                ENDIF 

                ! deallocate arrays
                DEALLOCATE(BoundaryAndGhostNodes)
                DEALLOCATE(AdjacentDomains)
                DEALLOCATE(SendRequestHandle)
                DEALLOCATE(ReceiveRequestHandle)
                DEALLOCATE(NumberInteriorAndLocalNodesOnRank)
                DEALLOCATE(LocalAndAdjacentDomains)
                DEALLOCATE(BoundaryAndGhostNodesDomain)
                DEALLOCATE(NumberNodesInDomain)
                DEALLOCATE(GhostNodes)
                DEALLOCATE(InternalNodes)
                DEALLOCATE(BoundaryNodes)
                DEALLOCATE(SendRequestHandle0)
                DEALLOCATE(SendRequestHandle1)
                DEALLOCATE(GhostNodesDomains)
                DEALLOCATE(SendBuffer)
                DEALLOCATE(NumberToSendToDomain)
                      
                !PRINT *, "Stop in DOMAIN_MAPPINGS_NODES_CALCULATE"
                !STOP
                
              ELSE
                CALL FlagError("Domain mesh is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Elements mapping is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Domain mappings nodes is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Domain mappings is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_MAPPINGS_NODES_CALCULATE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_NODES_CALCULATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_NODES_CALCULATE
  
  
  !
  !================================================================================================================================
  !

  !>Calculates the local/global node and dof mappings for a domain decomposition.
  SUBROUTINE DOMAIN_MAPPINGS_NODES_DOFS_CALCULATE(DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<A pointer to the domain to calculate the node dofs for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,no_adjacent_element,no_computational_node,no_ghost_node,adjacent_element,ghost_node, &
      & NUMBER_OF_NODES_PER_DOMAIN,domain_idx,domain_idx2,domain_no,node_idx,derivative_idx,version_idx,ny,NUMBER_OF_DOMAINS, &
      & MAX_NUMBER_DOMAINS,NUMBER_OF_GHOST_NODES,myComputationalNodeNumber,numberOfComputationalNodes,component_idx
    INTEGER(INTG), ALLOCATABLE :: LOCAL_NODE_NUMBERS(:),LOCAL_DOF_NUMBERS(:),NODE_COUNT(:),NUMBER_INTERNAL_NODES(:), &
      & NUMBER_BOUNDARY_NODES(:)
    INTEGER(INTG), ALLOCATABLE :: DOMAINS(:),ALL_DOMAINS(:),GHOST_NODES(:)
    LOGICAL :: BOUNDARY_DOMAIN
    TYPE(LIST_TYPE), POINTER :: ADJACENT_DOMAINS_LIST,ALL_ADJACENT_DOMAINS_LIST
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: GHOST_NODES_LIST(:)
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(MeshComponentTopologyType), POINTER :: MESH_TOPOLOGY
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: NODES_MAPPING
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOFS_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    ENTERS("DOMAIN_MAPPINGS_NODES_DOFS_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      IF(ASSOCIATED(DOMAIN%MAPPINGS)) THEN
        IF(ASSOCIATED(DOMAIN%MAPPINGS%NODES)) THEN
          NODES_MAPPING=>DOMAIN%MAPPINGS%NODES
          IF(ASSOCIATED(DOMAIN%MAPPINGS%DOFS)) THEN
            DOFS_MAPPING=>DOMAIN%MAPPINGS%DOFS
            IF(ASSOCIATED(DOMAIN%MAPPINGS%ELEMENTS)) THEN
              ELEMENTS_MAPPING=>DOMAIN%MAPPINGS%ELEMENTS
              IF(ASSOCIATED(DOMAIN%DECOMPOSITION)) THEN
                DECOMPOSITION=>DOMAIN%DECOMPOSITION
                IF(ASSOCIATED(DOMAIN%MESH)) THEN
                  MESH=>DOMAIN%MESH
                  component_idx=DOMAIN%MESH_COMPONENT_NUMBER
                  MESH_TOPOLOGY=>MESH%TOPOLOGY(component_idx)%PTR
                  
                  CALL WorkGroup_NumberOfGroupNodesGet(NODES_MAPPING%workGroup,numberOfComputationalNodes, &
                    & err,error,*999)
                  CALL WorkGroup_GroupNodeNumberGet(NODES_MAPPING%workGroup,myComputationalNodeNumber, &
                    & err,error,*999)
                  !numberOfComputationalNodes=ComputationalEnvironment_NumberOfNodesGet(ERR,ERROR)
                  !IF(ERR/=0) GOTO 999
                  !myComputationalNodeNumber=ComputationalEnvironment_NodeNumberGet(ERR,ERROR)
                  !IF(ERR/=0) GOTO 999
                  
                  !Calculate the local and global numbers and set up the mappings
                  ! allocate arrays 
                  ALLOCATE(NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(MESH_TOPOLOGY%NODES%numberOfNodes),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate node mapping global to local map.",ERR,ERROR,*999)
                  NODES_MAPPING%NUMBER_OF_GLOBAL=MESH_TOPOLOGY%NODES%numberOfNodes
                  ALLOCATE(LOCAL_NODE_NUMBERS(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate local node numbers.",ERR,ERROR,*999)
                  LOCAL_NODE_NUMBERS=0
                  ALLOCATE(DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(MESH_TOPOLOGY%dofs%numberOfDofs),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate dofs mapping global to local map.",ERR,ERROR,*999)
                  DOFS_MAPPING%NUMBER_OF_GLOBAL=MESH_TOPOLOGY%DOFS%numberOfDofs
                  ALLOCATE(LOCAL_DOF_NUMBERS(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate local dof numbers.",ERR,ERROR,*999)
                  LOCAL_DOF_NUMBERS=0
                  ALLOCATE(GHOST_NODES_LIST(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate ghost nodes list.",ERR,ERROR,*999)
                  DO domain_idx=0,DECOMPOSITION%NUMBER_OF_DOMAINS-1
                    NULLIFY(GHOST_NODES_LIST(domain_idx)%PTR)
                    CALL LIST_CREATE_START(GHOST_NODES_LIST(domain_idx)%PTR,ERR,ERROR,*999)
                    CALL LIST_DATA_TYPE_SET(GHOST_NODES_LIST(domain_idx)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                    CALL LIST_INITIAL_SIZE_SET(GHOST_NODES_LIST(domain_idx)%PTR,INT(MESH_TOPOLOGY%NODES%numberOfNodes/2), &
                      & ERR,ERROR,*999)
                    CALL LIST_CREATE_FINISH(GHOST_NODES_LIST(domain_idx)%PTR,ERR,ERROR,*999)
                  ENDDO !domain_idx
                  ALLOCATE(NUMBER_INTERNAL_NODES(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate number of internal nodes.",ERR,ERROR,*999)
                  NUMBER_INTERNAL_NODES=0
                  ALLOCATE(NUMBER_BOUNDARY_NODES(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate number of boundary nodes.",ERR,ERROR,*999)
                  NUMBER_BOUNDARY_NODES=0

                  IF (DIAGNOSTICS1) THEN
                    CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  n nodes: ",MESH_TOPOLOGY%NODES%numberOfNodes,ERR,ERROR,*999)
                  ENDIF 

                  !For the first pass just determine the internal and boundary nodes (to allocate space in lists)
                  ! loop over all nodes
                  DO node_idx = 1,MESH_TOPOLOGY%NODES%numberOfNodes
                  
                  ! ADJACENT_DOMAINS_LIST will contain the domains where an adjacent elements to this node is an internal or boundary element 
                    NULLIFY(ADJACENT_DOMAINS_LIST)
                    CALL LIST_CREATE_START(ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                    CALL LIST_DATA_TYPE_SET(ADJACENT_DOMAINS_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                    CALL LIST_INITIAL_SIZE_SET(ADJACENT_DOMAINS_LIST,DECOMPOSITION%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
                    CALL LIST_CREATE_FINISH(ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                    ! ALL_ADJACENT_DOMAINS_LIST will contain the domains where an adjacent elements to this node is stored, even as ghost 
                    NULLIFY(ALL_ADJACENT_DOMAINS_LIST)
                    CALL LIST_CREATE_START(ALL_ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                    CALL LIST_DATA_TYPE_SET(ALL_ADJACENT_DOMAINS_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                    CALL LIST_INITIAL_SIZE_SET(ALL_ADJACENT_DOMAINS_LIST,DECOMPOSITION%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
                    CALL LIST_CREATE_FINISH(ALL_ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)

                    ! loop over adjacent elements to that node
                    DO no_adjacent_element = 1,MESH_TOPOLOGY%NODES%NODES(node_idx)%numberOfSurroundingElements
                      adjacent_element = MESH_TOPOLOGY%NODES%NODES(node_idx)%surroundingElements(no_adjacent_element)
                      domain_no = DECOMPOSITION%ELEMENT_DOMAIN(adjacent_element) 
                      CALL LIST_ITEM_ADD(ADJACENT_DOMAINS_LIST,domain_no,ERR,ERROR,*999)
                      ! loop over all occurences of the element on different ranks 
                      DO domain_idx=1,ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%NUMBER_OF_DOMAINS
                        CALL LIST_ITEM_ADD(ALL_ADJACENT_DOMAINS_LIST,ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)% &
                          & DOMAIN_NUMBER(domain_idx),ERR,ERROR,*999)
                      ENDDO !domain_idx
                    ENDDO !no_adjacent_element
                    CALL LIST_REMOVE_DUPLICATES(ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                    CALL LIST_DETACH_AND_DESTROY(ADJACENT_DOMAINS_LIST,NUMBER_OF_DOMAINS,DOMAINS,ERR,ERROR,*999)
                    CALL LIST_REMOVE_DUPLICATES(ALL_ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                    CALL LIST_DETACH_AND_DESTROY(ALL_ADJACENT_DOMAINS_LIST,MAX_NUMBER_DOMAINS,ALL_DOMAINS,ERR,ERROR,*999)
! if an adjacent element is a ghost element on the local rank
                    IF(NUMBER_OF_DOMAINS /= MAX_NUMBER_DOMAINS) THEN !Ghost node
                      
                      ! loop over all adjacent domains at the nodes, those that have a ghost node included
                      DO domain_idx = 1,MAX_NUMBER_DOMAINS
                        domain_no=ALL_DOMAINS(domain_idx)
                        ! determine if the current domain contains a boundary element, i.e. not only ghost elements 
                        BOUNDARY_DOMAIN=.FALSE.

                        DO domain_idx2 = 1,NUMBER_OF_DOMAINS
                          IF(domain_no == DOMAINS(domain_idx2)) THEN 
                            BOUNDARY_DOMAIN=.TRUE.
                            EXIT
                          ENDIF
                        ENDDO !domain_idx2
                        ! if the domain domain_no only contains ghost elements adjacent to this node, add node to ghost nodes list of that domain 
                        IF(.NOT.BOUNDARY_DOMAIN) CALL LIST_ITEM_ADD(GHOST_NODES_LIST(domain_no)%PTR,node_idx,ERR,ERROR,*999)
                      ENDDO !domain_idx
                    ENDIF

                    ! allocate arrays 
                    ALLOCATE(NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_NUMBER(MAX_NUMBER_DOMAINS),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate node global to local map local number.",ERR,ERROR,*999)
                    ALLOCATE(NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%DOMAIN_NUMBER(MAX_NUMBER_DOMAINS),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate node global to local map domain number.",ERR,ERROR,*999)
                    ALLOCATE(NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_TYPE(MAX_NUMBER_DOMAINS),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate node global to local map local type.",ERR,ERROR,*999)
                    DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%numberOfDerivatives
                      DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%numberOfVersions
                        ny = MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%dofIndex(version_idx) 
                        ALLOCATE(DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(MAX_NUMBER_DOMAINS),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate dof global to local map local number.",ERR,ERROR,*999)
                        ALLOCATE(DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER(MAX_NUMBER_DOMAINS),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate dof global to local map domain number.",ERR,ERROR,*999)
                        ALLOCATE(DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE(MAX_NUMBER_DOMAINS),STAT=ERR)
                        IF(ERR/=0) CALL FlagError("Could not allocate dof global to local map local type.",ERR,ERROR,*999)
                      ENDDO !version_idx
                    ENDDO !derivative_idx
                    ! if only one domain is adjacent to the node 
                    IF(NUMBER_OF_DOMAINS == 1) THEN
                      !Node is an internal node
                      ! increase counter of internal nodes on its domain 
                      domain_no=DOMAINS(1)
                      NUMBER_INTERNAL_NODES(domain_no)=NUMBER_INTERNAL_NODES(domain_no)+1
                      !LOCAL_NODE_NUMBERS(domain_no)=LOCAL_NODE_NUMBERS(domain_no)+1
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%NUMBER_OF_DOMAINS=1
                      !NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_NUMBER(1)=LOCAL_NODE_NUMBERS(DOMAINS(1))
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_NUMBER(1)=-1
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%DOMAIN_NUMBER(1)=DOMAINS(1) 
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_TYPE(1)=DOMAIN_LOCAL_INTERNAL
                      ! fill DOFS_MAPPING 
                      DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%numberOfDerivatives
                        DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%numberOfVersions
                          ny=MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%dofIndex(version_idx)
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS=1
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(1)=-1
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER(1)=domain_no
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE(1)=DOMAIN_LOCAL_INTERNAL
                        ENDDO !version_idx
                      ENDDO !derivative_idx
                    ELSE
                      !Node is on the boundary of computational domains
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%NUMBER_OF_DOMAINS=NUMBER_OF_DOMAINS
                      DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%numberOfDerivatives
                        DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%numberOfVersions
                          ny=MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%dofIndex(version_idx)
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS=NUMBER_OF_DOMAINS
                        ENDDO !version_idx
                      ENDDO !derivative_idx
                      DO domain_idx=1,NUMBER_OF_DOMAINS
                        domain_no=DOMAINS(domain_idx)
                        !LOCAL_NODE_NUMBERS(domain_no)=LOCAL_NODE_NUMBERS(domain_no)+1
                        !NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_NUMBER(domain_idx)=LOCAL_NODE_NUMBERS(domain_no)
                        NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_NUMBER(domain_idx)=-1
                        NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%DOMAIN_NUMBER(domain_idx)=domain_no
                        NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_TYPE(domain_idx)=DOMAIN_LOCAL_BOUNDARY
                        DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%numberOfDerivatives
                          DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%numberOfVersions
                            ny=MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%dofIndex(version_idx)
                            DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(domain_idx)=-1
                            DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER(domain_idx)=domain_no
                            DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE(domain_idx)=DOMAIN_LOCAL_BOUNDARY
                          ENDDO !version_idx
                        ENDDO !derivative_idx
                      ENDDO !domain_idx
                    ENDIF
                    DEALLOCATE(DOMAINS) 
                    DEALLOCATE(ALL_DOMAINS)
                  ENDDO !node_idx

                  !For the second pass assign boundary nodes to one domain on the boundary and set local node numbers.
                  ! compute average number of nodes per domain 
                  NUMBER_OF_NODES_PER_DOMAIN=FLOOR(REAL(MESH_TOPOLOGY%NODES%numberOfNodes,DP)/ &
                    & REAL(DECOMPOSITION%NUMBER_OF_DOMAINS,DP))
                  ALLOCATE(DOMAIN%NODE_DOMAIN(MESH_TOPOLOGY%NODES%numberOfNodes),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate node domain",ERR,ERROR,*999)
                  DOMAIN%NODE_DOMAIN=-1
                  ! loop over nodes 
                  DO node_idx=1,MESH_TOPOLOGY%NODES%numberOfNodes
                    ! internal node 
                    IF(NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%NUMBER_OF_DOMAINS==1) THEN !Internal node
                      domain_no=NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%DOMAIN_NUMBER(1)
                      DOMAIN%NODE_DOMAIN(node_idx)=domain_no
                      ! set local node number 
                      LOCAL_NODE_NUMBERS(domain_no)=LOCAL_NODE_NUMBERS(domain_no)+1
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_NUMBER(1)=LOCAL_NODE_NUMBERS(domain_no)
                      ! DOFS_MAPPING 
                      DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%numberOfDerivatives
                        DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%numberOfVersions
                          ny=MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%dofIndex(version_idx)
                          LOCAL_DOF_NUMBERS(domain_no)=LOCAL_DOF_NUMBERS(domain_no)+1
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(1)=LOCAL_DOF_NUMBERS(domain_no)
                        ENDDO !version_idx
                      ENDDO !derivative_idx
                    ELSE !Boundary node
                      NUMBER_OF_DOMAINS=NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%NUMBER_OF_DOMAINS
                      ! loop over domains that have an adjacent not-ghost element to this node 
                      DO domain_idx=1,NUMBER_OF_DOMAINS
                        domain_no=NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%DOMAIN_NUMBER(domain_idx)
                        ! if node is not yet assigned to a domain 
                        IF(DOMAIN%NODE_DOMAIN(node_idx)<0) THEN
                          IF((NUMBER_INTERNAL_NODES(domain_no) + NUMBER_BOUNDARY_NODES(domain_no) < NUMBER_OF_NODES_PER_DOMAIN).OR.& 
                            & (domain_idx==NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%NUMBER_OF_DOMAINS)) THEN
                            !Allocate the node to this domain
                            DOMAIN%NODE_DOMAIN(node_idx)=domain_no
                            ! increase counter of boundary nodes for the domain 
                            NUMBER_BOUNDARY_NODES(domain_no)=NUMBER_BOUNDARY_NODES(domain_no)+1
                            ! assign local node number 
                            LOCAL_NODE_NUMBERS(domain_no)=LOCAL_NODE_NUMBERS(domain_no)+1
                            !Reset the boundary information to be in the first domain index. The remaining domain indicies will
                            !be overwritten when the ghost nodes are calculated below. 
                            NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%NUMBER_OF_DOMAINS=1
                            NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_NUMBER(1)=LOCAL_NODE_NUMBERS(domain_no)
                            NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%DOMAIN_NUMBER(1)=domain_no
                            NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_TYPE(1)=DOMAIN_LOCAL_BOUNDARY
                            DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%numberOfDerivatives
                              DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%numberOfVersions
                                ny=MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%dofIndex(version_idx)
                                LOCAL_DOF_NUMBERS(domain_no)=LOCAL_DOF_NUMBERS(domain_no)+1
                                DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS=1
                                DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(1)=LOCAL_DOF_NUMBERS(domain_no)
                                DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER(1)=domain_no
                                DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE(1)=DOMAIN_LOCAL_BOUNDARY
                              ENDDO !version_idx
                            ENDDO !derivative_idx
                          ELSE
                            !This domain doesn't take this node, so make it a ghost node in this domain 
                            CALL LIST_ITEM_ADD(GHOST_NODES_LIST(domain_no)%PTR,node_idx,ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          !The node has already been assigned to a domain so it must be a ghost node in this domain
                          CALL LIST_ITEM_ADD(GHOST_NODES_LIST(domain_no)%PTR,node_idx,ERR,ERROR,*999)
                        ENDIF
                      ENDDO !domain_idx
                    ENDIF
                  ENDDO !node_idx
                  DEALLOCATE(NUMBER_INTERNAL_NODES)
                  
                  !Calculate ghost node and dof mappings

                 ! loop over all domains
                 DO domain_idx = 0,DECOMPOSITION%NUMBER_OF_DOMAINS-1
                  
                    ! determine all ghost nodes of this domain 
                   CALL LIST_REMOVE_DUPLICATES(GHOST_NODES_LIST(domain_idx)%PTR,ERR,ERROR,*999)
                   CALL LIST_DETACH_AND_DESTROY(GHOST_NODES_LIST(domain_idx)%PTR,NUMBER_OF_GHOST_NODES,GHOST_NODES,ERR,ERROR,*999)

                   ! loop over the ghost nodes of this domain
                   DO no_ghost_node = 1,NUMBER_OF_GHOST_NODES 
                      ghost_node=GHOST_NODES(no_ghost_node)
                      ! add the 'ghost in this domain' information to GLOBAL_TO_LOCAL_MAP of this node 
                      LOCAL_NODE_NUMBERS(domain_idx)=LOCAL_NODE_NUMBERS(domain_idx)+1
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_node)%NUMBER_OF_DOMAINS= &
                        & NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_node)%NUMBER_OF_DOMAINS+1
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_node)%LOCAL_NUMBER( &
                        & NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_node)%NUMBER_OF_DOMAINS)= &
                        & LOCAL_NODE_NUMBERS(domain_idx)
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_node)%DOMAIN_NUMBER( &
                        & NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_node)%NUMBER_OF_DOMAINS)=domain_idx
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_node)%LOCAL_TYPE( &
                        & NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_node)%NUMBER_OF_DOMAINS)= &
                        & DOMAIN_LOCAL_GHOST
                      ! update DOFS_MAPPING 
                      DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(ghost_node)%numberOfDerivatives
                        DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(ghost_node)%DERIVATIVES(derivative_idx)%numberOfVersions
                          ny=MESH_TOPOLOGY%NODES%NODES(ghost_node)%DERIVATIVES(derivative_idx)%dofIndex(version_idx)
                          LOCAL_DOF_NUMBERS(domain_idx)=LOCAL_DOF_NUMBERS(domain_idx)+1
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS= &
                            & DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS+1
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER( &
                            & DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS)= &
                            & LOCAL_DOF_NUMBERS(domain_idx)
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER( &
                            & DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS)=domain_idx
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE( &
                            & DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS)= &
                            & DOMAIN_LOCAL_GHOST
                        ENDDO !version_idx
                      ENDDO !derivative_idx
                    ENDDO !no_ghost_node
                    DEALLOCATE(GHOST_NODES)
                  ENDDO !domain_idx
                  
                    !Check decomposition and check that each domain has a node in it.
                  ALLOCATE(NODE_COUNT(0:numberOfComputationalNodes-1),STAT=ERR) 
                  IF(ERR/=0) CALL FlagError("Could not allocate node count.",ERR,ERROR,*999)
                  NODE_COUNT=0
                  DO node_idx=1,MESH_TOPOLOGY%NODES%numberOfNodes
                    no_computational_node=DOMAIN%NODE_DOMAIN(node_idx)
                    IF(no_computational_node>=0.AND.no_computational_node<numberOfComputationalNodes) THEN
                      NODE_COUNT(no_computational_node)=NODE_COUNT(no_computational_node)+1
                    ELSE
                      LOCAL_ERROR="The computational node number of "// &
                        & TRIM(NUMBER_TO_VSTRING(no_computational_node,"*",ERR,ERROR))// &
                        & " for node number "//TRIM(NUMBER_TO_VSTRING(node_idx,"*",ERR,ERROR))// &
                        & " is invalid. The computational node number must be between 0 and "// &
                        & TRIM(NUMBER_TO_VSTRING(numberOfComputationalNodes-1,"*",ERR,ERROR))//"."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ENDDO !node_idx
                  DO no_computational_node=0,numberOfComputationalNodes-1
                    IF(NODE_COUNT(no_computational_node)==0) THEN
                      LOCAL_ERROR="Invalid decomposition. There are no nodes in computational node "// &
                        & TRIM(NUMBER_TO_VSTRING(no_computational_node,"*",ERR,ERROR))//"."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ENDDO !no_computational_node
                  DEALLOCATE(NODE_COUNT)
!               END IF

                  DEALLOCATE(GHOST_NODES_LIST)
                  DEALLOCATE(LOCAL_NODE_NUMBERS)
                  
                  !Calculate node and dof local to global maps from global to local map
                  CALL DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(NODES_MAPPING,ERR,ERROR,*999)
                  CALL DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(DOFS_MAPPING,ERR,ERROR,*999)
                  
                ELSE
                  CALL FlagError("Domain mesh is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Domain decomposition is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Domain mappings elements is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Domain mappings dofs is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Domain mappings nodes is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Domain mappings is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain is not associated.",ERR,ERROR,*998)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Node decomposition :",ERR,ERROR,*999)
      DO node_idx=1,MESH_TOPOLOGY%NODES%numberOfNodes
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Node = ",node_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Domain = ",DOMAIN%NODE_DOMAIN(node_idx),ERR,ERROR,*999)
      ENDDO !node_idx
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Node mappings :",ERR,ERROR,*999)
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Global to local map :",ERR,ERROR,*999)
      DO node_idx=1,MESH_TOPOLOGY%NODES%numberOfNodes
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global node = ",node_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of domains  = ", &
          & NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)% &
          & NUMBER_OF_DOMAINS,8,8,NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_NUMBER, &
          & '("      Local number :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)% &
            & NUMBER_OF_DOMAINS,8,8,NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%DOMAIN_NUMBER, &
            & '("      Domain number:",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)% &
          & NUMBER_OF_DOMAINS,8,8,NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_TYPE, &
          & '("      Local type   :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
      ENDDO !node_idx     
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Local to global map :",ERR,ERROR,*999)
      DO node_idx=1,NODES_MAPPING%TOTAL_NUMBER_OF_LOCAL
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Local node = ",node_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Global node = ", &
          & NODES_MAPPING%LOCAL_TO_GLOBAL_MAP(node_idx),ERR,ERROR,*999)
      ENDDO !node_idx
      IF(DIAGNOSTICS2) THEN
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Internal nodes :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of internal nodes = ", &
          & NODES_MAPPING%NUMBER_OF_INTERNAL,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%NUMBER_OF_INTERNAL,8,8, &
          & NODES_MAPPING%DOMAIN_LIST(NODES_MAPPING%INTERNAL_START:NODES_MAPPING%INTERNAL_FINISH), &
          & '("    Internal nodes:",8(X,I7))','(19X,8(X,I7))',ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary nodes :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of boundary nodes = ", &
          & NODES_MAPPING%NUMBER_OF_BOUNDARY,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%NUMBER_OF_BOUNDARY,8,8, &
          & NODES_MAPPING%DOMAIN_LIST(NODES_MAPPING%BOUNDARY_START:NODES_MAPPING%BOUNDARY_FINISH), &
          & '("    Boundary nodes:",8(X,I7))','(19X,8(X,I7))',ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Ghost nodes :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of ghost nodes = ", &
          & NODES_MAPPING%NUMBER_OF_GHOST,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%NUMBER_OF_GHOST,8,8, &
          & NODES_MAPPING%DOMAIN_LIST(NODES_MAPPING%GHOST_START:NODES_MAPPING%GHOST_FINISH), &
          & '("    Ghost nodes   :",8(X,I7))','(19X,8(X,I7))',ERR,ERROR,*999)
      ENDIF
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Adjacent domains :",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of adjacent domains = ", &
        & NODES_MAPPING%NUMBER_OF_ADJACENT_DOMAINS,ERR,ERROR,*999)

      DO domain_idx=1,NODES_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Adjacent domain idx : ",domain_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Domain number = ", &
          & NODES_MAPPING%ADJACENT_DOMAINS(domain_idx)%DOMAIN_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of send ghosts    = ", &
          & NODES_MAPPING%ADJACENT_DOMAINS(domain_idx)%NUMBER_OF_SEND_GHOSTS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%ADJACENT_DOMAINS(domain_idx)% &
          & NUMBER_OF_SEND_GHOSTS,6,6,NODES_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_SEND_INDICES, &
          & '("      Local send ghost indicies       :",6(X,I7))','(39X,6(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of recieve ghosts = ", &
          & NODES_MAPPING%ADJACENT_DOMAINS(domain_idx)%NUMBER_OF_RECEIVE_GHOSTS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%ADJACENT_DOMAINS(domain_idx)% &
          & NUMBER_OF_RECEIVE_GHOSTS,6,6,NODES_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_RECEIVE_INDICES, &
          & '("      Local receive ghost indicies    :",6(X,I7))','(39X,6(X,I7))',ERR,ERROR,*999)              
      ENDDO !domain_idx
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Dofs mappings :",ERR,ERROR,*999)
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Global to local map :",ERR,ERROR,*999)
      DO ny=1,MESH_TOPOLOGY%DOFS%numberOfDofs
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global dof = ",ny,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of domains  = ", &
          & DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)% &
          & NUMBER_OF_DOMAINS,8,8,DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER, &
          & '("      Local number :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)% &
            & NUMBER_OF_DOMAINS,8,8,DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER, &
            & '("      Domain number:",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)% &
          & NUMBER_OF_DOMAINS,8,8,DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE, &
          & '("      Local type   :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
      ENDDO !ny   
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Local to global map :",ERR,ERROR,*999)
      DO ny=1,DOFS_MAPPING%TOTAL_NUMBER_OF_LOCAL
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Local dof = ",ny,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Global dof = ", &
          & DOFS_MAPPING%LOCAL_TO_GLOBAL_MAP(ny),ERR,ERROR,*999)
      ENDDO !node_idx
      IF(DIAGNOSTICS2) THEN
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Internal dofs :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of internal dofs = ", &
          & DOFS_MAPPING%NUMBER_OF_INTERNAL,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%NUMBER_OF_INTERNAL,8,8, &
          & DOFS_MAPPING%DOMAIN_LIST(DOFS_MAPPING%INTERNAL_START:DOFS_MAPPING%INTERNAL_FINISH), &
          & '("    Internal dofs:",8(X,I7))','(18X,8(X,I7))',ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary dofs :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of boundary dofs = ", &
          & DOFS_MAPPING%NUMBER_OF_BOUNDARY,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%NUMBER_OF_BOUNDARY,8,8, &
          & DOFS_MAPPING%DOMAIN_LIST(DOFS_MAPPING%BOUNDARY_START:DOFS_MAPPING%BOUNDARY_FINISH), &
          & '("    Boundary dofs:",8(X,I7))','(18X,8(X,I7))',ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Ghost dofs :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of ghost dofs = ", &
          & DOFS_MAPPING%NUMBER_OF_GHOST,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%NUMBER_OF_GHOST,8,8, &
          & DOFS_MAPPING%DOMAIN_LIST(DOFS_MAPPING%GHOST_START:DOFS_MAPPING%GHOST_FINISH), &
          & '("    Ghost dofs   :",8(X,I7))','(18X,8(X,I7))',ERR,ERROR,*999)
      ENDIF
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Adjacent domains :",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of adjacent domains = ", &
        & DOFS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS,ERR,ERROR,*999)

      DO domain_idx=1,DOFS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Adjacent domain idx : ",domain_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Domain number = ", &
          & DOFS_MAPPING%ADJACENT_DOMAINS(domain_idx)%DOMAIN_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of send ghosts    = ", &
          & DOFS_MAPPING%ADJACENT_DOMAINS(domain_idx)%NUMBER_OF_SEND_GHOSTS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%ADJACENT_DOMAINS(domain_idx)% &
          & NUMBER_OF_SEND_GHOSTS,6,6,DOFS_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_SEND_INDICES, &
          & '("      Local send ghost indicies       :",6(X,I7))','(39X,6(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of recieve ghosts = ", &
          & DOFS_MAPPING%ADJACENT_DOMAINS(domain_idx)%NUMBER_OF_RECEIVE_GHOSTS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%ADJACENT_DOMAINS(domain_idx)% &
          & NUMBER_OF_RECEIVE_GHOSTS,6,6,DOFS_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_RECEIVE_INDICES, &
          & '("      Local receive ghost indicies    :",6(X,I7))','(39X,6(X,I7))',ERR,ERROR,*999)              
      ENDDO !domain_idx
    ENDIF
    
    EXITS("DOMAIN_MAPPINGS_NODES_DOFS_CALCULATE")
    RETURN
999 IF(ALLOCATED(DOMAINS)) DEALLOCATE(DOMAINS)
    IF(ALLOCATED(ALL_DOMAINS)) DEALLOCATE(ALL_DOMAINS)
    IF(ALLOCATED(GHOST_NODES)) DEALLOCATE(GHOST_NODES)
    IF(ALLOCATED(NUMBER_INTERNAL_NODES)) DEALLOCATE(NUMBER_INTERNAL_NODES)
    IF(ALLOCATED(NUMBER_BOUNDARY_NODES)) DEALLOCATE(NUMBER_BOUNDARY_NODES)
    IF(ASSOCIATED(DOMAIN%MAPPINGS%NODES)) CALL DOMAIN_MAPPINGS_NODES_FINALISE(DOMAIN%MAPPINGS,DUMMY_ERR,DUMMY_ERROR,*998)
998 IF(ASSOCIATED(DOMAIN%MAPPINGS%DOFS)) CALL DOMAIN_MAPPINGS_DOFS_FINALISE(DOMAIN%MAPPINGS,DUMMY_ERR,DUMMY_ERROR,*997)
997 ERRORSEXITS("DOMAIN_MAPPINGS_NODES_DOFS_CALCULATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_NODES_DOFS_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Finalises the node mapping in the given domain mappings. \todo pass in the nodes mapping
  SUBROUTINE DOMAIN_MAPPINGS_NODES_FINALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS !<A pointer to the domain mapping to finalise the nodes for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("DOMAIN_MAPPINGS_NODES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%NODES)) THEN
        CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(DOMAIN_MAPPINGS%NODES,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_MAPPINGS_NODES_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_NODES_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_NODES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the node mapping in the given domain mapping. \todo finalise on error
  SUBROUTINE DOMAIN_MAPPINGS_NODES_INITIALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS !<A pointer to the domain mappings to initialise the nodes for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("DOMAIN_MAPPINGS_NODES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%NODES)) THEN
        CALL FlagError("Domain nodes mappings are already associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(DOMAIN_MAPPINGS%NODES,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate domain mappings nodes.",ERR,ERROR,*999)
        CALL DOMAIN_MAPPINGS_MAPPING_INITIALISE(DOMAIN_MAPPINGS%NODES,DOMAIN_MAPPINGS%DOMAIN%DECOMPOSITION%workGroup, &
          & ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_MAPPINGS_NODES_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_MAPPINGS_NODES_INITIALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_NODES_INITIALISE

  !
  !================================================================================================================================
  !

  !>Calculates the domain topology.
  SUBROUTINE DOMAIN_TOPOLOGY_CALCULATE(TOPOLOGY,ERR,ERROR,*)

   !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to calculate.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ne,np
    TYPE(BASIS_TYPE), POINTER :: BASIS

    ENTERS("DOMAIN_TOPOLOGY_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      !Find maximum number of element parameters for all elements in the domain toplogy
      TOPOLOGY%ELEMENTS%MAXIMUM_NUMBER_OF_ELEMENT_PARAMETERS=-1
      DO ne=1,TOPOLOGY%ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
        BASIS=>TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS
        IF(ASSOCIATED(BASIS)) THEN
          IF(BASIS%NUMBER_OF_ELEMENT_PARAMETERS>TOPOLOGY%ELEMENTS%MAXIMUM_NUMBER_OF_ELEMENT_PARAMETERS) &
            & TOPOLOGY%ELEMENTS%MAXIMUM_NUMBER_OF_ELEMENT_PARAMETERS=BASIS%NUMBER_OF_ELEMENT_PARAMETERS
        ELSE
          CALL FlagError("Basis is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDDO !ne
      !Find maximum number of derivatives for all nodes in the domain toplogy
      TOPOLOGY%NODES%MAXIMUM_NUMBER_OF_DERIVATIVES=-1
      DO np=1,TOPOLOGY%NODES%TOTAL_NUMBER_OF_NODES
        IF(TOPOLOGY%NODES%NODES(np)%NUMBER_OF_DERIVATIVES>TOPOLOGY%NODES%MAXIMUM_NUMBER_OF_DERIVATIVES) &
            & TOPOLOGY%NODES%MAXIMUM_NUMBER_OF_DERIVATIVES=TOPOLOGY%NODES%NODES(np)%NUMBER_OF_DERIVATIVES
      ENDDO !np
      !Calculate the elements surrounding the nodes in the domain topology
      CALL DomainTopology_NodesSurroundingElementsCalculate(TOPOLOGY,ERR,ERROR,*999)
    ELSE
      CALL FlagError("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_CALCULATE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_CALCULATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Initialises the local domain topology from the mesh topology.
  SUBROUTINE DOMAIN_TOPOLOGY_INITIALISE_FROM_MESH(DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<A pointer to the domain to initialise the domain topology from the mesh topology for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: local_element,global_element,local_node,global_node,version_idx,derivative_idx,node_idx,dof_idx, &
      & component_idx
    INTEGER(INTG) :: ne,nn,nkk,INSERT_STATUS
    LOGICAL :: FOUND
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS
    TYPE(MeshElementsType), POINTER :: MESH_ELEMENTS
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(MeshNodesType), POINTER :: MESH_NODES
    TYPE(DOMAIN_DOFS_TYPE), POINTER :: DOMAIN_DOFS

    ENTERS("DOMAIN_TOPOLOGY_INITIALISE_FROM_MESH",ERR,ERROR,*999)
    
    IF(ASSOCIATED(DOMAIN)) THEN
      IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
        IF(ASSOCIATED(DOMAIN%MAPPINGS)) THEN
          IF(ASSOCIATED(DOMAIN%MESH)) THEN
            MESH=>DOMAIN%MESH
            component_idx=DOMAIN%MESH_COMPONENT_NUMBER
            IF(ASSOCIATED(MESH%TOPOLOGY(component_idx)%PTR)) THEN
              MESH_ELEMENTS=>MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS
              DOMAIN_ELEMENTS=>DOMAIN%TOPOLOGY%ELEMENTS
              MESH_NODES=>MESH%TOPOLOGY(component_idx)%PTR%NODES
              DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
              DOMAIN_DOFS=>DOMAIN%TOPOLOGY%DOFS
              ALLOCATE(DOMAIN_ELEMENTS%ELEMENTS(DOMAIN%MAPPINGS%ELEMENTS%TOTAL_NUMBER_OF_LOCAL),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate domain elements elements.",ERR,ERROR,*999)
              DOMAIN_ELEMENTS%NUMBER_OF_ELEMENTS=DOMAIN%MAPPINGS%ELEMENTS%NUMBER_OF_LOCAL
              DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS=DOMAIN%MAPPINGS%ELEMENTS%TOTAL_NUMBER_OF_LOCAL
              DOMAIN_ELEMENTS%NUMBER_OF_GLOBAL_ELEMENTS=DOMAIN%MAPPINGS%ELEMENTS%NUMBER_OF_GLOBAL
              ALLOCATE(DOMAIN_NODES%NODES(DOMAIN%MAPPINGS%NODES%TOTAL_NUMBER_OF_LOCAL),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate domain nodes nodes.",ERR,ERROR,*999)
              DOMAIN_NODES%NUMBER_OF_NODES=DOMAIN%MAPPINGS%NODES%NUMBER_OF_LOCAL
              DOMAIN_NODES%TOTAL_NUMBER_OF_NODES=DOMAIN%MAPPINGS%NODES%TOTAL_NUMBER_OF_LOCAL
              DOMAIN_NODES%NUMBER_OF_GLOBAL_NODES=DOMAIN%MAPPINGS%NODES%NUMBER_OF_GLOBAL
              ALLOCATE(DOMAIN_DOFS%DOF_INDEX(3,DOMAIN%MAPPINGS%DOFS%TOTAL_NUMBER_OF_LOCAL),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate domain dofs dof index.",ERR,ERROR,*999)
              DOMAIN_DOFS%NUMBER_OF_DOFS=DOMAIN%MAPPINGS%DOFS%NUMBER_OF_LOCAL
              DOMAIN_DOFS%TOTAL_NUMBER_OF_DOFS=DOMAIN%MAPPINGS%DOFS%TOTAL_NUMBER_OF_LOCAL
              DOMAIN_DOFS%NUMBER_OF_GLOBAL_DOFS=DOMAIN%MAPPINGS%DOFS%NUMBER_OF_GLOBAL
              !Loop over the domain nodes and calculate the parameters from the mesh nodes
              CALL TREE_CREATE_START(DOMAIN_NODES%NODES_TREE,ERR,ERROR,*999)
              CALL TREE_INSERT_TYPE_SET(DOMAIN_NODES%NODES_TREE,TREE_NO_DUPLICATES_ALLOWED,ERR,ERROR,*999)
              CALL TREE_CREATE_FINISH(DOMAIN_NODES%NODES_TREE,ERR,ERROR,*999)
              dof_idx=0
              DO local_node=1,DOMAIN_NODES%TOTAL_NUMBER_OF_NODES
                CALL DOMAIN_TOPOLOGY_NODE_INITIALISE(DOMAIN_NODES%NODES(local_node),ERR,ERROR,*999)
                global_node=DOMAIN%MAPPINGS%NODES%LOCAL_TO_GLOBAL_MAP(local_node)
                DOMAIN_NODES%NODES(local_node)%LOCAL_NUMBER=local_node
                DOMAIN_NODES%NODES(local_node)%MESH_NUMBER=global_node
                DOMAIN_NODES%NODES(local_node)%GLOBAL_NUMBER=MESH_NODES%NODES(global_node)%globalNumber
                DOMAIN_NODES%NODES(local_node)%USER_NUMBER=MESH_NODES%NODES(global_node)%userNumber
                CALL TREE_ITEM_INSERT(DOMAIN_NODES%NODES_TREE,DOMAIN_NODES%NODES(local_node)%USER_NUMBER,local_node, &
                  & INSERT_STATUS,ERR,ERROR,*999)
                DOMAIN_NODES%NODES(local_node)%NUMBER_OF_SURROUNDING_ELEMENTS=0
                NULLIFY(DOMAIN_NODES%NODES(local_node)%SURROUNDING_ELEMENTS)
                DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES=MESH_NODES%NODES(global_node)%numberOfDerivatives
                ALLOCATE(DOMAIN_NODES%NODES(local_node)%DERIVATIVES(MESH_NODES%NODES(global_node)%numberOfDerivatives),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate domain node derivatives.",ERR,ERROR,*999)
                DO derivative_idx=1,DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES
                  CALL DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE( &
                    & DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx),ERR,ERROR,*999)
                  DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx)%GLOBAL_DERIVATIVE_INDEX= & 
                    & MESH_NODES%NODES(global_node)%DERIVATIVES(derivative_idx)%globalDerivativeIndex
                  DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx)%PARTIAL_DERIVATIVE_INDEX= &
                    & MESH_NODES%NODES(global_node)%DERIVATIVES(derivative_idx)%partialDerivativeIndex
                  DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx)%numberOfVersions= & 
                    & MESH_NODES%NODES(global_node)%DERIVATIVES(derivative_idx)%numberOfVersions
                  ALLOCATE(DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx)%userVersionNumbers( &
                    & MESH_NODES%NODES(global_node)%DERIVATIVES(derivative_idx)%numberOfVersions),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate node derivative version numbers.",ERR,ERROR,*999)
                  DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx)%userVersionNumbers(1: &
                    & MESH_NODES%NODES(global_node)%DERIVATIVES(derivative_idx)%numberOfVersions)= &
                    & MESH_NODES%NODES(global_node)%DERIVATIVES(derivative_idx)%userVersionNumbers(1: &
                    & MESH_NODES%NODES(global_node)%DERIVATIVES(derivative_idx)%numberOfVersions)
                  ALLOCATE(DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx)%DOF_INDEX( &
                    & MESH_NODES%NODES(global_node)%DERIVATIVES(derivative_idx)%numberOfVersions),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate node dervative versions dof index.",ERR,ERROR,*999)
                  DO version_idx=1,DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx)%numberOfVersions
                    dof_idx=dof_idx+1
                    DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx)%DOF_INDEX(version_idx)=dof_idx
                    DOMAIN_DOFS%DOF_INDEX(1,dof_idx)=version_idx ! Fortran runtime error: Index '1' of dimension 2 of array 'domain_dofs%dof_index' above upper bound of 0
                    DOMAIN_DOFS%DOF_INDEX(2,dof_idx)=derivative_idx
                    DOMAIN_DOFS%DOF_INDEX(3,dof_idx)=local_node
                  ENDDO !version_idx
                ENDDO !derivative_idx
                DOMAIN_NODES%NODES(local_node)%BOUNDARY_NODE=MESH_NODES%NODES(global_node)%boundaryNode
              ENDDO !local_node
              !Loop over the domain elements and renumber from the mesh elements
              DO local_element=1,DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                CALL DOMAIN_TOPOLOGY_ELEMENT_INITIALISE(DOMAIN_ELEMENTS%ELEMENTS(local_element),ERR,ERROR,*999)
                global_element=DOMAIN%MAPPINGS%ELEMENTS%LOCAL_TO_GLOBAL_MAP(local_element)               
                BASIS=>MESH_ELEMENTS%ELEMENTS(global_element)%BASIS
                DOMAIN_ELEMENTS%ELEMENTS(local_element)%NUMBER=local_element
                DOMAIN_ELEMENTS%ELEMENTS(local_element)%BASIS=>BASIS
                ALLOCATE(DOMAIN_ELEMENTS%ELEMENTS(local_element)%ELEMENT_NODES(BASIS%NUMBER_OF_NODES),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate domain elements element nodes.",ERR,ERROR,*999)
                ALLOCATE(DOMAIN_ELEMENTS%ELEMENTS(local_element)%ELEMENT_DERIVATIVES(BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                  & BASIS%NUMBER_OF_NODES),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate domain elements element derivatives.",ERR,ERROR,*999)
                ALLOCATE(DOMAIN_ELEMENTS%ELEMENTS(local_element)%elementVersions(BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                  & BASIS%NUMBER_OF_NODES),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate domain elements element versions.",ERR,ERROR,*999)
                DO nn=1,BASIS%NUMBER_OF_NODES
                  global_node=MESH_ELEMENTS%ELEMENTS(global_element)%MESH_ELEMENT_NODES(nn)
                  local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(global_node)%LOCAL_NUMBER(1)
                  DOMAIN_ELEMENTS%ELEMENTS(local_element)%ELEMENT_NODES(nn)=local_node
                  DO derivative_idx=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
                    !Find equivalent node derivative by matching partial derivative index
                    !/todo Take a look at this - is it needed any more?
                    FOUND=.FALSE.
                    DO nkk=1,DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES
                      IF(DOMAIN_NODES%NODES(local_node)%DERIVATIVES(nkk)%PARTIAL_DERIVATIVE_INDEX == &
                        & BASIS%PARTIAL_DERIVATIVE_INDEX(derivative_idx,nn)) THEN
                        FOUND=.TRUE.
                        EXIT
                      ENDIF
                    ENDDO !nkk
                    IF(FOUND) THEN
                      DOMAIN_ELEMENTS%ELEMENTS(local_element)%ELEMENT_DERIVATIVES(derivative_idx,nn)=nkk
                      DOMAIN_ELEMENTS%ELEMENTS(local_element)%elementVersions(derivative_idx,nn) = & 
                        & MESH_ELEMENTS%ELEMENTS(global_element)%USER_ELEMENT_NODE_VERSIONS(derivative_idx,nn)
                    ELSE
                      CALL FlagError("Could not find equivalent node derivative",ERR,ERROR,*999)
                    ENDIF
                  ENDDO !derivative_idx
                ENDDO !nn
              ENDDO !local_element                       
            ELSE
              CALL FlagError("Mesh topology is not associated",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Mesh is not associated",ERR,ERROR,*999)

          ENDIF
        ELSE
          CALL FlagError("Domain mapping is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Domain topology is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain is not associated",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Initialised domain topology :",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Total number of domain nodes = ",DOMAIN_NODES%TOTAL_NUMBER_OF_NODES, &
        & ERR,ERROR,*999)
      DO node_idx=1,DOMAIN_NODES%TOTAL_NUMBER_OF_NODES
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Node number = ",DOMAIN_NODES%NODES(node_idx)%LOCAL_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Node mesh number = ",DOMAIN_NODES%NODES(node_idx)%MESH_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Node global number = ",DOMAIN_NODES%NODES(node_idx)%GLOBAL_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Node user number = ",DOMAIN_NODES%NODES(node_idx)%USER_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of derivatives = ", &
          & DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES,ERR,ERROR,*999)
        DO derivative_idx=1,DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Node local derivative number = ",derivative_idx,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Global derivative index = ", &
            & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%GLOBAL_DERIVATIVE_INDEX,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Partial derivative index = ", &
            & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%PARTIAL_DERIVATIVE_INDEX,ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
            & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%numberOfVersions,4,4, &
            & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%DOF_INDEX, &
            & '("        Degree-of-freedom index(version_idx)  :",4(X,I9))','(36X,4(X,I9))',ERR,ERROR,*999)
        ENDDO !derivative_idx
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Boundary node = ", &
          & DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE,ERR,ERROR,*999)
     ENDDO !node_idx
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Total Number of domain dofs = ",DOMAIN_DOFS%TOTAL_NUMBER_OF_DOFS, &
        & ERR,ERROR,*999)
      DO dof_idx=1,DOMAIN_DOFS%TOTAL_NUMBER_OF_DOFS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Dof number = ",dof_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,3,3, &
          & DOMAIN_DOFS%DOF_INDEX(:,dof_idx),'("    Degree-of-freedom index :",3(X,I9))','(29X,3(X,I9))', &
          & ERR,ERROR,*999)
      ENDDO !dof_idx
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Total number of domain elements = ", &
        & DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS,ERR,ERROR,*999)
      DO ne=1,DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Element number = ",DOMAIN_ELEMENTS%ELEMENTS(ne)%NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Basis user number = ", &
          & DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS%USER_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of local nodes = ", &
          & DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS%NUMBER_OF_NODES,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS%NUMBER_OF_NODES,8,8, &
          & DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES,'("    Element nodes(nn) :",8(X,I9))','(23X,8(X,I9))', &
          & ERR,ERROR,*999)
        DO nn=1,DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS%NUMBER_OF_NODES
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Local node number : ",nn,ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS%NUMBER_OF_DERIVATIVES(nn),8,8, &
            & DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_DERIVATIVES(:,nn), &
            & '("        Element derivatives :",8(X,I2))','(29X,8(X,I2))',ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS%NUMBER_OF_DERIVATIVES(nn),8,8, &
            & DOMAIN_ELEMENTS%ELEMENTS(ne)%elementVersions(:,nn), &
            & '("        Element versions    :",8(X,I2))','(29X,8(X,I2))',ERR,ERROR,*999)
        ENDDO !nn
      ENDDO !ne
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_INITIALISE_FROM_MESH")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_INITIALISE_FROM_MESH",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE DOMAIN_TOPOLOGY_INITIALISE_FROM_MESH
  
  !
  !================================================================================================================================
  !

  !>Finalises the dofs in the given domain topology. \todo pass in the dofs topolgy
  SUBROUTINE DOMAIN_TOPOLOGY_DOFS_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to finalise the dofs for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_DOFS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%DOFS)) THEN
        IF(ALLOCATED(TOPOLOGY%DOFS%DOF_INDEX)) DEALLOCATE(TOPOLOGY%DOFS%DOF_INDEX)
        DEALLOCATE(TOPOLOGY%DOFS)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_TOPOLOGY_DOFS_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_DOFS_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_DOFS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the dofs data structures for a domain topology. \todo finalise on exit
  SUBROUTINE DOMAIN_TOPOLOGY_DOFS_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to initialise the dofs for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_DOFS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%DOFS)) THEN
        CALL FlagError("Domain already has topology dofs associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%DOFS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology dofs",ERR,ERROR,*999)
        TOPOLOGY%DOFS%NUMBER_OF_DOFS=0
        TOPOLOGY%DOFS%TOTAL_NUMBER_OF_DOFS=0
        TOPOLOGY%DOFS%NUMBER_OF_GLOBAL_DOFS=0
        TOPOLOGY%DOFS%DOMAIN=>TOPOLOGY%DOMAIN
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_DOFS_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_DOFS_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_DOFS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the given domain topology element.
  SUBROUTINE DOMAIN_TOPOLOGY_ELEMENT_FINALISE(ELEMENT,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_ELEMENT_TYPE) :: ELEMENT !<The domain element to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(ELEMENT%ELEMENT_NODES)) DEALLOCATE(ELEMENT%ELEMENT_NODES)
    IF(ALLOCATED(ELEMENT%ELEMENT_DERIVATIVES)) DEALLOCATE(ELEMENT%ELEMENT_DERIVATIVES)
    IF(ALLOCATED(ELEMENT%elementVersions)) DEALLOCATE(ELEMENT%elementVersions)
 
    EXITS("DOMAIN_TOPOLOGY_ELEMENT_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_ELEMENT_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the given domain topology element.
  SUBROUTINE DOMAIN_TOPOLOGY_ELEMENT_INITIALISE(ELEMENT,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_ELEMENT_TYPE) :: ELEMENT !<The domain element to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_ELEMENT_INITIALISE",ERR,ERROR,*999)

    ELEMENT%NUMBER=0
    NULLIFY(ELEMENT%BASIS)
  
    EXITS("DOMAIN_TOPOLOGY_ELEMENT_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_ELEMENT_INITIALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_ELEMENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the elements in the given domain topology. \todo pass in the domain elements
  SUBROUTINE DOMAIN_TOPOLOGY_ELEMENTS_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to finalise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ne

    ENTERS("DOMAIN_TOPOLOGY_ELEMENTS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        DO ne=1,TOPOLOGY%ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
          CALL DOMAIN_TOPOLOGY_ELEMENT_FINALISE(TOPOLOGY%ELEMENTS%ELEMENTS(ne),ERR,ERROR,*999)
        ENDDO !ne
        IF(ASSOCIATED(TOPOLOGY%ELEMENTS%ELEMENTS)) DEALLOCATE(TOPOLOGY%ELEMENTS%ELEMENTS)
        DEALLOCATE(TOPOLOGY%ELEMENTS)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_TOPOLOGY_ELEMENTS_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_ELEMENTS_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_ELEMENTS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the element data structures for a domain topology. \todo finalise on error
  SUBROUTINE DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to initialise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        CALL FlagError("Decomposition already has topology elements associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%ELEMENTS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology elements",ERR,ERROR,*999)
        TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS=0
        TOPOLOGY%ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS=0
        TOPOLOGY%ELEMENTS%NUMBER_OF_GLOBAL_ELEMENTS=0
        TOPOLOGY%ELEMENTS%DOMAIN=>TOPOLOGY%DOMAIN
        NULLIFY(TOPOLOGY%ELEMENTS%ELEMENTS)
        TOPOLOGY%ELEMENTS%MAXIMUM_NUMBER_OF_ELEMENT_PARAMETERS=0
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE
  
  
  !
  !================================================================================================================================
  !

  !>Finalises the topology in the given domain. \todo pass in domain topology
  SUBROUTINE DOMAIN_TOPOLOGY_FINALISE(DOMAIN,ERR,ERROR,*)

   !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<A pointer to the domain to finalise the topology for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      CALL DOMAIN_TOPOLOGY_NODES_FINALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      CALL DOMAIN_TOPOLOGY_DOFS_FINALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      CALL DOMAIN_TOPOLOGY_ELEMENTS_FINALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      CALL DOMAIN_TOPOLOGY_LINES_FINALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      CALL DOMAIN_TOPOLOGY_FACES_FINALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      DEALLOCATE(DOMAIN%TOPOLOGY)
    ELSE
      CALL FlagError("Domain is not associated",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_TOPOLOGY_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the topology for a given domain. \todo finalise on error
  SUBROUTINE DOMAIN_TOPOLOGY_INITIALISE(DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !A pointer to the domain to initialise the topology for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
        CALL FlagError("Domain already has topology associated",ERR,ERROR,*999)
      ELSE
        !Allocate domain topology
        ALLOCATE(DOMAIN%TOPOLOGY,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Domain topology could not be allocated",ERR,ERROR,*999)
        DOMAIN%TOPOLOGY%DOMAIN=>DOMAIN
        NULLIFY(DOMAIN%TOPOLOGY%ELEMENTS)
        NULLIFY(DOMAIN%TOPOLOGY%NODES)
        NULLIFY(DOMAIN%TOPOLOGY%DOFS)
        NULLIFY(DOMAIN%TOPOLOGY%LINES)
        NULLIFY(DOMAIN%TOPOLOGY%FACES)
        !Initialise the topology components
        CALL DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
        CALL DOMAIN_TOPOLOGY_NODES_INITIALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
        CALL DOMAIN_TOPOLOGY_DOFS_INITIALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
        CALL DOMAIN_TOPOLOGY_LINES_INITIALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
        CALL DOMAIN_TOPOLOGY_FACES_INITIALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
        !Initialise the domain topology from the domain mappings and the mesh it came from
        CALL DOMAIN_TOPOLOGY_INITIALISE_FROM_MESH(DOMAIN,ERR,ERROR,*999)
        !Calculate the topological information.
        CALL DOMAIN_TOPOLOGY_CALCULATE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Get the basis for an element in the domain identified by its user number
  SUBROUTINE DomainTopology_ElementBasisGet(domainTopology,userElementNumber,basis,err,error,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: domainTopology !<A pointer to the domain topology to get the element basis for
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The element user number to get the basis for
    TYPE(BASIS_TYPE), POINTER, INTENT(OUT) :: basis !<On return, a pointer to the basis for the element.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    LOGICAL :: userElementExists,ghostElement
    INTEGER(INTG) :: localElementNumber
    TYPE(VARYING_STRING) :: localError

    ENTERS("DomainTopology_ElementBasisGet",err,error,*999)

    IF(.NOT.ASSOCIATED(domainTopology)) CALL FlagError("Domain topology is not associated.",err,error,*999)
    IF(ASSOCIATED(basis)) CALL FlagError("Basis is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(domainTopology%elements)) CALL FlagError("Domain topology elements is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(domainTopology%elements%elements)) &
      & CALL FlagError("Domain topology elements elements is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(domainTopology%domain)) CALL FlagError("Domain topology domain is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(domainTopology%domain%decomposition)) &
      & CALL FlagError("Domain topology domain decomposition is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(domainTopology%domain%decomposition%topology)) &
      & CALL FlagError("Domain topology domain decomposition topology is not associated.",err,error,*999)
    
    CALL DecompositionTopology_ElementCheckExists(domainTopology%domain%decomposition%topology,userElementNumber, &
      & userElementExists,localElementNumber,ghostElement,err,error,*999)
    IF(.NOT.userElementExists) THEN
      CALL FlagError("The specified user element number of "// &
        & TRIM(NumberToVstring(userElementNumber,"*",err,error))// &
        & " does not exist in the domain decomposition.",err,error,*999)
    END IF
    
    basis=>domainTopology%elements%elements(localElementNumber)%basis
    IF(.NOT.ASSOCIATED(basis)) THEN
      localError="The basis for user element number "//TRIM(NumberToVString(userElementNumber,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("DomainTopology_ElementBasisGet")
    RETURN
999 ERRORSEXITS("DomainTopology_ElementBasisGet",err,error)
    RETURN 1
    
  END SUBROUTINE DomainTopology_ElementBasisGet

  !
  !================================================================================================================================
  !

  !>Finalises a line in the given domain topology and deallocates all memory.
  SUBROUTINE DOMAIN_TOPOLOGY_LINE_FINALISE(LINE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_LINE_TYPE) :: LINE !<The domain line to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_LINE_FINALISE",ERR,ERROR,*999)

    LINE%NUMBER=0
    NULLIFY(LINE%BASIS)
    IF(ALLOCATED(LINE%NODES_IN_LINE)) DEALLOCATE(LINE%NODES_IN_LINE)
    IF(ALLOCATED(LINE%DERIVATIVES_IN_LINE)) DEALLOCATE(LINE%DERIVATIVES_IN_LINE)
 
    EXITS("DOMAIN_TOPOLOGY_LINE_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_LINE_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_LINE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the line data structure for a domain topology line.
  SUBROUTINE DOMAIN_TOPOLOGY_LINE_INITIALISE(LINE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_LINE_TYPE) :: LINE !<The domain line to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_LINE_INITIALISE",ERR,ERROR,*999)

    LINE%NUMBER=0
    NULLIFY(LINE%BASIS)
    LINE%BOUNDARY_LINE=.FALSE.
    
    EXITS("DOMAIN_TOPOLOGY_LINE_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_LINE_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_LINE_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the lines in the given domain topology. \todo pass in domain lines
  SUBROUTINE DOMAIN_TOPOLOGY_LINES_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to finalise the lines for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nl
    
    ENTERS("DOMAIN_TOPOLOGY_LINES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%LINES)) THEN
        DO nl=1,TOPOLOGY%LINES%NUMBER_OF_LINES
          CALL DOMAIN_TOPOLOGY_LINE_FINALISE(TOPOLOGY%LINES%LINES(nl),ERR,ERROR,*999)
        ENDDO !nl
        IF(ALLOCATED(TOPOLOGY%LINES%LINES)) DEALLOCATE(TOPOLOGY%LINES%LINES)
        DEALLOCATE(TOPOLOGY%LINES)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_TOPOLOGY_LINES_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_LINES_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_LINES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the line data structures for a domain topology. \todo finalise on error
  SUBROUTINE DOMAIN_TOPOLOGY_LINES_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to initialise the lines for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_LINES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%LINES)) THEN
        CALL FlagError("Decomposition already has topology lines associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%LINES,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology lines",ERR,ERROR,*999)
        TOPOLOGY%LINES%NUMBER_OF_LINES=0
        TOPOLOGY%LINES%DOMAIN=>TOPOLOGY%DOMAIN
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_LINES_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_LINES_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_LINES_INITIALISE
  
  !
  !================================================================================================================================
  !
  !>Finalises a face in the given domain topology and deallocates all memory.
  SUBROUTINE DOMAIN_TOPOLOGY_FACE_FINALISE(FACE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_FACE_TYPE) :: FACE !<The domain face to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_FACE_FINALISE",ERR,ERROR,*999)

    FACE%NUMBER=0
    NULLIFY(FACE%BASIS)
    IF(ALLOCATED(FACE%NODES_IN_FACE)) DEALLOCATE(FACE%NODES_IN_FACE)
    IF(ALLOCATED(FACE%DERIVATIVES_IN_FACE)) DEALLOCATE(FACE%DERIVATIVES_IN_FACE)
 
    EXITS("DOMAIN_TOPOLOGY_FACE_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_FACE_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_FACE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the face data structure for a domain topology face.
  SUBROUTINE DOMAIN_TOPOLOGY_FACE_INITIALISE(FACE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_FACE_TYPE) :: FACE !<The domain face to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_FACE_INITIALISE",ERR,ERROR,*999)

    FACE%NUMBER=0
    NULLIFY(FACE%BASIS)
    FACE%BOUNDARY_FACE=.FALSE.
    
    EXITS("DOMAIN_TOPOLOGY_FACE_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_FACE_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_FACE_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the faces in the given domain topology. \todo pass in domain faces
  SUBROUTINE DOMAIN_TOPOLOGY_FACES_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to finalise the faces for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nf
    
    ENTERS("DOMAIN_TOPOLOGY_FACES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%FACES)) THEN
        DO nf=1,TOPOLOGY%FACES%NUMBER_OF_FACES
          CALL DOMAIN_TOPOLOGY_FACE_FINALISE(TOPOLOGY%FACES%FACES(nf),ERR,ERROR,*999)
        ENDDO !nf
        IF(ALLOCATED(TOPOLOGY%FACES%FACES)) DEALLOCATE(TOPOLOGY%FACES%FACES)
        DEALLOCATE(TOPOLOGY%FACES)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_TOPOLOGY_FACES_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_FACES_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_FACES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the face data structures for a domain topology. \todo finalise on error
  SUBROUTINE DOMAIN_TOPOLOGY_FACES_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to initialise the faces for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_FACES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%FACES)) THEN
        CALL FlagError("Decomposition already has topology faces associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%FACES,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology faces",ERR,ERROR,*999)
        TOPOLOGY%FACES%NUMBER_OF_FACES=0
        TOPOLOGY%FACES%DOMAIN=>TOPOLOGY%DOMAIN
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_FACES_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_FACES_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_FACES_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Checks that a user node number exists in a domain nodes topolgoy. 
  SUBROUTINE DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS(DOMAIN_TOPOLOGY,USER_NODE_NUMBER,NODE_EXISTS,DOMAIN_LOCAL_NODE_NUMBER, &
    & GHOST_NODE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: DOMAIN_TOPOLOGY !<A pointer to the domain topology to check the node exists on
    INTEGER(INTG), INTENT(IN) :: USER_NODE_NUMBER !<The user node number to check if it exists
    LOGICAL, INTENT(OUT) :: NODE_EXISTS !<On exit, is .TRUE. if the node user number exists in the domain nodes topolgoy (even if it is a ghost node), .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: DOMAIN_LOCAL_NODE_NUMBER !<On exit, if the node exists the local number corresponding to the user node number. If the node does not exist then global number will be 0.
    LOGICAL, INTENT(OUT) :: GHOST_NODE !<On exit, is .TRUE. if the local node (if it exists) is a ghost node, .FALSE. if not. 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE
    
    ENTERS("DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS",ERR,ERROR,*999)

    NODE_EXISTS=.FALSE.
    DOMAIN_LOCAL_NODE_NUMBER=0
    GHOST_NODE=.FALSE.
    IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
      DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
      IF(ASSOCIATED(DOMAIN_NODES)) THEN
        NULLIFY(TREE_NODE)
        CALL TREE_SEARCH(DOMAIN_NODES%NODES_TREE,USER_NODE_NUMBER,TREE_NODE,ERR,ERROR,*999)
        IF(ASSOCIATED(TREE_NODE)) THEN
          CALL TREE_NODE_VALUE_GET(DOMAIN_NODES%NODES_TREE,TREE_NODE,DOMAIN_LOCAL_NODE_NUMBER,ERR,ERROR,*999)
          NODE_EXISTS=.TRUE.
          GHOST_NODE=DOMAIN_LOCAL_NODE_NUMBER>DOMAIN_NODES%NUMBER_OF_NODES
        ENDIF
      ELSE
        CALL FlagError("Domain topology nodes is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS
  
  !
  !================================================================================================================================
  !

  !>Finalises the given domain topology node derivative and deallocates all memory.
  SUBROUTINE DOMAIN_TOPOLOGY_NODE_DERIVATIVE_FINALISE(NODE_DERIVATIVE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_NODE_DERIVATIVE_TYPE) :: NODE_DERIVATIVE !<The domain node to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(NODE_DERIVATIVE%userVersionNumbers)) DEALLOCATE(NODE_DERIVATIVE%userVersionNumbers)
    IF(ALLOCATED(NODE_DERIVATIVE%DOF_INDEX)) DEALLOCATE(NODE_DERIVATIVE%DOF_INDEX)

    EXITS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_NODE_DERIVATIVE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the given mesh topology node.
  SUBROUTINE DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE(NODE_DERIVATIVE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_NODE_DERIVATIVE_TYPE) :: NODE_DERIVATIVE !<The domain node to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE",ERR,ERROR,*999)

    NODE_DERIVATIVE%numberOfVersions=0
    NODE_DERIVATIVE%GLOBAL_DERIVATIVE_INDEX=0
    NODE_DERIVATIVE%PARTIAL_DERIVATIVE_INDEX=0

    EXITS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the given domain topology node and deallocates all memory.
  SUBROUTINE DOMAIN_TOPOLOGY_NODE_FINALISE(NODE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_NODE_TYPE) :: NODE !<The domain node to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: derivative_idx

    ENTERS("DOMAIN_TOPOLOGY_NODE_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(NODE%DERIVATIVES)) THEN
      DO derivative_idx=1,NODE%NUMBER_OF_DERIVATIVES
        CALL DOMAIN_TOPOLOGY_NODE_DERIVATIVE_FINALISE(NODE%DERIVATIVES(derivative_idx),ERR,ERROR,*999)
      ENDDO !derivative_idx
      DEALLOCATE(NODE%DERIVATIVES)
    ENDIF
    IF(ASSOCIATED(NODE%SURROUNDING_ELEMENTS)) DEALLOCATE(NODE%SURROUNDING_ELEMENTS)
    IF(ALLOCATED(NODE%NODE_LINES)) DEALLOCATE(NODE%NODE_LINES)
 
    EXITS("DOMAIN_TOPOLOGY_NODE_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_NODE_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_NODE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the given domain topology node.
  SUBROUTINE DOMAIN_TOPOLOGY_NODE_INITIALISE(NODE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_NODE_TYPE) :: NODE !<The domain node to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_NODE_INITIALISE",ERR,ERROR,*999)

    NODE%LOCAL_NUMBER=0
    NODE%MESH_NUMBER=0
    NODE%GLOBAL_NUMBER=0
    NODE%USER_NUMBER=0
    NODE%NUMBER_OF_SURROUNDING_ELEMENTS=0
    NODE%NUMBER_OF_NODE_LINES=0
    NODE%BOUNDARY_NODE=.FALSE.
    
    EXITS("DOMAIN_TOPOLOGY_NODE_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_NODE_INITIALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_NODE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the nodees in the given domain topology. \todo pass in domain nodes
  SUBROUTINE DOMAIN_TOPOLOGY_NODES_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to initialise the nodes for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: np

    ENTERS("DOMAIN_TOPOLOGY_NODES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%NODES)) THEN
        DO np=1,TOPOLOGY%NODES%TOTAL_NUMBER_OF_NODES
          CALL DOMAIN_TOPOLOGY_NODE_FINALISE(TOPOLOGY%NODES%NODES(np),ERR,ERROR,*999)
        ENDDO !np
        IF(ASSOCIATED(TOPOLOGY%NODES%NODES)) DEALLOCATE(TOPOLOGY%NODES%NODES)
        IF(ASSOCIATED(TOPOLOGY%NODES%NODES_TREE)) CALL TREE_DESTROY(TOPOLOGY%NODES%NODES_TREE,ERR,ERROR,*999)
        DEALLOCATE(TOPOLOGY%NODES)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    EXITS("DOMAIN_TOPOLOGY_NODES_FINALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_NODES_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_NODES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the nodes data structures for a domain topology. \todo finalise on error
  SUBROUTINE DOMAIN_TOPOLOGY_NODES_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to initialise the nodes for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DOMAIN_TOPOLOGY_NODES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%NODES)) THEN
        CALL FlagError("Decomposition already has topology nodes associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%NODES,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology nodes",ERR,ERROR,*999)
        TOPOLOGY%NODES%NUMBER_OF_NODES=0
        TOPOLOGY%NODES%TOTAL_NUMBER_OF_NODES=0
        TOPOLOGY%NODES%NUMBER_OF_GLOBAL_NODES=0
        TOPOLOGY%NODES%MAXIMUM_NUMBER_OF_DERIVATIVES=0
        TOPOLOGY%NODES%DOMAIN=>TOPOLOGY%DOMAIN
        NULLIFY(TOPOLOGY%NODES%NODES)
        NULLIFY(TOPOLOGY%NODES%NODES_TREE)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("DOMAIN_TOPOLOGY_NODES_INITIALISE")
    RETURN
999 ERRORSEXITS("DOMAIN_TOPOLOGY_NODES_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_NODES_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Calculates the element numbers surrounding a node for a domain.
  SUBROUTINE DomainTopology_NodesSurroundingElementsCalculate(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to calculate the elements surrounding the nodes for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: element_no,insert_position,ne,nn,np,surrounding_elem_no
    INTEGER(INTG), POINTER :: NEW_SURROUNDING_ELEMENTS(:)
    LOGICAL :: FOUND_ELEMENT
    TYPE(BASIS_TYPE), POINTER :: BASIS

    NULLIFY(NEW_SURROUNDING_ELEMENTS)

    ENTERS("DomainTopology_NodesSurroundingElementsCalculate",ERR,ERROR,*999)
    
    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        IF(ASSOCIATED(TOPOLOGY%NODES)) THEN
          IF(ASSOCIATED(TOPOLOGY%NODES%NODES)) THEN
            DO np=1,TOPOLOGY%NODES%TOTAL_NUMBER_OF_NODES
              TOPOLOGY%NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS=0
              IF(ASSOCIATED(TOPOLOGY%NODES%NODES(np)%SURROUNDING_ELEMENTS)) &
                & DEALLOCATE(TOPOLOGY%NODES%NODES(np)%SURROUNDING_ELEMENTS)
            ENDDO !np
            DO ne=1,TOPOLOGY%ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
              BASIS=>TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS
              DO nn=1,BASIS%NUMBER_OF_NODES
                np=TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn)
                FOUND_ELEMENT=.FALSE.
                element_no=1
                insert_position=1
                DO WHILE(element_no<=TOPOLOGY%NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS.AND..NOT.FOUND_ELEMENT)
                  surrounding_elem_no=TOPOLOGY%NODES%NODES(np)%SURROUNDING_ELEMENTS(element_no)
                  IF(surrounding_elem_no==ne) THEN
                    FOUND_ELEMENT=.TRUE.
                  ENDIF
                  element_no=element_no+1
                  IF(ne>=surrounding_elem_no) THEN
                    insert_position=element_no
                  ENDIF
                ENDDO
                IF(.NOT.FOUND_ELEMENT) THEN
                  !Insert element into surrounding elements
                  ALLOCATE(NEW_SURROUNDING_ELEMENTS(TOPOLOGY%NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS+1),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate new surrounding elements",ERR,ERROR,*999)
                  IF(ASSOCIATED(TOPOLOGY%NODES%NODES(np)%SURROUNDING_ELEMENTS)) THEN
                    NEW_SURROUNDING_ELEMENTS(1:insert_position-1)=TOPOLOGY%NODES%NODES(np)% &
                      & SURROUNDING_ELEMENTS(1:insert_position-1)
                    NEW_SURROUNDING_ELEMENTS(insert_position)=ne
                    NEW_SURROUNDING_ELEMENTS(insert_position+1:TOPOLOGY%NODES%NODES(np)% &
                      & NUMBER_OF_SURROUNDING_ELEMENTS+1)=TOPOLOGY%NODES%NODES(np)%SURROUNDING_ELEMENTS(insert_position: &
                      & TOPOLOGY%NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS)
                    DEALLOCATE(TOPOLOGY%NODES%NODES(np)%SURROUNDING_ELEMENTS)
                  ELSE
                    NEW_SURROUNDING_ELEMENTS(1)=ne
                  ENDIF
                  TOPOLOGY%NODES%NODES(np)%SURROUNDING_ELEMENTS=>NEW_SURROUNDING_ELEMENTS
                  TOPOLOGY%NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS= &
                    & TOPOLOGY%NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS+1
                ENDIF
              ENDDO !nn
            ENDDO !ne
          ELSE
            CALL FlagError("Domain topology nodes nodes are not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Domain topology nodes are not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Domain topology elements is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Domain topology is not associated",ERR,ERROR,*999)
    ENDIF

    EXITS("DomainTopology_NodesSurroundingElementsCalculate")
    RETURN
999 IF(ASSOCIATED(NEW_SURROUNDING_ELEMENTS)) DEALLOCATE(NEW_SURROUNDING_ELEMENTS)
    ERRORS("DomainTopology_NodesSurroundingElementsCalculate",ERR,ERROR)
    EXITS("DomainTopology_NodesSurroundingElementsCalculate")
    RETURN 1   
  END SUBROUTINE DomainTopology_NodesSurroundingElementsCalculate
  
  !
  !================================================================================================================================
  !

  !>Finialises the mesh adjacent elements information and deallocates all memory
  SUBROUTINE MESH_ADJACENT_ELEMENT_FINALISE(MESH_ADJACENT_ELEMENT,ERR,ERROR,*)
    
    !Argument variables
    TYPE(MESH_ADJACENT_ELEMENT_TYPE) :: MESH_ADJACENT_ELEMENT !<The mesh adjacent element to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_ADJACENT_ELEMENT_FINALISE",ERR,ERROR,*999)

    MESH_ADJACENT_ELEMENT%NUMBER_OF_ADJACENT_ELEMENTS=0
    IF(ALLOCATED(MESH_ADJACENT_ELEMENT%ADJACENT_ELEMENTS)) DEALLOCATE(MESH_ADJACENT_ELEMENT%ADJACENT_ELEMENTS)
       
    EXITS("MESH_ADJACENT_ELEMENT_FINALISE")
    RETURN
999 ERRORSEXITS("MESH_ADJACENT_ELEMENT_FINALISE",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_ADJACENT_ELEMENT_FINALISE

  !
  !================================================================================================================================
  !
  !>Initalises the mesh adjacent elements information.
  SUBROUTINE MESH_ADJACENT_ELEMENT_INITIALISE(MESH_ADJACENT_ELEMENT,ERR,ERROR,*)
    
    !Argument variables
    TYPE(MESH_ADJACENT_ELEMENT_TYPE) :: MESH_ADJACENT_ELEMENT !<The mesh adjacent element to initialise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_ADJACENT_ELEMENT_INITIALISE",ERR,ERROR,*999)

    MESH_ADJACENT_ELEMENT%NUMBER_OF_ADJACENT_ELEMENTS=0
       
    EXITS("MESH_ADJACENT_ELEMENT_INITIALISE")
    RETURN
999 ERRORSEXITS("MESH_ADJACENT_ELEMENT_INITIALISE",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_ADJACENT_ELEMENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a mesh. \see OPENCMISS::Iron::cmfe_MeshCreateFinish
  SUBROUTINE MESH_CREATE_FINISH(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to finish creating
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx
    LOGICAL :: FINISHED
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(ASSOCIATED(MESH%TOPOLOGY)) THEN
        !Check that the mesh component elements have been finished
        FINISHED=.TRUE.
        DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
          IF(ASSOCIATED(MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS)) THEN
            IF(.NOT.MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%ELEMENTS_FINISHED) THEN
              LOCAL_ERROR="The elements for mesh component "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                & " have not been finished"
              FINISHED=.FALSE.
              EXIT
            ENDIF
          ELSE
            LOCAL_ERROR="The elements for mesh topology component "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
              & " are not associated"
            FINISHED=.FALSE.
            EXIT
          ENDIF
        ENDDO !component_idx
        IF(.NOT.FINISHED) CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        MESH%MESH_FINISHED=.TRUE.
        !Calulcate the mesh topology
        DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
          CALL MeshTopology_Calculate(MESH%TOPOLOGY(component_idx)%PTR,ERR,ERROR,*999)
        ENDDO !component_idx
      ELSE
        CALL FlagError("Mesh topology is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Mesh user number       = ",MESH%USER_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global number        = ",MESH%GLOBAL_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of dimensions = ",MESH%NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("MESH_CREATE_FINISH",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_CREATE_FINISH
        
  !
  !================================================================================================================================
  !

  !>Starts the process of creating a mesh 
  SUBROUTINE MESH_CREATE_START_GENERIC(MESHES,USER_NUMBER,NUMBER_OF_DIMENSIONS,MESH,ERR,ERROR,*)
    
    !Argument variables
    TYPE(MESHES_TYPE), POINTER :: MESHES !<The pointer to the meshes
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to create
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DIMENSIONS !<The number of dimensions in the mesh.
    TYPE(MESH_TYPE), POINTER :: MESH !<On return, a pointer to the mesh. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,mesh_idx
    TYPE(MESH_TYPE), POINTER :: NEW_MESH
    TYPE(MESH_PTR_TYPE), POINTER :: NEW_MESHES(:)
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    NULLIFY(NEW_MESH)
    NULLIFY(NEW_MESHES)

    ENTERS("MESH_CREATE_START_GENERIC",ERR,ERROR,*997)

    IF(ASSOCIATED(MESHES)) THEN
      IF(ASSOCIATED(MESH)) THEN
        CALL FlagError("Mesh is already associated.",ERR,ERROR,*997)
      ELSE
        CALL MESH_INITIALISE(NEW_MESH,ERR,ERROR,*999)
        !Set default mesh values
        NEW_MESH%USER_NUMBER=USER_NUMBER
        NEW_MESH%GLOBAL_NUMBER=MESHES%NUMBER_OF_MESHES+1
        NEW_MESH%MESHES=>MESHES
        NEW_MESH%NUMBER_OF_DIMENSIONS=NUMBER_OF_DIMENSIONS
        NEW_MESH%NUMBER_OF_COMPONENTS=1
        NEW_MESH%SURROUNDING_ELEMENTS_CALCULATE=.true. !default true
        !Initialise mesh topology and decompositions
        CALL MeshTopology_Initialise(NEW_MESH,ERR,ERROR,*999)
        CALL DECOMPOSITIONS_INITIALISE(NEW_MESH,ERR,ERROR,*999)
        !Add new mesh into list of meshes 
        ALLOCATE(NEW_MESHES(MESHES%NUMBER_OF_MESHES+1),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate new meshes",ERR,ERROR,*999)
        DO mesh_idx=1,MESHES%NUMBER_OF_MESHES
          NEW_MESHES(mesh_idx)%PTR=>MESHES%MESHES(mesh_idx)%PTR
        ENDDO !mesh_idx
        NEW_MESHES(MESHES%NUMBER_OF_MESHES+1)%PTR=>NEW_MESH
        IF(ASSOCIATED(MESHES%MESHES)) DEALLOCATE(MESHES%MESHES)
        MESHES%MESHES=>NEW_MESHES
        MESHES%NUMBER_OF_MESHES=MESHES%NUMBER_OF_MESHES+1
        MESH=>NEW_MESH
      ENDIF
    ELSE
      CALL FlagError("Meshes is not associated.",ERR,ERROR,*997)
    ENDIF
      
    EXITS("MESH_CREATE_START_GENERIC")
    RETURN
999 CALL MESH_FINALISE(NEW_MESH,DUMMY_ERR,DUMMY_ERROR,*998)
998 IF(ASSOCIATED(NEW_MESHES)) DEALLOCATE(NEW_MESHES)
    NULLIFY(MESH)    
997 ERRORSEXITS("MESH_CREATE_START_GENERIC",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_CREATE_START_GENERIC

  !
  !================================================================================================================================
  !

  !>Starts the process of creating a mesh defined by a user number with the specified NUMBER_OF_DIMENSIONS in an interface. \see OPENCMISS::Iron::cmfe_MeshCreateStart
  !>Default values set for the MESH's attributes are:
  !>- NUMBER_OF_COMPONENTS: 1
  SUBROUTINE MESH_CREATE_START_INTERFACE(USER_NUMBER,INTERFACE,NUMBER_OF_DIMENSIONS,MESH,ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to create
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to create the mesh on
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DIMENSIONS !<The number of dimensions in the mesh.
    TYPE(MESH_TYPE), POINTER :: MESH !<On exit, a pointer to the created mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(REGION_TYPE), POINTER :: PARENT_REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_CREATE_START_INTERFACE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(MESH)) THEN
        CALL FlagError("Mesh is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(MESH)
        IF(ASSOCIATED(INTERFACE%MESHES)) THEN
          CALL Mesh_UserNumberFindGeneric(USER_NUMBER,INTERFACE%MESHES,MESH,ERR,ERROR,*999)
          IF(ASSOCIATED(MESH)) THEN
            LOCAL_ERROR="Mesh number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
              & " has already been created on interface number "// &
              & TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            IF(ASSOCIATED(INTERFACE%INTERFACES)) THEN
              PARENT_REGION=>INTERFACE%INTERFACES%PARENT_REGION
              IF(ASSOCIATED(PARENT_REGION)) THEN
                IF(ASSOCIATED(PARENT_REGION%COORDINATE_SYSTEM)) THEN                  
                  IF(NUMBER_OF_DIMENSIONS>0) THEN
                    IF(NUMBER_OF_DIMENSIONS<=PARENT_REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS) THEN
                      CALL MESH_CREATE_START_GENERIC(INTERFACE%MESHES,USER_NUMBER,NUMBER_OF_DIMENSIONS,MESH,ERR,ERROR,*999)
                      MESH%INTERFACE=>INTERFACE
                    ELSE
                      LOCAL_ERROR="Number of mesh dimensions ("//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                        & ") must be <= number of parent region dimensions ("// &
                        & TRIM(NUMBER_TO_VSTRING(PARENT_REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))//")."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Number of mesh dimensions must be > 0.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Parent region coordinate system is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Interfaces parent region is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Interface interfaces is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDIF
        ELSE
          LOCAL_ERROR="The meshes on interface number "//TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))// &
            & " are not associated."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Interface is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_CREATE_START_INTERFACE")
    RETURN
999 ERRORSEXITS("MESH_CREATE_START_INTERFACE",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_CREATE_START_INTERFACE

  !
  !================================================================================================================================
  !

  !>Starts the process of creating a mesh defined by a user number with the specified NUMBER_OF_DIMENSIONS in the region identified by REGION. \see OPENCMISS::Iron::cmfe_MeshCreateStart
  !>Default values set for the MESH's attributes are:
  !>- NUMBER_OF_COMPONENTS: 1
  SUBROUTINE MESH_CREATE_START_REGION(USER_NUMBER,REGION,NUMBER_OF_DIMENSIONS,MESH,ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to create
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to create the mesh on
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DIMENSIONS !<The number of dimensions in the mesh.
    TYPE(MESH_TYPE), POINTER :: MESH !<On exit, a pointer to the created mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_CREATE_START_REGION",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(MESH)) THEN
        CALL FlagError("Mesh is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(MESH)
        IF(ASSOCIATED(REGION%MESHES)) THEN
          CALL Mesh_UserNumberFindGeneric(USER_NUMBER,REGION%MESHES,MESH,ERR,ERROR,*999)
          IF(ASSOCIATED(MESH)) THEN
            LOCAL_ERROR="Mesh number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
              & " has already been created on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            IF(ASSOCIATED(REGION%COORDINATE_SYSTEM)) THEN
              IF(NUMBER_OF_DIMENSIONS>0) THEN
                IF(NUMBER_OF_DIMENSIONS<=REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS) THEN
                  CALL MESH_CREATE_START_GENERIC(REGION%MESHES,USER_NUMBER,NUMBER_OF_DIMENSIONS,MESH,ERR,ERROR,*999)
                  MESH%REGION=>REGION
                ELSE
                  LOCAL_ERROR="Number of mesh dimensions ("//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                    & ") must be <= number of region dimensions ("// &
                    & TRIM(NUMBER_TO_VSTRING(REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))//")."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Number of mesh dimensions must be > 0.",ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The coordinate system on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
                & " are not associated."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ENDIF
        ELSE
          LOCAL_ERROR="The meshes on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
            & " are not associated."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Region is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("MESH_CREATE_START_REGION")
    RETURN
999 ERRORSEXITS("MESH_CREATE_START_REGION",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_CREATE_START_REGION

  !
  !================================================================================================================================
  !

  !>Destroys the mesh identified by a user number on the given region and deallocates all memory. \see OPENCMISS::Iron::cmfe_MeshDestroy
  SUBROUTINE MESH_DESTROY_NUMBER(USER_NUMBER,REGION,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to destroy
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region containing the mesh
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: mesh_idx,mesh_position
    LOGICAL :: FOUND
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(MESH_PTR_TYPE), POINTER :: NEW_MESHES(:)

    NULLIFY(NEW_MESHES)

    ENTERS("MESH_DESTROY_NUMBER",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%MESHES)) THEN

!!TODO: have a mesh_destory_ptr and mesh_destroy_number
        
        !Find the problem identified by the user number
        FOUND=.FALSE.
        mesh_position=0
        DO WHILE(mesh_position<REGION%MESHES%NUMBER_OF_MESHES.AND..NOT.FOUND)
          mesh_position=mesh_position+1
          IF(REGION%MESHES%MESHES(mesh_position)%PTR%USER_NUMBER==USER_NUMBER) FOUND=.TRUE.
        ENDDO
        
        IF(FOUND) THEN
          
          MESH=>REGION%MESHES%MESHES(mesh_position)%PTR

          CALL MESH_FINALISE(MESH,ERR,ERROR,*999)

          !Remove the mesh from the list of meshes
          IF(REGION%MESHES%NUMBER_OF_MESHES>1) THEN
            ALLOCATE(NEW_MESHES(REGION%MESHES%NUMBER_OF_MESHES-1),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate new meshes",ERR,ERROR,*999)
            DO mesh_idx=1,REGION%MESHES%NUMBER_OF_MESHES
              IF(mesh_idx<mesh_position) THEN
                NEW_MESHES(mesh_idx)%PTR=>REGION%MESHES%MESHES(mesh_idx)%PTR
              ELSE IF(mesh_idx>mesh_position) THEN
                REGION%MESHES%MESHES(mesh_idx)%PTR%GLOBAL_NUMBER=REGION%MESHES%MESHES(mesh_idx)%PTR%GLOBAL_NUMBER-1
                NEW_MESHES(mesh_idx-1)%PTR=>REGION%MESHES%MESHES(mesh_idx)%PTR
              ENDIF
            ENDDO !mesh_idx
            DEALLOCATE(REGION%MESHES%MESHES)
            REGION%MESHES%MESHES=>NEW_MESHES
            REGION%MESHES%NUMBER_OF_MESHES=REGION%MESHES%NUMBER_OF_MESHES-1
          ELSE
            DEALLOCATE(REGION%MESHES%MESHES)
            REGION%MESHES%NUMBER_OF_MESHES=0
          ENDIF
          
        ELSE
          LOCAL_ERROR="Mesh number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
            & " has not been created on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The meshes on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
          & " are not associated"
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Region is not associated",ERR,ERROR,*999)
    ENDIF

    EXITS("MESH_DESTROY_NUMBER")
    RETURN
999 IF(ASSOCIATED(NEW_MESHES)) DEALLOCATE(NEW_MESHES)
    ERRORSEXITS("MESH_DESTROY_NUMBER",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE MESH_DESTROY_NUMBER

  !
  !================================================================================================================================
  !

  !>Destroys the mesh and deallocates all memory. \see OPENCMISS::Iron::cmfe_MeshDestroy
  SUBROUTINE MESH_DESTROY(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to destroy.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: mesh_idx,mesh_position
    TYPE(MESHES_TYPE), POINTER :: MESHES
    TYPE(MESH_PTR_TYPE), POINTER :: NEW_MESHES(:)

    NULLIFY(NEW_MESHES)

    ENTERS("MESH_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      MESHES=>MESH%MESHES
      IF(ASSOCIATED(MESHES)) THEN
        mesh_position=MESH%GLOBAL_NUMBER
          
        CALL MESH_FINALISE(MESH,ERR,ERROR,*999)

        !Remove the mesh from the list of meshes
        IF(MESHES%NUMBER_OF_MESHES>1) THEN
          ALLOCATE(NEW_MESHES(MESHES%NUMBER_OF_MESHES-1),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate new meshes.",ERR,ERROR,*999)
          DO mesh_idx=1,MESHES%NUMBER_OF_MESHES
            IF(mesh_idx<mesh_position) THEN
              NEW_MESHES(mesh_idx)%PTR=>MESHES%MESHES(mesh_idx)%PTR
            ELSE IF(mesh_idx>mesh_position) THEN
              MESHES%MESHES(mesh_idx)%PTR%GLOBAL_NUMBER=MESHES%MESHES(mesh_idx)%PTR%GLOBAL_NUMBER-1
              NEW_MESHES(mesh_idx-1)%PTR=>MESHES%MESHES(mesh_idx)%PTR
            ENDIF
          ENDDO !mesh_idx
          DEALLOCATE(MESHES%MESHES)
          MESHES%MESHES=>NEW_MESHES
          MESHES%NUMBER_OF_MESHES=MESHES%NUMBER_OF_MESHES-1
        ELSE
          DEALLOCATE(MESHES%MESHES)
          MESHES%NUMBER_OF_MESHES=0
        ENDIF
      ELSE
        CALL FlagError("The mesh meshes is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("MESH_DESTROY")
    RETURN
999 IF(ASSOCIATED(NEW_MESHES)) DEALLOCATE(NEW_MESHES)
    ERRORSEXITS("MESH_DESTROY",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE MESH_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalises a mesh and deallocates all memory.
  SUBROUTINE MESH_FINALISE(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      CALL MeshTopology_Finalise(MESH,ERR,ERROR,*999)
      CALL DECOMPOSITIONS_FINALISE(MESH,ERR,ERROR,*999)
!      IF(ASSOCIATED(MESH%INTF)) CALL INTERFACE_MESH_FINALISE(MESH,ERR,ERROR,*999)  ! <<??>>
      DEALLOCATE(MESH)
    ENDIF
 
    EXITS("MESH_FINALISE")
    RETURN
999 ERRORSEXITS("MESH_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE MESH_FINALISE

  !
  !================================================================================================================================
  !

  !>Returns a region nodes pointer corresponding to the mesh global nodes accounting for interfaces.
  SUBROUTINE MeshGlobalNodesGet(mesh,nodes,err,error,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: mesh !<A pointer to the mesh to get the global nodes for
    TYPE(NODES_TYPE), POINTER :: nodes !<On return, the nodes pointer corresponding to the global nodes for the mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(INTERFACE_TYPE), POINTER :: interface
    TYPE(REGION_TYPE), POINTER :: region
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("MeshGlobalNodesGet",err,error,*999)
    
    IF(ASSOCIATED(mesh)) THEN
      IF(ASSOCIATED(nodes)) THEN
        CALL FlagError("Nodes is already associated.",err,error,*999)
      ELSE
        NULLIFY(nodes)
        region=>mesh%region
        IF(ASSOCIATED(region)) THEN
          nodes=>region%nodes
        ELSE
          INTERFACE=>mesh%INTERFACE
          IF(ASSOCIATED(interface)) THEN
            nodes=>interface%nodes
          ELSE
            localError="Mesh number "//TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))// &
              & " does not have an associated region or interface."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDIF
        IF(.NOT.ASSOCIATED(nodes)) THEN
          IF(ASSOCIATED(region)) THEN
            localError="Mesh number "//TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))// &
              & " does not have any nodes associated with the mesh region."
          ELSE
            localError="Mesh number "//TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))// &
              & " does not have any nodes associated with the mesh interface."
          ENDIF
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",err,error,*999)
    ENDIF
        
    EXITS("MeshGlobalNodesGet")
    RETURN
999 ERRORSEXITS("MeshGlobalNodesGet",err,error)
    RETURN 1
    
  END SUBROUTINE MeshGlobalNodesGet
  
  !
  !================================================================================================================================
  !

  !>Initialises a mesh.
  SUBROUTINE MESH_INITIALISE(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      CALL FlagError("Mesh is already associated.",ERR,ERROR,*999)
    ELSE
      ALLOCATE(MESH,STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate new mesh.",ERR,ERROR,*999)
      MESH%USER_NUMBER=0
      MESH%GLOBAL_NUMBER=0
      MESH%MESH_FINISHED=.FALSE.
      NULLIFY(MESH%MESHES)
      NULLIFY(MESH%REGION)
      NULLIFY(MESH%INTERFACE)
      NULLIFY(MESH%GENERATED_MESH)
      MESH%NUMBER_OF_DIMENSIONS=0
      MESH%NUMBER_OF_COMPONENTS=0
      MESH%MESH_EMBEDDED=.FALSE.
      NULLIFY(MESH%EMBEDDING_MESH)
      MESH%NUMBER_OF_EMBEDDED_MESHES=0
      NULLIFY(MESH%EMBEDDED_MESHES)
      MESH%NUMBER_OF_ELEMENTS=0
      NULLIFY(MESH%TOPOLOGY)
      NULLIFY(MESH%LOCAL_TOPOLOGY)
      NULLIFY(MESH%DECOMPOSITIONS)
    ENDIF
    
    EXITS("MESH_INITIALISE")
    RETURN
999 ERRORSEXITS("MESH_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_INITIALISE

  !
  !================================================================================================================================
  !
  
  !>Gets the number of mesh components for a mesh identified by a pointer. \see OPENCMISS::Iron::cmfe_MeshNumberOfComponentsGet
  SUBROUTINE MESH_NUMBER_OF_COMPONENTS_GET(MESH,NUMBER_OF_COMPONENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to get the number of components for
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_COMPONENTS !<On return, the number of components in the specified mesh.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("MESH_NUMBER_OF_COMPONENTS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(MESH%MESH_FINISHED) THEN
        NUMBER_OF_COMPONENTS=MESH%NUMBER_OF_COMPONENTS
      ELSE
        CALL FlagError("Mesh has not finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_NUMBER_OF_COMPONENTS_GET")
    RETURN
999 ERRORSEXITS("MESH_NUMBER_OF_COMPONENTS_GET",ERR,ERROR)    
    RETURN
  END SUBROUTINE MESH_NUMBER_OF_COMPONENTS_GET

  !
  !================================================================================================================================
  !

  !>Changes/sets the number of mesh components for a mesh. \see OPENCMISS::Iron::cmfe_MeshNumberOfComponentsSet
  SUBROUTINE MESH_NUMBER_OF_COMPONENTS_SET(MESH,NUMBER_OF_COMPONENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to set the number of components for
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_COMPONENTS !<The number of components to set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(MeshComponentTopologyPtrType), POINTER :: NEW_TOPOLOGY(:)

    NULLIFY(NEW_TOPOLOGY)
    
    ENTERS("MESH_NUMBER_OF_COMPONENTS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(NUMBER_OF_COMPONENTS>0) THEN
        IF(MESH%MESH_FINISHED) THEN
          CALL FlagError("Mesh has been finished",ERR,ERROR,*999)
        ELSE
          IF(NUMBER_OF_COMPONENTS/=MESH%NUMBER_OF_COMPONENTS) THEN
            ALLOCATE(NEW_TOPOLOGY(NUMBER_OF_COMPONENTS),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate new topology",ERR,ERROR,*999)
            IF(NUMBER_OF_COMPONENTS<MESH%NUMBER_OF_COMPONENTS) THEN
              DO component_idx=1,NUMBER_OF_COMPONENTS
                NEW_TOPOLOGY(component_idx)%PTR=>MESH%TOPOLOGY(component_idx)%PTR
              ENDDO !component_idx
            ELSE !NUMBER_OF_COMPONENTS>MESH%NUMBER_OF_COMPONENTS
              DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
                NEW_TOPOLOGY(component_idx)%PTR=>MESH%TOPOLOGY(component_idx)%PTR
              ENDDO !component_idx
!!TODO \todo sort out mesh_topology initialise/finalise so that they allocate and deal with this below then call that routine
              DO component_idx=MESH%NUMBER_OF_COMPONENTS+1,NUMBER_OF_COMPONENTS
                ALLOCATE(NEW_TOPOLOGY(component_idx)%PTR,STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate new topology component",ERR,ERROR,*999)
                NEW_TOPOLOGY(component_idx)%PTR%mesh=>mesh
                NEW_TOPOLOGY(component_idx)%PTR%meshComponentNumber=component_idx
                NULLIFY(NEW_TOPOLOGY(component_idx)%PTR%elements)
                NULLIFY(NEW_TOPOLOGY(component_idx)%PTR%nodes)
                NULLIFY(NEW_TOPOLOGY(component_idx)%PTR%dofs)
                NULLIFY(NEW_TOPOLOGY(component_idx)%PTR%dataPoints)
                !Initialise the topology components
                CALL MESH_TOPOLOGY_ELEMENTS_INITIALISE(NEW_TOPOLOGY(component_idx)%PTR,ERR,ERROR,*999)
                CALL MeshTopology_NodesInitialise(NEW_TOPOLOGY(component_idx)%PTR,ERR,ERROR,*999)
                CALL MeshTopology_DofsInitialise(NEW_TOPOLOGY(component_idx)%PTR,ERR,ERROR,*999)
                CALL MESH_TOPOLOGY_DATA_POINTS_INITIALISE(NEW_TOPOLOGY(component_idx)%PTR,ERR,ERROR,*999)
              ENDDO !component_idx
            ENDIF
            IF(ASSOCIATED(MESH%TOPOLOGY)) DEALLOCATE(MESH%TOPOLOGY)
            MESH%TOPOLOGY=>NEW_TOPOLOGY
            MESH%NUMBER_OF_COMPONENTS=NUMBER_OF_COMPONENTS
          ENDIF
        ENDIF
      ELSE
        LOCAL_ERROR="The specified number of mesh components ("//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
          & ") is illegal. You must have >0 mesh components"
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF      
    ELSE
      CALL FlagError("Mesh is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_NUMBER_OF_COMPONENTS_SET")
    RETURN
!!TODO: tidy up memory deallocation on error
999 ERRORSEXITS("MESH_NUMBER_OF_COMPONENTS_SET",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_NUMBER_OF_COMPONENTS_SET

  !
  !================================================================================================================================
  !
  
  !>Gets the number of elements for a mesh identified by a pointer. \see OPENCMISS::Iron::cmfe_MeshNumberOfElementsGet
  SUBROUTINE MESH_NUMBER_OF_ELEMENTS_GET(MESH,NUMBER_OF_ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to get the number of elements for
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_ELEMENTS !<On return, the number of elements in the specified mesh
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_NUMBER_OF_ELEMENTS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(MESH%MESH_FINISHED) THEN
        NUMBER_OF_ELEMENTS=MESH%NUMBER_OF_ELEMENTS
      ELSE
        CALL FlagError("Mesh has not been finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_NUMBER_OF_ELEMENTS_GET")
    RETURN
999 ERRORSEXITS("MESH_NUMBER_OF_ELEMENTS_GET",ERR,ERROR)    
    RETURN
  END SUBROUTINE MESH_NUMBER_OF_ELEMENTS_GET

  !
  !================================================================================================================================
  !

  !>Changes/sets the number of elements for a mesh. \see OPENCMISS::Iron::cmfe_MeshNumberOfElementsSet
  SUBROUTINE MESH_NUMBER_OF_ELEMENTS_SET(MESH,NUMBER_OF_ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to set the number of elements for
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_ELEMENTS !<The number of elements to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_NUMBER_OF_ELEMENTS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(NUMBER_OF_ELEMENTS>0) THEN
        IF(MESH%MESH_FINISHED) THEN
          CALL FlagError("Mesh has been finished.",ERR,ERROR,*999)
        ELSE
          IF(NUMBER_OF_ELEMENTS/=MESH%NUMBER_OF_ELEMENTS) THEN
            IF(ASSOCIATED(MESH%TOPOLOGY)) THEN
              DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
                IF(ASSOCIATED(MESH%TOPOLOGY(component_idx)%PTR)) THEN
                  IF(ASSOCIATED(MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS)) THEN
                    IF(MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%NUMBER_OF_ELEMENTS>0) THEN
!!TODO: Reallocate the elements and copy information. 
                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                    ENDIF
                  ENDIF
                ELSE
                  CALL FlagError("Mesh topology component pointer is not associated.",ERR,ERROR,*999)
                ENDIF
              ENDDO !component_idx
            ELSE
              CALL FlagError("Mesh topology is not associated.",ERR,ERROR,*999)
            ENDIF
            MESH%NUMBER_OF_ELEMENTS=NUMBER_OF_ELEMENTS
          ENDIF
        ENDIF
      ELSE
        LOCAL_ERROR="The specified number of elements ("//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_ELEMENTS,"*",ERR,ERROR))// &
          & ") is invalid. You must have > 0 elements."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF      
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_NUMBER_OF_ELEMENTS_SET")
    RETURN
999 ERRORSEXITS("MESH_NUMBER_OF_ELEMENTS_SET",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_NUMBER_OF_ELEMENTS_SET

  !
  !================================================================================================================================
  !

  !>Changes/sets the surrounding elements calculate flag. \see OPENCMISS::Iron::cmfe_MeshSurroundingElementsCalculateSet
  SUBROUTINE MESH_SURROUNDING_ELEMENTS_CALCULATE_SET(MESH,SURROUNDING_ELEMENTS_CALCULATE_FLAG,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to set the surrounding elements calculate flag for
    LOGICAL, INTENT(IN) :: SURROUNDING_ELEMENTS_CALCULATE_FLAG !<The surrounding elements calculate flag
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    ENTERS("MESH_SURROUNDING_ELEMENTS_CALCULATE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(MESH%MESH_FINISHED) THEN
        CALL FlagError("Mesh has been finished.",ERR,ERROR,*999)
      ELSE
        MESH%SURROUNDING_ELEMENTS_CALCULATE=SURROUNDING_ELEMENTS_CALCULATE_FLAG
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_SURROUNDING_ELEMENTS_CALCULATE_SET")
    RETURN
999 ERRORSEXITS("MESH_SURROUNDING_ELEMENTS_CALCULATE_SET",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_SURROUNDING_ELEMENTS_CALCULATE_SET

  !
  !================================================================================================================================
  !

  !>Calculates the mesh topology.
  SUBROUTINE MeshTopology_Calculate(topology,err,error,*)    

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: topology !<A pointer to the mesh topology to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopology_Calculate",err,error,*999)

    IF(ASSOCIATED(topology)) THEN
      !Calculate the nodes used in the mesh
      CALL MeshTopology_NodesCalculate(topology,err,error,*999)
      !Calculate the elements surrounding the nodes in a mesh
      CALL MeshTopology_SurroundingElementsCalculate(topology,err,error,*999)
      !Calculate the number of derivatives at each node in a mesh
      CALL MeshTopology_NodesDerivativesCalculate(topology,err,error,*999)
      !Calculate the number of versions for each derivative at each node in a mesh
      CALL MeshTopology_NodesVersionCalculate(topology,err,error,*999)
      !Calculate the elements surrounding the elements in the mesh
      CALL MeshTopology_ElementsAdjacentElementsCalculate(topology,err,error,*999)
      !Calculate the boundary nodes and elements in the mesh
      CALL MeshTopology_BoundaryCalculate(topology,err,error,*999)
      !Calculate the elements surrounding the elements in the mesh
      CALL MeshTopology_DofsCalculate(topology,err,error,*999)
    ELSE
      CALL FlagError("Topology is not associated",err,error,*999)
    ENDIF
    
    EXITS("MeshTopology_Calculate")
    RETURN
999 ERRORSEXITS("MeshTopology_Calculate",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_Calculate
  
  !
  !===============================================================================================================================
  !

  !>Calculates the boundary nodes and elements for a mesh topology. 
  SUBROUTINE MeshTopology_BoundaryCalculate(topology,err,error,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: topology !<A pointer to the mesh topology to calculate the boundary for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx,localNodeIdx,matchIndex,nodeIdx,xiCoordIdx,xiDirection
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(MeshElementsType), POINTER :: elements
    TYPE(MeshNodesType), POINTER :: nodes
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("MeshTopology_BoundaryCalculate",err,error,*999)

    IF(ASSOCIATED(topology)) THEN
      nodes=>topology%nodes
      IF(ASSOCIATED(nodes)) THEN
        elements=>topology%elements
        IF(ASSOCIATED(elements)) THEN
          DO elementIdx=1,elements%NUMBER_OF_ELEMENTS
            basis=>elements%elements(elementIdx)%basis
            SELECT CASE(basis%type)
            CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
              DO xiCoordIdx=-basis%NUMBER_OF_XI_COORDINATES,basis%NUMBER_OF_XI_COORDINATES
                IF(xiCoordIdx/=0) THEN
                  IF(elements%elements(elementIdx)%ADJACENT_ELEMENTS(xiCoordIdx)%NUMBER_OF_ADJACENT_ELEMENTS==0) THEN
                    elements%elements(elementIdx)%BOUNDARY_ELEMENT=.TRUE.
                    IF(xiCoordIdx<0) THEN
                      xiDirection=-xiCoordIdx
                      matchIndex=1
                    ELSE
                      xiDirection=xiCoordIdx
                      matchIndex=BASIS%NUMBER_OF_NODES_XIC(xiCoordIdx)
                    ENDIF                    
                    DO localNodeIdx=1,BASIS%NUMBER_OF_NODES
                      IF(basis%NODE_POSITION_INDEX(localNodeIdx,XIDIRECTION)==matchIndex) THEN
                        nodeIdx=elements%elements(elementIdx)%MESH_ELEMENT_NODES(localNodeIdx)
                        nodes%nodes(nodeIdx)%boundaryNode=.TRUE.
                      ENDIF
                    ENDDO !nn
                  ENDIF
                ENDIF
              ENDDO !xiCoordIdx            
            CASE(BASIS_SIMPLEX_TYPE)
              elements%elements(elementIdx)%BOUNDARY_ELEMENT=.FALSE.
              DO xiCoordIdx=1,basis%NUMBER_OF_XI_COORDINATES
                elements%elements(elementIdx)%BOUNDARY_ELEMENT=elements%elements(elementIdx)%BOUNDARY_ELEMENT.OR. &
                  & elements%elements(elementIdx)%ADJACENT_ELEMENTS(xiCoordIdx)%NUMBER_OF_ADJACENT_ELEMENTS==0
                IF(elements%elements(elementIdx)%ADJACENT_ELEMENTS(xiCoordIdx)%NUMBER_OF_ADJACENT_ELEMENTS==0) THEN
                  DO localNodeIdx=1,basis%NUMBER_OF_NODES
                    IF(basis%NODE_POSITION_INDEX(localNodeIdx,xiCoordIdx)==1) THEN
                      nodeIdx=elements%elements(elementIdx)%MESH_ELEMENT_NODES(localNodeIdx)
                      nodes%nodes(nodeIdx)%boundaryNode=.TRUE.
                    ENDIF
                  ENDDO !localNodeIdx                  
                ENDIF
              ENDDO !xiCoordIdx
            CASE(BASIS_SERENDIPITY_TYPE)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(BASIS_AUXILLIARY_TYPE)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(BASIS_B_SPLINE_TP_TYPE)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(BASIS_FOURIER_LAGRANGE_HERMITE_TP_TYPE)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(BASIS_EXTENDED_LAGRANGE_TP_TYPE)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE DEFAULT
              localError="The basis type of "//TRIM(NumberToVString(basis%TYPE,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDDO !elementIdx
        ELSE
          CALL FlagError("Topology elements is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Topology nodes is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",err,error,*999)
    ENDIF

    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Boundary elements:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number of elements = ",elements%NUMBER_OF_ELEMENTS,err,error,*999)
      DO elementIdx=1,elements%NUMBER_OF_ELEMENTS
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Element : ",elementIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary element = ",elements%elements(elementIdx)%BOUNDARY_ELEMENT, &
          & err,error,*999)        
      ENDDO !elementIdx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Boundary nodes:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number of nodes = ",nodes%numberOfNodes,err,error,*999)
      DO nodeIdx=1,nodes%numberOfNodes
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Node : ",nodeIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary node = ",nodes%nodes(nodeIdx)%boundaryNode,err,error,*999)        
      ENDDO !elementIdx            
    ENDIF
 
    EXITS("MeshTopology_BoundaryCalculate")
    RETURN
999 ERRORSEXITS("MeshTopology_BoundaryCalculate",err,error)
    RETURN 1
  END SUBROUTINE MeshTopology_BoundaryCalculate

  !
  !===============================================================================================================================
  !

  !>Calculates the degrees-of-freedom for a mesh topology. 
  SUBROUTINE MeshTopology_DofsCalculate(topology,err,error,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: topology !<A pointer to the mesh topology to calculate the dofs for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: derivativeIdx,nodeIdx,numberOfDofs,versionIdx
    TYPE(MeshDofsType), POINTER :: dofs
    TYPE(MeshNodesType), POINTER :: nodes

    ENTERS("MeshTopology_DofsCalculate",err,error,*999)

    IF(ASSOCIATED(topology)) THEN
      nodes=>topology%nodes
      IF(ASSOCIATED(nodes)) THEN
        dofs=>topology%dofs
        IF(ASSOCIATED(dofs)) THEN
          numberOfDofs=0
          DO nodeIdx=1,nodes%numberOfNodes
            DO derivativeIdx=1,nodes%nodes(nodeIdx)%numberOfDerivatives
              ALLOCATE(nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%dofIndex( &
                & nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate mesh topology node derivative version dof index.",err,error,*999)
              DO versionIdx=1,nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions
                numberOfDofs=numberOfDofs+1
                nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%dofIndex(versionIdx)=numberOfDofs
              ENDDO !versionIdx
            ENDDO !derivativeIdx
          ENDDO !nodeIdx
          dofs%numberOfDofs=numberOfDofs
        ELSE
          CALL FlagError("Topology dofs is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Topology nodes is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",err,error,*999)
    ENDIF
 
    EXITS("MeshTopology_DofsCalculate")
    RETURN
999 ERRORSEXITS("MeshTopology_DofsCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_DofsCalculate

  !
  !===============================================================================================================================
  !

  !>Finalises the dof data structures for a mesh topology and deallocates any memory. \todo pass in dofs
  SUBROUTINE MeshTopology_DofsFinalise(dofs,err,error,*)

    !Argument variables
    TYPE(MeshDofsType), POINTER :: dofs !<A pointer to the mesh topology to finalise the dofs for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopology_DofsFinalise",err,error,*999)

    IF(ASSOCIATED(dofs)) THEN
      DEALLOCATE(dofs)
    ENDIF
 
    EXITS("MeshTopology_DofsFinalise")
    RETURN
999 ERRORSEXITS("MeshTopology_DofsFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_DofsFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the dofs in a given mesh topology.
  SUBROUTINE MeshTopology_DofsInitialise(topology,err,error,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: topology !<A pointer to the mesh topology to initialise the dofs for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("MeshTopology_DofsInitialise",err,error,*998)

    IF(ASSOCIATED(topology)) THEN
      IF(ASSOCIATED(topology%dofs)) THEN
        CALL FlagError("Mesh already has topology dofs associated",err,error,*998)
      ELSE
        ALLOCATE(topology%dofs,STAT=err)
        IF(ERR/=0) CALL FlagError("Could not allocate topology dofs",err,error,*999)
        topology%dofs%numberOfDofs=0
        topology%dofs%meshComponentTopology=>topology
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",err,error,*998)
    ENDIF
    
    EXITS("MeshTopology_DofsInitialise")
    RETURN
999 CALL MeshTopology_DofsFinalise(topology%dofs,dummyErr,dummyError,*998)
998 ERRORSEXITS("MeshTopology_DofsInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_DofsInitialise

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating elements for a specified mesh component in a mesh topology. \see OPENCMISS::Iron::cmfe_MeshElementsCreateFinish
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the mesh elements to finish creating
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ne
    TYPE(MESH_TYPE), POINTER :: mesh
    TYPE(MeshComponentTopologyType), POINTER :: meshComponentTopology

    ENTERS("MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FlagError("Mesh elements have already been finished.",ERR,ERROR,*999)
      ELSE        
        ELEMENTS%ELEMENTS_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FlagError("Mesh elements is not associated.",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      meshComponentTopology=>elements%meshComponentTopology
      IF(ASSOCIATED(meshComponentTopology)) THEN
        mesh=>meshComponentTopology%mesh
        IF(ASSOCIATED(MESH)) THEN
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Number of global elements = ",MESH%NUMBER_OF_ELEMENTS, &
            & ERR,ERROR,*999)
          DO ne=1,MESH%NUMBER_OF_ELEMENTS
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Element = ",ne,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global number        = ",ELEMENTS%ELEMENTS(ne)%GLOBAL_NUMBER, &
              & ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    User number          = ",ELEMENTS%ELEMENTS(ne)%USER_NUMBER, &
              & ERR,ERROR,*999)
            IF(ASSOCIATED(ELEMENTS%ELEMENTS(ne)%BASIS)) THEN
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Basis number         = ",ELEMENTS%ELEMENTS(ne)%BASIS% &
                & USER_NUMBER,ERR,ERROR,*999)
            ELSE
              CALL FlagError("Basis is not associated.",ERR,ERROR,*999)
            ENDIF
            IF(ALLOCATED(ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODES)) THEN
              CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS%ELEMENTS(ne)% BASIS%NUMBER_OF_NODES,8,8, &
                & ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODES,'("    User element nodes   =",8(X,I6))','(26X,8(X,I6))', &
                & ERR,ERROR,*999)
            ELSE
              CALL FlagError("User element nodes are not associated.",ERR,ERROR,*999)
            ENDIF
            IF(ALLOCATED(ELEMENTS%ELEMENTS(ne)%GLOBAL_ELEMENT_NODES)) THEN
              CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS%ELEMENTS(ne)%BASIS%NUMBER_OF_NODES,8,8, &
                & ELEMENTS%ELEMENTS(ne)%GLOBAL_ELEMENT_NODES,'("    Global element nodes =",8(X,I6))','(26X,8(X,I6))', &
                & ERR,ERROR,*999)
            ELSE
              CALL FlagError("Global element nodes are not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !ne
        ELSE
          CALL FlagError("Mesh component topology mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Mesh elements mesh component topology is not associated.",ERR,ERROR,*999)
      ENDIF
    ENDIF
    
    EXITS("MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH
    
  !
  !================================================================================================================================
  !

  !>Starts the process of creating elements in the mesh component identified by MESH and component_idx. The elements will be created with a default basis of BASIS. ELEMENTS is the returned pointer to the MESH_ELEMENTS data structure. \see OPENCMISS::Iron::cmfe_MeshElementsCreateStart
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_CREATE_START(MESH,MESH_COMPONENT_NUMBER,BASIS,ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to start creating the elements on
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT_NUMBER !<The mesh component number
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the default basis to use
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<On return, a pointer to the created mesh elements
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,INSERT_STATUS,ne
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
 
    ENTERS("MESH_TOPOLOGY_ELEMENTS_CREATE_START",ERR,ERROR,*999)
    
    IF(ASSOCIATED(MESH)) THEN     
      IF(MESH_COMPONENT_NUMBER>0.AND.MESH_COMPONENT_NUMBER<=MESH%NUMBER_OF_COMPONENTS) THEN
        IF(ASSOCIATED(ELEMENTS)) THEN
          CALL FlagError("Elements is already associated.",ERR,ERROR,*999)
        ELSE
          IF(ASSOCIATED(MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR)) THEN
            IF(ASSOCIATED(MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%ELEMENTS)) THEN
              ELEMENTS=>MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%ELEMENTS
              IF(ASSOCIATED(ELEMENTS%ELEMENTS)) THEN
                CALL FlagError("Mesh topology already has elements associated",ERR,ERROR,*998)
              ELSE
                IF(ASSOCIATED(BASIS)) THEN
                  MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%meshComponentNumber=MESH_COMPONENT_NUMBER
                  ALLOCATE(ELEMENTS%ELEMENTS(MESH%NUMBER_OF_ELEMENTS),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate individual elements",ERR,ERROR,*999)
                  ELEMENTS%NUMBER_OF_ELEMENTS=MESH%NUMBER_OF_ELEMENTS !Psuedo inheritance of the number of elements
                  CALL TREE_CREATE_START(ELEMENTS%ELEMENTS_TREE,ERR,ERROR,*999)
                  CALL TREE_INSERT_TYPE_SET(ELEMENTS%ELEMENTS_TREE,TREE_NO_DUPLICATES_ALLOWED,ERR,ERROR,*999)
                  CALL TREE_CREATE_FINISH(ELEMENTS%ELEMENTS_TREE,ERR,ERROR,*999)
                  ELEMENTS%ELEMENTS_FINISHED=.FALSE.
                  !Set up the default values and allocate element structures
                  DO ne=1,ELEMENTS%NUMBER_OF_ELEMENTS
                    CALL MESH_TOPOLOGY_ELEMENT_INITIALISE(ELEMENTS%ELEMENTS(ne),ERR,ERROR,*999)
                    ELEMENTS%ELEMENTS(ne)%GLOBAL_NUMBER=ne
                    ELEMENTS%ELEMENTS(ne)%USER_NUMBER=ne
                    CALL TREE_ITEM_INSERT(ELEMENTS%ELEMENTS_TREE,ne,ne,INSERT_STATUS,ERR,ERROR,*999)
                    ELEMENTS%ELEMENTS(ne)%BASIS=>BASIS
                    ALLOCATE(ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODES(BASIS%NUMBER_OF_NODES),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate user element nodes",ERR,ERROR,*999)
                    ALLOCATE(ELEMENTS%ELEMENTS(ne)%GLOBAL_ELEMENT_NODES(BASIS%NUMBER_OF_NODES),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate global element nodes",ERR,ERROR,*999)
                    ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODES=1
                    ELEMENTS%ELEMENTS(ne)%GLOBAL_ELEMENT_NODES=1
                    ALLOCATE(ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODE_VERSIONS(BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                      & BASIS%NUMBER_OF_NODES),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate global element nodes versions",ERR,ERROR,*999)
                    ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODE_VERSIONS = 1
                  ENDDO !ne
                ELSE
                  CALL FlagError("Basis is not associated",ERR,ERROR,*999)
                ENDIF
              ENDIF
            ELSE
              CALL FlagError("Mesh topology elements is not associated",ERR,ERROR,*998)
            ENDIF
          ELSE
            CALL FlagError("Mesh topology is not associated",ERR,ERROR,*998)
          ENDIF
        ENDIF
      ELSE
        LOCAL_ERROR="The specified mesh component number of "//TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT_NUMBER,"*",ERR,ERROR))// &
          & " is invalid. The component number must be between 1 and "// &
          & TRIM(NUMBER_TO_VSTRING(MESH%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated",ERR,ERROR,*998)
    ENDIF

    EXITS("MESH_TOPOLOGY_ELEMENTS_CREATE_START")
    RETURN
999 CALL MESH_TOPOLOGY_ELEMENTS_FINALISE(ELEMENTS,DUMMY_ERR,DUMMY_ERROR,*998)
998 NULLIFY(ELEMENTS)
    ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_CREATE_START",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroys the elements in a mesh topology. \todo as this is a user routine it should take a mesh pointer like create start and finish? Split this into destroy and finalise?
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_DESTROY(ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the mesh elements to destroy 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_TOPOLOGY_ELEMENTS_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      CALL MESH_TOPOLOGY_ELEMENTS_FINALISE(ELEMENTS,ERR,ERROR,*999)
    ELSE
      CALL FlagError("Mesh topology is not associated",ERR,ERROR,*999)
    ENDIF

    EXITS("MESH_TOPOLOGY_ELEMENTS_DESTROY")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_DESTROY",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_DESTROY
  
  !
  !================================================================================================================================
  !

  !>Finalises the given mesh topology element.
  SUBROUTINE MESH_TOPOLOGY_ELEMENT_FINALISE(ELEMENT,ERR,ERROR,*)
    
    !Argument variables
    TYPE(MESH_ELEMENT_TYPE) :: ELEMENT !<The mesh element to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nic
    
    ENTERS("MESH_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR,*999)
    
    IF(ALLOCATED(ELEMENT%USER_ELEMENT_NODE_VERSIONS)) DEALLOCATE(ELEMENT%USER_ELEMENT_NODE_VERSIONS)
    IF(ALLOCATED(ELEMENT%USER_ELEMENT_NODES)) DEALLOCATE(ELEMENT%USER_ELEMENT_NODES)
    IF(ALLOCATED(ELEMENT%GLOBAL_ELEMENT_NODES)) DEALLOCATE(ELEMENT%GLOBAL_ELEMENT_NODES)
    IF(ALLOCATED(ELEMENT%MESH_ELEMENT_NODES)) DEALLOCATE(ELEMENT%MESH_ELEMENT_NODES)
    IF(ALLOCATED(ELEMENT%ADJACENT_ELEMENTS)) THEN
      DO nic=LBOUND(ELEMENT%ADJACENT_ELEMENTS,1),UBOUND(ELEMENT%ADJACENT_ELEMENTS,1)
        CALL MESH_ADJACENT_ELEMENT_FINALISE(ELEMENT%ADJACENT_ELEMENTS(nic),ERR,ERROR,*999)
      ENDDO !nic
      DEALLOCATE(ELEMENT%ADJACENT_ELEMENTS)
    ENDIF
    
    EXITS("MESH_TOPOLOGY_ELEMENT_FINALISE")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_ELEMENT_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the given mesh topology element.
  SUBROUTINE MESH_TOPOLOGY_ELEMENT_INITIALISE(ELEMENT,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_ELEMENT_TYPE) :: ELEMENT !<The mesh element to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_TOPOLOGY_ELEMENT_INITIALISE",ERR,ERROR,*999)

    ELEMENT%USER_NUMBER=0
    ELEMENT%GLOBAL_NUMBER=0
    NULLIFY(ELEMENT%BASIS)
    ELEMENT%BOUNDARY_ELEMENT=.FALSE.
    
    EXITS("MESH_TOPOLOGY_ELEMENT_INITIALISE")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENT_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_ELEMENT_INITIALISE
  
  !
  !================================================================================================================================
  !

!!MERGE: Take user number
  
  !>Gets the basis for a mesh element identified by a given global number. \todo should take user number \see OPENCMISS::Iron::cmfe_MeshElementsBasisGet
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET(GLOBAL_NUMBER,ELEMENTS,BASIS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the element to get the basis for
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the elements to get the basis for \todo before number?
    TYPE(BASIS_TYPE), POINTER :: BASIS !<On return, a pointer to the basis to get
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(MESH_ELEMENT_TYPE), POINTER :: ELEMENT
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(.NOT.ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FlagError("Elements have been finished",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          ELEMENT=>ELEMENTS%ELEMENTS(GLOBAL_NUMBER)
          BASIS=>ELEMENT%BASIS
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET",ERR,ERROR)    
    RETURN 1
    
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET

  !
  !================================================================================================================================
  !

  !>Changes/sets the basis for a mesh element identified by a given global number. \todo should take user number
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET(GLOBAL_NUMBER,ELEMENTS,BASIS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the element to set the basis for
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the elements to set the basis for \todo before number?
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG), ALLOCATABLE :: NEW_USER_ELEMENT_NODES(:),NEW_GLOBAL_ELEMENT_NODES(:),NEW_USER_ELEMENT_NODE_VERSIONS(:,:)
    INTEGER(INTG) :: OVERLAPPING_NUMBER_NODES,OVERLAPPING_NUMBER_DERIVATIVES
    TYPE(MESH_ELEMENT_TYPE), POINTER :: ELEMENT
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FlagError("Elements have been finished",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          IF(ASSOCIATED(BASIS)) THEN
            ELEMENT=>ELEMENTS%ELEMENTS(GLOBAL_NUMBER)
            IF(ELEMENT%BASIS%NUMBER_OF_NODES/=BASIS%NUMBER_OF_NODES.OR. &
                & ELEMENT%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES/=BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES) THEN
              !Allocate new user and global element nodes
              ALLOCATE(NEW_USER_ELEMENT_NODES(BASIS%NUMBER_OF_NODES),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate new user element nodes",ERR,ERROR,*999)
              ALLOCATE(NEW_GLOBAL_ELEMENT_NODES(BASIS%NUMBER_OF_NODES),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate new user element nodes",ERR,ERROR,*999)
              ALLOCATE(NEW_USER_ELEMENT_NODE_VERSIONS(BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                & BASIS%NUMBER_OF_NODES),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate element node versions",ERR,ERROR,*999)

              OVERLAPPING_NUMBER_NODES=MIN(BASIS%NUMBER_OF_NODES,ELEMENT%BASIS%NUMBER_OF_NODES)
              OVERLAPPING_NUMBER_DERIVATIVES=MIN(BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES,ELEMENT%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES)

              !Set default values
              NEW_USER_ELEMENT_NODE_VERSIONS=1
              NEW_USER_ELEMENT_NODES(OVERLAPPING_NUMBER_NODES+1:)=0
              NEW_GLOBAL_ELEMENT_NODES(OVERLAPPING_NUMBER_NODES+1:)=0
              !Copy previous values
              NEW_USER_ELEMENT_NODES(1:OVERLAPPING_NUMBER_NODES)=ELEMENT%USER_ELEMENT_NODES(1:OVERLAPPING_NUMBER_NODES)
              NEW_GLOBAL_ELEMENT_NODES(1:OVERLAPPING_NUMBER_NODES)=ELEMENT%GLOBAL_ELEMENT_NODES(1:OVERLAPPING_NUMBER_NODES)
              NEW_USER_ELEMENT_NODE_VERSIONS(1:OVERLAPPING_NUMBER_DERIVATIVES,1:OVERLAPPING_NUMBER_NODES)= &
                & ELEMENT%USER_ELEMENT_NODE_VERSIONS(1:OVERLAPPING_NUMBER_DERIVATIVES,1:OVERLAPPING_NUMBER_NODES)

              !Replace arrays with new ones
              CALL MOVE_ALLOC(NEW_USER_ELEMENT_NODE_VERSIONS,ELEMENT%USER_ELEMENT_NODE_VERSIONS)
              CALL MOVE_ALLOC(NEW_USER_ELEMENT_NODES,ELEMENT%USER_ELEMENT_NODES)
              CALL MOVE_ALLOC(NEW_GLOBAL_ELEMENT_NODES,ELEMENT%GLOBAL_ELEMENT_NODES)
            ENDIF            
            ELEMENT%BASIS=>BASIS
          ELSE
            CALL FlagError("Basis is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET",ERR,ERROR)    
    RETURN 1
    
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET
  
  !
  !================================================================================================================================
  !

  !>Returns the adjacent element number for a mesh element identified by a global number. \todo specify by user number not global number \see OPENCMISS::Iron::cmfe_MeshElementsNo
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENT_GET(GLOBAL_NUMBER,ELEMENTS,ADJACENT_ELEMENT_XI,ADJACENT_ELEMENT_NUMBER, &
    & ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the element to get the adjacent element for
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the elements of a mesh component from which to get the adjacent element from.
    INTEGER(INTG), INTENT(IN) :: ADJACENT_ELEMENT_XI !< The xi coordinate direction to get the adjacent element Note that -xiCoordinateDirection gives the adjacent element before the element in the xiCoordinateDirection'th direction and +xiCoordinateDirection gives the adjacent element after the element in the xiCoordinateDirection'th direction. The xiCoordinateDirection=0 index will give the information on the current element.
    INTEGER(INTG), INTENT(OUT) :: ADJACENT_ELEMENT_NUMBER !<On return, the adjacent element number in the specified xi coordinate direction. Return 0 if the specified element has no adjacent elements in the specified xi coordinate direction.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENT_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(.NOT.ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FlagError("Elements have not been finished",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          IF(ADJACENT_ELEMENT_XI>=-ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_XI .AND. &
            & ADJACENT_ELEMENT_XI<=ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_XI) THEN
            IF(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%ADJACENT_ELEMENTS(ADJACENT_ELEMENT_XI)%NUMBER_OF_ADJACENT_ELEMENTS > 0) THEN !\todo Currently returns only the first adjacent element for now as the python binding require the output array size of the adjacent element to be known a-prior. Add routine to first output number of adjacent elements and then loop over all adjacent elements
              ADJACENT_ELEMENT_NUMBER=ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%ADJACENT_ELEMENTS(ADJACENT_ELEMENT_XI)%ADJACENT_ELEMENTS(1)
            ELSE !Return 0 indicating the specified element has no adjacent elements in the specified xi coordinate direction.
              ADJACENT_ELEMENT_NUMBER=0
            ENDIF
          ELSE
            LOCAL_ERROR="The specified adjacent element xi is invalid. The supplied xi is "// &
            & TRIM(NUMBER_TO_VSTRING(ADJACENT_ELEMENT_XI,"*",ERR,ERROR))//" and needs to be >=-"// &
            & TRIM(NUMBER_TO_VSTRING(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_XI,"*",ERR,ERROR))//" and <="// &
            & TRIM(NUMBER_TO_VSTRING(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_XI,"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated",ERR,ERROR,*999)
    ENDIF

    EXITS("MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENT_GET")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENT_GET",ERR,ERROR)
    RETURN 1

  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENT_GET

  !
  !================================================================================================================================
  !

  !>Gets the element nodes for a mesh element identified by a given global number. \todo specify by user number not global number \see OPENCMISS::Iron::cmfe_MeshElementsNodesGet
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET(GLOBAL_NUMBER,ELEMENTS,USER_ELEMENT_NODES,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the element to set the nodes for
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the elements to set \todo before number?
    INTEGER(INTG), INTENT(OUT) :: USER_ELEMENT_NODES(:) !<On return, USER_ELEMENT_NODES(i). USER_ELEMENT_NODES(i) is the i'th user node number for the element
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(.NOT.ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FlagError("Elements have not been finished",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          IF(SIZE(USER_ELEMENT_NODES,1)>=SIZE(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_ELEMENT_NODES,1)) THEN
            USER_ELEMENT_NODES=ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_ELEMENT_NODES
          ELSE
            LOCAL_ERROR="The size of USER_ELEMENT_NODES is too small. The supplied size is "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(USER_ELEMENT_NODES,1),"*",ERR,ERROR))//" and it needs to be >= "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_ELEMENT_NODES,1),"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET",ERR,ERROR)    
    RETURN 1
    
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET

  !
  !================================================================================================================================
  !

  !>Changes/sets the element nodes for a mesh element identified by a given global number. \todo specify by user number not global number \see OPENCMISS::Iron::cmfe_MeshElementsNodesSet
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(GLOBAL_NUMBER,ELEMENTS,USER_ELEMENT_NODES,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the element to set the nodes for
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the elements to set \todo before number?
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NODES(:) !<USER_ELEMENT_NODES(i). USER_ELEMENT_NODES(i) is the i'th user node number for the element
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nn,NUMBER_OF_BAD_NODES,GLOBAL_NODE_NUMBER
    INTEGER(INTG), ALLOCATABLE :: GLOBAL_ELEMENT_NODES(:),BAD_NODES(:)
    LOGICAL :: ELEMENT_NODES_OK,NODE_EXISTS
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE
    TYPE(MESH_TYPE), POINTER :: mesh
    TYPE(MeshComponentTopologyType), POINTER :: meshComponentTopology
    TYPE(NODES_TYPE), POINTER :: NODES
    TYPE(REGION_TYPE), POINTER :: PARENT_REGION,REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FlagError("Elements have been finished.",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          IF(SIZE(USER_ELEMENT_NODES,1)==ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_NODES) THEN
            meshComponentTopology=>elements%meshComponentTopology
            IF(ASSOCIATED(meshComponentTopology)) THEN
              mesh=>meshComponentTopology%mesh
              IF(ASSOCIATED(mesh)) THEN
                REGION=>MESH%REGION
                IF(ASSOCIATED(REGION)) THEN
                  NODES=>REGION%NODES
                ELSE
                  INTERFACE=>MESH%INTERFACE
                  IF(ASSOCIATED(INTERFACE)) THEN
                    NODES=>INTERFACE%NODES
                    PARENT_REGION=>INTERFACE%PARENT_REGION
                    IF(.NOT.ASSOCIATED(PARENT_REGION)) CALL FlagError("Mesh interface has no parent region.",ERR,ERROR,*999)
                  ELSE
                    CALL FlagError("Elements mesh has no associated region or interface.",ERR,ERROR,*999)
                  ENDIF
                ENDIF
                IF(ASSOCIATED(NODES)) THEN
                  ELEMENT_NODES_OK=.TRUE.
                  ALLOCATE(GLOBAL_ELEMENT_NODES(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_NODES),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate global element nodes.",ERR,ERROR,*999)
                  ALLOCATE(BAD_NODES(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_NODES),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate bad nodes.",ERR,ERROR,*999)
                  NUMBER_OF_BAD_NODES=0
                  DO nn=1,ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_NODES
                    CALL NODE_CHECK_EXISTS(NODES,USER_ELEMENT_NODES(nn),NODE_EXISTS,GLOBAL_NODE_NUMBER,ERR,ERROR,*999)
                    IF(NODE_EXISTS) THEN
                      GLOBAL_ELEMENT_NODES(nn)=GLOBAL_NODE_NUMBER
                    ELSE
                      NUMBER_OF_BAD_NODES=NUMBER_OF_BAD_NODES+1
                      BAD_NODES(NUMBER_OF_BAD_NODES)=USER_ELEMENT_NODES(nn)
                      ELEMENT_NODES_OK=.FALSE.
                    ENDIF
                  ENDDO !nn
                  IF(ELEMENT_NODES_OK) THEN
                    ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_ELEMENT_NODES=USER_ELEMENT_NODES
                    ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%GLOBAL_ELEMENT_NODES=GLOBAL_ELEMENT_NODES
                  ELSE
                    IF(NUMBER_OF_BAD_NODES==1) THEN
                      IF(ASSOCIATED(REGION)) THEN
                        LOCAL_ERROR="The element user node number of "//TRIM(NUMBER_TO_VSTRING(BAD_NODES(1),"*",ERR,ERROR))// &
                          & " is not defined in region "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                      ELSE
                        LOCAL_ERROR="The element user node number of "//TRIM(NUMBER_TO_VSTRING(BAD_NODES(1),"*",ERR,ERROR))// &
                          & " is not defined in interface number "// &
                          & TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))// &
                          & " of parent region number "//TRIM(NUMBER_TO_VSTRING(PARENT_REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                      ENDIF
                    ELSE
                      LOCAL_ERROR="The element user node number of "//TRIM(NUMBER_TO_VSTRING(BAD_NODES(1),"*",ERR,ERROR))
                      DO nn=2,NUMBER_OF_BAD_NODES-1
                        LOCAL_ERROR=LOCAL_ERROR//","//TRIM(NUMBER_TO_VSTRING(BAD_NODES(nn),"*",ERR,ERROR))
                      ENDDO !nn
                      IF(ASSOCIATED(REGION)) THEN
                        LOCAL_ERROR=LOCAL_ERROR//" & "//TRIM(NUMBER_TO_VSTRING(BAD_NODES(NUMBER_OF_BAD_NODES),"*",ERR,ERROR))// &
                          & " are not defined in region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                      ELSE
                        LOCAL_ERROR=LOCAL_ERROR//" & "//TRIM(NUMBER_TO_VSTRING(BAD_NODES(NUMBER_OF_BAD_NODES),"*",ERR,ERROR))// &
                          & " are not defined in interface number "// &
                          & TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))//" of parent region number "// &
                          &  TRIM(NUMBER_TO_VSTRING(PARENT_REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                      ENDIF
                    ENDIF
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  IF(ASSOCIATED(REGION)) THEN                   
                    CALL FlagError("The elements mesh region does not have any associated nodes.",ERR,ERROR,*999)
                  ELSE
                    CALL FlagError("The elements mesh interface does not have any associated nodes.",ERR,ERROR,*999) 
                  ENDIF
                ENDIF
              ELSE
                CALL FlagError("The mesh component topology mesh is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("The elements mesh component topology is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Number of element nodes does not match number of basis nodes for this element.",ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified global element number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The global element number should be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET",ERR,ERROR)    
    RETURN 1
    
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET

  !
  !================================================================================================================================
  !

  !>Changes/sets an element node's version for a mesh element identified by a given global number. \todo specify by user number not global number \see OPENCMISS::Iron::cmfe_MeshElementsNodesSet
  SUBROUTINE MeshElements_ElementNodeVersionSet(GLOBAL_NUMBER,ELEMENTS,VERSION_NUMBER,DERIVATIVE_NUMBER, &
      & USER_ELEMENT_NODE_INDEX,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the element to set the nodes for
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the elements to set \todo before number?
    INTEGER(INTG), INTENT(IN) :: VERSION_NUMBER !<The version number of the specified element node to set.
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<The derivative number of the specified element node to set.
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NODE_INDEX !< The node index of the specified element node to set a version for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(MeshComponentTopologyType), POINTER :: meshComponentTopology
    TYPE(NODES_TYPE), POINTER :: NODES
    TYPE(REGION_TYPE), POINTER :: PARENT_REGION,REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MeshElements_ElementNodeVersionSet",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FlagError("Elements have been finished.",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          IF(USER_ELEMENT_NODE_INDEX>=1.AND.USER_ELEMENT_NODE_INDEX<=ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_NODES) THEN
            meshComponentTopology=>elements%meshComponentTopology
            IF(ASSOCIATED(meshComponentTopology)) THEN
              mesh=>meshComponentTopology%mesh
              IF(ASSOCIATED(mesh)) THEN
                REGION=>mesh%REGION
                IF(ASSOCIATED(REGION)) THEN
                  NODES=>REGION%NODES
                ELSE
                  INTERFACE=>mesh%INTERFACE
                  IF(ASSOCIATED(INTERFACE)) THEN
                    NODES=>INTERFACE%NODES
                    PARENT_REGION=>INTERFACE%PARENT_REGION
                    IF(.NOT.ASSOCIATED(PARENT_REGION)) CALL FlagError("Mesh interface has no parent region.",ERR,ERROR,*999)
                  ELSE
                    CALL FlagError("Elements mesh has no associated region or interface.",ERR,ERROR,*999)
                  ENDIF
                ENDIF
                IF(DERIVATIVE_NUMBER>=1.AND.DERIVATIVE_NUMBER<=ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS% &
                  & NUMBER_OF_DERIVATIVES(USER_ELEMENT_NODE_INDEX)) THEN !Check if the specified derivative exists
                  IF(VERSION_NUMBER>=1) THEN !Check if the specified version is greater than 1
                    ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_ELEMENT_NODE_VERSIONS(DERIVATIVE_NUMBER,USER_ELEMENT_NODE_INDEX) & 
                      & = VERSION_NUMBER
                    !\todo : There is redunancy in USER_ELEMENT_NODE_VERSIONS since it was allocated in MESH_TOPOLOGY_ELEMENTS_CREATE_START based on MAXIMUM_NUMBER_OF_DERIVATIVES for that elements basis:ALLOCATE(ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODE_VERSIONS(BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES,BASIS%NUMBER_OF_NODES),STAT=ERR)
                  ELSE
                    LOCAL_ERROR="The specified node version number of "//TRIM(NUMBER_TO_VSTRING(VERSION_NUMBER,"*", & 
                      & ERR,ERROR))//" is invalid. The element node index should be greater than 1."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The specified node derivative number of "//TRIM(NUMBER_TO_VSTRING(DERIVATIVE_NUMBER,"*", & 
                    & ERR,ERROR))//" is invalid. The element node derivative index should be between 1 and "// &
                    & TRIM(NUMBER_TO_VSTRING(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_DERIVATIVES( &
                    & USER_ELEMENT_NODE_INDEX),"*",ERR,ERROR))//"."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("The mesh component topology mesh is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("The elements mesh component topology is not associated.",err,error,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified element node index of "//TRIM(NUMBER_TO_VSTRING(USER_ELEMENT_NODE_INDEX,"*",ERR,ERROR))// &
              & " is invalid. The element node index should be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_NODES,"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified global element number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The global element number should be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MeshElements_ElementNodeVersionSet")
    RETURN
999 ERRORSEXITS("MeshElements_ElementNodeVersionSet",ERR,ERROR)    
    RETURN 1
    
  END SUBROUTINE MeshElements_ElementNodeVersionSet

  !
  !================================================================================================================================
  !

 !>Calculates the element numbers surrounding an element in a mesh topology.
  SUBROUTINE MeshTopology_ElementsAdjacentElementsCalculate(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: TOPOLOGY !<A pointer to the mesh topology to calculate the elements adjacent to elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: j,ne,ne1,nep1,ni,nic,nn,nn1,nn2,nn3,node_idx,np,np1,DUMMY_ERR,FACE_XI(2),FACE_XIC(3),NODE_POSITION_INDEX(4)
    INTEGER(INTG) :: xi_direction,direction_index,xi_dir_check,xi_dir_search,NUMBER_NODE_MATCHES
    INTEGER(INTG) :: NUMBER_SURROUNDING,NUMBER_OF_NODES_XIC(4)
    INTEGER(INTG), ALLOCATABLE :: NODE_MATCHES(:),ADJACENT_ELEMENTS(:)
    LOGICAL :: XI_COLLAPSED,FACE_COLLAPSED(-3:3),SUBSET
    TYPE(LIST_TYPE), POINTER :: NODE_MATCH_LIST
    TYPE(LIST_PTR_TYPE) :: ADJACENT_ELEMENTS_LIST(-4:4)
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    NULLIFY(NODE_MATCH_LIST)
    DO nic=-4,4
      NULLIFY(ADJACENT_ELEMENTS_LIST(nic)%PTR)
    ENDDO !nic
    
    ENTERS("MeshTopology_ElementsAdjacentElementsCalculate",ERR,ERROR,*999)
    
    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%NODES)) THEN
        IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
          !Loop over the global elements in the mesh
          DO ne=1,TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS
            !%%%% first we initialize lists that are required to find the adjacent elements list
            BASIS=>TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS
            DO nic=-BASIS%NUMBER_OF_XI_COORDINATES,BASIS%NUMBER_OF_XI_COORDINATES
              NULLIFY(ADJACENT_ELEMENTS_LIST(nic)%PTR)
              CALL LIST_CREATE_START(ADJACENT_ELEMENTS_LIST(nic)%PTR,ERR,ERROR,*999)
              CALL LIST_DATA_TYPE_SET(ADJACENT_ELEMENTS_LIST(nic)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
              CALL LIST_INITIAL_SIZE_SET(ADJACENT_ELEMENTS_LIST(nic)%PTR,5,ERR,ERROR,*999)
              CALL LIST_CREATE_FINISH(ADJACENT_ELEMENTS_LIST(nic)%PTR,ERR,ERROR,*999)
            ENDDO !ni
            NUMBER_OF_NODES_XIC=1
            NUMBER_OF_NODES_XIC(1:BASIS%NUMBER_OF_XI_COORDINATES)=BASIS%NUMBER_OF_NODES_XIC(1:BASIS%NUMBER_OF_XI_COORDINATES)
            !Place the current element in the surrounding list
            CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(0)%PTR,TOPOLOGY%ELEMENTS%ELEMENTS(ne)%GLOBAL_NUMBER,ERR,ERROR,*999)
            SELECT CASE(BASIS%TYPE)
            CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
              !Determine the collapsed "faces" if any
              NODE_POSITION_INDEX=1
              !Loop over the face normals of the element
              DO ni=1,BASIS%NUMBER_OF_XI
                !Determine the xi directions that lie in this xi direction
                FACE_XI(1)=OTHER_XI_DIRECTIONS3(ni,2,1)
                FACE_XI(2)=OTHER_XI_DIRECTIONS3(ni,3,1)
                !Reset the node_position_index in this xi direction
                NODE_POSITION_INDEX(ni)=1
                !Loop over the two faces with this normal
                DO direction_index=-1,1,2
                  xi_direction=direction_index*ni
                  FACE_COLLAPSED(xi_direction)=.FALSE.
                  DO j=1,2
                    xi_dir_check=FACE_XI(j)
                    IF(xi_dir_check<=BASIS%NUMBER_OF_XI) THEN
                      xi_dir_search=FACE_XI(3-j)
                      NODE_POSITION_INDEX(xi_dir_search)=1
                      XI_COLLAPSED=.TRUE.
                      DO WHILE(NODE_POSITION_INDEX(xi_dir_search)<=NUMBER_OF_NODES_XIC(xi_dir_search).AND.XI_COLLAPSED)
                        !Get the first local node along the xi check direction
                        NODE_POSITION_INDEX(xi_dir_check)=1
                        nn1=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2),NODE_POSITION_INDEX(3),1)
                        !Get the second local node along the xi check direction
                        NODE_POSITION_INDEX(xi_dir_check)=2
                        nn2=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2),NODE_POSITION_INDEX(3),1)
                        IF(nn1/=0.AND.nn2/=0) THEN
                          IF(TOPOLOGY%ELEMENTS%ELEMENTS(ne)%MESH_ELEMENT_NODES(nn1)/= &
                            & TOPOLOGY%ELEMENTS%ELEMENTS(ne)%MESH_ELEMENT_NODES(nn2)) XI_COLLAPSED=.TRUE.
                        ENDIF
                        NODE_POSITION_INDEX(xi_dir_search)=NODE_POSITION_INDEX(xi_dir_search)+1
                      ENDDO !xi_dir_search
                      IF(XI_COLLAPSED) FACE_COLLAPSED(xi_direction)=.TRUE.
                    ENDIF
                  ENDDO !j
                  NODE_POSITION_INDEX(ni)=NUMBER_OF_NODES_XIC(ni)
                ENDDO !direction_index
              ENDDO !ni
              !Loop over the xi directions and calculate the surrounding elements
              DO ni=1,BASIS%NUMBER_OF_XI
                !Determine the xi directions that lie in this xi direction
                FACE_XI(1)=OTHER_XI_DIRECTIONS3(ni,2,1)
                FACE_XI(2)=OTHER_XI_DIRECTIONS3(ni,3,1)
                !Loop over the two faces
                DO direction_index=-1,1,2
                  xi_direction=direction_index*ni
                  !Find nodes in the element on the appropriate face/line/point
                  NULLIFY(NODE_MATCH_LIST)
                  CALL LIST_CREATE_START(NODE_MATCH_LIST,ERR,ERROR,*999)
                  CALL LIST_DATA_TYPE_SET(NODE_MATCH_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)

                  CALL LIST_INITIAL_SIZE_SET(NODE_MATCH_LIST,16,ERR,ERROR,*999)
                  CALL LIST_CREATE_FINISH(NODE_MATCH_LIST,ERR,ERROR,*999)
                  IF(direction_index==-1) THEN
                    NODE_POSITION_INDEX(ni)=1
                  ELSE
                    NODE_POSITION_INDEX(ni)=NUMBER_OF_NODES_XIC(ni)
                  ENDIF
                  !If the face is collapsed then don't look in this xi direction. The exception is if the opposite face is also
                  !collpased. This may indicate that we have a funny element in non-rc coordinates that goes around the central
                  !axis back to itself
                  IF(FACE_COLLAPSED(xi_direction).AND..NOT.FACE_COLLAPSED(-xi_direction)) THEN
                    !Do nothing - the match lists are already empty
                  ELSE
                    !Find the nodes to match and add them to the node match list
                    DO nn1=1,NUMBER_OF_NODES_XIC(FACE_XI(1))
                      NODE_POSITION_INDEX(FACE_XI(1))=nn1
                      DO nn2=1,NUMBER_OF_NODES_XIC(FACE_XI(2))
                        NODE_POSITION_INDEX(FACE_XI(2))=nn2
                        nn=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2),NODE_POSITION_INDEX(3),1)
                        IF(nn/=0) THEN
                          np=TOPOLOGY%ELEMENTS%ELEMENTS(ne)%MESH_ELEMENT_NODES(nn)
                          CALL LIST_ITEM_ADD(NODE_MATCH_LIST,np,ERR,ERROR,*999)
                        ENDIF
                      ENDDO !nn2
                    ENDDO !nn1
                  ENDIF
                  CALL LIST_REMOVE_DUPLICATES(NODE_MATCH_LIST,ERR,ERROR,*999)
                  CALL LIST_DETACH_AND_DESTROY(NODE_MATCH_LIST,NUMBER_NODE_MATCHES,NODE_MATCHES,ERR,ERROR,*999)
                  NUMBER_SURROUNDING=0
                  IF(NUMBER_NODE_MATCHES>0) THEN
                    !Find list of elements surrounding those nodes
                    np1=NODE_MATCHES(1)
                    DO nep1=1,TOPOLOGY%NODES%NODES(np1)%numberOfSurroundingElements
                      ne1=TOPOLOGY%NODES%NODES(np1)%surroundingElements(nep1)
                      IF(ne1/=ne) THEN !Don't want the current element
                        ! grab the nodes list for current and this surrouding elements
                        ! current face : NODE_MATCHES
                        ! candidate elem : TOPOLOGY%ELEMENTS%ELEMENTS(ne1)%MESH_ELEMENT_NODES ! should this be GLOBAL_ELEMENT_NODES?
                        ! if all of current face belongs to the candidate element, we will have found the neighbour
                        CALL LIST_SUBSET_OF(NODE_MATCHES(1:NUMBER_NODE_MATCHES),TOPOLOGY%ELEMENTS%ELEMENTS(ne1)% &
                          & MESH_ELEMENT_NODES,SUBSET,ERR,ERROR,*999)
                        IF(SUBSET) THEN
                          CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(xi_direction)%PTR,ne1,ERR,ERROR,*999)
                          NUMBER_SURROUNDING=NUMBER_SURROUNDING+1
                        ENDIF
                      ENDIF
                    ENDDO !nep1
                  ENDIF
                  IF(ALLOCATED(NODE_MATCHES)) DEALLOCATE(NODE_MATCHES)
                ENDDO !direction_index
              ENDDO !ni
            CASE(BASIS_SIMPLEX_TYPE)
              !Loop over the xi coordinates and calculate the surrounding elements
              DO nic=1,BASIS%NUMBER_OF_XI_COORDINATES
                !Find the other coordinates of the face/line/point
                FACE_XIC(1)=OTHER_XI_DIRECTIONS4(nic,1)
                FACE_XIC(2)=OTHER_XI_DIRECTIONS4(nic,2)
                FACE_XIC(3)=OTHER_XI_DIRECTIONS4(nic,3)
                !Find nodes in the element on the appropriate face/line/point
                NULLIFY(NODE_MATCH_LIST)
                CALL LIST_CREATE_START(NODE_MATCH_LIST,ERR,ERROR,*999)
                CALL LIST_DATA_TYPE_SET(NODE_MATCH_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                CALL LIST_INITIAL_SIZE_SET(NODE_MATCH_LIST,16,ERR,ERROR,*999)
                CALL LIST_CREATE_FINISH(NODE_MATCH_LIST,ERR,ERROR,*999)
                NODE_POSITION_INDEX(nic)=1 !Furtherest away from node with the nic'th coordinate
                !Find the nodes to match and add them to the node match list
                DO nn1=1,NUMBER_OF_NODES_XIC(FACE_XIC(1))
                  NODE_POSITION_INDEX(FACE_XIC(1))=nn1
                  DO nn2=1,NUMBER_OF_NODES_XIC(FACE_XIC(2))
                    NODE_POSITION_INDEX(FACE_XIC(2))=nn2
                    DO nn3=1,NUMBER_OF_NODES_XIC(FACE_XIC(3))
                      NODE_POSITION_INDEX(FACE_XIC(3))=nn3
                      nn=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2),NODE_POSITION_INDEX(3), &
                        NODE_POSITION_INDEX(4))
                      IF(nn/=0) THEN
                        np=TOPOLOGY%ELEMENTS%ELEMENTS(ne)%MESH_ELEMENT_NODES(nn)
                        CALL LIST_ITEM_ADD(NODE_MATCH_LIST,np,ERR,ERROR,*999)
                      ENDIF
                    ENDDO !nn3
                  ENDDO !nn2
                ENDDO !nn1
                CALL LIST_REMOVE_DUPLICATES(NODE_MATCH_LIST,ERR,ERROR,*999)
                CALL LIST_DETACH_AND_DESTROY(NODE_MATCH_LIST,NUMBER_NODE_MATCHES,NODE_MATCHES,ERR,ERROR,*999)
                IF(NUMBER_NODE_MATCHES>0) THEN
                  !Find list of elements surrounding those nodes
                  DO node_idx=1,NUMBER_NODE_MATCHES
                    np1=NODE_MATCHES(node_idx)
                    DO nep1=1,TOPOLOGY%NODES%NODES(np1)%numberOfSurroundingElements
                      ne1=TOPOLOGY%NODES%NODES(np1)%surroundingElements(nep1)
                      IF(ne1/=ne) THEN !Don't want the current element
                        ! grab the nodes list for current and this surrouding elements
                        ! current face : NODE_MATCHES
                        ! candidate elem : TOPOLOGY%ELEMENTS%ELEMENTS(ne1)%MESH_ELEMENT_NODES 
                        ! if all of current face belongs to the candidate element, we will have found the neighbour
                        CALL LIST_SUBSET_OF(NODE_MATCHES(1:NUMBER_NODE_MATCHES),TOPOLOGY%ELEMENTS%ELEMENTS(ne1)% &
                          & MESH_ELEMENT_NODES,SUBSET,ERR,ERROR,*999)
                        IF(SUBSET) THEN
                          CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(nic)%PTR,ne1,ERR,ERROR,*999)
                        ENDIF
                      ENDIF
                    ENDDO !nep1
                  ENDDO !node_idx
                ENDIF
                IF(ALLOCATED(NODE_MATCHES)) DEALLOCATE(NODE_MATCHES)
              ENDDO !nic
            CASE(BASIS_SERENDIPITY_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE(BASIS_AUXILLIARY_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE(BASIS_B_SPLINE_TP_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE(BASIS_FOURIER_LAGRANGE_HERMITE_TP_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE(BASIS_EXTENDED_LAGRANGE_TP_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The basis type of "//TRIM(NUMBER_TO_VSTRING(BASIS%TYPE,"*",ERR,ERROR))// &
                & " is invalid."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
            !Set the surrounding elements for this element
            ALLOCATE(TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(-BASIS%NUMBER_OF_XI_COORDINATES: &
              & BASIS%NUMBER_OF_XI_COORDINATES),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate adjacent elements.",ERR,ERROR,*999)
            DO nic=-BASIS%NUMBER_OF_XI_COORDINATES,BASIS%NUMBER_OF_XI_COORDINATES
              CALL MESH_ADJACENT_ELEMENT_INITIALISE(TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic),ERR,ERROR,*999)
              CALL LIST_REMOVE_DUPLICATES(ADJACENT_ELEMENTS_LIST(nic)%PTR,ERR,ERROR,*999)
              CALL LIST_DETACH_AND_DESTROY(ADJACENT_ELEMENTS_LIST(nic)%PTR,TOPOLOGY%ELEMENTS%ELEMENTS(ne)% &
                & ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS,ADJACENT_ELEMENTS,ERR,ERROR,*999)
              ALLOCATE(TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%ADJACENT_ELEMENTS(TOPOLOGY%ELEMENTS%ELEMENTS(ne)% &
                ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate element adjacent elements.",ERR,ERROR,*999)
              TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%ADJACENT_ELEMENTS(1:TOPOLOGY%ELEMENTS%ELEMENTS(ne)% &
                ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS) = ADJACENT_ELEMENTS(1:TOPOLOGY%ELEMENTS% &
                & ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS)
              IF(ALLOCATED(ADJACENT_ELEMENTS)) DEALLOCATE(ADJACENT_ELEMENTS)
            ENDDO !nic
          ENDDO !ne           
        ELSE
          CALL FlagError("Mesh topology elements is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Mesh topology nodes is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh topology is not allocated.",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Number of elements = ",TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS,ERR,ERROR,*999)
      DO ne=1,TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS
        BASIS=>TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Global element number : ",ne,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of xi coordinates = ",BASIS%NUMBER_OF_XI_COORDINATES, &
          & ERR,ERROR,*999)
        DO nic=-BASIS%NUMBER_OF_XI_COORDINATES,BASIS%NUMBER_OF_XI_COORDINATES
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Xi coordinate : ",nic,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of adjacent elements = ", &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS,ERR,ERROR,*999)
          IF(TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS>0) THEN
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,TOPOLOGY%ELEMENTS%ELEMENTS(ne)% &
              & ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS,8,8,TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)% &
              & ADJACENT_ELEMENTS,'("        Adjacent elements :",8(X,I8))','(30x,8(X,I8))',ERR,ERROR,*999)
          ENDIF
        ENDDO !nic
      ENDDO !ne
    ENDIF
    
    EXITS("MeshTopology_ElementsAdjacentElementsCalculate")
    RETURN
999 IF(ALLOCATED(NODE_MATCHES)) DEALLOCATE(NODE_MATCHES)
    IF(ALLOCATED(ADJACENT_ELEMENTS)) DEALLOCATE(ADJACENT_ELEMENTS)
    IF(ASSOCIATED(NODE_MATCH_LIST)) CALL LIST_DESTROY(NODE_MATCH_LIST,DUMMY_ERR,DUMMY_ERROR,*998)
998 DO nic=-4,4
      IF(ASSOCIATED(ADJACENT_ELEMENTS_LIST(nic)%PTR)) CALL LIST_DESTROY(ADJACENT_ELEMENTS_LIST(nic)%PTR,DUMMY_ERR,DUMMY_ERROR,*997)
    ENDDO !ni
997 ERRORS("MeshTopology_ElementsAdjacentElementsCalculate",ERR,ERROR)
    EXITS("MeshTopology_ElementsAdjacentElementsCalculate")
    RETURN 1
    
  END SUBROUTINE MeshTopology_ElementsAdjacentElementsCalculate
  
  !
  !================================================================================================================================
  !

  !>Finalises the elements data structures for a mesh topology and deallocates any memory. \todo pass in elements
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_FINALISE(ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the mesh topology to finalise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ne

    ENTERS("MESH_TOPOLOGY_ELEMENTS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      DO ne=1,ELEMENTS%NUMBER_OF_ELEMENTS
        CALL MESH_TOPOLOGY_ELEMENT_FINALISE(ELEMENTS%ELEMENTS(ne),ERR,ERROR,*999)
      ENDDO !ne
      DEALLOCATE(ELEMENTS%ELEMENTS)
      IF(ASSOCIATED(ELEMENTS%ELEMENTS_TREE)) CALL TREE_DESTROY(ELEMENTS%ELEMENTS_TREE,ERR,ERROR,*999)
      DEALLOCATE(ELEMENTS)
    ENDIF
 
    EXITS("MESH_TOPOLOGY_ELEMENTS_FINALISE")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the elements in a given mesh topology. \todo finalise on error
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: TOPOLOGY !<A pointer to the mesh topology to initialise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_TOPOLOGY_ELEMENTS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        CALL FlagError("Mesh already has topology elements associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%ELEMENTS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology elements",ERR,ERROR,*999)
        TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS=0
        TOPOLOGY%ELEMENTS%meshComponentTopology=>TOPOLOGY
        NULLIFY(TOPOLOGY%ELEMENTS%ELEMENTS)
        NULLIFY(TOPOLOGY%ELEMENTS%ELEMENTS_TREE)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_TOPOLOGY_ELEMENTS_INITIALISE")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the elements in a given mesh topology. \todo finalise on error
  SUBROUTINE MESH_TOPOLOGY_DATA_POINTS_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: TOPOLOGY !<A pointer to the mesh topology to initialise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_TOPOLOGY_DATA_POINTS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%dataPoints)) THEN
        CALL FlagError("Mesh already has topology data points associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%dataPoints,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology data points",ERR,ERROR,*999)
        TOPOLOGY%dataPoints%totalNumberOfProjectedData=0
        TOPOLOGY%dataPoints%meshComponentTopology=>topology
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_TOPOLOGY_DATA_POINTS_INITIALISE")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_DATA_POINTS_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_DATA_POINTS_INITIALISE
  
  !
  !================================================================================================================================
  !

!!MERGE: ditto.
  
  !>Gets the user number for a global element identified by a given global number. \todo Check that the user number doesn't already exist. \see OPENCMISS::Iron::cmfe_MeshElementsUserNumberGet
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_NUMBER_GET(GLOBAL_NUMBER,USER_NUMBER,ELEMENTS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the elements to get.
    INTEGER(INTG), INTENT(OUT) :: USER_NUMBER !<The user number of the element to get
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<On return, a pointer to the elements to get the user number for \todo This should be the first parameter.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_TOPOLOGY_ELEMENTS_NUMBER_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FlagError("Elements have been finished",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          USER_NUMBER=ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_NUMBER
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_TOPOLOGY_ELEMENTS_NUMBER_GET")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_NUMBER_GET",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_NUMBER_GET

  !
  !================================================================================================================================
  !

  !>Returns the user number for a global element identified by a given global number. \see OPENCMISS::Iron::cmfe_MeshElementsUserNumberGet
  SUBROUTINE MeshElements_ElementUserNumberGet(GLOBAL_NUMBER,USER_NUMBER,ELEMENTS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the elements to get.
    INTEGER(INTG), INTENT(OUT) :: USER_NUMBER !<The user number of the element to get
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the elements to set the user number for \todo This should be the first parameter.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MeshElements_ElementUserNumberGet",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%ELEMENTS_FINISHED) THEN
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN          
          USER_NUMBER=ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_NUMBER
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Elements have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MeshElements_ElementUserNumberGet")
    RETURN
999 ERRORSEXITS("MeshElements_ElementUserNumberGet",ERR,ERROR)
    RETURN 1

   
  END SUBROUTINE MeshElements_ElementUserNumberGet

  !
  !================================================================================================================================
  !

  !>Changes/sets the user number for a global element identified by a given global number. \see OPENCMISS::Iron::cmfe_MeshElementsUserNumberSet
  SUBROUTINE MeshElements_ElementUserNumberSet(GLOBAL_NUMBER,USER_NUMBER,ELEMENTS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the elements to set.
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the element to set
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the elements to set the user number for \todo This should be the first parameter.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code

    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: GLOBAL_ELEMENT_NUMBER,INSERT_STATUS
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MeshElements_ElementUserNumberSet",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FlagError("Elements have been finished.",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          NULLIFY(TREE_NODE)
          CALL TREE_SEARCH(ELEMENTS%ELEMENTS_TREE,USER_NUMBER,TREE_NODE,ERR,ERROR,*999)
          IF(ASSOCIATED(TREE_NODE)) THEN
            CALL TREE_NODE_VALUE_GET(ELEMENTS%ELEMENTS_TREE,TREE_NODE,GLOBAL_ELEMENT_NUMBER,ERR,ERROR,*999)
            LOCAL_ERROR="Element user number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
              & " is already used by global element number "// &
              & TRIM(NUMBER_TO_VSTRING(GLOBAL_ELEMENT_NUMBER,"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            CALL TREE_ITEM_DELETE(ELEMENTS%ELEMENTS_TREE,ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_NUMBER,ERR,ERROR,*999)
            CALL TREE_ITEM_INSERT(ELEMENTS%ELEMENTS_TREE,USER_NUMBER,GLOBAL_NUMBER,INSERT_STATUS,ERR,ERROR,*999)
            ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_NUMBER=USER_NUMBER
          ENDIF
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MeshElements_ElementUserNumberSet")
    RETURN
999 ERRORSEXITS("MeshElements_ElementUserNumberSet",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MeshElements_ElementUserNumberSet
  
  !
  !================================================================================================================================
  !

  !>Changes/sets the user numbers for all elements.
  SUBROUTINE MeshTopology_ElementsUserNumbersAllSet(elements,userNumbers,err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: elements !<A pointer to the elements to set all the user numbers for 
    INTEGER(INTG), INTENT(IN) :: userNumbers(:) !<The user numbers to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx,insertStatus
    TYPE(TREE_TYPE), POINTER :: newElementsTree
    TYPE(VARYING_STRING) :: localError

    NULLIFY(newElementsTree)

    ENTERS("MeshTopology_ElementsUserNumbersAllSet",err,error,*999)

    IF(ASSOCIATED(elements)) THEN
      IF(elements%ELEMENTS_FINISHED) THEN
        CALL FlagError("Elements have been finished.",err,error,*999)
      ELSE
        IF(elements%NUMBER_OF_ELEMENTS==SIZE(userNumbers,1)) THEN
          !Check the users numbers to ensure that there are no duplicates          
          CALL TREE_CREATE_START(newElementsTree,err,error,*999)
          CALL TREE_INSERT_TYPE_SET(newElementsTree,TREE_NO_DUPLICATES_ALLOWED,err,error,*999)
          CALL TREE_CREATE_FINISH(newElementsTree,err,error,*999)
          DO elementIdx=1,elements%NUMBER_OF_ELEMENTS
            CALL TREE_ITEM_INSERT(newElementsTree,userNumbers(elementIdx),elementIdx,insertStatus,err,error,*999)
            IF(insertStatus/=TREE_NODE_INSERT_SUCESSFUL) THEN
              localError="The specified user number of "//TRIM(NumberToVstring(userNumbers(elementIdx),"*",err,error))// &
                & " for global element number "//TRIM(NUMBER_TO_VSTRING(elementIdx,"*",err,error))// &
                & " is a duplicate. The user element numbers must be unique."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ENDDO !elementIdx
          CALL TREE_DESTROY(elements%ELEMENTS_TREE,err,error,*999)
          elements%ELEMENTS_TREE=>newElementsTree
          NULLIFY(newElementsTree)
          DO elementIdx=1,elements%NUMBER_OF_ELEMENTS
            elements%ELEMENTS(elementIdx)%GLOBAL_NUMBER=elementIdx
            elements%ELEMENTS(elementIdx)%USER_NUMBER=userNumbers(elementIdx)
          ENDDO !elementIdx
        ELSE
          localError="The number of specified element user numbers ("// &
            TRIM(NumberToVstring(SIZE(userNumbers,1),"*",err,error))// &
            ") does not match number of elements ("// &
            TRIM(NumberToVstring(elements%NUMBER_OF_ELEMENTS,"*",err,error))//")."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopology_ElementsUserNumbersAllSet")
    RETURN
999 IF(ASSOCIATED(newElementsTree)) CALL TREE_DESTROY(newElementsTree,err,error,*998)
998 ERRORSEXITS("MeshTopology_ElementsUserNumbersAllSet",err,error)    
    RETURN 1
   
  END SUBROUTINE MeshTopology_ElementsUserNumbersAllSet
  
  !
  !================================================================================================================================
  !

  !>Calculates the data points in the given mesh topology.
  SUBROUTINE MeshTopology_DataPointsCalculateProjection(mesh,dataProjection,err,error,*)
  
    !Argument variables
    TYPE(MESH_TYPE), POINTER :: mesh !<A pointer to the mesh topology to calcualte the data projection for
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points
    TYPE(MeshDataPointsType), POINTER :: dataPointsTopology
    TYPE(DataProjectionResultType), POINTER :: dataProjectionResult
    TYPE(MeshElementsType), POINTER :: elements
    INTEGER(INTG) :: dataPointIdx,elementIdx,exitTag,countIdx,projectionNumber,globalCountIdx,elementNumber

    ENTERS("MeshTopology_DataPointsCalculateProjection",ERR,ERROR,*999)

    IF(ASSOCIATED(mesh)) THEN
      IF(dataProjection%dataProjectionFinished) THEN 
        dataPoints=>dataProjection%dataPoints
        !Default the first mesh component topology to contain data points ! \TODO: need to be changed once the data points topology is moved under meshTopologyType.
        dataPointsTopology=>mesh%topology(1)%PTR%dataPoints
        !Extract the global number of the data projection 
        projectionNumber=dataProjection%globalNumber
        !Hard code the first mesh component since element topology is the same for all mesh components
        !\TODO: need to be changed once the elements topology is moved under meshTopologyType.
        elements=>mesh%TOPOLOGY(1)%PTR%ELEMENTS
        ALLOCATE(dataPointsTopology%elementDataPoint(elements%NUMBER_OF_ELEMENTS),STAT=ERR)     
        IF(ERR/=0) CALL FlagError("Could not allocate data points topology element.",ERR,ERROR,*999)
        DO elementIdx=1,elements%NUMBER_OF_ELEMENTS
          dataPointsTopology%elementDataPoint(elementIdx)%elementNumber=elements%ELEMENTS(elementIdx)%GLOBAL_NUMBER
          dataPointsTopology%elementDataPoint(elementIdx)%numberOfProjectedData=0
        ENDDO        
        !Calculate number of projected data points on an element
        DO dataPointIdx=1,dataPoints%numberOfDataPoints
          dataProjectionResult=>dataProjection%dataProjectionResults(dataPointIdx)
          elementNumber=dataProjectionResult%elementNumber
          exitTag=dataProjectionResult%exitTag
          IF(exitTag/=DATA_PROJECTION_CANCELLED) THEN
            DO elementIdx=1,elements%NUMBER_OF_ELEMENTS
              IF(dataPointsTopology%elementDataPoint(elementIdx)%elementNumber==elementNumber) THEN
                dataPointsTopology%elementDataPoint(elementIdx)%numberOfProjectedData= &
                  & dataPointsTopology%elementDataPoint(elementIdx)%numberOfProjectedData+1;
              ENDIF
            ENDDO !elementIdx
          ENDIF
        ENDDO !dataPointIdx      
        !Allocate memory to store data indices and initialise them to be zero   
        DO elementIdx=1,elements%NUMBER_OF_ELEMENTS
          ALLOCATE(dataPointsTopology%elementDataPoint(elementIdx)%dataIndices(dataPointsTopology% &
            & elementDataPoint(elementIdx)%numberOfProjectedData),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate data points topology element data points.",ERR,ERROR,*999)
          DO countIdx=1,dataPointsTopology%elementDataPoint(elementIdx)%numberOfProjectedData
            dataPointsTopology%elementDataPoint(elementIdx)%dataIndices(countIdx)%userNumber=0
            dataPointsTopology%elementDataPoint(elementIdx)%dataIndices(countIdx)%globalNumber=0
          ENDDO
        ENDDO     
        !Record the indices of the data that projected on the elements 
        globalCountIdx=0
        dataPointsTopology%totalNumberOfProjectedData=0
        DO dataPointIdx=1,dataPoints%numberOfDataPoints 
          dataProjectionResult=>dataProjection%dataProjectionResults(dataPointIdx)
          elementNumber=dataProjectionResult%elementNumber
          exitTag=dataProjectionResult%exitTag
          IF(exitTag/=DATA_PROJECTION_CANCELLED) THEN
            DO elementIdx=1,elements%NUMBER_OF_ELEMENTS
              countIdx=1         
              IF(dataPointsTopology%elementDataPoint(elementIdx)%elementNumber==elementNumber) THEN
                globalCountIdx=globalCountIdx+1
                !Find the next data point index in this element
                DO WHILE(dataPointsTopology%elementDataPoint(elementIdx)%dataIndices(countIdx)%globalNumber/=0)
                  countIdx=countIdx+1
                ENDDO
                dataPointsTopology%elementDataPoint(elementIdx)%dataIndices(countIdx)%userNumber=dataPointIdx
                dataPointsTopology%elementDataPoint(elementIdx)%dataIndices(countIdx)%globalNumber=dataPointIdx!globalCountIdx (used this if only projected data are taken into account)
                dataPointsTopology%totalNumberOfProjectedData=dataPointsTopology%totalNumberOfProjectedData+1
              ENDIF
            ENDDO !elementIdx
          ENDIF
        ENDDO !dataPointIdx
        !Allocate memory to store total data indices in ascending order and element map
        ALLOCATE(dataPointsTopology%dataPoints(dataPointsTopology%totalNumberOfProjectedData),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate data points topology data points.",ERR,ERROR,*999)
        !The global number for the data points will be looping through elements.
        countIdx=1  
        DO elementIdx=1,elements%NUMBER_OF_ELEMENTS
          DO dataPointIdx=1,dataPointsTopology%elementDataPoint(elementIdx)%numberOfProjectedData
            dataPointsTopology%dataPoints(countIdx)%userNumber=dataPointsTopology%elementDataPoint(elementIdx)% &
              & dataIndices(dataPointIdx)%userNumber
             dataPointsTopology%dataPoints(countIdx)%globalNumber=dataPointsTopology%elementDataPoint(elementIdx)% &
              & dataIndices(dataPointIdx)%globalNumber
             dataPointsTopology%dataPoints(countIdx)%elementNumber=dataPointsTopology%elementDataPoint(elementIdx)% &
               & elementNumber
             countIdx=countIdx+1
          ENDDO !dataPointIdx
        ENDDO !elementIdx                      
      ELSE
        CALL FlagError("Data projection is not finished.",err,error,*999)
      ENDIF     
    ELSE
      CALL FlagError("Mesh is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopology_DataPointsCalculateProjection")
    RETURN
999 ERRORSEXITS("MeshTopology_DataPointsCalculateProjection",err,error)
    RETURN 1
  END SUBROUTINE MeshTopology_DataPointsCalculateProjection

  !
  !================================================================================================================================
  !

  !>Finalises the topology in the given mesh. \todo pass in the mesh topology
  SUBROUTINE MeshTopology_Finalise(mesh,err,error,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: mesh !<A pointer to the mesh to finalise the topology for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx

    ENTERS("MeshTopology_Finalise",err,error,*999)

    IF(ASSOCIATED(mesh)) THEN
      DO componentIdx=1,mesh%NUMBER_OF_COMPONENTS
        CALL MeshTopology_ComponentFinalise(mesh%topology(componentIdx)%ptr,err,error,*999)
      ENDDO !componentIdx
      DEALLOCATE(mesh%topology)
    ELSE
      CALL FlagError("Mesh is not associated.",err,error,*999)
    ENDIF
 
    EXITS("MeshTopology_Finalise")
    RETURN
999 ERRORSEXITS("MeshTopology_Finalise",err,error)
    RETURN 1
   
  END SUBROUTINE MeshTopology_Finalise

  !
  !================================================================================================================================
  !

  !>Finalises the topology in the given mesh. 
  SUBROUTINE MeshTopology_ComponentFinalise(meshComponent,err,error,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: meshComponent !<A pointer to the mesh component to finalise the topology for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopology_ComponentFinalise",err,error,*999)

    IF(ASSOCIATED(meshComponent)) THEN
      CALL MeshTopology_NodesFinalise(meshComponent%nodes,err,error,*999)
      CALL MESH_TOPOLOGY_ELEMENTS_FINALISE(meshComponent%elements,err,error,*999)
      CALL MeshTopology_DofsFinalise(meshComponent%dofs,err,error,*999)
      DEALLOCATE(meshComponent)
    ENDIF
 
    EXITS("MeshTopology_ComponentFinalise")
    RETURN
999 ERRORSEXITS("MeshTopology_ComponentFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE MeshTopology_ComponentFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the topology for a given mesh. \todo finalise on error
  SUBROUTINE MeshTopology_Initialise(mesh,err,error,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: mesh !<A pointer to the mesh to initialise the mesh topology for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx
    
    ENTERS("MeshTopology_Initialise",err,error,*999)

    IF(ASSOCIATED(mesh)) THEN
      IF(ASSOCIATED(mesh%topology)) THEN
        CALL FlagError("Mesh already has topology associated.",err,error,*999)
      ELSE
        !Allocate mesh topology
        ALLOCATE(mesh%topology(mesh%NUMBER_OF_COMPONENTS),STAT=err)
        IF(err/=0) CALL FlagError("Mesh topology could not be allocated.",err,error,*999)
        DO componentIdx=1,mesh%NUMBER_OF_COMPONENTS
          ALLOCATE(mesh%topology(componentIdx)%ptr,STAT=err)
          IF(err/=0) CALL FlagError("Mesh topology component could not be allocated.",err,error,*999)
          mesh%topology(componentIdx)%ptr%mesh=>mesh
          NULLIFY(mesh%topology(componentIdx)%ptr%elements)
          NULLIFY(mesh%topology(componentIdx)%ptr%nodes)
          NULLIFY(mesh%topology(componentIdx)%ptr%dofs)
          NULLIFY(mesh%topology(componentIdx)%ptr%dataPoints)
          !Initialise the topology components
          CALL MESH_TOPOLOGY_ELEMENTS_INITIALISE(mesh%topology(componentIdx)%ptr,err,error,*999)
          CALL MeshTopology_NodesInitialise(mesh%topology(componentIdx)%ptr,err,error,*999)
          CALL MeshTopology_DofsInitialise(mesh%topology(componentIdx)%ptr,err,error,*999)
          CALL MESH_TOPOLOGY_DATA_POINTS_INITIALISE(mesh%topology(componentIdx)%ptr,err,error,*999)
        ENDDO !componentIdx
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopology_Initialise")
    RETURN
999 ERRORSEXITS("MeshTopology_Initialise",err,error)
    RETURN 1
  END SUBROUTINE MeshTopology_Initialise
  
  !
  !================================================================================================================================
  !

  !>Checks that a user element number exists in a mesh component. 
  SUBROUTINE MeshTopology_ElementCheckExistsMesh(mesh,meshComponentNumber,userElementNumber,elementExists,globalElementNumber, &
    & err,error,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: mesh !<A pointer to the mesh to check the element exists on
    INTEGER(INTG), INTENT(IN) :: meshComponentNumber !<The mesh component to check the element exits on
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number to check if it exists
    LOGICAL, INTENT(OUT) :: elementExists !<On exit, is .TRUE. if the element user number exists in the mesh component, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: globalElementNumber !<On exit, if the element exists the global number corresponding to the user element number. If the element does not exist then global number will be 0.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(MeshElementsType), POINTER :: elements

    NULLIFY(elements)
    
    ENTERS("MeshTopology_ElementCheckExistsMesh",err,error,*999)

    IF(ASSOCIATED(mesh)) THEN
      IF(mesh%MESH_FINISHED) THEN
        CALL MESH_TOPOLOGY_ELEMENTS_GET(mesh,meshComponentNumber,elements,err,error,*999)
        CALL MeshTopology_ElementCheckExistsMeshElements(elements,userElementNumber,elementExists,globalElementNumber, &
          & err,error,*999)
      ELSE
        CALL FlagError("Mesh has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopology_ElementCheckExistsMesh")
    RETURN
999 ERRORSEXITS("MeshTopology_ElementCheckExistsMesh",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_ElementCheckExistsMesh
  
  !
  !================================================================================================================================
  !

  !>Checks that a user element number exists in a mesh elements. 
  SUBROUTINE MeshTopology_ElementCheckExistsMeshElements(meshElements,userElementNumber,elementExists,globalElementNumber, &
    & err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: meshElements !<A pointer to the mesh elements to check the element exists on
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number to check if it exists
    LOGICAL, INTENT(OUT) :: elementExists !<On exit, is .TRUE. if the element user number exists in the mesh elements, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: globalElementNumber !<On exit, if the element exists the global number corresponding to the user element number. If the element does not exist then global number will be 0.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(TREE_NODE_TYPE), POINTER :: treeNode
    
    ENTERS("MeshTopology_ElementCheckExistsMeshElements",err,error,*999)

    elementExists=.FALSE.
    globalElementNumber=0
    IF(ASSOCIATED(meshElements)) THEN
      NULLIFY(treeNode)
      CALL TREE_SEARCH(meshElements%ELEMENTS_TREE,userElementNumber,treeNode,err,error,*999)
      IF(ASSOCIATED(treeNode)) THEN
        CALL TREE_NODE_VALUE_GET(meshElements%ELEMENTS_TREE,treeNode,globalElementNumber,err,error,*999)
        elementExists=.TRUE.
      ENDIF
    ELSE
      CALL FlagError("Mesh elements is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopology_ElementCheckExistsMeshElements")
    RETURN
999 ERRORSEXITS("MeshTopology_ElementCheckExistsMeshElements",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_ElementCheckExistsMeshElements
  
  !
  !================================================================================================================================
  !

  !>Gets a mesh element number that corresponds to a user element number from mesh element. An error will be raised if the user element number does not exist.
  SUBROUTINE MeshTopology_ElementGetMeshElements(meshElements,userElementNumber,globalElementNumber,err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: meshElements !<A pointer to the mesh elements to get the element for
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number to get
    INTEGER(INTG), INTENT(OUT) :: globalElementNumber !<On exit, the global number corresponding to the user element number. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    LOGICAL :: elementExists
    TYPE(MESH_TYPE), POINTER :: mesh
    TYPE(MeshComponentTopologyType), POINTER :: meshComponentTopology
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("MeshTopology_ElementGetMeshElements",err,error,*999)

    CALL MeshTopology_ElementCheckExistsMeshElements(meshElements,userElementNumber,elementExists,globalElementNumber, &
      & err,error,*999)
    IF(.NOT.elementExists) THEN
      meshComponentTopology=>meshElements%meshComponentTopology
      IF(ASSOCIATED(meshComponentTopology)) THEN
        mesh=>meshComponentTopology%mesh
        IF(ASSOCIATED(mesh)) THEN
          localError="The user element number "//TRIM(NumberToVString(userElementNumber,"*",err,error))// &
            & " does not exist in mesh number "//TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ELSE
          localError="The user element number "//TRIM(NumberToVString(userElementNumber,"*",err,error))//" does not exist."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        localError="The user element number "//TRIM(NumberToVString(userElementNumber,"*",err,error))//" does not exist."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDIF
   
    EXITS("MeshTopology_ElementGetMeshElements")
    RETURN
999 globalElementNumber=0
    ERRORSEXITS("MeshTopology_ElementGetMeshElements",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_ElementGetMeshElements
  
  !
  !================================================================================================================================
  !
  
  !>Returns if the element in a mesh is on the boundary or not
  SUBROUTINE MeshTopology_ElementOnBoundaryGet(meshElements,userNumber,onBoundary,err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: meshElements !<A pointer to the mesh element containing the element to get the boundary type for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user element number to get the boundary type for
    INTEGER(INTG), INTENT(OUT) :: onBoundary !<On return, the boundary type of the specified user element number.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalNumber

    ENTERS("MeshTopology_ElementOnBoundaryGet",err,error,*999)

    IF(.NOT.ASSOCIATED(meshElements)) CALL FlagError("Mesh elements is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(meshElements%elements)) CALL FlagError("Mesh elements elements is not associated.",err,error,*999)
    
    CALL MeshTopology_ElementGet(meshElements,userNumber,globalNumber,err,error,*999)
    IF(meshElements%elements(globalNumber)%boundary_Element) THEN
      onBoundary=MESH_ON_DOMAIN_BOUNDARY
    ELSE
      onBoundary=MESH_OFF_DOMAIN_BOUNDARY
    ENDIF
    
    EXITS("MeshTopology_ElementOnBoundaryGet")
    RETURN
999 onBoundary=MESH_OFF_DOMAIN_BOUNDARY
    ERRORSEXITS("MeshTopology_ElementOnBoundaryGet",err,error)    
    RETURN 1
   
  END SUBROUTINE MeshTopology_ElementOnBoundaryGet
  
  !
  !================================================================================================================================
  !

  !>Checks that a user node number exists in a mesh component. 
  SUBROUTINE MeshTopology_NodeCheckExistsMesh(mesh,meshComponentNumber,userNodeNumber,nodeExists,meshNodeNumber,err,error,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: mesh !<A pointer to the mesh to check the node exists on
    INTEGER(INTG), INTENT(IN) :: meshComponentNumber !<The mesh component to check the node exits on
    INTEGER(INTG), INTENT(IN) :: userNodeNumber !<The user node number to check if it exists
    LOGICAL, INTENT(OUT) :: nodeExists !<On exit, is .TRUE. if the node user number exists in the mesh component, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: meshNodeNumber !<On exit, if the node exists the mesh number corresponding to the user node number. If the node does not exist then mesh number will be 0.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalNodeNumber
    TYPE(MeshNodesType), POINTER :: meshNodes
    TYPE(NODES_TYPE), POINTER :: nodes
    TYPE(REGION_TYPE), POINTER :: region
    TYPE(TREE_NODE_TYPE), POINTER :: treeNode

    NULLIFY(meshNodes)
    NULLIFY(nodes)
    NULLIFY(region)
    
    ENTERS("MeshTopology_NodeCheckExistsMesh",err,error,*999)

    nodeExists=.FALSE.
    meshNodeNumber=0
    IF(ASSOCIATED(mesh)) THEN
      IF(mesh%MESH_FINISHED) THEN
        CALL MeshTopology_NodesGet(mesh,meshComponentNumber,meshNodes,err,error,*999)
        CALL Mesh_RegionGet(mesh,region,err,error,*999)
        nodes=>region%nodes
        IF(ASSOCIATED(nodes)) THEN
          CALL NODE_CHECK_EXISTS(nodes,userNodeNumber,nodeExists,globalNodeNumber,err,error,*999)
          NULLIFY(treeNode)
          CALL TREE_SEARCH(meshNodes%nodesTree,globalNodeNumber,treeNode,err,error,*999)
          IF(ASSOCIATED(treeNode)) THEN
            CALL TREE_NODE_VALUE_GET(meshNodes%nodesTree,treeNode,meshNodeNumber,err,error,*999)
            nodeExists=.TRUE.
          ENDIF
        ELSE
          CALL FlagError("Region nodes is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Mesh has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopology_NodeCheckExistsMesh")
    RETURN
999 ERRORSEXITS("MeshTopology_NodeCheckExistsMesh",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_NodeCheckExistsMesh
  
  !
  !================================================================================================================================
  !

  !>Checks that a user node number exists in a mesh nodes. 
  SUBROUTINE MeshTopology_NodeCheckExistsMeshNodes(meshNodes,userNodeNumber,nodeExists,meshNodeNumber,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: meshNodes !<A pointer to the mesh nodes to check the node exists on
    INTEGER(INTG), INTENT(IN) :: userNodeNumber !<The user node number to check if it exists
    LOGICAL, INTENT(OUT) :: nodeExists !<On exit, is .TRUE. if the node user number exists in the mesh nodes, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: meshNodeNumber !<On exit, if the node exists the mesh number corresponding to the user node number. If the node does not exist then mesh number will be 0.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalNodeNumber
    TYPE(MESH_TYPE), POINTER :: mesh
    TYPE(MeshComponentTopologyType), POINTER :: meshComponentTopology
    TYPE(NODES_TYPE), POINTER :: nodes
    TYPE(REGION_TYPE), POINTER :: region
    TYPE(TREE_NODE_TYPE), POINTER :: treeNode

    NULLIFY(nodes)
    NULLIFY(region)
    
    ENTERS("MeshTopology_NodeCheckExistsMeshNodes",err,error,*999)

    nodeExists=.FALSE.
    meshNodeNumber=0
    IF(ASSOCIATED(meshNodes)) THEN
      meshComponentTopology=>meshNodes%meshComponentTopology
      IF(ASSOCIATED(meshComponentTopology)) THEN
        mesh=>meshComponentTopology%mesh
        IF(ASSOCIATED(mesh)) THEN
          IF(mesh%MESH_FINISHED) THEN
            CALL Mesh_RegionGet(mesh,region,err,error,*999)
            CALL Region_NodesGet(region,nodes,err,error,*999)
            CALL NODE_CHECK_EXISTS(nodes,userNodeNumber,nodeExists,globalNodeNumber,err,error,*999)
            NULLIFY(treeNode)
            CALL TREE_SEARCH(meshNodes%nodesTree,globalNodeNumber,treeNode,err,error,*999)
            IF(ASSOCIATED(treeNode)) THEN
              CALL TREE_NODE_VALUE_GET(meshNodes%nodesTree,treeNode,meshNodeNumber,err,error,*999)
              nodeExists=.TRUE.
            ENDIF
         ELSE
            CALL FlagError("Mesh has not been finished.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Mesh component topology mesh is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Mesh nodes mesh component topology is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh nodes is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopology_NodeCheckExistsMeshNodes")
    RETURN
999 ERRORSEXITS("MeshTopology_NodeCheckExistsMeshNodes",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_NodeCheckExistsMeshNodes
  
  !
  !================================================================================================================================
  !

  !>Gets a mesh node number that corresponds to a user node number from mesh nodes. An error will be raised if the user node number does not exist.
  SUBROUTINE MeshTopology_NodeGetMeshNodes(meshNodes,userNodeNumber,meshNodeNumber,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: meshNodes !<A pointer to the mesh nodes to get the node from
    INTEGER(INTG), INTENT(IN) :: userNodeNumber !<The user node number to get
    INTEGER(INTG), INTENT(OUT) :: meshNodeNumber !<On exit, the mesh number corresponding to the user node number. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: meshComponentNumber
    LOGICAL :: nodeExists
    TYPE(MESH_TYPE), POINTER :: mesh
    TYPE(MeshComponentTopologyType), POINTER :: meshComponentTopology
    TYPE(VARYING_STRING) :: localError

    ENTERS("MeshTopology_NodeGetMeshNodes",err,error,*999)

    CALL MeshTopology_NodeCheckExistsMeshNodes(meshNodes,userNodeNumber,nodeExists,meshNodeNumber,err,error,*999)
    IF(.NOT.nodeExists) THEN
      meshComponentTopology=>meshNodes%meshComponentTopology
      IF(ASSOCIATED(meshComponentTopology)) THEN
        mesh=>meshComponentTopology%mesh
        IF(ASSOCIATED(mesh)) THEN
          meshComponentNumber=meshComponentTopology%meshComponentNumber
          localError="The user node number "//TRIM(NumberToVString(userNodeNumber,"*",err,error))// &
            & " does not exist in mesh component number "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))// &
            & " of mesh number "//TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ELSE
          localError="The user node number "//TRIM(NumberToVString(userNodeNumber,"*",err,error))// &
            & " does not exist in mesh component number "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        localError="The user node number "//TRIM(NumberToVString(userNodeNumber,"*",err,error))//" does not exist."
        CALL FlagError(localError,err,error,*999)        
      ENDIF
    ENDIF
    
    EXITS("MeshTopology_NodeGetMeshNodes")
    RETURN
999 meshNodeNumber=0
    ERRORSEXITS("MeshTopology_NodeGetMeshNodes",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_NodeGetMeshNodes
  
  !
  !================================================================================================================================
  !

  !>Finalises the given mesh topology node. 
  SUBROUTINE MeshTopology_NodeFinalise(node,err,error,*)

    !Argument variables
    TYPE(MeshNodeType) :: node !<The mesh node to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: derivativeIdx

    ENTERS("MeshTopology_NodeFinalise",err,error,*999)

    IF(ALLOCATED(node%derivatives)) THEN
      DO derivativeIdx=1,node%numberOfDerivatives
        CALL MeshTopology_NodeDerivativeFinalise(node%derivatives(derivativeIdx),err,error,*999)
      ENDDO !derivativeIdx
      DEALLOCATE(node%derivatives)
    ENDIF
    IF(ASSOCIATED(node%surroundingElements)) DEALLOCATE(node%surroundingElements)
  
    EXITS("MeshTopology_NodeFinalise")
    RETURN
999 ERRORSEXITS("MeshTopology_NodeFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_NodeFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the given mesh topology node.
  SUBROUTINE MeshTopology_NodeInitialise(node,err,error,*)

    !Argument variables
    TYPE(MeshNodeType) :: node !<The mesh node to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopology_NodeInitialise",err,error,*999)

    node%userNumber=0
    node%globalNumber=0
    node%numberOfSurroundingElements=0
    NULLIFY(node%surroundingElements)
    node%numberOfDerivatives=0
    node%boundaryNode=.FALSE.
    
    EXITS("MeshTopology_NodeInitialise")
    RETURN
999 ERRORSEXITS("MeshTopology_NodeInitialise",err,error)
    RETURN 1
  END SUBROUTINE MeshTopology_NodeInitialise

  !
  !================================================================================================================================
  !

  !>Calculates the nodes used the mesh identified by a given mesh topology.
  SUBROUTINE MeshTopology_NodesCalculate(topology,err,error,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: topology !<A pointer to the mesh topology
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,elementIdx,insertStatus,localNodeIdx,globalNode,meshNodeIdx,meshNode,numberOfNodes 
    INTEGER(INTG), POINTER :: globalNodeNumbers(:)
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(MESH_TYPE), POINTER :: mesh
    TYPE(MeshElementsType), POINTER :: elements
    TYPE(MeshNodesType), POINTER :: meshNodes
    TYPE(NODES_TYPE), POINTER :: nodes
    TYPE(TREE_TYPE), POINTER :: globalNodesTree
    TYPE(TREE_NODE_TYPE), POINTER :: treeNode
    TYPE(VARYING_STRING) :: dummyError,localError

    NULLIFY(globalNodeNumbers)
    NULLIFY(globalNodesTree)
    
    ENTERS("MeshTopology_NodesCalculate",err,error,*998)

    IF(ASSOCIATED(topology)) THEN
      elements=>topology%elements
      IF(ASSOCIATED(elements)) THEN
        meshNodes=>topology%nodes
        IF(ASSOCIATED(meshNodes)) THEN
          mesh=>topology%mesh
          IF(ASSOCIATED(mesh)) THEN
            NULLIFY(nodes)
            CALL MeshGlobalNodesGet(mesh,nodes,err,error,*999)
            IF(ALLOCATED(meshNodes%nodes)) THEN
              localError="Mesh number "//TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))// &
                & " already has allocated mesh topology nodes."
              CALL FlagError(localError,err,error,*998)
            ELSE
              !Work out what nodes are in the mesh
              CALL TREE_CREATE_START(globalNodesTree,err,error,*999)
              CALL TREE_INSERT_TYPE_SET(globalNodesTree,TREE_NO_DUPLICATES_ALLOWED,err,error,*999)
              CALL TREE_CREATE_FINISH(globalNodesTree,err,error,*999)
              DO elementIdx=1,elements%NUMBER_OF_ELEMENTS
                basis=>elements%elements(elementIdx)%basis
                DO localNodeIdx=1,basis%NUMBER_OF_NODES
                  globalNode=elements%elements(elementIdx)%GLOBAL_ELEMENT_NODES(localNodeIdx)
                  CALL TREE_ITEM_INSERT(globalNodesTree,globalNode,globalNode,insertStatus,err,error,*999)
                ENDDO !localNodeIdx
              ENDDO !elementIdx
              CALL TREE_DETACH_AND_DESTROY(globalNodesTree,numberOfNodes,globalNodeNumbers,err,error,*999)
              !Set up the mesh nodes.
              ALLOCATE(meshNodes%nodes(numberOfNodes),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate mesh topology nodes nodes.",err,error,*999)
              CALL TREE_CREATE_START(meshNodes%nodesTree,err,error,*999)
              CALL TREE_INSERT_TYPE_SET(meshNodes%nodesTree,TREE_NO_DUPLICATES_ALLOWED,err,error,*999)
              CALL TREE_CREATE_FINISH(meshNodes%nodesTree,err,error,*999) 
              DO meshNodeIdx=1,numberOfNodes
                CALL MeshTopology_NodeInitialise(meshNodes%nodes(meshNodeIdx),err,error,*999)
                meshNodes%nodes(meshNodeIdx)%meshNumber=meshNodeIdx
                meshNodes%nodes(meshNodeIdx)%globalNumber=globalNodeNumbers(meshNodeIdx)
                meshNodes%nodes(meshNodeIdx)%userNumber=nodes%nodes(globalNodeNumbers(meshNodeIdx))%USER_NUMBER
                CALL TREE_ITEM_INSERT(meshNodes%nodesTree,globalNodeNumbers(meshNodeIdx),meshNodeIdx,insertStatus,err,error,*999)
              ENDDO !nodeIdx
              meshNodes%numberOfNodes=numberOfNodes
              IF(ASSOCIATED(globalNodeNumbers)) DEALLOCATE(globalNodeNumbers)
              !Now recalculate the mesh element nodes
              DO elementIdx=1,elements%NUMBER_OF_ELEMENTS
                basis=>elements%elements(elementIdx)%basis
                ALLOCATE(elements%elements(elementIdx)%MESH_ELEMENT_NODES(basis%NUMBER_OF_NODES),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate mesh topology elements mesh element nodes.",err,error,*999)
                DO localNodeIdx=1,basis%NUMBER_OF_NODES
                  globalNode=elements%elements(elementIdx)%GLOBAL_ELEMENT_NODES(localNodeIdx)
                  NULLIFY(treeNode)
                  CALL TREE_SEARCH(meshNodes%nodesTree,globalNode,treeNode,err,error,*999)
                  IF(ASSOCIATED(treeNode)) THEN
                    CALL TREE_NODE_VALUE_GET(meshNodes%nodesTree,treeNode,meshNode,err,error,*999)
                    elements%elements(elementIdx)%MESH_ELEMENT_NODES(localNodeIdx)=meshNode
                  ELSE
                    localError="Could not find global node "//TRIM(NumberToVString(globalNode,"*",err,error))//" (user node "// &
                      & TRIM(NumberToVString(nodes%nodes(globalNode)%USER_NUMBER,"*",err,error))//") in the mesh nodes."
                    CALL FlagError(localError,err,error,*999)
                  ENDIF
                ENDDO !localNodeIdx
              ENDDO !elementIdx
            ENDIF
          ELSE
            CALL FlagError("Mesh topology mesh is not associated.",err,error,*998)
          ENDIF
        ELSE
          CALL FlagError("Mesh topology nodes is not associated.",err,error,*998)
        ENDIF
      ELSE
        CALL FlagError("Mesh topology elements is not associated.",err,error,*998)
      ENDIF
    ELSE
      CALL FlagError("Mesh topology is not associated.",err,error,*998)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number of mesh nodes = ",meshNodes%numberOfNodes,err,error,*999)
      DO meshNodeIdx=1,meshNodes%numberOfNodes
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Mesh node number = ",meshNodeIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Global node number = ",meshNodes%nodes(meshNodeIdx)%globalNumber, &
          & err,error,*999)        
      ENDDO !meshNodeIdx
    ENDIF
    
    EXITS("MeshTopology_NodesCalculate")
    RETURN
999 IF(ASSOCIATED(globalNodeNumbers)) DEALLOCATE(globalNodeNumbers)
    IF(ASSOCIATED(globalNodesTree)) CALL TREE_DESTROY(globalNodesTree,dummyErr,dummyError,*998)
998 ERRORSEXITS("MeshTopology_NodesCalculate",err,error)
    RETURN 1
   
  END SUBROUTINE MeshTopology_NodesCalculate

  !
  !================================================================================================================================
  !

  !>Destroys the nodes in a mesh topology.
  SUBROUTINE MeshTopology_NodesDestroy(nodes,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: nodes !<A pointer to the mesh nodes to destroy 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopology_NodesDestroy",err,error,*999)

    IF(ASSOCIATED(nodes)) THEN
      CALL MeshTopology_NodesFinalise(nodes,err,error,*999)
    ELSE
      CALL FlagError("Mesh topology nodes is not associated",err,error,*999)
    ENDIF

    EXITS("MeshTopology_NodesDestroy")
    RETURN
999 ERRORSEXITS("MeshTopology_NodesDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_NodesDestroy
  
  !
  !================================================================================================================================
  !

  !>Finalises the given mesh topology node. 
  SUBROUTINE MeshTopology_NodeDerivativeFinalise(nodeDerivative,err,error,*)

    !Argument variables
    TYPE(MeshNodeDerivativeType) :: nodeDerivative !<The mesh node derivative to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopology_NodeDerivativeFinalise",err,error,*999)

    IF(ALLOCATED(nodeDerivative%userVersionNumbers)) DEALLOCATE(nodeDerivative%userVersionNumbers)
    IF(ALLOCATED(nodeDerivative%dofIndex)) DEALLOCATE(nodeDerivative%dofIndex)

    EXITS("MeshTopology_NodeDerivativeFinalise")
    RETURN
999 ERRORSEXITS("MeshTopology_NodeDerivativeFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_NodeDerivativeFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the given mesh topology node.
  SUBROUTINE MeshTopology_NodeDerivativeInitialise(nodeDerivative,err,error,*)

    !Argument variables
    TYPE(MeshNodeDerivativeType) :: nodeDerivative !<The mesh node derivative to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopology_NodeDerivativeInitialise",err,error,*999)

    nodeDerivative%numberOfVersions=0
    nodeDerivative%globalDerivativeIndex=0
    nodederivative%partialDerivativeIndex=0

    EXITS("MeshTopology_NodeDerivativeInitialise")
    RETURN
999 ERRORSEXITS("MeshTopology_NodeDerivativeInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_NodeDerivativeInitialise

  !
  !================================================================================================================================
  !

  !>Calculates the number of derivatives at each node in a topology.
  SUBROUTINE MeshTopology_NodesDerivativesCalculate(topology,err,error,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: topology !<A pointer to the mesh topology to calculate the derivates at each node for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: derivativeIdx,element,elementIdx,globalDerivative,localNodeIdx,maxNumberOfDerivatives,nodeIdx, &
      & numberOfDerivatives
    INTEGER(INTG), ALLOCATABLE :: derivatives(:)
    LOGICAL :: found
    TYPE(LIST_TYPE), POINTER :: nodeDerivativeList
    TYPE(MeshElementsType), POINTER :: elements
    TYPE(MeshNodesType), POINTER :: nodes
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("MeshTopology_NodesDerivativesCalculate",err,ERROR,*999)

     IF(ASSOCIATED(topology)) THEN
       elements=>topology%elements
       IF(ASSOCIATED(elements)) THEN
         nodes=>topology%nodes
         IF(ASSOCIATED(nodes)) THEN
          !Loop over the mesh nodes
          DO nodeIdx=1,nodes%numberOfNodes
            !Calculate the number of derivatives and versions at each node. This needs to be calculated by looking at the
            !mesh elements as we may have an adjacent element in another domain with a higher order basis also with versions.
            NULLIFY(nodeDerivativeList)
            CALL LIST_CREATE_START(nodeDerivativeList,err,error,*999)
            CALL LIST_DATA_TYPE_SET(nodeDerivativeList,LIST_INTG_TYPE,err,error,*999)
            CALL LIST_INITIAL_SIZE_SET(nodeDerivativeList,8,err,error,*999)
            CALL LIST_CREATE_FINISH(nodeDerivativeList,err,error,*999)
            maxNumberOfDerivatives=-1
            DO elementIdx=1,nodes%nodes(nodeIdx)%numberOfSurroundingElements
              element=nodes%nodes(nodeIdx)%surroundingElements(elementIdx)
              basis=>elements%elements(element)%basis
              !Find the local node corresponding to this node
              found=.FALSE.
              DO localNodeIdx=1,basis%NUMBER_OF_NODES
                IF(elements%elements(element)%MESH_ELEMENT_NODES(localNodeIdx)==nodeIdx) THEN
                  found=.TRUE.
                  EXIT
                ENDIF
              ENDDO !nn
              IF(found) THEN
                DO derivativeIdx=1,basis%NUMBER_OF_DERIVATIVES(localNodeIdx)
                  CALL LIST_ITEM_ADD(nodeDerivativeList,basis%PARTIAL_DERIVATIVE_INDEX(derivativeIdx,localNodeIdx),err,error,*999)
                ENDDO !derivativeIdx
                IF(basis%NUMBER_OF_DERIVATIVES(localNodeIdx)>maxNumberOfDerivatives) &
                  & maxNumberOfDerivatives=basis%NUMBER_OF_DERIVATIVES(localNodeidx)
              ELSE
                CALL FlagError("Could not find local node.",err,error,*999)
              ENDIF
            ENDDO !elem_idx
            CALL LIST_REMOVE_DUPLICATES(nodeDerivativeList,err,error,*999)
            CALL LIST_DETACH_AND_DESTROY(nodeDerivativeList,numberOfDerivatives,derivatives,err,error,*999)
            IF(numberOfDerivatives==maxNumberOfDerivatives) THEN
              !Set up the node derivatives.
              ALLOCATE(nodes%nodes(nodeIdx)%derivatives(maxNumberOfDerivatives),STAT=err)
              nodes%nodes(nodeIdx)%numberOfDerivatives=maxNumberOfDerivatives
              DO derivativeIdx=1,numberOfDerivatives
                CALL MeshTopology_NodeDerivativeInitialise(nodes%nodes(nodeIdx)%derivatives(derivativeIdx),err,error,*999)
                nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%partialDerivativeIndex = derivatives(derivativeIdx)
                globalDerivative=PARTIAL_DERIVATIVE_GLOBAL_DERIVATIVE_MAP(derivatives(derivativeIdx))
                IF(globalDerivative/=0) THEN
                  nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%globalDerivativeIndex=globalDerivative
                ELSE
                  localError="The partial derivative index of "//TRIM(NumberToVstring(derivatives(derivativeIdx),"*", &
                    & err,error))//" for derivative number "//TRIM(NumberToVstring(derivativeIdx,"*",err,error))// &
                    & " does not have a corresponding global derivative."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
              ENDDO !derivativeIdx
              DEALLOCATE(derivatives)
            ELSE
              localError="Invalid mesh configuration. User node "// &
                & TRIM(NumberToVstring(nodes%nodes(nodeIdx)%userNumber,"*",err,error))// &
                & " has inconsistent derivative directions."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ENDDO !nodeIdx
        ELSE
          CALL FlagError("Mesh topology nodes is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Mesh topology elements is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh topology is not associated.",err,error,*999)
    ENDIF
   
    IF(diagnostics1) THEN
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number of mesh global nodes = ",nodes%numberOfNodes,err,error,*999)
      DO nodeIdx=1,nodes%numberOfNodes
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Mesh global node number = ",nodeIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of derivatives = ",nodes%nodes(nodeIdx)%numberOfDerivatives, &
          & err,error,*999)
        DO derivativeIdx=1,nodes%nodes(nodeIdx)%numberOfDerivatives
          !TODO: change output string below so that it writes out derivativeIdx index as well
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Global derivative index(derivativeIdx) = ", &
            & nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%globalDerivativeIndex,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Partial derivative index(derivativeIdx) = ", &
            & nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%partialDerivativeIndex,err,error,*999)
        ENDDO !derivativeIdx
      ENDDO !node_idx
    ENDIF
    
    EXITS("MeshTopology_NodesDerivativesCalculate")
    RETURN
999 IF(ALLOCATED(derivatives)) DEALLOCATE(derivatives)
    IF(ASSOCIATED(nodeDerivativeList)) CALL LIST_DESTROY(nodeDerivativeList,err,error,*998)
998 ERRORSEXITS("MeshTopology_NodesDerivativesCalculate",err,error)
    RETURN 1
   
  END SUBROUTINE MeshTopology_NodesDerivativesCalculate

  !
  !================================================================================================================================
  !
  
  !>Returns the number of derivatives for a node in a mesh
  SUBROUTINE MeshTopology_NodeNumberOfDerivativesGet(meshNodes,userNumber,numberOfDerivatives,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: meshNodes !<A pointer to the mesh nodes containing the nodes to get the number of derivatives for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user node number to get the number of derivatives for
    INTEGER(INTG), INTENT(OUT) :: numberOfDerivatives !<On return, the number of global derivatives at the specified user node number.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: meshComponentNumber,meshNumber
    LOGICAL :: nodeExists
    TYPE(MESH_TYPE), POINTER :: mesh
    TYPE(MeshComponentTopologyType), POINTER :: meshComponentTopology
    TYPE(VARYING_STRING) :: localError

    ENTERS("MeshTopology_NodeNumberOfDerivativesGet",err,error,*999)

    IF(ASSOCIATED(meshNodes)) THEN
      CALL MeshTopology_NodeCheckExists(meshNodes,userNumber,nodeExists,meshNumber,err,error,*999)
      IF(nodeExists) THEN
        numberOfDerivatives=meshNodes%nodes(meshNumber)%numberOfDerivatives
      ELSE
        meshComponentTopology=>meshNodes%meshComponentTopology
        IF(ASSOCIATED(meshComponentTopology)) THEN
          mesh=>meshComponentTopology%mesh
          IF(ASSOCIATED(mesh)) THEN
            meshComponentNumber=meshComponentTopology%meshComponentNumber
            localError="The user node number "//TRIM(NumberToVString(userNumber,"*",err,error))// &
              & " does not exist in mesh component number "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))// &
              & " of mesh number "//TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ELSE
            CALL FlagError("Mesh component topology mesh is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Mesh nodes mesh component topology is not associated.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Mesh nodes is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopology_NodeNumberOfDerivativesGet")
    RETURN
999 ERRORSEXITS("MeshTopology_NodeNumberOfDerivativesGet",err,error)    
    RETURN 1
   
  END SUBROUTINE MeshTopology_NodeNumberOfDerivativesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the global derivative numbers for a node in mesh nodes
  SUBROUTINE MeshTopology_NodeDerivativesGet(meshNodes,userNumber,derivatives,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: meshNodes !<A pointer to the mesh nodes containing the node to get the derivatives for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the node to get the derivatives for
    INTEGER(INTG), INTENT(OUT) :: derivatives(:) !<On return, the global derivatives at the specified node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: derivativeIdx,meshComponentNumber,meshNumber,numberOfDerivatives
    LOGICAL :: nodeExists
    TYPE(MESH_TYPE), POINTER :: mesh
    TYPE(MeshComponentTopologyType), POINTER :: meshComponentTopology
    TYPE(VARYING_STRING) :: localError

    ENTERS("MeshTopology_NodeDerivativesGet",err,error,*999)

    IF(ASSOCIATED(meshNodes)) THEN
      CALL MeshTopology_NodeCheckExists(meshNodes,userNumber,nodeExists,meshNumber,err,error,*999)
      IF(nodeExists) THEN         
        numberOfDerivatives=meshNodes%nodes(meshNumber)%numberOfDerivatives
        IF(SIZE(derivatives,1)>=numberOfDerivatives) THEN
          DO derivativeIdx=1,numberOfDerivatives
            derivatives(derivativeIdx)=meshNodes%nodes(meshNumber)%derivatives(derivativeIdx)%globalDerivativeIndex
          ENDDO !derivativeIdx
        ELSE
          localError="The size of the supplied derivatives array of "// &
            & TRIM(NumberToVString(SIZE(derivatives,1),"*",err,error))// &
            & " is too small. The size should be >= "// &
            & TRIM(NumberToVString(numberOfDerivatives,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        meshComponentTopology=>meshNodes%meshComponentTopology
        IF(ASSOCIATED(meshComponentTopology)) THEN
          mesh=>meshComponentTopology%mesh
          IF(ASSOCIATED(mesh)) THEN
            meshComponentNumber=meshComponentTopology%meshComponentNumber
            localError="The user node number "//TRIM(NumberToVString(userNumber,"*",err,error))// &
              & " does not exist in mesh component number "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))// &
              & " of mesh number "//TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ELSE
            CALL FlagError("Mesh component topology mesh is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Mesh nodes mesh component topology is not associated.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Mesh nodes is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopology_NodeDerivativesGet")
    RETURN
999 ERRORSEXITS("MeshTopology_NodeDerivativesGet",err,error)    
    RETURN 1
   
  END SUBROUTINE MeshTopology_NodeDerivativesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of versions for a derivative of a node in mesh nodes
  SUBROUTINE MeshTopology_NodeNumberOfVersionsGet(meshNodes,derivativeNumber,userNumber,numberOfVersions,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: meshNodes !<A pointer to the mesh nodes containing the node to get the number of versions for
    INTEGER(INTG), INTENT(IN) :: derivativeNumber !<The derivative number of the node to get the number of versions for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the node to get the number of versions for
    INTEGER(INTG), INTENT(OUT) :: numberOfVersions !<On return, the number of versions for the specified derivative of the specified node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: meshComponentNumber,meshNumber
    LOGICAL :: nodeExists
    TYPE(MESH_TYPE), POINTER :: mesh
    TYPE(MeshComponentTopologyType), POINTER :: meshComponentTopology
    TYPE(VARYING_STRING) :: localError

    ENTERS("MeshTopology_NodeNumberOfVersionsGet",err,error,*999)

    IF(ASSOCIATED(meshNodes)) THEN
      CALL MeshTopology_NodeCheckExists(meshNodes,userNumber,nodeExists,meshNumber,err,error,*999)
      IF(nodeExists) THEN
        IF(derivativeNumber>=1.AND.derivativeNumber<=meshNodes%nodes(meshNumber)%numberOfDerivatives) THEN
          numberOfVersions=meshNodes%nodes(meshNumber)%derivatives(derivativeNumber)%numberOfVersions
        ELSE
          localError="The specified derivative index of "// &
            & TRIM(NumberToVString(derivativeNumber,"*",err,error))// &
            & " is invalid. The derivative index must be >= 1 and <= "// &
            & TRIM(NumberToVString(meshNodes%nodes(meshNumber)%numberOfDerivatives,"*",err,error))// &
            & " for user node number "//TRIM(NumberToVString(userNumber,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        meshComponentTopology=>meshNodes%meshComponentTopology
        IF(ASSOCIATED(meshComponentTopology)) THEN
          mesh=>meshComponentTopology%mesh
          IF(ASSOCIATED(mesh)) THEN
            meshComponentNumber=meshComponentTopology%meshComponentNumber
            localError="The user node number "//TRIM(NumberToVString(userNumber,"*",err,error))// &
              & " does not exist in mesh component number "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))// &
              & " of mesh number "//TRIM(NumberToVString(mesh%USER_NUMBER,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ELSE
            CALL FlagError("Mesh component topology mesh is not associated.",err,error,*999)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Mesh nodes is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopology_NodeNumberOfVersionsGet")
    RETURN
999 ERRORSEXITS("MeshTopology_NodeNumberOfVersionsGet",err,error)    
    RETURN 1
   
  END SUBROUTINE MeshTopology_NodeNumberOfVersionsGet
  
  !
  !================================================================================================================================
  !

  !>Returns the number of nodes for a node in a mesh
  SUBROUTINE MeshTopology_NodesNumberOfNodesGet(meshNodes,numberOfNodes,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: meshNodes !<A pointer to the mesh nodes containing the nodes to get the number of nodes for
    INTEGER(INTG), INTENT(OUT) :: numberOfNodes !<On return, the number of nodes in the mesh.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
  
    ENTERS("MeshTopology_NodesNumberOfNodesGet",err,error,*999)

    IF(ASSOCIATED(meshNodes)) THEN
      numberOfNodes=meshNodes%numberOfNodes
    ELSE
      CALL FlagError("Mesh nodes is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopology_NodesNumberOfNodesGet")
    RETURN
999 ERRORSEXITS("MeshTopology_NodesNumberOfNodesGet",err,error)    
    RETURN 1
   
  END SUBROUTINE MeshTopology_NodesNumberOfNodesGet
  
  !
  !================================================================================================================================
  !

  !>Calculates the number of versions at each node in a topology.
  SUBROUTINE MeshTopology_NodesVersionCalculate(topology,err,error,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: topology !<A pointer to the mesh topology to calculate the versions at each node for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: element,localNodeIdx,derivativeIdx,nodeIdx,numberOfVersions,versionIdx
    INTEGER(INTG), ALLOCATABLE :: versions(:)
    TYPE(LIST_PTR_TYPE), POINTER :: nodeVersionList(:,:)
    TYPE(MeshElementsType), POINTER :: elements
    TYPE(MeshNodesType), POINTER :: nodes
    TYPE(BASIS_TYPE), POINTER :: basis
    
    ENTERS("MeshTopology_NodesVersionCalculate",err,error,*999)

    IF(ASSOCIATED(topology)) THEN
      elements=>topology%elements
      IF(ASSOCIATED(elements)) THEN
        nodes=>topology%nodes
        IF(ASSOCIATED(nodes)) THEN
          !Loop over the mesh elements
          !Calculate the number of versions at each node. This needs to be calculated by looking at all the mesh elements
          !as we may have an adjacent elements in another domain with a higher order basis along with different versions
          !being assigned to its derivatives.
          !\todo : See if there are any constraints that can be applied to restrict the amount of memory being allocated here
          ALLOCATE(nodeVersionList(MAXIMUM_GLOBAL_DERIV_NUMBER,nodes%numberOfNodes),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate node version list.",err,error,*999)
          DO nodeIdx=1,nodes%numberOfNodes
            DO derivativeIdx=1,nodes%nodes(nodeIdx)%numberOfDerivatives
              NULLIFY(nodeVersionList(derivativeIdx,nodeIdx)%ptr)
              CALL LIST_CREATE_START(nodeVersionList(derivativeIdx,nodeIdx)%ptr,err,error,*999)
              CALL LIST_DATA_TYPE_SET(nodeVersionList(derivativeIdx,nodeIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
              CALL LIST_INITIAL_SIZE_SET(nodeVersionList(derivativeIdx,nodeIdx)%ptr,8,err,error,*999)
              CALL LIST_CREATE_FINISH(nodeVersionList(derivativeIdx,nodeIdx)%ptr,err,error,*999)
            ENDDO!derivativeIdx
          ENDDO!nodeIdx
          DO element=1,elements%NUMBER_OF_ELEMENTS
            basis=>elements%elements(element)%basis
            DO localNodeIdx=1,BASIS%NUMBER_OF_NODES
              DO derivativeIdx=1,basis%NUMBER_OF_DERIVATIVES(localNodeIdx)
                CALL LIST_ITEM_ADD(nodeVersionList(derivativeIdx,elements%elements(element)% &
                  & MESH_ELEMENT_NODES(localNodeIdx))%ptr,elements%elements(element)%USER_ELEMENT_NODE_VERSIONS( &
                  & derivativeIdx,localNodeIdx),err,error,*999)
              ENDDO!derivativeIdx
            ENDDO!localNodeIdx
          ENDDO!element
          DO nodeIdx=1,nodes%numberOfNodes
            DO derivativeIdx=1,nodes%nodes(nodeIdx)%numberOfDerivatives
              CALL LIST_REMOVE_DUPLICATES(nodeVersionList(derivativeIdx,nodeIdx)%ptr,err,error,*999)
              CALL LIST_DETACH_AND_DESTROY(nodeVersionList(derivativeIdx,nodeIdx)%ptr,numberOfVersions,versions, &
                & err,error,*999)
              nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions = MAXVAL(versions(1:numberOfVersions))
              ALLOCATE(nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%userVersionNumbers(nodes%nodes(nodeIdx)% &
                & derivatives(derivativeIdx)%numberOfVersions),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate node global derivative index.",err,error,*999)
              DO versionIdx=1,nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions 
                nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%userVersionNumbers(versionIdx) = versionIdx
              ENDDO !versionIdx
              DEALLOCATE(versions)
            ENDDO!derivativeIdx
          ENDDO!nodeIdx
          DEALLOCATE(nodeVersionList)
          NULLIFY(nodeVersionList)
        ELSE
          CALL FlagError("Mesh topology nodes is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Mesh topology elements is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh topology is not associated.",err,error,*999)
    ENDIF
   
    IF(diagnostics1) THEN
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number of mesh global nodes = ",nodes%numberOfNodes,err,error,*999)
      DO nodeIdx=1,nodes%numberOfNodes
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Mesh global node number = ",nodeIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of derivatives = ", &
          & nodes%nodes(nodeIdx)%numberOfDerivatives,err,error,*999)
        DO derivativeIdx=1,nodes%nodes(nodeIdx)%numberOfDerivatives
          !\todo : change output string below so that it writes out derivativeIdx index as well
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Global derivative index(derivativeIdx) = ", &
            & nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%globalDerivativeIndex,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Partial derivative index(derivativeIdx) = ", &
            & nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%partialDerivativeIndex,err,error,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
            & nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions,8,8, &
            & nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%userVersionNumbers, & 
            & '("    User Version index(derivativeIdx,:) :",8(X,I2))','(36X,8(X,I2))',err,error,*999)
        ENDDO!derivativeIdx
      ENDDO !nodeIdx
    ENDIF
    
    EXITS("MeshTopology_NodesVersionCalculate")
    RETURN
999 IF(ALLOCATED(versions)) DEALLOCATE(versions)
    IF(ASSOCIATED(nodeVersionList)) THEN
      DO nodeIdx=1,nodes%numberOfNodes
        DO derivativeIdx=1,nodes%nodes(nodeIdx)%numberOfDerivatives
          CALL LIST_DESTROY(nodeVersionList(derivativeIdx,nodeIdx)%ptr,err,error,*998)
        ENDDO !derivativeIdx
      ENDDO !nodeIdx
      DEALLOCATE(nodeVersionList)
    ENDIF
998 ERRORSEXITS("MeshTopology_NodesVersionCalculate",err,error)
    RETURN 1
   
  END SUBROUTINE MeshTopology_NodesVersionCalculate

  !
  !================================================================================================================================
  !
  
  !>Returns if the node in a mesh is on the boundary or not
  SUBROUTINE MeshTopology_NodeOnBoundaryGet(meshNodes,userNumber,onBoundary,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: meshNodes !<A pointer to the mesh nodes containing the nodes to get the boundary type for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user node number to get the boundary type for
    INTEGER(INTG), INTENT(OUT) :: onBoundary !<On return, the boundary type of the specified user node number.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: meshNumber
 
    ENTERS("MeshTopology_NodeOnBoundaryGet",err,error,*999)

    IF(.NOT.ASSOCIATED(meshNodes)) CALL FlagError("Mesh nodes is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(meshNodes%nodes)) CALL FlagError("Mesh nodes nodes is not allocated.",err,error,*999)
    
    CALL MeshTopology_NodeGet(meshNodes,userNumber,meshNumber,err,error,*999)
    IF(meshNodes%nodes(meshNumber)%boundaryNode) THEN
      onBoundary=MESH_ON_DOMAIN_BOUNDARY
    ELSE
      onBoundary=MESH_OFF_DOMAIN_BOUNDARY
    ENDIF
    
    EXITS("MeshTopology_NodeOnBoundaryGet")
    RETURN
999 ERRORSEXITS("MeshTopology_NodeOnBoundaryGet",err,error)    
    RETURN 1
   
  END SUBROUTINE MeshTopology_NodeOnBoundaryGet
  
  !
  !================================================================================================================================
  !

  !>Calculates the element numbers surrounding a node for a mesh.
  SUBROUTINE MeshTopology_SurroundingElementsCalculate(topology,err,error,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: topology !<A pointer to the mesh topology to calculate the elements surrounding each node for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: element,elementIdx,insertPosition,localNodeIdx,node,surroundingElementNumber
    INTEGER(INTG), POINTER :: newSurroundingElements(:)
    LOGICAL :: foundElement
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(MeshElementsType), POINTER :: elements
    TYPE(MeshNodesType), POINTER :: nodes

    NULLIFY(newSurroundingElements)

    ENTERS("MeshTopology_SurroundingElementsCalculate",err,error,*999)
    
    IF(ASSOCIATED(topology)) THEN
      elements=>topology%elements      
      IF(ASSOCIATED(elements)) THEN
        nodes=>topology%nodes       
        IF(ASSOCIATED(nodes)) THEN
          IF(ALLOCATED(nodes%nodes)) THEN
            DO elementIdx=1,elements%NUMBER_OF_ELEMENTS
              basis=>elements%elements(elementIdx)%basis
              DO localNodeIdx=1,basis%NUMBER_OF_NODES
                node=elements%elements(elementIdx)%MESH_ELEMENT_NODES(localNodeIdx)
                foundElement=.FALSE.
                element=1
                insertPosition=1
                DO WHILE(element<=nodes%nodes(node)%numberOfSurroundingElements.AND..NOT.foundElement)
                  surroundingElementNumber=nodes%nodes(node)%surroundingElements(element)
                  IF(surroundingElementNumber==elementIdx) THEN
                    foundElement=.TRUE.
                  ENDIF
                  element=element+1
                  IF(elementIdx>=surroundingElementNumber) THEN
                    insertPosition=element
                  ENDIF
                ENDDO
                IF(.NOT.foundElement) THEN
                  !Insert element into surrounding elements
                  ALLOCATE(newSurroundingElements(nodes%nodes(node)%numberOfSurroundingElements+1),STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate new surrounding elements.",err,error,*999)
                  IF(ASSOCIATED(nodes%nodes(node)%surroundingElements)) THEN
                    newSurroundingElements(1:insertPosition-1)=nodes%nodes(node)%surroundingElements(1:insertPosition-1)
                    newSurroundingElements(insertPosition)=elementIdx
                    newSurroundingElements(insertPosition+1:nodes%nodes(node)%numberOfSurroundingElements+1)= &
                      & nodes%nodes(node)%surroundingElements(insertPosition:nodes%nodes(node)%numberOfSurroundingElements)
                    DEALLOCATE(nodes%nodes(node)%surroundingElements)
                  ELSE
                    newSurroundingElements(1)=elementIdx
                  ENDIF
                  nodes%nodes(node)%surroundingElements=>newSurroundingElements
                  nodes%nodes(node)%numberOfSurroundingElements=nodes%nodes(node)%numberOfSurroundingElements+1
                ENDIF
              ENDDO !localNodeIdx
            ENDDO !elementIdx
          ELSE
            CALL FlagError("Mesh topology nodes nodes have not been allocated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Mesh topology nodes are not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Mesh topology elements is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh topology not associated.",err,error,*999)
    ENDIF

    EXITS("MeshTopology_SurroundingElementsCalculate")
    RETURN
999 IF(ASSOCIATED(newSurroundingElements)) DEALLOCATE(newSurroundingElements)
    ERRORSEXITS("MeshTopology_SurroundingElementsCalculate",err,error)
    RETURN 1   
  END SUBROUTINE MeshTopology_SurroundingElementsCalculate
  
  !
  !===============================================================================================================================
  !

  !>Finalises the nodes data structures for a mesh topology and deallocates any memory.
  SUBROUTINE MeshTopology_NodesFinalise(nodes,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: nodes !<A pointer to the mesh topology nodes to finalise the nodes for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: nodeIdx

    ENTERS("MeshTopology_NodesFinalise",err,error,*999)

    IF(ASSOCIATED(nodes)) THEN
      IF(ALLOCATED(nodes%nodes)) THEN
        DO nodeIdx=1,SIZE(nodes%nodes,1)
          CALL MeshTopology_NodeFinalise(nodes%nodes(nodeIdx),err,error,*999)
        ENDDO !nodesIdx
        DEALLOCATE(nodes%nodes)
      ENDIF
      IF(ASSOCIATED(nodes%nodesTree)) CALL TREE_DESTROY(nodes%nodesTree,err,error,*999)
      DEALLOCATE(nodes)
    ENDIF
 
    EXITS("MeshTopology_NodesFinalise")
    RETURN
999 ERRORSEXITS("MeshTopology_NodesFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_NodesFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the nodes in a given mesh topology. \todo finalise on errors
  SUBROUTINE MeshTopology_NodesInitialise(topology,err,error,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: topology !<A pointer to the mesh topology to initialise the nodes for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopology_NodesInitialise",err,error,*999)

    IF(ASSOCIATED(topology)) THEN
      IF(ASSOCIATED(topology%nodes)) THEN
        CALL FlagError("Mesh already has topology nodes associated.",err,error,*999)
      ELSE
        ALLOCATE(topology%nodes,STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate topology nodes.",err,error,*999)
        topology%nodes%numberOfNodes=0
        topology%nodes%meshComponentTopology=>topology
        NULLIFY(topology%nodes%nodesTree)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopology_NodesInitialise")
    RETURN
999 ERRORSEXITS("MeshTopology_NodesInitialise",err,error)
    RETURN 1
  END SUBROUTINE MeshTopology_NodesInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the meshes and deallocates all memory
  SUBROUTINE MESHES_FINALISE(MESHES,ERR,ERROR,*)

   !Argument variables
    TYPE(MESHES_TYPE), POINTER :: MESHES !<A pointer to the meshes to finalise the meshes for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(MESH_TYPE), POINTER :: MESH
 
    ENTERS("MESHES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESHES)) THEN
      DO WHILE(MESHES%NUMBER_OF_MESHES>0)
        MESH=>MESHES%MESHES(1)%PTR
        CALL MESH_DESTROY(MESH,ERR,ERROR,*999)
      ENDDO !mesh_idx
      DEALLOCATE(MESHES)
    ELSE
      CALL FlagError("Meshes is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("MESHES_FINALISE")
    RETURN
999 ERRORSEXITS("MESHES_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE MESHES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the generic meshes.
  SUBROUTINE MESHES_INITIALISE_GENERIC(MESHES,ERR,ERROR,*)

    !Argument variables
    TYPE(MESHES_TYPE), POINTER :: MESHES !<A pointer to the meshes to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    ENTERS("MESHES_INITIALISE_GENERIC",ERR,ERROR,*998)

    IF(ASSOCIATED(MESHES)) THEN
      CALL FlagError("Meshes is already associated.",ERR,ERROR,*998)
    ELSE
      ALLOCATE(MESHES,STAT=ERR)
      IF(ERR/=0) CALL FlagError("Meshes could not be allocated",ERR,ERROR,*999)
      NULLIFY(MESHES%REGION)
      NULLIFY(MESHES%INTERFACE)
      MESHES%NUMBER_OF_MESHES=0
      NULLIFY(MESHES%MESHES)
    ENDIF
    
    EXITS("MESHES_INITIALISE_GENERIC")
    RETURN
999 CALL MESHES_FINALISE(MESHES,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("MESHES_INITIALISE_GENERIC",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESHES_INITIALISE_GENERIC

  !
  !================================================================================================================================
  !

  !>Initialises the meshes for the given interface.
  SUBROUTINE MESHES_INITIALISE_INTERFACE(INTERFACE,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to initialise the meshes for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESHES_INITIALISE_INTERFACE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(INTERFACE%MESHES)) THEN
        LOCAL_ERROR="Interface number "//TRIM(NumberToVString(INTERFACE%USER_NUMBER,"*",ERR,ERROR))// &
          & " already has a mesh associated"
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        CALL MESHES_INITIALISE_GENERIC(INTERFACE%MESHES,ERR,ERROR,*999)
        INTERFACE%MESHES%INTERFACE=>INTERFACE
      ENDIF
    ELSE
      CALL FlagError("Interface is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESHES_INITIALISE_INTERFACE")
    RETURN
999 ERRORSEXITS("MESHES_INITIALISE_INTERFACE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESHES_INITIALISE_INTERFACE

  !
  !================================================================================================================================
  !

  !>Initialises the meshes for the given region.
  SUBROUTINE MESHES_INITIALISE_REGION(REGION,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to initialise the meshes for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESHES_INITIALISE_REGION",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%MESHES)) THEN
        LOCAL_ERROR="Region number "//TRIM(NumberToVString(REGION%USER_NUMBER,"*",ERR,ERROR))// &
          & " already has a mesh associated"
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        CALL MESHES_INITIALISE_GENERIC(REGION%MESHES,ERR,ERROR,*999)
        REGION%MESHES%REGION=>REGION
      ENDIF
    ELSE
      CALL FlagError("Region is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESHES_INITIALISE_REGION")
    RETURN
999 ERRORSEXITS("MESHES_INITIALISE_REGION",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESHES_INITIALISE_REGION

  !
  !================================================================================================================================
  !

  !>Gets the domain for a given node in a decomposition of a mesh. \todo should be able to specify lists of elements. \see OPENCMISS::Iron::cmfe_DecompositionNodeDomainGet
  SUBROUTINE DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,USER_NODE_NUMBER,MESH_COMPONENT_NUMBER,DOMAIN_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to set the element domain for
    INTEGER(INTG), INTENT(IN) :: USER_NODE_NUMBER !<The global element number to set the domain for.
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT_NUMBER !<The mesh component number to set the domain for.
    INTEGER(INTG), INTENT(OUT) :: DOMAIN_NUMBER !<On return, the domain of the global element.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(MeshComponentTopologyType), POINTER :: MESH_TOPOLOGY
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: GLOBAL_NODE_NUMBER
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE
    TYPE(MeshNodesType), POINTER :: MESH_NODES
    TYPE(DOMAIN_TYPE), POINTER :: MESH_DOMAIN

    ENTERS("DECOMPOSITION_NODE_DOMAIN_GET",ERR,ERROR,*999)

!!TODO: interface should specify user element number ???
    GLOBAL_NODE_NUMBER=0
    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        MESH=>DECOMPOSITION%MESH
        IF(ASSOCIATED(MESH)) THEN
          MESH_TOPOLOGY=>MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR
          IF(ASSOCIATED(MESH_TOPOLOGY)) THEN
            MESH_NODES=>MESH_TOPOLOGY%nodes
            IF(ASSOCIATED(MESH_NODES)) THEN
              NULLIFY(TREE_NODE)
              CALL TREE_SEARCH(MESH_NODES%nodesTree,USER_NODE_NUMBER,TREE_NODE,ERR,ERROR,*999)
              IF(ASSOCIATED(TREE_NODE)) THEN
                CALL TREE_NODE_VALUE_GET(MESH_NODES%nodesTree,TREE_NODE,GLOBAL_NODE_NUMBER,ERR,ERROR,*999)
                IF(GLOBAL_NODE_NUMBER>0.AND.GLOBAL_NODE_NUMBER<=MESH_TOPOLOGY%NODES%numberOfNodes) THEN
                  IF(MESH_COMPONENT_NUMBER>0.AND.MESH_COMPONENT_NUMBER<=MESH%NUMBER_OF_COMPONENTS) THEN
                    MESH_DOMAIN=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR
                    IF(ASSOCIATED(MESH_DOMAIN)) THEN
                      DOMAIN_NUMBER=MESH_DOMAIN%NODE_DOMAIN(GLOBAL_NODE_NUMBER)
                    ELSE
                      CALL FlagError("Decomposition domain is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="Mesh Component number "//TRIM(NumberToVString(MESH_COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid. The limits are 1 to "// &
                      & TRIM(NumberToVString(MESH%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))//"."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="Global element number found "//TRIM(NumberToVString(GLOBAL_NODE_NUMBER,"*",ERR,ERROR))// &
                    & " is invalid. The limits are 1 to "// &
                    & TRIM(NumberToVString(MESH_TOPOLOGY%NODES%numberOfNodes,"*",ERR,ERROR))//"."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="Decomposition mesh node corresponding to user node number "// &
                  & TRIM(NumberToVString(USER_NODE_NUMBER,"*",err,error))//" is not found."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Decomposition mesh nodes are not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Decomposition mesh topology is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Decomposition mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Decomposition has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("DECOMPOSITION_NODE_DOMAIN_GET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_NODE_DOMAIN_GET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_NODE_DOMAIN_GET

  !
  !================================================================================================================================
  !

  !>Initialises the embedded meshes.
  SUBROUTINE EMBEDDED_MESH_INITIALISE(MESH_EMBEDDING,ERR,ERROR,*)

    !Argument variables
    !TYPE(MESH_EMBEDDING_TYPE), INTENT(INOUT) :: MESH_EMBEDDING !<Mesh embedding to initialise
    TYPE(MESH_EMBEDDING_TYPE), POINTER :: MESH_EMBEDDING !<Mesh embedding to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("EMBEDDED_MESH_INITIALISE",ERR,ERROR,*998)
    
    ALLOCATE(MESH_EMBEDDING,STAT=ERR)
    NULLIFY(MESH_EMBEDDING%PARENT_MESH)
    NULLIFY(MESH_EMBEDDING%CHILD_MESH)
    
    EXITS("EMBEDDED_MESH_INITIALISE")
    RETURN
!999 CALL EMBEDDED_MESH_FINALISE(MESH_EMBEDDING,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("EMBEDDED_MESH_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE EMBEDDED_MESH_INITIALISE

  !
  !================================================================================================================================
  !

  !>Creates an embedding of one mesh in another
  SUBROUTINE MESH_EMBEDDING_CREATE(MESH_EMBEDDING, PARENT_MESH, CHILD_MESH,ERR,ERROR,*)
!    TYPE(MESH_EMBEDDING_TYPE), INTENT(INOUT) :: MESH_EMBEDDING !<Mesh embedding to create.
    TYPE(MESH_EMBEDDING_TYPE), POINTER :: MESH_EMBEDDING !<Mesh embedding to create.
    TYPE(MESH_TYPE), POINTER, INTENT(IN) :: PARENT_MESH !<The parent mesh
    TYPE(MESH_TYPE), POINTER, INTENT(IN) :: CHILD_MESH  !<The child mesh
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: NGP = 0, ne

    ENTERS("MESH_EMBEDDING_CREATE",ERR,ERROR,*999) 
    
    WRITE(*,*) 'parent mesh', PARENT_MESH%NUMBER_OF_ELEMENTS
    WRITE(*,*) 'child mesh', child_MESH%NUMBER_OF_ELEMENTS
    CALL EMBEDDED_MESH_INITIALISE(MESH_EMBEDDING,ERR,ERROR,*999)

    DO ne=1,PARENT_MESH%NUMBER_OF_ELEMENTS
      NGP = MAX(NGP,PARENT_MESH%TOPOLOGY(1)%PTR%ELEMENTS%ELEMENTS(ne)%BASIS%QUADRATURE%&
        & QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR%NUMBER_OF_GAUSS)
    ENDDO !ne

    MESH_EMBEDDING%PARENT_MESH => PARENT_MESH
    MESH_EMBEDDING%CHILD_MESH  => CHILD_MESH
    ALLOCATE(MESH_EMBEDDING%CHILD_NODE_XI_POSITION(PARENT_MESH%NUMBER_OF_ELEMENTS),STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate child node positions.",ERR,ERROR,*999)
    ALLOCATE(MESH_EMBEDDING%GAUSS_POINT_XI_POSITION(NGP,PARENT_MESH%NUMBER_OF_ELEMENTS),STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate gauss point positions.",ERR,ERROR,*999)
    
    EXITS("MESH_EMBEDDING_CREATE")
    RETURN 

999 ERRORSEXITS("MESH_EMBEDDING_CREATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_EMBEDDING_CREATE

  !
  !================================================================================================================================
  !

  !>Sets the positions of nodes in the child mesh for one element in the parent mesh
  SUBROUTINE MESH_EMBEDDING_SET_CHILD_NODE_POSITION(MESH_EMBEDDING, ELEMENT_NUMBER, NODE_NUMBERS, XI_COORDS,ERR,ERROR,*)
    TYPE(MESH_EMBEDDING_TYPE), INTENT(INOUT) :: MESH_EMBEDDING !<The mesh embedding object
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER  !<Element number in the parent mesh
    INTEGER(INTG), INTENT(IN) :: NODE_NUMBERS(:) !<NODE_NUMBERS(node_idx) Node numbers in child mesh for the node_idx'th embedded node in the ELEMENT_NUMBER'th element of the parent mesh
    REAL(DP), INTENT(IN)      :: XI_COORDS(:,:)  !<XI_COORDS(:,node_idx) Xi coordinates of the node_idx'th embedded node in the ELEMENT_NUMBER'th

    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    ENTERS("MESH_EMBEDDING_SET_CHILD_NODE_POSITION",ERR,ERROR,*999)

    IF(ELEMENT_NUMBER<1 .OR. ELEMENT_NUMBER > MESH_EMBEDDING%PARENT_MESH%NUMBER_OF_ELEMENTS) THEN
      CALL FlagError("Element number out of range",ERR,ERROR,*999)
    ENDIF

    MESH_EMBEDDING%CHILD_NODE_XI_POSITION(ELEMENT_NUMBER)%NUMBER_OF_NODES = SIZE(NODE_NUMBERS)

    ALLOCATE(MESH_EMBEDDING%CHILD_NODE_XI_POSITION(ELEMENT_NUMBER)%NODE_NUMBERS(SIZE(NODE_NUMBERS)))
    MESH_EMBEDDING%CHILD_NODE_XI_POSITION(ELEMENT_NUMBER)%NODE_NUMBERS(1:SIZE(NODE_NUMBERS)) = NODE_NUMBERS(1:SIZE(NODE_NUMBERS))

    ALLOCATE(MESH_EMBEDDING%CHILD_NODE_XI_POSITION(ELEMENT_NUMBER)%XI_COORDS(SIZE(XI_COORDS,1),SIZE(XI_COORDS,2)))
    MESH_EMBEDDING%CHILD_NODE_XI_POSITION(ELEMENT_NUMBER)%XI_COORDS(1:SIZE(XI_COORDS,1),1:SIZE(XI_COORDS,2)) = &
      & XI_COORDS(1:SIZE(XI_COORDS,1),1:SIZE(XI_COORDS,2))
    
    EXITS("MESH_EMBEDDING_SET_CHILD_NODE_POSITION")
    RETURN
999 ERRORSEXITS("MESH_EMBEDDING_SET_CHILD_NODE_POSITION",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_EMBEDDING_SET_CHILD_NODE_POSITION

  !
  !================================================================================================================================
  !

  !>Sets the positions of a Gauss point of the parent mesh in terms of element/xi coordinate in the child mesh
  SUBROUTINE MESH_EMBEDDING_SET_GAUSS_POINT_DATA(MESH_EMBEDDING, PARENT_ELEMENT_NUMBER, GAUSSPT_NUMBER,&
    & PARENT_XI_COORD, CHILD_ELEMENT_NUMBER, CHILD_XI_COORD,ERR,ERROR,*)
    TYPE(MESH_EMBEDDING_TYPE), INTENT(INOUT) :: MESH_EMBEDDING   !<The mesh embedding object
    INTEGER(INTG), INTENT(IN) :: PARENT_ELEMENT_NUMBER           !<Element number in the parent mesh
    INTEGER(INTG), INTENT(IN) :: GAUSSPT_NUMBER                  !<Gauss point number in this element
    REAL(DP), INTENT(IN) :: PARENT_XI_COORD(:)              !<Xi coordinate in parent element

    INTEGER(INTG), INTENT(IN) :: CHILD_ELEMENT_NUMBER !<Element number in the child mesh
    REAL(DP), INTENT(IN) :: CHILD_XI_COORD(:)    !<Xi coordinate in child element

    INTEGER(INTG), INTENT(OUT) :: ERR           !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR  !<The error string

    ENTERS("MESH_EMBEDDING_SET_GAUSS_POINT_DATA",ERR,ERROR,*999)

    IF(PARENT_ELEMENT_NUMBER<1 .OR. PARENT_ELEMENT_NUMBER > MESH_EMBEDDING%PARENT_MESH%NUMBER_OF_ELEMENTS) THEN
      CALL FlagError("Parent element number out of range",ERR,ERROR,*999)
    ENDIF
    IF(CHILD_ELEMENT_NUMBER<1 .OR. CHILD_ELEMENT_NUMBER > MESH_EMBEDDING%CHILD_MESH%NUMBER_OF_ELEMENTS) THEN
      CALL FlagError("Child element number out of range",ERR,ERROR,*999)
    ENDIF
    IF(GAUSSPT_NUMBER<1 .OR. GAUSSPT_NUMBER > SIZE(MESH_EMBEDDING%GAUSS_POINT_XI_POSITION,1)) THEN
      CALL FlagError("Gauss point number out of range",ERR,ERROR,*999)
    ENDIF

    ALLOCATE(MESH_EMBEDDING%GAUSS_POINT_XI_POSITION(GAUSSPT_NUMBER,PARENT_ELEMENT_NUMBER)&
     & %PARENT_XI_COORD(SIZE(PARENT_XI_COORD)))
    ALLOCATE(MESH_EMBEDDING%GAUSS_POINT_XI_POSITION(GAUSSPT_NUMBER,PARENT_ELEMENT_NUMBER)&
     & %CHILD_XI_COORD(SIZE(CHILD_XI_COORD)))


    MESH_EMBEDDING%GAUSS_POINT_XI_POSITION(GAUSSPT_NUMBER,PARENT_ELEMENT_NUMBER)%PARENT_XI_COORD(1:SIZE(PARENT_XI_COORD)) = &
      & PARENT_XI_COORD(1:SIZE(PARENT_XI_COORD))
    MESH_EMBEDDING%GAUSS_POINT_XI_POSITION(GAUSSPT_NUMBER,PARENT_ELEMENT_NUMBER)%CHILD_XI_COORD(1:SIZE(CHILD_XI_COORD)) = &
      & CHILD_XI_COORD(1:SIZE(CHILD_XI_COORD))
    MESH_EMBEDDING%GAUSS_POINT_XI_POSITION(GAUSSPT_NUMBER,PARENT_ELEMENT_NUMBER)%ELEMENT_NUMBER = CHILD_ELEMENT_NUMBER

    EXITS("MESH_EMBEDDING_SET_GAUSS_POINT_DATA")
    RETURN
999 ERRORSEXITS("MESH_EMBEDDING_SET_GAUSS_POINT_DATA",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_EMBEDDING_SET_GAUSS_POINT_DATA


  !
  !================================================================================================================================
  !

END MODULE MESH_ROUTINES

