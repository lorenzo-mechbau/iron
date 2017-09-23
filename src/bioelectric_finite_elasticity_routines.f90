!> \file
!> \authors Thomas Heidlauf
!> \brief This module handles all routines pertaining to bioelectrics coupled with finite elasticity.
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
!> Contributor(s): Thomas Klotz
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

!>This module handles all routines pertaining to bioelectrics coupled with finite elasticity.


MODULE BIOELECTRIC_FINITE_ELASTICITY_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE BIOELECTRIC_ROUTINES
  USE BIODOMAIN_EQUATION_ROUTINES
  USE COMP_ENVIRONMENT
  USE CONSTANTS
  USE CONTROL_LOOP_ROUTINES
  USE EQUATIONS_ROUTINES
  USE EQUATIONS_MAPPING_ROUTINES
  USE EQUATIONS_MATRICES_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_IO_ROUTINES
  USE FIELD_ROUTINES
  USE FINITE_ELASTICITY_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE LISTS
  USE MATHS
  USE MESH_ROUTINES
#ifndef NOMPIMOD
  USE MPI
#endif
  USE PRINT_TYPES_ROUTINES
  USE PROBLEM_CONSTANTS
  USE STRINGS
  USE SOLVER_ROUTINES
  USE TIMER
  USE TYPES
#include "macros.h"

  IMPLICIT NONE

  PUBLIC BioelectricFiniteElasticity_EquationsSetSetup
  PUBLIC BioelectricFiniteElasticity_EquationsSetSolutionMethodSet

  PUBLIC BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SETUP
  PUBLIC BioelectricFiniteElasticity_ProblemSpecificationSet
 
  PUBLIC BioelectricFiniteElasticity_FiniteElementCalculate

  PUBLIC BIOELECTRIC_FINITE_ELASTICITY_PRE_SOLVE
  PUBLIC BIOELECTRIC_FINITE_ELASTICITY_POST_SOLVE

  PUBLIC BioelectricFiniteElasticity_ControlLoopPreLoop
  PUBLIC BioelectricFiniteElasticity_ControlLoopPostLoop
  
  PUBLIC BioelectricFiniteElasticity_UpdateGeometricField

  REAL(DP), PUBLIC :: TIMING_FILE_OUTPUT_USER = 0_DP
  REAL(DP), PUBLIC :: TIMING_FILE_OUTPUT_SYSTEM = 0_DP
  
CONTAINS
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a bioelectrics finite elasticity equation type of a multi physics equations set class.
  SUBROUTINE BioelectricFiniteElasticity_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("BioelectricFiniteElasticity_EquationsSetSolutionMethodSet",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a "// &
          & "Bioelectric-finite elasticity type equations set.",err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE, &
        & EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE, &
        & EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
        & EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE, &
        & EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE)
        SELECT CASE(SOLUTION_METHOD)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
        CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The specified solution method of "//TRIM(NUMBER_TO_VSTRING(SOLUTION_METHOD,"*",ERR,ERROR))//" is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " is not valid for a bioelectrics finite elasticity equation type of a multi physics equations set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("BioelectricFiniteElasticity_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_EquationsSetSolutionMethodSet",ERR,ERROR)
    EXITS("BioelectricFiniteElasticity_EquationsSetSolutionMethodSet")
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Sets up the bioelectrics finite elasticity equation.
  SUBROUTINE BioelectricFiniteElasticity_EquationsSetSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    ENTERS("BioelectricFiniteElasticity_EquationsSetSetup",ERR,ERROR,*999)

    CALL FlagError("BioelectricFiniteElasticity_EquationsSetSetup is not implemented.",ERR,ERROR,*999)

    EXITS("BioelectricFiniteElasticity_EquationsSetSetup")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_EquationsSetSetup",ERR,ERROR)
    EXITS("BioelectricFiniteElasticity_EquationsSetSetup")
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a bioelectrics finite elasticity equation finite element equations set.
  SUBROUTINE BioelectricFiniteElasticity_FiniteElementCalculate(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    ENTERS("BioelectricFiniteElasticity_FiniteElementCalculate",ERR,ERROR,*999)

    CALL FlagError("BioelectricFiniteElasticity_FiniteElementCalculate is not implemented.",ERR,ERROR,*999)

    EXITS("BioelectricFiniteElasticity_FiniteElementCalculate")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_FiniteElementCalculate",ERR,ERROR)
    EXITS("BioelectricFiniteElasticity_FiniteElementCalculate")
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a bioelectric finite elasticity problem type .
  SUBROUTINE BioelectricFiniteElasticity_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype

    ENTERS("BioelectricFiniteElasticity_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(SIZE(problemSpecification,1)==3) THEN
        problemSubtype=problemSpecification(3)
        SELECT CASE(problemSubtype)
        CASE(PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE, &
            & PROBLEM_STRANG_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE, &
            & PROBLEM_GUDUNOV_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
            & PROBLEM_STRANG_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
            & PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE, &
            & PROBLEM_STRANG_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE)
          
        CASE DEFAULT
          localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
            & " is not valid for a bioelectric finite elasticity type of a multi physics problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        IF(ALLOCATED(problem%specification)) THEN
          CALL FlagError("Problem specification is already allocated.",err,error,*999)
        ELSE
          ALLOCATE(problem%specification(3),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
        END IF
        problem%specification(1:3)=[PROBLEM_MULTI_PHYSICS_CLASS,PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE, &
          & problemSubtype]
      ELSE
        CALL FlagError("Bioelectric finite elasticity problem specification must have three entries.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    END IF

    EXITS("BioelectricFiniteElasticity_ProblemSpecificationSet")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_ProblemSpecificationSet",err,error)
    EXITS("BioelectricFiniteElasticity_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the bioelectric finite elasticity problem.
  SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to setup
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(CONTROL_LOOP_TYPE), POINTER :: MONODOMAIN_SUB_LOOP,ELASTICITY_SUB_LOOP
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(CELLML_EQUATIONS_TYPE), POINTER :: CELLML_EQUATIONS
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS,MONODOMAIN_SOLVERS,ELASTICITY_SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SETUP",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(MONODOMAIN_SUB_LOOP)
    NULLIFY(ELASTICITY_SUB_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVERS)
    NULLIFY(MONODOMAIN_SOLVERS)
    NULLIFY(ELASTICITY_SOLVERS)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(CELLML_EQUATIONS)
    
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(problem%specification)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(problem%specification,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a bioelectric-finite elasticity problem.",err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(3))

      !--------------------------------------------------------------------
      !   Transient monodomain, simple finite elasticity
      !--------------------------------------------------------------------
     CASE(PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE, &
        & PROBLEM_STRANG_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE, &
        & PROBLEM_GUDUNOV_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
        & PROBLEM_STRANG_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
        & PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE, &
        & PROBLEM_STRANG_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE)
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a bioelectrics finite elasticity equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a time control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,ERR,ERROR,*999)
            CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(CONTROL_LOOP,2,ERR,ERROR,*999)
            CALL CONTROL_LOOP_OUTPUT_TYPE_SET(CONTROL_LOOP,CONTROL_LOOP_PROGRESS_OUTPUT,ERR,ERROR,*999)

            !Set up the control sub loop for monodomain
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,MONODOMAIN_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_LABEL_SET(MONODOMAIN_SUB_LOOP,'MONODOMAIN_TIME_LOOP',ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(MONODOMAIN_SUB_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,ERR,ERROR,*999)
            CALL CONTROL_LOOP_OUTPUT_TYPE_SET(MONODOMAIN_SUB_LOOP,CONTROL_LOOP_PROGRESS_OUTPUT,ERR,ERROR,*999)

            !Set up the control sub loop for finite elasicity
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,ELASTICITY_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_LABEL_SET(ELASTICITY_SUB_LOOP,'ELASTICITY_LOAD_INCREMENT_LOOP',ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(ELASTICITY_SUB_LOOP,PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE,ERR,ERROR,*999)
            CALL CONTROL_LOOP_OUTPUT_TYPE_SET(ELASTICITY_SUB_LOOP,CONTROL_LOOP_PROGRESS_OUTPUT,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,ERR,ERROR,*999)
            !Sub-loops are finished when parent is finished
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a bioelectrics finite elasticity equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the monodomain sub loop
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,MONODOMAIN_SUB_LOOP,ERR,ERROR,*999)
            !Start the solvers creation
            CALL SOLVERS_CREATE_START(MONODOMAIN_SUB_LOOP,MONODOMAIN_SOLVERS,ERR,ERROR,*999)
            
            SELECT CASE(PROBLEM%specification(3))
            CASE(PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE, &
              & PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE, &
              & PROBLEM_GUDUNOV_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE)
            ! ------------------------------------
            ! Godunov Splitting   
              CALL SOLVERS_NUMBER_SET(MONODOMAIN_SOLVERS,2,ERR,ERROR,*999)
              !Set the first solver to be a differential-algebraic equations solver
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(MONODOMAIN_SOLVERS,1,SOLVER,ERR,ERROR,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DAE_TYPE,ERR,ERROR,*999)
              CALL SOLVER_LABEL_SET(SOLVER,"ODE Solver",ERR,ERROR,*999)
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
              !Set the second solver to be a dynamic solver 
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(MONODOMAIN_SOLVERS,2,SOLVER,ERR,ERROR,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
              CALL SOLVER_LABEL_SET(SOLVER,"Parabolic solver",ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_RESTART_SET(SOLVER,.TRUE.,ERR,ERROR,*999)
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
            
            CASE(PROBLEM_STRANG_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE, &
              & PROBLEM_STRANG_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE, &
              & PROBLEM_STRANG_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE)
            
              ! ------------------------------------
              ! Strang Splitting
              CALL SOLVERS_NUMBER_SET(MONODOMAIN_SOLVERS,3,ERR,ERROR,*999)
              
              !Set the first solver to be a differential-algebraic equations solver
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(MONODOMAIN_SOLVERS,1,SOLVER,ERR,ERROR,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DAE_TYPE,ERR,ERROR,*999)
              CALL SOLVER_LABEL_SET(SOLVER,"ODE Solver",ERR,ERROR,*999)
              CALL SOLVER_DAE_EULER_SOLVER_TYPE_SET(SOLVER,SOLVER_DAE_EULER_IMPROVED,ERR,ERROR,*999)
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)

              !Set the second solver to be a dynamic solver 
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(MONODOMAIN_SOLVERS,2,SOLVER,ERR,ERROR,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
              CALL SOLVER_LABEL_SET(SOLVER,"Parabolic solver",ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_RESTART_SET(SOLVER,.TRUE.,ERR,ERROR,*999)
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
              
              !Set the third solver to be a differential-algebraic equations solver
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(MONODOMAIN_SOLVERS,3,SOLVER,ERR,ERROR,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DAE_TYPE,ERR,ERROR,*999)
              CALL SOLVER_LABEL_SET(SOLVER,"ODE Solver2",ERR,ERROR,*999)
              CALL SOLVER_DAE_EULER_SOLVER_TYPE_SET(SOLVER,SOLVER_DAE_EULER_IMPROVED,ERR,ERROR,*999)
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
              
            CASE DEFAULT
              LOCAL_ERROR="The problem specification(3) of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%specification(3),"*",ERR,ERROR))// &
                & " is invalid for a bioelectrics finite elasticity equation."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
            ! ------------------------------------
            
            !Get the finite elasticity sub loop
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,ELASTICITY_SUB_LOOP,ERR,ERROR,*999)
            !Start the solvers creation
            CALL SOLVERS_CREATE_START(ELASTICITY_SUB_LOOP,ELASTICITY_SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_NUMBER_SET(ELASTICITY_SOLVERS,1,ERR,ERROR,*999)
            !Set the finite elasticity solver to be a nonlinear solver
            NULLIFY(SOLVER)
            CALL SOLVERS_SOLVER_GET(ELASTICITY_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the monodomain solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,MONODOMAIN_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(MONODOMAIN_SUB_LOOP,MONODOMAIN_SOLVERS,ERR,ERROR,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(MONODOMAIN_SOLVERS,ERR,ERROR,*999)

            !Get the finite elasticity solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,ELASTICITY_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(ELASTICITY_SUB_LOOP,ELASTICITY_SOLVERS,ERR,ERROR,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(ELASTICITY_SOLVERS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a bioelectrics finite elasticity equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop and solvers
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)

            !Get the monodomain sub loop and solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,MONODOMAIN_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(MONODOMAIN_SUB_LOOP,MONODOMAIN_SOLVERS,ERR,ERROR,*999)
            !Create the solver equations for the second (parabolic) solver
            NULLIFY(SOLVER)
            CALL SOLVERS_SOLVER_GET(MONODOMAIN_SOLVERS,2,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)

            !Get the finite elasticity sub loop and solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,ELASTICITY_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(ELASTICITY_SUB_LOOP,ELASTICITY_SOLVERS,ERR,ERROR,*999)
            !Get the finite elasticity solver and create the finite elasticity solver equations
            NULLIFY(SOLVER)
            NULLIFY(SOLVER_EQUATIONS)
            CALL SOLVERS_SOLVER_GET(ELASTICITY_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            
            !Get the monodomain sub loop and solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,MONODOMAIN_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(MONODOMAIN_SUB_LOOP,MONODOMAIN_SOLVERS,ERR,ERROR,*999)
            !Get the solver equations for the second (parabolic) solver
            NULLIFY(SOLVER)
            NULLIFY(SOLVER_EQUATIONS)
            CALL SOLVERS_SOLVER_GET(MONODOMAIN_SOLVERS,2,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)             
            
            !Get the finite elasticity sub loop and solvers
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,ELASTICITY_SUB_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(ELASTICITY_SUB_LOOP,ELASTICITY_SOLVERS,ERR,ERROR,*999)
            !Finish the creation of the finite elasticity solver equations
            NULLIFY(SOLVER)
            NULLIFY(SOLVER_EQUATIONS)
            CALL SOLVERS_SOLVER_GET(ELASTICITY_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)             
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a bioelectrics finite elasticity equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CELLML_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,MONODOMAIN_SUB_LOOP,ERR,ERROR,*999)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(MONODOMAIN_SUB_LOOP,MONODOMAIN_SOLVERS,ERR,ERROR,*999)
            !Create the CellML equations for the first DAE solver
            CALL SOLVERS_SOLVER_GET(MONODOMAIN_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL CELLML_EQUATIONS_CREATE_START(SOLVER,CELLML_EQUATIONS,ERR,ERROR,*999)
            
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,MONODOMAIN_SUB_LOOP,ERR,ERROR,*999)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(MONODOMAIN_SUB_LOOP,MONODOMAIN_SOLVERS,ERR,ERROR,*999)
            !Get the CellML equations for the first DAE solver
            CALL SOLVERS_SOLVER_GET(MONODOMAIN_SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_CELLML_EQUATIONS_GET(SOLVER,CELLML_EQUATIONS,ERR,ERROR,*999)
            !Finish the CellML equations creation
            CALL CELLML_EQUATIONS_CREATE_FINISH(CELLML_EQUATIONS,ERR,ERROR,*999)
            
            
            SELECT CASE(PROBLEM%specification(3))
            CASE(PROBLEM_STRANG_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE, &
              & PROBLEM_STRANG_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE, &
              & PROBLEM_STRANG_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE)
                
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(MONODOMAIN_SOLVERS,3,SOLVER,ERR,ERROR,*999)
              
              ! set the cellml equations to the second DAE_solver
              SOLVER%CELLML_EQUATIONS=>CELLML_EQUATIONS
              
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a bioelectrics finite elasticity equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a bioelectrics finite elasticity equation."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " does not equal a transient monodomain quasistatic finite elasticity equation subtype."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SETUP")
    RETURN
999 ERRORSEXITS("BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SETUP",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SETUP

  !
  !================================================================================================================================
  !
 
  !>Sets up the bioelectrics finite elasticity problem pre-solve.
  SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("BIOELECTRIC_FINITE_ELASTICITY_PRE_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(control_loop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a bioelectric-finite elasticity problem.", &
              & err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
          CASE(PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE, &
            & PROBLEM_STRANG_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE, &
            & PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE, &
            & PROBLEM_STRANG_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE, &
            & PROBLEM_GUDUNOV_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
            & PROBLEM_STRANG_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE)
            SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
            CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
              CALL BIODOMAIN_PRE_SOLVE(SOLVER,ERR,ERROR,*999)
            CASE(PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE)
              CALL FINITE_ELASTICITY_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="Control loop loop type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%LOOP_TYPE,"*",ERR,ERROR))// &
                & " is not valid for a bioelectrics finite elasticity type of a multi physics problem class."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
              & " is not valid for a bioelectrics finite elasticity type of a multi physics problem class."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("BIOELECTRIC_FINITE_ELASTICITY_PRE_SOLVE")
    RETURN
999 ERRORSEXITS("BIOELECTRIC_FINITE_ELASTICITY_PRE_SOLVE",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_PRE_SOLVE
      
  !   
  !================================================================================================================================
  !

  !>Sets up the bioelectrics finite elasticity problem post solve.
  SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("BIOELECTRIC_FINITE_ELASTICITY_POST_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN 
          IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(control_loop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a bioelectric-finite elasticity problem.", &
              & err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
          CASE(PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE, &
            & PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE, &
            & PROBLEM_GUDUNOV_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
            & PROBLEM_STRANG_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE, &
            & PROBLEM_STRANG_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE, &
            & PROBLEM_STRANG_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE)
            SELECT CASE(SOLVER%SOLVE_TYPE)
            CASE(SOLVER_DAE_TYPE)
              CALL BIOELECTRIC_POST_SOLVE(SOLVER,ERR,ERROR,*999)
            CASE(SOLVER_DYNAMIC_TYPE)
              CALL BIOELECTRIC_POST_SOLVE(SOLVER,ERR,ERROR,*999)
            CASE(SOLVER_NONLINEAR_TYPE)
              CALL FINITE_ELASTICITY_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="Solver solve type "//TRIM(NUMBER_TO_VSTRING(SOLVER%SOLVE_TYPE,"*",ERR,ERROR))// &
                & " is not valid for a bioelectrics finite elasticity type of a multi physics problem class."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
              & " is not valid for a bioelectrics finite elasticity type of a multi physics problem class."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("BIOELECTRIC_FINITE_ELASTICITY_POST_SOLVE")
    RETURN
999 ERRORSEXITS("BIOELECTRIC_FINITE_ELASTICITY_POST_SOLVE",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_POST_SOLVE

  !
  !================================================================================================================================
  !

  !>Sets up the bioelectrics finite elasticity problem pre-control loop.
  SUBROUTINE BioelectricFiniteElasticity_ControlLoopPreLoop(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("BioelectricFiniteElasticity_ControlLoopPreLoop",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
        PROBLEM=>CONTROL_LOOP%PROBLEM
        IF(ASSOCIATED(PROBLEM)) THEN
          IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(control_loop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a bioelectric-finite elasticity problem.", &
              & err,error,*999)
          END IF
          SELECT CASE(PROBLEM%SPECIFICATION(2))
          CASE(PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE)
            SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
            CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
              !do nothing ???
            CASE(PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE)
              SELECT CASE(PROBLEM%SPECIFICATION(3))
              CASE(PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE, &
                & PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE, &
                & PROBLEM_STRANG_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE, &
                & PROBLEM_STRANG_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE)
                CALL BioelectricFiniteElasticity_IndependentFieldInterpolate(CONTROL_LOOP,ERR,ERROR,*999)
                
              CASE(PROBLEM_GUDUNOV_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
                & PROBLEM_STRANG_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE)
                CALL BIOELECTRIC_FINITE_ELASTICITY_COMPUTE_TITIN(CONTROL_LOOP,ERR,ERROR,*999)
                CALL BioelectricFiniteElasticity_IndependentFieldInterpolate(CONTROL_LOOP,ERR,ERROR,*999)
                
              CASE DEFAULT
                LOCAL_ERROR="The third problem specification of "// &
                  & TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
                  & " is not valid for bioelectric finite elasticity problem."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
              CALL FiniteElasticity_ControlTimeLoopPreLoop(CONTROL_LOOP,ERR,ERROR,*999)
            CASE(PROBLEM_CONTROL_WHILE_LOOP_TYPE)
              CALL FiniteElasticity_ControlTimeLoopPreLoop(CONTROL_LOOP,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="Control loop loop type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%LOOP_TYPE,"*",ERR,ERROR))// &
                & " is not valid for bioelectric finite elasticity problem type."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The second problem specification of "// &
              & TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
              & " is not valid for a multi physics problem."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Control loop problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        !the main time loop - do nothing!
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("BioelectricFiniteElasticity_ControlLoopPreLoop")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_ControlLoopPreLoop",ERR,ERROR)
    EXITS("BioelectricFiniteElasticity_ControlLoopPreLoop")
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_ControlLoopPreLoop

  !
  !================================================================================================================================
  !

  !>Computes the fibre stretch for bioelectrics finite elasticity problems.
  SUBROUTINE BioelectricFiniteElasticity_ComputeFibreStretch(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,FIBRE_FIELD,GEOMETRIC_FIELD,INDEPENDENT_FIELD
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: DEPENDENT_QUADRATURE_SCHEME
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: GEOMETRIC_INTERPOLATION_PARAMETERS, &
      & FIBRE_INTERPOLATION_PARAMETERS,DEPENDENT_INTERPOLATION_PARAMETERS
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: GEOMETRIC_INTERPOLATED_POINT,FIBRE_INTERPOLATED_POINT, &
      & DEPENDENT_INTERPOLATED_POINT
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: GEOMETRIC_INTERPOLATED_POINT_METRICS, &
      & DEPENDENT_INTERPOLATED_POINT_METRICS
    TYPE(BASIS_TYPE), POINTER :: GEOMETRIC_BASIS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VAR_U1
    INTEGER(INTG) :: DEPENDENT_NUMBER_OF_GAUSS_POINTS
    INTEGER(INTG) :: MESH_COMPONENT_NUMBER,NUMBER_OF_ELEMENTS,FIELD_VAR_TYPE
    INTEGER(INTG) :: equations_set_idx,gauss_idx,DofIdx,FEElementGlobalNumber,idx
    REAL(DP) :: DZDNU(3,3),DZDNUT(3,3),dZdXi(3,3),AZL(3,3)

    ENTERS("BioelectricFiniteElasticity_ComputeFibreStretch",ERR,ERROR,*999)

    NULLIFY(SOLVERS)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(FIELD_VAR_U1)
    NULLIFY(DEPENDENT_BASIS)
    NULLIFY(EQUATIONS)
    NULLIFY(DEPENDENT_FIELD,FIBRE_FIELD,GEOMETRIC_FIELD,INDEPENDENT_FIELD)
    NULLIFY(DEPENDENT_QUADRATURE_SCHEME)
    NULLIFY(GEOMETRIC_INTERPOLATION_PARAMETERS,FIBRE_INTERPOLATION_PARAMETERS,DEPENDENT_INTERPOLATION_PARAMETERS)
    NULLIFY(GEOMETRIC_INTERPOLATED_POINT,FIBRE_INTERPOLATED_POINT,DEPENDENT_INTERPOLATED_POINT)
    NULLIFY(GEOMETRIC_INTERPOLATED_POINT_METRICS,DEPENDENT_INTERPOLATED_POINT_METRICS)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
        SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
        CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
          !do nothing
        CASE(PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE)
          !do nothing
        CASE(PROBLEM_CONTROL_WHILE_LOOP_TYPE)
          CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
          CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
          CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
          !Loop over the equations sets associated with the solver
          IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
            SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
            IF(ASSOCIATED(SOLVER_MAPPING)) THEN
              DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                IF(ASSOCIATED(EQUATIONS_SET)) THEN
                  GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                  IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD)) THEN
                    LOCAL_ERROR="Geometric field is not associated for equations set index "// &
                      & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))//" in the solver mapping."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                  DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                  IF(.NOT.ASSOCIATED(DEPENDENT_FIELD)) THEN
                    LOCAL_ERROR="Dependent field is not associated for equations set index "// &
                      & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))//" in the solver mapping."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                  INDEPENDENT_FIELD=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
                  IF(.NOT.ASSOCIATED(INDEPENDENT_FIELD)) THEN
                    LOCAL_ERROR="Independent field is not associated for equations set index "// &
                      & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))//" in the solver mapping."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                  EQUATIONS=>EQUATIONS_SET%EQUATIONS
                  IF(.NOT.ASSOCIATED(EQUATIONS)) THEN
                    LOCAL_ERROR="Equations is not associated for equations set index "// &
                      & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))//" in the solver mapping."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                  FIBRE_FIELD=>EQUATIONS%INTERPOLATION%FIBRE_FIELD
                  IF(.NOT.ASSOCIATED(FIBRE_FIELD)) THEN
                    LOCAL_ERROR="Fibre field is not associated for equations set index "// &
                      & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))//" in the solver mapping."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="Equations set is not associated for equations set index "// &
                    & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))//" in the solver mapping."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDDO !equations_set_idx
            ELSE
              CALL FlagError("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Solver solver equations are not associated.",ERR,ERROR,*999)
          ENDIF

          CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_VAR_U1,ERR,ERROR,*999)

          DECOMPOSITION=>DEPENDENT_FIELD%DECOMPOSITION
          MESH_COMPONENT_NUMBER=DECOMPOSITION%MESH_COMPONENT_NUMBER

          IF(CONTROL_LOOP%WHILE_LOOP%ITERATION_NUMBER==1) THEN
            !copy the old fibre stretch to the previous values parameter set
            CALL FIELD_PARAMETER_SETS_COPY(INDEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,1.0_DP,ERR,ERROR,*999)
          ENDIF

          NUMBER_OF_ELEMENTS=DEPENDENT_FIELD%GEOMETRIC_FIELD%DECOMPOSITION%TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS

          !loop over the elements of the finite elasticity mesh (internal and boundary elements)
          DO FEElementGlobalNumber=1,NUMBER_OF_ELEMENTS

            DEPENDENT_BASIS=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(FEElementGlobalNumber)%BASIS       
            DEPENDENT_QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
            DEPENDENT_NUMBER_OF_GAUSS_POINTS=DEPENDENT_QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
              & TOPOLOGY%ELEMENTS%ELEMENTS(FEElementGlobalNumber)%BASIS

            !Initialise tensors and matrices
            DZDNU=0.0_DP
            DO idx=1,3
              DZDNU(idx,idx)=1.0_DP
            ENDDO

            !Grab interpolation parameters
            FIELD_VAR_TYPE=EQUATIONS_SET%EQUATIONS%EQUATIONS_MAPPING%NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR%VARIABLE_TYPE
            DEPENDENT_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR
            GEOMETRIC_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR
            FIBRE_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%FIBRE_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR

            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,FEElementGlobalNumber, &
              & GEOMETRIC_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,FEElementGlobalNumber, &
              & FIBRE_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,FEElementGlobalNumber, &
              & DEPENDENT_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)

            !Point interpolation pointer
            GEOMETRIC_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR
            GEOMETRIC_INTERPOLATED_POINT_METRICS=>EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR
            FIBRE_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%FIBRE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR
            DEPENDENT_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR
            DEPENDENT_INTERPOLATED_POINT_METRICS=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT_METRICS(FIELD_VAR_TYPE)%PTR

            !Loop over gauss points
            DO gauss_idx=1,DEPENDENT_NUMBER_OF_GAUSS_POINTS
            
              !Interpolate dependent, geometric, fibre and materials fields
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                & DEPENDENT_INTERPOLATED_POINT,ERR,ERROR,*999)
              CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(DEPENDENT_BASIS%NUMBER_OF_XI, &
                & DEPENDENT_INTERPOLATED_POINT_METRICS,ERR,ERROR,*999)
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                & GEOMETRIC_INTERPOLATED_POINT,ERR,ERROR,*999)
              CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI, &
                & GEOMETRIC_INTERPOLATED_POINT_METRICS,ERR,ERROR,*999)
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                & FIBRE_INTERPOLATED_POINT,ERR,ERROR,*999)

              !Calculate F=dZ/dNU, the deformation gradient tensor at the gauss point
              CALL FiniteElasticity_GaussDeformationGradientTensor(DEPENDENT_INTERPOLATED_POINT_METRICS, &
                & GEOMETRIC_INTERPOLATED_POINT_METRICS,FIBRE_INTERPOLATED_POINT,dZdXi,DZDNU,ERR,ERROR,*999)
              
              !compute C=F^T F
              CALL MATRIX_TRANSPOSE(DZDNU,DZDNUT,ERR,ERROR,*999)
              CALL MATRIX_PRODUCT(DZDNUT,DZDNU,AZL,ERR,ERROR,*999)

              !store the fibre stretch lambda_f, i.e., sqrt(C_11) or AZL(1,1)
              DofIdx=FIELD_VAR_U1%COMPONENTS(1)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx, &
               & FEElementGlobalNumber)
              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & DofIdx,SQRT(AZL(1,1)),ERR,ERROR,*999)

            ENDDO !gauss_idx
          ENDDO !FEElementGlobalNumber

          !now the ghost elements -- get the relevant info from the other computational nodes
          CALL FIELD_PARAMETER_SET_UPDATE_START(INDEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
          CALL FIELD_PARAMETER_SET_UPDATE_FINISH(INDEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)

        CASE DEFAULT
          LOCAL_ERROR="Control loop type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%LOOP_TYPE,"*",ERR,ERROR))// &
            & " is not valid for a bioelectrics finite elasticity type of a multi physics problem class."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE !the control loop contains subloops
        !do nothing
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("BioelectricFiniteElasticity_ComputeFibreStretch")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_ComputeFibreStretch",ERR,ERROR)
    EXITS("BioelectricFiniteElasticity_ComputeFibreStretch")
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_ComputeFibreStretch

  !
  !================================================================================================================================
  !

  !>Computes the bioelectrics finite elasticity force-length and force_velocity relations.
  SUBROUTINE BioelectricFiniteElasticity_ForceLengthVelocityRelation(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP_PARENT
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,FIBRE_FIELD,GEOMETRIC_FIELD,INDEPENDENT_FIELD
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: DEPENDENT_QUADRATURE_SCHEME
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: GEOMETRIC_INTERPOLATION_PARAMETERS, &
      & FIBRE_INTERPOLATION_PARAMETERS,DEPENDENT_INTERPOLATION_PARAMETERS
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: GEOMETRIC_INTERPOLATED_POINT,FIBRE_INTERPOLATED_POINT, &
      & DEPENDENT_INTERPOLATED_POINT
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: GEOMETRIC_INTERPOLATED_POINT_METRICS, &
      & DEPENDENT_INTERPOLATED_POINT_METRICS
    TYPE(BASIS_TYPE), POINTER :: GEOMETRIC_BASIS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VAR_U,FIELD_VAR_U1
    INTEGER(INTG) :: DEPENDENT_NUMBER_OF_GAUSS_POINTS
    INTEGER(INTG) :: MESH_COMPONENT_NUMBER,NUMBER_OF_ELEMENTS
    INTEGER(INTG) :: equations_set_idx,gauss_idx,DofIdx,FEElementGlobalNumber
    INTEGER(INTG) :: ITERATION_NUMBER,MAXIMUM_NUMBER_OF_ITERATIONS
    REAL(DP) :: LENGTH_HS,LENGTH_HS_0,ACTIVE_STRESS,FIBRE_STRETCH,FIBRE_STRETCH_OLD
    REAL(DP) :: FACTOR_LENGTH,FACTOR_VELO,SARCO_LENGTH,VELOCITY,VELOCITY_MAX,TIME_STEP,kappa,A,S,d,c

!tomo
    REAL(DP) :: VELOCITY_AVERAGE,STRETCH_AVERAGE,OLD_STRETCH_AVERAGE
    INTEGER(INTG) :: counter

    ENTERS("BioelectricFiniteElasticity_ForceLengthVelocityRelation",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP_PARENT)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(SOLVER_MAPPING)
    NULLIFY(GEOMETRIC_BASIS)
    NULLIFY(DECOMPOSITION)
    NULLIFY(SOLVERS)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(FIELD_VAR_U,FIELD_VAR_U1)
    NULLIFY(DEPENDENT_BASIS)
    NULLIFY(EQUATIONS)
    NULLIFY(DEPENDENT_FIELD,FIBRE_FIELD,GEOMETRIC_FIELD,INDEPENDENT_FIELD)
    NULLIFY(DEPENDENT_QUADRATURE_SCHEME)
    NULLIFY(GEOMETRIC_INTERPOLATION_PARAMETERS,FIBRE_INTERPOLATION_PARAMETERS,DEPENDENT_INTERPOLATION_PARAMETERS)
    NULLIFY(GEOMETRIC_INTERPOLATED_POINT,FIBRE_INTERPOLATED_POINT,DEPENDENT_INTERPOLATED_POINT)
    NULLIFY(GEOMETRIC_INTERPOLATED_POINT_METRICS,DEPENDENT_INTERPOLATED_POINT_METRICS)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
        SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
        CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
          !do nothing
        CASE(PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE)
          !do nothing
        CASE(PROBLEM_CONTROL_WHILE_LOOP_TYPE)
          CALL CONTROL_LOOP_GET(CONTROL_LOOP%PROBLEM%CONTROL_LOOP,CONTROL_LOOP_NODE,CONTROL_LOOP_PARENT,ERR,ERROR,*999)
          CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
          CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
          CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
          !Loop over the equations sets associated with the solver
          IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
            SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
            IF(ASSOCIATED(SOLVER_MAPPING)) THEN
              DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                IF(ASSOCIATED(EQUATIONS_SET)) THEN
                  DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                  IF(.NOT.ASSOCIATED(DEPENDENT_FIELD)) THEN
                    LOCAL_ERROR="Dependent field is not associated for equations set index "// &
                      & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))//" in the solver mapping."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                  INDEPENDENT_FIELD=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
                  IF(.NOT.ASSOCIATED(INDEPENDENT_FIELD)) THEN
                    LOCAL_ERROR="Independent field is not associated for equations set index "// &
                      & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))//" in the solver mapping."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="Equations set is not associated for equations set index "// &
                    & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))//" in the solver mapping."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDDO !equations_set_idx
            ELSE
              CALL FlagError("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Solver solver equations are not associated.",ERR,ERROR,*999)
          ENDIF

          CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VAR_U,ERR,ERROR,*999)
          CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_VAR_U1,ERR,ERROR,*999)

          !get the inital half-sarcomere length
          DofIdx=FIELD_VAR_U1%COMPONENTS(2)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,DofIdx,LENGTH_HS_0,ERR,ERROR,*999)

          !get the maximum contraction velocity
          DofIdx=FIELD_VAR_U1%COMPONENTS(3)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,DofIdx,VELOCITY_MAX,ERR,ERROR,*999)

          ITERATION_NUMBER=CONTROL_LOOP%WHILE_LOOP%ITERATION_NUMBER
          MAXIMUM_NUMBER_OF_ITERATIONS=CONTROL_LOOP%WHILE_LOOP%MAXIMUM_NUMBER_OF_ITERATIONS
          !in the first iteration store the unaltered homogenized active stress field 
          IF(ITERATION_NUMBER==1) THEN
            CALL FIELD_PARAMETER_SETS_COPY(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,1.0_DP,ERR,ERROR,*999)
          ELSE
            !restore the solution of the previous time step
            CALL FIELD_PARAMETER_SETS_COPY(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
              & FIELD_PREVIOUS_VALUES_SET_TYPE,FIELD_VALUES_SET_TYPE,1.0_DP,ERR,ERROR,*999)
          ENDIF

          TIME_STEP=CONTROL_LOOP_PARENT%TIME_LOOP%TIME_INCREMENT

          DECOMPOSITION=>DEPENDENT_FIELD%DECOMPOSITION
          MESH_COMPONENT_NUMBER=DECOMPOSITION%MESH_COMPONENT_NUMBER

!tomo start
          VELOCITY_AVERAGE=0.0_DP
          STRETCH_AVERAGE=0.0_DP
          OLD_STRETCH_AVERAGE=0.0_DP
          counter=0
!tomo end

          NUMBER_OF_ELEMENTS=DEPENDENT_FIELD%GEOMETRIC_FIELD%DECOMPOSITION%TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS

          !loop over the elements of the finite elasticity mesh (internal and boundary elements)
          DO FEElementGlobalNumber=1,NUMBER_OF_ELEMENTS

            DEPENDENT_BASIS=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(FEElementGlobalNumber)%BASIS       
            DEPENDENT_QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
            DEPENDENT_NUMBER_OF_GAUSS_POINTS=DEPENDENT_QUADRATURE_SCHEME%NUMBER_OF_GAUSS

            !Loop over gauss points
            DO gauss_idx=1,DEPENDENT_NUMBER_OF_GAUSS_POINTS

              !get the unaltered ACTIVE_STRESS at the GP
              DofIdx=FIELD_VAR_U%COMPONENTS(1)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx, &
               & FEElementGlobalNumber)
              CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE, &
                & DofIdx,ACTIVE_STRESS,ERR,ERROR,*999)

              ! FORCE-LENGTH RELATION -------------------------------------------------------------------------------
              
              !get the current fibre stretch at the GP
              DofIdx=FIELD_VAR_U1%COMPONENTS(1)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx, &
               & FEElementGlobalNumber)
              CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,DofIdx,FIBRE_STRETCH,ERR,ERROR,*999)

              !compute the current half-sarcomere length at the GP: l_hs = lambda_f * l_hs_0
              LENGTH_HS=LENGTH_HS_0*FIBRE_STRETCH
              
              !compute the scale factor (0,1) due to sarcomere F-l relation of Gordon, A. M., A.F. Huxley, and F.J. Julian. 
              !The variation in isometric tension with sarcomere length in vertebrate muscle fibres. 
              !The Journal of Physiology 184.1 (1966): 170-192.
              SARCO_LENGTH=2.0_DP*LENGTH_HS
              IF(SARCO_LENGTH.LE.1.27_DP) THEN
                FACTOR_LENGTH=0.0_DP
              ELSEIF(SARCO_LENGTH.LE.1.7_DP) THEN
                FACTOR_LENGTH=1.6047_DP*SARCO_LENGTH-2.0379_DP
              ELSEIF(SARCO_LENGTH.LE.2.34_DP) THEN
                FACTOR_LENGTH=0.4844_DP*SARCO_LENGTH-0.1334_DP
              ELSEIF(SARCO_LENGTH.LE.2.51_DP) THEN
                FACTOR_LENGTH=1.0_DP
              ELSEIF(SARCO_LENGTH.LE.3.94_DP) THEN
                FACTOR_LENGTH=-0.6993_DP*SARCO_LENGTH+2.7552_DP
              ELSE
                FACTOR_LENGTH=0.0_DP
              ENDIF

              !multiply the ACTIVE_STRESS with the scale factor
              ACTIVE_STRESS=ACTIVE_STRESS*FACTOR_LENGTH

              !update the ACTIVE_STRESS at GP
              DofIdx=FIELD_VAR_U%COMPONENTS(1)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx, &
               & FEElementGlobalNumber)
              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & DofIdx,ACTIVE_STRESS,ERR,ERROR,*999)





              ! FORCE-VELOCITY RELATION -------------------------------------------------------------------------------
              
              !get fibre stretch at the GP of the previous time step
              DofIdx=FIELD_VAR_U1%COMPONENTS(1)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx, &
               & FEElementGlobalNumber)
              CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                & FIELD_PREVIOUS_VALUES_SET_TYPE,DofIdx,FIBRE_STRETCH_OLD,ERR,ERROR,*999)

              !compute the contraction velocity
              VELOCITY=(FIBRE_STRETCH-FIBRE_STRETCH_OLD)/TIME_STEP

!!!! original
              !NOTE: VELOCITY_MAX is the max shortening velocity, and hence negative
              IF(VELOCITY<VELOCITY_MAX) THEN
                CALL FLAG_WARNING('Exceeded maximum contraction velocity (shortening).',ERR,ERROR,*999)
!                VELOCITY=VELOCITY_MAX
                !damping
                IF(ITERATION_NUMBER<(MAXIMUM_NUMBER_OF_ITERATIONS/2)) THEN
                  VELOCITY=VELOCITY*1.0_DP/DBLE((MAXIMUM_NUMBER_OF_ITERATIONS/2)-ITERATION_NUMBER)
                ENDIF
              ELSEIF(VELOCITY>(ABS(VELOCITY_MAX))) THEN
                ! warning disabled
                !CALL FLAG_WARNING('Exceeded maximum contraction velocity (lengthening).',ERR,ERROR,*999)
!!!                VELOCITY=-VELOCITY_MAX
!                !damping
!                IF(ITERATION_NUMBER<(MAXIMUM_NUMBER_OF_ITERATIONS/2)) THEN
!                  VELOCITY=VELOCITY*1.0_DP/DBLE((MAXIMUM_NUMBER_OF_ITERATIONS/2)-ITERATION_NUMBER)
!                ENDIF
              ENDIF

              !compute scale factor (FACTOR_VELO)
              kappa=0.24_DP !make mat param TODO
              A=1.6_DP !make mat param TODO
              S=2.5_DP !make mat param TODO
              IF(VELOCITY<0.0_DP) THEN
                !shortening contraction              
                FACTOR_VELO=(1-VELOCITY/VELOCITY_MAX)/(1+VELOCITY/VELOCITY_MAX/kappa)
              ELSE
                !lengthening contraction
                d=kappa*(1-A)
                c=VELOCITY/VELOCITY_MAX*S*(kappa+1)
                FACTOR_VELO=1+c*(A-1)/(d+c)
              ENDIF
              
              !multiply the ACTIVE_STRESS with the scale factor
              ACTIVE_STRESS=ACTIVE_STRESS*FACTOR_VELO

              !update the ACTIVE_STRESS at GP
              DofIdx=FIELD_VAR_U%COMPONENTS(1)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx, & 
               & FEElementGlobalNumber)
              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & DofIdx,ACTIVE_STRESS,ERR,ERROR,*999)

            ENDDO !gauss_idx
          ENDDO !FEElementGlobalNumber
!!!!end original




!!!!tomo --- try averaging              
!!!              VELOCITY_AVERAGE=VELOCITY_AVERAGE+VELOCITY
!!!              STRETCH_AVERAGE=STRETCH_AVERAGE+FIBRE_STRETCH
!!!              OLD_STRETCH_AVERAGE=OLD_STRETCH_AVERAGE+FIBRE_STRETCH_OLD
!!!              
!!!              counter=counter+1

!!!            ENDDO !gauss_idx
!!!          ENDDO !FEElementGlobalNumber

!!!          VELOCITY=VELOCITY_AVERAGE/REAL(counter)
!!!          FIBRE_STRETCH=STRETCH_AVERAGE/REAL(counter)
!!!          FIBRE_STRETCH_OLD=OLD_STRETCH_AVERAGE/REAL(counter)

!!!          VELOCITY_AVERAGE=(FIBRE_STRETCH-FIBRE_STRETCH_OLD)/TIME_STEP

!!!!          VELOCITY=0.0_DP

!!!!          !damping
!!!!          IF(ITERATION_NUMBER<(MAXIMUM_NUMBER_OF_ITERATIONS/2)) THEN
!!!!            VELOCITY=VELOCITY*1.0_DP/DBLE((MAXIMUM_NUMBER_OF_ITERATIONS/2)-ITERATION_NUMBER)
!!!!          ENDIF

!!!          LOCAL_ERROR="######### VELOCITY: "//TRIM(NUMBER_TO_VSTRING(VELOCITY,"*",ERR,ERROR))//" #########"
!!!          CALL FLAG_WARNING(LOCAL_ERROR,ERR,ERROR,*999)
!!!          LOCAL_ERROR="######### VELOCITY AVERAGE: "//TRIM(NUMBER_TO_VSTRING(VELOCITY_AVERAGE,"*",ERR,ERROR))//" #########"
!!!          CALL FLAG_WARNING(LOCAL_ERROR,ERR,ERROR,*999)
!!!          LOCAL_ERROR="######### STRETCH: "//TRIM(NUMBER_TO_VSTRING(FIBRE_STRETCH,"*",ERR,ERROR))//" #########"
!!!          CALL FLAG_WARNING(LOCAL_ERROR,ERR,ERROR,*999)
!!!          LOCAL_ERROR="######### OLD STRETCH: "//TRIM(NUMBER_TO_VSTRING(FIBRE_STRETCH_OLD,"*",ERR,ERROR))//" #########"
!!!          CALL FLAG_WARNING(LOCAL_ERROR,ERR,ERROR,*999)


!!!          
!!!          
!!!          VELOCITY=VELOCITY_AVERAGE
!!!          
!!!          
!!!          

!!!          !compute scale factor (FACTOR_VELO)
!!!          kappa=0.24_DP !make mat param TODO
!!!          A=1.6_DP !make mat param TODO
!!!          S=2.5_DP !make mat param TODO
!!!          IF(VELOCITY<0.0_DP) THEN
!!!            !shortening contraction              
!!!            FACTOR_VELO=(1-VELOCITY/VELOCITY_MAX)/(1+VELOCITY/VELOCITY_MAX/kappa)
!!!            IF(VELOCITY<VELOCITY_MAX) THEN
!!!              CALL FLAG_WARNING('Exceeded maximum contraction velocity (shortening).',ERR,ERROR,*999)
!!!            ENDIF
!!!          ELSE
!!!            !lengthening contraction
!!!            d=kappa*(1-A)
!!!            c=VELOCITY/VELOCITY_MAX*S*(kappa+1)
!!!            FACTOR_VELO=1+c*(A-1)/(d+c)
!!!            IF(VELOCITY>ABS(VELOCITY_MAX)) THEN
!!!              CALL FLAG_WARNING('Exceeded maximum contraction velocity (lengthening).',ERR,ERROR,*999)
!!!            ENDIF
!!!          ENDIF

!!!          DO FEElementGlobalNumber=1,NUMBER_OF_ELEMENTS

!!!            DEPENDENT_BASIS=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(FEElementGlobalNumber)%BASIS       
!!!            DEPENDENT_QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
!!!            DEPENDENT_NUMBER_OF_GAUSS_POINTS=DEPENDENT_QUADRATURE_SCHEME%NUMBER_OF_GAUSS

!!!            !Loop over gauss points
!!!            DO gauss_idx=1,DEPENDENT_NUMBER_OF_GAUSS_POINTS

!!!              !get the ACTIVE_STRESS at the GP
!!!              DofIdx=FIELD_VAR_U%COMPONENTS(1)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx,FEElementGlobalNumber)
!!!              CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
!!!                & DofIdx,ACTIVE_STRESS,ERR,ERROR,*999)

!!!              !multiply the ACTIVE_STRESS with the scale factor
!!!              ACTIVE_STRESS=ACTIVE_STRESS*FACTOR_VELO

!!!              !update the ACTIVE_STRESS at GP
!!!              DofIdx=FIELD_VAR_U%COMPONENTS(1)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx,FEElementGlobalNumber)
!!!              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
!!!                & DofIdx,ACTIVE_STRESS,ERR,ERROR,*999)
!!!              
!!!            ENDDO !gauss_idx
!!!          ENDDO !FEElementGlobalNumber
!tomo end

          !now the ghost elements -- get the relevant info from the other computational nodes
          CALL FIELD_PARAMETER_SET_UPDATE_START(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
          CALL FIELD_PARAMETER_SET_UPDATE_FINISH(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)

        CASE DEFAULT
          LOCAL_ERROR="Control loop type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%LOOP_TYPE,"*",ERR,ERROR))// &
            & " is not valid for a bioelectrics finite elasticity type of a multi physics problem class."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE !the control loop contains subloops
        !do nothing
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF


    EXITS("BioelectricFiniteElasticity_ForceLengthVelocityRelation")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_ForceLengthVelocityRelation",ERR,ERROR)
    EXITS("BioelectricFiniteElasticity_ForceLengthVelocityRelation")
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_ForceLengthVelocityRelation

  !
  !================================================================================================================================
  !

  !>Sets up the bioelectrics finite elasticity problem post-control loop.
  SUBROUTINE BioelectricFiniteElasticity_ControlLoopPostLoop(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG) :: equations_set_idx
    TYPE(CONTROL_LOOP_TIME_TYPE), POINTER :: TIME_LOOP
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(REGION_TYPE), POINTER :: DEPENDENT_REGION   
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: FILENAME,LOCAL_ERROR,METHOD
    TYPE(CONTROL_LOOP_TYPE), POINTER :: ELASTICITY_SUB_LOOP,BIOELECTRIC_SUB_LOOP
    REAL(SP) :: TIME_USER_START(1), TIME_USER_STOP(1), TIME_SYSTEM_START(1), TIME_SYSTEM_STOP(1)

    ENTERS("BioelectricFiniteElasticity_ControlLoopPostLoop",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
        PROBLEM=>CONTROL_LOOP%PROBLEM
        IF(ASSOCIATED(PROBLEM)) THEN
          IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(control_loop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a bioelectric-finite elasticity problem.", &
              & err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
          CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
            SELECT CASE(PROBLEM%SPECIFICATION(2))
            CASE(PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE)
              !the monodomain time loop - output of the monodomain fields
              CALL BIODOMAIN_CONTROL_LOOP_POST_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
                & " is not valid for a multi physics problem class."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE)
            CALL BioelectricFiniteElasticity_UpdateGeometricField(CONTROL_LOOP,.FALSE.,ERR,ERROR,*999)
          CASE(PROBLEM_CONTROL_WHILE_LOOP_TYPE)
            CALL BioelectricFiniteElasticity_ConvergenceCheck(CONTROL_LOOP,ERR,ERROR,*999)
          CASE DEFAULT
            !do nothing
          END SELECT
        ELSE
          CALL FlagError("Control loop problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        !the main time loop - output the finite elasticity fields MainTime_1_i and MainTime_M_2_<i>
        IF(CONTROL_LOOP%OUTPUT_TYPE>=CONTROL_LOOP_PROGRESS_OUTPUT .OR. CONTROL_LOOP%OUTPUT_TYPE == CONTROL_LOOP_FILE_OUTPUT) THEN
          TIME_LOOP=>CONTROL_LOOP%TIME_LOOP
          
          ! Only export in every nth time step, n=TIME_LOOP%OUTPUT_NUMBER, can be set via cmfe_ControlLoop_TimeOutputSet(ControlLoopMain,OUTPUT_FREQUENCY,Err)
          IF(MOD(TIME_LOOP%GLOBAL_ITERATION_NUMBER, TIME_LOOP%OUTPUT_NUMBER) == 0) THEN
            
            CALL CPU_TIMER(USER_CPU, TIME_USER_START, ERR,ERROR,*999)
            CALL CPU_TIMER(SYSTEM_CPU, TIME_SYSTEM_START, ERR,ERROR,*999)
            
            !Export the dependent field for this time step
            IF(ASSOCIATED(TIME_LOOP)) THEN
              PROBLEM=>CONTROL_LOOP%PROBLEM
              IF(ASSOCIATED(PROBLEM)) THEN
              
                NULLIFY(SOLVERS)
                NULLIFY(SOLVER)
                NULLIFY(SOLVER_EQUATIONS)
                NULLIFY(ELASTICITY_SUB_LOOP)
                !Get the solver. The first solver of the second sub loop will contain the finite elasticity dependent field equation set
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,ELASTICITY_SUB_LOOP,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(ELASTICITY_SUB_LOOP,SOLVERS,ERR,ERROR,*999)
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                !Loop over the equations sets associated with the solver
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                    DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                      EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                      IF(ASSOCIATED(EQUATIONS_SET)) THEN
                        DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                        NULLIFY(DEPENDENT_REGION)
                        CALL FIELD_REGION_GET(DEPENDENT_FIELD,DEPENDENT_REGION,ERR,ERROR,*999)
                        FILENAME="MainTime_"//TRIM(NUMBER_TO_VSTRING(DEPENDENT_REGION%USER_NUMBER,"*",ERR,ERROR))// &
                          & "_"//TRIM(NUMBER_TO_VSTRING(TIME_LOOP%GLOBAL_ITERATION_NUMBER,"*",ERR,ERROR))
                        METHOD="FORTRAN"
                        CALL FIELD_IO_NODES_EXPORT(DEPENDENT_REGION%FIELDS,FILENAME,METHOD,ERR,ERROR,*999)
                      ELSE
                        LOCAL_ERROR="Equations set is not associated for equations set index "// &
                          & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))// &
                          & " in the solver mapping."
                        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ENDDO !equations_set_idx
                  ELSE
                    CALL FlagError("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Solver solver equations are not associated.",ERR,ERROR,*999)
                ENDIF
                IF((PROBLEM%SPECIFICATION(3)==PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE).OR. &
                 & (PROBLEM%SPECIFICATION(3)==PROBLEM_GUDUNOV_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE).OR. &
                 & (PROBLEM%SPECIFICATION(3)==PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE).OR. &
                 & (PROBLEM%SPECIFICATION(3)==PROBLEM_STRANG_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE).OR. &
                 & (PROBLEM%SPECIFICATION(3)==PROBLEM_STRANG_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE).OR. &
                 & (PROBLEM%SPECIFICATION(3)==PROBLEM_STRANG_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE)) THEN
                  NULLIFY(SOLVERS)
                  NULLIFY(SOLVER)
                  NULLIFY(SOLVER_EQUATIONS)
                  NULLIFY(BIOELECTRIC_SUB_LOOP)
                  !Get the solver. The second solver of the first sub loop will contain the bioelectrics equation set
                  CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,BIOELECTRIC_SUB_LOOP,ERR,ERROR,*999)
                  CALL CONTROL_LOOP_SOLVERS_GET(BIOELECTRIC_SUB_LOOP,SOLVERS,ERR,ERROR,*999)
                  CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
                  CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                  !Loop over the equations sets associated with the solver
                  IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                    SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                    IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                      DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                        EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                        IF(ASSOCIATED(EQUATIONS_SET)) THEN
                          DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                          NULLIFY(DEPENDENT_REGION)
                          CALL FIELD_REGION_GET(DEPENDENT_FIELD,DEPENDENT_REGION,ERR,ERROR,*999)
                          FILENAME="MainTime_M_"//TRIM(NUMBER_TO_VSTRING(DEPENDENT_REGION%USER_NUMBER,"*",ERR,ERROR))// &
                            & "_"//TRIM(NUMBER_TO_VSTRING(TIME_LOOP%GLOBAL_ITERATION_NUMBER,"*",ERR,ERROR))
                          METHOD="FORTRAN"
                          CALL FIELD_IO_NODES_EXPORT(DEPENDENT_REGION%FIELDS,FILENAME,METHOD,ERR,ERROR,*999)
                          
                          !WRITE(*,*) TIME_LOOP%ITERATION_NUMBER
                          
                        ELSE
                          LOCAL_ERROR="Equations set is not associated for equations set index "// &
                            & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))// &
                            & " in the solver mapping."
                          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                        ENDIF
                      ENDDO !equations_set_idx
                    ELSE
                      CALL FlagError("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Solver solver equations are not associated.",ERR,ERROR,*999)
                  ENDIF
                ENDIF !PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE,PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE
              ELSE
                CALL FlagError("Control loop problem is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Time loop is not associated.",ERR,ERROR,*999)
            ENDIF
          
            CALL CPU_TIMER(USER_CPU, TIME_USER_STOP, ERR,ERROR,*999)
            CALL CPU_TIMER(SYSTEM_CPU, TIME_SYSTEM_STOP, ERR,ERROR,*999)
            TIMING_FILE_OUTPUT_USER = TIMING_FILE_OUTPUT_USER + (TIME_USER_STOP(1) - TIME_USER_START(1))
            TIMING_FILE_OUTPUT_SYSTEM = TIMING_FILE_OUTPUT_SYSTEM + (TIME_SYSTEM_STOP(1) - TIME_SYSTEM_START(1))
            !PRINT *, "duration file output: user: ", (TIME_USER_STOP(1) - TIME_USER_START(1))
            !PRINT *, "                      system: ",(TIME_SYSTEM_STOP(1) - TIME_SYSTEM_START(1)), &
            ! & ", new total duration: ",TIMING_FILE_OUTPUT_USER,",",TIMING_FILE_OUTPUT_SYSTEM
          
          ENDIF
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("BioelectricFiniteElasticity_ControlLoopPostLoop")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_ControlLoopPostLoop",ERR,ERROR)
    EXITS("BioelectricFiniteElasticity_ControlLoopPostLoop")
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_ControlLoopPostLoop


  !
  !================================================================================================================================
  !

  !>Check for the convergence of the bioelectric finite elasticity while loop, i.e.,
  !> if the force-length and force-velocity relations yielded a different actual configuration.
  !> This method is not called for cuboid.
  SUBROUTINE BioelectricFiniteElasticity_ConvergenceCheck(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    INTEGER(INTG) :: equations_set_idx,NUMBER_OF_NODES,DofIdx,node_idx
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VAR
    REAL(DP) :: x1,x2,x3,y1,y2,y3,my_sum

    ENTERS("BioelectricFiniteElasticity_ConvergenceCheck",ERR,ERROR,*999)

    NULLIFY(SOLVERS)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(FIELD_VAR)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
        PROBLEM=>CONTROL_LOOP%PROBLEM
        IF(ASSOCIATED(PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
          CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
            !do nothing
          CASE(PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE)
            !do nothing
          CASE(PROBLEM_CONTROL_WHILE_LOOP_TYPE)
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_SOLUTION_UPDATE(SOLVER,ERR,ERROR,*999) !tomo added this
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            !Loop over the equations sets associated with the solver
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
              IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                  EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                  IF(ASSOCIATED(EQUATIONS_SET)) THEN
                    DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                  ELSE
                    LOCAL_ERROR="Equations set is not associated for equations set index "// &
                      & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))// &
                      & " in the solver mapping."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ENDDO !equations_set_idx
              ELSE
                CALL FlagError("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Solver solver equations are not associated.",ERR,ERROR,*999)
            ENDIF

            IF(CONTROL_LOOP%WHILE_LOOP%ITERATION_NUMBER==1) THEN
              !
            ELSE
              CALL FIELD_VARIABLE_GET(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VAR,ERR,ERROR,*999)

              NUMBER_OF_NODES=DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION% &
                & MESH_COMPONENT_NUMBER)%PTR%MAPPINGS%NODES%NUMBER_OF_LOCAL

              my_sum=0.0_DP
              DO node_idx=1,NUMBER_OF_NODES

                !get the current node position (x) and the node position of the last iteration (y)
                DofIdx=FIELD_VAR%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)% &
                  & DERIVATIVES(1)%VERSIONS(1)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,DofIdx,x1,ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,DofIdx,y1,ERR,ERROR,*999)
                DofIdx=FIELD_VAR%COMPONENTS(2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)% &
                  & DERIVATIVES(1)%VERSIONS(1)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,DofIdx,x2,ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,DofIdx,y2,ERR,ERROR,*999)
                DofIdx=FIELD_VAR%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)% &
                  & DERIVATIVES(1)%VERSIONS(1)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,DofIdx,x3,ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,DofIdx,y3,ERR,ERROR,*999)

                !sum up
                my_sum=my_sum+SQRT((x1-y1)**2+(x2-y2)**2+(x3-y3)**2)
              ENDDO
              
              IF(my_sum<1.0E-06_DP) THEN !if converged then:
              
                PRINT*, "In bioelectric_finite_elasticity_routines.f90:1650: while loop converged after", &
                  CONTROL_LOOP%WHILE_LOOP%ITERATION_NUMBER, "iterations."
              
                CONTROL_LOOP%WHILE_LOOP%CONTINUE_LOOP=.FALSE.
                CALL BioelectricFiniteElasticity_UpdateGeometricField(CONTROL_LOOP,.FALSE.,ERR,ERROR,*999)
                !copy the current solution to the previous solution
                CALL FIELD_PARAMETER_SETS_COPY(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,1.0_DP,ERR,ERROR,*999)
              ELSEIF(CONTROL_LOOP%WHILE_LOOP%ITERATION_NUMBER==CONTROL_LOOP%WHILE_LOOP%MAXIMUM_NUMBER_OF_ITERATIONS) THEN
                CALL BioelectricFiniteElasticity_UpdateGeometricField(CONTROL_LOOP,.FALSE.,ERR,ERROR,*999)
                CALL FLAG_WARNING('----------- Maximum number of iterations in while loop reached. -----------',ERR,ERROR,*999)
                !copy the current solution to the previous solution
                CALL FIELD_PARAMETER_SETS_COPY(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,1.0_DP,ERR,ERROR,*999)
              ENDIF

            ENDIF

            !copy the current solution to the previous iteration solution
            CALL FIELD_PARAMETER_SETS_COPY(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,1.0_DP,ERR,ERROR,*999)

          CASE DEFAULT
            !do nothing
          END SELECT
        ELSE
          CALL FlagError("Control loop problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE !the main time loop
        !do nothing
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("BioelectricFiniteElasticity_ConvergenceCheck")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_ConvergenceCheck",ERR,ERROR)
    EXITS("BioelectricFiniteElasticity_ConvergenceCheck")
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_ConvergenceCheck

  !
  !================================================================================================================================
  !
  
  !> Get the local element number for an element with its global number
  SUBROUTINE BioelectricFiniteElasticity_GetLocalElementNumber(GEOMETRIC_FIELD, ElementGlobalNumber, ElementLocalNumber, &
   & ERR, ERROR, *)
    TYPE(FIELD_TYPE), POINTER, INTENT(IN) :: GEOMETRIC_FIELD  !< the geometric field of the elements
    INTEGER(INTG), INTENT(IN) :: ElementGlobalNumber !< the global element number of the element for which the local number is seeked
    INTEGER(INTG), INTENT(OUT) :: ElementLocalNumber !< the local number of the element with the global number (or 0 if the element is not on the local domain)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    ! local variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: DomainIdx, J
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: ELEMENTS_TOPOLOGY
    TYPE(DOMAIN_GLOBAL_MAPPING_TYPE), POINTER :: GLOBAL_TO_LOCAL_MAP
    LOGICAL, PARAMETER :: DEBUGGING = .FALSE.
    
    ENTERS("BioelectricFiniteElasticity_GetLocalElementNumber",ERR,ERROR,*999)

    ELEMENTS_TOPOLOGY=>GEOMETRIC_FIELD%DECOMPOSITION%TOPOLOGY%ELEMENTS
    
    ! get local element number from global element number
    GLOBAL_TO_LOCAL_MAP=>ELEMENTS_TOPOLOGY%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION% &
      & MESH_COMPONENT_NUMBER)%PTR%MAPPINGS%ELEMENTS%GLOBAL_TO_LOCAL_MAP(ElementGlobalNumber)
    DomainIdx = 0
    DO J = 1, GLOBAL_TO_LOCAL_MAP%NUMBER_OF_DOMAINS
      IF (GLOBAL_TO_LOCAL_MAP%DOMAIN_NUMBER(j) == COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)) THEN
        DomainIdx = J
        EXIT
      ENDIF
    ENDDO
    
    IF (DomainIdx == 0) THEN
      LOCAL_ERROR="Element with global number "//TRIM(NUMBER_TO_VSTRING(ElementGlobalNumber,"*",ERR,ERROR))// &
        & " is not on local domain of process "//TRIM(NUMBER_TO_VSTRING(COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR),"*",ERR,ERROR))
      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
    ENDIF

    ElementLocalNumber = GLOBAL_TO_LOCAL_MAP%LOCAL_NUMBER(DomainIdx)
    
    IF (DEBUGGING) &
      PRINT "(I1.1,2(A,I6),2(A,I2))", COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR),": global->local: ", &
       & ElementGlobalNumber,"->",ElementLocalNumber, & 
       & ", domain",DomainIdx,"/",&
       & GLOBAL_TO_LOCAL_MAP%NUMBER_OF_DOMAINS
    
    EXITS("BioelectricFiniteElasticity_GetLocalElementNumber")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_GetLocalElementNumber",ERR,ERROR)
    EXITS("BioelectricFiniteElasticity_GetLocalElementNumber")
    RETURN 1
  END SUBROUTINE BioelectricFiniteElasticity_GetLocalElementNumber

  !
  !================================================================================================================================
  !
  
  !> Advance to the next monodomain node and set the corresponding FE element
  SUBROUTINE BioelectricFiniteElasticity_IterateNextMonodomainNode(SOLVER_MAPPING_MONODOMAIN, GEOMETRIC_FIELD_ELASTICITY, &
   & INDEPENDENT_FIELD_ELASTICITY, &
   & Initialize, GEOMETRIC_FIELD_MONODOMAIN, DEPENDENT_FIELD_MONODOMAIN, INDEPENDENT_FIELD_MONODOMAIN, &
   & FIELD_VAR_DEP_M, FIELD_VAR_GEO_M, FIELD_VAR_IND_M_U1, FIELD_VAR_IND_M_U2, FIELD_VAR_IND_M_V, &
   & IsFinished, FEElementGlobalNumber, FEElementLocalNumber, BioelectricNodeGlobalNumber, BioelectricNodeLocalNumber, &
   & BioelectricNodeInFibreNumber, LeftNeighbourBioelectricNodeLocalNumber, FibreNumber, XI, DELTA_XI1, &
   & IsFirstBioelectricNodeOfFibre, LeftNeighbourIsGhost, PreviousNodeHadNoLeftNeighbour, ERR, ERROR, *)
    TYPE(SOLVER_MAPPING_TYPE), POINTER, INTENT(IN) :: SOLVER_MAPPING_MONODOMAIN
    TYPE(FIELD_TYPE), POINTER, INTENT(IN) :: GEOMETRIC_FIELD_ELASTICITY
    TYPE(FIELD_TYPE), POINTER, INTENT(IN) :: INDEPENDENT_FIELD_ELASTICITY
    LOGICAL, INTENT(IN) :: Initialize
    
    TYPE(FIELD_TYPE), POINTER, INTENT(OUT) :: DEPENDENT_FIELD_MONODOMAIN
    TYPE(FIELD_TYPE), POINTER, INTENT(OUT) :: INDEPENDENT_FIELD_MONODOMAIN
    TYPE(FIELD_TYPE), POINTER, INTENT(OUT) :: GEOMETRIC_FIELD_MONODOMAIN
    TYPE(FIELD_VARIABLE_TYPE), POINTER, INTENT(OUT) :: FIELD_VAR_DEP_M
    TYPE(FIELD_VARIABLE_TYPE), POINTER, INTENT(OUT) :: FIELD_VAR_GEO_M
    TYPE(FIELD_VARIABLE_TYPE), POINTER, INTENT(OUT) :: FIELD_VAR_IND_M_U1
    TYPE(FIELD_VARIABLE_TYPE), POINTER, INTENT(OUT) :: FIELD_VAR_IND_M_U2
    TYPE(FIELD_VARIABLE_TYPE), POINTER, INTENT(OUT) :: FIELD_VAR_IND_M_V
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    LOGICAL, INTENT(OUT) :: IsFinished    !< if the iteration over all monodomain nodes is finished
    INTEGER(INTG), INTENT(OUT) :: FEElementGlobalNumber  !< current FE element user number
    INTEGER(INTG), INTENT(OUT) :: FEElementLocalNumber     !< local element number of FE element
    INTEGER(INTG), INTENT(OUT) :: BioelectricNodeGlobalNumber   !< global number of the bioelectric node, per fibre, i.e. there can be multiple nodes with the same global number
    INTEGER(INTG), INTENT(OUT) :: BioelectricNodeLocalNumber    !< local number of the bioelectric node
    INTEGER(INTG), INTENT(OUT) :: BioelectricNodeInFibreNumber    !< contiguous number of the node starting with 1 for the first node of the fibre
    INTEGER(INTG), INTENT(OUT) :: LeftNeighbourBioelectricNodeLocalNumber    !< local number of the bioelectrics node to the left of the current, if it exists, else 0
    INTEGER(INTG), INTENT(OUT) :: FibreNumber  !< index of the current fibre
    REAL(DP), INTENT(OUT) :: XI(3)        !< position of the current bioelectric node in the current FE element, in [0,1]^3
    REAL(DP), INTENT(OUT) :: DELTA_XI1    !< the increment of XI(1) between two bioelectric nodes of a fibre
    LOGICAL, INTENT(OUT) :: IsFirstBioelectricNodeOfFibre    !< if the current node is the first node of the fibre
    LOGICAL, INTENT(OUT) :: LeftNeighbourIsGhost   !< If LeftNeighbourBioelectricNodeLocalNumber is the number of a ghost node to the local domain
    LOGICAL, INTENT(OUT) :: PreviousNodeHadNoLeftNeighbour   !< for the previous node LeftNeighbourBioelectricNodeLocalNumber was 0, i.e. the current node is next to a node that had no left neighbour (useful to update the values of the left node from the current)
    
    ! local variables
    INTEGER(INTG), SAVE :: PreviousFibreNumber
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET_MONODOMAIN
    TYPE(DOMAIN_TYPE), POINTER, SAVE :: MONODOMAIN_DOMAIN
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER, SAVE :: MONODOMAIN_MAPPINGS
    TYPE(DOMAIN_MAPPING_TYPE), POINTER, SAVE :: M_NODES_MAPPING
    TYPE(DOMAIN_MAPPING_TYPE), POINTER, SAVE :: M_ELEMENTS_MAPPING
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER, SAVE :: M_DOMAIN_TOPOLOGY
    TYPE(DOMAIN_NODES_TYPE), POINTER, SAVE :: M_DOMAIN_TOPOLOGY_NODES
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER, SAVE :: M_DOMAIN_TOPOLOGY_ELEMENTS
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER, SAVE :: M_ELEMENTS_TOPOLOGY
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER, SAVE :: FE_ELEMENTS_TOPOLOGY
    TYPE(FIELD_VARIABLE_TYPE), POINTER, SAVE  :: FIELD_VAR_IND_FE
    
    !INTEGER(INTG), SAVE :: NumberOfElementsMPerElementFE
    !INTEGER(INTG), SAVE :: BioelectricElementIdx
    INTEGER(INTG), SAVE :: BioelectricElementLocalNumber
    LOGICAL, SAVE :: CurrentBioelectricNodeIsLeftNode               !< if the currently considered monodomain node is the left node of the current monodomain element (only true at the beginning of a fibre domain, subsequently always the right node of the element is considered)
    INTEGER(INTG), SAVE :: ComputationalNodeNumber
    INTEGER(INTG), SAVE :: LastBioelectricElementLocalNumber
    INTEGER(INTG), SAVE :: LastFEElementGlobalNumber
    INTEGER(INTG), SAVE :: BioelectricElementInFEElementNumber
    INTEGER(INTG), SAVE :: BioelectricNodeInFEElementNumber
    INTEGER(INTG), SAVE :: NumberOfBioelectricElementsPerElementFE
    INTEGER(INTG), SAVE :: NumberOfNodesInXi2PerElementFE
    INTEGER(INTG), SAVE :: NumberOfNodesInXi3PerElementFE
    INTEGER(INTG), SAVE :: NumberGlobalYElements
    INTEGER(INTG), SAVE :: LastLeftNeighbourBioelectricNodeLocalNumber     !< the value of LastLeftNeighbourBioelectricNodeLocalNumber after the last call to this function
    INTEGER(INTG) :: TemporaryFEElementLocalNumber, TemporaryPreviousFEElementLocalNumber
    INTEGER(INTG) :: TemporaryBiolectricNodeLocalNumber
    INTEGER(INTG) :: DofIdx
    INTEGER(INTG) :: LeftNodeFibreNumber
    INTEGER(INTG) :: NodeDomain
    INTEGER(INTG) :: I    
    INTEGER(INTG) :: PreviousBioelectricElementLocalNumber
    LOGICAL :: DEBUGGING = .FALSE.
    !LOGICAL, PARAMETER :: DEBUGGING = .FALSE.
    
    ENTERS("BioelectricFiniteElasticity_IterateNextMonodomainNode",ERR,ERROR,*999)

#if 0     
    IF (COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR) == 0) THEN
      DEBUGGING = .TRUE.
    ENDIF
#endif    
    
    IsFinished = .FALSE.
    IsFirstBioelectricNodeOfFibre = .FALSE.
  
    ! in the first IF-ELSE statement the bioelectric element number gets set/advanced
    ! afterwards all needed numbers get derived
  
    ! Initialize values
    IF (Initialize) THEN
      
      IF (DEBUGGING) THEN
        PRINT *, "'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''"
        PRINT *, "BioelectricFiniteElasticity_IterateNextMonodomainNode, initialize"
      ENDIF
      
      ! field variables
      EQUATIONS_SET_MONODOMAIN=>SOLVER_MAPPING_MONODOMAIN%EQUATIONS_SETS(1)%PTR       ! Select the 1st equation set. If there were multiple fibre meshes, this would be affected.
      GEOMETRIC_FIELD_MONODOMAIN=>EQUATIONS_SET_MONODOMAIN%GEOMETRY%GEOMETRIC_FIELD
      DEPENDENT_FIELD_MONODOMAIN=>EQUATIONS_SET_MONODOMAIN%DEPENDENT%DEPENDENT_FIELD
      INDEPENDENT_FIELD_MONODOMAIN=>EQUATIONS_SET_MONODOMAIN%INDEPENDENT%INDEPENDENT_FIELD
      
      NULLIFY(FIELD_VAR_GEO_M)
      NULLIFY(FIELD_VAR_DEP_M)
      NULLIFY(FIELD_VAR_IND_M_U1)
      NULLIFY(FIELD_VAR_IND_M_U2)
      NULLIFY(FIELD_VAR_IND_M_V)
      NULLIFY(FIELD_VAR_IND_FE)
    
      CALL FIELD_VARIABLE_GET(GEOMETRIC_FIELD_MONODOMAIN,  FIELD_U_VARIABLE_TYPE, FIELD_VAR_GEO_M,   ERR,ERROR,*999)
      CALL FIELD_VARIABLE_GET(DEPENDENT_FIELD_MONODOMAIN,  FIELD_V_VARIABLE_TYPE, FIELD_VAR_DEP_M,   ERR,ERROR,*999)
      CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U1_VARIABLE_TYPE,FIELD_VAR_IND_M_U1,ERR,ERROR,*999)
      CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE,FIELD_VAR_IND_M_U2,ERR,ERROR,*999)
      CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, FIELD_VAR_IND_M_V, ERR,ERROR,*999)
      CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD_ELASTICITY,FIELD_V_VARIABLE_TYPE, FIELD_VAR_IND_FE,  ERR,ERROR,*999)

      ! helper variables
      MONODOMAIN_DOMAIN=>GEOMETRIC_FIELD_MONODOMAIN%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD_MONODOMAIN%DECOMPOSITION% &
       & MESH_COMPONENT_NUMBER)%PTR              ! DOMAIN_TYPE    
      MONODOMAIN_MAPPINGS=>MONODOMAIN_DOMAIN%MAPPINGS
      M_ELEMENTS_MAPPING=>MONODOMAIN_MAPPINGS%ELEMENTS
      M_NODES_MAPPING=>MONODOMAIN_MAPPINGS%NODES

      M_DOMAIN_TOPOLOGY=>MONODOMAIN_DOMAIN%TOPOLOGY     ! TYPE DOMAIN_TOPOLOGY_TYPE
      M_DOMAIN_TOPOLOGY_ELEMENTS=>M_DOMAIN_TOPOLOGY%ELEMENTS  ! TYPE DOMAIN_ELEMENTS_TYPE
              
      M_ELEMENTS_TOPOLOGY=>GEOMETRIC_FIELD_MONODOMAIN%DECOMPOSITION%TOPOLOGY%ELEMENTS
      FE_ELEMENTS_TOPOLOGY=>GEOMETRIC_FIELD_ELASTICITY%DECOMPOSITION%TOPOLOGY%ELEMENTS
              
      ComputationalNodeNumber = COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
            
      IF (SIZE(M_ELEMENTS_MAPPING%DOMAIN_LIST) >= M_ELEMENTS_MAPPING%GHOST_START) THEN
        LastBioelectricElementLocalNumber = M_ELEMENTS_MAPPING%DOMAIN_LIST(M_ELEMENTS_MAPPING%GHOST_START)
      ELSE
        LastBioelectricElementLocalNumber = M_ELEMENTS_MAPPING%DOMAIN_LIST(M_ELEMENTS_MAPPING%GHOST_START-1)+1
      ENDIF
      
      ! indices
      CurrentBioelectricNodeIsLeftNode = .TRUE.
      PreviousNodeHadNoLeftNeighbour = .FALSE.
      
      PreviousFibreNumber = 0
      BioelectricNodeInFEElementNumber = 1
      LastFEElementGlobalNumber = 0
      LastLeftNeighbourBioelectricNodeLocalNumber = -1 ! this may not be initialized to 0, because then PreviousNodeHadNoLeftNeighbour would be .TRUE. in the first call to this function
      
      !BioelectricElementIdx = M_ELEMENTS_MAPPING%INTERNAL_START
      ! directly initialize local element number to 1
      BioelectricElementLocalNumber = 1
      
      NumberGlobalYElements = GEOMETRIC_FIELD_ELASTICITY%DECOMPOSITION%MESH%GENERATED_MESH%REGULAR_MESH%NUMBER_OF_ELEMENTS_XI(2)
    
    ELSE
      IF (CurrentBioelectricNodeIsLeftNode) THEN
        ! if the last node was the left node of its element now use the right node of this element
        CurrentBioelectricNodeIsLeftNode = .FALSE.
        IF (DEBUGGING) PRINT *," at el. ", BioelectricElementLocalNumber,", going from left node to right node of element"
      ELSE
        ! advance to next element
        !BioelectricElementIdx = BioelectricElementIdx + 1
        BioelectricElementLocalNumber = BioelectricElementLocalNumber + 1
        IF (DEBUGGING) PRINT *," going to el. ", BioelectricElementLocalNumber
        
        BioelectricElementInFEElementNumber = BioelectricElementInFEElementNumber + 1
        
        IF (BioelectricElementLocalNumber >= LastBioelectricElementLocalNumber) THEN
          IF (DEBUGGING) PRINT "(2(A,I3))", " el ",BioelectricElementLocalNumber," has reached finish ", &
            & LastBioelectricElementLocalNumber
          IsFinished = .TRUE.
          RETURN
        ENDIF
      ENDIF
      
      ! also increment the counter of bioelectric node in the current FE element
      BioelectricNodeInFEElementNumber = BioelectricNodeInFEElementNumber + 1
      
    ENDIF ! Initialize            
    
    ! get the local number of bioelectric element from the mapping index
    !BioelectricElementLocalNumber = M_ELEMENTS_MAPPING%DOMAIN_LIST(BioelectricElementIdx)
    !IF (DEBUGGING) PRINT *, " el idx. ", BioelectricElementIdx," is el no ", BioelectricElementLocalNumber
    
    IF (DEBUGGING) THEN
      PRINT "(A,I3,A,I3,A)", " element ", BioelectricElementLocalNumber, " has", &
        & SIZE(M_DOMAIN_TOPOLOGY_ELEMENTS%ELEMENTS(BioelectricElementLocalNumber)%ELEMENT_NODES), &
        & " nodes: "
      DO I = 1,SIZE(M_DOMAIN_TOPOLOGY_ELEMENTS%ELEMENTS(BioelectricElementLocalNumber)%ELEMENT_NODES)
        PRINT "(A,I3,2(A,I3))", "I=",I, &
          & ": local ",M_DOMAIN_TOPOLOGY_ELEMENTS%ELEMENTS(BioelectricElementLocalNumber)%ELEMENT_NODES(I), &
          & ", global ",M_NODES_MAPPING%LOCAL_TO_GLOBAL_MAP(M_DOMAIN_TOPOLOGY_ELEMENTS%ELEMENTS(BioelectricElementLocalNumber)% &
            & ELEMENT_NODES(I))
      ENDDO
    ENDIF
    
    ! save the previous bioelectric node number as left neighbour node
    LeftNeighbourBioelectricNodeLocalNumber = BioelectricNodeLocalNumber
    LeftNeighbourIsGhost = .FALSE.
    
    DO I = 1,2
      ! get the current bioelectric node (local number) from the current local element (local number)
      IF (CurrentBioelectricNodeIsLeftNode) THEN
        BioelectricNodeLocalNumber = M_DOMAIN_TOPOLOGY_ELEMENTS%ELEMENTS(BioelectricElementLocalNumber)%ELEMENT_NODES(1)
        IF (DEBUGGING) PRINT *, " take left node, local",BioelectricNodeLocalNumber
      ELSE
        IF (SIZE(M_DOMAIN_TOPOLOGY_ELEMENTS%ELEMENTS(BioelectricElementLocalNumber)%ELEMENT_NODES) == 2) THEN
          BioelectricNodeLocalNumber = M_DOMAIN_TOPOLOGY_ELEMENTS%ELEMENTS(BioelectricElementLocalNumber)%ELEMENT_NODES(2)
          IF (DEBUGGING) PRINT *, " take right node, local",BioelectricNodeLocalNumber
        ELSE        
          ! if the element only contains one node, take this (should not happen)
          BioelectricNodeLocalNumber = M_DOMAIN_TOPOLOGY_ELEMENTS%ELEMENTS(BioelectricElementLocalNumber)%ELEMENT_NODES(1)
          IF (DEBUGGING) PRINT *, " take right and only node, local",BioelectricNodeLocalNumber
        ENDIF
      ENDIF
      
      ! get global bioelectric node number from local bioelectric node number
      BioelectricNodeGlobalNumber = M_NODES_MAPPING%LOCAL_TO_GLOBAL_MAP(BioelectricNodeLocalNumber)
      
      ! determine if this node is on the local domain (it may not, despite the element is)
      CALL DECOMPOSITION_NODE_DOMAIN_GET(GEOMETRIC_FIELD_MONODOMAIN%DECOMPOSITION,BioelectricNodeGlobalNumber,1,NodeDomain, &
        & ERR,ERROR,*999)
      IF (DEBUGGING) PRINT "(2(A,I3))", " Node global", BioelectricNodeGlobalNumber, " is on domain ", NodeDomain
      
      ! if node is not on local domain, iterate to next node
      IF (NodeDomain /= ComputationalNodeNumber) THEN
        IF (DEBUGGING) PRINT *," Node is not on local domain, go to next node"
        LeftNeighbourIsGhost = .TRUE.
        
        IF (CurrentBioelectricNodeIsLeftNode) THEN
          ! if the last node was the left node of its element, now use the right node of this element
          CurrentBioelectricNodeIsLeftNode = .FALSE.
          IF (DEBUGGING) PRINT *," at element local", BioelectricElementLocalNumber, &
            & ", going from left node to right node of element"
            
          LeftNeighbourBioelectricNodeLocalNumber = BioelectricNodeLocalNumber
        ELSE
          ! advance to next element
          !BioelectricElementIdx = BioelectricElementIdx + 1
          BioelectricElementLocalNumber = BioelectricElementLocalNumber + 1
          IF (DEBUGGING) PRINT *," going to element local", BioelectricElementLocalNumber
                  
          IF (BioelectricElementLocalNumber >= LastBioelectricElementLocalNumber) THEN
            IF (DEBUGGING) PRINT "(2(A,I3))", " element local",BioelectricElementLocalNumber," has reached finish ", &
              & LastBioelectricElementLocalNumber
            IsFinished = .TRUE.
            RETURN
          ENDIF
          
          CurrentBioelectricNodeIsLeftNode = .TRUE.
        ENDIF
        
        ! get the current bioelectric node (local number) from the current local element (local number)
        IF (SIZE(M_DOMAIN_TOPOLOGY_ELEMENTS%ELEMENTS(BioelectricElementLocalNumber)%ELEMENT_NODES) == 2) THEN
          BioelectricNodeLocalNumber = M_DOMAIN_TOPOLOGY_ELEMENTS%ELEMENTS(BioelectricElementLocalNumber)%ELEMENT_NODES(2)
          IF (DEBUGGING) PRINT *, " take right node local",BioelectricNodeLocalNumber
        ELSE        
          ! if the element only contains one node, take this (should not happen)
          BioelectricNodeLocalNumber = M_DOMAIN_TOPOLOGY_ELEMENTS%ELEMENTS(BioelectricElementLocalNumber)%ELEMENT_NODES(1)
          IF (DEBUGGING) PRINT *, " take right and only node local",BioelectricNodeLocalNumber
        ENDIF
        
        ! get global bioelectric node number from local bioelectric node number
        BioelectricNodeGlobalNumber = M_NODES_MAPPING%LOCAL_TO_GLOBAL_MAP(BioelectricNodeLocalNumber)

        ! if something in the current endif branch fails, try to use a recursive call like the following instead (tested and works)
        !CALL BioelectricFiniteElasticity_IterateNextMonodomainNode(SOLVER_MAPPING_MONODOMAIN, FE_ELEMENTS_TOPOLOGY, .FALSE., &
        ! & GEOMETRIC_FIELD_MONODOMAIN, DEPENDENT_FIELD_MONODOMAIN, INDEPENDENT_FIELD_MONODOMAIN, &
        ! & FIELD_VAR_DEP_M, FIELD_VAR_GEO_M, FIELD_VAR_IND_M_U1, FIELD_VAR_IND_M_U2, FIELD_VAR_IND_M_V, &
        ! & IsFinished, FEElementGlobalNumber, FEElementLocalNumber, BioelectricNodeGlobalNumber, BioelectricNodeLocalNumber, &
        ! & BioelectricNodeInFibreNumber, LeftNeighbourBioelectricNodeLocalNumber, FibreNumber, XI, IsFirstBioelectricNodeOfFibre, &
        ! & ERR, ERROR, *999)
        !RETURN
      ENDIF ! node not on local domain
      
      ! get fibre number
      DofIdx=FIELD_VAR_IND_M_V%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%&
        & NODES(BioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
        & FIELD_VALUES_SET_TYPE,DofIdx,FibreNumber,ERR,ERROR,*999)

      IF (DEBUGGING) PRINT "(2(A,I3))", " Node ", BioelectricNodeLocalNumber, " is of fibre ", FibreNumber
        
      ! determine if this is a new fibre
      IF (FibreNumber /= PreviousFibreNumber) THEN
        PreviousFibreNumber = FibreNumber
        
        !IF (PreviousFibreNumber /= 0) THEN
        IF (DEBUGGING) PRINT *, " new fibre"

        ! because this is the first element of the fibre in the local domain, update the node to be the left node of the element (not the right one as usual), this is done in the next iteration of the I loop
        CurrentBioelectricNodeIsLeftNode = .TRUE.
        
      ELSE    ! this is the same fibre as the last node, so the obtained right node of the current element is correct
        EXIT
      ENDIF
    ENDDO
    
    ! increment the counter of bioelectric nodes in the current FE element
    !!BioelectricNodeInFEElementNumber = BioelectricNodeInFEElementNumber + 1
        
    ! set the left neighbour node
    IF (CurrentBioelectricNodeIsLeftNode) THEN
      ! if there is no left neighbour element that would contain the neighbour node
      IF (M_ELEMENTS_TOPOLOGY%ELEMENTS(BioelectricElementLocalNumber)%ADJACENT_ELEMENTS(-1)%NUMBER_OF_ADJACENT_ELEMENTS == 0) THEN
        LeftNeighbourBioelectricNodeLocalNumber = 0
        LeftNeighbourIsGhost = .FALSE.
        IF (DEBUGGING) PRINT *, "no previous node"
      ELSE    ! there is a left neighbour element
        PreviousBioelectricElementLocalNumber = M_ELEMENTS_TOPOLOGY%ELEMENTS(BioelectricElementLocalNumber)% &
         & ADJACENT_ELEMENTS(-1)%ADJACENT_ELEMENTS(1)
        ! take left node of element
        LeftNeighbourBioelectricNodeLocalNumber = M_DOMAIN_TOPOLOGY_ELEMENTS%ELEMENTS(PreviousBioelectricElementLocalNumber)% &
          & ELEMENT_NODES(1)
        LeftNeighbourIsGhost = .TRUE.
        IF (DEBUGGING) &
          &  PRINT *, "previous node is of ElementM loc.", PreviousBioelectricElementLocalNumber, " node loc.", &
          & LeftNeighbourBioelectricNodeLocalNumber
      ENDIF
    ENDIF
        
    ! determine if this is the first node of a fibre or if the fibre already started more left on a different processors domain
    IsFirstBioelectricNodeOfFibre = .FALSE.
    IF (CurrentBioelectricNodeIsLeftNode) THEN    ! if this is another fibre than the previous element
      
      ! if the current node has no left neighbour
      IF (LeftNeighbourBioelectricNodeLocalNumber == 0) THEN  ! if there is no left neighbour, it is the first node of the 1D mesh and therefore the first node of a fibre
        IsFirstBioelectricNodeOfFibre = .TRUE.
      ELSE
        ! this does not work because there is no communication that sets fibre numbers in ghost nodes
        ! but it is not required when there are no in-series fibres, so no problem here
#if 0        
        ! get the fibre number of the left node, maybe it is from a different in series fibre
        DofIdx=FIELD_VAR_IND_M_V%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%&
          & NODES(LeftNeighbourBioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
        CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,DofIdx,LeftNodeFibreNumber,ERR,ERROR,*999)

        IF (DEBUGGING) PRINT *, "LeftNodeFibreNumber: ", LeftNodeFibreNumber,", current FibreNumber", FibreNumber
        IF (LeftNodeFibreNumber /= FibreNumber) THEN
          IsFirstBioelectricNodeOfFibre = .TRUE.
        ENDIF
#endif        
      ENDIF
    ENDIF
      
    ! get corresponding finite elasticity element (local element number)
    DofIdx=FIELD_VAR_IND_M_V%COMPONENTS(5)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%&
      & NODES(BioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
    CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
      & FIELD_VALUES_SET_TYPE,DofIdx,FEElementLocalNumber,ERR,ERROR,*999)
  
    IF (DEBUGGING) PRINT *, "bioelectric node local", BioelectricNodeLocalNumber, " is in FEElement: local", FEElementLocalNumber
  
    ! get global number of FE element
    LastFEElementGlobalNumber = FEElementGlobalNumber
    FEElementGlobalNumber = FE_ELEMENTS_TOPOLOGY%ELEMENTS(FEElementLocalNumber)%GLOBAL_NUMBER

    ! if new FE element is reached, reset BioelectricElementInFEElementNumber counter
    IF (CurrentBioelectricNodeIsLeftNode) THEN    ! new element is reached and left node is considered
      IF (IsFirstBioelectricNodeOfFibre) THEN     ! the fibre starts here, so the current node is on the left boundary of the domain and the first node of the element
        BioelectricNodeInFEElementNumber = 1
      ELSE      ! the fibre started ealier, so the current node is owned by the left FE element and the last node of that element
        BioelectricNodeInFEElementNumber = NumberOfBioelectricElementsPerElementFE+1
      ENDIF
    ELSEIF (LastFEElementGlobalNumber /= FEElementGlobalNumber) THEN ! new element is reached
      BioelectricNodeInFEElementNumber = 2
    ENDIF
    
    IF (DEBUGGING) PRINT *, "CurrentBioelectricNodeIsLeftNode: ", CurrentBioelectricNodeIsLeftNode, &
      & ", IsFirstBioelectricNodeOfFibre: ", IsFirstBioelectricNodeOfFibre, &
      & ", BioelectricNodeInFEElementNumber: ", BioelectricNodeInFEElementNumber
      
    ! set flag if last node had no previous neighbour
    PreviousNodeHadNoLeftNeighbour = (LastLeftNeighbourBioelectricNodeLocalNumber == 0)
    LastLeftNeighbourBioelectricNodeLocalNumber = LeftNeighbourBioelectricNodeLocalNumber
    
    IF (Initialize) THEN
      ! determine number of bioelectric elements per FE element
      DofIdx=FIELD_VAR_IND_FE%COMPONENTS(1)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%ELEMENTS(FEElementLocalNumber)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
        & DofIdx,NumberOfBioelectricElementsPerElementFE,ERR,ERROR,*999)
      
      DofIdx=FIELD_VAR_IND_FE%COMPONENTS(2)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%ELEMENTS(FEElementLocalNumber)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
        & DofIdx,NumberOfNodesInXi2PerElementFE,ERR,ERROR,*999)
      
      DofIdx=FIELD_VAR_IND_FE%COMPONENTS(3)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%ELEMENTS(FEElementLocalNumber)
      CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
        & DofIdx,NumberOfNodesInXi3PerElementFE,ERR,ERROR,*999)
        
      ! compute the increment for XI(1) between nodes
      DELTA_XI1 = 1.0 / NumberOfBioelectricElementsPerElementFE
      
      IF (DEBUGGING) THEN 
        PRINT *, "NumberOfBioelectricElementsPerElementFE: ", NumberOfBioelectricElementsPerElementFE
        PRINT *, "NumberOfNodesInXi2PerElementFE: ", NumberOfNodesInXi2PerElementFE
        PRINT *, "NumberOfNodesInXi3PerElementFE: ", NumberOfNodesInXi3PerElementFE
        PRINT *, "NumberGlobalYElements: ", NumberGlobalYElements
      ENDIF
    ENDIF
    
    ! compute value of XI
    ! XI is the position inside the containing FE element, in a coordinate frame [0,1]^3
    XI(1) = FLOAT(BioelectricNodeInFEElementNumber-1) / NumberOfBioelectricElementsPerElementFE     ! in direction of fibre
    XI(2) = (MOD(FibreNumber-1, NumberOfNodesInXi2PerElementFE)+0.5) / NumberOfNodesInXi2PerElementFE  ! in y-direction
    XI(3) = (MOD(INT((FibreNumber-1) / (NumberOfNodesInXi2PerElementFE * NumberGlobalYElements)), &    ! in z-direction
      & NumberOfNodesInXi3PerElementFE)+0.5) / NumberOfNodesInXi3PerElementFE
      
    IF (DEBUGGING) PRINT "(I1, A, I5, A, I5)", ComputationalNodeNumber, ": global ",BioelectricNodeGlobalNumber, &
       & ", BioelectricNodeInFEElementNumber=",BioelectricNodeInFEElementNumber
      
    EXITS("BioelectricFiniteElasticity_IterateNextMonodomainNode")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_IterateNextMonodomainNode",ERR,ERROR)
    EXITS("BioelectricFiniteElasticity_IterateNextMonodomainNode")
    RETURN 1
  END SUBROUTINE BioelectricFiniteElasticity_IterateNextMonodomainNode
  
  !
  !================================================================================================================================
  !
  !> Make all processes wait until one sets the variable gdb_resume to 1 via a debugger (e.g. gdb)
  SUBROUTINE gdbParallelDebuggingBarrier()
    INTEGER(Intg) :: Gdb_Resume
    INTEGER(Intg) :: MPI_IERROR, MPI_COMM_WORLD
    INTEGER(Intg) :: ComputationalNodeNumber, NumberOfComputationalNodes
    Gdb_Resume = 0

    CALL MPI_COMM_RANK(MPI_COMM_WORLD,ComputationalNodeNumber,MPI_IERROR)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NumberOfComputationalNodes,MPI_IERROR)

    IF (NumberOfComputationalNodes > 0) THEN
      PRINT*, "Node ", ComputationalNodeNumber, ", UID ",GETPID()," is waiting for Gdb_Resume=", Gdb_Resume &
        & , " to become 1 " // NEW_LINE('A') // "sudo gdb cuboid ",GETPID(), NEW_LINE('A') //"select-frame 2" // &
        & NEW_LINE('A') // "set var gdb_resume = 1" // NEW_LINE('A') // &
        & "info locals" // NEW_LINE('A') // "next"
      DO WHILE (Gdb_Resume == 0)
        CALL Sleep(1)
      ENDDO
      PRINT*, "Node ", ComputationalNodeNumber, " resumes because gdb_resume=", Gdb_Resume, "."
    ENDIF

  END SUBROUTINE gdbParallelDebuggingBarrier

  !
  !================================================================================================================================
  !
  !> Make a consistent state such that also in bioelectric node ghost elements the GeometricFieldM FIELD_U_VARIABLE_TYPE GeometryM 
  !> field has the correct values. This is done by exchanging node distances over adjacent domains.
  !> This method can be considered a helper method for BioelectricFiniteElasticity_UpdateGeometricField.
  SUBROUTINE BioelectricFiniteElasticity_GhostElements(INDEPENDENT_FIELD_MONODOMAIN,GEOMETRIC_FIELD_MONODOMAIN,&
    & M_NODES_MAPPING,M_ELEMENTS_MAPPING,M_DOMAIN_TOPOLOGY_ELEMENTS,FIELD_VAR_GEO_M,FIELD_VAR_IND_M_U2,M_ELEMENTS_TOPOLOGY,&
    & ERR,ERROR,*)
    TYPE(FIELD_TYPE), POINTER, INTENT(IN) :: INDEPENDENT_FIELD_MONODOMAIN
    TYPE(FIELD_TYPE), POINTER, INTENT(IN) :: GEOMETRIC_FIELD_MONODOMAIN
    TYPE(DOMAIN_MAPPING_TYPE), POINTER, INTENT(IN) :: M_NODES_MAPPING
    TYPE(DOMAIN_MAPPING_TYPE), POINTER, INTENT(IN) :: M_ELEMENTS_MAPPING
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER, INTENT(IN) :: M_DOMAIN_TOPOLOGY_ELEMENTS
    TYPE(FIELD_VARIABLE_TYPE), POINTER, INTENT(IN) :: FIELD_VAR_GEO_M
    TYPE(FIELD_VARIABLE_TYPE), POINTER, INTENT(IN) :: FIELD_VAR_IND_M_U2
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER, INTENT(IN) :: M_ELEMENTS_TOPOLOGY
    
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    INTEGER(INTG) :: BioelectricElementIdx, BioelectricElementLocalNumber, BioelectricElementGlobalNumber
    INTEGER(INTG) :: BioelectricNodeIdx, BioelectricNodeLocalNumber, BioelectricNodeGlobalNumber
    INTEGER(INTG) :: LeftBioelectricNodeLocalNumber, RightBioelectricNodeLocalNumber
    INTEGER(INTG) :: LeftBioelectricNodeGlobalNumber, RightBioelectricNodeGlobalNumber
    INTEGER(INTG) :: SecondLeftBioelectricNodeLocalNumber, SecondRightBioelectricNodeLocalNumber
    INTEGER(INTG) :: LeftBioelectricElementLocalNumber, RightBioelectricElementLocalNumber
    REAL(DP) :: Position1DLeftNode, Position1DRightNode, DistanceToLeftNode, DistanceToRightNode, DistanceNodes
    REAL(DP) :: Position1DSecondLeftNode, Position1DSecondRightNode
    REAL(DP) :: TemporaryValue
    INTEGER(INTG) :: DofIdx, I
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: DISTRIBUTED_VECTOR
    INTEGER(INTG) :: DOMAIN_IDX, RecvIndex, SendIndex
    !LOGICAL, PARAMETER :: DEBUGGING = .TRUE.   ! enable debugging output with this parameter
    LOGICAL :: DEBUGGING = .FALSE.   ! enable debugging output with this parameter

    INTEGER(INTG), SAVE :: CALL_NO = 0
    
    TYPE(LIST_TYPE), POINTER :: List
    
    ENTERS("BioelectricFiniteElasticity_GhostElements",ERR,ERROR,*999)

    CALL_NO = CALL_NO + 1
    
    ! test case for LIST_REMOVE_DUPLICATES_WITHOUT_SORTING
#if 0
    NULLIFY(List)
    CALL LIST_CREATE_START(List,ERR,ERROR,*999)
    CALL LIST_DATA_TYPE_SET(List,LIST_INTG_TYPE,ERR,ERROR,*999)
    CALL LIST_INITIAL_SIZE_SET(List,2,ERR,ERROR,*999)
    CALL LIST_CREATE_FINISH(List,ERR,ERROR,*999)
    
    CALL LIST_ITEM_ADD(List,5,ERR,ERROR,*999)
    CALL LIST_ITEM_ADD(List,5,ERR,ERROR,*999)
    CALL LIST_ITEM_ADD(List,5,ERR,ERROR,*999)
    CALL LIST_ITEM_ADD(List,5,ERR,ERROR,*999)
    CALL LIST_ITEM_ADD(List,5,ERR,ERROR,*999)
    !CALL LIST_ITEM_ADD(List,3,ERR,ERROR,*999)
    !CALL LIST_ITEM_ADD(List,7,ERR,ERROR,*999)
    !CALL LIST_ITEM_ADD(List,5,ERR,ERROR,*999)
    !CALL LIST_ITEM_ADD(List,7,ERR,ERROR,*999)
    !CALL LIST_ITEM_ADD(List,7,ERR,ERROR,*999)
    !CALL LIST_ITEM_ADD(List,8,ERR,ERROR,*999)
    !CALL LIST_ITEM_ADD(List,8,ERR,ERROR,*999)
    !CALL LIST_ITEM_ADD(List,1,ERR,ERROR,*999)
    !CALL LIST_ITEM_ADD(List,2,ERR,ERROR,*999)
    !CALL LIST_ITEM_ADD(List,3,ERR,ERROR,*999)
    
    PRINT *, "list"
    DO I=1,List%NUMBER_IN_LIST
      PRINT *, List%LIST_INTG(I)
    ENDDO
    
    CALL LIST_REMOVE_DUPLICATES_WITHOUT_SORTING(List, ERR,ERROR,*999)
    !CALL LIST_REMOVE_DUPLICATES(List, ERR,ERROR,*999)
    
    PRINT *, "list w/o duplicates"
    DO I=1,List%NUMBER_IN_LIST
      PRINT *, List%LIST_INTG(I)
    ENDDO
    
    PRINT *, "stop in bioelectric_finite_elasticity_routines.f90:2244"
    STOP
#endif
    
#if 0
    DEBUGGING = .FALSE.
    IF (COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR) == 0) DEBUGGING = .TRUE.
#endif    
    
    IF (DEBUGGING) THEN
      ! output all values
      PRINT *, "process ",COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR), " in BioelectricFiniteElasticity_GhostElements, " // &
        & "before communication"
      PRINT *, "    process, glob. node, dof no.(l),        distance (l),              distance (r)"
      PRINT *, "internal"
      DO BioelectricNodeIdx = M_NODES_MAPPING%INTERNAL_START, M_NODES_MAPPING%INTERNAL_FINISH
        BioelectricNodeLocalNumber = M_NODES_MAPPING%DOMAIN_LIST(BioelectricNodeIdx)
        BioelectricNodeGlobalNumber = M_NODES_MAPPING%LOCAL_TO_GLOBAL_MAP(BioelectricNodeLocalNumber)
        
        ! store global number in U(1)
        !CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
        !  & FIELD_VALUES_SET_TYPE,1,1,BioelectricNodeLocalNumber,1,DBLE(BioelectricNodeGlobalNumber),ERR,ERROR,*999)
        
        DofIdx=FIELD_VAR_IND_M_U2%COMPONENTS(7)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%&
          & NODES(BioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
        CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,DofIdx,DistanceToRightNode,ERR,ERROR,*999)
          
        DofIdx=FIELD_VAR_IND_M_U2%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%&
          & NODES(BioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
        CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,DofIdx,DistanceToLeftNode,ERR,ERROR,*999)
        
        PRINT *, COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR), BioelectricNodeGlobalNumber, DofIdx, &
          & DistanceToLeftNode, DistanceToRightNode
      ENDDO
      PRINT *, "boundary"
      DO BioelectricNodeIdx = M_NODES_MAPPING%BOUNDARY_START, M_NODES_MAPPING%BOUNDARY_FINISH
        BioelectricNodeLocalNumber = M_NODES_MAPPING%DOMAIN_LIST(BioelectricNodeIdx)
        BioelectricNodeGlobalNumber = M_NODES_MAPPING%LOCAL_TO_GLOBAL_MAP(BioelectricNodeLocalNumber)
        
        ! store global number in U(1)
        !CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
        !  & FIELD_VALUES_SET_TYPE,1,1,BioelectricNodeLocalNumber,1,DBLE(BioelectricNodeGlobalNumber),ERR,ERROR,*999)
          
        DofIdx=FIELD_VAR_IND_M_U2%COMPONENTS(7)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%&
          & NODES(BioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
        CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,DofIdx,DistanceToRightNode,ERR,ERROR,*999)
        
        DofIdx=FIELD_VAR_IND_M_U2%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%&
          & NODES(BioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
        CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,DofIdx,DistanceToLeftNode,ERR,ERROR,*999)
        
        PRINT *, COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR), BioelectricNodeGlobalNumber, DofIdx, &
          & DistanceToLeftNode, DistanceToRightNode
      ENDDO
      PRINT *, "ghost"
      DO BioelectricNodeIdx = M_NODES_MAPPING%GHOST_START, M_NODES_MAPPING%GHOST_FINISH
        BioelectricNodeLocalNumber = M_NODES_MAPPING%DOMAIN_LIST(BioelectricNodeIdx)
        BioelectricNodeGlobalNumber = M_NODES_MAPPING%LOCAL_TO_GLOBAL_MAP(BioelectricNodeLocalNumber)
        
        ! store global number in U(1)
        !CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
        !  & FIELD_VALUES_SET_TYPE,1,1,BioelectricNodeLocalNumber,1,DBLE(BioelectricNodeGlobalNumber),ERR,ERROR,*999)
          
        DofIdx=FIELD_VAR_IND_M_U2%COMPONENTS(7)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%&
          & NODES(BioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
        CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,DofIdx,DistanceToRightNode,ERR,ERROR,*999)
        
        DofIdx=FIELD_VAR_IND_M_U2%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%&
          & NODES(BioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
        CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,DofIdx,DistanceToLeftNode,ERR,ERROR,*999)
        
        PRINT *, COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR), BioelectricNodeGlobalNumber, DofIdx, &
          & DistanceToLeftNode, DistanceToRightNode
      ENDDO
    ENDIF
      
    IF (DEBUGGING.AND..FALSE.) THEN
      DISTRIBUTED_VECTOR=>INDEPENDENT_FIELD_MONODOMAIN%VARIABLE_TYPE_MAP(FIELD_U2_VARIABLE_TYPE)%PTR%PARAMETER_SETS% &
        & SET_TYPE(FIELD_VALUES_SET_TYPE)%PTR%PARAMETERS
      
      DO domain_idx = 1,DISTRIBUTED_VECTOR%DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
        DO i = 1,DISTRIBUTED_VECTOR%DOMAIN_MAPPING%ADJACENT_DOMAINS(domain_idx)%NUMBER_OF_RECEIVE_GHOSTS
          RecvIndex = DISTRIBUTED_VECTOR%DOMAIN_MAPPING%ADJACENT_DOMAINS(domain_idx)% &
            & LOCAL_GHOST_RECEIVE_INDICES(i)
          SendIndex = DISTRIBUTED_VECTOR%DOMAIN_MAPPING%ADJACENT_DOMAINS(domain_idx)% &
            & LOCAL_GHOST_SEND_INDICES(i)
          PRINT "(I1,(A,I3),4(A,I4))", COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR), &
            & ": Domain ", domain_idx, ", RecvIndex ", RecvIndex, " <- ", I, ", SendIndex ", SendIndex," -> ",I
        ENDDO
      ENDDO
      
    ENDIF
    
    ! bioelectric elements: transfer "distance between nodes" to ghost elements
    CALL FIELD_PARAMETER_SET_UPDATE_START(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
      & ERR,ERROR,*999)
    CALL FIELD_PARAMETER_SET_UPDATE_FINISH(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, & 
      & ERR,ERROR,*999)
    
    
    IF (DEBUGGING) THEN
      ! output all values
      ! distance (r) = physical distance to right neighouring node
      ! distance (l) = physical distance to left neighbouring node
      PRINT *, "process ",COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR), " after communication"
      PRINT *, "    process, glob. node, dof no.(l),        distance (l),              distance (r)"
      PRINT *, "internal"
      DO BioelectricNodeIdx = M_NODES_MAPPING%INTERNAL_START, M_NODES_MAPPING%INTERNAL_FINISH
        BioelectricNodeLocalNumber = M_NODES_MAPPING%DOMAIN_LIST(BioelectricNodeIdx)
        BioelectricNodeGlobalNumber = M_NODES_MAPPING%LOCAL_TO_GLOBAL_MAP(BioelectricNodeLocalNumber)
        
        DofIdx=FIELD_VAR_IND_M_U2%COMPONENTS(7)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%&
          & NODES(BioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
        CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,DofIdx,DistanceToRightNode,ERR,ERROR,*999)
        
        DofIdx=FIELD_VAR_IND_M_U2%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%&
          & NODES(BioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
        CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,DofIdx,DistanceToLeftNode,ERR,ERROR,*999)
        
        PRINT *, COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR),BioelectricNodeGlobalNumber, DofIdx, &
          & DistanceToLeftNode, DistanceToRightNode
      ENDDO
      PRINT *, "boundary"
      DO BioelectricNodeIdx = M_NODES_MAPPING%BOUNDARY_START, M_NODES_MAPPING%BOUNDARY_FINISH
        BioelectricNodeLocalNumber = M_NODES_MAPPING%DOMAIN_LIST(BioelectricNodeIdx)
        BioelectricNodeGlobalNumber = M_NODES_MAPPING%LOCAL_TO_GLOBAL_MAP(BioelectricNodeLocalNumber)
        
        DofIdx=FIELD_VAR_IND_M_U2%COMPONENTS(7)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%&
          & NODES(BioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
        CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,DofIdx,DistanceToRightNode,ERR,ERROR,*999)
        
        DofIdx=FIELD_VAR_IND_M_U2%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%&
          & NODES(BioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
        CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,DofIdx,DistanceToLeftNode,ERR,ERROR,*999)
        
        PRINT *, COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR),BioelectricNodeGlobalNumber, DofIdx, &
          & DistanceToLeftNode, DistanceToRightNode
      ENDDO
      PRINT *, "ghost"
      DO BioelectricNodeIdx = M_NODES_MAPPING%GHOST_START, M_NODES_MAPPING%GHOST_FINISH
        BioelectricNodeLocalNumber = M_NODES_MAPPING%DOMAIN_LIST(BioelectricNodeIdx)
        BioelectricNodeGlobalNumber = M_NODES_MAPPING%LOCAL_TO_GLOBAL_MAP(BioelectricNodeLocalNumber)
        
        DofIdx=FIELD_VAR_IND_M_U2%COMPONENTS(7)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%&
          & NODES(BioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
        CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,DofIdx,DistanceToRightNode,ERR,ERROR,*999)
        
        DofIdx=FIELD_VAR_IND_M_U2%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%&
          & NODES(BioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
        CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,DofIdx,DistanceToLeftNode,ERR,ERROR,*999)
        
        PRINT *, COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR),BioelectricNodeGlobalNumber, DofIdx, &
          & DistanceToLeftNode, DistanceToRightNode
      ENDDO
    ENDIF
    
    !PRINT *, "stop in bioelectric_finite_elasticity_routines.f90:2426"
    !STOP
    
    ! transfer geometric field of biolectric nodes 
    !CALL FIELD_PARAMETER_SET_UPDATE_START(GEOMETRIC_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
    !  & ERR,ERROR,*999)
    !CALL FIELD_PARAMETER_SET_UPDATE_FINISH(GEOMETRIC_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, & 
    !  & ERR,ERROR,*999)
     
    
    ! check on which side of fibre part ghost element is located
    DO BioelectricElementIdx = M_ELEMENTS_MAPPING%GHOST_START, M_ELEMENTS_MAPPING%GHOST_FINISH
      BioelectricElementLocalNumber = M_ELEMENTS_MAPPING%DOMAIN_LIST(BioelectricElementIdx)
      BioelectricElementGlobalNumber = M_ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(BioelectricElementLocalNumber)
      
      IF (SIZE(M_DOMAIN_TOPOLOGY_ELEMENTS%ELEMENTS(BioelectricElementLocalNumber)%ELEMENT_NODES) /= 2) THEN
        LOCAL_ERROR="Bioelectric Element with global number"//TRIM(NUMBER_TO_VSTRING(BioelectricElementGlobalNumber,"*",ERR, &
          & ERROR))//" does not have 2 nodes!"
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
      
      ! get local numbers of 2 adjacent nodes
      LeftBioelectricNodeLocalNumber = M_DOMAIN_TOPOLOGY_ELEMENTS%ELEMENTS(BioelectricElementLocalNumber)%ELEMENT_NODES(1)
      RightBioelectricNodeLocalNumber = M_DOMAIN_TOPOLOGY_ELEMENTS%ELEMENTS(BioelectricElementLocalNumber)%ELEMENT_NODES(2)
      
      LeftBioelectricNodeGlobalNumber = M_NODES_MAPPING%LOCAL_TO_GLOBAL_MAP(LeftBioelectricNodeLocalNumber)
      RightBioelectricNodeGlobalNumber = M_NODES_MAPPING%LOCAL_TO_GLOBAL_MAP(RightBioelectricNodeLocalNumber)
    
      ! if fibre continues on left side, update value of right node
      IF (M_ELEMENTS_TOPOLOGY%ELEMENTS(BioelectricElementLocalNumber)%ADJACENT_ELEMENTS(-1)%NUMBER_OF_ADJACENT_ELEMENTS == 1) THEN
        
        ! get the left adjacent element
        LeftBioelectricElementLocalNumber = M_ELEMENTS_TOPOLOGY%ELEMENTS(BioelectricElementLocalNumber)%ADJACENT_ELEMENTS(-1)%&
          & ADJACENT_ELEMENTS(1)
          
        ! get local number of left node
        SecondLeftBioelectricNodeLocalNumber = M_DOMAIN_TOPOLOGY_ELEMENTS%ELEMENTS(LeftBioelectricElementLocalNumber)%&
          & ELEMENT_NODES(1)
          
        ! get the 1D position of the left node of the left element
        DofIdx=FIELD_VAR_GEO_M%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%&
          & NODES(SecondLeftBioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
        CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(GEOMETRIC_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,DofIdx,Position1DSecondLeftNode,ERR,ERROR,*999)
        
        ! get distance between nodes, which is stored at that second left node
        DofIdx=FIELD_VAR_IND_M_U2%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%&
          & NODES(SecondLeftBioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
        CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,DofIdx,DistanceNodes,ERR,ERROR,*999)
        
        ! compute the 1D position of the left node
        Position1DLeftNode = Position1DSecondLeftNode + DistanceNodes
        
        ! store the computed position
        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(GEOMETRIC_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,1,1,LeftBioelectricNodeLocalNumber,1,Position1DLeftNode,ERR,ERROR,*999)
        
        ! get 1D position of left node (should already been computed by the local process [no, therefore the previous computation])
        !DofIdx=FIELD_VAR_GEO_M%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%&
        !  & NODES(LeftBioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
        !CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(GEOMETRIC_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
        !  & FIELD_VALUES_SET_TYPE,DofIdx,Position1DLeftNode,ERR,ERROR,*999)
    
        ! get distance between nodes, which is stored at the right node
        DofIdx=FIELD_VAR_IND_M_U2%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%&
          & NODES(RightBioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
        CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,DofIdx,DistanceNodes,ERR,ERROR,*999)
        
        IF (ABS(DistanceNodes) < 1e-15) THEN
            
          ! get distance between nodes, which is stored at the left node
          DofIdx=FIELD_VAR_IND_M_U2%COMPONENTS(7)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%&
            & NODES(LeftBioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,DofIdx,DistanceNodes,ERR,ERROR,*999)
          
          IF (DEBUGGING) PRINT *, "get d(r) from left node (global ", LeftBioelectricNodeGlobalNumber,&
            & ") DistanceNodes=",DistanceNodes
            
        ELSEIF (DEBUGGING) THEN 
          PRINT *, "get d(l) from right node, DistanceNodes=",DistanceNodes
        
        ENDIF
        
        ! compute position of right node
        Position1DRightNode = Position1DLeftNode + DistanceNodes
        
        ! store the new position of the right node to ghost node
        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(GEOMETRIC_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,1,1,RightBioelectricNodeLocalNumber,1,Position1DRightNode,ERR,ERROR,*999)

        IF (DEBUGGING) PRINT "(I2,3(A,F8.3),A,I3,A)", COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR), &
          & ": fibre continues left, set r=l+d(at r) ", Position1DRightNode,"=",&
          & Position1DLeftNode,"+",DistanceNodes,"(node global ",RightBioelectricNodeGlobalNumber,")"
          
      ELSE        ! fibre is on the right side
      
        ! get the right adjacent element
        RightBioelectricElementLocalNumber = M_ELEMENTS_TOPOLOGY%ELEMENTS(BioelectricElementLocalNumber)%ADJACENT_ELEMENTS(1)%&
          & ADJACENT_ELEMENTS(1)
          
       ! ! get local number of right node
       ! SecondRightBioelectricNodeLocalNumber = M_DOMAIN_TOPOLOGY_ELEMENTS%ELEMENTS(RightBioelectricElementLocalNumber)%&
       !   & ELEMENT_NODES(2)
       !   
       ! ! get the 1D position of the right node of the right element
       ! DofIdx=FIELD_VAR_GEO_M%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%&
       !   & NODES(SecondRightBioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
       ! CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(GEOMETRIC_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
       !   & FIELD_VALUES_SET_TYPE,DofIdx,Position1DSecondRightNode,ERR,ERROR,*999)
       ! 
       ! ! get distance between nodes, which is stored at the right node
       ! DofIdx=FIELD_VAR_IND_M_U2%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%&
       !   & NODES(RightBioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
       ! CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
       !   & FIELD_VALUES_SET_TYPE,DofIdx,DistanceNodes,ERR,ERROR,*999)
       ! 
       ! ! compute the 1D position of the left node
       ! Position1DRightNode = Position1DSecondRightNode - DistanceNodes
        
          ! get 1D position of right node (should already been computed by the local process [not sure therefore the computation]))
          DofIdx=FIELD_VAR_GEO_M%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%&
            & NODES(RightBioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(GEOMETRIC_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,DofIdx,Position1DRightNode,ERR,ERROR,*999)
      
        ! store the computed position
       ! CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(GEOMETRIC_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
       !   & FIELD_VALUES_SET_TYPE,1,1,RightBioelectricNodeLocalNumber,1,Position1DRightNode,ERR,ERROR,*999)
        
    
        !IF (DEBUGGING) PRINT *, "position 1D right node: r=rr-d(at r)=", Position1DRightNode,"=",Position1DSecondRightNode,&
        !  & "-",DistanceNodes," (was: ", TemporaryValue,")","(node global ",RightBioelectricNodeGlobalNumber,")"
    
        ! get distance between nodes, which is stored at the right node
        DofIdx=FIELD_VAR_IND_M_U2%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%&
          & NODES(RightBioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
        CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,DofIdx,DistanceNodes,ERR,ERROR,*999)
        
        IF (ABS(DistanceNodes) > 1e-5) THEN
        
          ! compute position of left node
          Position1DLeftNode = Position1DRightNode - DistanceNodes
          
          IF (DEBUGGING) PRINT "(I2,3(A,F8.3),A,I3,A)", COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR), &
            & ": fibre continues right, set l=r-d(at r) ", Position1DLeftNode, "=",&
            & Position1DRightNode,"-",DistanceNodes,"(node global ",RightBioelectricNodeGlobalNumber,")"
        ELSE
          
          ! get distance between nodes, which is stored at the left node
          DofIdx=FIELD_VAR_IND_M_U2%COMPONENTS(7)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%&
            & NODES(LeftBioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,DofIdx,DistanceNodes,ERR,ERROR,*999)
            
          ! compute position of left node
          Position1DLeftNode = Position1DRightNode - DistanceNodes
          
          IF (DEBUGGING) PRINT "(I2,3(A,F8.3),A,I3,A)", COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR), &
            & ": fibre continues right, d(at r)=0, l=r-d(at l) ",Position1DLeftNode, "=",& 
            & Position1DRightNode,"-",DistanceNodes,"(node global ",LeftBioelectricNodeGlobalNumber,")"
       ENDIF
          
        ! store the new position of the left node to ghost node
        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(GEOMETRIC_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,1,1,LeftBioelectricNodeLocalNumber,1,Position1DLeftNode,ERR,ERROR,*999)

      ENDIF
    ENDDO
    
    IF (DEBUGGING) THEN
      PRINT "(I1,A)", COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR), ": Updated values at boundary and ghost elements"
        
      ! show distance between nodes on ghost nodes
      DO BioelectricElementIdx = M_ELEMENTS_MAPPING%BOUNDARY_START, M_ELEMENTS_MAPPING%GHOST_FINISH
        BioelectricElementLocalNumber = M_ELEMENTS_MAPPING%DOMAIN_LIST(BioelectricElementIdx)
        BioelectricElementGlobalNumber = M_ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(BioelectricElementLocalNumber)
        
        PRINT *, "   Bioelectric element local", BioelectricElementLocalNumber, "global", BioelectricElementGlobalNumber
        
        IF (SIZE(M_DOMAIN_TOPOLOGY_ELEMENTS%ELEMENTS(BioelectricElementLocalNumber)%ELEMENT_NODES) == 2) THEN
          LeftBioelectricNodeLocalNumber = M_DOMAIN_TOPOLOGY_ELEMENTS%ELEMENTS(BioelectricElementLocalNumber)%ELEMENT_NODES(1)
          RightBioelectricNodeLocalNumber = M_DOMAIN_TOPOLOGY_ELEMENTS%ELEMENTS(BioelectricElementLocalNumber)%ELEMENT_NODES(2)
          
          LeftBioelectricNodeGlobalNumber = M_NODES_MAPPING%LOCAL_TO_GLOBAL_MAP(LeftBioelectricNodeLocalNumber)
          RightBioelectricNodeGlobalNumber = M_NODES_MAPPING%LOCAL_TO_GLOBAL_MAP(RightBioelectricNodeLocalNumber)

          PRINT *, "     left node, local",LeftBioelectricNodeLocalNumber, "global", LeftBioelectricNodeGlobalNumber
          PRINT *, "     right node, local",RightBioelectricNodeLocalNumber, "global", RightBioelectricNodeGlobalNumber
          
          ! get the position in 1D
          DofIdx=FIELD_VAR_GEO_M%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%&
            & NODES(LeftBioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(GEOMETRIC_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,DofIdx,Position1DLeftNode,ERR,ERROR,*999)                              
      
          ! get the position in 1D
          DofIdx=FIELD_VAR_GEO_M%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%&
            & NODES(RightBioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(GEOMETRIC_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,DofIdx,Position1DRightNode,ERR,ERROR,*999)    
      
          ! get the distance to the left node from the left node
          DofIdx=FIELD_VAR_IND_M_U2%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%&
            & NODES(LeftBioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,DofIdx,DistanceToLeftNode,ERR,ERROR,*999)                              
      
          ! get the distance to the left node from the right node
          DofIdx=FIELD_VAR_IND_M_U2%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%&
            & NODES(RightBioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,DofIdx,DistanceToRightNode,ERR,ERROR,*999)                           
      
          PRINT *, "      Position1DLeftNode:", Position1DLeftNode
          PRINT *, "      Position1DRightNode:", Position1DRightNode                             
          PRINT *, "      DistanceToLeftNode distance (l) (at l):", DistanceToLeftNode                             
          PRINT *, "      DistanceToRightNode distance (l) (at r):", DistanceToRightNode                             
      
          
        ELSE        
          ! if the element only contains one node, take this (should not happen)
          LeftBioelectricNodeLocalNumber = M_DOMAIN_TOPOLOGY_ELEMENTS%ELEMENTS(BioelectricElementLocalNumber)%ELEMENT_NODES(1)
          PRINT *, "      only 1 node, local",LeftBioelectricNodeLocalNumber
        ENDIF
      ENDDO
      
      PRINT *, "node  local       global   1D position"
      DO BioelectricNodeIdx = M_NODES_MAPPING%INTERNAL_START, M_NODES_MAPPING%GHOST_FINISH
      
        LeftBioelectricNodeLocalNumber = M_NODES_MAPPING%DOMAIN_LIST(BioelectricNodeIdx)
        LeftBioelectricNodeGlobalNumber = M_NODES_MAPPING%LOCAL_TO_GLOBAL_MAP(LeftBioelectricNodeLocalNumber)
          
        ! get the position in 1D
        DofIdx=FIELD_VAR_GEO_M%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%&
          & NODES(LeftBioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
        CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(GEOMETRIC_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,DofIdx,Position1DLeftNode,ERR,ERROR,*999)   
          
          PRINT *, LeftBioelectricNodeLocalNumber, LeftBioelectricNodeGlobalNumber, Position1DLeftNode
      ENDDO
      
    ENDIF

#if 0
    IF (CALL_NO == 2) THEN
      PRINT *, "stop in bioelectric_finite_elastcity_routines.f90:2704"
      STOP
    ENDIF
#endif
    
    EXITS("BioelectricFiniteElasticity_GhostElements")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_GhostElements",ERR,ERROR)
    EXITS("BioelectricFiniteElasticity_GhostElements")
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_GhostElements
  
  !
  !================================================================================================================================
  !
  
  !> Check if value is NAN. This is not used so far.
  FUNCTION IS_NAN(Value)
    REAL(DP), INTENT(IN) :: Value
    LOGICAL :: IS_NAN
    CHARACTER(LEN=3) :: STR
    
    ! normally it should work like this, but is does not! (tested with GCC 5.4)
    IF (Value /= Value) THEN
      PRINT *, "Is nan (by comparison)!"
      IS_NAN = .TRUE.
      RETURN
    ENDIF
    ! this is a GCC extension, but does not work either
    IF (ISNAN(Value)) THEN
      PRINT *, "Is nan (by gcc)!"
      IS_NAN = .TRUE.
      RETURN
    ENDIF
    
    ! this works (but not with the flags "-ffpe-trap=invalid,zero" that make floating point exceptions fatal, which is good in debug mode)
    Write(Str, "(F3.5)") Value
    !Print*, "STR=[", TRIM(STR), "]"
    IF (Str == "NaN") THEN
      PRINT *, "Is nan (by string comparison)"
      IS_NAN = .TRUE.
      RETURN
    ENDIF
    IS_NAN = .FALSE.
    
  END FUNCTION IS_NAN
  
  !
  !================================================================================================================================
  !

  !>Update the bioelectric equation geometric field from the finite elasticity dependent field (deformed geometry)
  !>NOTE: this is only temporary - will be replaced once embedded meshes are available (if ever)
  !>This methods gets called FortranExample before the simulation starts with CALC_CLOSEST_GAUSS_POINT=.TRUE.
  !> Then it gets called in BioelectricFiniteElasticity_ControlLoopPostLoop and again after the while loop converged with CALC_CLOSEST_GAUSS_POINT=.FALSE.
  SUBROUTINE BioelectricFiniteElasticity_UpdateGeometricField(CONTROL_LOOP,CALC_CLOSEST_GAUSS_POINT,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    LOGICAL, INTENT(IN) :: CALC_CLOSEST_GAUSS_POINT !<If true then the closest finite elasticity Gauss point for each bioelectrics node is calculated. This is to be set .TRUE. only once (in the first call)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP_ROOT,CONTROL_LOOP_PARENT,CONTROL_LOOP_ELASTICITY,CONTROL_LOOP_MONODOMAIN
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(FIELD_TYPE), POINTER :: INDEPENDENT_FIELD_ELASTICITY,GEOMETRIC_FIELD_MONODOMAIN,GEOMETRIC_FIELD_ELASTICITY
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD_MONODOMAIN,INDEPENDENT_FIELD_MONODOMAIN,DEPENDENT_FIELD_ELASTICITY
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING_ELASTICITY, SOLVER_MAPPING_MONODOMAIN
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: INTERPOLATION_PARAMETERS
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: M_NODES_MAPPING, M_ELEMENTS_MAPPING
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: FE_ELEMENTS_TOPOLOGY
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: M_ELEMENTS_TOPOLOGY
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VAR_DEP_M,FIELD_VAR_GEO_M,FIELD_VAR_IND_FE,FIELD_VAR_IND_M_U1,FIELD_VAR_IND_M_U2
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VAR_IND_M_V
    TYPE(DOMAIN_GLOBAL_MAPPING_TYPE),POINTER :: GLOBAL_TO_LOCAL_MAP
    INTEGER(INTG) :: FEElementGlobalNumber,FEElementLocalNumber,FibreStartsInCurrentElement
    INTEGER(INTG) :: i, j, ElementGlobalNumber
    INTEGER(INTG) :: GAUSS_POINT,gauss_idx,FibreIdx,DomainIdx
    INTEGER(INTG) :: nodes_in_Xi_1_small
    INTEGER(INTG) :: DofIdx
    REAL(DP) :: VALUE
    REAL(DP) :: DISTANCE,ContractionVelocity,MaximumAllowedShorteningVelocity
    REAL(DP) :: OLD_DIST,OLD_DIST_2,OLD_DIST_3,OLD_DIST_4
    REAL(DP) :: HalfSarcomereLength
    REAL(DP) :: DELTA_XI1
    REAL(DP) :: XI(3),XI_DEBUG(3),InitialNodeMDistance,HalfSarcomereInitialLength,TIME_STEP,DIST
    REAL(DP) :: Position3DCurrentNode(3), Position3DLeftNode(3), Position3DRightNode(3)
    REAL(DP) :: RelativeContractionVelocity
    REAL(DP) :: Position1DCurrentNode, Position1DLeftNode
    REAL(DP), POINTER :: GAUSS_POSITIONS(:,:)
    LOGICAL :: IsFinished, IsFirstBioelectricNodeOfFibre
    LOGICAL :: LeftNeighbourIsGhost, PreviousNodeHadNoLeftNeighbour
    INTEGER(INTG) :: BioelectricNodeLocalNumber, BioelectricNodeGlobalNumber, BioelectricNodeIdx
    TYPE(DOMAIN_TYPE), POINTER :: MONODOMAIN_DOMAIN
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: MONODOMAIN_MAPPINGS
    INTEGER(INTG) :: NumberNodesInCurrentFEElement
    REAL(DP) :: RegionWidth, FibrePhysicalLength
    INTEGER(INTG) :: BioelectricNodeInFibreNumber
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: M_DOMAIN_TOPOLOGY
    TYPE(DOMAIN_NODES_TYPE), POINTER :: M_DOMAIN_TOPOLOGY_NODES
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: M_DOMAIN_TOPOLOGY_ELEMENTS

    !LOGICAL, PARAMETER :: DEBUGGING = .FALSE.   ! enable debugging output with this parameter
    LOGICAL :: DEBUGGING = .FALSE.   ! enable debugging output with this parameter
    INTEGER(Intg) :: MPI_IERROR, LeftNeighbourBioelectricNodeLocalNumber
    INTEGER(Intg) :: ComputationalNodeNumber, NumberOfComputationalNodes

    ENTERS("BioelectricFiniteElasticity_UpdateGeometricField",ERR,ERROR,*999)

    ! this needs to be called
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,ComputationalNodeNumber,MPI_IERROR)
    !CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NumberOfComputationalNodes,MPI_IERROR)
#if 0
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,ComputationalNodeNumber,MPI_IERROR)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NumberOfComputationalNodes,MPI_IERROR)

    ! turn on debugging depending on computational node number
    IF (NumberOfComputationalNodes == 4 .AND. ComputationalNodeNumber == 0) THEN
      DEBUGGING = .TRUE.      ! note when line searching for 'debugging = .true.': this line is eventually commented out
      PRINT*, "ComputationalNodeNumber=",ComputationalNodeNumber,"of",NumberOfComputationalNodes
    ENDIF
#endif
    
    NULLIFY(CONTROL_LOOP_ROOT)
    NULLIFY(CONTROL_LOOP_PARENT)
    NULLIFY(CONTROL_LOOP_ELASTICITY)
    NULLIFY(CONTROL_LOOP_MONODOMAIN)
    NULLIFY(PROBLEM)
    NULLIFY(SOLVERS)
    NULLIFY(SOLVER)
    NULLIFY(INDEPENDENT_FIELD_ELASTICITY)
    NULLIFY(DEPENDENT_FIELD_MONODOMAIN)
    NULLIFY(INDEPENDENT_FIELD_MONODOMAIN)
    NULLIFY(DEPENDENT_FIELD_ELASTICITY)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVER_MAPPING_MONODOMAIN)
    NULLIFY(SOLVER_MAPPING_ELASTICITY)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(GEOMETRIC_FIELD_MONODOMAIN)
    NULLIFY(GEOMETRIC_FIELD_ELASTICITY)
    NULLIFY(FE_ELEMENTS_TOPOLOGY)
    NULLIFY(INTERPOLATED_POINT)
    NULLIFY(INTERPOLATION_PARAMETERS)
    NULLIFY(FIELD_VAR_DEP_M)
    NULLIFY(FIELD_VAR_GEO_M)
    NULLIFY(FIELD_VAR_IND_FE)
    NULLIFY(FIELD_VAR_IND_M_U1)
    NULLIFY(FIELD_VAR_IND_M_U2)
    NULLIFY(FIELD_VAR_IND_M_V)
    NULLIFY(GAUSS_POSITIONS)
    
    IF (DEBUGGING) THEN
      PRINT*, "BioelectricFiniteElasticity_UpdateGeometricField, CALC_CLOSEST_GAUSS_POINT=",CALC_CLOSEST_GAUSS_POINT
    ENDIF
    
    !CALL gdbParallelDebuggingBarrier()
    
    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
        PROBLEM=>CONTROL_LOOP%PROBLEM
        IF(ASSOCIATED(PROBLEM)) THEN
          IF(.NOT.ALLOCATED(problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a bioelectric-finite elasticity problem.", &
              & err,error,*999)
          END IF
          SELECT CASE(PROBLEM%SPECIFICATION(2))
          CASE(PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE)
            SELECT CASE(PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE, &
              & PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE, &
              & PROBLEM_GUDUNOV_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
              & PROBLEM_STRANG_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE, &
              & PROBLEM_STRANG_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE, &
              & PROBLEM_STRANG_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE)
              
              ! this case
              CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
              CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP_PARENT,ERR,ERROR,*999)
              !get the monodomain sub loop, solvers, solver, and finally geometric field and dependent field
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP_PARENT,1,CONTROL_LOOP_MONODOMAIN,ERR,ERROR,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP_MONODOMAIN,SOLVERS,ERR,ERROR,*999)
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING_MONODOMAIN=>SOLVER_EQUATIONS%SOLVER_MAPPING
                
                IF (DEBUGGING) THEN
                  PRINT *, "SolverMapping has ", SOLVER_MAPPING_MONODOMAIN%NUMBER_OF_EQUATIONS_SETS, "equation sets"
                ENDIF
                
                IF(ASSOCIATED(SOLVER_MAPPING_MONODOMAIN)) THEN
                  EQUATIONS_SET=>SOLVER_MAPPING_MONODOMAIN%EQUATIONS_SETS(1)%PTR
                  IF(ASSOCIATED(EQUATIONS_SET)) THEN
                    GEOMETRIC_FIELD_MONODOMAIN=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                    IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD_MONODOMAIN)) THEN
                      CALL FlagError("Geometric field is not associated.",ERR,ERROR,*999)
                    ENDIF
                    ! the Field_V_Variable_Type contains the 3D nodal positions
                    DEPENDENT_FIELD_MONODOMAIN=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                    IF(.NOT.ASSOCIATED(DEPENDENT_FIELD_MONODOMAIN)) THEN
                      CALL FlagError("Dependent field is not associated.",ERR,ERROR,*999)
                    ENDIF
                    INDEPENDENT_FIELD_MONODOMAIN=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
                    IF(.NOT.ASSOCIATED(INDEPENDENT_FIELD_MONODOMAIN)) THEN
                      CALL FlagError("Independent field is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Solver mapping is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Solver equations is not associated.",ERR,ERROR,*999)
              ENDIF
              
              NULLIFY(SOLVERS)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_MAPPING_ELASTICITY)
              NULLIFY(EQUATIONS_SET)
              NULLIFY(SOLVER_EQUATIONS)
              
              !get the finite elasticity sub loop, solvers, solver, and finally the dependent and independent fields
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP_PARENT,2,CONTROL_LOOP_ELASTICITY,ERR,ERROR,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP_ELASTICITY,SOLVERS,ERR,ERROR,*999)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING_ELASTICITY=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING_ELASTICITY)) THEN
                  EQUATIONS_SET=>SOLVER_MAPPING_ELASTICITY%EQUATIONS_SETS(1)%PTR
                  IF(ASSOCIATED(EQUATIONS_SET)) THEN
                    INDEPENDENT_FIELD_ELASTICITY=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
                    IF(.NOT.ASSOCIATED(INDEPENDENT_FIELD_ELASTICITY)) THEN
                      CALL FlagError("Independent field is not associated.",ERR,ERROR,*999)
                    ENDIF
                    GEOMETRIC_FIELD_ELASTICITY=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                    IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD_ELASTICITY)) THEN
                      CALL FlagError("Dependent field is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Solver mapping is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Solver equations is not associated.",ERR,ERROR,*999)
              ENDIF

              ! get field variables
              CALL FIELD_VARIABLE_GET(GEOMETRIC_FIELD_MONODOMAIN,  FIELD_U_VARIABLE_TYPE, FIELD_VAR_GEO_M,   ERR,ERROR,*999)
              CALL FIELD_VARIABLE_GET(DEPENDENT_FIELD_MONODOMAIN,  FIELD_V_VARIABLE_TYPE, FIELD_VAR_DEP_M,   ERR,ERROR,*999)
              CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U1_VARIABLE_TYPE,FIELD_VAR_IND_M_U1,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE,FIELD_VAR_IND_M_U2,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, FIELD_VAR_IND_M_V, ERR,ERROR,*999)
              CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD_ELASTICITY,FIELD_V_VARIABLE_TYPE, FIELD_VAR_IND_FE,  ERR,ERROR,*999)

              MONODOMAIN_DOMAIN=>GEOMETRIC_FIELD_MONODOMAIN%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD_MONODOMAIN%DECOMPOSITION% &
                & MESH_COMPONENT_NUMBER)%PTR              ! DOMAIN_TYPE    
              MONODOMAIN_MAPPINGS=>MONODOMAIN_DOMAIN%MAPPINGS
              
              M_NODES_MAPPING=>MONODOMAIN_MAPPINGS%NODES
              M_ELEMENTS_MAPPING=>MONODOMAIN_MAPPINGS%ELEMENTS
              
              FE_ELEMENTS_TOPOLOGY=>GEOMETRIC_FIELD_ELASTICITY%DECOMPOSITION%TOPOLOGY%ELEMENTS
              M_ELEMENTS_TOPOLOGY=>GEOMETRIC_FIELD_MONODOMAIN%DECOMPOSITION%TOPOLOGY%ELEMENTS

              M_DOMAIN_TOPOLOGY=>MONODOMAIN_DOMAIN%TOPOLOGY     ! TYPE DOMAIN_TOPOLOGY_TYPE
              M_DOMAIN_TOPOLOGY_NODES=>M_DOMAIN_TOPOLOGY%NODES  ! TYPE DOMAIN_NODES_TYPE
              M_DOMAIN_TOPOLOGY_ELEMENTS=>M_DOMAIN_TOPOLOGY%ELEMENTS  ! TYPE DOMAIN_ELEMENTS_TYPE
              
              ! --------- get constants --------------------
              ! get number of in series fibres (not used anymore, set NumberInSeriesFibres to 1)
              !DofIdx=FIELD_VAR_IND_FE%COMPONENTS(5)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
              !CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              !  & DofIdx,NumberInSeriesFibres,ERR,ERROR,*999)
  
              !get the initial sarcomere half length
              DofIdx=FIELD_VAR_IND_M_U1%COMPONENTS(2)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP                          
              CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U1_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,DofIdx,HalfSarcomereInitialLength,ERR,ERROR,*999)
              
              !get initial node distance
              DofIdx=FIELD_VAR_IND_M_U1%COMPONENTS(3)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP                         
              CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U1_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,DofIdx,InitialNodeMDistance,ERR,ERROR,*999)
                        
              ! get physical length of fibre
              RegionWidth = FIELD_VAR_IND_FE%REGION%GENERATED_MESHES%GENERATED_MESHES(1)%PTR%&
                & REGULAR_MESH%MAXIMUM_EXTENT(1)
              !FibrePhysicalLength = RegionWidth / NumberInSeriesFibres
              FibrePhysicalLength = RegionWidth
              
              
              IF (DEBUGGING) THEN
                PRINT*, "FibrePhysicalLength: ", FibrePhysicalLength
              ENDIF
              
              !get the maximum contraction velocity 
              DofIdx=FIELD_VAR_IND_M_U2%COMPONENTS(2)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
              CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,DofIdx,MaximumAllowedShorteningVelocity,ERR,ERROR,*999)
              !NOTE: MaximumAllowedShorteningVelocity is the max shortening velocity, and hence negative!!!
              !The max lengthening velocity is assumed to be   abs(MaximumAllowedShorteningVelocity)/2.0
              
              !get the time step of the elasticity problem
              TIME_STEP=CONTROL_LOOP_PARENT%TIME_LOOP%TIME_INCREMENT
              
              ! initialize counters
              BioelectricNodeIdx=M_NODES_MAPPING%INTERNAL_START 
              nodes_in_Xi_1_small = 0
                
              ! debugging output
              IF (DEBUGGING) THEN
                PRINT*, "========== fibre beginning output, Process ",ComputationalNodeNumber," =========="
                  
                DO i = 1, FE_ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS
                      
                  ! get local element number from global element number
                  GLOBAL_TO_LOCAL_MAP=>FE_ELEMENTS_TOPOLOGY%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD_ELASTICITY%DECOMPOSITION% &
                    & MESH_COMPONENT_NUMBER)%PTR%MAPPINGS%ELEMENTS%GLOBAL_TO_LOCAL_MAP(i)
                  
                  DomainIdx = 0
                  DO j = 1, GLOBAL_TO_LOCAL_MAP%NUMBER_OF_DOMAINS
                    IF (GLOBAL_TO_LOCAL_MAP%DOMAIN_NUMBER(j) == ComputationalNodeNumber) THEN
                      DomainIdx = j
                      EXIT
                    ENDIF
                  ENDDO
                  
                  IF (DomainIdx == 0) CYCLE
                  
                  FEElementLocalNumber = GLOBAL_TO_LOCAL_MAP%LOCAL_NUMBER(DomainIdx)
                  ElementGlobalNumber = FE_ELEMENTS_TOPOLOGY%ELEMENTS(FEElementLocalNumber)%GLOBAL_NUMBER
                    
                  !cycle if element is not on local domain
                  IF (FE_ELEMENTS_TOPOLOGY%DECOMPOSITION%ELEMENT_DOMAIN(ElementGlobalNumber) /= ComputationalNodeNumber) THEN
                    IF (DEBUGGING) THEN
                      PRINT*, "Element idx ",i," local number ",FEElementLocalNumber, & 
                       & " is not on local domain of computational node, skip"
                    ENDIF
                    CYCLE
                  ENDIF
                
                  DofIdx=FIELD_VAR_IND_FE%COMPONENTS(4)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%ELEMENTS(FEElementLocalNumber)
                  CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                    & DofIdx,FibreStartsInCurrentElement,ERR,ERROR,*999)
                  
                    PRINT*, "Element idx ", i, " local number ",FEElementLocalNumber, &
                    & "contains fibre beginning: ",FibreStartsInCurrentElement
                ENDDO
                PRINT*, "==========================================================="                  
              ENDIF
              
              IF (DEBUGGING) THEN
                PRINT*, "topology statistics: NUMBER_OF_ELEMENTS: ",FE_ELEMENTS_TOPOLOGY%NUMBER_OF_ELEMENTS, &
                  & ", TOTAL_NUMBER_OF_ELEMENTS:",FE_ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS, &
                  & ", NUMBER_OF_GLOBAL_ELEMENTS", FE_ELEMENTS_TOPOLOGY%NUMBER_OF_GLOBAL_ELEMENTS               
                PRINT*, "loop over ",FE_ELEMENTS_TOPOLOGY%NUMBER_OF_ELEMENTS," 3D finite elasticity elements"
                !CALL Print_DECOMPOSITION_ELEMENTS(FE_ELEMENTS_TOPOLOGY, 2, 10)
              ENDIF
              
              !EQUATIONS_SET_MONODOMAIN%geometry%geometric_field%decomposition%domain(1)%ptr%mappings%nodes
              ! Initialize Monodomain iteration routine
              CALL BioelectricFiniteElasticity_IterateNextMonodomainNode(SOLVER_MAPPING_MONODOMAIN, GEOMETRIC_FIELD_ELASTICITY, &
                 & INDEPENDENT_FIELD_ELASTICITY, .TRUE., &
                 & GEOMETRIC_FIELD_MONODOMAIN, DEPENDENT_FIELD_MONODOMAIN, INDEPENDENT_FIELD_MONODOMAIN, &
                 & FIELD_VAR_DEP_M, FIELD_VAR_GEO_M, FIELD_VAR_IND_M_U1, FIELD_VAR_IND_M_U2, FIELD_VAR_IND_M_V, &
                 & IsFinished, FEElementGlobalNumber, FEElementLocalNumber, BioelectricNodeGlobalNumber, &
                 & BioelectricNodeLocalNumber, BioelectricNodeInFibreNumber, LeftNeighbourBioelectricNodeLocalNumber, FibreIdx, &
                 & XI, DELTA_XI1, IsFirstBioelectricNodeOfFibre, LeftNeighbourIsGhost, PreviousNodeHadNoLeftNeighbour, &
                 & ERR, ERROR, *999)

              ! synchronously iterate over monodomain elements and 3D finite elasticity elements (only internal and boundary elements, not ghost elements)
              DO
                ! debugging output
                IF (DEBUGGING) THEN
                  PRINT *, ""
                  PRINT "(I1,6(A,I4))", ComputationalNodeNumber, ": FE element, global:", FEElementGlobalNumber, &
                   & ", local: ", FEElementLocalNumber, &
                   & ", Bioelectric node, global: ", BioelectricNodeGlobalNumber, ", local: ", BioelectricNodeLocalNumber, &
                   & ", previous: ", LeftNeighbourBioelectricNodeLocalNumber, &
                   & ", Fibre index: ", FibreIdx
                  PRINT "(3(A,F10.5),A)", "    Xi: (", XI(1), ",", XI(2), ",", XI(3), ")"
                  !PRINT *, "       previous bioelectric node, local: ",LeftNeighbourBioelectricNodeLocalNumber 
                  IF (IsFirstBioelectricNodeOfFibre) THEN
                    PRINT *, "(first bioelectric node of fibre)"
                  ENDIF
                ENDIF
          
                ! --------------------- determine 1D and 3D position of current node ------------------
                ! get the finite elasticity dependent field interpolation parameters of this element
                INTERPOLATION_PARAMETERS=>EQUATIONS_SET%EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS &
                  & (FIELD_U_VARIABLE_TYPE)%PTR
                CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,FEElementLocalNumber, &
                  & INTERPOLATION_PARAMETERS, ERR,ERROR,*999)
                INTERPOLATED_POINT=>EQUATIONS_SET%EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR

                ! find the interpolated position of the bioelectric grid node from the finite elasticity FE dependent field
                ! XI goes from 0 to 1 per FE element
                CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,XI,INTERPOLATED_POINT,ERR,ERROR,*999)
                Position3DCurrentNode(:) = INTERPOLATED_POINT%VALUES(1:3,1)
                                    
                ! store the spatial position of the bioelectric node to DEPENDENT_FIELD_MONODOMAIN, V ("GeometryM3D")
                ! it is used for visualization only
                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(DEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,1,1,BioelectricNodeLocalNumber,1,Position3DCurrentNode(1),ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(DEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,1,1,BioelectricNodeLocalNumber,2,Position3DCurrentNode(2),ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(DEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,1,1,BioelectricNodeLocalNumber,3,Position3DCurrentNode(3),ERR,ERROR,*999)
       
                IF (DEBUGGING) PRINT *, "LeftNeighbourBioelectricNodeLocalNumber: ", LeftNeighbourBioelectricNodeLocalNumber, &
                  & "LeftNeighbourIsGhost: ", LeftNeighbourIsGhost, "DELTA_XI1:", DELTA_XI1, " XI(1)-DELTA_XI1=", (XI(1)-DELTA_XI1)
                  
                ! if a previous node is available, which was processed earlier
                IF (LeftNeighbourBioelectricNodeLocalNumber /= 0 .AND..NOT. LeftNeighbourIsGhost) THEN
                 
                  IF(DEBUGGING) PRINT *, "   Get stored position of left node"
                 
                  ! get the stored position in 3D of the previous node
                  DofIdx=FIELD_VAR_DEP_M%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                    & NODES(LeftNeighbourBioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
                  CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,DofIdx,Position3DLeftNode(1),ERR,ERROR,*999)
                  
                  DofIdx=FIELD_VAR_DEP_M%COMPONENTS(2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                    & NODES(LeftNeighbourBioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
                  CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,DofIdx,Position3DLeftNode(2),ERR,ERROR,*999)
                  
                  DofIdx=FIELD_VAR_DEP_M%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                    & NODES(LeftNeighbourBioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
                  CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(DEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,DofIdx,Position3DLeftNode(3),ERR,ERROR,*999)
                  
                  ! get the stored position in 1D of the previous node
                  DofIdx=FIELD_VAR_GEO_M%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%&
                    & NODES(LeftNeighbourBioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
                  CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(GEOMETRIC_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,DofIdx,Position1DLeftNode,ERR,ERROR,*999)                              
              
                ! if the left neighbour node is still positioned inside the current FE element (but yet on another domain)
                ELSEIF (XI(1)-DELTA_XI1 >= -1.0E-13_DP) THEN
                  
                  IF(DEBUGGING) PRINT *, "   Compute position of left node"
                  
                  ! also compute the position of that node
                  XI(1) = XI(1)-DELTA_XI1
                  
                  ! find the interpolated position of the left neighbour bioelectric grid node from the finite elasticity FE dependent field
                  ! XI goes from 0 to 1 per FE element
                  CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,XI,INTERPOLATED_POINT,ERR,ERROR,*999)
                  Position3DLeftNode(:) = INTERPOLATED_POINT%VALUES(1:3,1)
                  
                  ! set the 1D position of the left node to 0
                  Position1DLeftNode = 0.0_DP
                
                ! if the left neighbour node really is on a different domain, set its position to the current position
                ELSE  
                
                  IF(DEBUGGING) PRINT *, "   Set position of left node to current position"
                
                  Position3DLeftNode(:) = Position3DCurrentNode(:)
                  Position1DLeftNode = 0.0_DP
                ENDIF
                    
                IF (DEBUGGING) THEN
                  PRINT *, "Position3DCurrentNode :", Position3DCurrentNode
                  PRINT *, "Position3DLeftNode:    ", Position3DLeftNode
                ENDIF
  
                  
                ! compute the distance between the previous node and the current node
                DIST=SQRT( &
                  & (Position3DCurrentNode(1)-Position3DLeftNode(1)) * (Position3DCurrentNode(1)-Position3DLeftNode(1))+ &
                  & (Position3DCurrentNode(2)-Position3DLeftNode(2)) * (Position3DCurrentNode(2)-Position3DLeftNode(2))+ &
                  & (Position3DCurrentNode(3)-Position3DLeftNode(3)) * (Position3DCurrentNode(3)-Position3DLeftNode(3)))
  
                ! update the current sarcomere half length
                HalfSarcomereLength = HalfSarcomereInitialLength * DIST / InitialNodeMDistance
                
                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,1,1,BioelectricNodeLocalNumber,1,HalfSarcomereLength,ERR,ERROR,*999)

                ! ---------------- compute contraction velocity -----------------
                !get the distance between current and left node in the previous time step
                DofIdx=FIELD_VAR_IND_M_U2%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                  & NODES(BioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,DofIdx,OLD_DIST,ERR,ERROR,*999)

                !get the distance between current and left node nodes before 2 time steps
                DofIdx=FIELD_VAR_IND_M_U2%COMPONENTS(4)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                  & NODES(BioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,DofIdx,OLD_DIST_2,ERR,ERROR,*999)

                !get the distance between current and left node nodes before 3 time steps
                DofIdx=FIELD_VAR_IND_M_U2%COMPONENTS(5)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                  & NODES(BioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,DofIdx,OLD_DIST_3,ERR,ERROR,*999)

                !get the distance between current and left node nodes before 4 time steps
                DofIdx=FIELD_VAR_IND_M_U2%COMPONENTS(6)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                  & NODES(BioelectricNodeLocalNumber)%DERIVATIVES(1)%VERSIONS(1)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,DofIdx,OLD_DIST_4,ERR,ERROR,*999)
  
                IF (DEBUGGING) PRINT "(5(A, F6.3))", "Previous node distances: ", DIST, ", ", OLD_DIST, ", ", OLD_DIST_2, ", ", &
                  & OLD_DIST_3, ", ", OLD_DIST_4
  
                ! compute the new contraction velocity
                ! ContractionVelocity=(DIST-OLD_DIST)/TIME_STEP
                ! v = 1/4 * [(dist-dist1)/dt + (dist-dist2)/(2dt) + (dist-dist3)/(3dt) + (dist-dist4)/(4dt)]
                ContractionVelocity=0.25_DP*((DIST-OLD_DIST)/TIME_STEP+(DIST-OLD_DIST_2)/(2.0_DP*TIME_STEP)+ &
                  & (DIST-OLD_DIST_3)/(3.0_DP*TIME_STEP)+(DIST-OLD_DIST_4)/(4.0_DP*TIME_STEP))
        
                ! CALC_CLOSEST_GAUSS_POINT is only the first time .TRUE. and subsequent times .FALSE.
                IF (.NOT. CALC_CLOSEST_GAUSS_POINT) THEN
                
                  ! MaximumAllowedShorteningVelocity is the maximum shortening velocity and hence negative
                  IF (ContractionVelocity<MaximumAllowedShorteningVelocity) THEN
                    !CALL FLAG_WARNING('Exceeded maximum contraction velocity (shortening).',ERR,ERROR,*999)
                    !CALL FLAG_WARNING('  ContractionVelocity: '// &
                    !  TRIM(NUMBER_TO_VSTRING(ContractionVelocity,"*",ERR,ERROR))//", MaximumAllowedShorteningVelocity:"// &
                    !  TRIM(NUMBER_TO_VSTRING(MaximumAllowedShorteningVelocity,"*",ERR,ERROR)),ERR,ERROR,*999)
                    ContractionVelocity=MaximumAllowedShorteningVelocity
                  !The max lengthening velocity is assumed to be MaximumAllowedShorteningVelocity/2.0
                  ELSEIF(ContractionVelocity>(-MaximumAllowedShorteningVelocity/2.0_DP)) THEN
                    ! warning disabled
                    !CALL FLAG_WARNING('Exceeded maximum contraction velocity (lengthening).',ERR,ERROR,*999)
                    !CALL FLAG_WARNING('Exceeded maximum contraction velocity (lengthening).',ERR,ERROR,*999)
                    !CALL FLAG_WARNING('  ContractionVelocity: '// &
                    !  TRIM(NUMBER_TO_VSTRING(ContractionVelocity,"*",ERR,ERROR))//", MaximumAllowedLengtheningVelocity:"// &
                    !  TRIM(NUMBER_TO_VSTRING(-MaximumAllowedShorteningVelocity/2.0_DP,"*",ERR,ERROR)),ERR,ERROR,*999)
                    ContractionVelocity=-MaximumAllowedShorteningVelocity/2.0_DP
                  ENDIF
                ENDIF
    
                RelativeContractionVelocity = ContractionVelocity/ABS(MaximumAllowedShorteningVelocity)
    
                ! store node distance to previous node
                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,1,1,BioelectricNodeLocalNumber,1,DIST,ERR,ERROR,*999)

                IF (DEBUGGING) PRINT *, "Store d(l)=",DIST," to node global ",BioelectricNodeGlobalNumber
                
                ! update distances for old timesteps
                ! store the relative contraction velocity in component 3 of the U2 variable of the monodomain independent field
                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,1,1,BioelectricNodeLocalNumber,3,RelativeContractionVelocity,ERR,ERROR,*999)

                ! store node distance to previous node before 1 time step
                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,1,1,BioelectricNodeLocalNumber,4,OLD_DIST,ERR,ERROR,*999)

                ! store node distance to previous node before 2 time step
                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,1,1,BioelectricNodeLocalNumber,5,OLD_DIST_2,ERR,ERROR,*999)

                ! store node distance to previous node before 3 time step
                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,1,1,BioelectricNodeLocalNumber,6,OLD_DIST_3,ERR,ERROR,*999)

                ! update the current 1D node position
                Position1DCurrentNode = Position1DLeftNode+DIST
                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(GEOMETRIC_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,1,1,BioelectricNodeLocalNumber,1,Position1DCurrentNode,ERR,ERROR,*999)

                IF (DEBUGGING) THEN
                  PRINT "(I1,2(A,F8.5))", ComputationalNodeNumber, ": prev node has pos ", Position1DLeftNode, ", DIST=", DIST
                  PRINT "(I1,A,I4,3(A,F8.5))", ComputationalNodeNumber, ": node global ", BioelectricNodeGlobalNumber, " pos: ", &
                    & Position1DCurrentNode, " dist", DIST, " |xi(1)+delta_xi1-1|=",ABS(XI(1) + DELTA_XI1 - 1.0)
                ENDIF
                  
                ! if we are currently at the second node where left previous node had no neighbour, also store the current quantities to the left node
                IF (PreviousNodeHadNoLeftNeighbour) THEN
                  ! current sarcomere half length
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U1_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,1,1,LeftNeighbourBioelectricNodeLocalNumber,1,HalfSarcomereLength,ERR,ERROR,*999)                            
                  ! old node distance
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,1,1,LeftNeighbourBioelectricNodeLocalNumber,1,DIST,ERR,ERROR,*999)
                  ! relative contraction velocity
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,1,1,LeftNeighbourBioelectricNodeLocalNumber,3,RelativeContractionVelocity, &
                    & ERR,ERROR,*999)
                ENDIF
                
                ! if we are currently at the last node on the right of the domain where the next node lies on the element boundary, store distance to right node
                IF (ABS(XI(1) + DELTA_XI1 - 1.0) < 1e-5) THEN
                  ! compute distance to right node
                  
                  XI(1) = 1.0
                  
                  ! find the interpolated position of the right neighbour bioelectric grid node from the finite elasticity FE dependent field
                  ! XI goes from 0 to 1 per FE element
                  CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,XI,INTERPOLATED_POINT,ERR,ERROR,*999)
                  Position3DRightNode(:) = INTERPOLATED_POINT%VALUES(1:3,1)
                    
                  ! compute the distance between the current node and the right node
                  DIST=SQRT( &
                    & (Position3DCurrentNode(1)-Position3DRightNode(1)) * (Position3DCurrentNode(1)-Position3DRightNode(1))+ &
                    & (Position3DCurrentNode(2)-Position3DRightNode(2)) * (Position3DCurrentNode(2)-Position3DRightNode(2))+ &
                    & (Position3DCurrentNode(3)-Position3DRightNode(3)) * (Position3DCurrentNode(3)-Position3DRightNode(3)))
    
                  ! store node distance
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U2_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,1,1,BioelectricNodeLocalNumber,7,DIST,ERR,ERROR,*999)

                  IF (DEBUGGING) PRINT "(I1,A,I4,A,F8.5)", ComputationalNodeNumber, ": node global ", &
                    & BioelectricNodeGlobalNumber, &
                    & " distance to right node d(r): ", &
                    & DIST   
                ENDIF
                  
                ! ------------- for the bioelectric element store the position of the nearest FE element Gauss point ------------------
                IF(CALC_CLOSEST_GAUSS_POINT) THEN       ! only in the first call to this function
                
                  ! get the positions of the Gauss points of the Finite Elasticity element, GAUSS_POSITIONS(components,number_of_Gauss_points)
                  GAUSS_POSITIONS=>GEOMETRIC_FIELD_ELASTICITY%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD_ELASTICITY%DECOMPOSITION% &
                    & MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(FEElementLocalNumber)%BASIS%QUADRATURE% & 
                    & QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR%GAUSS_POSITIONS
                  
                  ! calculate the closest finite elasticity Gauss point of each bioelectrics node
                  DISTANCE = 1000000.0_DP ! huge value
                  GAUSS_POINT = 0
                  IF (DEBUGGING) THEN
                    PRINT*, "        perform search for nearest gauss point in list of ",SIZE(GAUSS_POSITIONS,2),"points"
                  ENDIF
                  DO gauss_idx = 1,SIZE(GAUSS_POSITIONS,2)
                    ! compute the squared distance between the bioelectrics node and the Gauss point
                    VALUE = ( &
                      & (Xi(1)-GAUSS_POSITIONS(1,gauss_idx))*(Xi(1)-GAUSS_POSITIONS(1,gauss_idx))+ &
                      & (Xi(2)-GAUSS_POSITIONS(2,gauss_idx))*(Xi(2)-GAUSS_POSITIONS(2,gauss_idx))+ &
                      & (Xi(3)-GAUSS_POSITIONS(3,gauss_idx))*(Xi(3)-GAUSS_POSITIONS(3,gauss_idx)))
                    
                    IF (VALUE < DISTANCE) THEN
                      DISTANCE = VALUE
                      GAUSS_POINT = gauss_idx
                    ENDIF
                  ENDDO !gauss_idx
                  IF (GAUSS_POINT == 0) CALL FLAG_WARNING("Closest Gauss Point not found",ERR,ERROR,*999)
                  
                  !store the number of the nearest Gauss point
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,1,1,BioelectricNodeLocalNumber,4,GAUSS_POINT,ERR,ERROR,*999)
                ENDIF ! CALC_CLOSEST_GAUSS_POINT

                ! go to next monodomain node
                CALL BioelectricFiniteElasticity_IterateNextMonodomainNode(SOLVER_MAPPING_MONODOMAIN, GEOMETRIC_FIELD_ELASTICITY, &
                  & INDEPENDENT_FIELD_ELASTICITY, .FALSE., &
                  & GEOMETRIC_FIELD_MONODOMAIN, DEPENDENT_FIELD_MONODOMAIN, INDEPENDENT_FIELD_MONODOMAIN, &
                  & FIELD_VAR_DEP_M, FIELD_VAR_GEO_M, FIELD_VAR_IND_M_U1, FIELD_VAR_IND_M_U2, FIELD_VAR_IND_M_V, &
                  & IsFinished, FEElementGlobalNumber, FEElementLocalNumber, BioelectricNodeGlobalNumber, &
                  & BioelectricNodeLocalNumber, BioelectricNodeInFibreNumber, LeftNeighbourBioelectricNodeLocalNumber, FibreIdx, &
                  & XI, DELTA_XI1, IsFirstBioelectricNodeOfFibre, LeftNeighbourIsGhost, PreviousNodeHadNoLeftNeighbour, &
                  & ERR,ERROR,*999)
    
                ! exit loop if all the elements were traversed
                IF (IsFinished) EXIT
                
              ENDDO
              
                
              ! handle ghost elements, update geometric information
              CALL BioelectricFiniteElasticity_GhostElements(INDEPENDENT_FIELD_MONODOMAIN, GEOMETRIC_FIELD_MONODOMAIN, &
                & M_NODES_MAPPING,M_ELEMENTS_MAPPING,M_DOMAIN_TOPOLOGY_ELEMENTS,FIELD_VAR_GEO_M,FIELD_VAR_IND_M_U2, &
                & M_ELEMENTS_TOPOLOGY,ERR,ERROR,*999)
                
              !IF (CONTROL_LOOP_ELASTICITY%LOAD_INCREMENT_LOOP%ITERATION_NUMBER > 0) THEN
                !PRINT *, "Exit program in bioelectric_finite_elasticity_routines.f90:3084"
                !STOP
              !ENDIF
              
            CASE DEFAULT
              LOCAL_ERROR="The third problem specification of "// &
                & TRIM(NUMBER_TO_VSTRING(PROBLEM%specification(2),"*",ERR,ERROR))// &
                & " is not valid for a bioelectrics finite elasticity of a multi physics problem."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The second problem specification of "// &
              & TRIM(NUMBER_TO_VSTRING(PROBLEM%specification(2),"*",ERR,ERROR))// &
              & " is not valid for a multi physics problem."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Control loop problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        !the main time loop - do nothing!
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("BioelectricFiniteElasticity_UpdateGeometricField")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_UpdateGeometricField",ERR,ERROR)
    EXITS("BioelectricFiniteElasticity_UpdateGeometricField")
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_UpdateGeometricField

  !
  !================================================================================================================================
  !

  !>Interpolates the finite elasticity independent field from the biolectrics independent field. (1D->3D)
  !>NOTE: this is only temporary - will be replaced once embedded meshes are available
  SUBROUTINE BioelectricFiniteElasticity_IndependentFieldInterpolate(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the time control loop
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP_ROOT,CONTROL_LOOP_PARENT,CONTROL_LOOP_ELASTICITY,CONTROL_LOOP_MONODOMAIN
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(FIELD_TYPE), POINTER :: INDEPENDENT_FIELD_MONODOMAIN,INDEPENDENT_FIELD_ELASTICITY
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: M_ELEMENTS_MAPPING,M_NODES_MAPPING
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE_U,FIELD_VARIABLE_V,FIELD_VARIABLE_FE
    INTEGER(INTG) :: node_idx,FEElementGlobalNumber,gauss_idx,FEElementLocalNumber
    INTEGER(INTG) :: nearestGP,DofIdx
    INTEGER(INTG) :: NUMBER_OF_GAUSS_POINTS
    REAL(DP) :: ACTIVE_STRESS
    REAL(DP) :: TITIN_STRESS_UNBOUND,TITIN_STRESS_BOUND,TITIN_STRESS_CROSS_FIBRE_UNBOUND,TITIN_STRESS_CROSS_FIBRE_BOUND,ACTIVATION
    INTEGER(INTG), PARAMETER :: MAX_NUMBER_OF_GAUSS_POINTS=64
    
    REAL(DP) :: A_1,A_2,x_1,x_2
    INTEGER(INTG) :: NumberOfBioelectricNodesAtGaussPoint
    
    INTEGER(INTG), SAVE :: MaximumFEElementLocalNumber
    INTEGER(INTG), ALLOCATABLE, DIMENSION(:,:), SAVE :: NUMBER_OF_NODES
    REAL(INTG), ALLOCATABLE, DIMENSION(:,:), SAVE :: ACTIVE_STRESS_VALUES
    REAL(INTG), ALLOCATABLE, DIMENSION(:,:), SAVE :: TITIN_STRESS_VALUES_UNBOUND, TITIN_STRESS_VALUES_BOUND
    REAL(INTG), ALLOCATABLE, DIMENSION(:,:), SAVE :: TITIN_STRESS_VALUES_CROSS_FIBRE_UNBOUND
    REAL(INTG), ALLOCATABLE, DIMENSION(:,:), SAVE :: TITIN_STRESS_VALUES_CROSS_FIBRE_BOUND
    REAL(INTG), ALLOCATABLE, DIMENSION(:,:), SAVE :: ACTIVATION_VALUES
    REAL(INTG), ALLOCATABLE, DIMENSION(:,:), SAVE :: A_1_VALUES, A_2_VALUES, x_1_VALUES, x_2_VALUES
    
    !LOGICAL :: DEBUGGING = .FALSE.
    LOGICAL, PARAMETER :: DEBUGGING = .FALSE.
    
    ENTERS("BioelectricFiniteElasticity_IndependentFieldInterpolate",ERR,ERROR,*999)

#if 0
    IF (COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR) == 1) THEN
      DEBUGGING = .TRUE.
    ENDIF
#endif    
    
    NULLIFY(CONTROL_LOOP_PARENT)
    NULLIFY(CONTROL_LOOP_MONODOMAIN)
    NULLIFY(CONTROL_LOOP_ELASTICITY)
    NULLIFY(SOLVERS)
    NULLIFY(SOLVER)
    NULLIFY(FIELD_VARIABLE_U)
    NULLIFY(FIELD_VARIABLE_V)
    NULLIFY(FIELD_VARIABLE_FE)
    
    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      PROBLEM=>CONTROL_LOOP%PROBLEM
      IF(ASSOCIATED(PROBLEM)) THEN
        IF(.NOT.ALLOCATED(problem%specification)) THEN
          CALL FlagError("Problem specification is not allocated.",err,error,*999)
        ELSE IF(SIZE(problem%specification,1)<3) THEN
          CALL FlagError("Problem specification must have three entries for a bioelectric-finite elasticity problem.", &
            & err,error,*999)
        END IF
        SELECT CASE(PROBLEM%SPECIFICATION(3))
        CASE(PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE, &
          & PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE, &
          & PROBLEM_GUDUNOV_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
          & PROBLEM_STRANG_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE, &
          & PROBLEM_STRANG_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE, &
          & PROBLEM_STRANG_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE)
          IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP_PARENT,ERR,ERROR,*999)
            !--- MONODOMAIN ---
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP_PARENT,1,CONTROL_LOOP_MONODOMAIN,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP_MONODOMAIN,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
            SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
              IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR
                IF(ASSOCIATED(EQUATIONS_SET)) THEN
                  INDEPENDENT_FIELD_MONODOMAIN=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
                  IF(.NOT.ASSOCIATED(INDEPENDENT_FIELD_MONODOMAIN)) CALL FlagError("Independent field is not associated.", &
                    & ERR,ERROR,*999)
                ELSE
                  CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Solver mapping is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Solver equations is not associated.",ERR,ERROR,*999)
            ENDIF

            !--- FINITE ELASTICITY ---
            NULLIFY(SOLVERS)
            NULLIFY(SOLVER)
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP_PARENT,2,CONTROL_LOOP_ELASTICITY,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP_ELASTICITY,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
              IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR
                IF(ASSOCIATED(EQUATIONS_SET)) THEN
                  INDEPENDENT_FIELD_ELASTICITY=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
                  IF(.NOT.ASSOCIATED(INDEPENDENT_FIELD_ELASTICITY)) CALL FlagError("Independent field is not associated.",ERR, &
                    & ERROR,*999)
                ELSE
                  CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Solver mapping is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Solver equations is not associated.",ERR,ERROR,*999)
            ENDIF
                    
            !--- ALLOCATE AND INITIALIZE ---
            M_ELEMENTS_MAPPING=>INDEPENDENT_FIELD_ELASTICITY%DECOMPOSITION%DOMAIN(INDEPENDENT_FIELD_ELASTICITY%DECOMPOSITION% &
              & MESH_COMPONENT_NUMBER)%PTR%MAPPINGS%ELEMENTS
            M_NODES_MAPPING=>INDEPENDENT_FIELD_MONODOMAIN%DECOMPOSITION%DOMAIN(INDEPENDENT_FIELD_MONODOMAIN%DECOMPOSITION% &
              & MESH_COMPONENT_NUMBER)%PTR%MAPPINGS%NODES

            CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE_U,ERR,ERROR,*999)
            CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE,FIELD_VARIABLE_V,ERR,ERROR,*999)
            CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD_ELASTICITY,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE_FE,ERR,ERROR,*999)

            ! allocate fields at gauss points that are only needed for multiple subtypes
            IF (.NOT. ALLOCATED(NUMBER_OF_NODES)) THEN
            
              ! check if number of gauss points is below 64 in every element, also determine the maximum local number of a FE element
              MaximumFEElementLocalNumber = 0
              
              IF (DEBUGGING) PRINT *, "FEElementGlobalNumber,FEElementLocalNumber,MaximumFEElementLocalNumber"
              
              ! loop over the finite elasticity elements
              DO FEElementGlobalNumber=M_ELEMENTS_MAPPING%INTERNAL_START,M_ELEMENTS_MAPPING%BOUNDARY_FINISH
                FEElementLocalNumber=M_ELEMENTS_MAPPING%DOMAIN_LIST(FEElementGlobalNumber)
                
                NUMBER_OF_GAUSS_POINTS=INDEPENDENT_FIELD_ELASTICITY%DECOMPOSITION%DOMAIN(INDEPENDENT_FIELD_ELASTICITY% &
                  & DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(FEElementLocalNumber)%BASIS%QUADRATURE% & 
                  & QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR%NUMBER_OF_GAUSS

                IF(NUMBER_OF_GAUSS_POINTS>MAX_NUMBER_OF_GAUSS_POINTS) CALL FlagError( & 
                  & "NUMBER_OF_GAUSS_POINTS is greater than MAX_NUMBER_OF_GAUSS_POINTS.",ERR,ERROR,*999)
                MaximumFEElementLocalNumber = MAX(MaximumFEElementLocalNumber, FEElementLocalNumber)
                  
                IF (DEBUGGING) PRINT *, FEElementGlobalNumber,FEElementLocalNumber,MaximumFEElementLocalNumber
                  
              ENDDO
              
              ! allocate variable for number of bioelectric nodes for which the gauss point is the nearest
              ALLOCATE(NUMBER_OF_NODES(MaximumFEElementLocalNumber, MAX_NUMBER_OF_GAUSS_POINTS), Stat=Err)
              IF(Err /= 0) CALL FlagError("Could not allocate NUMBER_OF_NODES.",err,error,*999)
              
            ENDIF
            
            SELECT CASE(PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE, &
             & PROBLEM_STRANG_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE, &
             & PROBLEM_GUDUNOV_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
             & PROBLEM_STRANG_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE)
              IF (.NOT. ALLOCATED(ACTIVE_STRESS_VALUES)) THEN
                ALLOCATE(ACTIVE_STRESS_VALUES(MaximumFEElementLocalNumber, MAX_NUMBER_OF_GAUSS_POINTS), Stat=Err)
                IF(Err /= 0) CALL FlagError("Could not allocate ACTIVE_STRESS_VALUES.",err,error,*999)
                
                ACTIVE_STRESS_VALUES = 0.0_DP
              ENDIF
            END SELECT
              
            ! allocate fields at gauss points that are only needed for a single subtype
            SELECT CASE(PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_GUDUNOV_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
              & PROBLEM_STRANG_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE)
              
              IF (.NOT. ALLOCATED(TITIN_STRESS_VALUES_UNBOUND)) THEN
                ALLOCATE(TITIN_STRESS_VALUES_UNBOUND(MaximumFEElementLocalNumber, MAX_NUMBER_OF_GAUSS_POINTS), Stat=Err)
                IF(Err /= 0) CALL FlagError("Could not allocate TITIN_STRESS_VALUES_UNBOUND.",err,error,*999)
                
              ENDIF
              IF (.NOT. ALLOCATED(TITIN_STRESS_VALUES_BOUND)) THEN
                ALLOCATE(TITIN_STRESS_VALUES_BOUND(MaximumFEElementLocalNumber, MAX_NUMBER_OF_GAUSS_POINTS), Stat=Err)
                IF(Err /= 0) CALL FlagError("Could not allocate TITIN_STRESS_VALUES_BOUND.",err,error,*999)
                
              ENDIF
              IF (.NOT. ALLOCATED(TITIN_STRESS_VALUES_CROSS_FIBRE_UNBOUND)) THEN
                ALLOCATE(TITIN_STRESS_VALUES_CROSS_FIBRE_UNBOUND(MaximumFEElementLocalNumber, MAX_NUMBER_OF_GAUSS_POINTS), Stat=Err)
                IF(Err /= 0) CALL FlagError("Could not allocate TITIN_STRESS_VALUES_CROSS_FIBRE_UNBOUND.",err,error,*999)
                
              ENDIF
              IF (.NOT. ALLOCATED(TITIN_STRESS_VALUES_CROSS_FIBRE_BOUND)) THEN
                ALLOCATE(TITIN_STRESS_VALUES_CROSS_FIBRE_BOUND(MaximumFEElementLocalNumber, MAX_NUMBER_OF_GAUSS_POINTS), Stat=Err)
                IF(Err /= 0) CALL FlagError("Could not allocate TITIN_STRESS_VALUES_CROSS_FIBRE_BOUND.",err,error,*999)
                
              ENDIF
              
              IF (.NOT. ALLOCATED(ACTIVATION_VALUES)) THEN
                ALLOCATE(ACTIVATION_VALUES(MaximumFEElementLocalNumber, MAX_NUMBER_OF_GAUSS_POINTS), Stat=Err)
                IF(Err /= 0) CALL FlagError("Could not allocate ACTIVATION_VALUES.",err,error,*999)
                
              ENDIF
              
            CASE(PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE, PROBLEM_STRANG_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE)
              IF (.NOT. ALLOCATED(A_1_VALUES)) THEN
                ALLOCATE(A_1_VALUES(MaximumFEElementLocalNumber, MAX_NUMBER_OF_GAUSS_POINTS), Stat=Err)
                IF(Err /= 0) CALL FlagError("Could not allocate A_1_VALUES.",err,error,*999)
                
              ENDIF
              
              IF (.NOT. ALLOCATED(A_2_VALUES)) THEN
                ALLOCATE(A_2_VALUES(MaximumFEElementLocalNumber, MAX_NUMBER_OF_GAUSS_POINTS), Stat=Err)
                IF(Err /= 0) CALL FlagError("Could not allocate A_2_VALUES.",err,error,*999)
                
              ENDIF
              
              
              IF (.NOT. ALLOCATED(x_1_VALUES)) THEN
                ALLOCATE(x_1_VALUES(MaximumFEElementLocalNumber, MAX_NUMBER_OF_GAUSS_POINTS), Stat=Err)
                IF(Err /= 0) CALL FlagError("Could not allocate x_1_VALUES.",err,error,*999)
                
              ENDIF
              
              IF (.NOT. ALLOCATED(x_2_VALUES)) THEN
                ALLOCATE(x_2_VALUES(MaximumFEElementLocalNumber, MAX_NUMBER_OF_GAUSS_POINTS), Stat=Err)
                IF(Err /= 0) CALL FlagError("Could not allocate x_2_VALUES.",err,error,*999)
                
              ENDIF
              
            END SELECT
            
            ! initialize values to 0
            SELECT CASE(PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE,PROBLEM_STRANG_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE)

              NUMBER_OF_NODES = 0
              ACTIVE_STRESS_VALUES = 0.0_DP
              
            CASE(PROBLEM_GUDUNOV_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE,PROBLEM_STRANG_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE)
            
              NUMBER_OF_NODES = 0
              ACTIVE_STRESS_VALUES = 0.0_DP
              
              TITIN_STRESS_VALUES_UNBOUND = 0.0_DP
              TITIN_STRESS_VALUES_BOUND = 0.0_DP
              TITIN_STRESS_VALUES_CROSS_FIBRE_UNBOUND = 0.0_DP
              TITIN_STRESS_VALUES_CROSS_FIBRE_BOUND = 0.0_DP
              ACTIVATION_VALUES = 0.0_DP
              
            CASE(PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE,PROBLEM_STRANG_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE)

              NUMBER_OF_NODES = 0
              A_1_VALUES = 0.0_DP
              A_2_VALUES = 0.0_DP
              x_1_VALUES = 0.0_DP
              x_2_VALUES = 0.0_DP
              
            END SELECT
            
            IF (DEBUGGING) THEN
              PRINT *, "FE elements associated with bioelectric nodes:"
              DO node_idx=1,M_NODES_MAPPING%NUMBER_OF_LOCAL
                  
                ! get the local FE element local number
                DofIdx=FIELD_VARIABLE_V%COMPONENTS(5)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                  & VERSIONS(1)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & DofIdx,FEElementLocalNumber,ERR,ERROR,*999) !component 5 of variable V contains inElem info (LOCAL NUMBERING!!!)
                  
                PRINT *, "node ", node_idx,", FE element local", FEElementLocalNumber
                  
              ENDDO
              
            ENDIF
              
            IF (DEBUGGING) PRINT*, "Interpolate"
            !--- NOW INTERPOLATE ---
            ! ------------ accumulate stress values at gauss points ----------------------              
            !loop over the bioelectrics nodes
            DO node_idx=1,M_NODES_MAPPING%NUMBER_OF_LOCAL
              
              ! get the local FE element local number
              DofIdx=FIELD_VARIABLE_V%COMPONENTS(5)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                & VERSIONS(1)
              CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & DofIdx,FEElementLocalNumber,ERR,ERROR,*999) !component 5 of variable V contains inElem info (LOCAL NUMBERING!!!)
                
              ! do not consider nodes for which the containing FE element is on a different processors subdomain
              IF (FEElementLocalNumber > MaximumFEElementLocalNumber) CYCLE
                
              ! get the number of the nearest Gauss Point
              !component 4 of variable V contains Nearest Gauss Point info
              DofIdx=FIELD_VARIABLE_V%COMPONENTS(4)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                & VERSIONS(1)
              CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & DofIdx,nearestGP,ERR,ERROR,*999)
              IF(nearestGP>MAX_NUMBER_OF_GAUSS_POINTS) CALL FlagError( &
                & "Nearest Gauss Point is greater than MAX_NUMBER_OF_GAUSS_POINTS.",ERR,ERROR,*999)

              IF (DEBUGGING) THEN
                WRITE(*,"(3(A, I3))",advance='no') "node ", node_idx, ", FE element: ", FEElementLocalNumber, & 
                  & ", nearestGP:", nearestGP
                WRITE(*,"(1(A,I3))",advance='no') ", NUMBER_OF_NODES:", NUMBER_OF_NODES(FEElementLocalNumber, nearestGP)
                PRINT *, ""
              ENDIF
              
              ! depending on problem subtype accumulate stress values
              SELECT CASE(PROBLEM%SPECIFICATION(3))
              CASE(PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE, PROBLEM_STRANG_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE)

                !component 1 of variable U contains the active stress
                DofIdx=FIELD_VARIABLE_U%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                  & VERSIONS(1)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,DofIdx,ACTIVE_STRESS,ERR,ERROR,*999)
                
                !count the number of bioelectrics nodes that are closest to each finite elasticity Gauss point
                NUMBER_OF_NODES(FEElementLocalNumber, nearestGP)=NUMBER_OF_NODES(FEElementLocalNumber, nearestGP)+1
                !add up the active stress value
                ACTIVE_STRESS_VALUES(FEElementLocalNumber, nearestGP) = &
                  & ACTIVE_STRESS_VALUES(FEElementLocalNumber, nearestGP)+ACTIVE_STRESS
                
              CASE(PROBLEM_GUDUNOV_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, PROBLEM_STRANG_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE)

                !component 1 of variable U contains the active stress
                DofIdx=FIELD_VARIABLE_U%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                  & VERSIONS(1)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,DofIdx,ACTIVE_STRESS,ERR,ERROR,*999)
                
                !count the number of bioelectrics nodes that are closest to each finite elasticity Gauss point
                NUMBER_OF_NODES(FEElementLocalNumber, nearestGP)=NUMBER_OF_NODES(FEElementLocalNumber, nearestGP)+1
                !add up the active stress value
                ACTIVE_STRESS_VALUES(FEElementLocalNumber, nearestGP) = &
                  & ACTIVE_STRESS_VALUES(FEElementLocalNumber, nearestGP)+ACTIVE_STRESS
                
                ! get stress values from bioelectric nodes
                !component 2 of variable U contains the titin stress unbound
                DofIdx=FIELD_VARIABLE_U%COMPONENTS(2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                  & VERSIONS(1)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,DofIdx,TITIN_STRESS_UNBOUND,ERR,ERROR,*999)
                !component 3 of variable U contains the titin stress bound
                DofIdx=FIELD_VARIABLE_U%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                  & VERSIONS(1)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,DofIdx,TITIN_STRESS_BOUND,ERR,ERROR,*999)
                !component 4 of variable U contains the titin XF-stress (cross-fibre directions) unbound
                DofIdx=FIELD_VARIABLE_U%COMPONENTS(4)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                  & VERSIONS(1)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,DofIdx,TITIN_STRESS_CROSS_FIBRE_UNBOUND,ERR,ERROR,*999)
                !component 5 of variable U contains the titin XF-stress (cross-fibre directions) bound
                DofIdx=FIELD_VARIABLE_U%COMPONENTS(5)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                  & VERSIONS(1)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,DofIdx,TITIN_STRESS_CROSS_FIBRE_BOUND,ERR,ERROR,*999)
                !component 6 of variable U contains the titin activation
                DofIdx=FIELD_VARIABLE_U%COMPONENTS(6)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                  & VERSIONS(1)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,DofIdx,ACTIVATION,ERR,ERROR,*999)

                TITIN_STRESS_VALUES_UNBOUND(FEElementLocalNumber, nearestGP) = &
                  & TITIN_STRESS_VALUES_UNBOUND(FEElementLocalNumber, nearestGP)+TITIN_STRESS_UNBOUND
                TITIN_STRESS_VALUES_BOUND(FEElementLocalNumber, nearestGP) = &
                  & TITIN_STRESS_VALUES_BOUND(FEElementLocalNumber, nearestGP)+TITIN_STRESS_BOUND
                TITIN_STRESS_VALUES_CROSS_FIBRE_UNBOUND(FEElementLocalNumber, nearestGP) = &
                  & TITIN_STRESS_VALUES_CROSS_FIBRE_UNBOUND(FEElementLocalNumber, nearestGP) + TITIN_STRESS_CROSS_FIBRE_UNBOUND
                TITIN_STRESS_VALUES_CROSS_FIBRE_BOUND(FEElementLocalNumber, nearestGP) = &
                  & TITIN_STRESS_VALUES_CROSS_FIBRE_BOUND(FEElementLocalNumber, nearestGP) + TITIN_STRESS_CROSS_FIBRE_BOUND
                ACTIVATION_VALUES(FEElementLocalNumber, nearestGP) = ACTIVATION_VALUES(FEElementLocalNumber, nearestGP)+ACTIVATION
                
              CASE(PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE, PROBLEM_STRANG_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE)

                !count the number of bioelectrics nodes that are closest to each finite elasticity Gauss point
                NUMBER_OF_NODES(FEElementLocalNumber, nearestGP)=NUMBER_OF_NODES(FEElementLocalNumber, nearestGP)+1

                !component 1 of variable U contains A_1
                DofIdx=FIELD_VARIABLE_U%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                  & VERSIONS(1)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,DofIdx,A_1,ERR,ERROR,*999)
                !component 2 of variable U contains A_2
                DofIdx=FIELD_VARIABLE_U%COMPONENTS(2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                  & VERSIONS(1)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,DofIdx,A_2,ERR,ERROR,*999)
                !component 3 of variable U contains x_1
                DofIdx=FIELD_VARIABLE_U%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                  & VERSIONS(1)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,DofIdx,x_1,ERR,ERROR,*999)
                !component 4 of variable U contains x_2
                DofIdx=FIELD_VARIABLE_U%COMPONENTS(4)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                  & VERSIONS(1)
                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,DofIdx,x_2,ERR,ERROR,*999)

                A_1_VALUES(FEElementLocalNumber, nearestGP) = A_1_VALUES(FEElementLocalNumber, nearestGP) + A_1
                A_2_VALUES(FEElementLocalNumber, nearestGP) = A_2_VALUES(FEElementLocalNumber, nearestGP) + A_2
                x_1_VALUES(FEElementLocalNumber, nearestGP) = x_1_VALUES(FEElementLocalNumber, nearestGP) + x_1
                x_2_VALUES(FEElementLocalNumber, nearestGP) = x_2_VALUES(FEElementLocalNumber, nearestGP) + x_2
              END SELECT
            ENDDO
          
            IF (DEBUGGING) PRINT *, "compute averages"
          
            ! ------------ compute average values at gauss points and store them to FE elements --------
            ! loop over internal and boundary finite elasticity elements
            DO FEElementGlobalNumber=M_ELEMENTS_MAPPING%INTERNAL_START,M_ELEMENTS_MAPPING%BOUNDARY_FINISH
              FEElementLocalNumber=M_ELEMENTS_MAPPING%DOMAIN_LIST(FEElementGlobalNumber)
              
              NUMBER_OF_GAUSS_POINTS=INDEPENDENT_FIELD_ELASTICITY%DECOMPOSITION%DOMAIN(INDEPENDENT_FIELD_ELASTICITY% &
                & DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(FEElementLocalNumber)%BASIS%QUADRATURE% & 
                & QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR%NUMBER_OF_GAUSS

              IF (DEBUGGING) PRINT *, "FE element local ", FEElementLocalNumber,", NUMBER_OF_GAUSS_POINTS:", NUMBER_OF_GAUSS_POINTS
                
              !loop over the finite elasticity Gauss points
              DO gauss_idx=1,NUMBER_OF_GAUSS_POINTS
              
                NumberOfBioelectricNodesAtGaussPoint = NUMBER_OF_NODES(FEElementLocalNumber, gauss_idx)
                IF (DEBUGGING) PRINT *, "FE element local",FEElementLocalNumber,", Gauss Point",gauss_idx, &
                  & ", number of nodes: ",NumberOfBioelectricNodesAtGaussPoint
              
                SELECT CASE(PROBLEM%SPECIFICATION(3))
                CASE(PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE,PROBLEM_STRANG_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE)

                  ! average values at gauss point
                  IF(NumberOfBioelectricNodesAtGaussPoint == 0) THEN
                    ACTIVE_STRESS = 0.0_DP
                  ELSE
                    ACTIVE_STRESS = ACTIVE_STRESS_VALUES(FEElementLocalNumber, gauss_idx) / NumberOfBioelectricNodesAtGaussPoint
                  ENDIF

                  ! store values in independent field
                  DofIdx=FIELD_VARIABLE_FE%COMPONENTS(1)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx, & 
                   & FEElementLocalNumber)
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,DofIdx,ACTIVE_STRESS,ERR,ERROR,*999)

                CASE(PROBLEM_GUDUNOV_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, PROBLEM_STRANG_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE)
                  ! average values at gauss point
                  IF(NumberOfBioelectricNodesAtGaussPoint == 0) THEN
                    
                    ACTIVE_STRESS=0.0_DP
                    ACTIVATION=0.0_DP
                    TITIN_STRESS_UNBOUND=0.0_DP
                    TITIN_STRESS_BOUND=0.0_DP
                    TITIN_STRESS_CROSS_FIBRE_UNBOUND=0.0_DP
                    TITIN_STRESS_CROSS_FIBRE_BOUND=0.0_DP

                  ELSE
                    ACTIVE_STRESS = ACTIVE_STRESS_VALUES(FEElementLocalNumber, gauss_idx) / NumberOfBioelectricNodesAtGaussPoint
                    ACTIVATION = ACTIVATION_VALUES(FEElementLocalNumber, gauss_idx) / NumberOfBioelectricNodesAtGaussPoint
                    TITIN_STRESS_UNBOUND = TITIN_STRESS_VALUES_UNBOUND(FEElementLocalNumber, gauss_idx) / &
                      & NumberOfBioelectricNodesAtGaussPoint
                    TITIN_STRESS_BOUND = TITIN_STRESS_VALUES_BOUND(FEElementLocalNumber, gauss_idx) / &
                      & NumberOfBioelectricNodesAtGaussPoint
                    TITIN_STRESS_CROSS_FIBRE_UNBOUND = TITIN_STRESS_VALUES_CROSS_FIBRE_UNBOUND(FEElementLocalNumber, gauss_idx) / &
                      & NumberOfBioelectricNodesAtGaussPoint
                    TITIN_STRESS_CROSS_FIBRE_BOUND = TITIN_STRESS_VALUES_CROSS_FIBRE_BOUND(FEElementLocalNumber, gauss_idx) / &
                      & NumberOfBioelectricNodesAtGaussPoint
                  ENDIF
                  
                  ! store values in independent field
                  DofIdx=FIELD_VARIABLE_FE%COMPONENTS(1)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx, & 
                    & FEElementLocalNumber)
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,DofIdx,ACTIVE_STRESS,ERR,ERROR,*999)
                  DofIdx=FIELD_VARIABLE_FE%COMPONENTS(2)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx, & 
                    & FEElementLocalNumber)
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,DofIdx,TITIN_STRESS_UNBOUND,ERR,ERROR,*999)
                  DofIdx=FIELD_VARIABLE_FE%COMPONENTS(3)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx, & 
                    & FEElementLocalNumber)
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,DofIdx,TITIN_STRESS_BOUND,ERR,ERROR,*999)
                  DofIdx=FIELD_VARIABLE_FE%COMPONENTS(4)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx, & 
                    & FEElementLocalNumber)
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,DofIdx,TITIN_STRESS_CROSS_FIBRE_UNBOUND,ERR,ERROR,*999)
                  DofIdx=FIELD_VARIABLE_FE%COMPONENTS(5)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx, & 
                    & FEElementLocalNumber)
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,DofIdx,TITIN_STRESS_CROSS_FIBRE_BOUND,ERR,ERROR,*999)
                  DofIdx=FIELD_VARIABLE_FE%COMPONENTS(6)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx, & 
                    & FEElementLocalNumber)
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,DofIdx,ACTIVATION,ERR,ERROR,*999)

                CASE(PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE,PROBLEM_STRANG_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE)
                  ! average values at gauss points
                  IF(NumberOfBioelectricNodesAtGaussPoint == 0) THEN
                    A_1=0.0_DP
                    A_2=0.0_DP
                    x_1=0.0_DP
                    x_2=0.0_DP
                  ELSE
                    A_1 = A_1_VALUES(FEElementLocalNumber, gauss_idx) / NumberOfBioelectricNodesAtGaussPoint
                    A_2 = A_2_VALUES(FEElementLocalNumber, gauss_idx) / NumberOfBioelectricNodesAtGaussPoint
                    x_1 = x_1_VALUES(FEElementLocalNumber, gauss_idx) / NumberOfBioelectricNodesAtGaussPoint
                    x_2 = x_2_VALUES(FEElementLocalNumber, gauss_idx) / NumberOfBioelectricNodesAtGaussPoint
                  
                  ENDIF
                  
                  ! store values in independent field
                  DofIdx=FIELD_VARIABLE_FE%COMPONENTS(1)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx, & 
                    & FEElementLocalNumber)
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,DofIdx,A_1,ERR,ERROR,*999)
                  DofIdx=FIELD_VARIABLE_FE%COMPONENTS(2)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx, & 
                    & FEElementLocalNumber)
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,DofIdx,A_2,ERR,ERROR,*999)
                  DofIdx=FIELD_VARIABLE_FE%COMPONENTS(3)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx, & 
                    & FEElementLocalNumber)
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,DofIdx,x_1,ERR,ERROR,*999)
                  DofIdx=FIELD_VARIABLE_FE%COMPONENTS(4)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gauss_idx, & 
                    & FEElementLocalNumber)
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,DofIdx,x_2,ERR,ERROR,*999)

                END SELECT

              ENDDO !gauss_idx
            ENDDO !FEElementGlobalNumber

            !now the ghost elements -- get the relevant info from the other computational nodes
            CALL FIELD_PARAMETER_SET_UPDATE_START(INDEPENDENT_FIELD_ELASTICITY, & 
              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
            CALL FIELD_PARAMETER_SET_UPDATE_FINISH(INDEPENDENT_FIELD_ELASTICITY, & 
              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)

          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="Independent field interpolation is not implemented for problem subtype " &
            & //TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("BioelectricFiniteElasticity_IndependentFieldInterpolate")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_IndependentFieldInterpolate",ERR,ERROR)
    EXITS("BioelectricFiniteElasticity_IndependentFieldInterpolate")
    RETURN 1

  END SUBROUTINE BioelectricFiniteElasticity_IndependentFieldInterpolate

  !
  !================================================================================================================================
  !

  !>Computes force enhancement based on the titin model of C Rode et al. (2009). Force depression is not yet implemented.
  !>Titin-induced force enhancement and force depression: A 'sticky-spring' mechanism in muscle contractions? C Rode, T Siebert, R Blickhan - Journal of theoretical biology, 2009.
  SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_COMPUTE_TITIN(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !A pointer to the time control loop
    INTEGER(INTG), INTENT(OUT) :: ERR !The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP_ROOT,CONTROL_LOOP_PARENT,CONTROL_LOOP_MONODOMAIN
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(FIELD_TYPE), POINTER :: INDEPENDENT_FIELD_MONODOMAIN
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VAR_IND_M_U1
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: M_NODES_MAPPING
    INTEGER(INTG) :: node_idx,DofIdx
    REAL(DP), PARAMETER :: LENGTH_ACTIN=1.04_DP,LENGTH_MBAND=0.0625_DP,LENGTH_MYOSIN=0.7375_DP
    REAL(DP), PARAMETER :: LENGTH_ZERO=0.635_DP,LENGTH_ZDISC=0.05_DP
    REAL(DP) :: F0,FORCE,TITIN_BOUND,TITIN_XF_BOUND,TITIN_UNBOUND,TITIN_XF_UNBOUND
    REAL(DP) :: ELONGATION,ELONGATION_NEW
    REAL(DP) :: ELONGATION_DIST_IG,ELONGATION_PEVK
    REAL(DP) :: LENGTH_TITIN,LENGTH_INIT_TITIN,LENGTH_DIST_IG_F0,LENGTH_DIST_IG
    REAL(DP) :: STIFFNESS_PEVK
    REAL(DP) :: SARCO_LENGTH_AT_ACTIVATION,SARCO_LENGTH
    REAL(DP) :: ACTIN_MYOSIN_DISTANCE,d10
    REAL(DP) :: SLOPE,STIFFNESS_DIST,FORCE_DISTAL_IG,DELTA_F,DIFF_QUOT
    REAL(DP), PARAMETER, DIMENSION(5) :: COEFF_MATRIX=[5.0239_DP,-0.6717_DP,-2.5841_DP,-5.0128_DP,-5.0239_DP]
    REAL(DP), DIMENSION(250) :: LENGTHS_DIST_IG,FORCES_DIST_IG
    REAL(DP), PARAMETER :: DX=0.001_DP
    REAL(DP), PARAMETER :: FORCE_INCREMENT=1.e-5_DP,TOL=1.e-5_DP
    INTEGER(INTG) :: INDEX_REF,INDEX_PSEUDO,INDEX_I
    INTEGER(INTG), PARAMETER :: DIM_DATA=250
    INTEGER(INTG) :: SWITCH_MODEL
    
    NULLIFY(CONTROL_LOOP_ROOT)
    NULLIFY(CONTROL_LOOP_PARENT)
    NULLIFY(CONTROL_LOOP_MONODOMAIN)
    NULLIFY(PROBLEM)
    NULLIFY(SOLVERS)
    NULLIFY(SOLVER)
    NULLIFY(INDEPENDENT_FIELD_MONODOMAIN)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVER_MAPPING)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(FIELD_VAR_IND_M_U1)
    
    ENTERS("BIOELECTRIC_FINITE_ELASTICITY_COMPUTE_TITIN",ERR,ERROR,*999)

    !the relation between length and force of the distal_Ig region is very nonlinear. Linear interpolation of Rode's data is used here.
    FORCES_DIST_IG= &
      & [0.0_DP,0.001_DP,0.002_DP,0.003_DP,0.004_DP,0.005_DP,0.006_DP,0.007_DP,0.008_DP,0.009_DP,0.01_DP,0.011_DP,0.012_DP, &
      & 0.013_DP,0.014_DP,0.015_DP,0.016_DP,0.017_DP,0.018_DP,0.019_DP,0.02_DP,0.021_DP,0.022_DP,0.023_DP,0.024_DP, &
      & 0.025_DP,0.026_DP,0.027_DP,0.028_DP,0.029_DP,0.03_DP,0.031_DP,0.032_DP,0.033_DP,0.034_DP,0.035_DP,0.036_DP, &
      & 0.037_DP,0.038_DP,0.039_DP,0.04_DP,0.041_DP,0.042_DP,0.043_DP,0.044_DP,0.045_DP,0.046_DP,0.047_DP,0.048_DP, &
      & 0.049_DP,0.05_DP,0.051_DP,0.052_DP,0.053_DP,0.054_DP,0.055_DP,0.056_DP,0.057_DP,0.058_DP,0.059_DP,0.06_DP, &
      & 0.061_DP,0.062_DP,0.063_DP,0.064_DP,0.065_DP,0.066_DP,0.067_DP,0.068_DP,0.069_DP,0.070_DP,0.071_DP,0.072_DP, &
      & 0.073_DP,0.074_DP,0.075_DP,0.076_DP,0.077_DP,0.078_DP,0.079_DP,0.080_DP,0.081_DP,0.082_DP,0.083_DP,0.084_DP, &
      & 0.085_DP,0.086_DP,0.087_DP,0.088_DP,0.089_DP,0.090_DP,0.091_DP,0.092_DP,0.093_DP,0.094_DP,0.095_DP,0.096_DP, &
      & 0.097_DP,0.098_DP,0.099_DP,0.1_DP,0.101_DP,0.102_DP,0.103_DP,0.104_DP,0.105_DP,0.106_DP,0.107_DP,0.108_DP, &
      & 0.109_DP,0.11_DP,0.111_DP,0.112_DP,0.113_DP,0.114_DP,0.115_DP,0.116_DP,0.117_DP,0.118_DP,0.119_DP,0.12_DP, &
      & 0.121_DP,0.122_DP,0.123_DP,0.124_DP,0.125_DP,0.126_DP,0.127_DP,0.128_DP,0.129_DP,0.13_DP,0.131_DP,0.132_DP, &
      & 0.133_DP,0.134_DP,0.135_DP,0.136_DP,0.137_DP,0.138_DP,0.139_DP,0.14_DP,0.141_DP,0.142_DP,0.143_DP,0.144_DP, &
      & 0.145_DP,0.146_DP,0.147_DP,0.148_DP,0.149_DP,0.15_DP,0.151_DP,0.152_DP,0.153_DP,0.154_DP,0.155_DP,0.156_DP, &
      & 0.157_DP,0.158_DP,0.159_DP,0.16_DP,0.161_DP,0.162_DP,0.163_DP,0.164_DP,0.165_DP,0.166_DP,0.167_DP, 0.168_DP, &
      & 0.169_DP,0.17_DP,0.171_DP,0.172_DP,0.173_DP,0.174_DP,0.175_DP,0.176_DP,0.177_DP,0.178_DP,0.179_DP,0.18_DP, &
      & 0.181_DP,0.182_DP,0.183_DP,0.184_DP,0.185_DP,0.186_DP,0.187_DP,0.188_DP,0.189_DP,0.19_DP,0.191_DP,0.192_DP, &
      & 0.193_DP,0.194_DP,0.195_DP,0.196_DP,0.197_DP,0.198_DP,0.199_DP,0.2_DP,0.201_DP,0.202_DP,0.203_DP,0.204_DP, &
      & 0.205_DP,0.206_DP,0.207_DP,0.208_DP,0.209_DP,0.21_DP,0.211_DP,0.212_DP,0.213_DP,0.214_DP,0.215_DP,0.216_DP, &
      & 0.217_DP,0.218_DP,0.219_DP,0.22_DP,0.221_DP,0.222_DP,0.223_DP,0.224_DP,0.225_DP,0.226_DP,0.227_DP,0.228_DP, &
      & 0.229_DP,0.23_DP,0.231_DP,0.232_DP,0.233_DP,0.234_DP,0.235_DP,0.236_DP,0.237_DP,0.238_DP,0.239_DP,0.24_DP, &
      & 0.241_DP,0.242_DP,0.243_DP,0.244_DP,0.245_DP,0.246_DP,0.247_DP,0.248_DP,0.249_DP]
    LENGTHS_DIST_IG= &
      & [0.0_DP,0.03461753545561_DP,0.049729169766010_DP,0.058506390323323_DP,0.064606296848594_DP,0.06922519775133_DP, &
      & 0.0729080120998386_DP,0.0759458896446241_DP,0.0785230355395668_DP,0.0807314335143191_DP,0.0826674161660979_DP, &
      & 0.0843819721302_DP,0.0859161360822_DP,0.087300738288_DP,0.0885510536196_DP,0.08970061165_DP,0.090751366_DP, &
      & 0.0917285714001_DP,0.0926271710799_DP,0.093467026018_DP,0.094254010845_DP,0.094992014919_DP,0.09568255451_DP, &
      & 0.0963346932312_DP,0.0969518155718_DP,0.097537000419_DP,0.098093047899_DP,0.098622503780_DP,0.09912768169_DP, &
      & 0.0996103583026_DP,0.1000734170008_DP,0.100517614158_DP,0.100942907967_DP,0.101351270601_DP,0.101745244913_DP, &
      & 0.1021260518375_DP,0.1024947934706_DP,0.102851281365_DP,0.103194317381_DP,0.103528086731_DP,0.103853341524_DP, &
      & 0.1041686065999_DP,0.1044739635997_DP,0.104772842609_DP,0.105064758806_DP,0.105347078318_DP,0.105624524436_DP, &
      & 0.1058959740402_DP,0.1061596374590_DP,0.106419674124_DP,0.106673222977_DP,0.106921779223_DP,0.107167251003_DP, &
      & 0.1074055929211_DP,0.1076418938850_DP,0.107872798106_DP,0.108100580266_DP,0.108324935479_DP,0.108545154704_DP, &
      & 0.1087634145225_DP,0.1089769308376_DP,0.109189523632_DP,0.109397109133_DP,0.109604335645_DP,0.109806785604_DP, &
      & 0.1100090293896_DP,0.1102069598668_DP,0.110404790711_DP,0.110598542577_DP,0.110792294444_DP,0.110982362315_DP, &
      & 0.1111722575399_DP,0.1113591719379_DP,0.111545514634_DP,0.111729654458_DP,0.111912731875_DP,0.112094428489_DP, &
      & 0.1122745119339_DP,0.1124540532647_DP,0.112631398979_DP,0.112808744694_DP,0.112983883284_DP,0.113158733268_DP, &
      & 0.1133324054841_DP,0.1135049882670_DP,0.113677360515_DP,0.113847891889_DP,0.114018423263_DP,0.114187784974_DP, &
      & 0.1143564686814_DP,0.1145249703398_DP,0.114691998719_DP,0.114859027099_DP,0.115025270473_DP,0.115190825075_DP, &
      & 0.1153563796784_DP,0.1155207617063_DP,0.115685013868_DP,0.115849024045_DP,0.116012135436_DP,0.116175246828_DP, &
      & 0.1163378963393_DP,0.1165000194746_DP,0.116662142610_DP,0.116823704015_DP,0.116984982740_DP,0.117146261465_DP, &
      & 0.1173069696799_DP,0.1174675396300_DP,0.117628109580_DP,0.117788165045_DP,0.117948154078_DP,0.118108143111_DP, &
      & 0.1182677147299_DP,0.1184272433357_DP,0.118586771941_DP,0.118745999791_DP,0.118905181478_DP,0.119064363164_DP, &
      & 0.1192233610196_DP,0.1193823026790_DP,0.119541244338_DP,0.119700101992_DP,0.119858904246_DP,0.120017706500_DP, &
      & 0.1201764919257_DP,0.1203352494533_DP,0.120494006981_DP,0.120652768322_DP,0.120811570168_DP,0.120970372014_DP, &
      & 0.1211291738603_DP,0.1212880693042_DP,0.121446999172_DP,0.121605929041_DP,0.121764923098_DP,0.121924059629_DP, &
      & 0.1220831961613_DP,0.1222423326927_DP,0.122401700265_DP,0.122561117298_DP,0.122720534331_DP,0.122880045624_DP, &
      & 0.1230398124447_DP,0.1231995792652_DP,0.123359346085_DP,0.123519381317_DP,0.123679562894_DP,0.123839744470_DP, &
      & 0.1239999260467_DP,0.1241605622478_DP,0.124321219453_DP,0.124481876659_DP,0.124642637543_DP,0.124803827368_DP, &
      & 0.1249650171945_DP,0.1251262070202_DP,0.125287609380_DP,0.125449385133_DP,0.125611160886_DP,0.125772936638_DP, &
      & 0.1259350043157_DP,0.1260974158090_DP,0.126259827302_DP,0.126422238795_DP,0.126584979901_DP,0.126748073636_DP, &
      & 0.1269111673704_DP,0.1270742611047_DP,0.127237669513_DP,0.127401488846_DP,0.127565308178_DP,0.127729127511_DP, &
      & 0.1278931842348_DP,0.1280577695422_DP,0.128222354849_DP,0.128386940157_DP,0.128551614623_DP,0.128717003455_DP, &
      & 0.1288823922862_DP,0.1290477811173_DP,0.129213169948_DP,0.129379259581_DP,0.129545486803_DP,0.129711714025_DP, &
      & 0.1298779412472_DP,0.1300445897140_DP,0.130211687650_DP,0.130378785586_DP,0.130545883522_DP,0.130713029866_DP, &
      & 0.1308810284274_DP,0.1310490269884_DP,0.131217025549_DP,0.131385024110_DP,0.131553528414_DP,0.131722455223_DP, &
      & 0.1318913820322_DP,0.1320603088411_DP,0.132229235650_DP,0.132399074312_DP,0.132568954822_DP,0.132738835332_DP, &
      & 0.1329087158420_DP,0.1330788764544_DP,0.133249734060_DP,0.133420591667_DP,0.133591449273_DP,0.133762306879_DP, &
      & 0.1339336992411_DP,0.1341055553883_DP,0.134277411535_DP,0.134449267682_DP,0.134621123829_DP,0.134793694483_DP, &
      & 0.1349665687658_DP,0.1351394430487_DP,0.135312317331_DP,0.135485191614_DP,0.135658878561_DP,0.135832788820_DP, &
      & 0.1360066990800_DP,0.1361806093395_DP,0.136354519599_DP,0.136529253243_DP,0.136704215658_DP,0.136879178072_DP, &
      & 0.1370541404878_DP,0.1372291029027_DP,0.137404806924_DP,0.137580836098_DP,0.137756865271_DP,0.137932894445_DP, &
      & 0.1381089236188_DP,0.1382855157880_DP,0.138462624830_DP,0.138639733872_DP,0.138816842915_DP,0.138993951957_DP, &
      & 0.1391713448843_DP,0.1393495454909_DP,0.139527746097_DP,0.139705946704_DP,0.139884147310_DP,0.140062347917_DP, &
      & 0.1402415516727_DP,0.1404208541990_DP,0.140600156725270_DP,0.140779459251523_DP,0.140958761777775_DP]

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      PROBLEM=>CONTROL_LOOP%PROBLEM
      IF(ASSOCIATED(PROBLEM)) THEN
        SELECT CASE(PROBLEM%SPECIFICATION(3))
        CASE(PROBLEM_GUDUNOV_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, PROBLEM_STRANG_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE)
          IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP_PARENT,ERR,ERROR,*999)
            !The first control_loop is the one for monodomain
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP_PARENT,1,CONTROL_LOOP_MONODOMAIN,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP_MONODOMAIN,SOLVERS,ERR,ERROR,*999)
            !The second solver is associated with the diffusion part of the monodomain equation
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
            SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
              IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR
                IF(ASSOCIATED(EQUATIONS_SET)) THEN
                  INDEPENDENT_FIELD_MONODOMAIN=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
                  IF(.NOT.ASSOCIATED(INDEPENDENT_FIELD_MONODOMAIN)) THEN
                    CALL FlagError("Independent field is not associated.",ERR,ERROR,*999)
                  ENDIF

                  CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U1_VARIABLE_TYPE,FIELD_VAR_IND_M_U1,ERR,ERROR,*999)

                  M_NODES_MAPPING=>INDEPENDENT_FIELD_MONODOMAIN%DECOMPOSITION%DOMAIN(INDEPENDENT_FIELD_MONODOMAIN%DECOMPOSITION% &
                    & MESH_COMPONENT_NUMBER)%PTR%MAPPINGS%NODES

                  ! Initialization
                  INDEX_REF=1

                  DO node_idx=1,M_NODES_MAPPING%NUMBER_OF_LOCAL
                  
                    !the fourth component of the U1 variable contains the half sarcomere length at activation
                    DofIdx=FIELD_VAR_IND_M_U1%COMPONENTS(4)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                      & VERSIONS(1)
                    CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U1_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,DofIdx,SARCO_LENGTH_AT_ACTIVATION,ERR,ERROR,*999)

                    !the first component of the U1 variable contains the actual half sarcomere length
                    DofIdx=FIELD_VAR_IND_M_U1%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)% &
                      & VERSIONS(1)
                    CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(INDEPENDENT_FIELD_MONODOMAIN,FIELD_U1_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,DofIdx,SARCO_LENGTH,ERR,ERROR,*999)
                    
                    ELONGATION=SARCO_LENGTH-SARCO_LENGTH_AT_ACTIVATION
                    
                    IF(ELONGATION.LT.0) THEN
                      LENGTH_TITIN=SARCO_LENGTH-LENGTH_MYOSIN-LENGTH_MBAND-LENGTH_ZDISC
                      !function to approximate the relation between the initial titin length and the initial passive force F0	
                      TITIN_UNBOUND=COEFF_MATRIX(1)*EXP(LENGTH_TITIN)+COEFF_MATRIX(2)*LENGTH_TITIN**3+COEFF_MATRIX(3)* &
                        & LENGTH_TITIN**2+COEFF_MATRIX(4)*LENGTH_TITIN+COEFF_MATRIX(5) 
                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN, &
                        & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,node_idx,2,TITIN_UNBOUND,ERR,ERROR,*999)
                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN, &
                        & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,node_idx,3,TITIN_UNBOUND,ERR,ERROR,*999)

                      ! calculate x-fibre titin stress (trigonometry):                      
                      ! fitting function to calculate the filament-lattice parameter d10 in nanometer (from Elliot 1963 -mammalian muscle)
                      d10=-13.39_DP*SARCO_LENGTH+58.37_DP
                      ! d10=-13.39_DP*1.0_DP+58.37_DP
                      ! calculate the distance between actin and myosin filament in micro meter (geometrical relationship)
                      ACTIN_MYOSIN_DISTANCE=0.001_DP*(2.0_DP/3.0_DP*d10)  
                      ! calculate x-fibre stress with tangens-function (unbound titin)
                      TITIN_XF_UNBOUND=0.5_DP*TITIN_UNBOUND*ACTIN_MYOSIN_DISTANCE/LENGTH_TITIN
                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN, &
                        & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,node_idx,4,TITIN_XF_UNBOUND,ERR,ERROR,*999)
                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN, &
                        & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,node_idx,5,TITIN_XF_UNBOUND,ERR,ERROR,*999)

                    ELSE  ! Force enhancement  

                      ! Calculate Titin-Force for unbound titin filaments
                      LENGTH_TITIN=SARCO_LENGTH-LENGTH_MYOSIN-LENGTH_MBAND-LENGTH_ZDISC
                      !function to approximate the relation between the initial titin length and the initial passive force F0	
                      TITIN_UNBOUND=COEFF_MATRIX(1)*EXP(LENGTH_TITIN)+COEFF_MATRIX(2)*LENGTH_TITIN**3+COEFF_MATRIX(3)* &
                        & LENGTH_TITIN**2+COEFF_MATRIX(4)*LENGTH_TITIN+COEFF_MATRIX(5)

                      ! Calculate Titin-Force for bound titin filaments
                      ! Switch between different force-enhancent models/implementations 
                      ! 1=linear approximation 2=square approximation 3=simple iterative solution of Rode's model
                      ! 4=Newton's-method to solve Rode's model
                      SWITCH_MODEL=4
                      !linear approximation
                      IF(SWITCH_MODEL.EQ.1) THEN
                        LENGTH_INIT_TITIN=SARCO_LENGTH_AT_ACTIVATION-LENGTH_MYOSIN-LENGTH_MBAND-LENGTH_ZDISC
                        !function to approximate the relation between the initial titin length and the initial passive force F0	
                        F0=COEFF_MATRIX(1)*EXP(LENGTH_INIT_TITIN)+COEFF_MATRIX(2)*LENGTH_INIT_TITIN**3+COEFF_MATRIX(3)* &
                          & LENGTH_INIT_TITIN**2+COEFF_MATRIX(4)*LENGTH_INIT_TITIN+COEFF_MATRIX(5)
                        !function to approximate the relation between the initial sarcomere length and the stiffness of the PEVK region.
                        STIFFNESS_PEVK=1000.0_DP*(0.1880_DP*SARCO_LENGTH_AT_ACTIVATION**4-0.8694_DP*SARCO_LENGTH_AT_ACTIVATION**3+ &
                          & 1.5084_DP*SARCO_LENGTH_AT_ACTIVATION**2-1.1577_DP*SARCO_LENGTH_AT_ACTIVATION+0.3345_DP)
                        IF(ELONGATION.LT.0.02) THEN ! 0.02 -> offset value 
                          FORCE=0.0_DP
                        ELSE 
                          FORCE=(ELONGATION-0.02_DP)*STIFFNESS_PEVK
                        ENDIF
                        TITIN_BOUND=F0+FORCE

                      !square approximation
                      ELSEIF(SWITCH_MODEL.EQ.2) THEN
                        LENGTH_INIT_TITIN=SARCO_LENGTH_AT_ACTIVATION-LENGTH_MYOSIN-LENGTH_MBAND-LENGTH_ZDISC
                        !function to approximate the relation between the initial titin length and the initial passive force F0	
                        F0=COEFF_MATRIX(1)*EXP(LENGTH_INIT_TITIN)+COEFF_MATRIX(2)*LENGTH_INIT_TITIN**3+COEFF_MATRIX(3)* &
                          & LENGTH_INIT_TITIN**2+COEFF_MATRIX(4)*LENGTH_INIT_TITIN+COEFF_MATRIX(5)
                        FORCE=2.0_DP*ELONGATION**2
                        TITIN_BOUND=F0+FORCE

                      !Rode's model -> Ekin's implementation
                      ELSEIF(SWITCH_MODEL.EQ.3) THEN
                        LENGTH_INIT_TITIN=SARCO_LENGTH_AT_ACTIVATION-LENGTH_MYOSIN-LENGTH_MBAND-LENGTH_ZDISC
                        !function to approximate the relation between the initial titin length and the initial passive force F0	
                        F0=COEFF_MATRIX(1)*EXP(LENGTH_INIT_TITIN)+COEFF_MATRIX(2)*LENGTH_INIT_TITIN**3+COEFF_MATRIX(3)* &
                          & LENGTH_INIT_TITIN**2+COEFF_MATRIX(4)*LENGTH_INIT_TITIN+COEFF_MATRIX(5)
                        !function to approximate the relation between the initial sarcomere length and the stiffness of the PEVK region.
                        STIFFNESS_PEVK=1000.0_DP*(0.1880_DP*SARCO_LENGTH_AT_ACTIVATION**4-0.8694_DP*SARCO_LENGTH_AT_ACTIVATION**3+ &
                          & 1.5084_DP*SARCO_LENGTH_AT_ACTIVATION**2-1.1577_DP*SARCO_LENGTH_AT_ACTIVATION+0.3345_DP)
                      
                        IF(F0.LE.0) THEN
                          LENGTH_DIST_IG_F0=-35.63_DP+39.58889_DP*(F0+0.9_DP)
                        ELSEIF(F0.GE.0.24_DP) THEN
                          LENGTH_DIST_IG_F0=0.1411_DP+0.196576763485477_DP*(F0-0.2495_DP)
                        ELSE                    
                          INDEX_PSEUDO=CEILING(F0/DX)
                          INDEX_I=INDEX_REF+INDEX_PSEUDO-1
                          LENGTH_DIST_IG_F0=LENGTHS_DIST_IG(INDEX_I)-(LENGTHS_DIST_IG(INDEX_I+1)-LENGTHS_DIST_IG(INDEX_I))* &
                            & (FORCES_DIST_IG(INDEX_I)-F0)/(FORCES_DIST_IG(INDEX_I+1)-FORCES_DIST_IG(INDEX_I))
                        ENDIF          
	               
                        ELONGATION_NEW=-1.0_DP
                        FORCE=0.0_DP                  
                        DO WHILE (ELONGATION_NEW.LT.ELONGATION)
                          FORCE=FORCE+FORCE_INCREMENT
                          TITIN_BOUND=FORCE+F0
                          ELONGATION_PEVK=FORCE/STIFFNESS_PEVK  
                          IF (TITIN_BOUND.LE.0.0_DP) THEN
                            LENGTH_DIST_IG=-35.63_DP+39.58889_DP*(TITIN_BOUND+0.9_DP)
                          ELSE IF (TITIN_BOUND.GE.0.24_DP) THEN
                            LENGTH_DIST_IG=0.1411_DP+0.196576763485477_DP*(TITIN_BOUND-0.2495_DP)
                          ELSE                  
                            INDEX_PSEUDO=CEILING(TITIN_BOUND/DX)
                            INDEX_I=INDEX_REF+INDEX_PSEUDO-1
                            LENGTH_DIST_IG=LENGTHS_DIST_IG(INDEX_I)-(LENGTHS_DIST_IG(INDEX_I+1)-LENGTHS_DIST_IG( &
                              & INDEX_I))*(FORCES_DIST_IG(INDEX_I)-TITIN_BOUND)/(FORCES_DIST_IG(INDEX_I+1)-FORCES_DIST_IG(INDEX_I))
                          END IF
                          ELONGATION_DIST_IG=LENGTH_DIST_IG-LENGTH_DIST_IG_F0                  
                          ELONGATION_NEW=ELONGATION_PEVK+ELONGATION_DIST_IG
                        ENDDO

                      !Rode's titin model -> solve with Newton's method
                      ELSEIF(SWITCH_MODEL.EQ.4) THEN
                        LENGTH_INIT_TITIN=SARCO_LENGTH_AT_ACTIVATION-LENGTH_MYOSIN-LENGTH_MBAND-LENGTH_ZDISC
                        !function to approximate the relation between the initial titin length and the initial passive force F0	
                        F0=COEFF_MATRIX(1)*EXP(LENGTH_INIT_TITIN)+COEFF_MATRIX(2)*LENGTH_INIT_TITIN**3+COEFF_MATRIX(3)* &
                          & LENGTH_INIT_TITIN**2+COEFF_MATRIX(4)*LENGTH_INIT_TITIN+COEFF_MATRIX(5)
                        !function to approximate the relation between the initial sarcomere length and the stiffness of the PEVK region.
                        STIFFNESS_PEVK=1000.0_DP*(0.1880_DP*SARCO_LENGTH_AT_ACTIVATION**4-0.8694_DP*SARCO_LENGTH_AT_ACTIVATION**3+ &
                          & 1.5084_DP*SARCO_LENGTH_AT_ACTIVATION**2-1.1577_DP*SARCO_LENGTH_AT_ACTIVATION+0.3345_DP)

                        !calculate LENGTH_DIST_IG_F0 with linear inter- or extrapolation
                        INDEX_I=2
                        SLOPE=0.0_DP
                        LENGTH_DIST_IG_F0=0.0_DP
                        IF(F0.GE.FORCES_DIST_IG(DIM_DATA)) THEN
                          INDEX_I=DIM_DATA
                        ELSE
                          DO WHILE (F0.GT.FORCES_DIST_IG(INDEX_I))
                            INDEX_I=INDEX_I+1
                          ENDDO
                        ENDIF
                        SLOPE=(FORCES_DIST_IG(INDEX_I)-FORCES_DIST_IG(INDEX_I-1))/ &
                          & (LENGTHS_DIST_IG(INDEX_I)-LENGTHS_DIST_IG(INDEX_I-1))
                        LENGTH_DIST_IG_F0=LENGTHS_DIST_IG(INDEX_I-1)+SLOPE*(F0-FORCES_DIST_IG(INDEX_I-1))

                        !initialize Newton-method to calculate the titin force
                        INDEX_I=2                
                        STIFFNESS_DIST=1.0_DP   
                        FORCE_DISTAL_IG=F0      
                        FORCE=0.0_DP                  ! delta P
                        TITIN_BOUND=0.0_DP            ! total titin force (= P_PEVK)
                        DELTA_F=10.0_DP               ! start value for iteration
                        DIFF_QUOT=1.0_DP              ! numerical derivative              
                        ELONGATION_PEVK=ELONGATION/2.0_DP                               
                        LENGTH_DIST_IG=LENGTH_DIST_IG_F0+ELONGATION-ELONGATION_PEVK     

                        DO WHILE(ABS(DELTA_F).GT.TOL) !Newton-method (solve FORCE_DISTAL_IG-FORCE_0=0)

                          IF(LENGTH_DIST_IG.GE.LENGTHS_DIST_IG(DIM_DATA)) THEN  !Extrapolation if LENGTHS_DIST_IG(end)>LENGTH_DIST_IG
                            INDEX_I=DIM_DATA
                          ELSE 
                            DO WHILE(LENGTH_DIST_IG.GT.LENGTHS_DIST_IG(INDEX_I))
                              INDEX_I=INDEX_I+1
                            ENDDO
                          ENDIF
                          STIFFNESS_DIST=(FORCES_DIST_IG(INDEX_I)-FORCES_DIST_IG(INDEX_I-1))/ &
                            & (LENGTHS_DIST_IG(INDEX_I)-LENGTHS_DIST_IG(INDEX_I-1))
                          FORCE_DISTAL_IG=FORCES_DIST_IG(INDEX_I-1)+STIFFNESS_DIST* &
                            & (LENGTH_DIST_IG-LENGTHS_DIST_IG(INDEX_I-1))                

                          FORCE=STIFFNESS_PEVK*ELONGATION_PEVK
                          TITIN_BOUND=FORCE+F0

                          DELTA_F=TITIN_BOUND-FORCE_DISTAL_IG !new iteration for FORCE_DISTAL_IG-TITIN_BOUND
                          DIFF_QUOT=STIFFNESS_PEVK+STIFFNESS_DIST                 
                          ELONGATION_PEVK=ELONGATION_PEVK-DELTA_F/DIFF_QUOT                 
                          LENGTH_DIST_IG=LENGTH_DIST_IG_F0+ELONGATION-ELONGATION_PEVK       
                        ENDDO
                      ENDIF ! switch model

                      ! calculate x-fibre titin stress                       
                      ! fitting function to calculate the filament-lattice parameter d10 in nanometer (from Elliot 1963 -mammalian muscle)
                      d10=-13.39_DP*SARCO_LENGTH+58.37_DP
                      ! d10=-13.39_DP*1.0_DP+58.37_DP
                      ! calculate the distance between actin and myosin filament in micro meter (geometrical relationship)
                      ACTIN_MYOSIN_DISTANCE=0.001_DP*(2.0_DP/3.0_DP*d10) 
                      ! calculate x-fibre stress with tangent (unbound titin)
                      TITIN_XF_UNBOUND=0.5_DP*TITIN_UNBOUND*ACTIN_MYOSIN_DISTANCE/LENGTH_TITIN 
                      ! calculate x-fibre stress with tangent (for bound titin filaments)
                      TITIN_XF_BOUND=0.5_DP*TITIN_BOUND*ACTIN_MYOSIN_DISTANCE/LENGTH_DIST_IG  
                      
                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN, &
                        & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,node_idx,2,TITIN_UNBOUND,ERR,ERROR,*999)
                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN, &
                        & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,node_idx,3,TITIN_BOUND,ERR,ERROR,*999)
                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN, &
                        & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,node_idx,4,TITIN_XF_UNBOUND,ERR,ERROR,*999)
                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(INDEPENDENT_FIELD_MONODOMAIN, &
                        & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,node_idx,5,TITIN_XF_BOUND,ERR,ERROR,*999)

                    ENDIF ! Check if elongation is positive or not
                  ENDDO ! Over the nodes

                  !now the ghost elements -- get the relevant info from the other computational nodes
                  CALL FIELD_PARAMETER_SET_UPDATE_START(INDEPENDENT_FIELD_MONODOMAIN, & 
                    & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                  CALL FIELD_PARAMETER_SET_UPDATE_FINISH(INDEPENDENT_FIELD_MONODOMAIN, & 
                    & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)

                ELSE
                  CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Solver mapping is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Solver equations is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDIF
        CASE DEFAULT
          CALL FlagError("Problem subtype not implemented for titin",ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FlagError("Problem not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("BIOELECTRIC_FINITE_ELASTICITY_COMPUTE_TITIN")
    RETURN
999 ERRORSEXITS("BIOELECTRIC_FINITE_ELASTICITY_COMPUTE_TITIN",ERR,ERROR)
    RETURN 1

  END SUBROUTINE BIOELECTRIC_FINITE_ELASTICITY_COMPUTE_TITIN

  !
  !================================================================================================================================
  !

END MODULE BIOELECTRIC_FINITE_ELASTICITY_ROUTINES

