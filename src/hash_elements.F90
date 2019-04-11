
  !>Calculates the local/global element mappings for a domain decomposition.
  SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_CALCULATE(DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<A pointer to the domain to calculate the element mappings for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,no_adjacent_element,adjacent_element,domain_no,domain_idx,ne,nn,np,NUMBER_OF_DOMAINS, &
      & NUMBER_OF_ADJACENT_ELEMENTS,myComputationalNodeNumber,component_idx, domainLocalType, &
      & numberOfHashKeys
    INTEGER(INTG), ALLOCATABLE :: ADJACENT_ELEMENTS(:),DOMAINS(:),LOCAL_ELEMENT_NUMBERS(:), &
      & hashKeysArray(:), hashValuesArray(:), hashIntegerArray(:)
    INTEGER(INTG), ALLOCATABLE :: hashValuesMatrix(:,:),hashValuesSubMatrix(:,:)
    TYPE(LIST_TYPE), POINTER :: ADJACENT_DOMAINS_LIST
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: ADJACENT_ELEMENTS_LIST(:), hashKeysList(:)!, hashValuesList(:)
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    ENTERS("DOMAIN_MAPPINGS_ELEMENTS_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      IF(ASSOCIATED(DOMAIN%MAPPINGS)) THEN
        IF(ASSOCIATED(DOMAIN%MAPPINGS%ELEMENTS)) THEN
          ELEMENTS_MAPPING=>DOMAIN%MAPPINGS%ELEMENTS
          IF(ASSOCIATED(DOMAIN%DECOMPOSITION)) THEN
            DECOMPOSITION=>DOMAIN%DECOMPOSITION
            IF(ASSOCIATED(DOMAIN%MESH)) THEN
              MESH=>DOMAIN%MESH
              component_idx=DOMAIN%MESH_COMPONENT_NUMBER
              myComputationalNodeNumber=ComputationalEnvironment_NodeNumberGet(ERR,ERROR)
              IF(ERR/=0) GOTO 999

              !Calculate the local and global numbers and set up the mappings
              ALLOCATE(ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(MESH%NUMBER_OF_ELEMENTS),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate element mapping global to local map.",ERR,ERROR,*999)
              !ALLOCATE(ELEMENTS_MAPPING%domainMappingHashes(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
              !Should be only local!
              ALLOCATE(ELEMENTS_MAPPING%domainMappingHashes(1),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate hash tables.",ERR,ERROR,*999)

              ELEMENTS_MAPPING%NUMBER_OF_GLOBAL=MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%NUMBER_OF_ELEMENTS
              !Loop over the global elements and calculate local numbers
              ALLOCATE(LOCAL_ELEMENT_NUMBERS(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate local element numbers.",ERR,ERROR,*999)
              LOCAL_ELEMENT_NUMBERS=0
              ALLOCATE(ADJACENT_ELEMENTS_LIST(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate adjacent elements list.",ERR,ERROR,*999)

              !ALLOCATE(hashKeysList(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
              ALLOCATE(hashKeysList(1),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate hash keys list.",ERR,ERROR,*999)
!              ALLOCATE(hashValuesList(ELEMENTS_MAPPING%NUMBER_OF_GLOBAL),STAT=ERR)
!              IF(ERR/=0) CALL FlagError("Could not allocate hashValuesList.",ERR,ERROR,*999)
              ALLOCATE(hashValuesMatrix(DECOMPOSITION%NUMBER_OF_DOMAINS*3+1,MESH%NUMBER_OF_ELEMENTS),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate hashValuesMatrix.",ERR,ERROR,*999)
              hashValuesMatrix=0

              NULLIFY(hashKeysList(1)%PTR)
              CALL LIST_CREATE_START(hashKeysList(1)%PTR,ERR,ERROR,*999)
              CALL LIST_DATA_TYPE_SET(hashKeysList(1)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
              CALL LIST_CREATE_FINISH(hashKeysList(1)%PTR,ERR,ERROR,*999)

              DO domain_idx=0,DECOMPOSITION%NUMBER_OF_DOMAINS-1

                NULLIFY(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR)
                CALL LIST_CREATE_START(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,ERR,ERROR,*999)
                CALL LIST_DATA_TYPE_SET(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                CALL LIST_INITIAL_SIZE_SET(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,MAX(INT(MESH%NUMBER_OF_ELEMENTS/2),1), &
                  & ERR,ERROR,*999)
                CALL LIST_CREATE_FINISH(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,ERR,ERROR,*999)

                !NULLIFY(hashKeysList(domain_idx)%PTR)
                !CALL LIST_CREATE_START(hashKeysList(domain_idx)%PTR,ERR,ERROR,*999)
                !CALL LIST_DATA_TYPE_SET(hashKeysList(domain_idx)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                !CALL LIST_INITIAL_SIZE_SET(hashKeysList(domain_idx)%PTR,MAX(INT(MESH%NUMBER_OF_ELEMENTS/2),1), &
                !  & ERR,ERROR,*999)
                !CALL LIST_CREATE_FINISH(hashKeysList(domain_idx)%PTR,ERR,ERROR,*999)

              ENDDO !domain_idx

              DO ne=1,MESH%NUMBER_OF_ELEMENTS
                !Calculate the local numbers
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

                !NULLIFY(hashValuesList(ne)%PTR)
                !CALL LIST_CREATE_START(hashValuesList(ne)%PTR,ERR,ERROR,*999)
                !CALL LIST_DATA_TYPE_SET(hashValuesList(ne)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                !CALL List_MutableSet(hashValuesList(ne)%PTR,.TRUE., err, error, *999)
                !CALL LIST_INITIAL_SIZE_SET(hashValuesList(ne)%PTR,4,ERR,ERROR,*999)
                !CALL LIST_CREATE_FINISH(hashValuesList(ne)%PTR,ERR,ERROR,*999)

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
                IF(NUMBER_OF_DOMAINS==1) THEN
                  !Element is an internal element
                  ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_TYPE(1)=DOMAIN_LOCAL_INTERNAL
                  domainLocalType = DOMAIN_LOCAL_INTERNAL
                ELSE
                  !Element is on the boundary of computational domains
                  ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_TYPE(1)=DOMAIN_LOCAL_BOUNDARY
                  domainLocalType = DOMAIN_LOCAL_BOUNDARY
                ENDIF
                !Add local element
                !This info should be collected in a different way once gtl does not exist any more!
                hashValuesMatrix(2:4,ne)= &
                    & [domain_no,LOCAL_ELEMENT_NUMBERS(domain_no),domainLocalType]
                hashValuesMatrix(1,ne)=1

                !CALL LIST_ITEM_ADD (hashValuesList(ne)%PTR, &
                !  & 1,err,error,*999) ! First element is number of domains
                !CALL LIST_ITEM_ADD (hashValuesList(ne)%PTR, &
                !  & domain_no,err,error,*999)
                !CALL LIST_ITEM_ADD (hashValuesList(ne)%PTR, &
                !  & LOCAL_ELEMENT_NUMBERS(domain_no),err,error,*999)
                !CALL LIST_ITEM_ADD (hashValuesList(ne)%PTR, &
                !  & domainLocalType,err,error,*999)
                !CALL LIST_ITEM_ADD (hashKeysList(domain_no)%PTR, ne, err, error,*999)
                IF(domain_no==myComputationalNodeNumber) THEN
                  CALL LIST_ITEM_ADD (hashKeysList(1)%PTR, ne, err, error,*999)
                END IF
                !This info should be collected in a different way once gtl does not exist any more!
                !hashValuesMatrix(2:4,sizeOfKeyList)= &
                !  & [domain_no,LOCAL_ELEMENT_NUMBERS(domain_no),domainLocalType]
                !hashValuesMatrix(1,sizeOfKeyList)=1
                !END IF
              ENDDO !ne

              !Compute ghost element mappings
              !Hash: we can keep adding keys because new elements are ghosts => no risk of doubles
              DO domain_idx=0,DECOMPOSITION%NUMBER_OF_DOMAINS-1
                CALL LIST_REMOVE_DUPLICATES(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,ERR,ERROR,*999)
                CALL LIST_DETACH_AND_DESTROY(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,NUMBER_OF_ADJACENT_ELEMENTS, &
                  & ADJACENT_ELEMENTS,ERR,ERROR,*999)
                DO no_adjacent_element=1,NUMBER_OF_ADJACENT_ELEMENTS
                  adjacent_element=ADJACENT_ELEMENTS(no_adjacent_element)
                  LOCAL_ELEMENT_NUMBERS(domain_idx)=LOCAL_ELEMENT_NUMBERS(domain_idx)+1
                  ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%NUMBER_OF_DOMAINS= &
                    & ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%NUMBER_OF_DOMAINS+1
                  ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%LOCAL_NUMBER( &
                    & ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%NUMBER_OF_DOMAINS)=LOCAL_ELEMENT_NUMBERS(domain_idx)
                  ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%DOMAIN_NUMBER( &
                    & ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%NUMBER_OF_DOMAINS)=domain_idx
                  ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%LOCAL_TYPE( &
                    & ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%NUMBER_OF_DOMAINS)= &
                    & DOMAIN_LOCAL_GHOST
                  !Add keys and values to the table
                  !hash_t(domain_idx)_Skey(idx)=adjacent_element !check for double entry? (no)
                  !CALL LIST_ITEM_ADD(hashKeysList(domain_idx)%PTR,adjacent_element,ERR,ERROR,*999)
                  !This info should be collected in a different way once gtl does not exist any more!
                  hashValuesMatrix(hashValuesMatrix(1,adjacent_element)*3+2:4,adjacent_element)= &
                    & [domain_idx,LOCAL_ELEMENT_NUMBERS(domain_idx),DOMAIN_LOCAL_GHOST]
                  hashValuesMatrix(1,adjacent_element)=hashValuesMatrix(1,adjacent_element)+1
                  IF(domain_idx==myComputationalNodeNumber) THEN
                    CALL LIST_ITEM_ADD (hashKeysList(1)%PTR, adjacent_element, err, error,*999)
                  END IF
                  !hashValuesMatrix(hashValuesMatrix(1,sizeOfKeyList)*3+2:4,sizeOfKeyList)= &
                  !    & [domain_idx,LOCAL_ELEMENT_NUMBERS(domain_idx),DOMAIN_LOCAL_GHOST]
                  !  hashValuesMatrix(1,sizeOfKeyList)=hashValuesMatrix(1,adjacent_element)+1
                  !END IF

!                  CALL LIST_ITEM_ADD (hashValuesList(adjacent_element)%PTR, &
!                    & domain_idx,err,error,*999)
!                  CALL LIST_ITEM_ADD (hashValuesList(adjacent_element)%PTR, &
!                    & LOCAL_ELEMENT_NUMBERS(domain_idx),err,error,*999)
!                  CALL LIST_ITEM_ADD (hashValuesList(adjacent_element)%PTR, &
!                    & DOMAIN_LOCAL_GHOST,err,error,*999)
                  !Change number of domains
!                  CALL LIST_ITEM_GET (hashValuesList(adjacent_element)%PTR, &
!                    & 1,numberOfDomainsElement,err,error,*999)
!                  CALL LIST_ITEM_SET (hashValuesList(adjacent_element)%PTR, &
!                    & 1,numberOfDomainsElement+1,err,error,*999)

                ENDDO !no_adjacent_element
                IF(ALLOCATED(ADJACENT_ELEMENTS)) DEALLOCATE(ADJACENT_ELEMENTS)

              ENDDO !domain_idx

              DEALLOCATE(ADJACENT_ELEMENTS_LIST)
              DEALLOCATE(LOCAL_ELEMENT_NUMBERS)

              !Transfer the lists to a matrix
              !DO ne = 1, MESH%NUMBER_OF_ELEMENTS
              !  CALL LIST_DETACH_AND_DESTROY(hashValuesList(ne)%PTR,numberOfValues, &
              !    & hashIntegerArray,ERR,ERROR,*999)
              !  ALLOCATE(hashValuesArray(numberOfValues))
              !  hashValuesArray=hashIntegerArray(1:numberOfValues)
              !  IF(ALLOCATED(hashIntegerArray)) DEALLOCATE(hashIntegerArray)
              !  DO rowMat = 1,numberOfValues
              !    hashValuesMatrix(rowMat,ne)=hashValuesArray(rowMat)
              !  END DO
              !  DEALLOCATE(hashValuesArray)
              !END DO

              CALL LIST_DETACH_AND_DESTROY(hashKeysList(1)%PTR,numberOfHashKeys, &
                & hashIntegerArray,ERR,ERROR,*999)

              ALLOCATE(hashKeysArray(numberOfHashKeys))
              hashKeysArray=hashIntegerArray(1:numberOfHashKeys)
              IF(ALLOCATED(hashIntegerArray)) DEALLOCATE(hashIntegerArray)

              WRITE(*,*) "NofKeys"
              WRITE(*,*) numberOfHashKeys
              WRITE(*,*) "HashKarray"
              WRITE(*,*) hashKeysArray

              ALLOCATE(hashValuesSubMatrix(DECOMPOSITION%NUMBER_OF_DOMAINS*3+1,numberOfHashKeys),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate hashValuesMatrix.",ERR,ERROR,*999)
              hashValuesSubMatrix=0

              hashValuesSubMatrix = hashValuesMatrix(:,hashKeysArray)
              WRITE(*,*) "HashVmatrix1"
              WRITE(*,*) hashValuesSubMatrix(:,1)
              WRITE(*,*) "HashVmatrix2"
              WRITE(*,*) hashValuesSubMatrix(:,2)

              !Finally compute the table
              NULLIFY(ELEMENTS_MAPPING%domainMappingHashes(1)%PTR)
              CALL HashTable_CreateStart(ELEMENTS_MAPPING%domainMappingHashes(1)%PTR,ERR,ERROR,*999)
              ! define some parameters here if needed
              CALL HashTable_CreateFinish(ELEMENTS_MAPPING%domainMappingHashes(1)%PTR,ERR,ERROR,*999)

              CALL HashTable_ValuesSetAndInsert(ELEMENTS_MAPPING%domainMappingHashes(1)%PTR, &
                & hashKeysArray,hashValuesSubMatrix, .FALSE., ERR, ERROR, *999)

              IF(ALLOCATED(hashValuesSubMatrix)) DEALLOCATE(hashValuesSubMatrix)
              IF(ALLOCATED(hashKeysArray)) DEALLOCATE(hashKeysArray)
              IF(ALLOCATED(hashValuesMatrix)) DEALLOCATE(hashValuesMatrix)
              IF(ALLOCATED(hashKeysList)) DEALLOCATE(hashKeysList)

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
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%NUMBER_OF_DOMAINS+1,8,8, &
        & ELEMENTS_MAPPING%ADJACENT_DOMAINS_PTR,'("    Adjacent domains ptr  :",8(X,I7))','(27X,8(X,I7))',ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%ADJACENT_DOMAINS_PTR( &
        & ELEMENTS_MAPPING%NUMBER_OF_DOMAINS)-1,8,8,ELEMENTS_MAPPING%ADJACENT_DOMAINS_LIST, &
        '("    Adjacent domains list :",8(X,I7))','(27X,8(X,I7))',ERR,ERROR,*999)
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
    IF(ALLOCATED(LOCAL_ELEMENT_NUMBERS)) DEALLOCATE(LOCAL_ELEMENT_NUMBERS)
    IF(ALLOCATED(hashKeysArray)) DEALLOCATE(hashKeysArray)
    IF(ALLOCATED(hashValuesArray)) DEALLOCATE(hashValuesArray)
    IF(ALLOCATED(hashIntegerArray)) DEALLOCATE(hashIntegerArray)
    IF(ALLOCATED(hashValuesSubMatrix)) DEALLOCATE(hashValuesSubMatrix)
    IF(ALLOCATED(hashValuesMatrix)) DEALLOCATE(hashValuesMatrix)
    IF(ASSOCIATED(DOMAIN%MAPPINGS%ELEMENTS)) CALL DOMAIN_MAPPINGS_ELEMENTS_FINALISE(DOMAIN%MAPPINGS,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("DOMAIN_MAPPINGS_ELEMENTS_CALCULATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_CALCULATE
