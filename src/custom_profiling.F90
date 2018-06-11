
MODULE Custom_Profiling

#include "macros.h"
  USE KINDS
  IMPLICIT NONE
  PRIVATE
#include "mpif.h"

  PUBLIC :: CustomProfilingStart
  PUBLIC :: CustomProfilingStop
  PUBLIC :: CustomProfilingMemory
  PUBLIC :: CustomProfilingGetInfo
  PUBLIC :: CustomProfilingGetDuration
  PUBLIC :: CustomProfilingGetMemory
  PUBLIC :: CustomProfilingGetSizePerElement
  PUBLIC :: CustomProfilingGetNumberObjects
  PUBLIC :: CustomProfilingReset
  PRIVATE :: GetDurationIndex
  PRIVATE :: GetMemoryIndex
  PRIVATE :: PrintWarningDuration
  PRIVATE :: PrintWarningMemory

  INTEGER, PARAMETER :: IDENTIFIER_LENGTH = 80
  INTEGER, PARAMETER :: NUMBER_OF_RECORDS = 200
  LOGICAL, PARAMETER ::  WITH_MEMORY_CONSUMPTION_MEASUREMENT = .FALSE.   ! enable this for measuring of memory consumption increase between timing markers

  CHARACTER(LEN=IDENTIFIER_LENGTH), DIMENSION(NUMBER_OF_RECORDS) :: DurationIdentifiers    !< identifiers for duration records
  CHARACTER(LEN=IDENTIFIER_LENGTH), DIMENSION(NUMBER_OF_RECORDS) :: MemoryIdentifiers      !< identifier for memory records
  REAL(DP), DIMENSION(NUMBER_OF_RECORDS) :: Durations
  REAL(DP), DIMENSION(NUMBER_OF_RECORDS) :: StartTime
  INTEGER(LINTG), DIMENSION(NUMBER_OF_RECORDS) :: MemoryConsumptions   !< Memory Consumption at call to CustomProfilingMemory
  INTEGER(INTG), DIMENSION(NUMBER_OF_RECORDS) :: SizesPerElement       !< size of an element at call to CustomProfilingMemory
  INTEGER(INTG), DIMENSION(NUMBER_OF_RECORDS) :: NumberOfObjects       !< number of objects
  INTEGER(INTG), DIMENSION(NUMBER_OF_RECORDS) :: TimeCount             !< how often a CustomProfilingStart was called for each identifier
  INTEGER(LINTG), DIMENSION(NUMBER_OF_RECORDS) :: StartMemory     !<   memory (resident set size) at call to last CustomProfilingStart, in multiples of pagesize (=2048)
  INTEGER(LINTG), DIMENSION(NUMBER_OF_RECORDS) :: TotalMemory     !<   total new memory consumption (RSS) between call to CustomProfilingStart and CustomProfilingEnd, in multiples of pagesize (=2048)
  INTEGER :: SizeDuration = 0
  INTEGER :: SizeMemory = 0
  INTEGER(LINTG) :: CurrentMemoryConsumption

CONTAINS

  !
  !================================================================================================================================
  !
  SUBROUTINE CustomProfilingStart(Identifier)
    ! PARAMETERS
    CHARACTER(LEN=*), INTENT(IN)  :: Identifier !< A custom Identifier that describes the timer

    ! LOCAL VARIABLES
    INTEGER :: CurrentIndex

    ! find index of identifier
    CurrentIndex = GetDurationIndex(Identifier)

    ! If record with this identifier does not yet exist, create new
    IF (CurrentIndex == 0) THEN
      SizeDuration = SizeDuration + 1
      CurrentIndex = SizeDuration
      DurationIdentifiers(CurrentIndex) = Identifier
      Durations(CurrentIndex) = 0.0_8
      TimeCount(CurrentIndex) = 0
      
      IF (WITH_MEMORY_CONSUMPTION_MEASUREMENT) THEN
        StartMemory(CurrentIndex) = GetCurrentMemoryConsumption()
        TotalMemory(CurrentIndex) = 0
      ENDIF
    ENDIF

    !CALL CPU_TIME(StartTime(CurrentIndex))
    StartTime(CurrentIndex) = MPI_WTIME()
  END SUBROUTINE
  !
  !================================================================================================================================
  !
  SUBROUTINE CustomProfilingStop(Identifier)
    ! PARAMETERS
    CHARACTER(LEN=*), INTENT(IN)  :: Identifier !< A custom Identifier that describes the timer

    ! LOCAL VARIABLES
    INTEGER :: CurrentIndex
    REAL(8) :: EndTime, Duration

    !CALL CPU_TIME(EndTime)
    EndTime = MPI_WTIME()
    
    ! find index of identifier
    CurrentIndex = GetDurationIndex(Identifier)
    
    IF (CurrentIndex == 0) THEN
      CALL PrintWarningDuration(Identifier)
      RETURN
    ENDIF
    
    Duration = EndTime - StartTime(CurrentIndex)
    Durations(CurrentIndex) = Durations(CurrentIndex) + Duration
    TimeCount(CurrentIndex) = TimeCount(CurrentIndex) + 1
    
    IF (WITH_MEMORY_CONSUMPTION_MEASUREMENT) THEN
      CurrentMemoryConsumption = GetCurrentMemoryConsumption()
      TotalMemory(CurrentIndex) = TotalMemory(CurrentIndex) + (CurrentMemoryConsumption - StartMemory(CurrentIndex))
      StartMemory(CurrentIndex) = CurrentMemoryConsumption
    ENDIF
  END SUBROUTINE

  !
  !================================================================================================================================
  !
  SUBROUTINE CustomProfilingMemory(Identifier, NumberOfElements, TotalSize)
    ! PARAMETERS
    CHARACTER(LEN=*), INTENT(IN)  :: Identifier !< A custom Identifier that describes the timer
    INTEGER(INTG), INTENT(IN) :: NumberOfElements  !< number of elements
    INTEGER(INTG), INTENT(IN) :: TotalSize  !< number of bytes in total (e.g. by call to SIZEOF(<array>))
    INTEGER(LINTG) :: MemoryConsumption

    ! LOCAL VARIABLES
    INTEGER :: CurrentIndex
    INTEGER(INTG) :: SizePerElement  !< number of bytes of one element

    MemoryConsumption = INT8(TotalSize)
    
    IF (NumberOfElements > 0) THEN
      SizePerElement = TotalSize / NumberOfElements
    ELSE
      SizePerElement = 1
    ENDIF
    
    ! find index of identifier
    CurrentIndex = GetMemoryIndex(Identifier)
    
    ! If record with this identifier does not yet exist, create new
    IF (CurrentIndex == 0) THEN
      SizeMemory = SizeMemory + 1
      CurrentIndex = SizeMemory
      MemoryIdentifiers(CurrentIndex) = Identifier
      MemoryConsumptions(CurrentIndex) = 0
      SizesPerElement(CurrentIndex) = 0
      NumberOfObjects(CurrentIndex) = 0
    ENDIF
    
    ! Add value to record of memory consumption
    MemoryConsumptions(CurrentIndex) = MemoryConsumptions(CurrentIndex) + MemoryConsumption
    SizesPerElement(CurrentIndex) = SizePerElement
    NumberOfObjects(CurrentIndex) = NumberOfObjects(CurrentIndex) + 1

  END SUBROUTINE
  !
  !================================================================================================================================
  !
  SUBROUTINE CustomProfilingReset()
    ! reset is done by dropping the values in Duration* and Memory* fields
    SizeDuration = 0
    SizeMemory = 0
  END SUBROUTINE CustomProfilingReset
  !
  !================================================================================================================================
  !
  !>Checks if a given Identifier is contained in the records. If yes, return the matching index, if no returns 0
  FUNCTION GetDurationIndex(Identifier)
    CHARACTER(LEN=*), INTENT(IN) :: Identifier
    INTEGER :: GetDurationIndex    !< return value
    INTEGER :: I

    GetDurationIndex = 0

    DO I=SizeDuration,1,-1
      IF (TRIM(DurationIdentifiers(I)) == Identifier) THEN
      GetDurationIndex = I
      EXIT
      ENDIF
    END DO
  END FUNCTION GetDurationIndex
  !
  !================================================================================================================================
  !
  !>Checks if a given Identifier is contained in the records. If yes, return the matching index, if no returns 0
  FUNCTION GetMemoryIndex(Identifier)
    CHARACTER(LEN=*), INTENT(IN) :: Identifier
    INTEGER :: GetMemoryIndex    !< return value
    INTEGER :: I

    GetMemoryIndex = 0

    DO I=SizeMemory,1,-1
      IF (TRIM(MemoryIdentifiers(I)) == Identifier) THEN
      GetMemoryIndex = I
      EXIT
      ENDIF
    END DO
  END FUNCTION GetMemoryIndex
  !
  !================================================================================================================================
  !
  FUNCTION GetCurrentMemoryConsumption()
    INTEGER(LIntg) :: GetCurrentMemoryConsumption
    INTEGER(LIntg) :: VmSize, VmRSS, Shared, Text, Lib, Data, Dt, RssAnon
    INTEGER(Intg) :: Stat

    ! read from /proc/self/statm
    OPEN(UNIT=10, FILE="/proc/self/statm", ACTION="read", IOSTAT=stat)
    IF (STAT /= 0) THEN
      PRINT*, "Could not read memory consumption from /proc/self/statm."
      GetCurrentMemoryConsumption = 0
    ELSE
      READ(10,*, IOSTAT=stat) VmSize, VmRSS, Shared, Text, Lib, Data, Dt
      CLOSE(UNIT=10)

      ! VmRSS = RssAnon + Shared (all as number of pages)
      ! RssAnon is the resident set size of anonymous memory (real memory in RAM, not laid out in files)
      RssAnon = VmRSS - Shared
      GetCurrentMemoryConsumption = RssAnon
    ENDIF

  END FUNCTION GetCurrentMemoryConsumption
  !
  !================================================================================================================================
  !
  FUNCTION CustomProfilingGetInfo()
    CHARACTER(LEN=200000) :: CustomProfilingGetInfo
    CHARACTER(LEN=10000) :: Line
    INTEGER :: I, TotalNumberOfObjects = 0
    INTEGER(LINTG) :: TotalMemoryConsumption = 0

    ! collect timing
    CustomProfilingGetInfo = NEW_LINE('A') // '---------------------Profiling Info-----------------------' // NEW_LINE('A')
    CustomProfilingGetInfo = TRIM(CustomProfilingGetInfo) // "Timing" // NEW_LINE('A')
    DO I = 1,SizeDuration
      WRITE(Line,"(3A,F25.13,A,I10,2A)") "   ", (DurationIdentifiers(I)), ": ", Durations(I), ' s', &
        & TimeCount(I), 'x', NEW_LINE('A')
      CustomProfilingGetInfo = TRIM(CustomProfilingGetInfo) // TRIM(Line)
    ENDDO

    IF (WITH_MEMORY_CONSUMPTION_MEASUREMENT) THEN
      CustomProfilingGetInfo = TRIM(CustomProfilingGetInfo) // "Memory Consumption while Timing" // NEW_LINE('A')
      DO I = 1,SizeDuration
        WRITE(Line,"(3A,I17,A,I17,2A)") "   ", (DurationIdentifiers(I)), ": ", TotalMemory(I), ' B, avg.', &
          & TotalMemory(I)/TimeCount(I), ' B per profiling interval, ', NEW_LINE('A')
        CustomProfilingGetInfo = TRIM(CustomProfilingGetInfo) // TRIM(Line)
      ENDDO
    ENDIF
    
    ! Collect memory consumption
    CustomProfilingGetInfo = TRIM(CustomProfilingGetInfo) // "Static Memory Consumption" // NEW_LINE('A')
    TotalNumberOfObjects = 0
    TotalMemoryConsumption = 0
    DO I = 1,SizeMemory
      TotalMemoryConsumption = TotalMemoryConsumption + MemoryConsumptions(I)
      TotalNumberOfObjects = TotalNumberOfObjects + NumberOfObjects(I)
      WRITE(Line,"(3A,I17,A,I5,A,I8,A,I5,2A)") "   ", (MemoryIdentifiers(I)), ": ", MemoryConsumptions(I), ' B,', &
        & SizesPerElement(I), ' B,', MemoryConsumptions(I) / SizesPerElement(I), ' entries, ', &
        & NumberOfObjects(I), ' objects', NEW_LINE('A')
      CustomProfilingGetInfo = TRIM(CustomProfilingGetInfo) // TRIM(Line)
    ENDDO

    WRITE(Line,"(A,I17,A,I17,2A)") " Total: ", TotalMemoryConsumption, " B, ", TotalNumberOfObjects, " objects", NEW_LINE('A')

    CustomProfilingGetInfo = TRIM(CustomProfilingGetInfo) // TRIM(Line) // &
      & "----------------------------------------------------------" // NEW_LINE('A')
  END FUNCTION
  !
  !================================================================================================================================
  !
  FUNCTION CustomProfilingGetDuration(Identifier)
    CHARACTER(LEN=*), INTENT(IN)  :: Identifier !< A custom Identifier that describes the timer
    REAL(8) :: CustomProfilingGetDuration

    INTEGER :: CurrentIndex

    ! find index of identifier
    CurrentIndex = GetDurationIndex(Identifier)

    IF (CurrentIndex == 0) THEN
      CALL PrintWarningDuration(Identifier)
      CustomProfilingGetDuration = 0
      RETURN
    ENDIF

    CustomProfilingGetDuration = Durations(CurrentIndex)
  END FUNCTION
  !
  !================================================================================================================================
  !
  FUNCTION CustomProfilingGetMemory(Identifier)
    CHARACTER(LEN=*), INTENT(IN)  :: Identifier !< A custom Identifier that describes the timer
    INTEGER(LINTG) :: CustomProfilingGetMemory

    INTEGER :: CurrentIndex

    ! find index of identifier
    CurrentIndex = GetMemoryIndex(Identifier)

    IF (CurrentIndex == 0) THEN
      CALL PrintWarningMemory(Identifier)
      CustomProfilingGetMemory = 0
      RETURN
    ENDIF

    CustomProfilingGetMemory = MemoryConsumptions(CurrentIndex)
  END FUNCTION
  !
  !================================================================================================================================
  !
  FUNCTION CustomProfilingGetSizePerElement(Identifier)
    CHARACTER(LEN=*), INTENT(IN)  :: Identifier !< A custom Identifier that describes the timer
    INTEGER(INTG) :: CustomProfilingGetSizePerElement

    INTEGER :: CurrentIndex

    ! find index of identifier
    CurrentIndex = GetMemoryIndex(Identifier)

    IF (CurrentIndex == 0) THEN
      CALL PrintWarningMemory(Identifier)
      CustomProfilingGetSizePerElement = 0
      RETURN
    ENDIF

    CustomProfilingGetSizePerElement = SizesPerElement(CurrentIndex)
  END FUNCTION
  !
  !================================================================================================================================
  !
  FUNCTION CustomProfilingGetNumberObjects(Identifier)
    CHARACTER(LEN=*), INTENT(IN)  :: Identifier !< A custom Identifier that describes the timer
    INTEGER :: CustomProfilingGetNumberObjects

    INTEGER :: CurrentIndex

    ! find index of identifier
    CurrentIndex = GetMemoryIndex(Identifier)

    IF (CurrentIndex == 0) THEN
      CALL PrintWarningMemory(Identifier)
      CustomProfilingGetNumberObjects = 0
      RETURN
    ENDIF

    CustomProfilingGetNumberObjects = NumberOfObjects(CurrentIndex)
  END FUNCTION
  !
  !================================================================================================================================
  !
  SUBROUTINE PrintWarningDuration(Identifier)
    CHARACTER(LEN=*), INTENT(IN)  :: Identifier
    INTEGER :: I

    PRINT *, "Warning! Identifier for Duration, '", Identifier, "' does not exist. "
    PRINT *, SizeDuration, " existing records: "
    DO I = 1,SizeDuration
      PRINT *, "  '", DurationIdentifiers(I), "'"
    ENDDO
  END SUBROUTINE
  !
  !================================================================================================================================
  !
  SUBROUTINE PrintWarningMemory(Identifier)
    CHARACTER(LEN=*), INTENT(IN)  :: Identifier
    INTEGER :: I

    PRINT *, "Warning! Identifier for Memory, '", Identifier, "' does not exist. "
    PRINT *, SizeMemory, " existing records: "
    DO I = 1,SizeMemory
      PRINT *, "  '", MemoryIdentifiers(I), "'"
    ENDDO
  END SUBROUTINE

END MODULE Custom_Profiling
