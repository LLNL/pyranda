!===================================================================================================
! Copyright (c) 2018, Lawrence Livemore National Security, LLC.
! Produced at the Lawrence Livermore National Laboratory.
!
! LLNL-CODE-749864
! This file is part of pyranda
! For details about use and distribution, please read: pyranda/LICENSE
!
! Written by: Britton J. Olson, olson45@llnl.gov
!===================================================================================================
!===================================================================================================
 MODULE LES_quest
!===================================================================================================  
  USE iso_c_binding
  USE MPI
  USE LES_comm
!#ifdef USE_QUEST
  !  USE axom_quest, ONLY : quest_initialize,quest_finalize,quest_distance
    USE axom_quest
!#endif
  
!----------------------------------------------------------------------------------------------------
  INTEGER, PARAMETER :: QUEST_SIGNED_DISTANCE_MODE = 1
  INTEGER, PARAMETER :: QUEST_IN_OUT_MODE          = 2
  INTEGER, PARAMETER :: QUEST_WATERTIGHT_GEOMETRY  = 0
  INTEGER, PARAMETER :: QUEST_TOPOLOGICALLY_PLANAR = 1
  INTEGER, PARAMETER :: QUEST_MAX_ELEMENTS  = 0
  INTEGER, PARAMETER :: QUEST_MAX_LEVELS    = 1
  INTEGER, PARAMETER :: QUEST_QUERY_MODE    = 2
  INTEGER, PARAMETER :: QUEST_VERBOSE       = 3
  INTEGER, PARAMETER :: QUEST_GEOMETRY_TYPE = 4
  INTEGER, PARAMETER :: QUEST_NUM_OPTIONS   = 5
  
  
!===================================================================================================
  CONTAINS
!===================================================================================================

    SUBROUTINE setup_quest(surfile)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: surfile
      INTEGER, dimension(0:QUEST_NUM_OPTIONS-1) :: quest_options
      LOGICAL(c_bool), PARAMETER :: useQuestdistance = .TRUE.
      INTEGER(c_int), PARAMETER :: nQdims = 3
      INTEGER(c_int), PARAMETER :: Q_max_elements = 25
      INTEGER(c_int), PARAMETER :: Q_max_levels = 20
      
!#ifdef USE_QUEST
         CALL quest_signed_distance_set_dimension( 3 )
         CALL quest_signed_distance_set_closed_surface( .True. )
         CALL quest_signed_distance_set_max_levels( QUEST_MAX_LEVELS )
         CALL quest_signed_distance_set_max_levels( QUEST_MAX_ELEMENTS )
         !CALL quest_signed_distance_set_verbose( .False. )
         CALL quest_signed_distance_set_shared_memory( .False. )

         !CALL quest_signed_distance_init_mpi(LES_comm_world, TRIM(surfile) )

!#endif
         
    END SUBROUTINE setup_quest
      
    REAL(KIND=c_double) FUNCTION distance(x,y,z) BIND(c)
      USE, INTRINSIC :: iso_c_binding, only : c_ptr, c_double, c_int
      REAL(c_double), INTENT(IN) :: x,y,z
       distance = 0.0
!#ifdef USE_QUEST
       distance = quest_signed_distance_evaluate(x,y,z)
!#endif
    END FUNCTION distance

    
   SUBROUTINE remove_quest()
     IMPLICIT NONE
     
!#ifdef USE_QUEST
     CALL quest_signed_distance_finalize()
!#endif
   END SUBROUTINE remove_quest


   
!===================================================================================================
 END MODULE LES_quest
!===================================================================================================
