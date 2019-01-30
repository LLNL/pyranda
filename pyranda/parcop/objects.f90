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
MODULE LES_objects
!===================================================================================================  
  USE iso_c_binding
  USE LES_patch,   ONLY :   patch_type                                    
  USE LES_comm,    ONLY :   comm_type
  USE LES_compact, ONLY :   compact_type
  USE LES_mesh,    ONLY :   mesh_type

  IMPLICIT NONE
  
  TYPE(patch_type),   TARGET  ::   patch_data(0:100,0:10) 
  TYPE(comm_type),    TARGET  ::    comm_data(0:100,0:10) 
  TYPE(compact_type), TARGET  :: compact_data(0:100,0:10) 
  TYPE(mesh_type),    TARGET  :: mesh_data(0:100,0:10) 

  TYPE(patch_type),   POINTER :: patch_ptr                    ! Points to an element of patch_data
  TYPE(comm_type),    POINTER :: comm_ptr                     ! Points to an element of comm_data
  TYPE(compact_type), POINTER :: compact_ptr                  ! Points to an element of compact_data
  TYPE(mesh_type),    POINTER :: mesh_ptr                     ! Points to an element of compact_data
  
  
CONTAINS
  
  SUBROUTINE setup_objects(patch,level,color,key,coordsys,nx,ny,nz,px,py,pz,x1,xn,y1,yn,z1,zn,bx1,bxn,by1,byn,bz1,bzn,simtime)
    IMPLICIT NONE
    INTEGER(c_int),               INTENT(IN) :: patch,level
    INTEGER(c_int),               INTENT(IN) :: color,key,coordsys,nx,ny,nz,px,py,pz
    REAL(c_double),               INTENT(IN) :: x1,xn,y1,yn,z1,zn
    CHARACTER(KIND=c_char,LEN=*), INTENT(IN) :: bx1,bxn,by1,byn,bz1,bzn
    REAL(c_double),               INTENT(IN) :: simtime
    CALL   patch_data(patch,level)%setup(color,key,coordsys,nx,ny,nz,px,py,pz,x1,xn,y1,yn,z1,zn,bx1,bxn,by1,byn,bz1,bzn,simtime)
    CALL    comm_data(patch,level)%setup(patch_data(patch,level))
    CALL compact_data(patch,level)%setup(patch_data(patch,level),comm_data(patch,level))
    

  END SUBROUTINE setup_objects

  SUBROUTINE setup_mesh_data(patch,level)
    IMPLICIT NONE 
    INTEGER(c_int), INTENT(IN) :: patch,level
    CALL mesh_data(patch,level)%setup(patch_data(patch,level),&
         comm_data(patch,level),&
         compact_data(patch,level))
  END SUBROUTINE setup_mesh_data

  SUBROUTINE setup_mesh_data_x3(patch,level,x1,x2,x3,mesh_per) !mesh_perX,mesh_perY,mesh_perZ)
    IMPLICIT NONE 
    INTEGER(c_int), INTENT(IN) :: patch,level
    REAL(kind=8), DIMENSION(:,:,:), INTENT(IN) :: x1,x2,x3
    LOGICAL, INTENT(IN) :: mesh_per
    !LOGICAL, INTENT(IN) :: mesh_perX,mesh_perY,mesh_perZ
    CALL mesh_data(patch,level)%setup(patch_data(patch,level),&
         comm_data(patch,level),&
         compact_data(patch,level),&
         xpy=x1,ypy=x2,zpy=x3,custom_per=mesh_per) !,&
         !custom_periodicX=mesh_perX,&
         !custom_periodicY=mesh_perY,&
         !custom_periodicZ=mesh_perZ)
  END SUBROUTINE setup_mesh_data_x3
  
  


  SUBROUTINE point_to_objects(patch,level)
    IMPLICIT NONE 
    INTEGER(c_int), INTENT(IN) :: patch,level
    !TYPE(patch_type),   INTENT(IN)  ::   patch_data 
    !TYPE(comm_type),    INTENT(IN)  ::    comm_data 
    !TYPE(compact_type), INTENT(IN)  :: compact_data 

    patch_ptr   => patch_data(patch,level)
    comm_ptr    => comm_data(patch,level)
    compact_ptr => compact_data(patch,level)
    mesh_ptr    => mesh_data(patch,level)
        
  END SUBROUTINE point_to_objects


   SUBROUTINE remove_objects(patch,level)
    IMPLICIT NONE
    INTEGER(c_int), INTENT(IN) :: patch,level
     CALL    mesh_data(patch,level)%remove()
     CALL compact_data(patch,level)%remove()
     CALL    comm_data(patch,level)%remove()
     CALL   patch_data(patch,level)%remove()
   END SUBROUTINE remove_objects


END MODULE LES_objects
