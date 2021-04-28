!===================================================================================================
 MODULE LES_operators
!===================================================================================================
  USE iso_c_binding
  USE LES_objects, ONLY : patch_ptr,mesh_ptr,compact_ptr
  USE LES_compact_operators, ONLY : ddx=>d1x,ddy=>d1y,ddz=>d1z,d2x,d2y,d2z,d8x,d8y,d8z, &
       bppfx=>filterx,bppfy=>filtery,bppfz=>filterz
!  USE LES_explicit, ONLY : ddx=>ddx,ddy=>ddy,ddz=>ddz,d2x=>ddx,d2y=>ddy,d2z=>ddz, &
!       d8x=>dd4x,d8y=>dd4y,d8z=>dd4z,bppfx=>xFilter,bppfy=>yFilter,bppfz=>zFilter

  !USE LES_FFTs, ONLY : sfilterx,sfiltery,sfilterz
  !USE mapp_exosim_annotation, ONLY : exosim_annotation_begin,exosim_annotation_end
#ifdef fexlpool
  USE LES_ompsync
#else
  USE LES_ompsync, ONLY : sync_var
#endif
  
  INTERFACE div
   MODULE PROCEDURE divV, divT
  END INTERFACE
  INTERFACE grad
   MODULE PROCEDURE gradS, gradV, gradVc, gradTc
  END INTERFACE
  INTERFACE curl
   MODULE PROCEDURE curlV, curlT
  END INTERFACE
  INTERFACE cross
   MODULE PROCEDURE VcrossV, VcrossT
  END INTERFACE
  INTERFACE Laplacian
   MODULE PROCEDURE LapS, LapV
  END INTERFACE
  INTERFACE ring
   MODULE PROCEDURE ringS, ringV
  END INTERFACE
  
!===================================================================================================
  CONTAINS
!===================================================================================================

! DIVERGENCE OF A VECTOR ===========================================================================
   SUBROUTINE divV(fx,fy,fz,df,ax,ay,az)
    IMPLICIT NONE
    INTEGER(c_int), INTENT(IN) :: ax,ay,az 
    DOUBLE PRECISION, DIMENSION(ax,ay,az), INTENT(IN) :: fx,fy,fz
    DOUBLE PRECISION, DIMENSION(ax,ay,az), INTENT(OUT) :: df
    DOUBLE PRECISION, DIMENSION(ax,ay,az) :: fA,fB,fC,tmp
    INTEGER :: error
    integer :: ifunr,jfunr,kfunr,nfunr 

#ifdef omppool
    !$omp target enter data map(alloc:fA,fB,fC,tmp) nowait depend(inout:sync_var)
    !! $omp taskwait
#endif
!!$#ifdef fexlpool
!!$    ASSOCIATE(fA=>fexlPool_get(ax,ay,az),fB=>fexlPool_get(ax,ay,az),&
!!$         fC=>fexlPool_get(ax,ay,az),tmp=>fexlPool_get(ax,ay,az) )
!!$#endif      
!!$     SELECT CASE(patch_ptr%coordsys)
!!$     CASE(0) ! Cartesian
!!$      CALL ddx(fx,fA,patch_ptr%isymX)
!!$      CALL ddy(fy,fB,patch_ptr%isymY)
!!$      CALL ddz(fz,fC,patch_ptr%isymZ)
!!$
!!$!$omp target teams distribute parallel do collapse(3) nowait depend(inout:sync_var)
!!$      do kfunr=1,size(df,3)
!!$        do jfunr=1,size(df,2)
!!$          do ifunr=1,size(df,1)
!!$            df(ifunr,jfunr,kfunr)= fA(ifunr,jfunr,kfunr)+ fB(ifunr,jfunr,kfunr)+ fC(ifunr,jfunr,kfunr)
!!$           end do
!!$        end do
!!$      end do
!!$!$omp end target teams distribute parallel do
!!$
!!$     CASE(1)
!!$
!!$!$omp target teams distribute parallel do collapse(3)
!!$      do kfunr=1,size(tmp,3)
!!$        do jfunr=1,size(tmp,2)
!!$          do ifunr=1,size(tmp,1)
!!$            tmp(ifunr,jfunr,kfunr)= mesh_ptr%xgrid(ifunr,jfunr,kfunr)* fx(ifunr,jfunr,kfunr)! r*v_r
!!$          end do
!!$        end do
!!$      end do
!!$!$omp end target teams distribute parallel do
!!$
!!$      CALL ddx(tmp,fA,patch_ptr%isymX**2)
!!$      CALL ddy(fy,fB, patch_ptr%isymY)
!!$      CALL ddz(fz,fC, patch_ptr%isymZ)
!!$
!!$!$omp target teams distribute parallel do collapse(3)
!!$      do kfunr=1,size(df,3)
!!$        do jfunr=1,size(df,2)
!!$          do ifunr=1,size(df,1)
!!$            df(ifunr,jfunr,kfunr)=(fA(ifunr,jfunr,kfunr)+ fB(ifunr,jfunr,kfunr))/mesh_ptr%xgrid(ifunr,jfunr,kfunr)+ fC(ifunr,jfunr,kfunr)
!!$           end do
!!$        end do
!!$      end do
!!$!$omp end target teams distribute parallel do
!!$
!!$     CASE(2)
!!$      ! *** does sin(mesh_ptr%ygrid) ever generate patch_ptr%isymY factor? ***
!!$
!!$!$omp target teams distribute parallel do collapse(3)
!!$      do kfunr=1,size(tmp,3)
!!$        do jfunr=1,size(tmp,2)
!!$          do ifunr=1,size(tmp,1)
!!$            tmp(ifunr,jfunr,kfunr)= mesh_ptr%xgrid(ifunr,jfunr,kfunr)**2*fx(ifunr,jfunr,kfunr)! r** 2*v_r
!!$            df(ifunr,jfunr,kfunr)= fy(ifunr,jfunr,kfunr)* SIN(mesh_ptr%ygrid(ifunr,jfunr,kfunr)) ! v_theta*SIN(theta )
!!$          end do
!!$        end do
!!$      end do
!!$!$omp end target teams distribute parallel do
!!$
!!$      CALL ddx(tmp,fA,patch_ptr%isymX**3)
!!$      CALL ddy(df,fB, patch_ptr%isymY**2)  ! Assumes theta has some symmetry here... move outside derivative
!!$      CALL ddz(fz,fC, patch_ptr%isymZ)
!!$
!!$!$omp target teams distribute parallel do collapse(3)
!!$      do kfunr=1,size(fA,3)
!!$        do jfunr=1,size(fA,2)
!!$          do ifunr=1,size(fA,1)
!!$            df(ifunr,jfunr,kfunr)= fA(ifunr,jfunr,kfunr)/ mesh_ptr%xgrid(ifunr,jfunr,kfunr)**2+( fB(ifunr,jfunr,kfunr)+ fC(ifunr,jfunr,kfunr))/( mesh_ptr%xgrid(ifunr,jfunr,kfunr)* SIN(mesh_ptr%ygrid(ifunr,jfunr,kfunr)) )
!!$          end do
!!$        end do
!!$      end do
!!$!$omp end target teams distribute parallel do
!!$
!!$     CASE(3) ! General curvilinear
!!$
!!$!$omp target teams distribute parallel do collapse(3)
!!$      do kfunr=1,size(fA,3)
!!$        do jfunr=1,size(fA,2)
!!$          do ifunr=1,size(fA,1)
!!$            fA(ifunr,jfunr,kfunr)=(fx(ifunr,jfunr,kfunr)* mesh_ptr%dAdx(ifunr,jfunr,kfunr)+ fy(ifunr,jfunr,kfunr)* mesh_ptr%dAdy(ifunr,jfunr,kfunr)+ fz(ifunr,jfunr,kfunr)* mesh_ptr%dAdz(ifunr,jfunr,kfunr))*mesh_ptr%detxyz(ifunr,jfunr,kfunr)
!!$             fB(ifunr,jfunr,kfunr)=(fx(ifunr,jfunr,kfunr)* mesh_ptr%dBdx(ifunr,jfunr,kfunr)+ fy(ifunr,jfunr,kfunr)* mesh_ptr%dBdy(ifunr,jfunr,kfunr)+ fz(ifunr,jfunr,kfunr)* mesh_ptr%dBdz(ifunr,jfunr,kfunr))*mesh_ptr%detxyz(ifunr,jfunr,kfunr)
!!$             fC(ifunr,jfunr,kfunr)=(fx(ifunr,jfunr,kfunr)* mesh_ptr%dCdx(ifunr,jfunr,kfunr)+ fy(ifunr,jfunr,kfunr)* mesh_ptr%dCdy(ifunr,jfunr,kfunr)+ fz(ifunr,jfunr,kfunr)* mesh_ptr%dCdz(ifunr,jfunr,kfunr))*mesh_ptr%detxyz(ifunr,jfunr,kfunr)
!!$           end do
!!$        end do
!!$      end do
!!$!$omp end target teams distribute parallel do
!!$
!!$      CALL ddx(fA,df)
!!$      CALL ddy(fB,tmp)
!!$
!!$!$omp target teams distribute parallel do collapse(3)
!!$      do kfunr=1,size(df,3)
!!$        do jfunr=1,size(df,2)
!!$          do ifunr=1,size(df,1)
!!$            df(ifunr,jfunr,kfunr)= df(ifunr,jfunr,kfunr)+ tmp(ifunr,jfunr,kfunr)
!!$           end do
!!$        end do
!!$      end do
!!$!$omp end target teams distribute parallel do
!!$
!!$      CALL ddz(fC,tmp)
!!$
!!$!$omp target teams distribute parallel do collapse(3)
!!$      do kfunr=1,size(df,3)
!!$        do jfunr=1,size(df,2)
!!$          do ifunr=1,size(df,1)
!!$            df(ifunr,jfunr,kfunr)=(df(ifunr,jfunr,kfunr)+ tmp(ifunr,jfunr,kfunr))/mesh_ptr%detxyz(ifunr,jfunr,kfunr)
!!$           end do
!!$        end do
!!$      end do
!!$!$omp end target teams distribute parallel do
!!$
!!$   END SELECT

#ifdef fexlpool
 END ASSOCIATE
 error = fexlPool_free(4)
#endif
#ifdef omppool
 !$omp target exit data map(release:fA,fB,fC,tmp) nowait depend(inout:sync_var)
 !$! omp taskwait
#endif
   END SUBROUTINE divV
 
! DIVERGENCE OF A TENSOR ===========================================================================
   SUBROUTINE divT(fxx,fxy,fxz,fyx,fyy,fyz,fzx,fzy,fzz,dfx,dfy,dfz)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: fxx,fxy,fxz,fyx,fyy,fyz,fzx,fzy,fzz
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: dfx,dfy,dfz
    DOUBLE PRECISION, DIMENSION(SIZE(fxx,1),SIZE(fxx,2),SIZE(fxx,3)) :: fA,fB,fC,tmp
    integer :: ifunr,jfunr,kfunr,nfunr 

    !$omp target data map(alloc:fA,fB,fC,tmp)
     SELECT CASE(patch_ptr%coordsys)
     CASE(0) ! Cartesian
      CALL ddx(fxx,fA,patch_ptr%isymX**2)
      CALL ddy(fyx,fB,patch_ptr%isymY)
      CALL ddz(fzx,fC,patch_ptr%isymZ)

!$omp target teams distribute parallel do collapse(3)
      do kfunr=1,size(dfx,3)
        do jfunr=1,size(dfx,2)
          do ifunr=1,size(dfx,1)
            dfx(ifunr,jfunr,kfunr)= fA(ifunr,jfunr,kfunr)+ fB(ifunr,jfunr,kfunr)+ fC(ifunr,jfunr,kfunr)
           end do
        end do
      end do
!$omp end target teams distribute parallel do

      CALL ddx(fxy,fA,patch_ptr%isymX)
      CALL ddy(fyy,fB,patch_ptr%isymY**2)
      CALL ddz(fzy,fC,patch_ptr%isymZ)

!$omp target teams distribute parallel do collapse(3)
      do kfunr=1,size(dfy,3)
        do jfunr=1,size(dfy,2)
          do ifunr=1,size(dfy,1)
            dfy(ifunr,jfunr,kfunr)= fA(ifunr,jfunr,kfunr)+ fB(ifunr,jfunr,kfunr)+ fC(ifunr,jfunr,kfunr)
           end do
        end do
      end do
!$omp end target teams distribute parallel do

      CALL ddx(fxz,fA,patch_ptr%isymX)
      CALL ddy(fyz,fB,patch_ptr%isymY)
      CALL ddz(fzz,fC,patch_ptr%isymZ**2)

!$omp target teams distribute parallel do collapse(3)
      do kfunr=1,size(dfz,3)
        do jfunr=1,size(dfz,2)
          do ifunr=1,size(dfz,1)
            dfz(ifunr,jfunr,kfunr)= fA(ifunr,jfunr,kfunr)+ fB(ifunr,jfunr,kfunr)+ fC(ifunr,jfunr,kfunr)
           end do
        end do
      end do
!$omp end target teams distribute parallel do

     CASE(1)

!$omp target teams distribute parallel do collapse(3)
      do kfunr=1,size(tmp,3)
        do jfunr=1,size(tmp,2)
          do ifunr=1,size(tmp,1)
            tmp(ifunr,jfunr,kfunr)= mesh_ptr%xgrid(ifunr,jfunr,kfunr)* fxx(ifunr,jfunr,kfunr)
           end do
        end do
      end do
!$omp end target teams distribute parallel do

      CALL ddx(tmp,fA,patch_ptr%isymX**3)
      CALL ddy(fyx,fB,patch_ptr%isymY)
      CALL ddz(fzx,fC,patch_ptr%isymZ)

!$omp target teams distribute parallel do collapse(3)
      do kfunr=1,size(dfx,3)
        do jfunr=1,size(dfx,2)
          do ifunr=1,size(dfx,1)
            dfx(ifunr,jfunr,kfunr)=(fA(ifunr,jfunr,kfunr)+ fB(ifunr,jfunr,kfunr)- fyy(ifunr,jfunr,kfunr))/mesh_ptr%xgrid(ifunr,jfunr,kfunr)+ fC(ifunr,jfunr,kfunr)
             tmp(ifunr,jfunr,kfunr)= mesh_ptr%xgrid(ifunr,jfunr,kfunr)**2*fxy(ifunr,jfunr,kfunr)
           end do
        end do
      end do
!$omp end target teams distribute parallel do

      CALL ddx(tmp,fA,patch_ptr%isymX**3)
      CALL ddy(fyy,fB,patch_ptr%isymY**2)
      CALL ddz(fzy,fC,patch_ptr%isymZ)

!$omp target teams distribute parallel do collapse(3)
      do kfunr=1,size(dfy,3)
        do jfunr=1,size(dfy,2)
          do ifunr=1,size(dfy,1)
            dfy(ifunr,jfunr,kfunr)= fA(ifunr,jfunr,kfunr)/ mesh_ptr%xgrid(ifunr,jfunr,kfunr)**2+( fB(ifunr,jfunr,kfunr)+ fyx(ifunr,jfunr,kfunr)- fxy(ifunr,jfunr,kfunr))/mesh_ptr%xgrid(ifunr,jfunr,kfunr)+ fC(ifunr,jfunr,kfunr)
             tmp(ifunr,jfunr,kfunr)= mesh_ptr%xgrid(ifunr,jfunr,kfunr)* fxz(ifunr,jfunr,kfunr)
           end do
        end do
      end do
!$omp end target teams distribute parallel do

      CALL ddx(tmp,fA,patch_ptr%isymX**2)
      CALL ddy(fyz,fB,patch_ptr%isymY)
      CALL ddz(fzz,fC,patch_ptr%isymZ**2)

!$omp target teams distribute parallel do collapse(3)
      do kfunr=1,size(dfz,3)
        do jfunr=1,size(dfz,2)
          do ifunr=1,size(dfz,1)
            dfz(ifunr,jfunr,kfunr)=(fA(ifunr,jfunr,kfunr)+ fB(ifunr,jfunr,kfunr))/mesh_ptr%xgrid(ifunr,jfunr,kfunr)+ fC(ifunr,jfunr,kfunr)
           end do
        end do
      end do
!$omp end target teams distribute parallel do

    CASE(2)
    ! *** does sin(mesh_ptr%ygrid) ever generate patch_ptr%isymY factor? ***
      tmp = mesh_ptr%xgrid**2*fxx
      CALL ddx(tmp,fA,patch_ptr%isymX**4)
      tmp = fyx*SIN(mesh_ptr%ygrid)
      CALL ddy(tmp,fB,patch_ptr%isymY**2)   !?
      CALL ddz(fzx,fC,patch_ptr%isymZ)
      dfx = fA/mesh_ptr%xgrid**2 + (fB+fC)/(mesh_ptr%xgrid*SIN(mesh_ptr%ygrid))-(fyy+fzz)/mesh_ptr%xgrid
      tmp = mesh_ptr%xgrid**3*fxy
      CALL ddx(tmp,fA,patch_ptr%isymX**4)
      tmp = fyy*SIN(mesh_ptr%ygrid)
      CALL ddy(tmp,fB,patch_ptr%isymY**3)   !?
      CALL ddz(fzy,fC,patch_ptr%isymZ)
      dfy = fA/mesh_ptr%xgrid**3 + (fB+fC)/(mesh_ptr%xgrid*SIN(mesh_ptr%ygrid)) + (fyx-fxy-fzz/TAN(mesh_ptr%ygrid))/mesh_ptr%xgrid
      tmp = mesh_ptr%xgrid**3*fxz
      CALL ddx(tmp,fA,patch_ptr%isymX**4)
      tmp = fyz*SIN(mesh_ptr%ygrid)
      CALL ddy(tmp,fB,patch_ptr%isymY**2)   !?
      CALL ddz(fzz,fC,patch_ptr%isymZ**2)
      dfz = fA/mesh_ptr%xgrid**3 + (fB+fC)/(mesh_ptr%xgrid*SIN(mesh_ptr%ygrid)) + (fzx-fxz+fzy/TAN(mesh_ptr%ygrid))/mesh_ptr%xgrid
     CASE(3)
      CALL divV(fxx,fxy,fxz,dfx,patch_ptr%ax,patch_ptr%ay,patch_ptr%az)
      CALL divV(fyx,fyy,fyz,dfy,patch_ptr%ax,patch_ptr%ay,patch_ptr%az)
      CALL divV(fzx,fzy,fzz,dfz,patch_ptr%ax,patch_ptr%ay,patch_ptr%az)
     END SELECT
     !$omp end target data
   END SUBROUTINE divT
 
! GRADIENT OF A SCALAR =============================================================================
   SUBROUTINE gradS(f,dfdx,dfdy,dfdz)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: f
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: dfdx,dfdy,dfdz
    DOUBLE PRECISION, DIMENSION(SIZE(f,1),SIZE(f,2),SIZE(f,3)) :: dfdA,dfdB,dfdC
    integer :: ifunr,jfunr,kfunr,nfunr 

    !$omp target data map(alloc:dfdA,dfdB,dfdC)
     CALL ddx(f,dfdx)
     CALL ddy(f,dfdy)
     CALL ddz(f,dfdz)
     SELECT CASE(patch_ptr%coordsys)
     CASE(1) ! Cylindrical

!$omp target teams distribute parallel do collapse(3)
      do kfunr=1,size(dfdy,3)
        do jfunr=1,size(dfdy,2)
          do ifunr=1,size(dfdy,1)
            dfdy(ifunr,jfunr,kfunr)= dfdy(ifunr,jfunr,kfunr)/ mesh_ptr%xgrid(ifunr,jfunr,kfunr)
           end do
        end do
      end do
!$omp end target teams distribute parallel do

     CASE(2) ! Spherical
      dfdy = dfdy/mesh_ptr%xgrid
      dfdz = dfdz/(mesh_ptr%xgrid*SIN(mesh_ptr%ygrid))
     CASE(3) ! General
      dfdA = dfdx
      dfdB = dfdy
      dfdC = dfdz
      dfdx = dfdA*mesh_ptr%dAdx + dfdB*mesh_ptr%dBdx + dfdC*mesh_ptr%dCdx
      dfdy = dfdA*mesh_ptr%dAdy + dfdB*mesh_ptr%dBdy + dfdC*mesh_ptr%dCdy
      dfdz = dfdA*mesh_ptr%dAdz + dfdB*mesh_ptr%dBdz + dfdC*mesh_ptr%dCdz
     END SELECT
     !$omp end target data
   END SUBROUTINE gradS
 
! GRADIENT OF A VECTOR COMPONENT =============================================================================
   SUBROUTINE gradVc(f,dfdx,dfdy,dfdz,vc)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: f
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: dfdx,dfdy,dfdz
    INTEGER, INTENT(IN) :: vc
    DOUBLE PRECISION, DIMENSION(SIZE(f,1),SIZE(f,2),SIZE(f,3)) :: dfdA,dfdB,dfdC
    integer :: ifunr,jfunr,kfunr,nfunr 

    !$omp target data map(alloc:dfdA,dfdB,dfdC)
     SELECT CASE( vc )
     CASE( 1 )
       CALL ddx(f,dfdx,patch_ptr%isymX); CALL ddy(f,dfdy);                 CALL ddz(f,dfdz)
     CASE( 2 )
       CALL ddx(f,dfdx);                 CALL ddy(f,dfdy,patch_ptr%isymY); CALL ddz(f,dfdz)
     CASE( 3 )
       CALL ddx(f,dfdx);                 CALL ddy(f,dfdy);                 CALL ddz(f,dfdz,patch_ptr%isymZ)
     CASE DEFAULT
       CALL ddx(f,dfdx);                 CALL ddy(f,dfdy);                 CALL ddz(f,dfdz)
     END SELECT
     SELECT CASE(patch_ptr%coordsys)
     CASE(1) ! Cylindrical

!$omp target teams distribute parallel do collapse(3)
      do kfunr=1,size(dfdy,3)
        do jfunr=1,size(dfdy,2)
          do ifunr=1,size(dfdy,1)
            dfdy(ifunr,jfunr,kfunr)= dfdy(ifunr,jfunr,kfunr)/ mesh_ptr%xgrid(ifunr,jfunr,kfunr)
           end do
        end do
      end do
!$omp end target teams distribute parallel do

     CASE(2) ! Spherical
      dfdy = dfdy/mesh_ptr%xgrid
      dfdz = dfdz/(mesh_ptr%xgrid*SIN(mesh_ptr%ygrid))
     CASE(3) ! General
      dfdA = dfdx
      dfdB = dfdy
      dfdC = dfdz
      dfdx = dfdA*mesh_ptr%dAdx + dfdB*mesh_ptr%dBdx + dfdC*mesh_ptr%dCdx
      dfdy = dfdA*mesh_ptr%dAdy + dfdB*mesh_ptr%dBdy + dfdC*mesh_ptr%dCdy
      dfdz = dfdA*mesh_ptr%dAdz + dfdB*mesh_ptr%dBdz + dfdC*mesh_ptr%dCdz
     END SELECT
     !$omp end target data
   END SUBROUTINE gradVc
 
! GRADIENT OF A TENSOR COMPONENT =============================================================================
   SUBROUTINE gradTc(f,dfdx,dfdy,dfdz,tc)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: f
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: dfdx,dfdy,dfdz
    CHARACTER(LEN=2), INTENT(IN) :: tc
    DOUBLE PRECISION, DIMENSION(SIZE(f,1),SIZE(f,2),SIZE(f,3)) :: dfdA,dfdB,dfdC
    integer :: ifunr,jfunr,kfunr,nfunr 

     SELECT CASE( tc )
     CASE( 'xx', 'yy', 'zz' )
       CALL ddx(f,dfdx);                 CALL ddy(f,dfdy);                 CALL ddz(f,dfdz)
     CASE( 'xy', 'yx' )
       CALL ddx(f,dfdx,patch_ptr%isymX); CALL ddy(f,dfdy,patch_ptr%isymY); CALL ddz(f,dfdz)
     CASE( 'xz', 'zx' )
       CALL ddx(f,dfdx,patch_ptr%isymX); CALL ddy(f,dfdy);                 CALL ddz(f,dfdz,patch_ptr%isymZ)
     CASE( 'yz', 'zy' )
       CALL ddx(f,dfdx);                 CALL ddy(f,dfdy,patch_ptr%isymY); CALL ddz(f,dfdz,patch_ptr%isymZ)
     CASE DEFAULT
       CALL ddx(f,dfdx);                 CALL ddy(f,dfdy);                 CALL ddz(f,dfdz)
     END SELECT
     SELECT CASE(patch_ptr%coordsys) !-----switches for non-Cartesian to be implemented---
     CASE(1) ! Cylindrical

!$omp target teams distribute parallel do collapse(3)
      do kfunr=1,size(dfdy,3)
        do jfunr=1,size(dfdy,2)
          do ifunr=1,size(dfdy,1)
            dfdy(ifunr,jfunr,kfunr)= dfdy(ifunr,jfunr,kfunr)/ mesh_ptr%xgrid(ifunr,jfunr,kfunr)
           end do
        end do
      end do
!$omp end target teams distribute parallel do

     CASE(2) ! Spherical
      dfdy = dfdy/mesh_ptr%xgrid
      dfdz = dfdz/(mesh_ptr%xgrid*SIN(mesh_ptr%ygrid))
     CASE(3) ! General
      dfdA = dfdx
      dfdB = dfdy
      dfdC = dfdz
      dfdx = dfdA*mesh_ptr%dAdx + dfdB*mesh_ptr%dBdx + dfdC*mesh_ptr%dCdx
      dfdy = dfdA*mesh_ptr%dAdy + dfdB*mesh_ptr%dBdy + dfdC*mesh_ptr%dCdy
      dfdz = dfdA*mesh_ptr%dAdz + dfdB*mesh_ptr%dBdz + dfdC*mesh_ptr%dCdz
     END SELECT
   END SUBROUTINE gradTc
 
! GRADIENT OF A VECTOR =============================================================================
   SUBROUTINE gradV(fx,fy,fz,dfxx,dfxy,dfxz,dfyx,dfyy,dfyz,dfzx,dfzy,dfzz)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: fx,fy,fz
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: dfxx,dfxy,dfxz,dfyx,dfyy,dfyz,dfzx,dfzy,dfzz
    DOUBLE PRECISION, DIMENSION(SIZE(fx,1),SIZE(fx,2),SIZE(fx,3)) :: dfdA,dfdB,dfdC
    integer :: ifunr,jfunr,kfunr,nfunr 

     CALL ddx(fx,dfxx,patch_ptr%isymX)
     CALL ddx(fy,dfxy)
     CALL ddx(fz,dfxz)
     CALL ddy(fx,dfyx)
     CALL ddy(fy,dfyy,patch_ptr%isymY)
     CALL ddy(fz,dfyz)
     CALL ddz(fx,dfzx)
     CALL ddz(fy,dfzy)
     CALL ddz(fz,dfzz,patch_ptr%isymZ)
     SELECT CASE(patch_ptr%coordsys)
     CASE(1) ! Cylindrical

!$omp target teams distribute parallel do collapse(3)
      do kfunr=1,size(dfyx,3)
        do jfunr=1,size(dfyx,2)
          do ifunr=1,size(dfyx,1)
            dfyx(ifunr,jfunr,kfunr)=(dfyx(ifunr,jfunr,kfunr)- fy(ifunr,jfunr,kfunr))/mesh_ptr%xgrid(ifunr,jfunr,kfunr)
             dfyy(ifunr,jfunr,kfunr)=(dfyy(ifunr,jfunr,kfunr)+ fx(ifunr,jfunr,kfunr))/mesh_ptr%xgrid(ifunr,jfunr,kfunr)
             dfyz(ifunr,jfunr,kfunr)= dfyz(ifunr,jfunr,kfunr)/ mesh_ptr%xgrid(ifunr,jfunr,kfunr)
           end do
        end do
      end do
!$omp end target teams distribute parallel do

     CASE(2) ! Spherical
      dfyx = (dfyx-fy)/mesh_ptr%xgrid
      dfyy = (dfyy+fx)/mesh_ptr%xgrid
      dfyz = dfyz/mesh_ptr%xgrid
      dfzx = (dfzx/SIN(mesh_ptr%ygrid) - fz)/mesh_ptr%xgrid
      dfzy = (dfzy/SIN(mesh_ptr%ygrid) - fz/TAN(mesh_ptr%ygrid))/mesh_ptr%xgrid
      dfzz = (dfzz/SIN(mesh_ptr%ygrid) + fx + fy/TAN(mesh_ptr%ygrid))/mesh_ptr%xgrid
     CASE(3) ! General
      dfdA = dfxx; dfdB = dfyx; dfdC = dfzx
      dfxx = dfdA*mesh_ptr%dAdx + dfdB*mesh_ptr%dBdx + dfdC*mesh_ptr%dCdx
      dfyx = dfdA*mesh_ptr%dAdy + dfdB*mesh_ptr%dBdy + dfdC*mesh_ptr%dCdy
      dfzx = dfdA*mesh_ptr%dAdz + dfdB*mesh_ptr%dBdz + dfdC*mesh_ptr%dCdz
      dfdA = dfxy; dfdB = dfyy; dfdC = dfzy
      dfxy = dfdA*mesh_ptr%dAdx + dfdB*mesh_ptr%dBdx + dfdC*mesh_ptr%dCdx
      dfyy = dfdA*mesh_ptr%dAdy + dfdB*mesh_ptr%dBdy + dfdC*mesh_ptr%dCdy
      dfzy = dfdA*mesh_ptr%dAdz + dfdB*mesh_ptr%dBdz + dfdC*mesh_ptr%dCdz
      dfdA = dfxz; dfdB = dfyz; dfdC = dfzz
      dfxz = dfdA*mesh_ptr%dAdx + dfdB*mesh_ptr%dBdx + dfdC*mesh_ptr%dCdx
      dfyz = dfdA*mesh_ptr%dAdy + dfdB*mesh_ptr%dBdy + dfdC*mesh_ptr%dCdy
      dfzz = dfdA*mesh_ptr%dAdz + dfdB*mesh_ptr%dBdz + dfdC*mesh_ptr%dCdz
     END SELECT
   END SUBROUTINE gradV
 
! CURL OF A VECTOR =================================================================================
!   SUBROUTINE curlV(fx,fy,fz,cx,cy,cz)
!    IMPLICIT NONE
!    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: fx,fy,fz
!    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: cx,cy,cz
!    DOUBLE PRECISION, DIMENSION(SIZE(fx,1),SIZE(fx,2),SIZE(fx,3)) :: tx,ty,tz,tmp
!     SELECT CASE(patch_ptr%coordsys)
!     CASE(0) ! Cartesian
!      CALL ddy(fz,cx)
!      CALL ddz(fy,tx)
!      cx = cx - tx
!      CALL ddz(fx,cy)
!      CALL ddx(fz,ty)
!      cy = cy - ty
!      CALL ddx(fy,cz)
!      CALL ddy(fx,tz)
!      cz = cz - tz
!     CASE(1) ! Cylindrical
!      CALL ddy(fz,tmp)
!      cx = tmp/mesh_ptr%xgrid
!      CALL ddz(fy,tx)
!      cx = cx - tx
!      CALL ddz(fx,cy)
!      CALL ddx(fz,ty)
!      cy = cy - ty
!      cz = mesh_ptr%xgrid*fy
!      CALL ddx(cz,tmp,patch_ptr%isymX)    !?
!      cz = tmp/mesh_ptr%xgrid
!      CALL ddy(fx,tmp)
!      tz = tmp/mesh_ptr%xgrid
!      cz = cz - tz
!     CASE(2) ! Spherical
!      cx = fz*SIN(mesh_ptr%ygrid)
!      CALL ddy(cx,tmp,patch_ptr%isymY)    !?
!      cx = tmp/(mesh_ptr%xgrid*SIN(mesh_ptr%ygrid))
!      CALL ddz(fy,tmp)
!      tx = tmp/(mesh_ptr%xgrid*SIN(mesh_ptr%ygrid))
!      cx = cx - tx
!      CALL ddz(fx,tmp)
!      cy = tmp/(mesh_ptr%xgrid*SIN(mesh_ptr%ygrid))
!      ty = mesh_ptr%xgrid*fz
!      CALL ddx(ty,tmp,patch_ptr%isymX)    !?
!      ty = tmp/mesh_ptr%xgrid
!      cy = cy - ty
!      cz = mesh_ptr%xgrid*fy
!      CALL ddx(cz,tmp,patch_ptr%isymX)
!      cz = tmp/mesh_ptr%xgrid
!      CALL ddy(fx,tmp)
!      tz = tmp/mesh_ptr%xgrid
!      cz = cz - tz
!     CASE(3) ! Hack to avoid NANs in cylinder_curvy.csv file
!      cx = 0.0D0
!      cy = 0.0D0
!      cz = 0.0D0      
!     END SELECT
!   END SUBROUTINE curlV
 
! CURL OF A SYMMETRIC (e.g., for velocity, vsym=1) OR ANTISYMMETRIC (e.g., for magnetic field, vsym=-1) VECTOR
   SUBROUTINE curlV(vsym,fx,fy,fz,cx,cy,cz)
    IMPLICIT NONE
    INTEGER(c_int), INTENT(IN) :: vsym
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: fx,fy,fz
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: cx,cy,cz
    DOUBLE PRECISION, DIMENSION(SIZE(fx,1),SIZE(fx,2),SIZE(fx,3)) :: tx,ty,tz,tmp
    INTEGER(c_int) :: isymX,isymY,isymZ
     isymX = -vsym*patch_ptr%isymX
     isymY = -vsym*patch_ptr%isymY    
     isymZ = -vsym*patch_ptr%isymZ        
     SELECT CASE(patch_ptr%coordsys)
     CASE(0) ! Cartesian
      CALL ddy(fz,cx,isymY)
      CALL ddz(fy,tx,isymZ)
      cx = cx - tx
      CALL ddz(fx,cy,isymZ)
      CALL ddx(fz,ty,isymX)
      cy = cy - ty
      CALL ddx(fy,cz,isymX)
      CALL ddy(fx,tz,isymY)
      cz = cz - tz
! NEED TO CHECK THESE SYMMETRY SWITCHES!------------------------------------------------------------      
     CASE(1) ! Cylindrical
      CALL ddy(fz,tmp,isymY)                      !?  
      cx = tmp/mesh_ptr%xgrid
      CALL ddz(fy,tx,isymZ)                       !?
      cx = cx - tx
      CALL ddz(fx,cy,isymZ)                       !?
      CALL ddx(fz,ty,isymX)                       !?
      cy = cy - ty
      cz = mesh_ptr%xgrid*fy
      CALL ddx(cz,tmp,isymX**2)                   !??
      cz = tmp/mesh_ptr%xgrid
      CALL ddy(fx,tmp,isymY)                      !?
      tz = tmp/mesh_ptr%xgrid
      cz = cz - tz
     CASE(2) ! Spherical
      ! *** Does sin(mesh_ptr%ygrid) generate isymY factor? ***
      cx = fz*SIN(mesh_ptr%ygrid)
      CALL ddy(cx,tmp,isymY**2)                   !??
      cx = tmp/(mesh_ptr%xgrid*SIN(mesh_ptr%ygrid))
      CALL ddz(fy,tmp,isymZ)                      !?
      tx = tmp/(mesh_ptr%xgrid*SIN(mesh_ptr%ygrid))
      cx = cx - tx
      CALL ddz(fx,tmp,isymZ)                      !?
      cy = tmp/(mesh_ptr%xgrid*SIN(mesh_ptr%ygrid))
      ty = mesh_ptr%xgrid*fz
      CALL ddx(ty,tmp,isymX**2)                   !??
      ty = tmp/mesh_ptr%xgrid
      cy = cy - ty
      cz = mesh_ptr%xgrid*fy
      CALL ddx(cz,tmp,isymX**2)                   !??
      cz = tmp/mesh_ptr%xgrid
      CALL ddy(fx,tmp,isymY)                      !?
      tz = tmp/mesh_ptr%xgrid
      cz = cz - tz
!---------------------------------------------------------------------------------------------------      
     CASE(3) ! Hack to avoid NANs in cylinder_curvy.csv file
      cx = 0.0D0
      cy = 0.0D0
      cz = 0.0D0      
     END SELECT
   END SUBROUTINE curlV
 
! CURL OF A TENSOR =================================================================================
   ! This needs symmetry switch for symmetric/anti-symmetric tensors.
   SUBROUTINE curlT(fxx,fxy,fxz,fyx,fyy,fyz,fzx,fzy,fzz,ttxx,ttxy,ttxz,ttyx,ttyy,ttyz,ttzx,ttzy,ttzz)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: fxx,fxy,fxz,fyx,fyy,fyz,fzx,fzy,fzz
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: ttxx,ttxy,ttxz,ttyx,ttyy,ttyz,ttzx,ttzy,ttzz
    DOUBLE PRECISION, DIMENSION(SIZE(fxx,1),SIZE(fxx,2),SIZE(fxx,3)) :: tmp!tx,ty,tz
     SELECT CASE(patch_ptr%coordsys)
     CASE(0) ! Cartesian
      CALL ddy(fzx,ttxx);                  CALL ddz(fyx,tmp);                 ttxx = ttxx - tmp
      CALL ddy(fzy,ttxy,patch_ptr%isymY);  CALL ddz(fyy,tmp);                 ttxy = ttxy - tmp
      CALL ddy(fzz,ttxz);                  CALL ddz(fyz,tmp,patch_ptr%isymZ); ttxz = ttxz - tmp

      CALL ddz(fxx,ttyx);                  CALL ddx(fzx,tmp,patch_ptr%isymX); ttyx = ttyx - tmp
      CALL ddz(fxy,ttyy);                  CALL ddx(fzy,tmp);                 ttyy = ttyy - tmp
      CALL ddz(fxz,ttyz,patch_ptr%isymZ);  CALL ddx(fzz,tmp);                 ttyz = ttyz - tmp

      CALL ddx(fyx,ttzx,patch_ptr%isymX);  CALL ddy(fxx,tmp);                 ttzx = ttzx - tmp
      CALL ddx(fyy,ttzy);                  CALL ddy(fxy,tmp,patch_ptr%isymY); ttzy = ttzy - tmp
      CALL ddx(fyz,ttzz);                  CALL ddy(fxz,tmp);                 ttzz = ttzz - tmp
! SYMMETRY SWITCHES FOR NON-CARTESIAN NEED TO BE FIXED----------------------------------------------      
     CASE(1) ! Cylindrical
      CALL ddy(fzx,ttxx);   CALL ddz(fyx,tmp); ttxx = (ttxx - fzy)/mesh_ptr%xgrid - tmp
      CALL ddy(fzy,ttxy);   CALL ddz(fyy,tmp); ttxy = (ttxy + fzx)/mesh_ptr%xgrid - tmp
      CALL ddy(fzz,ttxz);   CALL ddz(fyz,tmp); ttxz =  ttxz       /mesh_ptr%xgrid - tmp

      CALL ddz(fxx,ttyx);   CALL ddx(fzx,tmp); ttyx = ttyx - tmp
      CALL ddz(fxy,ttyy);   CALL ddx(fzy,tmp); ttyy = ttyy - tmp
      CALL ddz(fxz,ttyz);   CALL ddx(fzz,tmp); ttyz = ttyz - tmp

      CALL ddx(fyx,ttzx);   CALL ddy(fxx,tmp); ttzx = ttzx - (tmp - fyx - fxy)/mesh_ptr%xgrid
      CALL ddx(fyy,ttzy);   CALL ddy(fxy,tmp); ttzy = ttzy - (tmp - fyy + fxx)/mesh_ptr%xgrid
      CALL ddx(fyz,ttzz);   CALL ddy(fxz,tmp); ttzz = ttzz - (tmp - fyz      )/mesh_ptr%xgrid
     CASE(2) ! Spherical
      CALL ddy(fzx,ttxx);   CALL ddz(fyx,tmp); ttxx = (ttxx - tmp/SIN(mesh_ptr%ygrid) - fzy + fyz + fzx  /TAN(mesh_ptr%ygrid))/mesh_ptr%xgrid
      CALL ddy(fzy,ttxy);   CALL ddz(fyy,tmp); ttxy = (ttxy - tmp/SIN(mesh_ptr%ygrid) + fzx + (fyz + fzy)/TAN(mesh_ptr%ygrid))/mesh_ptr%xgrid
      CALL ddy(fzz,ttxz);   CALL ddz(fyz,tmp); ttxz = (ttxz - tmp/SIN(mesh_ptr%ygrid) - fyx + (fzz - fyy)/TAN(mesh_ptr%ygrid))/mesh_ptr%xgrid
 
      CALL ddz(fxx,ttyx);   CALL ddx(fzx,tmp); ttyx = (ttyx/SIN(mesh_ptr%ygrid) - tmp - fzx - fxz                          )/mesh_ptr%xgrid
      CALL ddz(fxy,ttyy);   CALL ddx(fzy,tmp); ttyy = (ttyy/SIN(mesh_ptr%ygrid) - tmp - fzy       - fxz/TAN(mesh_ptr%ygrid))/mesh_ptr%xgrid
      CALL ddz(fxz,ttyz);   CALL ddx(fzz,tmp); ttyz = (ttyz/SIN(mesh_ptr%ygrid) - tmp + fxx - fzz + fxy/TAN(mesh_ptr%ygrid))/mesh_ptr%xgrid

      CALL ddx(fyx,ttzx);   CALL ddy(fxx,tmp); ttzx = ttzx - (tmp - fyx - fxy)/mesh_ptr%xgrid
      CALL ddx(fyy,ttzy);   CALL ddy(fxy,tmp); ttzy = ttzy - (tmp - fyy + fxx)/mesh_ptr%xgrid
      CALL ddx(fyz,ttzz);   CALL ddy(fxz,tmp); ttzz = ttzz - (tmp - fyz      )/mesh_ptr%xgrid
!---------------------------------------------------------------------------------------------------      
     END SELECT
   END SUBROUTINE curlT
 
! LAPLACIAN OF A SCALAR ============================================================================ 
   SUBROUTINE LapS(f,Lapf)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: f
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: Lapf
    DOUBLE PRECISION, DIMENSION(SIZE(f,1),SIZE(f,2),SIZE(f,3)) :: tmp,dum
    integer :: ifunr,jfunr,kfunr,nfunr 

    !$omp target data map(alloc:tmp,dum)
! SYMMETRY SWITCHES?--------------------------------------------------------------------------------    
     SELECT CASE(patch_ptr%coordsys)
     CASE(0) ! Cartesian
      CALL d2x(f,Lapf)
      CALL d2y(f,tmp)
      CALL d2z(f,dum)
      Lapf = Lapf+tmp+dum
     CASE(1) ! Cylindrical
      CALL ddx(f,tmp)

!$omp target teams distribute parallel do collapse(3)
      do kfunr=1,size(dum,3)
        do jfunr=1,size(dum,2)
          do ifunr=1,size(dum,1)
            dum(ifunr,jfunr,kfunr)= mesh_ptr%xgrid(ifunr,jfunr,kfunr)* tmp(ifunr,jfunr,kfunr)
           end do
        end do
      end do
!$omp end target teams distribute parallel do

      CALL ddx(dum,tmp)

!$omp target teams distribute parallel do collapse(3)
      do kfunr=1,size(Lapf,3)
        do jfunr=1,size(Lapf,2)
          do ifunr=1,size(Lapf,1)
            Lapf(ifunr,jfunr,kfunr)= tmp(ifunr,jfunr,kfunr)/ mesh_ptr%xgrid(ifunr,jfunr,kfunr)
           end do
        end do
      end do
!$omp end target teams distribute parallel do

      CALL d2y(f,tmp)

!$omp target teams distribute parallel do collapse(3)
      do kfunr=1,size(Lapf,3)
        do jfunr=1,size(Lapf,2)
          do ifunr=1,size(Lapf,1)
            dum(ifunr,jfunr,kfunr)= tmp(ifunr,jfunr,kfunr)/ mesh_ptr%xgrid(ifunr,jfunr,kfunr)**2
            Lapf(ifunr,jfunr,kfunr)= Lapf(ifunr,jfunr,kfunr)+ dum(ifunr,jfunr,kfunr)
           end do
        end do
      end do
!$omp end target teams distribute parallel do

      CALL d2z(f,tmp)

!$omp target teams distribute parallel do collapse(3)
      do kfunr=1,size(Lapf,3)
        do jfunr=1,size(Lapf,2)
          do ifunr=1,size(Lapf,1)
            Lapf(ifunr,jfunr,kfunr)= Lapf(ifunr,jfunr,kfunr)+ tmp(ifunr,jfunr,kfunr)
           end do
        end do
      end do
!$omp end target teams distribute parallel do

     CASE(2) ! Spherical
      CALL ddx(f,tmp)
      dum = mesh_ptr%xgrid**2*tmp
      CALL ddx(dum,tmp)
      Lapf = tmp/mesh_ptr%xgrid**2
     
      CALL ddy(f,tmp)
      dum = SIN(mesh_ptr%ygrid)*tmp
      CALL ddy(dum,tmp)
      dum = tmp/(mesh_ptr%xgrid**2*SIN(mesh_ptr%ygrid))
      Lapf = Lapf + dum

      CALL d2z(f,tmp)
      dum = tmp/(mesh_ptr%xgrid*SIN(mesh_ptr%ygrid))**2
      Lapf = Lapf+dum
   END SELECT
     !$omp end target data
   END SUBROUTINE LapS
 
! LAPLACIAN OF A VECTOR ============================================================================ 
   SUBROUTINE LapV(fx,fy,fz,Lx,Ly,Lz)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: fx,fy,fz
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: Lx,Ly,Lz
     SELECT CASE(patch_ptr%coordsys)
     CASE(0) ! Cartesian
     CASE(1) ! Cylindrical
     CASE(2) ! Spherical
     CASE(3) ! General
     END SELECT
     Lx = 0.0D0
     Ly = 0.0D0
     Lz = 0.0D0
   END SUBROUTINE LapV

! CROSS PRODUCT OF TWO VECTORS =================================================================================
   SUBROUTINE VcrossV(ux,uy,uz,cx,cy,cz)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: ux,uy,uz    ! first vector
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: cx,cy,cz   ! second vector on input; result on output
    DOUBLE PRECISION, DIMENSION(SIZE(ux,1),SIZE(ux,2),SIZE(ux,3)) :: tx,ty,tz
     !SELECT CASE(patch_ptr%coordsys)
     !CASE(0) ! Cartesian
     !CASE(1) ! Cylindrical
     !CASE(2) ! Spherical
      tx = cx; ty = cy; tz = cz
      cx = uy*tz - uz*ty
      cy = uz*tx - ux*tz
      cz = ux*ty - uy*tx
!     END SELECT
   END SUBROUTINE VcrossV

! CROSS PRODUCT OF A VECTOR AND A TENSOR =================================================================================
   SUBROUTINE VcrossT(ux,uy,uz,fxx,fxy,fxz,fyx,fyy,fyz,fzx,fzy,fzz,ttxx,ttxy,ttxz,ttyx,ttyy,ttyz,ttzx,ttzy,ttzz)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: ux,uy,uz                                        ! input vector
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: fxx,fxy,fxz,fyx,fyy,fyz,fzx,fzy,fzz             ! input tensor
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: ttxx,ttxy,ttxz,ttyx,ttyy,ttyz,ttzx,ttzy,ttzz   ! output tensor
     !SELECT CASE(patch_ptr%coordsys)
     !CASE(0) ! Cartesian
     !CASE(1) ! Cylindrical
     !CASE(2) ! Spherical
      ttxx = uy*fzx-uz*fyx;    ttxy = uy*fzy-uz*fyy;    ttxz = uy*fzz-uz*fyz
      ttyx = uz*fxx-ux*fzx;    ttyy = uz*fxy-ux*fzy;    ttyz = uz*fxz-ux*fzz
      ttzx = ux*fyx-uy*fxx;    ttzy = ux*fyy-uy*fxy;    ttzz = ux*fyz-uy*fxz
     !END SELECT
   END SUBROUTINE VcrossT
   
! RINGING DETECTOR FOR ARTIFICIAL VISCOSITY ETC.====================================================

   SUBROUTINE ringS(f,fbar,L) ! adds cm^L to the dimensions of f
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: f
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: fbar
    INTEGER,                            INTENT(IN) :: L
    DOUBLE PRECISION, DIMENSION(SIZE(f,1),SIZE(f,2),SIZE(f,3)) :: ring1,ring2,ring3
    integer :: ifunr,jfunr,kfunr,nfunr 

    !CALL exosim_annotation_begin("operators.ringS")
    !$omp target data map(alloc:ring1,ring2,ring3)    
     CALL ringx(f,ring1)
     CALL ringy(f,ring2)
     CALL ringz(f,ring3)
     SELECT CASE(L)
     CASE(0)   

!$omp target teams distribute parallel do collapse(3)
        do kfunr=1,size(fbar,3)
          do jfunr=1,size(fbar,2)
            do ifunr=1,size(fbar,1)
              fbar(ifunr,jfunr,kfunr)=MAX( ABS(ring1(ifunr,jfunr,kfunr)),ABS(ring2(ifunr,jfunr,kfunr)),ABS(ring3(ifunr,jfunr,kfunr)) )
            end do
          end do
        end do
!$omp end target teams distribute parallel do

     CASE(1)                        

!$omp target teams distribute parallel do collapse(3)
        do kfunr=1,size(fbar,3)
          do jfunr=1,size(fbar,2)
            do ifunr=1,size(fbar,1)
              fbar(ifunr,jfunr,kfunr)=MAX( ABS(ring1(ifunr,jfunr,kfunr))*mesh_ptr%d1(ifunr,jfunr,kfunr), ABS(ring2(ifunr,jfunr,kfunr))*mesh_ptr%d2(ifunr,jfunr,kfunr), ABS(ring3(ifunr,jfunr,kfunr))*mesh_ptr%d3(ifunr,jfunr,kfunr))
            end do
          end do
        end do
!$omp end target teams distribute parallel do

     CASE(2)                        

!$omp target teams distribute parallel do collapse(3)
        do kfunr=1,size(fbar,3)
          do jfunr=1,size(fbar,2)
            do ifunr=1,size(fbar,1)
              fbar(ifunr,jfunr,kfunr)=MAX( ABS(ring1(ifunr,jfunr,kfunr))*mesh_ptr%d1(ifunr,jfunr,kfunr)**2,ABS(ring2(ifunr,jfunr,kfunr))*mesh_ptr%d2(ifunr,jfunr,kfunr)**2,ABS(ring3(ifunr,jfunr,kfunr))*mesh_ptr%d3(ifunr,jfunr,kfunr)**2 )
            end do
          end do
        end do
!$omp end target teams distribute parallel do

     END SELECT
     !$omp end target data 
     !CALL exosim_annotation_end("operators.ringS")
   END SUBROUTINE ringS
     
   SUBROUTINE ringV(f,g,h,fbar,L) ! adds cm^L to the dimensions of (f,g,h)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: f,g,h
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: fbar
    INTEGER,                            INTENT(IN) :: L
    !DOUBLE PRECISION, DIMENSION(SIZE(f,1),SIZE(f,2),SIZE(f,3)) :: ringV
    DOUBLE PRECISION, DIMENSION(SIZE(f,1),SIZE(f,2),SIZE(f,3)) :: ringxL,ringyL,ringzL
    DOUBLE PRECISION, DIMENSION(SIZE(f,1),SIZE(f,2),SIZE(f,3)) :: ringx1,ringy1,ringz1
    DOUBLE PRECISION, DIMENSION(SIZE(f,1),SIZE(f,2),SIZE(f,3)) :: ringx2,ringy2,ringz2
    DOUBLE PRECISION, DIMENSION(SIZE(f,1),SIZE(f,2),SIZE(f,3)) :: ringx3,ringy3,ringz3
    integer :: ifunr,jfunr,kfunr,nfunr 

    !CALL exosim_annotation_begin("operators.ringV")
    !$omp target data map(alloc:ringx1,ringx2,ringx3) & 
    !$omp             map(alloc:ringy1,ringy2,ringy3) &  
    !$omp             map(alloc:ringz1,ringz2,ringz3) &
    !$omp             map(alloc:ringxL,ringyL,ringzL)

     CALL ringx(f,ringx1,patch_ptr%isymX)
     CALL ringy(f,ringy1)
     CALL ringz(f,ringz1)

     CALL ringx(g,ringx2)
     CALL ringy(g,ringy2,patch_ptr%isymY)
     CALL ringz(g,ringz2)

     CALL ringx(h,ringx3)
     CALL ringy(h,ringy3)
     CALL ringz(h,ringz3,patch_ptr%isymZ)

     SELECT CASE(L)
     CASE(0)

!$omp target teams distribute parallel do collapse(3)
        do kfunr=1,size(fbar,3)
          do jfunr=1,size(fbar,2)
            do ifunr=1,size(fbar,1)
              ringxL(ifunr,jfunr,kfunr)=MAX( ABS(ringx1(ifunr,jfunr,kfunr)),ABS(ringx2(ifunr,jfunr,kfunr)),ABS(ringx3(ifunr,jfunr,kfunr)) )
              ringyL(ifunr,jfunr,kfunr)=MAX( ABS(ringy1(ifunr,jfunr,kfunr)),ABS(ringy2(ifunr,jfunr,kfunr)),ABS(ringy3(ifunr,jfunr,kfunr)) )
              ringzL(ifunr,jfunr,kfunr)=MAX( ABS(ringz1(ifunr,jfunr,kfunr)),ABS(ringz2(ifunr,jfunr,kfunr)),ABS(ringz3(ifunr,jfunr,kfunr)) )
              fbar(ifunr,jfunr,kfunr)=MAX( ringxL(ifunr,jfunr,kfunr), ringyL(ifunr,jfunr,kfunr), ringzL(ifunr,jfunr,kfunr))
            end do
          end do
        end do
!$omp end target teams distribute parallel do

     CASE(1)

!$omp target teams distribute parallel do collapse(3)
        do kfunr=1,size(fbar,3)
          do jfunr=1,size(fbar,2)
            do ifunr=1,size(fbar,1)
              ringxL(ifunr,jfunr,kfunr)=MAX( ABS(ringx1(ifunr,jfunr,kfunr)),ABS(ringx2(ifunr,jfunr,kfunr)),ABS(ringx3(ifunr,jfunr,kfunr)) )*mesh_ptr%d1(ifunr,jfunr,kfunr)
               ringyL(ifunr,jfunr,kfunr)=MAX( ABS(ringy1(ifunr,jfunr,kfunr)),ABS(ringy2(ifunr,jfunr,kfunr)),ABS(ringy3(ifunr,jfunr,kfunr)) )*mesh_ptr%d2(ifunr,jfunr,kfunr)
               ringzL(ifunr,jfunr,kfunr)=MAX( ABS(ringz1(ifunr,jfunr,kfunr)),ABS(ringz2(ifunr,jfunr,kfunr)),ABS(ringz3(ifunr,jfunr,kfunr)) )*mesh_ptr%d3(ifunr,jfunr,kfunr)
               fbar(ifunr,jfunr,kfunr)=MAX( ringxL(ifunr,jfunr,kfunr), ringyL(ifunr,jfunr,kfunr), ringzL(ifunr,jfunr,kfunr))
            end do
          end do
        end do
!$omp end target teams distribute parallel do

     CASE(2)

!$omp target teams distribute parallel do collapse(3)
        do kfunr=1,size(fbar,3)
          do jfunr=1,size(fbar,2)
            do ifunr=1,size(fbar,1)
              ringxL(ifunr,jfunr,kfunr)=MAX( ABS(ringx1(ifunr,jfunr,kfunr)),ABS(ringx2(ifunr,jfunr,kfunr)),ABS(ringx3(ifunr,jfunr,kfunr)) )*mesh_ptr%d1(ifunr,jfunr,kfunr)**2
              ringyL(ifunr,jfunr,kfunr)=MAX( ABS(ringy1(ifunr,jfunr,kfunr)),ABS(ringy2(ifunr,jfunr,kfunr)),ABS(ringy3(ifunr,jfunr,kfunr)) )*mesh_ptr%d2(ifunr,jfunr,kfunr)**2
              ringzL(ifunr,jfunr,kfunr)=MAX( ABS(ringz1(ifunr,jfunr,kfunr)),ABS(ringz2(ifunr,jfunr,kfunr)),ABS(ringz3(ifunr,jfunr,kfunr)) )*mesh_ptr%d3(ifunr,jfunr,kfunr)**2
              fbar(ifunr,jfunr,kfunr)=MAX( ringxL(ifunr,jfunr,kfunr), ringyL(ifunr,jfunr,kfunr), ringzL(ifunr,jfunr,kfunr))
            end do
          end do
        end do
!$omp end target teams distribute parallel do

     END SELECT
     !$omp end target data
     !CALL exosim_annotation_end("operators.ringV")
   END SUBROUTINE ringV

   SUBROUTINE ringx(f,fbar,bc)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: f
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: fbar
    INTEGER, INTENT(IN), OPTIONAL :: bc
    INTEGER :: n
    integer :: ifunr,jfunr,kfunr,nfunr 

    !CALL exosim_annotation_begin("operators.ringx")
     IF (SIZE(f,1) == 1) THEN

!$omp target teams distribute parallel do collapse(3)
        do kfunr=1,size(fbar,3)
          do jfunr=1,size(fbar,2)
            do ifunr=1,size(fbar,1)
              fbar(ifunr,jfunr,kfunr)= 0.0D0
            end do
          end do
        end do
!$omp end target teams distribute parallel do

     ELSE
        CALL d8x(f,fbar,bc)
     ENDIF
     !CALL exosim_annotation_end("operators.ringx")
   END SUBROUTINE ringx

   SUBROUTINE ringy(f,fbar,bc)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: f
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: fbar
    INTEGER, INTENT(IN), OPTIONAL :: bc
    INTEGER :: n
    integer :: ifunr,jfunr,kfunr,nfunr 

    !CALL exosim_annotation_begin("operators.ringy")
     IF (SIZE(f,2) == 1) THEN

!$omp target teams distribute parallel do collapse(3)
        do kfunr=1,size(fbar,3)
          do jfunr=1,size(fbar,2)
            do ifunr=1,size(fbar,1)
              fbar(ifunr,jfunr,kfunr)= 0.0D0
            end do
          end do
        end do
!$omp end target teams distribute parallel do

     ELSE
       CALL d8y(f,fbar,bc)
     ENDIF
     !CALL exosim_annotation_end("operators.ringy")
   END SUBROUTINE ringy

   SUBROUTINE ringz(f,fbar,bc)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: f
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: fbar
    INTEGER, INTENT(IN), OPTIONAL :: bc
    INTEGER :: n
    integer :: ifunr,jfunr,kfunr,nfunr 

    !CALL exosim_annotation_begin("operators.ringz")
     IF (SIZE(f,3) == 1) THEN

!$omp target teams distribute parallel do collapse(3)
        do kfunr=1,size(fbar,3)
          do jfunr=1,size(fbar,2)
            do ifunr=1,size(fbar,1)
              fbar(ifunr,jfunr,kfunr)= 0.0D0
            end do
          end do
        end do
!$omp end target teams distribute parallel do

     ELSE
       CALL d8z(f,fbar,bc)
     ENDIF
     !CALL exosim_annotation_end("operators.ringz")
   END SUBROUTINE ringz

   SUBROUTINE filterGdir(filtype,fun,bar,direction)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: filtype
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: fun
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: bar
    INTEGER, INTENT(IN)  :: direction ! 1-x,2-y,3-z
    DOUBLE PRECISION, DIMENSION(SIZE(fun,1),SIZE(fun,2),SIZE(fun,3)) :: tmp
    REAL(c_double), DIMENSION(:,:,:), POINTER :: CellBar
    INTEGER :: filnum,xasym,yasym,zasym

    ASSOCIATE( Gfilter=>compact_ptr%control%gfspec )

      IF (direction == 1)  THEN
         CALL bppfx(fun,bar,Gfilter,1)
      END IF
      IF (direction == 2)  THEN
         CALL bppfy(fun,bar,Gfilter,1)
      END IF
      IF (direction == 3)  THEN
         CALL bppfz(fun,bar,Gfilter,1)
      END IF
      
    END ASSOCIATE
  END SUBROUTINE filterGdir
   
! CONSERVATIVE FILTER (except for 'smooth' filtype)=================================================
   SUBROUTINE filter(filtype,fun,bar,component,vsym)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: filtype
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: fun
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: bar
    INTEGER(c_int), INTENT(IN), OPTIONAL :: component ! -1=no filter, 0=scalar, 1=x, 2=y, 3=z
    INTEGER(c_int), INTENT(IN), OPTIONAL :: vsym ! symmetry flag for anti-symmetric magnetic field
    LOGICAL, PARAMETER :: sharp = .TRUE.
    DOUBLE PRECISION, DIMENSION(SIZE(fun,1),SIZE(fun,2),SIZE(fun,3)) :: tmp
    REAL(c_double), DIMENSION(:,:,:), POINTER :: CellBar
    INTEGER :: filnum,xasym,yasym,zasym
    integer :: ifunr,jfunr,kfunr,nfunr 


    !$omp target enter data map(alloc:tmp)
    

!$omp target teams distribute parallel do collapse(3)
    do kfunr=1,size(tmp,3)
      do jfunr=1,size(tmp,2)
        do ifunr=1,size(tmp,1)
          tmp(ifunr,jfunr,kfunr)= 0.0D0
        end do
      end do
    end do
!$omp end target teams distribute parallel do


    
     IF (PRESENT(component)) THEN ! it's a vector
       SELECT CASE(component)
       CASE(-1)     ! Do not filter.
         bar = fun
         !$omp target exit data map(delete:tmp)
         RETURN
       CASE(0)      ! symmetric scalar
         xasym =  1
         yasym =  1
         zasym =  1
       CASE(1)
         xasym = -1 ! antisymmetric normal vector component
         yasym =  1
         zasym =  1
       CASE(2)
         xasym =  1
         yasym = -1 ! antisymmetric normal vector component
         zasym =  1
       CASE(3)
         xasym =  1
         yasym =  1
         zasym = -1 ! antisymmetric normal vector component
       CASE(12)
         xasym = -1 ! antisymmetric in x, y directions
         yasym = -1 ! antisymmetric in x, y directions
         zasym =  1
       CASE(13)
         xasym = -1 ! antisymmetric in x, z directions
         yasym =  1
         zasym = -1 ! antisymmetric in x, z directions
       CASE(23)
         xasym =  1
         yasym = -1 ! antisymmetric in y, z directions
         zasym = -1 ! antisymmetric in y, z directions
       END SELECT     
       IF (PRESENT(vsym)) THEN
         xasym = vsym*xasym
         yasym = vsym*yasym
         zasym = vsym*zasym
       END IF	  
     ELSE ! it's a scalar
       xasym =  1
       yasym =  1
       zasym =  1 
     END IF
 
     SELECT CASE(TRIM(filtype))
     CASE('smooth','Smooth','SMOOTH')
       CALL bppfx(fun,bar,compact_ptr%control%gfspec,1)
       CALL bppfy(bar,tmp,compact_ptr%control%gfspec,1)
       CALL bppfz(tmp,bar,compact_ptr%control%gfspec,1)
       !$omp target exit data map(delete:tmp)
       RETURN
     !CASE('fourier','Fourier','FOURIER')
     !  CALL sfilterx(fun,bar,1,sharp)
     !  CALL sfiltery(bar,tmp,1,sharp)
     !  CALL sfilterz(tmp,bar,1,sharp)
     !  !$omp target exit data map(delete:tmp)
     !  RETURN
     CASE('spectral','Spectral','SPECTRAL','shrpspct','eightord')
       filnum = compact_ptr%control%sfspec
       CellBar => mesh_ptr%CellVolS
     CASE('gaussian','Gaussian','GAUSSIAN')
       filnum = compact_ptr%control%gfspec
       CellBar => mesh_ptr%CellVolG
     CASE('com6ten')
       filnum = 4
       CellBar => mesh_ptr%CellVolS
     END SELECT

     SELECT CASE(patch_ptr%coordsys)
     CASE(0) ! Cartesian
       CALL bppfx(fun,bar,filnum,xasym)
       CALL bppfy(bar,tmp,filnum,yasym)
       CALL bppfz(tmp,bar,filnum,zasym)
     CASE(1) ! Cylindrical
       xasym = -xasym ! mesh_ptr%CellVol is an odd function in radial direction.

!$omp target teams distribute parallel do collapse(3)
       do kfunr=1,size(fun,3)
         do jfunr=1,size(fun,2)
           do ifunr=1,size(fun,1)
             tmp(ifunr,jfunr,kfunr)= fun(ifunr,jfunr,kfunr)* mesh_ptr%CellVol(ifunr,jfunr,kfunr)
            end do
         end do
       end do
!$omp end target teams distribute parallel do

       CALL bppfx(tmp,bar,filnum,xasym)
       CALL bppfy(bar,tmp,filnum,yasym)
       CALL bppfz(tmp,bar,filnum,zasym)

!$omp target teams distribute parallel do collapse(3)
       do kfunr=1,size(bar,3)
         do jfunr=1,size(bar,2)
           do ifunr=1,size(bar,1)
             bar(ifunr,jfunr,kfunr)= bar(ifunr,jfunr,kfunr)/ mesh_ptr%CellVol(ifunr,jfunr,kfunr)
            end do
         end do
       end do
!$omp end target teams distribute parallel do

     CASE(2,3) ! Spherical, Curvilinear
       tmp = fun*mesh_ptr%CellVol
       CALL bppfx(tmp,bar,filnum,xasym)
       CALL bppfy(bar,tmp,filnum,yasym)
       CALL bppfz(tmp,bar,filnum,zasym)
       bar = bar/CellBar
     END SELECT

     !$omp target exit data map(delete:tmp)
   END SUBROUTINE filter

  SUBROUTINE filtRands(Nspan,Ni,No,Nbuff, &
        bmnI,bmnO,rands,vfilt, &
        ny, nz, ay, az, iy1, iz1,y_r )        
     !USE inputs, ONLY: ny,nz
     !USE globals, ONLY: iy,iz,ay,az
     !USE nozfull_data, ONLY: y_r
     IMPLICIT NONE

     INTEGER, INTENT(IN) :: Nspan,Ni,No,Nbuff
     INTEGER, INTENT(IN) :: ny,nz,iy1,iz1,ay,az
     DOUBLE PRECISION, INTENT(IN) :: y_r(ay)
     DOUBLE PRECISION, INTENT(IN) ::  rands(ny+Nbuff,nz+Nbuff)
     DOUBLE PRECISION, INTENT(IN) ::  bmnI(2*Ni+1,2*Nspan+1), bmnO(2*No+1,2*Nspan+1)
     DOUBLE PRECISION, DIMENSION(ay,az), INTENT(OUT) :: vfilt(ay,az)
     DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: filt

     INTEGER :: N1,N2
     INTEGER :: j,k,m,n,mm,nn
     INTEGER :: mF,mG,nF,nG
     
     N2 = Nspan
     DO j=1,ay

        ! Inner Layer
        IF( y_r(j) < 1.0D0 ) THEN
           N1 = Ni
           IF(.NOT. ALLOCATED(filt)) ALLOCATE(filt(SIZE(bmnI,1),SIZE(bmnI,2) )  )
           filt = bmnI
           ! Outer Layer
        ELSE
           N1 = No
           IF(.NOT. ALLOCATED(filt)) ALLOCATE(filt(SIZE(bmnO,1),SIZE(bmnO,2) ) )
           filt = bmnO
        END IF
        
        DO k=1,az
           vfilt(j,k) = 0.0D0
           
           DO m=-N1,N1,1
              mG = (iy1 ) + j + (m + N1)
              mF = m + N1 + 1
              
              DO n=-N2,N2,1
                 nG = (iz1 ) + k + (n + N2)
                 nF = n + N2 + 1
                 !vfilt(j,k) = vfilt(j,k) + filt(mF,nF)*rands(mG,nG)
                 vfilt(j,k) = vfilt(j,k) + filt(mF,nF)*rands( 1+MOD(mG-1,ny) , 1+MOD(nG-1,nz) )
                 
              END DO
           END DO
        END DO
        
        DEALLOCATE(filt)  ! Every new y, compute a new filter
     END DO
     
   END SUBROUTINE filtRands

   SUBROUTINE get_rands_normal(rands,ny,nz,Nbuff,time_seed)
     IMPLICIT NONE
     DOUBLE PRECISION, DIMENSION(4,ny+Nbuff,nz+Nbuff),INTENT(OUT) :: rands
     INTEGER, INTENT(IN) :: ny,nz,Nbuff
     INTEGER, INTENT(IN) :: time_seed 
     integer, allocatable :: seed(:)
     DOUBLE PRECISION, DIMENSION(4,ny+Nbuff,nz+Nbuff) :: rtmp
     integer :: n
     DOUBLE PRECISION :: two=2.0D0,pi=3.14159265359

     call random_seed(size = n)
     allocate(seed(n))
     seed(:) = time_seed
     call random_seed(put=seed)


     ! Get the RNs... same for all procs     
     CALL random_number(rands)

     ! Make them have normal distribution (Box-Mueller theorem)
     rtmp(1:2,:,:) = sqrt( -two*LOG(rands((/1,3/),:,:))) * cos(two*pi*rands((/2,4/),:,:))
     rtmp(3:4,:,:) = sqrt( -two*LOG(rands((/1,3/),:,:))) * sin(two*pi*rands((/2,4/),:,:))
     rands = rtmp
     
   END SUBROUTINE get_rands_normal

   
!===================================================================================================
 END MODULE LES_operators
!=================================================================================================== 


