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
 MODULE LES_operators
!===================================================================================================
  USE iso_c_binding
  USE LES_objects, ONLY : patch_ptr,mesh_ptr,compact_ptr
!  USE LES_matrices, ONLY : ddx,ddy,ddz,d2x,d2y,d2z,d8x,d8y,d8z,bppfx,bppfy,bppfz,Sfilter,Gfilter,Tfilter
  USE LES_compact_operators, ONLY : ddx=>d1x,ddy=>d1y,ddz=>d1z,d2x,d2y,d2z,d8x,d8y,d8z, &
    bppfx=>filterx,bppfy=>filtery,bppfz=>filterz,isrx,isry,isrz,islx,isly,islz
!  USE LES_FFTs, ONLY : sfilterx,sfiltery,sfilterz
 
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
  
!===================================================================================================
  CONTAINS
!===================================================================================================

! DIVERGENCE OF A VECTOR ===========================================================================
   SUBROUTINE divV(fx,fy,fz,df)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: fx,fy,fz
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: df
    DOUBLE PRECISION, DIMENSION(SIZE(df,1),SIZE(df,2),SIZE(df,3)) :: fA,fB,fC,tmp
     SELECT CASE(patch_ptr%coordsys)
    CASE(0) ! Cartesian
      CALL ddx(fx,fA,patch_ptr%isymX)
      CALL ddy(fy,fB,patch_ptr%isymY)
      CALL ddz(fz,fC,patch_ptr%isymZ)
      df = fA + fB + fC
     CASE(1)
      tmp = mesh_ptr%xgrid*fx             ! r*v_r
      CALL ddx(tmp,fA,patch_ptr%isymX**2)
      CALL ddy(fy,fB, patch_ptr%isymY)
      CALL ddz(fz,fC, patch_ptr%isymZ)
      df = (fA+fB)/mesh_ptr%xgrid + fC
     CASE(2)
     ! *** does sin(mesh_ptr%ygrid) ever generate patch_ptr%isymY factor? ***
      tmp = mesh_ptr%xgrid**2*fx          ! r**2*v_r
      df = fy*SIN(mesh_ptr%ygrid)          ! v_theta*SIN(theta)
      CALL ddx(tmp,fA,patch_ptr%isymX**3)
      CALL ddy(df,fB, patch_ptr%isymY)      ! Assumes theta has some symmetry here... move outside derivative
      CALL ddz(fz,fC, patch_ptr%isymZ)
      df = fA/mesh_ptr%xgrid**2 + (fB + fC)/(mesh_ptr%xgrid*SIN(mesh_ptr%ygrid))
     CASE(3) ! General curvilinear
      fA = (fx*mesh_ptr%dAdx + fy*mesh_ptr%dAdy + fz*mesh_ptr%dAdz)*mesh_ptr%detxyz
      fB = (fx*mesh_ptr%dBdx + fy*mesh_ptr%dBdy + fz*mesh_ptr%dBdz)*mesh_ptr%detxyz
      fC = (fx*mesh_ptr%dCdx + fy*mesh_ptr%dCdy + fz*mesh_ptr%dCdz)*mesh_ptr%detxyz
      CALL ddx(fA,df)
      CALL ddy(fB,tmp)
      df = df+tmp
      CALL ddz(fC,tmp)
      df = (df+tmp)/mesh_ptr%detxyz
     END SELECT
   END SUBROUTINE divV
 
! DIVERGENCE OF A TENSOR ===========================================================================
   SUBROUTINE divT(fxx,fxy,fxz,fyx,fyy,fyz,fzx,fzy,fzz,dfx,dfy,dfz)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: fxx,fxy,fxz,fyx,fyy,fyz,fzx,fzy,fzz
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: dfx,dfy,dfz
    DOUBLE PRECISION, DIMENSION(SIZE(fxx,1),SIZE(fxx,2),SIZE(fxx,3)) :: fA,fB,fC,tmp
     SELECT CASE(patch_ptr%coordsys)
     CASE(0) ! Cartesian
      CALL ddx(fxx,fA,patch_ptr%isymX**2)
      CALL ddy(fyx,fB,patch_ptr%isymY)
      CALL ddz(fzx,fC,patch_ptr%isymZ)
      dfx = fA + fB + fC   
      CALL ddx(fxy,fA,patch_ptr%isymX)
      CALL ddy(fyy,fB,patch_ptr%isymY**2)
      CALL ddz(fzy,fC,patch_ptr%isymZ)
      dfy = fA + fB + fC     
      CALL ddx(fxz,fA,patch_ptr%isymX)
      CALL ddy(fyz,fB,patch_ptr%isymY)
      CALL ddz(fzz,fC,patch_ptr%isymZ**2)
      dfz = fA + fB + fC   
     CASE(1)
      tmp = mesh_ptr%xgrid*fxx
      CALL ddx(tmp,fA,patch_ptr%isymX**3)
      CALL ddy(fyx,fB,patch_ptr%isymY)
      CALL ddz(fzx,fC,patch_ptr%isymZ)
      dfx = (fA+fB-fyy)/mesh_ptr%xgrid+fC
      tmp = mesh_ptr%xgrid**2*fxy
      CALL ddx(tmp,fA,patch_ptr%isymX**3)
      CALL ddy(fyy,fB,patch_ptr%isymY**2)
      CALL ddz(fzy,fC,patch_ptr%isymZ)
      dfy = fA/mesh_ptr%xgrid**2 + (fB+fyx-fxy)/mesh_ptr%xgrid + fC
      tmp = mesh_ptr%xgrid*fxz
      CALL ddx(tmp,fA,patch_ptr%isymX**2)
      CALL ddy(fyz,fB,patch_ptr%isymY)
      CALL ddz(fzz,fC,patch_ptr%isymZ**2)
      dfz = (fA+fB)/mesh_ptr%xgrid + fC
    CASE(2)
    ! *** does sin(mesh_ptr%ygrid) ever generate patch_ptr%isymY factor? ***
      tmp = mesh_ptr%xgrid**2*fxx
      CALL ddx(tmp,fA,patch_ptr%isymX**4)
      tmp = fyx*SIN(mesh_ptr%ygrid)
      CALL ddy(tmp,fB,patch_ptr%isymY)
      CALL ddz(fzx,fC,patch_ptr%isymZ)
      dfx = fA/mesh_ptr%xgrid**2 + (fB+fC)/(mesh_ptr%xgrid*SIN(mesh_ptr%ygrid))-(fyy+fzz)/mesh_ptr%xgrid
      tmp = mesh_ptr%xgrid**3*fxy
      CALL ddx(tmp,fA,patch_ptr%isymX**4)
      tmp = fyy*SIN(mesh_ptr%ygrid)
      CALL ddy(tmp,fB,patch_ptr%isymY**2)
      CALL ddz(fzy,fC,patch_ptr%isymZ)
      dfy = fA/mesh_ptr%xgrid**3 + (fB+fC)/(mesh_ptr%xgrid*SIN(mesh_ptr%ygrid)) + (fyx-fxy-fzz/TAN(mesh_ptr%ygrid))/mesh_ptr%xgrid
      tmp = mesh_ptr%xgrid**3*fxz
      CALL ddx(tmp,fA,patch_ptr%isymX**4)
      tmp = fyz*SIN(mesh_ptr%ygrid)
      CALL ddy(tmp,fB,patch_ptr%isymY)
      CALL ddz(fzz,fC,patch_ptr%isymZ**2)
      dfz = fA/mesh_ptr%xgrid**3 + (fB+fC)/(mesh_ptr%xgrid*SIN(mesh_ptr%ygrid)) + (fzx-fxz+fzy/TAN(mesh_ptr%ygrid))/mesh_ptr%xgrid
     CASE(3)
      CALL divV(fxx,fxy,fxz,dfx)
      CALL divV(fyx,fyy,fyz,dfy)
      CALL divV(fzx,fzy,fzz,dfz)
     END SELECT
   END SUBROUTINE divT
 
! GRADIENT OF A SCALAR =============================================================================
   SUBROUTINE gradS(f,dfdx,dfdy,dfdz)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: f
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: dfdx,dfdy,dfdz
    DOUBLE PRECISION, DIMENSION(SIZE(f,1),SIZE(f,2),SIZE(f,3)) :: dfdA,dfdB,dfdC
     CALL ddx(f,dfdx)
     CALL ddy(f,dfdy)
     CALL ddz(f,dfdz)
     SELECT CASE(patch_ptr%coordsys)
     CASE(1) ! Cylindrical
      dfdy = dfdy/mesh_ptr%xgrid
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
   END SUBROUTINE gradS
 
! GRADIENT OF A VECTOR COMPONENT =============================================================================
   SUBROUTINE gradVc(f,dfdx,dfdy,dfdz,vc)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: f
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: dfdx,dfdy,dfdz
    INTEGER, INTENT(IN) :: vc
    DOUBLE PRECISION, DIMENSION(SIZE(f,1),SIZE(f,2),SIZE(f,3)) :: dfdA,dfdB,dfdC
     SELECT CASE( vc )
     CASE( 1 )
       CALL ddx(f,dfdx,patch_ptr%isymX) ; CALL ddy(f,dfdy) ; CALL ddz(f,dfdz)
     CASE( 2 )
       CALL ddx(f,dfdx) ; CALL ddy(f,dfdy,patch_ptr%isymY) ; CALL ddz(f,dfdz)
     CASE( 3 )
       CALL ddx(f,dfdx) ; CALL ddy(f,dfdy) ; CALL ddz(f,dfdz,patch_ptr%isymZ)
     CASE DEFAULT
       CALL ddx(f,dfdx) ; CALL ddy(f,dfdy) ; CALL ddz(f,dfdz)
     END SELECT
     SELECT CASE(patch_ptr%coordsys)
     CASE(1) ! Cylindrical
      dfdy = dfdy/mesh_ptr%xgrid
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
   END SUBROUTINE gradVc
 
! GRADIENT OF A TENSOR COMPONENT =============================================================================
   SUBROUTINE gradTc(f,dfdx,dfdy,dfdz,tc)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: f
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: dfdx,dfdy,dfdz
    CHARACTER(LEN=2), INTENT(IN) :: tc
    DOUBLE PRECISION, DIMENSION(SIZE(f,1),SIZE(f,2),SIZE(f,3)) :: dfdA,dfdB,dfdC
     SELECT CASE( tc )
     CASE( 'xx', 'yy', 'zz' )
       CALL ddx(f,dfdx) ; CALL ddy(f,dfdy) ; CALL ddz(f,dfdz)
     CASE( 'xy', 'yx' )
       CALL ddx(f,dfdx,patch_ptr%isymX) ; CALL ddy(f,dfdy,patch_ptr%isymY) ; CALL ddz(f,dfdz)
     CASE( 'xz', 'zx' )
       CALL ddx(f,dfdx,patch_ptr%isymX) ; CALL ddy(f,dfdy) ; CALL ddz(f,dfdz,patch_ptr%isymZ)
     CASE( 'yz', 'zy' )
       CALL ddx(f,dfdx) ; CALL ddy(f,dfdy,patch_ptr%isymY) ; CALL ddz(f,dfdz,patch_ptr%isymZ)
     CASE DEFAULT
       CALL ddx(f,dfdx) ; CALL ddy(f,dfdy) ; CALL ddz(f,dfdz)
     END SELECT
     SELECT CASE(patch_ptr%coordsys) !-----switches for non-Cartesian to be implemented---
     CASE(1) ! Cylindrical 
      dfdy = dfdy/mesh_ptr%xgrid
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
      dfyx = (dfyx-fy)/mesh_ptr%xgrid
      dfyy = (dfyy+fx)/mesh_ptr%xgrid
      dfyz = dfyz/mesh_ptr%xgrid
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
   SUBROUTINE curlV(fx,fy,fz,cx,cy,cz)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: fx,fy,fz
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: cx,cy,cz
    DOUBLE PRECISION, DIMENSION(SIZE(fx,1),SIZE(fx,2),SIZE(fx,3)) :: tx,ty,tz,tmp
     SELECT CASE(patch_ptr%coordsys)
     CASE(0) ! Cartesian
      CALL ddy(fz,cx)
      CALL ddz(fy,tx)
      cx = cx - tx
      CALL ddz(fx,cy)
      CALL ddx(fz,ty)
      cy = cy - ty
      CALL ddx(fy,cz)
      CALL ddy(fx,tz)
      cz = cz - tz
     CASE(1) ! Cylindrical
      CALL ddy(fz,tmp)
      cx = tmp/mesh_ptr%xgrid
      CALL ddz(fy,tx)
      cx = cx - tx
      CALL ddz(fx,cy)
      CALL ddx(fz,ty)
      cy = cy - ty
      cz = mesh_ptr%xgrid*fy
      CALL ddx(cz,tmp)
      cz = tmp/mesh_ptr%xgrid
      CALL ddy(fx,tmp)
      tz = tmp/mesh_ptr%xgrid
      cz = cz - tz
     CASE(2) ! Spherical
      cx = fz*SIN(mesh_ptr%ygrid)
      CALL ddy(cx,tmp)
      cx = tmp/(mesh_ptr%xgrid*SIN(mesh_ptr%ygrid))
      CALL ddz(fy,tmp)
      tx = tmp/(mesh_ptr%xgrid*SIN(mesh_ptr%ygrid))
      cx = cx - tx
      CALL ddz(fx,tmp)
      cy = tmp/(mesh_ptr%xgrid*SIN(mesh_ptr%ygrid))
      ty = mesh_ptr%xgrid*fz
      CALL ddx(ty,tmp)
      ty = tmp/mesh_ptr%xgrid
      cy = cy - ty
      cz = mesh_ptr%xgrid*fy
      CALL ddx(cz,tmp)
      cz = tmp/mesh_ptr%xgrid
      CALL ddy(fx,tmp)
      tz = tmp/mesh_ptr%xgrid
      cz = cz - tz
     CASE(3) ! Hack to avoid NANs in cylinder_curvy.csv file
      cx = 0.0D0
      cy = 0.0D0
      cz = 0.0D0      
     END SELECT
   END SUBROUTINE curlV
 
! CURL OF A TENSOR =================================================================================
   SUBROUTINE curlT(fxx,fxy,fxz,fyx,fyy,fyz,fzx,fzy,fzz,ttxx,ttxy,ttxz,ttyx,ttyy,ttyz,ttzx,ttzy,ttzz)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: fxx,fxy,fxz,fyx,fyy,fyz,fzx,fzy,fzz
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: ttxx,ttxy,ttxz,ttyx,ttyy,ttyz,ttzx,ttzy,ttzz
    DOUBLE PRECISION, DIMENSION(SIZE(fxx,1),SIZE(fxx,2),SIZE(fxx,3)) :: tmp!tx,ty,tz
     SELECT CASE(patch_ptr%coordsys)
     CASE(0) ! Cartesian
      CALL ddy(fzx,ttxx);                  CALL ddz(fyx,tmp);                 ttxx = ttxx - tmp
      CALL ddy(fzy,ttxy,patch_ptr%isymY);   CALL ddz(fyy,tmp);                ttxy = ttxy - tmp
      CALL ddy(fzz,ttxz);                  CALL ddz(fyz,tmp,patch_ptr%isymZ); ttxz = ttxz - tmp

      CALL ddz(fxx,ttyx);                  CALL ddx(fzx,tmp,patch_ptr%isymX); ttyx = ttyx - tmp
      CALL ddz(fxy,ttyy);                  CALL ddx(fzy,tmp);                 ttyy = ttyy - tmp
      CALL ddz(fxz,ttyz,patch_ptr%isymZ);   CALL ddx(fzz,tmp);                ttyz = ttyz - tmp

      CALL ddx(fyx,ttzx,patch_ptr%isymX);   CALL ddy(fxx,tmp);                ttzx = ttzx - tmp
      CALL ddx(fyy,ttzy);                  CALL ddy(fxy,tmp,patch_ptr%isymY); ttzy = ttzy - tmp
      CALL ddx(fyz,ttzz);                  CALL ddy(fxz,tmp);                 ttzz = ttzz - tmp
     CASE(1) ! Cylindrical ----symmetry switches for non-Cartesian to be fixed-----
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
     END SELECT
   END SUBROUTINE curlT
 
! LAPLACIAN OF A SCALAR ============================================================================ 
   SUBROUTINE LapS(f,Lapf)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: f
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: Lapf
    DOUBLE PRECISION, DIMENSION(SIZE(f,1),SIZE(f,2),SIZE(f,3)) :: tmp,dum
     SELECT CASE(patch_ptr%coordsys)
     CASE(0) ! Cartesian
      CALL d2x(f,Lapf)
      CALL d2y(f,tmp)
      CALL d2z(f,dum)
      Lapf = Lapf+tmp+dum
     CASE(1) ! Cylindrical
      CALL ddx(f,tmp)
      dum = mesh_ptr%xgrid*tmp
      CALL ddx(dum,tmp)
      Lapf = tmp/mesh_ptr%xgrid
     
      CALL d2y(f,tmp)
      dum = tmp/mesh_ptr%xgrid**2
      Lapf = Lapf+dum
     
      CALL d2z(f,tmp)
      Lapf = Lapf+tmp
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
   FUNCTION ring(f)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: f
    DOUBLE PRECISION, DIMENSION(SIZE(f,1),SIZE(f,2),SIZE(f,3)) :: ring
     ring = MAX(ABS(ringx(f))*mesh_ptr%d1**2,ABS(ringy(f))*mesh_ptr%d2**2,ABS(ringz(f))*mesh_ptr%d3**2)
   END FUNCTION ring
     
   FUNCTION ringV(f,g,h)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: f,g,h
    DOUBLE PRECISION, DIMENSION(SIZE(f,1),SIZE(f,2),SIZE(f,3)) :: ringV
    DOUBLE PRECISION, DIMENSION(SIZE(f,1),SIZE(f,2),SIZE(f,3)) :: ringxL,ringyL,ringzL
     ringxL = MAX(ABS(ringx(f,patch_ptr%isymX)),ABS(ringx(g))                ,ABS(ringx(h))                )*mesh_ptr%d1
     ringyL = MAX(ABS(ringy(f))                ,ABS(ringy(g,patch_ptr%isymY)),ABS(ringy(h))                )*mesh_ptr%d2
     ringzL = MAX(ABS(ringz(f))                ,ABS(ringz(g))                ,ABS(ringz(h,patch_ptr%isymZ)))*mesh_ptr%d3
     ringV = MAX(ringxL,ringyL,ringzL)
   END FUNCTION ringV

   FUNCTION ringx(f,bc)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: f
    INTEGER, INTENT(IN), OPTIONAL :: bc
    DOUBLE PRECISION, DIMENSION(SIZE(f,1),SIZE(f,2),SIZE(f,3)) :: ringx
    DOUBLE PRECISION, DIMENSION(-1:SIZE(f,1)+2,SIZE(f,2),SIZE(f,3)) :: g
    INTEGER :: n
     IF (SIZE(f,1) == 1) THEN
       ringx = 0.0D0
     ELSE
       CALL d8x(f,ringx,bc)
     ENDIF
   END FUNCTION ringx

   FUNCTION ringy(f,bc)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: f
    INTEGER, INTENT(IN), OPTIONAL :: bc
    DOUBLE PRECISION, DIMENSION(SIZE(f,1),SIZE(f,2),SIZE(f,3)) :: ringy
    DOUBLE PRECISION, DIMENSION(SIZE(f,1),-1:SIZE(f,2)+2,SIZE(f,3)) :: g
    INTEGER :: n
     IF (SIZE(f,2) == 1) THEN
       ringy = 0.0D0
     ELSE
       CALL d8y(f,ringy,bc)
     ENDIF
   END FUNCTION ringy

   FUNCTION ringz(f,bc)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: f
    INTEGER, INTENT(IN), OPTIONAL :: bc
    DOUBLE PRECISION, DIMENSION(SIZE(f,1),SIZE(f,2),SIZE(f,3)) :: ringz
    DOUBLE PRECISION, DIMENSION(SIZE(f,1),SIZE(f,2),-1:SIZE(f,3)+2) :: g
    INTEGER :: n
     IF (SIZE(f,3) == 1) THEN
       ringz = 0.0D0
     ELSE
       CALL d8z(f,ringz,bc)
     ENDIF
   END FUNCTION ringz


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
   SUBROUTINE filter(filtype,fun,bar,component)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: filtype
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: fun
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: bar
    INTEGER, INTENT(IN), OPTIONAL :: component ! -1=no filter, 0=scalar, 1=x, 2=y, 3=z
    LOGICAL, PARAMETER :: sharp = .TRUE.
    DOUBLE PRECISION, DIMENSION(SIZE(fun,1),SIZE(fun,2),SIZE(fun,3)) :: tmp
    REAL(c_double), DIMENSION(:,:,:), POINTER :: CellBar
    INTEGER :: filnum,xasym,yasym,zasym
    
     IF (PRESENT(component)) THEN
       SELECT CASE(component)
       CASE(-1)     ! Do not filter.
         bar = fun
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
     ELSE
       xasym =  1
       yasym =  1
       zasym =  1 
     END IF
 
     ASSOCIATE( Gfilter=>compact_ptr%control%gfspec, Sfilter=>compact_ptr%control%sfspec, &
                Tfilter=>compact_ptr%control%tfspec )
     
     SELECT CASE(TRIM(filtype))
     CASE('smooth','Smooth','SMOOTH')
       CALL bppfx(fun,bar,Gfilter,1)
       CALL bppfy(bar,tmp,Gfilter,1)
       CALL bppfz(tmp,bar,Gfilter,1)
       RETURN
!     CASE('fourier','Fourier','FOURIER')
!       CALL sfilterx(fun,bar,1,sharp)
!       CALL sfiltery(bar,tmp,1,sharp)
!       CALL sfilterz(tmp,bar,1,sharp)
!       RETURN
     CASE('spectral','Spectral','SPECTRAL','shrpspct','eightord')
       filnum = Sfilter
       CellBar => mesh_ptr%CellVolS
     CASE('gaussian','Gaussian','GAUSSIAN')
       filnum = Gfilter
       CellBar => mesh_ptr%CellVolG
     END SELECT
     
     END ASSOCIATE

     SELECT CASE(patch_ptr%coordsys)
     CASE(0) ! Cartesian
       CALL bppfx(fun,bar,filnum,xasym)
       CALL bppfy(bar,tmp,filnum,yasym)
       CALL bppfz(tmp,bar,filnum,zasym)
     CASE(1) ! Cylindrical
       xasym = -xasym ! mesh_ptr%CellVol is an odd function in radial direction.
       tmp = fun*mesh_ptr%CellVol
       CALL bppfx(tmp,bar,filnum,xasym)
       CALL bppfy(bar,tmp,filnum,yasym)
       CALL bppfz(tmp,bar,filnum,zasym)
       bar = bar/mesh_ptr%CellVol
     CASE(2,3) ! Spherical, Curvilinear
       tmp = fun*mesh_ptr%CellVol
       CALL bppfx(tmp,bar,filnum,xasym)
       CALL bppfy(bar,tmp,filnum,yasym)
       CALL bppfz(tmp,bar,filnum,zasym)
       bar = bar/CellBar
     END SELECT

   END SUBROUTINE filter

! refine by a factor of 3 interpolating left shift(1) = -1/3 and right shift(2) = 1/3
   SUBROUTINE interp_Rubik_compact(f,ff,component)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN)  :: f
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: ff
    INTEGER, INTENT(IN), OPTIONAL :: component ! -1=no refinement, 0=scalar, 1=x, 2=y, 3=z
    DOUBLE PRECISION, DIMENSION(:,:,:,:,:,:), allocatable :: f3
    DOUBLE PRECISION, DIMENSION(SIZE(f,1),SIZE(f,2),SIZE(f,3)) :: tmp
    DOUBLE PRECISION, parameter :: half=0.5D0, third=1.0D0/3.0D0
    INTEGER :: asymx,asymy,asymz,ax,ay,az,bx,by,bz,cx,cy,cz,i,j,k,l,m,n,comp
     ax = SIZE(f,1)  ; ay = SIZE(f,2)  ; az = SIZE(f,3)
     bx = SIZE(ff,1) ; by = SIZE(ff,2) ; bz = SIZE(ff,3)
     if( ax /= bx .and. bx /= 3*ax ) return
     if( ay /= by .and. by /= 3*ay ) return
     if( az /= bz .and. bz /= 3*az ) return
     if( ax == bx .and. ay == by .and. az == bz ) then
       ff = f ; return
     endif
     comp  = 0 ; if(present(component)) comp=component
     if( comp < 0 ) then
       ff = f ; return
     endif
     asymx = 1 ; if( comp == 1 ) asymx = -1
     asymy = 1 ; if( comp == 2 ) asymy = -1
     asymz = 1 ; if( comp == 3 ) asymz = -1
     cx = (bx/ax-1)/2 ; cy = (by/ay-1)/2 ; cz = (bz/az-1)/2 ! 0 or 1
     allocate( f3(ax,ay,az,-cx:cx,-cy:cy,-cz:cz) )
     ! lines, 6 ops max
     f3(:,:,:,0,0,0) = f
     if( bx > ax ) then
       call islx(f,f3(:,:,:,-1,0,0),asymx)
       call isrx(f,f3(:,:,:,+1,0,0),asymx)
       if( by==ay .and. bz==az ) then ! 1D
         forall(i=1:ax,j=1:ay,k=1:az,l=-cx:cx) ff(3*i+l-1,j,k) = f3(i,j,k,l,0,0)
         deallocate(f3)
         return
       endif
     endif
     if( by > ay ) then
       call isly(f,f3(:,:,:,0,-1,0),asymy)
       call isry(f,f3(:,:,:,0,+1,0),asymy)
       if( bx==ax .and. bz==az ) then ! 1D
         forall(i=1:ax,j=1:ay,k=1:az,m=-cy:cy) ff(i,3*j+m-1,k) = f3(i,j,k,0,m,0)
         deallocate(f3)
         return
       endif
     endif
     if( bz > az ) then
       call islz(f,f3(:,:,:,0,0,-1),asymz)
       call isrz(f,f3(:,:,:,0,0,+1),asymz)
       if( bx==ax .and. by==ay ) then ! 1D
         forall(i=1:ax,j=1:ay,k=1:az,n=-cz:cz) ff(i,j,3*k+n-1) = f3(i,j,k,0,0,n)
         deallocate(f3)
         return
       endif
     endif
     ! planes, 24 ops max
     if( bx > ax .and. by > ay ) then ! xy + yx
       call islx(f3(:,:,:,0,-1,0),f3(:,:,:,-1,-1,0),asymx)
       call islx(f3(:,:,:,0,+1,0),f3(:,:,:,-1,+1,0),asymx)
       call isrx(f3(:,:,:,0,-1,0),f3(:,:,:,+1,-1,0),asymx)
       call isrx(f3(:,:,:,0,+1,0),f3(:,:,:,+1,+1,0),asymx)
       call isly(f3(:,:,:,-1,0,0),tmp,asymy)
       f3(:,:,:,-1,-1,0) = half*(f3(:,:,:,-1,-1,0)+tmp)
       call isly(f3(:,:,:,+1,0,0),tmp,asymy)
       f3(:,:,:,+1,-1,0) = half*(f3(:,:,:,+1,-1,0)+tmp)
       call isry(f3(:,:,:,-1,0,0),tmp,asymy)
       f3(:,:,:,-1,+1,0) = half*(f3(:,:,:,-1,+1,0)+tmp)
       call isry(f3(:,:,:,+1,0,0),tmp,asymy)
       f3(:,:,:,+1,+1,0) = half*(f3(:,:,:,+1,+1,0)+tmp)
       if( bz == az ) then ! 2D
         forall(i=1:ax,j=1:ay,k=1:az,l=-cx:cx,m=-cy:cy) ff(3*i+l-1,3*j+m-1,k) = f3(i,j,k,l,m,0)
         deallocate(f3)
         return
       endif
     endif
     if( bx > ax .and. bz > az ) then ! xz + zx
       call islx(f3(:,:,:,0,0,-1),f3(:,:,:,-1,0,-1),asymx)
       call islx(f3(:,:,:,0,0,+1),f3(:,:,:,-1,0,+1),asymx)
       call isrx(f3(:,:,:,0,0,-1),f3(:,:,:,+1,0,-1),asymx)
       call isrx(f3(:,:,:,0,0,+1),f3(:,:,:,+1,0,+1),asymx)
       call islz(f3(:,:,:,-1,0,0),tmp,asymz)
       f3(:,:,:,-1,0,-1) = half*(f3(:,:,:,-1,0,-1)+tmp)
       call islz(f3(:,:,:,+1,0,0),tmp,asymz)
       f3(:,:,:,+1,0,-1) = half*(f3(:,:,:,+1,0,-1)+tmp)
       call isrz(f3(:,:,:,-1,0,0),tmp,asymz)
       f3(:,:,:,-1,0,+1) = half*(f3(:,:,:,-1,0,+1)+tmp)
       call isrz(f3(:,:,:,+1,0,0),tmp,asymz)
       f3(:,:,:,+1,0,+1) = half*(f3(:,:,:,+1,0,+1)+tmp)
       if( by == ay ) then ! 2D
         forall(i=1:ax,j=1:ay,k=1:az,l=-cx:cx,n=-cz:cz) ff(3*i+l-1,j,3*k+n-1) = f3(i,j,k,l,0,n)
         deallocate(f3)
         return
       endif
     endif
     if( by > ay .and. bz > az ) then ! yz + zy
       call isly(f3(:,:,:,0,0,-1),f3(:,:,:,0,-1,-1),asymy)
       call isly(f3(:,:,:,0,0,+1),f3(:,:,:,0,-1,+1),asymy)
       call isry(f3(:,:,:,0,0,-1),f3(:,:,:,0,+1,-1),asymy)
       call isry(f3(:,:,:,0,0,+1),f3(:,:,:,0,+1,+1),asymy)
       call islz(f3(:,:,:,0,-1,0),tmp,asymz)
       f3(:,:,:,0,-1,-1) = half*(f3(:,:,:,0,-1,-1)+tmp)
       call islz(f3(:,:,:,0,+1,0),tmp,asymz)
       f3(:,:,:,0,+1,-1) = half*(f3(:,:,:,0,+1,-1)+tmp)
       call isrz(f3(:,:,:,0,-1,0),tmp,asymz)
       f3(:,:,:,0,-1,+1) = half*(f3(:,:,:,0,-1,+1)+tmp)
       call isrz(f3(:,:,:,0,+1,0),tmp,asymz)
       f3(:,:,:,0,+1,+1) = half*(f3(:,:,:,0,+1,+1)+tmp)
       if( bx == ax ) then ! 2D
         forall(i=1:ax,j=1:ay,k=1:az,m=-cy:cy,n=-cz:cz) ff(i,3*j+m-1,3*k+n-1) = f3(i,j,k,0,m,n)
         deallocate(f3)
         return
       endif
     endif
     ! 3D corners, 24 ops max
     call islx(f3(:,:,:,0,-1,-1),f3(:,:,:,-1,-1,-1),asymx)
     call islx(f3(:,:,:,0,+1,-1),f3(:,:,:,-1,+1,-1),asymx)
     call islx(f3(:,:,:,0,-1,+1),f3(:,:,:,-1,-1,+1),asymx)
     call islx(f3(:,:,:,0,+1,+1),f3(:,:,:,-1,+1,+1),asymx)
     call isrx(f3(:,:,:,0,-1,-1),f3(:,:,:,+1,-1,-1),asymx)
     call isrx(f3(:,:,:,0,+1,-1),f3(:,:,:,+1,+1,-1),asymx)
     call isrx(f3(:,:,:,0,-1,+1),f3(:,:,:,+1,-1,+1),asymx)
     call isrx(f3(:,:,:,0,+1,+1),f3(:,:,:,+1,+1,+1),asymx)
     call isly(f3(:,:,:,-1,0,-1),tmp,asymy)
     f3(:,:,:,-1,-1,-1) = f3(:,:,:,-1,-1,-1)+tmp
     call isly(f3(:,:,:,+1,0,-1),tmp,asymy)
     f3(:,:,:,+1,-1,-1) = f3(:,:,:,+1,-1,-1)+tmp
     call isly(f3(:,:,:,-1,0,+1),tmp,asymy)
     f3(:,:,:,-1,-1,+1) = f3(:,:,:,-1,-1,+1)+tmp
     call isly(f3(:,:,:,+1,0,+1),tmp,asymy)
     f3(:,:,:,+1,-1,+1) = f3(:,:,:,+1,-1,+1)+tmp
     call isry(f3(:,:,:,-1,0,-1),tmp,asymy)
     f3(:,:,:,-1,+1,-1) = f3(:,:,:,-1,+1,-1)+tmp
     call isry(f3(:,:,:,+1,0,-1),tmp,asymy)
     f3(:,:,:,+1,+1,-1) = f3(:,:,:,+1,+1,-1)+tmp
     call isry(f3(:,:,:,-1,0,+1),tmp,asymy)
     f3(:,:,:,-1,+1,+1) = f3(:,:,:,-1,+1,+1)+tmp
     call isry(f3(:,:,:,+1,0,+1),tmp,asymy)
     f3(:,:,:,+1,+1,+1) = f3(:,:,:,+1,+1,+1)+tmp
     call islz(f3(:,:,:,-1,-1,0),tmp,asymz)
     f3(:,:,:,-1,-1,-1) = third*(f3(:,:,:,-1,-1,-1)+tmp)
     call islz(f3(:,:,:,+1,-1,0),tmp,asymz)
     f3(:,:,:,+1,-1,-1) = third*(f3(:,:,:,+1,-1,-1)+tmp)
     call islz(f3(:,:,:,-1,+1,0),tmp,asymz)
     f3(:,:,:,-1,+1,-1) = third*(f3(:,:,:,-1,+1,-1)+tmp)
     call islz(f3(:,:,:,+1,+1,0),tmp,asymz)
     f3(:,:,:,+1,+1,-1) = third*(f3(:,:,:,+1,+1,-1)+tmp)
     call isrz(f3(:,:,:,-1,-1,0),tmp,asymz)
     f3(:,:,:,-1,-1,+1) = third*(f3(:,:,:,-1,-1,+1)+tmp)
     call isrz(f3(:,:,:,+1,-1,0),tmp,asymz)
     f3(:,:,:,+1,-1,+1) = third*(f3(:,:,:,+1,-1,+1)+tmp)
     call isrz(f3(:,:,:,-1,+1,0),tmp,asymz)
     f3(:,:,:,-1,+1,+1) = third*(f3(:,:,:,-1,+1,+1)+tmp)
     call isrz(f3(:,:,:,+1,+1,0),tmp,asymz)
     f3(:,:,:,+1,+1,+1) = third*(f3(:,:,:,+1,+1,+1)+tmp)
     forall(i=1:ax,j=1:ay,k=1:az,l=-cx:cx,m=-cy:cy,n=-cz:cz) ff(3*i+l-1,3*j+m-1,3*k+n-1) = f3(i,j,k,l,m,n)
     deallocate(f3)
   END SUBROUTINE interp_Rubik_compact


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


