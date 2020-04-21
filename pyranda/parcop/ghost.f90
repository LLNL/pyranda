! *** Warning: FORALLs in SELECT CASE statements avoided here due to BGP compiler bug! ***
!===================================================================================================
 MODULE LES_ghost
!=================================================================================================== 
  USE MPI
  USE LES_objects, ONLY : mpistatus,mpierr,patch_ptr,comm_ptr,mesh_ptr,patch_data
  !USE mapp_exosim_annotation, ONLY : exosim_annotation_begin,exosim_annotation_end
  IMPLICIT NONE
  INTERFACE ghost
   MODULE PROCEDURE ghostdS, ghostdV, ghostdVc, ghostiS, ghostiV, ghostiVc
  END INTERFACE

!===================================================================================================
  CONTAINS
!===================================================================================================

   SUBROUTINE ghostdS(gstext,a) ! Fill ghost cells with scalar data from nearest neighbors.
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: gstext ! extrapolation for ghost cells past boundaries: 0=copy bndry, 1=linear
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(INOUT) :: a ! (1-np:ax+np,1-np:ay+np,1-np:az+np), np=number of overlapping planes
     CALL ghostx(gstext,a)
     CALL ghosty(gstext,a)
     CALL ghostz(gstext,a)
   END SUBROUTINE ghostdS

   SUBROUTINE ghostdV(gstext,a,b,c) ! Fill ghost cells with vector data from nearest neighbors.
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: gstext ! extrapolation for ghost cells past boundaries: 0=copy bndry, 1=linear
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(INOUT) :: a,b,c ! (1-np:ax+np,1-np:ay+np,1-np:az+np), np=number of overlapping planes
     CALL ghostx(gstext,a,-1); CALL ghosty(gstext,a)   ; CALL ghostz(gstext,a)
     CALL ghostx(gstext,b)   ; CALL ghosty(gstext,b,-1); CALL ghostz(gstext,b)
     CALL ghostx(gstext,c)   ; CALL ghosty(gstext,c)   ; CALL ghostz(gstext,c,-1)
   END SUBROUTINE ghostdV

   SUBROUTINE ghostdVc(gstext,a,vc) ! Fill ghost cells with vector component data from nearest neighbors.
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: gstext ! extrapolation for ghost cells past boundaries: 0=copy bndry, 1=linear
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(INOUT) :: a ! (1-np:ax+np,1-np:ay+np,1-np:az+np), np=number of overlapping planes
    INTEGER, INTENT(IN) :: vc
    select case( vc )
    case( 1 )
     CALL ghostx(gstext,a,-1) ; CALL ghosty(gstext,a) ; CALL ghostz(gstext,a)
    case( 2 )
     CALL ghostx(gstext,a) ; CALL ghosty(gstext,a,-1) ; CALL ghostz(gstext,a)
    case( 3 )
     CALL ghostx(gstext,a) ; CALL ghosty(gstext,a) ; CALL ghostz(gstext,a,-1)
    case default
     CALL ghostx(gstext,a) ; CALL ghosty(gstext,a) ; CALL ghostz(gstext,a)
    end select
   END SUBROUTINE ghostdVc

   SUBROUTINE ghostx(gstext,a,isym) ! Fill ghost cells in x.
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: gstext ! extrapolation for ghost cells past boundaries: 0=copy bndry, 1=linear
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(INOUT)  :: a ! (1:ax+2*np,:,:)
    DOUBLE PRECISION, DIMENSION((SIZE(a,1)-patch_ptr%ax)/2,size(a,2),size(a,3)) :: a1,a2
    INTEGER, INTENT(IN), OPTIONAL :: isym
    INTEGER :: bx,by,bz,np,n,i,gstx,sym,j,k
    integer :: ifunr,jfunr,kfunr,nfunr 

    
     sym = 1
     if( present(isym) ) sym = isym
     bx = SIZE(a,1)
     by = SIZE(a,2)
     bz = SIZE(a,3)
     np = (bx-patch_ptr%ax)/2
     SELECT CASE(np)
     CASE(:0)
       RETURN
     END SELECT
     n = np*by*bz

     !$omp target data map(from:a1)
     !$omp target teams distribute parallel do collapse(3)
     !a1 = a(   np+1:2*np,:,:)
     do k=1,bz
        do j=1,by
           do i=1,np
              a1(i,j,k) = a( i+np ,j,k)
           end do
        end do
     end do
     !$omp end target teams distribute parallel do
     !$omp end target data
     
     CALL MPI_SENDRECV(a1,n,MPI_DOUBLE_PRECISION,comm_ptr%xcom_lo,0, &
                     & a2,n,MPI_DOUBLE_PRECISION,comm_ptr%xcom_hi,0, &
                     & comm_ptr%xcom,mpistatus,mpierr)

     !$omp target data map(to:a2) map(from:a1)    
     !$omp target teams distribute parallel do collapse(3)
     !a(1:np,:,:) = a2
     !a1 = a(bx-2*np+1:bx-np,:,:)
     do k=1,bz
        do j=1,by
           do i=1,np
              a(bx-np+i,j,k) = a2(i,j,k)
              a1(i,j,k) = a(bx-2*np+i,j,k)
           end do
        end do
     end do
     !$omp end target teams distribute parallel do
     !$omp end target data
     
     CALL MPI_SENDRECV(a1,n,MPI_DOUBLE_PRECISION,comm_ptr%xcom_hi,1, &
                     & a2,n,MPI_DOUBLE_PRECISION,comm_ptr%xcom_lo,1, &
                     & comm_ptr%xcom,mpistatus,mpierr)

     !$omp target data map(to:a2)
     !$omp target teams distribute parallel do collapse(3)
     !a(1:np,:,:) = a2
     do k=1,bz
        do j=1,by
           do i=1,np
              a(i,j,k) = a2(i,j,k)
           end do
        end do
     end do
     !$omp end target teams distribute parallel do
     !$omp end target data
     
     IF (patch_ptr%periodicx) RETURN

     gstx = gstext
     if( patch_ptr%nx == 1 ) gstx = 0
   
     IF (mesh_ptr%x1proc .OR. (patch_ptr%nx .EQ. 1)) THEN
       SELECT CASE(patch_ptr%bx1)
       CASE('GSYM')

!$omp target teams distribute parallel do collapse(2)
          do kfunr=1,size(a,3)
            do jfunr=1,size(a,2)
              DO i=1,np
              a(i,jfunr,kfunr)= a(2*( np )+1-i,jfunr,kfunr)
               END DO
            end do
          end do
!$omp end target teams distribute parallel do

       CASE('SYMM')

!$omp target teams distribute parallel do collapse(2)
          do kfunr=1,size(a,3)
            do jfunr=1,size(a,2)
              DO i=1,np
              a(i,jfunr,kfunr)= dble(sym )*a(2*( np )+1-i,jfunr,kfunr)
               END DO
            end do
          end do
!$omp end target teams distribute parallel do

       CASE DEFAULT
         SELECT CASE( gstx )
         CASE( 0 )

!$omp target teams distribute parallel do collapse(2)
            do kfunr=1,size(a,3)
              do jfunr=1,size(a,2)
                DO i=1,np
                a(i,jfunr,kfunr)= a(np+1,jfunr,kfunr)
                 END DO
              end do
            end do
!$omp end target teams distribute parallel do

         CASE DEFAULT

!$omp target teams distribute parallel do collapse(2)
            do kfunr=1,size(a,3)
              do jfunr=1,size(a,2)
                DO i=1,np
                a(i,jfunr,kfunr)= DBLE(np+2-i )*a(np+1,jfunr,kfunr)- DBLE(np+1-i )*a(np+2,jfunr,kfunr)
                 END DO
              end do
            end do
!$omp end target teams distribute parallel do

         END SELECT
       END SELECT
     END IF
     
     IF (mesh_ptr%xnproc .OR. (patch_ptr%nx .EQ. 1)) THEN
       SELECT CASE(patch_ptr%bxn)
       CASE('GSYM')

!$omp target teams distribute parallel do collapse(2)
         do kfunr=1,size(a,3)
           do jfunr=1,size(a,2)
             DO i=1,np
             a(bx+1-i,jfunr,kfunr)= a(bx-2*( np )+i,jfunr,kfunr)
              END DO
           end do
         end do
!$omp end target teams distribute parallel do

       CASE('SYMM')

!$omp target teams distribute parallel do collapse(2)
         do kfunr=1,size(a,3)
           do jfunr=1,size(a,2)
             DO i=1,np
             a(bx+1-i,jfunr,kfunr)= dble(sym )*a(bx-2*( np )+i,jfunr,kfunr)
              END DO
           end do
         end do
!$omp end target teams distribute parallel do

       CASE DEFAULT
         SELECT CASE( gstx )
         CASE( 0 )

!$omp target teams distribute parallel do collapse(2)
            do kfunr=1,size(a,3)
              do jfunr=1,size(a,2)
                DO i=1,np
                a(bx-np+i,jfunr,kfunr)= a(bx-np,jfunr,kfunr)
                 END DO
              end do
            end do
!$omp end target teams distribute parallel do

         CASE DEFAULT

!$omp target teams distribute parallel do collapse(2)
            do kfunr=1,size(a,3)
              do jfunr=1,size(a,2)
                DO i=1,np
                a(bx-np+i,jfunr,kfunr)= DBLE(i+1 )*a(bx-np,jfunr,kfunr)- DBLE(i )*a(bx-np-1,jfunr,kfunr)
                 END DO
              end do
            end do
!$omp end target teams distribute parallel do

         END SELECT
       END SELECT
     END IF

     
   END SUBROUTINE ghostx
 
   SUBROUTINE ghosty(gstext,a,isym) ! Fill ghost cells in y.
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: gstext ! extrapolation for ghost cells past boundaries: 0=copy bndry, 1=linear
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(INOUT)  :: a ! (:,1:ay+2*np,:)
    DOUBLE PRECISION, DIMENSION(SIZE(a,1),(SIZE(a,2)-patch_ptr%ay)/2,size(a,3)) :: a1,a2
    INTEGER, INTENT(IN), OPTIONAL :: isym
    INTEGER :: bx,by,bz,np,n,i,gstx,sym,j,k
    integer :: ifunr,jfunr,kfunr,nfunr 

    
     sym = 1
     if( present(isym) ) sym = isym
     bx = SIZE(a,1)
     by = SIZE(a,2)
     bz = SIZE(a,3)
     np = (by-patch_ptr%ay)/2
     SELECT CASE(np)
     CASE(:0)
       RETURN
     END SELECT
     n = np*bx*bz
     
     !$omp target data map(from:a1)
     !$omp target teams distribute parallel do collapse(3)
     !a1 = a(   np+1:2*np,:,:)
     do k=1,bz
        do j=1,np
           do i=1,bx
              a1(i,j,k) = a(i,j+np,k)
           end do
        end do
     end do
     !$omp end target teams distribute parallel do
     !$omp end target data
     
     CALL MPI_SENDRECV(a1,n,MPI_DOUBLE_PRECISION,comm_ptr%ycom_lo,0, &
                     & a2,n,MPI_DOUBLE_PRECISION,comm_ptr%ycom_hi,0, &
                     & comm_ptr%ycom,mpistatus,mpierr)

     !$omp target data map(to:a2) map(from:a1)     
     !$omp target teams distribute parallel do collapse(3)
     do k=1,bz
        do j=1,np
           do i=1,bx
              a(i,by-np+j,k) = a2(i,j,k)
              a1(i,j,k) = a(i,by-2*np+j,k)
           end do
        end do
     end do
     !$omp end target teams distribute parallel do
     !$omp end target data
     
     CALL MPI_SENDRECV(a1,n,MPI_DOUBLE_PRECISION,comm_ptr%ycom_hi,1, &
                     & a2,n,MPI_DOUBLE_PRECISION,comm_ptr%ycom_lo,1, &
                     & comm_ptr%ycom,mpistatus,mpierr)

     !$omp target data map(to:a2)
     !$omp target teams distribute parallel do collapse(3)
     do k=1,bz
        do j=1,np
           do i=1,bx
              a(i,j,k) = a2(i,j,k)
           end do
        end do
     end do
     !$omp end target teams distribute parallel do
     !$omp end target data
     
     IF (patch_ptr%periodicy) return

     gstx = gstext
     if( patch_ptr%ny == 1 ) gstx = 0
   
     IF (mesh_ptr%y1proc .OR. (patch_ptr%ny .EQ. 1)) THEN
       SELECT CASE(patch_ptr%by1)
       CASE('GSYM')

!$omp target teams distribute parallel do collapse(2)
         do kfunr=1,size(a,3)
           do ifunr=1,size(a,1)
             DO i=1,np
             a(ifunr, i,kfunr)= a(ifunr, 2*( np )+1-i,kfunr)
              END DO
           end do
         end do
!$omp end target teams distribute parallel do

      CASE('SYMM')

!$omp target teams distribute parallel do collapse(2)
         do kfunr=1,size(a,3)
           do ifunr=1,size(a,1)
             DO i=1,np
             a(ifunr, i,kfunr)= dble(sym )*a(ifunr, 2*( np )+1-i,kfunr)
              END DO
           end do
         end do
!$omp end target teams distribute parallel do

       CASE DEFAULT
         SELECT CASE( gstx )
         CASE( 0 )

!$omp target teams distribute parallel do collapse(2)
            do kfunr=1,size(a,3)
              do ifunr=1,size(a,1)
                DO i=1,np
                a(ifunr, i,kfunr)= a(ifunr, np+1,kfunr)
                 END DO
              end do
            end do
!$omp end target teams distribute parallel do

         CASE DEFAULT

!$omp target teams distribute parallel do collapse(2)
            do kfunr=1,size(a,3)
              do ifunr=1,size(a,1)
                DO i=1,np
                a(ifunr, i,kfunr)= DBLE(np+2-i )*a(ifunr, np+1,kfunr)- DBLE(np+1-i )*a(ifunr, np+2,kfunr)
                 END DO
              end do
            end do
!$omp end target teams distribute parallel do

         END SELECT
       END SELECT
     END IF
     
     IF (mesh_ptr%ynproc .OR. (patch_ptr%ny .EQ. 1)) THEN
       SELECT CASE(patch_ptr%byn)
       CASE('GSYM')

!$omp target teams distribute parallel do collapse(2)
          do kfunr=1,size(a,3)
            do ifunr=1,size(a,1)
              DO i=1,np
              a(ifunr, by+1-i,kfunr)= a(ifunr, by-2*( np )+i,kfunr)
               END DO
            end do
          end do
!$omp end target teams distribute parallel do

       CASE('SYMM')

!$omp target teams distribute parallel do collapse(2)
          do kfunr=1,size(a,3)
            do ifunr=1,size(a,1)
              DO i=1,np
              a(ifunr, by+1-i,kfunr)= dble(sym )*a(ifunr, by-2*( np )+i,kfunr)
               END DO
            end do
          end do
!$omp end target teams distribute parallel do

       CASE DEFAULT
         SELECT CASE( gstx )
         CASE( 0 )

!$omp target teams distribute parallel do collapse(2)
            do kfunr=1,size(a,3)
              do ifunr=1,size(a,1)
                DO i=1,np
                a(ifunr, by-np+i,kfunr)= a(ifunr, by-np,kfunr)
                 END DO
              end do
            end do
!$omp end target teams distribute parallel do

         CASE DEFAULT

!$omp target teams distribute parallel do collapse(2)
            do kfunr=1,size(a,3)
              do ifunr=1,size(a,1)
                DO i=1,np
                a(ifunr, by-np+i,kfunr)= DBLE(i+1 )*a(ifunr, by-np,kfunr)- DBLE(i )*a(ifunr, by-np-1,kfunr)
                 END DO
              end do
            end do
!$omp end target teams distribute parallel do

         END SELECT
       END SELECT
     END IF

     
   END SUBROUTINE ghosty
 
   SUBROUTINE ghostz(gstext,a,isym) ! Fill ghost cells in z.
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: gstext ! extrapolation for ghost cells past boundaries: 0=copy bndry, 1=linear
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(INOUT)  :: a ! (:,:,1:az+2*np)
    DOUBLE PRECISION, DIMENSION(SIZE(a,1),SIZE(a,2),(SIZE(a,3)-patch_ptr%az)/2) :: a1,a2
    INTEGER, INTENT(IN), OPTIONAL :: isym
    INTEGER :: bx,by,bz,np,n,i,gstx,sym,j,k
    integer :: ifunr,jfunr,kfunr,nfunr 

    
     sym = 1
     if( present(isym) ) sym = isym
     bx = SIZE(a,1)
     by = SIZE(a,2)
     bz = SIZE(a,3)
     np = (bz-patch_ptr%az)/2
     SELECT CASE(np)
     CASE(:0)
       RETURN
     END SELECT
     n = np*bx*by

     !$omp target data map(from:a1)
     !$omp target teams distribute parallel do collapse(3)
     !a1 = a(   np+1:2*np,:,:)
     do k=1,np
        do j=1,by
           do i=1,bx
              a1(i,j,k) = a(i,j,k+np)
           end do
        end do
     end do
     !$omp end target teams distribute parallel do
     !$omp end target data
     
     CALL MPI_SENDRECV(a1,n,MPI_DOUBLE_PRECISION,comm_ptr%zcom_lo,0, &
                     & a2,n,MPI_DOUBLE_PRECISION,comm_ptr%zcom_hi,0, &
                     & comm_ptr%zcom,mpistatus,mpierr)

     !$omp target data map(to:a2) map(from:a1)
     !$omp target teams distribute parallel do collapse(3)
     do k=1,np
        do j=1,by
           do i=1,bx
              a(i,j,bz-np+k) = a2(i,j,k)
              a1(i,j,k) = a(i,j,bz-2*np+k)
           end do
        end do
     end do
     !$omp end target teams distribute parallel do
     !$omp end target data
     
     CALL MPI_SENDRECV(a1,n,MPI_DOUBLE_PRECISION,comm_ptr%zcom_hi,1, &
                     & a2,n,MPI_DOUBLE_PRECISION,comm_ptr%zcom_lo,1, &
                     & comm_ptr%zcom,mpistatus,mpierr)

     !$omp target data map(to:a2)
     !$omp target teams distribute parallel do collapse(3)
     do k=1,np
        do j=1,by
           do i=1,bx
              a(i,j,k) = a2(i,j,k)
           end do
        end do
     end do
     !$omp end target teams distribute parallel do
     !$omp end target data

     IF (patch_ptr%periodicz) return

     gstx = gstext
     if( patch_ptr%nz == 1 ) gstx = 0
   
     IF (mesh_ptr%z1proc .OR. (patch_ptr%nz .EQ. 1)) THEN
       SELECT CASE(patch_ptr%bz1)
       CASE('GSYM')

!$omp target teams distribute parallel do collapse(2)
          do jfunr=1,size(a,2)
            do ifunr=1,size(a,1)
              DO i=1,np
              a(ifunr,jfunr, i )=a(ifunr,jfunr, 2*( np )+1-i )
              END DO
            end do
          end do
!$omp end target teams distribute parallel do

       CASE('SYMM')

!$omp target teams distribute parallel do collapse(2)
          do jfunr=1,size(a,2)
            do ifunr=1,size(a,1)
              DO i=1,np
              a(ifunr,jfunr, i )=dble(sym )*a(ifunr,jfunr, 2*( np )+1-i )
              END DO
            end do
          end do
!$omp end target teams distribute parallel do

       CASE DEFAULT
         SELECT CASE( gstx )
         CASE( 0 )

!$omp target teams distribute parallel do collapse(2)
            do jfunr=1,size(a,2)
              do ifunr=1,size(a,1)
                DO i=1,np
                a(ifunr,jfunr, i )=a(ifunr,jfunr, np+1 )
                END DO
              end do
            end do
!$omp end target teams distribute parallel do

         CASE DEFAULT

!$omp target teams distribute parallel do collapse(2)
            do jfunr=1,size(a,2)
              do ifunr=1,size(a,1)
                DO i=1,np
                a(ifunr,jfunr, i )=DBLE(np+2-i )*a(ifunr,jfunr, np+1 )-DBLE(np+1-i )*a(ifunr,jfunr, np+2 )
                END DO
              end do
            end do
!$omp end target teams distribute parallel do

         END SELECT
       END SELECT
     END IF
     
     IF (mesh_ptr%znproc .OR. (patch_ptr%nz .EQ. 1)) THEN
       SELECT CASE(patch_ptr%bzn)
       CASE('GSYM')

!$omp target teams distribute parallel do collapse(2)
          do jfunr=1,size(a,2)
            do ifunr=1,size(a,1)
              DO i=1,np
              a(ifunr,jfunr, bz+1-i )=a(ifunr,jfunr, bz-2*( np )+i )
              END DO
            end do
          end do
!$omp end target teams distribute parallel do

       CASE('SYMM')

!$omp target teams distribute parallel do collapse(2)
          do jfunr=1,size(a,2)
            do ifunr=1,size(a,1)
              DO i=1,np
              a(ifunr,jfunr, bz+1-i )=dble(sym )*a(ifunr,jfunr, bz-2*( np )+i )
              END DO
            end do
          end do
!$omp end target teams distribute parallel do

       CASE DEFAULT
         SELECT CASE( gstx )
         CASE( 0 )

!$omp target teams distribute parallel do collapse(2)
            do jfunr=1,size(a,2)
              do ifunr=1,size(a,1)
                DO i=1,np
                a(ifunr,jfunr, bz-np+i )=a(ifunr,jfunr, bz-np )
                END DO
              end do
            end do
!$omp end target teams distribute parallel do

         CASE DEFAULT

!$omp target teams distribute parallel do collapse(2)
            do jfunr=1,size(a,2)
              do ifunr=1,size(a,1)
                DO i=1,np
                a(ifunr,jfunr, bz-np+i )=DBLE(i+1 )*a(ifunr,jfunr, bz-np )-DBLE(i )*a(ifunr,jfunr, bz-np-1 )
                END DO
              end do
            end do
!$omp end target teams distribute parallel do

         END SELECT
       END SELECT
     END IF

     
   END SUBROUTINE ghostz

   SUBROUTINE ghostiS(gstext,a) ! Fill ghost cells with (integer) scalar data from nearest neighbors.
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: gstext ! extrapolation for ghost cells past boundaries: 0=copy bndry, 1=linear
    INTEGER, DIMENSION(:,:,:), INTENT(INOUT)  :: a ! (1-np:ax+np,1-np:ay+np,1-np:az+np), np = number of overlap planes
     CALL ghostix(gstext,a)
     CALL ghostiy(gstext,a)
     CALL ghostiz(gstext,a)
   END SUBROUTINE ghostiS

   SUBROUTINE ghostiV(gstext,a,b,c) ! Fill ghost cells with vector data from nearest neighbors.
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: gstext ! extrapolation for ghost cells past boundaries: 0=copy bndry, 1=linear
    INTEGER, DIMENSION(:,:,:), INTENT(INOUT) :: a,b,c ! (1-np:ax+np,1-np:ay+np,1-np:az+np), np=number of overlapping planes
     CALL ghostix(gstext,a,-1) ; CALL ghostiy(gstext,a) ; CALL ghostiz(gstext,a)
     CALL ghostix(gstext,b) ; CALL ghostiy(gstext,b,-1) ; CALL ghostiz(gstext,b)
     CALL ghostix(gstext,c) ; CALL ghostiy(gstext,c) ; CALL ghostiz(gstext,c,-1)
   END SUBROUTINE ghostiV

   SUBROUTINE ghostiVc(gstext,a,vc) ! Fill ghost cells with vector component data from nearest neighbors.
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: gstext ! extrapolation for ghost cells past boundaries: 0=copy bndry, 1=linear
    INTEGER, DIMENSION(:,:,:), INTENT(INOUT) :: a ! (1-np:ax+np,1-np:ay+np,1-np:az+np), np=number of overlapping planes
    INTEGER, INTENT(IN) :: vc
    select case( vc )
    case( 1 )
     CALL ghostix(gstext,a,-1) ; CALL ghostiy(gstext,a) ; CALL ghostiz(gstext,a)
    case( 2 )
     CALL ghostix(gstext,a) ; CALL ghostiy(gstext,a,-1) ; CALL ghostiz(gstext,a)
    case( 3 )
     CALL ghostix(gstext,a) ; CALL ghostiy(gstext,a) ; CALL ghostiz(gstext,a,-1)
    case default
     CALL ghostix(gstext,a) ; CALL ghostiy(gstext,a) ; CALL ghostiz(gstext,a)
    end select
   END SUBROUTINE ghostiVc

   SUBROUTINE ghostix(gstext,a,isym) ! Fill ghost cells in x.
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: gstext ! extrapolation for ghost cells past boundaries: 0=copy bndry, 1=linear
    INTEGER, DIMENSION(:,:,:), INTENT(INOUT)  :: a ! (1:ax+2*np,:,:)
    INTEGER, INTENT(IN), OPTIONAL :: isym
    INTEGER :: bx,by,bz,np,n,i,gstx,sym
     sym = 1
     if( present(isym) ) sym = isym
     bx = SIZE(a,1)
     by = SIZE(a,2)
     bz = SIZE(a,3)
     np = (bx-patch_ptr%ax)/2
     SELECT CASE(np)
     CASE(:0)
       RETURN
     END SELECT
     n = np*by*bz
     CALL MPI_SENDRECV(a(   np+1:2*np,:,:),n,MPI_INTEGER,comm_ptr%xcom_lo,0, &
                     & a(bx-np+1:bx  ,:,:),n,MPI_INTEGER,comm_ptr%xcom_hi,0, &
                     & comm_ptr%xcom,mpistatus,mpierr)
     CALL MPI_SENDRECV(a(bx-2*np+1:bx-np,:,:),n,MPI_INTEGER,comm_ptr%xcom_hi,1, &
                     & a(1:np,:,:)           ,n,MPI_INTEGER,comm_ptr%xcom_lo,1, &
                     & comm_ptr%xcom,mpistatus,mpierr)
           
     IF (patch_ptr%periodicx) return

     gstx = gstext
     if( patch_ptr%nx == 1 ) gstx = 0
   
     IF (mesh_ptr%x1proc .OR. (patch_ptr%nx .EQ. 1)) THEN
       SELECT CASE(patch_ptr%bx1)
       CASE('GSYM')
         DO i=1,np
           a(i,:,:) = a(2*(np)+1-i,:,:)
         END DO
       CASE('SYMM')
         DO i=1,np
           a(i,:,:) = sym*a(2*(np)+1-i,:,:)
         END DO
       CASE DEFAULT
         SELECT CASE( gstx )
         CASE( 0 )
           DO i=1,np ; a(i,:,:) = a(np+1,:,:) ; END DO
         CASE DEFAULT
           DO i=1,np ; a(i,:,:) = (np+2-i)*a(np+1,:,:) - (np+1-i)*a(np+2,:,:) ; END DO
         END SELECT
       END SELECT
     END IF
     
     IF (mesh_ptr%xnproc .OR. (patch_ptr%nx .EQ. 1)) THEN
       SELECT CASE(patch_ptr%bxn)
       CASE('GSYM')
         DO i=1,np
           a(bx+1-i,:,:) = a(bx-2*(np)+i,:,:)
         END DO
       CASE('SYMM')
         DO i=1,np
           a(bx+1-i,:,:) = sym*a(bx-2*(np)+i,:,:)
         END DO
       CASE DEFAULT
         SELECT CASE( gstx )
         CASE( 0 )
           DO i=1,np ; a(bx-np+i,:,:) = a(bx-np,:,:) ; END DO
         CASE DEFAULT
           DO i=1,np ; a(bx-np+i,:,:) = (i+1)*a(bx-np,:,:) - (i)*a(bx-np-1,:,:) ; END DO
         END SELECT
       END SELECT
     END IF
   
   END SUBROUTINE ghostix
 
   SUBROUTINE ghostiy(gstext,a,isym) ! Fill ghost cells in y.
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: gstext ! extrapolation for ghost cells past boundaries: 0=copy bndry, 1=linear
    INTEGER, DIMENSION(:,:,:), INTENT(INOUT)  :: a ! (:,1:ay+2*np,:)
    INTEGER, INTENT(IN), OPTIONAL :: isym
    INTEGER :: bx,by,bz,np,n,i,gstx,sym
     sym = 1
     if( present(isym) ) sym = isym
     bx = SIZE(a,1)
     by = SIZE(a,2)
     bz = SIZE(a,3)
     np = (by-patch_ptr%ay)/2
     SELECT CASE(np)
     CASE(:0)
       RETURN
     END SELECT
     n = np*bx*bz
     CALL MPI_SENDRECV(a(:,   np+1:2*np,:),n,MPI_INTEGER,comm_ptr%ycom_lo,0, &
                     & a(:,by-np+1:by  ,:),n,MPI_INTEGER,comm_ptr%ycom_hi,0, &
                     & comm_ptr%ycom,mpistatus,mpierr)
     CALL MPI_SENDRECV(a(:,by-2*np+1:by-np,:),n,MPI_INTEGER,comm_ptr%ycom_hi,1, &
                     & a(:,1:np,:)           ,n,MPI_INTEGER,comm_ptr%ycom_lo,1, &
                     & comm_ptr%ycom,mpistatus,mpierr)
           
     IF (patch_ptr%periodicy) return

     gstx = gstext
     if( patch_ptr%ny == 1 ) gstx = 0
   
     IF (mesh_ptr%y1proc .OR. (patch_ptr%ny .EQ. 1)) THEN
       SELECT CASE(patch_ptr%by1)
       CASE('GSYM')
         DO i=1,np
           a(:,i,:) = a(:,2*(np)+1-i,:)
         END DO
       CASE('SYMM')
         DO i=1,np
           a(:,i,:) = sym*a(:,2*(np)+1-i,:)
         END DO
       CASE DEFAULT
         SELECT CASE( gstx )
         CASE( 0 )
           DO i=1,np ; a(:,i,:) = a(:,np+1,:) ; END DO
         CASE DEFAULT
           DO i=1,np ; a(:,i,:) = (np+2-i)*a(:,np+1,:) - (np+1-i)*a(:,np+2,:) ; END DO
         END SELECT
       END SELECT
     END IF
     
     IF (mesh_ptr%ynproc .OR. (patch_ptr%ny .EQ. 1)) THEN
       SELECT CASE(patch_ptr%byn)
       CASE('GSYM')
         DO i=1,np
           a(:,by+1-i,:) = a(:,by-2*(np)+i,:)
         END DO
       CASE('SYMM')
         DO i=1,np
           a(:,by+1-i,:) = sym*a(:,by-2*(np)+i,:)
         END DO
       CASE DEFAULT
         SELECT CASE( gstx )
         CASE( 0 )
           DO i=1,np ; a(:,by-np+i,:) = a(:,by-np,:) ; END DO
         CASE DEFAULT
           DO i=1,np ; a(:,by-np+i,:) = (i+1)*a(:,by-np,:) - (i)*a(:,by-np-1,:) ; END DO
         END SELECT
       END SELECT
     END IF

   END SUBROUTINE ghostiy
 
   SUBROUTINE ghostiz(gstext,a,isym) ! Fill ghost cells in z.
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: gstext ! extrapolation for ghost cells past boundaries: 0=copy bndry, 1=linear
    INTEGER, DIMENSION(:,:,:), INTENT(INOUT)  :: a ! (:,:,1:az+2*np)
    INTEGER, INTENT(IN), OPTIONAL :: isym
    INTEGER :: bx,by,bz,np,n,i,gstx,sym
     sym = 1
     if( present(isym) ) sym = isym
     bx = SIZE(a,1)
     by = SIZE(a,2)
     bz = SIZE(a,3)
     np = (bz-patch_ptr%az)/2
     SELECT CASE(np)
     CASE(:0)
       RETURN
     END SELECT
     n = np*bx*by
     CALL MPI_SENDRECV(a(:,:,   np+1:2*np),n,MPI_INTEGER,comm_ptr%zcom_lo,0, &
                     & a(:,:,bz-np+1:bz  ),n,MPI_INTEGER,comm_ptr%zcom_hi,0, &
                    & comm_ptr%zcom,mpistatus,mpierr)
     CALL MPI_SENDRECV(a(:,:,bz-2*np+1:bz-np),n,MPI_INTEGER,comm_ptr%zcom_hi,1, &
                     & a(:,:,1:np)           ,n,MPI_INTEGER,comm_ptr%zcom_lo,1, &
                     & comm_ptr%zcom,mpistatus,mpierr)
           
     IF (patch_ptr%periodicz) return

     gstx = gstext
     if( patch_ptr%nz == 1 ) gstx = 0
   
     IF (mesh_ptr%z1proc .OR. (patch_ptr%nz .EQ. 1)) THEN
       SELECT CASE(patch_ptr%bz1)
       CASE('GSYM')
         DO i=1,np
           a(:,:,i) = a(:,:,2*(np)+1-i)
         END DO
       CASE('SYMM')
         DO i=1,np
           a(:,:,i) = sym*a(:,:,2*(np)+1-i)
         END DO
       CASE DEFAULT
         SELECT CASE( gstx )
         CASE( 0 )
           DO i=1,np ; a(:,:,i) = a(:,:,np+1) ; END DO
         CASE DEFAULT
           DO i=1,np ; a(:,:,i) = (np+2-i)*a(:,:,np+1) - (np+1-i)*a(:,:,np+2) ; END DO
         END SELECT
       END SELECT
     END IF
     
     IF (mesh_ptr%znproc .OR. (patch_ptr%nz .EQ. 1)) THEN
       SELECT CASE(patch_ptr%bzn)
       CASE('GSYM')
         DO i=1,np
           a(:,:,bz+1-i) = a(:,:,bz-2*(np)+i)
         END DO
       CASE('SYMM')
         DO i=1,np
           a(:,:,bz+1-i) = sym*a(:,:,bz-2*(np)+i)
         END DO
       CASE DEFAULT
         SELECT CASE( gstx )
         CASE( 0 )
           DO i=1,np ; a(:,:,bz-np+i) = a(:,:,bz-np) ; END DO
         CASE DEFAULT
           DO i=1,np ; a(:,:,bz-np+i) = (i+1)*a(:,:,bz-np) - (i)*a(:,:,bz-np-1) ; END DO
         END SELECT
       END SELECT
     END IF

   END SUBROUTINE ghostiz
   
! FOR IMPLICIT CONDUCTION AND RADIATION ROUTINES----------------------------------------------------
   SUBROUTINE faces(D,Dip,Dim,Djp,Djm,Dkp,Dkm)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: D
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: Dip,Dim,Djp,Djm,Dkp,Dkm
    DOUBLE PRECISION, DIMENSION(0:patch_ptr%ax+1,0:patch_ptr%ay+1,0:patch_ptr%az+1) :: Dghost
    INTEGER :: i,j,k
     Dghost(1:patch_ptr%ax,1:patch_ptr%ay,1:patch_ptr%az) = D
     CALL ghost(0,Dghost)
     FORALL(i=1:patch_ptr%ax,j=1:patch_ptr%ay,k=1:patch_ptr%az) Dip(i,j,k) = 0.5D0*(Dghost(i,j,k) + Dghost(i+1,j  ,k  ))
     FORALL(i=1:patch_ptr%ax,j=1:patch_ptr%ay,k=1:patch_ptr%az) Dim(i,j,k) = 0.5D0*(Dghost(i,j,k) + Dghost(i-1,j  ,k  ))
     FORALL(i=1:patch_ptr%ax,j=1:patch_ptr%ay,k=1:patch_ptr%az) Djp(i,j,k) = 0.5D0*(Dghost(i,j,k) + Dghost(i  ,j+1,k  ))
     FORALL(i=1:patch_ptr%ax,j=1:patch_ptr%ay,k=1:patch_ptr%az) Djm(i,j,k) = 0.5D0*(Dghost(i,j,k) + Dghost(i  ,j-1,k  ))
     FORALL(i=1:patch_ptr%ax,j=1:patch_ptr%ay,k=1:patch_ptr%az) Dkp(i,j,k) = 0.5D0*(Dghost(i,j,k) + Dghost(i  ,j  ,k+1))
     FORALL(i=1:patch_ptr%ax,j=1:patch_ptr%ay,k=1:patch_ptr%az) Dkm(i,j,k) = 0.5D0*(Dghost(i,j,k) + Dghost(i  ,j  ,k-1))
   END SUBROUTINE faces

!---------------------------------------------------------------------------
! FOR IMPLICIT RADIATION ROUTINE: IMPLEMENTS AVERAGING AS IN HOWELL AND
! GREENOUGH (2003), FIGURE 4.  CONDUCTION ROUTINE HANDLES COEFFICIENTS
! DIFFERENTLY SO THIS FORMULA WOULD NOT BE APPROPRIATE WITHOUT MODIFICATION.
!---------------------------------------------------------------------------
   SUBROUTINE radfaces(D,Dip,Dim,Djp,Djm,Dkp,Dkm)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: D
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: Dip,Dim,Djp,Djm,Dkp,Dkm
    DOUBLE PRECISION, DIMENSION(0:patch_ptr%ax+1,0:patch_ptr%ay+1,0:patch_ptr%az+1) :: Dghost
    INTEGER :: i,j,k
     Dghost(1:patch_ptr%ax,1:patch_ptr%ay,1:patch_ptr%az) = D
     CALL ghost(0,Dghost)
     FORALL(i=1:patch_ptr%ax,j=1:patch_ptr%ay,k=1:patch_ptr%az)
       Dip(i,j,k) = MIN(0.5D0  *  (Dghost(i,j,k) + Dghost(i+1,j  ,k  )), &
                        MAX(2.D0 * Dghost(i,j,k) * Dghost(i+1,j  ,k  )   &
                                / (Dghost(i,j,k) + Dghost(i+1,j  ,k  )), &
                            4.D0 / (3.D0 * patch_ptr%dx)))
       Dim(i,j,k) = MIN(0.5D0  *  (Dghost(i,j,k) + Dghost(i-1,j  ,k  )), &
                        MAX(2.D0 * Dghost(i,j,k) * Dghost(i-1,j  ,k  )   &
                                / (Dghost(i,j,k) + Dghost(i-1,j  ,k  )), &
                            4.D0 / (3.D0 * patch_ptr%dx)))
       Djp(i,j,k) = MIN(0.5D0  *  (Dghost(i,j,k) + Dghost(i  ,j+1,k  )), &
                        MAX(2.D0 * Dghost(i,j,k) * Dghost(i  ,j+1,k  )   &
                                / (Dghost(i,j,k) + Dghost(i  ,j+1,k  )), &
                            4.D0 / (3.D0 * patch_ptr%dy)))
       Djm(i,j,k) = MIN(0.5D0  *  (Dghost(i,j,k) + Dghost(i  ,j-1,k  )), &
                        MAX(2.D0 * Dghost(i,j,k) * Dghost(i  ,j-1,k  )   &
                                / (Dghost(i,j,k) + Dghost(i  ,j-1,k  )), &
                            4.D0 / (3.D0 * patch_ptr%dy)))
       Dkp(i,j,k) = MIN(0.5D0  *  (Dghost(i,j,k) + Dghost(i  ,j  ,k+1)), &
                        MAX(2.D0 * Dghost(i,j,k) * Dghost(i  ,j  ,k+1)   &
                                / (Dghost(i,j,k) + Dghost(i  ,j  ,k+1)), &
                            4.D0 / (3.D0 * patch_ptr%dz)))
       Dkm(i,j,k) = MIN(0.5D0  *  (Dghost(i,j,k) + Dghost(i  ,j  ,k-1)), &
                        MAX(2.D0 * Dghost(i,j,k) * Dghost(i  ,j  ,k-1)   &
                                / (Dghost(i,j,k) + Dghost(i  ,j  ,k-1)), &
                            4.D0 / (3.D0 * patch_ptr%dz)))
     END FORALL
   END SUBROUTINE radfaces

   SUBROUTINE harmfaces(D,Dip,Dim,Djp,Djm,Dkp,Dkm)
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: D
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: Dip,Dim,Djp,Djm,Dkp,Dkm
    DOUBLE PRECISION, DIMENSION(0:patch_ptr%ax+1,0:patch_ptr%ay+1,0:patch_ptr%az+1) :: Dghost
    INTEGER :: i,j,k
     Dghost(1:patch_ptr%ax,1:patch_ptr%ay,1:patch_ptr%az) = D
     CALL ghost(0,Dghost)
     FORALL(i=1:patch_ptr%ax,j=1:patch_ptr%ay,k=1:patch_ptr%az)
       Dip(i,j,k) = (2.D0 * Dghost(i,j,k) * Dghost(i+1,j  ,k  ))   &
                         / (Dghost(i,j,k) + Dghost(i+1,j  ,k  ))
       Dim(i,j,k) = (2.D0 * Dghost(i,j,k) * Dghost(i-1,j  ,k  ))   &
                         / (Dghost(i,j,k) + Dghost(i-1,j  ,k  ))
       Djp(i,j,k) = (2.D0 * Dghost(i,j,k) * Dghost(i  ,j+1,k  ))   &
                         / (Dghost(i,j,k) + Dghost(i  ,j+1,k  ))
       Djm(i,j,k) = (2.D0 * Dghost(i,j,k) * Dghost(i  ,j-1,k  ))   &
                         / (Dghost(i,j,k) + Dghost(i  ,j-1,k  ))
       Dkp(i,j,k) = (2.D0 * Dghost(i,j,k) * Dghost(i  ,j  ,k+1))   &
                         / (Dghost(i,j,k) + Dghost(i  ,j  ,k+1))
       Dkm(i,j,k) = (2.D0 * Dghost(i,j,k) * Dghost(i  ,j  ,k-1))   &
                         / (Dghost(i,j,k) + Dghost(i  ,j  ,k-1))
     END FORALL
   END SUBROUTINE harmfaces
 
!===================================================================================================
 END MODULE LES_ghost
!===================================================================================================
