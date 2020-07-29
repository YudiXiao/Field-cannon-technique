MODULE oalloc

  IMPLICIT NONE

  INTEGER :: info
  INTEGER,PARAMETER :: MKS = KIND(1.0E0)
  INTEGER,PARAMETER :: MKD = KIND(1.0D0)
  REAL(MKS),DIMENSION(:,:),ALLOCATABLE :: swork
  REAL(MKD),DIMENSION(:,:),ALLOCATABLE :: dwork

  INTERFACE alloc
      MODULE PROCEDURE salloc,dalloc
  END INTERFACE alloc

  CONTAINS
      SUBROUTINE salloc(array,Dimx,Dimy)
      INTEGER :: Dimx,Dimy      
      REAL(MKS),DIMENSION(:,:),ALLOCATABLE :: array

      IF(ALLOCATED(array))THEN
         swork = array
      ELSE
         ALLOCATE(array(Dimx,Dimy),STAT=info)
         !print*,'status of allocate:',info,',return 0 on success'
      ENDIF
      END SUBROUTINE salloc

      SUBROUTINE dalloc(array,Dimx,Dimy)
      INTEGER :: Dimx,Dimy    
      REAL(MKD),DIMENSION(:,:),ALLOCATABLE :: array

      IF(ALLOCATED(array))THEN
         dwork = array
      ELSE
         ALLOCATE(array(Dimx,Dimy),STAT=info)
         !print*,'status of allocate:',info,',return 0 on success'
      ENDIF
      END SUBROUTINE dalloc

END MODULE oalloc
