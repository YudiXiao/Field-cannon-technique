MODULE odealloc

  IMPLICIT NONE

  INTEGER :: info
  INTEGER,PARAMETER :: MKS = KIND(1.0E0)
  INTEGER,PARAMETER :: MKD = KIND(1.0D0)

  INTERFACE dealloc
      MODULE PROCEDURE sdealloc,ddealloc
  END INTERFACE dealloc

  CONTAINS
      SUBROUTINE sdealloc(array)     
      REAL(MKS),DIMENSION(:,:),ALLOCATABLE :: array

      DEALLOCATE(array,STAT=info)
      !print*,'status of deallocate:',info,',return 0 on success'

      END SUBROUTINE sdealloc

      SUBROUTINE ddealloc(array)    
      REAL(MKD),DIMENSION(:,:),ALLOCATABLE :: array

      DEALLOCATE(array,STAT=info)
      !print*,'status of deallocate:',info,',return 0 on success'

      END SUBROUTINE ddealloc

END MODULE odealloc
