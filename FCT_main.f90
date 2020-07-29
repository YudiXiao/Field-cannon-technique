PROGRAM FCT

USE glodata
USE readin,ONLY:readnml
USE procgrid,ONLY:calgrid
USE oalloc,ONLY:alloc
USE odealloc,ONLY:dealloc
USE inij,ONLY:initij

IMPLICIT NONE

CALL readnml('external_input.txt') ! read input data (geometry)
permea_vacuum = 4.0*3.14/10000000.0
cond_copper = 5.96*10000000.0
relapermea_copper = 0.999994


CALL calgrid() ! discretize the geometry

CALL alloc(A_z_old,Nx+1,Ny+1)
CALL alloc(A_z_new,Nx+1,Ny+1)
CALL alloc(J_z,Nx+1,Ny+1)
CALL alloc(B_x,Nx+1,Ny+1)
CALL alloc(B_y,Nx+1,Ny+1)

CALL initij() ! initialize J field, which is uniform under static field

! initial conditions
A_z_old(2:Nx,2:Ny) = 0.0
! end of initialization

! boundary conditions
! Dirichlet
A_z_old(1,1:Ny+1) = 0.0
A_z_old(2:Nx,Ny+1) = 0.0
A_z_old(2:Nx,1) = 0.0
A_z_old(Nx+1,1:Ny+1) = 0.0
! end of assigning boundary conditions

A_z_new = A_z_old

! perform time stepping
DO istep = 1,Nstep

  DO j = 2,Ny-1
     DO i = 2,Nx-1

        IF((i.EQ.N_core_x+N_core_wind_x+1).AND.(j.EQ.N_dis_wind+1))THEN! calculate interface #4, copper air
           A_z_new(i,j) = (relapermea_copper/(3.0*relapermea_copper+1.0))&
                          *(A_z_old(i-1,j) + A_z_old(i,j-1))&
                          +((relapermea_copper+1.0)/(2.0*(3.0*relapermea_copper+1.0)))&
                          *(A_z_old(i+1,j) + A_z_old(i,j+1))&
                          +((relapermea_copper*permea_vacuum)&
                          /(4.0*(3.0*relapermea_copper+1.0)))&
                          *grid_size*grid_size*J_z(i,j)

        ELSEIF((i.EQ.N_core_x+1).AND.(j.EQ.N_win_width+1))THEN! calculate interface #3, air core
           A_z_new(i,j) = (1.0/(3.0+relapermea_core))&
                          *(A_z_old(i,j+1) + A_z_old(i-1,j))&
                          +((1.0+relapermea_core)/(2.0*(3.0+relapermea_core)))&
                          *(A_z_old(i,j-1) + A_z_old(i+1,j))

        ELSEIF((i.EQ.N_core_x+N_core_wind_x+1).AND.(j.EQ.N_dis_wind+N_wind_y+1))THEN! calculate interface #3, copper air
           A_z_new(i,j) = (relapermea_copper/(3.0*relapermea_copper+1.0))&
                          *(A_z_old(i,j+1) + A_z_old(i-1,j))&
                          +((relapermea_copper+1.0)/(2.0*(3.0*relapermea_copper+1.0)))&
                          *(A_z_old(i,j-1) + A_z_old(i+1,j))&
                          +((relapermea_copper*permea_vacuum)&
                          /(4.0*(3.0*relapermea_copper+1.0)))&
                          *grid_size*grid_size*J_z(i,j)

        ELSEIF((i.EQ.N_core_x+N_bump_x+1).AND.(j.EQ.N_bump_y+1))THEN! calculate interface #5, core air
           A_z_new(i,j) = (relapermea_core/(3.0*relapermea_core+1.0))&
                          *(A_z_old(i+1,j) + A_z_old(i,j+1))&
                          +((1.0+relapermea_core)/(2.0*(3.0*relapermea_core+1)))&
                          *(A_z_old(i-1,j) + A_z_old(i,j-1))

        ELSEIF((i.EQ.N_core_x+N_core_wind_x+N_wind_x+1).AND.(j.EQ.N_dis_wind+1))THEN! calculate interface #6, copper air
           A_z_new(i,j) = (relapermea_copper/(3.0*relapermea_copper+1.0))&
                          *(A_z_old(i+1,j) + A_z_old(i,j-1))&
                          +((relapermea_copper+1.0)/(2.0*(3.0*relapermea_copper+1.0)))&
                          *(A_z_old(i-1,j) + A_z_old(i,j+1))&
                          +((relapermea_copper*permea_vacuum)&
                          /(4.0*(3.0*relapermea_copper+1.0)))&
                          *grid_size*grid_size*J_z(i,j)

        ELSEIF((i.EQ.Nx-N_core_x+1).AND.(j.EQ.N_win_width+1))THEN! calculate interface #5, air core
           A_z_new(i,j) = (1.0/(3.0+relapermea_core))&
                          *(A_z_old(i+1,j) + A_z_old(i,j+1))&
                          +((1.0+relapermea_core)/(2.0*(3.0+relapermea_core)))&
                          *(A_z_old(i-1,j) + A_z_old(i,j-1))

        ELSEIF((i.EQ.N_core_x+N_core_wind_x+N_wind_x+1).AND.(j.EQ.N_dis_wind+N_wind_y+1))THEN! calculate interface #5, copper air
           A_z_new(i,j) = (1.0/(3.0+relapermea_core))&
                          *(A_z_old(i+1,j) + A_z_old(i,j+1))&
                          +((1.0+relapermea_core)/(2.0*(3.0+relapermea_core)))&
                          *(A_z_old(i-1,j) + A_z_old(i,j-1))&
                          +((relapermea_copper*permea_vacuum)&
                          /(4.0*(3.0*relapermea_copper+1.0)))&
                          *grid_size*grid_size*J_z(i,j)

        ELSEIF((i.EQ.Nx-N_core_x-N_bump_x+1).AND.(j.EQ.N_bump_y+1))THEN! calculate interface #3, core air
           A_z_new(i,j) = (relapermea_core/(3.0*relapermea_core+1.0))&
                          *(A_z_old(i,j+1) + A_z_old(i-1,j))&
                          +((relapermea_core+1.0)/(2.0*(3.0*relapermea_core+1.0)))&
                          *(A_z_old(i,j-1) + A_z_old(i+1,j))

        ELSEIF((i.EQ.N_core_x+1).AND.(j.GE.N_bump_y+2).AND.(j.LE.N_win_width))THEN! calculate A at interface #1, core and air
           A_z_new(i,j) = 0.25*(((2.0*relapermea_core)/(relapermea_core+1.0))&
                          *A_z_old(i+1,j) + A_z_old(i,j+1)&
                          + (2.0/(relapermea_core+1.0))*A_z_old(i-1,j)&
                          + A_z_old(i,j-1))

        ELSEIF((i.EQ.N_core_x+N_bump_x+1).AND.(j.GE.2).AND.(j.LE.N_bump_y))THEN! calculate A at interface #1, core and air
           A_z_new(i,j) = 0.25*(((2.0*relapermea_core)/(relapermea_core+1.0))&
                          *A_z_old(i+1,j) + A_z_old(i,j+1)&
                          + (2.0/(relapermea_core+1.0))*A_z_old(i-1,j)&
                          + A_z_old(i,j-1))

        ELSEIF((i.EQ.Nx-N_core_x+1).AND.(j.GE.N_bump_y+2).AND.(j.LE.N_win_width))THEN! calculate A at interface #1, air and core
           A_z_new(i,j) = 0.25*(((2.0)/(1.0+relapermea_core))&
                          *A_z_old(i+1,j) + A_z_old(i,j+1)&
                          + (2.0*relapermea_core/(1.0+relapermea_core))*A_z_old(i-1,j)&
                          + A_z_old(i,j-1))

        ELSEIF((i.EQ.Nx-N_core_x-N_bump_x+1).AND.(j.GE.2).AND.(j.LE.N_bump_y))THEN! calculate A at interface #1, air and core
           A_z_new(i,j) = 0.25*(((2.0)/(1.0+relapermea_core))&
                          *A_z_old(i+1,j) + A_z_old(i,j+1)&
                          + (2.0*relapermea_core/(1.0+relapermea_core))*A_z_old(i-1,j)&
                          + A_z_old(i,j-1))

        ELSEIF((i.EQ.N_core_x+N_core_wind_x+1).AND.(j.GE.N_dis_wind+2)&
              .AND.(j.LE.N_dis_wind+N_wind_y))THEN! calculate A at interface #1, air and copper
           A_z_new(i,j) = 0.25*((2.0/(1.0+relapermea_copper))*A_z_old(i+1,j)&
                          + A_z_old(i,j+1) + ((2.0*relapermea_copper)&
                          /(1.0+relapermea_copper))&
                          *A_z_old(i-1,j) + A_z_old(i,j-1)&
                          + (relapermea_copper*permea_vacuum&
                          /(1.0 + relapermea_copper))&
                          *grid_size*grid_size*J_z(i,j))

        ELSEIF((i.EQ.N_core_x+N_core_wind_x+N_wind_x+1).AND.(j.GE.N_dis_wind+2)&
              .AND.(j.LE.N_dis_wind+N_wind_y))THEN! calculate A at interface #1, copper and air
           A_z_new(i,j) = 0.25*((2.0*relapermea_copper/(relapermea_copper+1.0))&
                          *A_z_old(i+1,j)&
                          + A_z_old(i,j+1) + ((2.0)&
                          /(1.0+relapermea_copper))&
                          *A_z_old(i-1,j) + A_z_old(i,j-1)&
                          + (relapermea_copper*permea_vacuum&
                          /(1.0 + relapermea_copper))&
                          *grid_size*grid_size*J_z(i,j))

        ELSEIF((i.GE.N_core_x+2).AND.(i.LE.Nx-N_core_x).AND.(j.EQ.N_win_width+1))THEN! calculate A at interface #2, core and air
           A_z_new(i,j) = 0.25*(A_z_old(i+1,j) + (2.0/(relapermea_core+1))&
                           *A_z_old(i,j+1) + A_z_old(i-1,j)&
                           + ((2.0*relapermea_core)/(relapermea_core+1.0))&
                           *A_z_old(i,j-1))

        ELSEIF((i.GE.N_core_x+N_core_wind_x+2).AND.(i.LE.N_core_x+N_core_wind_x+N_wind_x).AND.&
               (j.EQ.N_dis_wind+N_wind_y+1))THEN! calculate A at interface #2, air and copper
           A_z_new(i,j) = 0.25*(A_z_old(i+1,j)&
                           + ((2.0*relapermea_copper)/(1.0 + relapermea_copper))&
                           *A_z_old(i,j+1) + A_z_old(i-1,j) + (2.0/(1.0 + &
                           relapermea_copper))*A_z_old(i,j-1) + (relapermea_copper&
                           *permea_vacuum/(1.0 + relapermea_copper))&
                           *grid_size*grid_size*J_z(i,j))

        ELSEIF((i.GE.N_core_x+2).AND.(i.LE.N_core_x+N_bump_x).AND.&
               (j.EQ.N_bump_y+1))THEN! calculate A at interface #2, air and core
           A_z_new(i,j) = 0.25*(A_z_old(i+1,j)&
                           + ((2.0*relapermea_core)/(1.0 + relapermea_core))&
                           *A_z_old(i,j+1) + A_z_old(i-1,j) + (2.0/(1.0 + &
                           relapermea_core))*A_z_old(i,j-1))

        ELSEIF((i.GE.Nx-N_core_x-N_bump_x+2).AND.(i.LE.Nx-N_core_x).AND.&
               (j.EQ.N_bump_y+1))THEN! calculate A at interface #2, air and core
           A_z_new(i,j) = 0.25*(A_z_old(i+1,j)&
                           + ((2.0*relapermea_core)/(1.0 + relapermea_core))&
                           *A_z_old(i,j+1) + A_z_old(i-1,j) + (2.0/(1.0 + &
                           relapermea_core))*A_z_old(i,j-1))

        ELSEIF((i.GE.N_core_x+N_core_wind_x+2).AND.(i.LE.N_core_x+N_core_wind_x+N_wind_x).AND.&
               (j.EQ.N_dis_wind+1))THEN! calculate A at interface #2, copper and air
           A_z_new(i,j) = 0.25*(A_z_old(i+1,j) + (2.0/(relapermea_copper+1.0))&
                          *A_z_old(i,j+1) + A_z_old(i-1,j)&
                           + ((2.0*relapermea_copper)/(relapermea_copper+1.0))&
                           *A_z_old(i,j-1) + (relapermea_copper*permea_vacuum/&
                           (relapermea_copper+1.0))&
                           *grid_size*grid_size*J_z(i,j))

        ELSE! calculate A within the same material
            IF((i.LE.N_core_x).OR.(j.GE.N_win_width+2).OR.(i.GE.Nx-N_core_x+2))THEN
              permea_cal = permea_vacuum*relapermea_core
            ELSEIF((i.GE.N_core_x+1).AND.(i.LE.N_core_x+N_bump_x).AND.&
                  (j.GE.2).AND.(j.LE.N_bump_y))THEN!
              permea_cal = permea_vacuum*relapermea_core
            ELSEIF((i.GE.Nx-N_core_x-N_bump_x+2).AND.(i.LE.Nx-N_core_x+1).AND.&
                  (j.GE.2).AND.(j.LE.N_bump_y))THEN!
              permea_cal = permea_vacuum*relapermea_core
            ELSEIF((i.GE.N_core_x+N_core_wind_x+2).AND.(i.LE.N_core_x+N_core_wind_x+N_wind_x)&
                   .AND.(j.GE.N_dis_wind+2).AND.(j.LE.N_dis_wind+N_wind_y))THEN
              permea_cal = permea_vacuum*relapermea_copper
            ELSE
              permea_cal = permea_vacuum
            ENDIF
            A_z_new(i,j) = 0.25*(A_z_old(i+1,j) + A_z_old(i,j+1) &
                                 + A_z_old(i-1,j) + A_z_old(i,j-1) &
                                 + grid_size*grid_size*permea_cal*J_z(i,j))
        ENDIF

     ENDDO
  ENDDO

  A_z_old = A_z_new

  IF(istep.EQ.Nstep)THEN
    OPEN(10,FILE = 'A_old.dat')
    DO j = 1,Ny+1
       DO i = 1,Nx+1
          WRITE(10,'(3E12.4)') REAL(i-1)*grid_size,REAL(j-1)*grid_size,A_z_old(i,j)
       ENDDO
       WRITE(10,'(A)')
    ENDDO
    CLOSE(10)
  ENDIF

ENDDO
! end of do loop for time stepping

PRINT*,'Nx = ',Nx
PRINT*,'Ny = ',Ny
PRINT*,'istep = ',istep
PRINT*,'N_core_x = ',N_core_x
PRINT*,'N_bump_x = ',N_bump_x
PRINT*,'N_bump_y = ',N_bump_y
PRINT*,'N_dis_wind = ',N_dis_wind
PRINT*,'N_core_wind_x = ',N_core_wind_x
PRINT*,'N_wind_x = ',N_wind_x
PRINT*,'N_wind_y = ',N_wind_y

!CALL dealloc(A_z_old)
!CALL dealloc(A_z_new)
!CALL dealloc(J_z)

OPEN(10,FILE = 'A_new.dat')
DO j = 1,Ny+1
   DO i = 1,Nx+1
      WRITE(10,'(3E12.4)') REAL(i-1)*grid_size,REAL(j-1)*grid_size,A_z_old(i,j)
   ENDDO
   WRITE(10,'(A)')
ENDDO
CLOSE(10)

! calculate B_x and B_y
DO j = 1,Ny+1
     DO i = 1,Nx+1
        IF((i.EQ.1).AND.(j.GE.1).AND.(j.LE.Ny))THEN
          B_x(i,j) = (A_z_old(i,j+1)-A_z_old(i,j))/(grid_size)
          B_y(i,j) = -(A_z_old(i+1,j)-A_z_old(i,j))/(grid_size)
        ELSEIF((j.EQ.Ny+1).AND.(i.GE.1).AND.(i.LE.Nx))THEN
          B_x(i,j) = (A_z_old(i,j)-A_z_old(i,j-1))/(grid_size)
          B_y(i,j) = -(A_z_old(i+1,j)-A_z_old(i,j))/(grid_size)
        ELSEIF((i.EQ.Nx+1).AND.(j.GE.2).AND.(j.LE.Ny+1))THEN
          B_x(i,j) = (A_z_old(i,j)-A_z_old(i,j-1))/(grid_size)
          B_y(i,j) = -(A_z_old(i,j)-A_z_old(i-1,j))/(grid_size)
        ELSEIF((j.EQ.1).AND.(i.GE.2).AND.(i.LE.Nx+1))THEN
          B_x(i,j) = (A_z_old(i,j+1)-A_z_old(i,j))/(grid_size)
          B_y(i,j) = -(A_z_old(i,j)-A_z_old(i-1,j))/(grid_size)
        ELSE
          B_x(i,j) = (A_z_old(i,j+1)-A_z_old(i,j-1))/(2.0*grid_size)
          B_y(i,j) = -(A_z_old(i+1,j)-A_z_old(i-1,j))/(2.0*grid_size)
        ENDIF
     ENDDO
ENDDO
! end of calculating B_x and B_y

OPEN(10,FILE = 'B_x.dat')
DO j = 1,Ny+1
   DO i = 1,Nx+1
      WRITE(10,'(3E12.4)') REAL(i-1)*grid_size,REAL(j-1)*grid_size,B_x(i,j)
   ENDDO
   WRITE(10,'(A)')
ENDDO
CLOSE(10)

OPEN(10,FILE = 'B_y.dat')
DO j = 1,Ny+1
   DO i = 1,Nx+1
      WRITE(10,'(3E12.4)') REAL(i-1)*grid_size,REAL(j-1)*grid_size,B_y(i,j)
   ENDDO
   WRITE(10,'(A)')
ENDDO
CLOSE(10)

END PROGRAM FCT
