MODULE glodata

IMPLICIT NONE
  
  INTEGER :: Nx,Ny,Nstep,i,j
  REAL :: grid_size,I_rms,J_rms
  REAL :: core_length,win_length,win_width,bump_length,bump_width
  REAL :: permea_vacuum,relapermea_core,cond_copper,relapermea_copper,permea_cal

  INTEGER :: N_core_x,N_core_y,N_core_wind_x,N_core_wind
  INTEGER :: N_win_width
  INTEGER :: N_wind_x,N_wind_y
  INTEGER :: N_bump_x,N_bump_y

  REAL :: wind_thickness,wind_width
  REAL :: distance_wind_to_core,distance_wind_to_wind

  INTEGER :: N_dis_wind

  REAL,DIMENSION(:,:),ALLOCATABLE :: A_z_old,A_z_new
  REAL,DIMENSION(:,:),ALLOCATABLE :: B_x,B_y
  REAL,DIMENSION(:,:),ALLOCATABLE :: J_z
  INTEGER :: istep

END MODULE glodata
