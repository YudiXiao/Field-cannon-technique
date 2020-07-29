MODULE inij

IMPLICIT NONE

CONTAINS

   SUBROUTINE initij
     USE glodata
     
     I_rms = 1.0
     J_rms = I_rms/(wind_thickness*wind_width)

     J_z(1:Nx+1,1:Ny+1) = 0.0

     J_z(N_core_x+N_core_wind_x+1:N_core_x+N_core_wind_x &
         +N_wind_x+1,N_dis_wind+1:N_dis_wind+1+N_wind_y) &
        = J_rms

   END SUBROUTINE initij

END MODULE inij
