MODULE procgrid

IMPLICIT NONE

CONTAINS

   SUBROUTINE calgrid
     USE glodata

     distance_wind_to_core = 3.0*(win_length - wind_width)/2.0
     win_width = distance_wind_to_core + wind_thickness & 
                 + distance_wind_to_wind/2.0

     Nx = INT(core_length/(1.0*grid_size))
     N_core_x = INT((core_length - win_length)/grid_size)/2
     N_core_y = N_core_x
     N_win_width = INT(win_width/grid_size)
     Ny = N_win_width + N_core_y
     
     !calculate displacment grids for core-to-winding distance
     N_core_wind_x = INT((win_length - wind_width)/grid_size)/2
     N_core_wind = INT(distance_wind_to_core/grid_size)
     !calculate grids for windings
     N_wind_x = INT(wind_width/(1.0*grid_size))
     N_wind_y = INT(wind_thickness/grid_size)
     !calculate grids for winding distance
     N_dis_wind = INT(distance_wind_to_wind/(2.0*grid_size))
     !calculate grids for bumps
     N_bump_x = INT(bump_length/grid_size)
     N_bump_y = INT(bump_width/(2.0*grid_size))

   END SUBROUTINE calgrid

END MODULE procgrid
