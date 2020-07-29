MODULE readin

IMPLICIT NONE

CONTAINS

   SUBROUTINE readnml(filename)
     USE glodata
     CHARACTER(LEN=*),INTENT(IN) :: filename
     INTEGER :: ios
     NAMELIST /list/grid_size,Nstep,core_length,&
                    win_length,bump_length,bump_width,&
                    wind_thickness,wind_width,&
                    distance_wind_to_wind,relapermea_core

     OPEN(10,FILE=filename)
     READ(10,nml=list,IOSTAT=ios)
     CLOSE(10)

   END SUBROUTINE readnml

END MODULE readin
