!
! file:   kind_param.f90
! author: aixa marã­a rivera rã­os
! created on 9 nov 2010, 1:07 pm
! module to define precision values
module kind_param
    integer, parameter:: single=selected_real_kind(6,37)  !kind(0.0) !single precision: 1 or 4 (processor dependent)
    integer, parameter:: double=selected_real_kind(13,200)  !kind(0.0d0) !double precision: 2 or 8 (processor dependent)
    integer, parameter:: short=selected_int_kind(3) !short precision for integer
    integer,parameter:: long=selected_int_kind(9)   !long precision for integer
end module kind_param
