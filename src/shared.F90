#include "fabm_driver.h"
module uvic_shared
   use fabm_types
   implicit none
   public
   
   real(rk), parameter :: d_per_s      = 1.0_rk/86400.0_rk
   real(rk), parameter :: year_length  = 365._rk
   real(rk), parameter :: redptn       = 1._rk/16.0_rk       ! redfield ratio P:N
   real(rk), parameter :: redctn       = 7._rk/1000._rk      ! C:N Redfield ratio, includes conversion from mmol m-3 to umol cm-3
   real(rk), parameter :: redotn       = 10.6_rk/1000._rk    ! O:N Redfield ratio, includes conversion from mmol m-3 to umol cm-3
   real(rk), parameter :: redotp       = redotn/redptn       ! O:P Redfield ratio, includes conversion from mmol m-3 to umol cm-3
   real(rk), parameter :: redctp       = redctn/redptn       ! C:P Redfield ratio, includes conversion from mmol m-3 to umol cm-3
   real(rk), parameter :: redntp       = 1._rk/redptn        ! N:P Redfield ratio
   real(rk), parameter :: trcmin       = 5e-12_rk            ! min tracer concentration
   real(rk), parameter :: xconv        = 25.1_rk/3.6e+05_rk  ! xconv is constant to convert piston_vel from cm/hr -> cm/s
   real(rk), parameter :: C2K          = 273.15_rk           ! Celsius to Kelvin constant

   
end module