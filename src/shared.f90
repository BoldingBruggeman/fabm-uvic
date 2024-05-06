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
   
   ! New stoichiometric model parameters formulated as in Paulmier et al. 2009 BG
   real(rk), parameter :: Paulmier_a   = 1.e3_rk*redctn      
   real(rk), parameter :: Paulmier_z   = 4._rk*1.e3_rk*redotp - 4._rk*Paulmier_a - 8._rk*redntp
   real(rk), parameter :: Paulmier_R0  = Paulmier_a + 0.25_rk*Paulmier_z
   real(rk), parameter :: Paulmier_con = (4._rk/5._rk*Paulmier_R0 + (3._rk/5._rk +1._rk)*redntp) * 1.e-3_rk !conversion between phosphorous sms and alkalinity sms
end module