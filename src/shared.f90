module uvic_shared
   use fabm_types, only: rk
   real(rk), parameter :: d_per_s     = 1.0_rk/86400.0_rk
   real(rk), parameter :: year_length = 365._rk
   real(rk), parameter :: redptn      = 1._rk/16.0_rk       ! redfield ratio P:N
   real(rk), parameter :: redctn      = 7._rk               ! C:N Redfield ratio
   real(rk), parameter :: redotn      = 10.6_rk             ! O2:N Redfield ratio
end module