#include "fabm_driver.h"

module uvic_solar
   use fabm_types
   use uvic_shared
   
   implicit none

   private

   type, extends(type_base_model), public :: type_uvic_solar
      type (type_global_dependency_id)           :: id_doy
      type (type_surface_dependency_id)          :: id_latitude
      type (type_surface_diagnostic_variable_id) :: id_rctheta, id_dayfrac
   contains
      procedure :: initialize
      procedure :: do_surface
   end type

contains

   subroutine initialize(self, configunit)
      class (type_uvic_solar), intent(inout), target :: self
      integer,                       intent(in)            :: configunit

      call self%register_dependency(self%id_doy, standard_variables%number_of_days_since_start_of_the_year)
      call self%register_dependency(self%id_latitude, standard_variables%latitude)
      call self%register_diagnostic_variable(self%id_rctheta, 'rctheta', '', 'incoming solar angle of incidence')
      call self%register_diagnostic_variable(self%id_dayfrac, 'dayfrac', '-', 'daylength fraction')
   end subroutine initialize

   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_uvic_solar), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_
      
      real(rk)            :: doy, relyr, declin, rctheta, dayfrac, lat
      real(rk), parameter :: pi      = 3.141592653589793_rk           !: pi
      real(rk), parameter :: radians = 180._rk/pi                     !: conversion from radians to degrees
      
      _GET_GLOBAL_(self%id_doy, doy)
      
      relyr = doy/year_length
      declin = sin((relyr - 0.22_rk)*2._rk*pi)*0.4_rk   ! declination
      
      _SURFACE_LOOP_BEGIN_
         _GET_SURFACE_(self%id_latitude, lat)
         
         rctheta = max(-1.5_rk, min(1.5_rk, lat/radians - declin))
         rctheta = sqrt(1._rk - (1._rk - cos(rctheta)**2._rk)/1.33_rk**2._rk)! NIC: Move this line to the light/phytoplankton model? kw is pure water light extinction coefficient
         dayfrac = min( 1._rk, -tan(lat/radians)*tan(declin))
         dayfrac = max(1e-12_rk, acos(max(-1._rk, dayfrac))/pi)
         
         _SET_SURFACE_DIAGNOSTIC_(self%id_rctheta, rctheta)
         _SET_SURFACE_DIAGNOSTIC_(self%id_dayfrac, dayfrac)
      _SURFACE_LOOP_END_
   end subroutine

end module