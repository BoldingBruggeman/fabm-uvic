#include "fabm_driver.h"

module uvic_light
   use fabm_types
   use fabm_particle

   use uvic_shared

   implicit none

   private

   type, extends(type_particle_model), public :: type_uvic_light
      type (type_state_variable_id)                            :: id_
      type (type_dependency_id)                                :: id_dzt, id_phyt, id_depth
      type (type_surface_dependency_id)                        :: id_swr, id_rctheta
      type (type_diagnostic_variable_id)                       :: id_sw_par, id_f1
      
      real(rk)                            :: kc, kw, par
   contains
      procedure :: initialize
      procedure :: do_column
   end type

contains

   subroutine initialize(self, configunit)
      class (type_uvic_light), intent(inout), target :: self
      integer,                 intent(in)            :: configunit

      ! parameters
      call self%get_parameter(self%kc,  'kc',  '1/(mmol m-2)',  'light attenuation by phytoplankton',  default=0.047_rk)
      call self%get_parameter(self%kw,  'kw',  '1/m',           'light attenuation due to water',      default=0.04_rk)
      call self%get_parameter(self%par, 'par', '-',             'PAR fraction of swr',                 default=0.43_rk)
      
      ! variable registrations      
      call self%register_diagnostic_variable(self%id_sw_par, 'sw_par', 'W m-2', 'average layer PAR',                     source=source_do_column)
      call self%register_diagnostic_variable(self%id_f1,     'f1',     '-',     'current layer light extinction coeff.', source=source_do_column)
      
      ! environmental dependencies      
      call self%register_dependency(self%id_rctheta, 'rctheta', '', 'incoming solar angle of incidence')
      call self%register_dependency(self%id_swr,     standard_variables%surface_downwelling_shortwave_flux)
      call self%register_dependency(self%id_depth,   standard_variables%depth)
      call self%register_dependency(self%id_dzt,     standard_variables%cell_thickness)
      call self%register_dependency(self%id_phyt,    type_interior_standard_variable(name='total_phytoplankton',units='mmol N m-3'))

   end subroutine initialize
   
   subroutine do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_uvic_light), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_
      
      real(rk) :: swr, dzt, phyt, rctheta, depth, f1, phin, ztt, sw_par
      _GET_SURFACE_(self%id_swr,swr)
      phin = 0.0_rk
      ztt  = 0.0_rk

      _DOWNWARD_LOOP_BEGIN_
         _GET_(self%id_dzt,dzt)
         _GET_(self%id_depth,depth)
         _GET_(self%id_phyt,phyt)
         _GET_SURFACE_(self%id_rctheta,rctheta)
         swr = swr*exp(-self%kc*phin)
         
         sw_par = self%par*swr*exp(ztt*(self%kw/rctheta))
         
         
         ztt  = -1._rk*(depth+dzt*0.5_rk) !check sign of depth
         phin = phyt*dzt
         f1 = (-self%kw - self%kc*phyt)
         
         _SET_DIAGNOSTIC_(self%id_sw_par,sw_par)
         _SET_DIAGNOSTIC_(self%id_f1,f1)
      _DOWNWARD_LOOP_END_
    end subroutine do_column
end module