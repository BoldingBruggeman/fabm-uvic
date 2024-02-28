#include "fabm_driver.h"

module uvic_tracer
   use fabm_types
   use fabm_particle

   use uvic_shared

   implicit none

   private

   type, extends(type_particle_model), public :: type_uvic_zooplankton
      type (type_state_variable_id) :: id_
   contains
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class (type_uvic_tracer), intent(inout), target :: self
      integer,                    intent(in)            :: configunit

      logical :: has_carbon, has_nitrogen, has_phosphorus

      call self%get_parameter(has_carbon,     'has_carbon',     '', 'tracer contains carbon',     default=.false.)
      call self%get_parameter(has_nitrogen,   'has_nitrogen',   '', 'tracer contains nitrogen',   default=.false.)
      call self%get_parameter(has_phosphorus, 'has_phosphorus', '', 'tracer contains phosphorus', default=.false.)
      
      if (has_carbon .or. has_nitrogen .or. has_phosphorus) then
         call self%register_state_variable(self%id_c, 'c', 'kmol P m-3', 'concentration in P units', minimum=0.0_rk)
         if (has_carbon)     call self%add_to_aggregate_variable(standard_variables%total_carbon,     self%id_c, scale_factor=rcar * 1e6_rk)
         if (has_nitrogen)   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,   self%id_c, scale_factor=rnit * 1e6_rk)
         if (has_phosphorus) call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_c, scale_factor=1e6_rk)
         if (has_iron)       call self%add_to_aggregate_variable(standard_variables%total_iron,       self%id_c, scale_factor=riron * 1e9_rk)
      end if
   end subroutine initialize
   
   subroutine do(self, _ARGUMENTS_DO_SURFACE_)
      class (type_ihamocc_zooplankton), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_
      
      real(rk) :: 
      
      _LOOP_BEGIN_
      
      _LOOP_END_
   end subroutine do

end module