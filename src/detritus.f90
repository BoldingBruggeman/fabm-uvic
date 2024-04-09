#include "fabm_driver.h"

module uvic_detritus
   use fabm_types
   use fabm_particle

   use uvic_shared

   implicit none

   private

   type, extends(type_particle_model), public :: type_uvic_detritus
      type (type_state_variable_id)                            :: id_det, id_o2, id_phosphorus, id_no3
      type (type_dependency_id)                                :: id_temp, id_depth
      
      real(rk)                            :: nu0, wd0
      logical                             :: o2_sens, nitrogen, no_temp_sens
   contains
      procedure :: initialize
      procedure :: do
      procedure :: get_vertical_movement
   end type

contains

   subroutine initialize(self, configunit)
      class (type_uvic_detritus), intent(inout), target :: self
      integer,                    intent(in)            :: configunit

      ! parameters
      call self%get_parameter(self%o2_sens,      '02_sens',      '',       'turn on remineralization oxygen dependence',  default=.true.)
      call self%get_parameter(self%nitrogen,     'nitrogen',     '',       'turn on nitrogen dependency',                 default=.true.)
      if (self%nitrogen) then
          call self%get_parameter(self%no3_sens, 'n03_sens',     '',       'turn on remineralization nitrate dependence', default=.true.)
      endif
      call self%get_parameter(self%no_temp_sens, 'no_temp_sens', '',       'turn off temperature sensitivity',            default=.false.)                                                                                                                               
      if (.not. self%no_temp_sens) then
          call self%get_parameter(self%bbio,     'bbio',         '-',      'temperature dependency intercept',            default=1.066_rk)
          call self%get_parameter(self%cbio,     'cbio',         'degC-1', 'temperature dependency exponent',             default=1.0_rk)
      endif
      call self%get_parameter(self%nud0,         'nud0',         'd-1',    'base remineralization rate',                  default=0.5_rk, scale_factor=d_per_s)
      call self%get_parameter(self%wd0,          'wd0',          'm d-1',  'base detritud sinking speed',                 default=6.0_rk, scale_factor=d_per_s)
      
      ! variable registrations
      call self%register_state_variable(self%id_det, 'detritus', 'mmol N m-3', 'detritus concentration', minimum=0.0_rk, )
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_det)
      call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_det, scale_factor=redptn)
      
      ! environmental dependencies
      call self%register_state_dependency(self%id_o2, 'o2',     'mmol O m-3',  'oxygen concentration')
      call self%register_state_dependency(self%id_phosphorus, 'p',     'mmol p m-3',  'dissolved inorganic phosphorous concentration')
      if (self%nitrogen) then
          call self%register_state_dependency(self%id_no3, 'no3',     'mmol n m-3',  'dissolved inorganic nitrogen concentration')
      endif
      call self%register_dependency(self%id_temp,      standard_variables%temperature)
      call self%register_dependency(self%id_depth,     standard_variables%depth)

   end subroutine initialize
   
   subroutine do(self, _ARGUMENTS_DO_)
      class (type_uvic_detritus), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_
      
      real(rk) :: o2, no3, nud, temp, bct, det
      
      _LOOP_BEGIN_
         _GET_(self%id_detritus, det)
         if (self%o2_sens) then
             _GET_(self%id_o2, o2)
             nud = self%nud0*(0.65_rk+0.35_rk*tanh(o2*1000._rk-6._rk)) ! decrease remineralisation rate in oxygen minimum zone
             _GET_(self%id_no3, no3)
             if (self%nitrogen .and. no3_sens .and. (no3.lt.0.0_rk)) then !!! The following implies that remineralisation is zero when no3 is zero, even in oxygen replete conditions !!!
                 nud = 0.0_rk
             endif
         else
             nud=self%nud0
         endif
         
         if (self%no_temp_sens) then
             remi = nud*2.99_rk*det
         else
             _GET_(self%id_temp, temp)
             bct = self%bbio**(self%cbio*temp)
             remi = nud*bct*det
         endif
         
         _ADD_SOURCE_(self%id_det, -remi)
         _ADD_SOURCE_(self%id_phosphorus, redptn*remi)
         if (self%nitrogen)
             _ADD_SOURCE_(self%id_no3, remi)
         endif
      _LOOP_END_
    end subroutine do
    
    subroutine get_vertical_movement(self, _ARGUMENTS_GET_VERTICAL_MOVEMENT_)
      class (type_uvic_detritus), intent(in) :: self
      _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_
      
      real(rk) :: wd, depth
      
      _LOOP_BEGIN_
         _GET_(self%id_depth,depth)
         wd = self%wd0+6.0e-2_rk*depth*d_per_s 
         
         _ADD_VERTICAL_VELOCITY_(self%id_det,  -wd)
      _LOOP_END_
   end subroutine get_vertical_movement
end module