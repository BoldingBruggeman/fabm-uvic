#include "fabm_driver.h"

module uvic_detritus
   use fabm_types
   use fabm_particle

   use uvic_shared

   implicit none

   private

   type, extends(type_particle_model), public :: type_uvic_detritus
      type (type_state_variable_id)                            :: id_det, id_o2, id_oxi, id_phosphorus, id_no3
      type (type_dependency_id)                                :: id_temp, id_depth, id_wd_in
      type (type_diagnostic_variable_id)                       :: id_remi_out, id_wd_out
      
      real(rk)                            :: nud0, wd0, bbio, cbio
      logical                             :: o2_sens, nitrogen, no_temp_sens, no3_sens
   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_bottom
      procedure :: get_vertical_movement
   end type

contains

   subroutine initialize(self, configunit)
      class (type_uvic_detritus), intent(inout), target :: self
      integer,                    intent(in)            :: configunit

      ! parameters
      call self%get_parameter(self%o2_sens,      'o2_sens',      '',       'turn on remineralization oxygen dependence',  default=.true.)
      call self%get_parameter(self%nitrogen,     'nitrogen',     '',       'turn on nitrogen dependency',                 default=.true.)
      if (self%nitrogen) then
          call self%get_parameter(self%no3_sens, 'no3_sens',     '',       'turn on remineralization nitrate dependence', default=.true.)
      endif
      call self%get_parameter(self%no_temp_sens, 'no_temp_sens', '',       'turn off temperature sensitivity',            default=.false.)                                                                                                                               
      if (.not. self%no_temp_sens) then
          call self%get_parameter(self%bbio,     'bbio',         '-',      'temperature dependency intercept',            default=1.066_rk)
          call self%get_parameter(self%cbio,     'cbio',         'degC-1', 'temperature dependency exponent',             default=1.0_rk)
      endif
      call self%get_parameter(self%nud0,         'nud0',         'd-1',    'base remineralization rate',                  default=0.5_rk, scale_factor=d_per_s)
      call self%get_parameter(self%wd0,          'wd0',          'm d-1',  'base detritus sinking speed',                 default=6.0_rk, scale_factor=d_per_s)
      
      ! variable registrations
      call self%register_state_variable(self%id_det, 'detritus', 'mmol N m-3', 'detritus concentration')
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_det)
      call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_det, scale_factor=redptn)
      
      ! register diagnostic variables for output
      call self%register_diagnostic_variable(self%id_remi_out,  'remi',  'mmol m-3 s-1', 'remineralization rate')
      call self%register_diagnostic_variable(self%id_wd_out,    'wd',    'm s-1',        'layer-specific sinking speed')
      
      ! environmental dependencies
      call self%register_state_dependency(self%id_o2,         'o2',      'umol O cm-3', 'oxygen concentration')
      call self%register_state_dependency(self%id_oxi,        'oxi',     'umol cm-3',   'oxidative demand')
      call self%register_state_dependency(self%id_phosphorus, 'p',       'mmol p m-3',  'dissolved inorganic phosphorous concentration')
      call self%register_dependency(self%id_wd_in,            'wd',      'm s-1',       'layer-specific sinking speed')
      
      if (self%nitrogen) then
          call self%register_state_dependency(self%id_no3, 'no3',     'mmol n m-3',  'dissolved inorganic nitrogen concentration')
      endif
      call self%register_dependency(self%id_temp,      standard_variables%temperature)
      call self%register_dependency(self%id_depth,     standard_variables%depth)

   end subroutine initialize
   
   subroutine do(self, _ARGUMENTS_DO_)
      class (type_uvic_detritus), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_
      
      real(rk) :: o2, no3, nud, temp, bct, det, remi, dflag, depth, wd
      
      _LOOP_BEGIN_
         _GET_(self%id_det, det)
         _GET_(self%id_o2, o2)
         if (self%o2_sens) then
             nud = self%nud0*(0.65_rk+0.35_rk*tanh(o2*1000._rk-6._rk)) ! decrease remineralisation rate in oxygen minimum zone
             _GET_(self%id_no3, no3)
             if (self%nitrogen .and. self%no3_sens .and. (no3.lt.0.0_rk)) then !!! The following implies that remineralisation is zero when no3 is zero, even in oxygen replete conditions !!!
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
         dflag = 0.5_rk + sign(0.5_rk,det - trcmin)
         remi = remi*dflag
         
         ! calculate layer-specific sinking rate and save as diagnostic for the do_bottom and get_vertical_movement subroutines
         _GET_(self%id_depth,depth)
         wd = self%wd0+6.0e-2_rk*depth*d_per_s*dflag
         _SET_DIAGNOSTIC_(self%id_wd_out,wd)
         
         ! update state variables
         _ADD_SOURCE_(self%id_det, -remi)
         _ADD_SOURCE_(self%id_phosphorus, redptn*remi)
         _ADD_SOURCE_(self%id_oxi,        -remi*redptn*redotp)
         if (self%nitrogen) then
             _ADD_SOURCE_(self%id_no3, remi)
         endif
         
         ! diagnostic output
         _SET_DIAGNOSTIC_(self%id_remi_out, remi)
      _LOOP_END_
   end subroutine do
    
   subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
      class (type_uvic_detritus), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_
      
      real(rk) :: wd, det
      
      _BOTTOM_LOOP_BEGIN_
         _GET_(self%id_wd_in,wd)
         _GET_(self%id_det, det)
         
         _ADD_BOTTOM_FLUX_(self%id_det,   -wd*det)
         _ADD_BOTTOM_FLUX_(self%id_phosphorus, redptn*wd*det)
         if (self%nitrogen) then
             _ADD_BOTTOM_FLUX_(self%id_no3, wd*det)
         endif
      _BOTTOM_LOOP_END_
   end subroutine do_bottom
   
   subroutine get_vertical_movement(self, _ARGUMENTS_GET_VERTICAL_MOVEMENT_)
      class (type_uvic_detritus), intent(in) :: self
      _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_
      
      real(rk) :: wd, det
      
      _LOOP_BEGIN_
         _GET_(self%id_wd_in,wd)
         _GET_(self%id_det, det)
                  
         _ADD_VERTICAL_VELOCITY_(self%id_det,  -wd)
      _LOOP_END_
   end subroutine get_vertical_movement
end module