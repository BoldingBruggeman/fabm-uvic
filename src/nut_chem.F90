#include "fabm_driver.h"

module uvic_nut_chem
   use fabm_types
   use fabm_particle
   use fabm_builtin_depth_integral, only: type_depth_integral
   
   use uvic_shared

   implicit none

   private

   type, extends(type_particle_model), public :: type_uvic_nut_chem
      type (type_state_variable_id)                            :: id_det, id_o2, id_phosphorus, id_no3, id_dic, id_c14!id_det, id_o2, id_phosphorus, id_no3
      type (type_dependency_id)                                :: id_o2_sms, id_phosphorus_sms, id_calc_sms!, id_temp, id_depth, id_wd_in
      type (type_surface_dependency_id)                        :: id_intcalc
      type (type_bottom_dependency_id)                         :: id_bdepth
      type (type_diagnostic_variable_id)                       :: id_o2_dummy, id_phosphorus_dummy, id_calc_dummy
      
      real(rk)                            :: dcaco3
      logical                             :: nitrogen
   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_column
      procedure :: do_bottom
      procedure :: get_vertical_movement
   end type

contains

   subroutine initialize(self, configunit)
      class (type_uvic_nut_chem), intent(inout), target :: self
      integer,                    intent(in)            :: configunit
      
      ! parameters
      call self%get_parameter(self%nitrogen,     'nitrogen',     '',            'turn on nitrogen dependency',                    default=.true.)
      call self%get_parameter(self%dcaco3,       'dcaco3',       'cm',          'calcite remineralization depth',                 default=650000.0_rk, scale_factor=1.e-2)
      call self%get_parameter(self%rstd,         'rstd',         '',            'standard c14/c12 ratio',                         default=1.176e-12_rk)

      call self%register_state_variable(self%id_o2,         'o2',       'umol O cm-3',   'oxygen concentration')
      call self%register_state_variable(self%id_phosphorus, 'p',        'mmol P m-3',    'dissolved inorganic phosphorous concentration')
      call self%register_state_variable(self%id_dic,        'DIC',      'umol C m-3',    'dissolved inorganic carbon concentration')
      call self%register_state_variable(self%id_c14,        'c14',      'umol c14 cm-3', 'carbon 14 concentration')
      if (self%nitrogen) then
          call self%register_state_variable(self%id_no3,    'no3',      'mmol n m-3',  'dissolved inorganic nitrogen concentration')
          call self%register_diagnostic_variable(self%id_deni, 'denitrification',  'mmol m-3 s-1', 'denitrification rate')
      endif


      call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_phosphorus)
      
      ! set up dummy state vars to pick up sources from the npzd model
      call self%register_diagnostic_variable(self%id_o2_dummy, 'o2_dummy','umol cm-3', 'target pool for net respiration',source=source_constant, output=output_none, act_as_state_variable=.true.)
      call self%request_coupling(self%id_o2_sms,'./o2_dummy_sms_tot')
      
      call self%register_diagnostic_variable(self%id_calc_dummy, 'calc_dummy','umol cm-3', 'target pool for calcite production',source=source_constant, output=output_none, act_as_state_variable=.true.)
      call self%request_coupling(self%id_calc_sms,'./calc_dummy_sms_tot')

      call self%register_diagnostic_variable(self%id_phosphorus_dummy, 'phosphorus_dummy','mmol P m-3', 'target pool for p release',source=source_constant, output=output_none, act_as_state_variable=.true.)
      call self%request_coupling(self%id_phosphorus_sms,'./phosphorus_dummy_sms_tot')


      ! set up depth integral of calcite production
      class (type_depth_integral), pointer :: calc_sms_integrator
      allocate(calc_sms_integrator)
      call self%add_child(calc_sms_integrator, 'calc_sms_integrator')
      call calc_sms_integrator%request_coupling('source','../calc_dummy_sms_tot')
      
      ! couple depth-integrated calcite production dependency to calc_sms_integrator 
      call self%register_dependency(self%id_intcalc, 'intcalc', 'umol cm-2', 'depth-integrated calcite production')
      call self%request_coupling(self%id_intcalc, '../calc_sms_integrator/result')
      
      ! environmental dependencies
      call self%register_dependency(self%id_depth,   standard_variables%depth)
      call self%register_dependency(self%id_dzt,     standard_variables%cell_thickness)
      call self%register_dependency(self%id_bdepth,   standard_variables%bottom_depth)

! declare p, calcite, optional n, alk and dic
! get roc of phosphorous and apply to P, alk, and dic
      
   end subroutine initialize
   
   subroutine do(self, _ARGUMENTS_DO_)
      class (type_uvic_nut_chem), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_
      
      real(rk) :: p, p_release, o2, no3, nud, temp, bct, det, o2_demand, no3, fo2, o2_roc, no3flag, deni, calc_prod, dic_roc, alk_roc, 
      
      _LOOP_BEGIN_
         ! apply phosphorus sources to actual phosphorus state var.
         _GET_(self%id_phosphorus_sms,    p_release)
         _ADD_SOURCE_(self%id_phosphorus, p_release) !apply phosphorus sms to real phosphorus state variable
         
         ! calculate and apply carbon sources from p_release and calcite production
         _GET_(self%id_calc_sms,    calc_prod)
         
         dic_roc = p_release*redctp-calc_prod
         
         _ADD_SOURCE_(self%id_dic, dic_roc)
         
         alk_roc = -p_release*redntp*1.e-3_rk-calc_prod*2._rk
         
         ! calculate and apply oxygen use and denitrification
         _GET_(self%id_o2, o2)
         _GET_(self%id_o2_sms, o2_demand)
         
         fo2 = 0.5_rk*tanh(o2*1000._rk - 5._rk) ! limit oxygen consumption below concentrations of 5umol/kg as recommended in OCMIP
         
         so2 = (p_release*redotp+o2_demand)*(0.5_rk + fo2)
         
         _ADD_SOURCE_(self%id_o2,         so2)
         if (self%nitrogen) then
             _GET_(self%id_no3, no3)
             no3flag = 0.5_rk+sign(0.5_rk,no3-trcmin)
             deni = max(0._rk, 800._rk*no3flag*so2*(0.5_rk - fo2))
             _ADD_SOURCE_(self%id_no3,          -deni)
             _SET_DIAGNOSTIC_(self%id_deni,      deni)
             
             !calculate alkalinity change from p_release
             Paulmier_a   = 1.e3_rk*redctn      
             Paulmier_z   = 4._rk*1.e3_rk*redotp - 4._rk*Paulmier_a - 8._rk*redntp
             Paulmier_R0  = Paulmier_a + 0.25_rk*Paulmier_z
                  
             alk_roc = alk_roc + no3flag*p_release*(0.5_rk-fo2)*(4._rk/5._rk*Paulmier_R0 + (3._rk/5._rk +1._rk)*redntp) * 1.e-3_rk
         endif
         _ADD_SOURCE_(self%id_alk, alk_roc)
         
         ! calculate c14 production and decay
         _GET_(self%id_c14, c14)
         c14_roc = dic_roc*self%rstd - 3.836e-12_rk*c14
      _LOOP_END_
   end subroutine do
   
   subroutine do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_uvic_nut_chem), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_
      
      real(rk) :: dzt, depth, rcak, rcab, zw, zw_prev, intcalc, bdepth
      _GET_SURFACE_(self%id_intcalc,intcalc)
      intcalc = intcalc*100._rk              ! we multiply with [cm m-1], to convert back to cm-2 (because of unit mismatch, type_depth_integral returns [umol cm-3 m])
      _GET_BOTTOM_(self%id_bdepth,bdepth

      rcak = 0.0_rk 
      rcab = 0.0_rk
      _DOWNWARD_LOOP_BEGIN_
         _GET_(self%id_dzt,dzt)
         _GET_(self%id_depth,depth)

         zw = depth+0.5_rk*dzt
         
         if (rcak == 0.0_rk) then
             rcak = -(exp(-zw/self%dcaco3)-1._rk)/dzt
             rcab = -1._rk/dzt
         else
             rcak = -(exp(-zw/self%dcaco3))/dzt
             rcab = (exp(-zw_prev/self%dcaco3))/dzt
         endif
         
         zw_prev = zw
         
         if (zw < bdepth) then
             dic_roc = intcalc*rcak
             alk_roc = intcalc*rcak*2._rk
         else
             dic_roc = intcalc*rcab
             alk_roc = intcalc*rcab*2._rk
         endif
         
         _ADD_SOURCE_(self%id_dic, dic_roc) !NIC: NOT FINISHED we still need to interact with the sediment (line 742 of tracer.f90) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         _ADD_SOURCE_(self%id_alk, alk_roc)
      _DOWNWARD_LOOP_END_
   end subroutine do_column
end module
    