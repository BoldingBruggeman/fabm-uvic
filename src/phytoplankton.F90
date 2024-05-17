#include "fabm_driver.h"

module uvic_phytoplankton
   use fabm_types
   use fabm_particle

   use uvic_shared

   implicit none

   private

   type, extends(type_particle_model), public :: type_uvic_phytoplankton
      type (type_state_variable_id)         :: id_phyt, id_phosphorus, id_no3, id_det, id_o2, id_oxi, id_alk, id_calc
      type (type_surface_dependency_id)     :: id_fe, id_dayfrac
      type (type_dependency_id)             :: id_temp, id_depth, id_sw_par, id_f1, id_dzt
      type (type_horizontal_dependency_id)  :: id_latitude
      type (type_diagnostic_variable_id)    :: id_morp_out, id_npp_out, id_morpt_out, id_avej_out, id_no3P_out, id_po4P_out
      
      real(rk)                            :: abio, bbio, cbio, kfe, alpha, k1n, nup, nupt0, jdiar, capr
      logical                             :: no_temp_sens, nitrogen, no3_sens, fe_limitation, so_fe_fer, cdom_attenuation, extra_diags, nfix
   contains
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class (type_uvic_phytoplankton), intent(inout), target :: self
      integer,                    intent(in)            :: configunit

      ! parameters
      call self%get_parameter(self%alpha,        'alpha',        '(W/m2)-1 d-1','Initial slope P-I curve',                        default=0.1_rk, scale_factor=d_per_s)
      call self%get_parameter(self%k1n,          'k1n',          'mmol N m-3',  'Half saturation constant for N uptake',          default=0.7_rk)
      call self%get_parameter(self%nup,          'nup',          'd-1',         'Specific mortality rate',                        default=0.025_rk, scale_factor=d_per_s)
      call self%get_parameter(self%nupt0,        'nupt0',        'd-1',         'temp. dependent specific mortality rate',        default=0.02_rk, scale_factor=d_per_s)
      call self%get_parameter(self%capr,         'capr',         '',            'carbonate to carbon production ratio',           default=0.018_rk)
      call self%get_parameter(self%nitrogen,     'nitrogen',     '',            'turn on nitrogen dependency',                    default=.true.)
      if (self%nitrogen) then                                                   
          call self%get_parameter(self%nfix,     'nfix',         '',            'turn on n-fixation',                             default=.false.)
          if (self%nfix) then                                                   
              call self%get_parameter(self%jdiar,'jdiar',        '',            'factor reducing the growth rate of diazotrophs', default=0.5_rk)
          endif
      endif                                                                     
      call self%get_parameter(self%extra_diags,  'extra_diags',  '',            'turn on extra diagnostic variables',             default=.false.)    
      call self%get_parameter(self%abio,         'abio',         'd-1',         'maximum growth rate',                            default=0.18_rk, scale_factor=d_per_s)
      call self%get_parameter(self%no_temp_sens, 'no_temp_sens', '',            'turn off temperature sensitivity',               default=.false.)                                                                                                                               
      if (.not. self%no_temp_sens) then                                                                                           
          call self%get_parameter(self%bbio,     'bbio',         '-',           'temperature dependency intercept',               default=1.066_rk)
          call self%get_parameter(self%cbio,     'cbio',         'degC-1',      'temperature dependency exponent',                default=1.0_rk)
      endif                                                                                                                       
      call self%get_parameter(self%fe_limitation,'fe_limitation','',            'turn on iron limitation',                        default=.true.)                                                                                                                               
      if (self%fe_limitation) then                                                                                                
          call self%get_parameter(self%kfe,      'kfe',          'mmol Fe m-3', 'iron half saturation constant',                  default=0.1_rk)
          call self%get_parameter(self%so_fe_fer,'so_fe_fer',    '',            'iron fertilization in the southern ocean',       default=.false.)
      endif                                                                                                                       
      call self%get_parameter(self%cdom_attenuation,'cdom_attenuation','',      'cdom_attenuation',                               default=.false.)
      
      
      ! variable registrations
      call self%register_state_variable(self%id_phyt, 'phytoplankton', 'mmol N m-3', 'phytoplankton concentration', minimum=0.0_rk)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_phyt)
      call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_phyt, scale_factor=redptn)
      call self%add_to_aggregate_variable(type_interior_standard_variable(name='total_phytoplankton',units='mmol N m-3',aggregate_variable=.true.), self%id_phyt)
      
      ! environmental dependencies
      call self%register_dependency(self%id_sw_par,           'sw_par',           'W m-2',       'average layer PAR')
      call self%register_dependency(self%id_f1,               'f1',               '-',           'current layer light extinction coeff.')
      call self%register_dependency(self%id_dayfrac,          'dayfrac',          '-',           'daylength fraction')
      call self%register_state_dependency(self%id_phosphorus, 'p',                'mmol p m-3',  'dissolved inorganic phosphorous concentration')
      call self%register_state_dependency(self%id_det,        'detritus',         'mmol N m-3',  'detritus concentration')
      call self%register_state_dependency(self%id_o2,         'oxygen',           'umol cm-3',   'oxygen concentration')
      call self%register_state_dependency(self%id_oxi,        'oxi',              'umol cm-3',   'oxidative demand')
      call self%register_state_dependency(self%id_alk,        'alkalinity',       'umol cm-3',   'ocean alkalinity')
      call self%register_state_dependency(self%id_calc,       'calcite_detritus', 'umol C m-3',  'detrital calcite concentration')
      if (self%nitrogen) then                                                            
          call self%register_state_dependency(self%id_no3,    'no3',      'mmol n m-3',  'dissolved inorganic nitrogen concentration')
      endif                                                                              
      if (self%fe_limitation) then                                                       
          call self%register_surface_dependency(self%id_fe,   'fe',       'mmol Fe m-3', 'photic zone iron concentration')
      endif
      call self%register_dependency(self%id_temp,      standard_variables%temperature)
      call self%register_dependency(self%id_depth,     standard_variables%depth)
      call self%register_dependency(self%id_dzt,       standard_variables%cell_thickness)
      call self%register_dependency(self%id_latitude,  standard_variables%latitude)
      
      ! register diagnostic variables for output
      call self%register_diagnostic_variable(self%id_morp_out,  'morp_out',  'mmol m-3 s-1', 'quadratic mortality of phytoplankton')
      call self%register_diagnostic_variable(self%id_npp_out,   'npp_out',   'mmol m-3 s-1', 'net primary production')
      call self%register_diagnostic_variable(self%id_morpt_out, 'morpt_out', 'mmol m-3 s-1', 'specific mortality of phytoplankton')
      if (self%extra_diags) then
          call self%register_diagnostic_variable(self%id_avej_out, 'avej_out', 's-1', 'light-depend phyt. growth rate')
          call self%register_diagnostic_variable(self%id_no3P_out, 'no3P_out', 's-1', 'no3 depend. phyt growth rate')
          call self%register_diagnostic_variable(self%id_po4P_out, 'po4P_out', 's-1', 'po4 depend. phyt growth rate')
      endif
   end subroutine initialize
   
   subroutine do(self, _ARGUMENTS_DO_SURFACE_)
      class (type_uvic_phytoplankton), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_
      
      real(rk) :: depth, lat, temp, felimit, fe_conc, jmax, sw_par, gl, f1, gd, u1, u2, phi1, phi2
      real(rk) :: dayfrac, p, cdom_factor, avej, k1p, u_P, po4P, no3, no3P, npp, phyt, no3upt, o2
      real(rk) :: Paulmier_a, Paulmier_z, Paulmier_R0, pflag, nflag, jdiar, nupt, bct, dzt, morp
      real(rk) :: morpt, no3flag, o2_roc, alk_roc
      
      _LOOP_BEGIN_
         _GET_(self%id_sw_par, sw_par)
         _GET_(self%id_f1, f1)
         _GET_(self%id_dzt,dzt)
         _GET_(self%id_phosphorus, p)
         _GET_(self%id_phyt, phyt)
         _GET_(self%id_temp, temp)
         _GET_(self%id_o2, o2)
         _GET_SURFACE_(self%id_dayfrac, dayfrac)
         pflag = 0.5_rk + sign(0.5_rk,phyt - trcmin)
         nflag = 0.5_rk + sign(0.5_rk,p - trcmin)
         
         gl = 2._rk*self%alpha*sw_par
         
         ! iron limitation
         if (self%fe_limitation) then
             _GET_(self%id_temp, depth)
             _GET_HORIZONTAL_(self%id_latitude, lat)
             if ((self%so_fe_fer .and. lat < -40.0_rk) .or. (depth .le. 225._rk)) then
                 felimit = 1.0_rk
             else
                 _GET_SURFACE_(self%id_fe, fe_conc)
                 felimit = fe_conc/(self%kfe+fe_conc)
             endif
         else
             felimit = 1.0_rk
         endif
                  
         ! temperature- and irondependent max growth rate:
         if (self%no_temp_sens) then
             if (self%nitrogen .and. self%nfix) then
                 jmax = max(0.0_rk,self%abio*2.72_rk*felimit)*self%jdiar
             else
                 jmax = self%abio*3.84_rk*felimit
             endif
             nupt = self%nupt0*4.01_rk
         else
             _GET_(self%id_temp, temp)
             bct = self%bbio**(self%cbio*temp)
             if (self%nitrogen .and. self%nfix) then
                 jmax = max(0.0_rk,self%abio*(bct-2.6_rk)*felimit)*self%jdiar
             else
                 jmax = self%abio*bct*felimit
             endif
             nupt = self%nupt0*bct
         endif
         
         gd = max(1.0e-14_rk,jmax*dayfrac)
         u1 = max(gl/gd,1.e-6_rk)
         u2 = u1*exp(f1*dzt)
         ! for the following approximation ensure that u1 < 20
         phi1 = log(u1+sqrt(1._rk+u1**2._rk))-(sqrt(1._rk+u1**2._rk)-1._rk)/u1
         phi2 = log(u2+sqrt(1._rk+u2**2._rk))-(sqrt(1._rk+u2**2._rk)-1._rk)/u2
         
         if (self%cdom_attenuation) then
             cdom_factor = 1.2_rk
         else
             cdom_factor = 1._rk
         endif
         avej = gd*(phi1 - phi2)/(cdom_factor*abs(f1)*dzt)
         
         ! mortality
         morp = self%nup*phyt * pflag ! if (self%nitrogen) morp = 0.0_rk
         morpt = nupt*phyt * pflag

         k1p = self%k1n*redptn
         po4P = jmax*p/(k1p + p) 
         u_P = min(avej, po4P)
         no3flag = 1.0_rk
         if (self%nitrogen) then
             _GET_(self%id_no3, no3)
             no3flag = 0.5_rk + sign(0.5_rk,no3 - trcmin)
             if (self%nfix) then
                 morp = 0.0_rk
                 npp = u_P*phyt * no3flag*nflag
                 no3upt = no3/(self%k1n + no3)*npp * no3flag*nflag ! nitrate uptake
                 o2_roc =    (npp - no3upt)*1.25e-3_rk    ! check units of oxygen tracer
                 alk_roc = - (npp - no3upt)*1e-3_rk
                 
                 _ADD_SOURCE_(self%id_oxi,  o2_roc)
                 _ADD_SOURCE_(self%id_alk, alk_roc)
                 _ADD_SOURCE_(self%id_no3, morpt-no3upt)
             else
                 u_P = min(u_P, jmax*no3/(self%k1n + no3))
                 npp = u_P*phyt * no3flag*nflag
                 no3P = jmax*no3/(self%k1n + no3)
                 _ADD_SOURCE_(self%id_no3, morpt-npp)
             endif             
         endif
         
         ! update state variables
         _ADD_SOURCE_(self%id_phosphorus, (morpt-npp)*redptn)
         _ADD_SOURCE_(self%id_phyt,        npp-morpt-morp)
         _ADD_SOURCE_(self%id_det,         morp)
         _ADD_SOURCE_(self%id_calc,        morp*self%capr*redctn)
         
         ! diagnostic output
         _SET_DIAGNOSTIC_(self%id_morp_out, morp)
         _SET_DIAGNOSTIC_(self%id_npp_out, npp)
         _SET_DIAGNOSTIC_(self%id_morpt_out, morpt)
         if (self%extra_diags) then
             _SET_DIAGNOSTIC_(self%id_avej_out, avej)
             _SET_DIAGNOSTIC_(self%id_po4P_out, po4P)
             if (self%nitrogen) then
                 _SET_DIAGNOSTIC_(self%id_no3P_out, no3P)
             endif
         endif
      _LOOP_END_
   end subroutine do

end module