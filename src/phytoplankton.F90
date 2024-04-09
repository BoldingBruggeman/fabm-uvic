#include "fabm_driver.h"

module uvic_phytoplankton
   use fabm_types
   use fabm_particle

   use uvic_shared

   implicit none

   private

   type, extends(type_particle_model), public :: type_uvic_phytoplankton
      type (type_state_variable_id)     :: id_phyt, id_phosphorus, id_no3, id_det
      type (type_surface_dependency_id) :: id_fe, id_dayfrac
      type (type_dependency_id)         :: id_temp, id_depth, id_sw_par, id_f1
      type (type_global_dependency_id)  :: id_latitude
      
      real(rk)                            :: abio, bbio, cbio, kfe, alpha, k1n, nup, nupt0
      logical                             :: no_temp_sens, nitrogen, no3_sens, fe_limitation, so_fe_fer, cdom_attenuation
   contains
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class (type_uvic_phytoplankton), intent(inout), target :: self
      integer,                    intent(in)            :: configunit

      ! parameters
      call self%get_parameter(self%alpha,        'alpha',        '(W/m2)-1 d-1','Initial slope P-I curve',                     default=0.1_rk, scale_factor=d_per_s)
      call self%get_parameter(self%k1n,          'k1n',          'mmol N m-3',  'Half saturation constant for N uptake',       default=0.7_rk)
      call self%get_parameter(self%nup,          'nup',          'd-1',         'Specific mortality rate',                     default=0.025_rk)
      call self%get_parameter(self%nupt0,        'nupt0',        'd-1',         'temp. dependent specific mortality rate',     default=0.02_rk)
      call self%get_parameter(self%nitrogen,     'nitrogen',     '',            'turn on nitrogen dependency',                 default=.true.)
      !if (self%nitrogen) then                                                   
      !    call self%get_parameter(self%no3_sens, 'n03_sens',     '',            'turn on remineralization nitrate dependence', default=.true.)
      !endif                                                                     
      call self%get_parameter(self%abio,         'abio',         'd-1',         'maximum growth rate',                         default=0.18_rk, scale_factor=d_per_s)
      call self%get_parameter(self%no_temp_sens, 'no_temp_sens', '',            'turn off temperature sensitivity',            default=.false.)                                                                                                                               
      if (.not. self%no_temp_sens) then                                         
          call self%get_parameter(self%bbio,     'bbio',         '-',           'temperature dependency intercept',            default=1.066_rk)
          call self%get_parameter(self%cbio,     'cbio',         'degC-1',      'temperature dependency exponent',             default=1.0_rk)
      endif                                                                     
      call self%get_parameter(self%fe_limitation,'fe_limitation','',            'turn on iron limitation',                     default=.true.)                                                                                                                               
      if (self%fe_limitation) then                                              
          call self%get_parameter(self%kfe,      'kfe',          'mmol Fe m-3', 'iron half saturation constant',               default=0.1_rk)
          call self%get_parameter(self%so_fe_fer,'so_fe_fer',    '',            'iron fertilization in the southern ocean',    default=.false.)
      endif
      call self%get_parameter(self%cdom_attenuation,'cdom_attenuation','',      'cdom_attenuation',                            default=.false.)
      
      
      ! variable registrations
      call self%register_state_variable(self%id_phyt, 'phytoplankton', 'mmol N m-3', 'phytoplankton concentration', minimum=0.0_rk, )
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_phyt)
      call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_phyt, scale_factor=redptn)
      call self%add_to_aggregate_variable(type_interior_standard_variable(name='total_phytoplankton',units='mmol N m-3'), self%id_phyt)
      
      ! environmental dependencies
      call self%register_dependency(self%id_sw_par, 'sw_par', 'W m-2', 'average layer PAR')
      call self%register_dependency(self%id_f1,     'f1',     '-',     'current layer light extinction coeff.')
      call self%register_dependency(self%id_dayfrac,'dayfrac','-',     'daylength fraction')
      call self%register_state_dependency(self%id_phosphorus, 'p',     'mmol p m-3',  'dissolved inorganic phosphorous concentration')
      call self%register_state_dependency(self%id_det,        'detritus', 'mmol N m-3',  'detritus concentration')
      if (self%nitrogen) then
          call self%register_state_dependency(self%id_no3, 'no3',     'mmol n m-3',  'dissolved inorganic nitrogen concentration')
      endif
      if (self%fe_limitation) then
          call self%register_surface_dependency(self%id_fe, 'fe',     'mmol Fe m-3',  'photic zone iron concentration')
      endif
      call self%register_dependency(self%id_temp,      standard_variables%temperature)
      call self%register_dependency(self%id_depth,     standard_variables%depth)
      call self%register_dependency(self%id_latitude,  standard_variables%latitude)
      
   end subroutine initialize
   
   subroutine do(self, _ARGUMENTS_DO_SURFACE_)
      class (type_uvic_phytoplankton), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_
      
      real(rk) :: depth, lat, temp, felimit, fe_conc, jmax, sw_par, gl, f1, gd, u1, u2, phi1, phi2, dayfrac, p, cdom_factor, avej, k1p, u_P, po4P, no3, no3P, npp, phyt, temp
      
      _LOOP_BEGIN_
         _GET_(self%id_sw_par, sw_par)
         _GET_(self%id_f1, f1)
         _GET_(self%id_dayfrac, dayfrac)
         _GET_(self%id_phosphorus, p)
         _GET_(self%id_phyt, phyt)
         _GET_(self%id_temp, temp)
         
         gl = 2._rk*self%alpha*sw_par
         
         ! iron limitation
         if (self%fe_limitation) then
             _GET_(self%id_temp, depth)
             _GET_(self%id_latitude, lat)
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
             jmax = self%abio*3.84_rk*felimit
             nupt = self%nupt0*4.01_rk
         else
             _GET_(self%id_temp, temp)
             bct = self%bbio**(self%cbio*temp)
             jmax = self%abio*bct*felimit
             nupt = self%nupt0*bct
         endif
         
         gd = jmax*dayfrac
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
         morp = self%nup*phyt
         morpt = nupt*phyt

         
         k1p = self%k1n*redptn
         po4P = jmax*p/(k1p + p) 
         u_P = min(avej, po4P)
         if (self%nitrogen) then
             _GET_(self%id_no3, no3)
             u_P = min(u_P, jmax*no3/(self%k1n + no3))
             no3P = jmax*no3/(self%k1n + no3)
             _ADD_SOURCE_(self%id_no3, morpt-u_P*phyt)
         endif
         npp = u_P*phyt
         
         _ADD_SOURCE_(self%id_phosphorus, (morpt-npp)*redptn)
         _ADD_SOURCE_(self%id_phyt,        npp-morpt-morp)
         _ADD_SOURCE_(self%id_det,         morp)
      _LOOP_END_
   end subroutine do

end module