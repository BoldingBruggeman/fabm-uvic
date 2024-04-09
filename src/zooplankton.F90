#include "fabm_driver.h"

module uvic_zooplankton
   use fabm_types
   use fabm_particle

   use uvic_shared

   implicit none

   private

   type, extends(type_particle_model), public :: type_uvic_zooplankton
      type (type_state_variable_id)                            :: id_zoop, id_o2, id_phosphorus, id_no3, id_det
      type (type_dependency_id)                                :: id_temp
      type (type_model_id),          allocatable, dimension(:) :: id_prey
      type (type_state_variable_id), allocatable, dimension(:) :: id_pbio

      
      real(rk)                            :: gbio, bbio, cbio, nu, gamma1, ge
      real(rk), allocatable, dimension(:) :: zpref
      integer                             :: nprey
      logical                             :: no_temp_sens, graz_upper_temp_limit, nexcr 
   contains
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class (type_uvic_zooplankton), intent(inout), target :: self
      integer,                    intent(in)            :: configunit

      ! parameters
      call self%get_parameter(self%no_temp_sens,          'no_temp_sens',          '',                    'turn off temperature sensitivity',  default=.false.)
      if (.not. self%no_temp_sens) then
          call self%get_parameter(self%bbio,                  'bbio',                  '-',                   'upper limit grazing intercept',     default=1.066_rk)
          call self%get_parameter(self%cbio,                  'cbio',                  'degC-1',              'upper limit grazing exponent',      default=1.0_rk)
      endif
      call self%get_parameter(self%graz_upper_temp_limit, 'graz_upper_temp_limit', '',                    'turn on upper temperature limit',   default=.true.)
      call self%get_parameter(self%nexcr,                 'nexcr',                 '',                    'turn on nitrogen excretion',        default=.false.)
                                                                                                                                               
      call self%get_parameter(self%gbio,                  'gbio',                  'd-1',                 'Maximum grazing rate at 0 deg C',   default=0.18_rk, scale_factor=d_per_s)
      call self%get_parameter(self%nprey,                 'nprey',                 '-',                   'number of prey types',              default=1)
      call self%get_parameter(self%nu,                    'nu',                    '(mmol m-3)^(-2) d-1', 'quadratic mortality',               default=0.34_rk, scale_factor=d_per_s)
      call self%get_parameter(self%gamma1,                'gamma1',                '-',                   'assimilation efficiency',           default=0.7_rk)
      call self%get_parameter(self%ge,                    'ge',                    '-',                   'growth efficiency',                 default=0.5_rk)
      
      ! prey coupling
      allocate(self%zpref(self%nprey))
      allocate(self%id_prey(self%nprey))
      allocate(self%id_pbio(self%nprey))
      do iprey = 1, self%nprey
          write (index, '(i0)') iprey
          call self%get_parameter(self%zpref(iprey),      'zpref'//trim(index),    '-',      'preference for prey '//trim(index), default=0.0_rk)
          
          call self%register_model_dependency(self%id_prey(iprey), 'prey'//trim(index))
          call self%register_state_dependency(self%id_pbio(i_prey), 'pbio'//trim(index), 'mmol N m^-3', 'prey '//trim(index)//' biomass')
          call self%request_coupling_to_model(self%id_pbio(i_prey), self%id_prey(i_prey), standard_variables%total_nitrogen)
      enddo
      
      ! variable registrations
      call self%register_state_variable(self%id_zoop, 'zoop', 'mmol N m-3', 'zooplankton concentration', minimum=0.0_rk, )
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_zoop)
      call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_zoop, scale_factor=redptn)
      
      ! environmental dependencies
      call self%register_dependency(self%id_temp,      standard_variables%temperature)
      call self%register_state_dependency(self%id_o2,         'o2',       'mmol O m-3',  'oxygen concentration')
      call self%register_state_dependency(self%id_det,        'detritus', 'mmol N m-3',  'detritus concentration')
      call self%register_state_dependency(self%id_phosphorus, 'p',        'mmol P m-3',  'dissolved inorganic phosphorous concentration')
      if (self%nexcr) then
          call self%register_state_dependency(self%id_no3, 'n03',     'mmol N m-3',  'dissolved inorganic nitrogen concentration')
      endif
   end subroutine initialize
   
   subroutine do(self, _ARGUMENTS_DO_)
      class (type_uvic_zooplankton), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_
      
      real(rk)                           :: bctz, zoop, o2, temp, gmax, thetaZ, excr, prey_s, zoop_roc
      real(rk),   dimension(self%nPrey)  :: pbio, ing, gs, graz
      
      _LOOP_BEGIN_
         _GET_(self%id_zoop, zoop)
         _GET_(self%id_o2, o2)
         _GET_(self%id_temp, temp)
         do iprey=1,self%nprey
             _GET_(self%id_pbio(iprey),pbio(iprey))
         enddo
         
         if (no_temp_sens) then
             if (graz_upper_temp_limit) then
                 bctz = (0.5_rk*(tanh(o2*1000._rk - 8._rk)+1._rk))*3.59_rk
             else
                 bctz = (0.5_rk*(tanh(o2*1000._rk - 8._rk)+1._rk))*3.92_rk
             endif
         else
             bctz = (0.5_rk*(tanh(o2*1000._rk - 8._rk)+1._rk))*self%bbio**(self%cbio*temp)
         endif
         ! Make the max grazing rate a function of temperature
         gmax = self%gbio*bctz
         
         if (sum(self%zpref) - 1._rk > 1.e-15_rk) then
             call self%fatal_error('zooplankton', 'prey preferences do not sum up to one')         
         endif
         thetaZ = sum(self%zpref*pbio) + kzoo
         ing = self%zpref/thetaZ
         gs = 
         graz = gmax*ing*pbio*zoop
         
         mor = self%nu*zoop*zoop
         
         excr = (self%gamma1-self%geZ)*sum(graz)
         
         zoop_roc = self%geZ*sum(graz) - mor
         det_roc = (1._rk-self%gamma1)*sum(graz) + mor
         
         !update state variable sources/sinks
         _ADD_SOURCE_(self%id_zoop, zoop_roc)
         
         do iprey = 1, self%nprey
             do istate = 1, size(self%id_prey(iprey)%state
                 _GET_(self%id_prey(iprey)%state(istate), prey_s)
                 _ADD_SOURCE_(self%id_prey(iprey)%state(istate), - gs(iprey)*prey_s)
             enddo
         enddo
                 
         _ADD_SOURCE_(self%id_det, det_roc)
         
         _ADD_SOURCE_(self%id_phosphorus, redptn*excr)
         if (self%nexcr) then
             _ADD_SOURCE_(self%id_no3, excr)
         endif  
      _LOOP_END_
   end subroutine do

end module