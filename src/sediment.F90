#include "fabm_driver.h"

module uvic_sediment
   use fabm_types
   use fabm_particle
   use fabm_expressions

   use uvic_shared

   implicit none

   private

   type, extends(type_particle_model), public :: type_uvic_sediment
      type (type_bottom_diagnostic_variable_id), allocatable, dimension(:) :: id_co2, id_hco3, id_co3, id_o2s, id_orgml, id_calgg, id_orggg, id_bmass, id_bcalfrac
      type (type_bottom_diagnostic_variable_id), allocatable, dimension(:) :: id_cal_c
      type (type_bottom_dependency_id), allocatable, dimension(:)          :: id_co2_in, id_hco3_in, id_co3_in, id_o2s_in, id_orgml_in, id_calgg_in, id_orggg_in
      type (type_bottom_dependency_id), allocatable, dimension(:)          :: id_bmass_in, id_bcalfrac_in, id_cal_c_in

      type (type_bottom_diagnostic_variable_id)                            :: id_zrct, id_weathflx, id_sed_ml_mass
      type (type_state_variable_id)                                        :: id_dic, id_alk, id_o2
      type (type_bottom_dependency_id)                                     :: id_zb, id_zrct_in
      type (type_bottom_dependency_id)                                     :: id_sed_ml_mass_in
      type (type_dependency_id)                                            :: id_sal, id_temp, id_dic_mean, id_alk_mean, id_sal_mean, id_temp_mean, id_o2_mean
      type (type_dependency_id)                                            :: id_dzt
      
      type (type_dependency_id)                                            :: id_rain_org_mean, id_rain_cal_mean
      type (type_dependency_id)                                            :: id_rain_org, id_rain_cal 
      
      type (type_bottom_diagnostic_variable_id)                            :: id_counter
      type (type_bottom_dependency_id)                                            :: id_counter_in

      
      real(rk)                            :: dtsed, dissc, dissn
      integer                             :: ibmax, kmax, nzmax
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

contains

   subroutine initialize(self, configunit)
      class (type_uvic_sediment), intent(inout), target :: self
      integer,                    intent(in)            :: configunit

      integer           :: i
      character(len=64) :: index
      
      call self%register_diagnostic_variable(self%id_counter, 'counter', '-', 'counter', output=output_instantaneous, missing_value=0.0_rk)
      call self%register_dependency(self%id_counter_in,  'counter', '-', 'counter')
      
      ! parameters
      call self%get_parameter(self%dtsed, 'dtsed', 'd', 'sediment time step',                  default=30.0_rk)
      call self%get_parameter(self%ibmax, 'ibmax', '-', 'max number of buried layers',         default=20)
      call self%get_parameter(self%kmax,  'kmax',  '-', 'number of mixed layers',              default=7)
      call self%get_parameter(self%nzmax, 'nzmax', '-', 'max number of mixed layers',          default=8)
      call self%get_parameter(self%dissc, 'dissc', '',  'calcite dissolution constant 1',      default=1.1574e-5_rk)
      call self%get_parameter(self%dissn, 'dissn', '',  'calcite dissolution constant 2',      default=4.5_rk)
      
      ! variable registrations
      call self%register_diagnostic_variable(self%id_zrct, 'zrct', 'cm', 'max respiration depth', output=output_instantaneous)
      call self%register_dependency(self%id_zrct_in,  'zrct', 'cm', 'max respiration depth')
      
      call self%register_diagnostic_variable(self%id_sed_ml_mass, 'sed_ml_mass', 'g cm-2', 'mixed layer mass', output=output_instantaneous)
      call self%register_dependency(self%id_sed_ml_mass_in,  'sed_ml_mass', 'g cm-2', 'mixed layer mass')

      allocate(self%id_co2  (self%nzmax))
      allocate(self%id_hco3 (self%nzmax))
      allocate(self%id_co3  (self%nzmax))
      allocate(self%id_o2s  (self%nzmax))
      allocate(self%id_orgml(self%nzmax))
      allocate(self%id_calgg(self%nzmax))
      allocate(self%id_orggg(self%nzmax))
      allocate(self%id_cal_c(self%nzmax))
      allocate(self%id_co2_in  (self%nzmax))
      allocate(self%id_hco3_in (self%nzmax))
      allocate(self%id_co3_in  (self%nzmax))
      allocate(self%id_o2s_in  (self%nzmax))
      allocate(self%id_orgml_in(self%nzmax))
      allocate(self%id_calgg_in(self%nzmax))
      allocate(self%id_orggg_in(self%nzmax))
      allocate(self%id_cal_c_in(self%nzmax))
      
      do i=1,self%nzmax
          write(index,'(i0)') i
          call self%register_diagnostic_variable(self%id_co2(i),     'co2_'//trim(index),     'mol l-1','co2 conc in layer '//trim(index),                output=output_instantaneous, missing_value=-1.e-300_rk)
          call self%register_diagnostic_variable(self%id_hco3(i),    'hco3_'//trim(index),    'mol l-1','hco3 conc in layer '//trim(index),               output=output_instantaneous)
          call self%register_diagnostic_variable(self%id_co3(i),     'co3_'//trim(index),     'mol l-1','co3 conc in layer '//trim(index),                output=output_instantaneous)
          call self%register_diagnostic_variable(self%id_o2s(i),     'o2s_'//trim(index),     'mol l-1','o2 conc in layer '//trim(index),                 output=output_instantaneous)
          call self%register_diagnostic_variable(self%id_orgml(i),   'orgml_'//trim(index),   'g cm-2', 'organic carbon mass in layer '//trim(index),     output=output_instantaneous)
          call self%register_diagnostic_variable(self%id_calgg(i),   'calgg_'//trim(index),   '-',      'calcite fraction in layer '//trim(index),        output=output_instantaneous)
          call self%register_diagnostic_variable(self%id_orggg(i),   'orggg_'//trim(index),   '-',      'organic carbon fraction in layer '//trim(index), output=output_instantaneous)
          call self%register_diagnostic_variable(self%id_cal_c(i),   'cal_c_'//trim(index),   '-',      'cal_c'//trim(index),                             output=output_instantaneous)

          call self%register_dependency(self%id_co2_in(i),  'co2_'//trim(index),  'mol l-1', 'co2 conc in layer '//trim(index))
          call self%register_dependency(self%id_hco3_in(i), 'hco3_'//trim(index), 'mol l-1', 'hco3 conc in layer '//trim(index))
          call self%register_dependency(self%id_co3_in(i),  'co3_'//trim(index),  'mol l-1', 'co3 conc in layer '//trim(index))
          call self%register_dependency(self%id_o2s_in(i),  'o2s_'//trim(index),  'mol l-1', 'o2 conc in layer '//trim(index))
          call self%register_dependency(self%id_orgml_in(i),'orgml_'//trim(index),'g cm-2',  'organic carbon mass in layer '//trim(index))
          call self%register_dependency(self%id_calgg_in(i),'calgg_'//trim(index),'-',       'calcite fraction in layer '//trim(index))
          call self%register_dependency(self%id_orggg_in(i),'orggg_'//trim(index),'-',       'organic carbon fraction in layer '//trim(index))
          call self%register_dependency(self%id_cal_c_in(i),'cal_c_'//trim(index),'-',       'cal_c '//trim(index))
      end do
      allocate(self%id_bmass(self%ibmax))
      allocate(self%id_bcalfrac(self%ibmax))
      allocate(self%id_bmass_in(self%ibmax))
      allocate(self%id_bcalfrac_in(self%ibmax))
      do i=1,self%ibmax
          write(index,'(i0)') i
          call self%register_diagnostic_variable(self%id_bmass(i),   'bmass_'//trim(index),   'g cm-2', 'buried CaCO3 mass in layer '//trim(index))
          call self%register_diagnostic_variable(self%id_bcalfrac(i),'bcalfrac_'//trim(index),'-',      'buried calcite fraction in layer '//trim(index))

          call self%register_dependency(self%id_bmass_in(i),   'bmass_'//trim(index),   'g cm-2', 'buried CaCO3 mass in layer '//trim(index))
          call self%register_dependency(self%id_bcalfrac_in(i),'bcalfrac_'//trim(index),'-',      'buried calcite fraction in layer '//trim(index))
      end do
      
      ! environmental dependencies
      call self%register_dependency(self%id_temp,standard_variables%temperature)
      call self%register_dependency(self%id_sal, standard_variables%practical_salinity)
      call self%register_dependency(self%id_zb,  standard_variables%bottom_depth)
      call self%register_dependency(self%id_dzt, standard_variables%cell_thickness)
      
      ! register dependencies
      call self%register_dependency(self%id_rain_org, 'rain_org', 'umol cm-2 s-1', 'detritus sedimentation flux')
      call self%register_dependency(self%id_rain_cal, 'rain_cal', 'umol cm-2 s-1', 'calcite sedimentation flux')
      
      call self%register_state_dependency(self%id_dic,  'DIC',        'umol C cm-3',   'dissolved inorganic carbon concentration')
      call self%register_state_dependency(self%id_alk,  'alkalinity', 'umol cm-3',     'ocean alkalinity')
      call self%register_state_dependency(self%id_o2,   'o2',         'umol O cm-3',   'oxygen concentration')
      
      call self%register_dependency(self%id_temp_mean,     temporal_mean(self%id_temp,     period=self%dtsed*86400._rk, resolution=1800._rk, missing_value=-1.e-300_rk))
      call self%register_dependency(self%id_sal_mean,      temporal_mean(self%id_sal,      period=self%dtsed*86400._rk, resolution=1800._rk, missing_value=-1.e-300_rk))
      call self%register_dependency(self%id_rain_org_mean, temporal_mean(self%id_rain_org, period=self%dtsed*86400._rk, resolution=1800._rk, missing_value=-1.e-300_rk))
      call self%register_dependency(self%id_rain_cal_mean, temporal_mean(self%id_rain_cal, period=self%dtsed*86400._rk, resolution=1800._rk, missing_value=-1.e-300_rk))
      call self%register_dependency(self%id_dic_mean,      temporal_mean(self%id_dic,      period=self%dtsed*86400._rk, resolution=1800._rk, missing_value=-1.e-300_rk))
      call self%register_dependency(self%id_alk_mean,      temporal_mean(self%id_alk,      period=self%dtsed*86400._rk, resolution=1800._rk, missing_value=-1.e-300_rk))
      call self%register_dependency(self%id_o2_mean,       temporal_mean(self%id_o2,       period=self%dtsed*86400._rk, resolution=1800._rk, missing_value=-1.e-300_rk))   
      
      
      !register output diagnostics
      call self%register_diagnostic_variable(self%id_weathflx,   'weathflux',   'umol cm-2 s-1', 'weathering flux')

   end subroutine initialize
    
   subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
      class (type_uvic_sediment), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_
      
!      integer  :: kmax=7, nzmax=8
      real(rk) :: rain_org, rain_org_p, temp, sal, alk, dic, o2_bw, rain_cal, rain_cal_p, rain_clay_p, water_z, k1, k2, k3, csat, co2, hco3
      real(rk) :: co3, carb(self%nzmax,3), delz(self%nzmax), zsed(self%nzmax), o2(self%nzmax), rc, orgml(self%nzmax), calgg(self%nzmax)
      real(rk) :: orggg(self%nzmax), buried_mass(self%ibmax), buried_calfrac(self%ibmax), dcpls(self%nzmax,3), dcmin(self%nzmax,3)
      real(rk) :: pore(self%nzmax), form(self%nzmax), ttrcal, zrct, dopls(self%nzmax), domin(self%nzmax), dbpls(self%nzmax), weathflx
      real(rk) :: dbmin(self%nzmax), tmp, dzt, resp_c(self%nzmax,3), ttrorg, cal_c(self%nzmax), ttrtc, ttral, diftc, difal=0._rk, sed_ml_mass
      integer  :: k, i, counter
      
      _BOTTOM_LOOP_BEGIN_
         _GET_BOTTOM_(self%id_counter_in,counter)
         counter = counter + 1._rk
         _SET_BOTTOM_DIAGNOSTIC_(self%id_counter, counter)
      
         ! get environmental variables
         _GET_(self%id_temp_mean,temp)                       ! I NEED AVERAGE CONCENTRATIONS SINCE LAST SEDIMENT CALL
         _GET_(self%id_sal_mean,sal)                         ! I NEED AVERAGE CONCENTRATIONS SINCE LAST SEDIMENT CALL
         _GET_(self%id_alk_mean,alk)                         ! I NEED AVERAGE CONCENTRATIONS SINCE LAST SEDIMENT CALL
         _GET_(self%id_dic_mean,dic)                         ! I NEED AVERAGE CONCENTRATIONS SINCE LAST SEDIMENT CALL
         _GET_(self%id_o2_mean,o2_bw)                        ! I NEED AVERAGE CONCENTRATIONS SINCE LAST SEDIMENT CALL
         _GET_(self%id_rain_org_mean,rain_org)     ! I NEED AVERAGE CONCENTRATIONS SINCE LAST SEDIMENT CALL
         _GET_(self%id_rain_cal_mean,rain_cal)   ! I NEED AVERAGE CONCENTRATIONS SINCE LAST SEDIMENT CALL
         if (temp == -1.e-300_rk) then
             _GET_(self%id_temp,temp)
             _GET_(self%id_sal,sal)           
             _GET_(self%id_alk,alk)           
             _GET_(self%id_dic,dic)           
             _GET_(self%id_o2,o2_bw)          
             _GET_(self%id_rain_org,rain_org) 
             _GET_(self%id_rain_cal,rain_cal) 
         endif
         
         !temp = max(1e-32_rk, temp)
         !sal = max(1e-32_rk, sal)
         !alk = max(1e-32_rk, alk)
         !dic = max(1e-32_rk, dic)
         !o2_bw = max(1e-32_rk, o2_bw)
         !rain_org = max(1e-32_rk, rain_org)
         !rain_cal = max(1e-32_rk, rain_cal)
         _GET_BOTTOM_(self%id_zb,water_z)
         _GET_(self%id_dzt,dzt)
         
         ! get saved internal variables
         _GET_BOTTOM_(self%id_zrct_in,zrct)
         _GET_BOTTOM_(self%id_sed_ml_mass_in,sed_ml_mass)
         do i=1,self%nzmax              
            _GET_BOTTOM_(self%id_co2_in(i),  carb(i,1))
            _GET_BOTTOM_(self%id_hco3_in(i), carb(i,2))
            _GET_BOTTOM_(self%id_co3_in(i),  carb(i,3))
            _GET_BOTTOM_(self%id_o2s_in(i),  o2(i))
            _GET_BOTTOM_(self%id_orgml_in(i),orgml(i))
            _GET_BOTTOM_(self%id_calgg_in(i),calgg(i))
            _GET_BOTTOM_(self%id_orggg_in(i),orggg(i))
            _GET_BOTTOM_(self%id_cal_c_in(i),cal_c(i))
         end do
         do i=1,self%ibmax              
            _GET_BOTTOM_(self%id_bmass_in(i),   buried_mass(i))
            _GET_BOTTOM_(self%id_bcalfrac_in(i),buried_calfrac(i))
         end do
         
         dcpls          = 0._rk
         dcmin          = 0._rk
         pore           = 0._rk
         form           = 0._rk
         dopls          = 0._rk
         domin          = 0._rk
         dbpls          = 0._rk
         dbmin          = 0._rk
         delz           = 3._rk
         delz(1:5)      = [0._rk, 0.5_rk, 0.5_rk, 1._rk, 2._rk]
         ttrcal=0._rk
         zsed(1)        = 0._rk
         do k=2,self%kmax
             zsed(k) = zsed(k-1) + delz(k)
         enddo
         
         if (carb(1,1) == -1.e-300_rk) then
             ! initialize internal variables
             carb(:,1)      = 20.e-6_rk
             carb(:,2)      = 2000.e-6_rk
             carb(:,3)      = 80.e-6_rk
             o2             = 150.e-6_rk
             calgg          = 0.1_rk
             orggg          = 0.0002_rk
             orgml          = 0._rk
             sed_ml_mass    = 0._rk
             do k=1,self%kmax
                 sed_ml_mass = sed_ml_mass + delz(k)*(1._rk-pore(k))*2.5_rk
             enddo
             zrct = zsed(self%kmax)
             buried_mass    = 0._rk
             buried_mass(1) = 500._rk
             buried_calfrac = 0.1_rk
             cal_c = 0.0_rk
         endif         
         
         ! convert units
         !sal = sal * 1000._rk + 35._rk ! No conversion needed - FABM standard is already in PSU
         alk = alk * 0.001_rk                        !convert concentrations from mol m-3 (or umol cm-3) to mol l-1
         dic = dic * 0.001_rk                        !convert concentrations from mol m-3 (or umol cm-3) to mol l-1
         o2_bw = max(o2_bw * 0.001_rk, 1.e-6_rk)     !convert concentrations from mol m-3 (or umol cm-3) to mol l-1
         rain_org_p = rain_org * 1.e-6_rk * self%dtsed*86400._rk ! convert from umol/cm2/s to mol/cm2/dtsed
         rain_cal_p = rain_cal * 1.e-6_rk * self%dtsed*86400._rk ! convert from umol/cm2/s to mol/cm2/dtsed
         rain_clay_p = 9.1e-6_rk * (1.e6_rk * rain_org_p)**1.41_rk
         
         call calc_k (temp, sal, water_z, k1, k2, k3, csat)
         
         call calc_buff (alk, dic, sal, k1, k2, k3, co2, hco3, co3)
         
         o2(1)     = o2_bw
         carb(1,1) = co2
         carb(1,2) = hco3
         carb(1,3) = co3
         rc = 2.e-9_rk
                  
         !-----------------------------------------------------------------------
         !         step ahead time-dependent calcite concentration
         !-----------------------------------------------------------------------
          call sed_const_cal (zsed, delz, form, pore, self%kmax, o2, zrct, carb, orgml, orggg, calgg, &
              rain_cal_p, rain_org_p, rain_clay_p, rc, self%dissc, self%dissn, csat, &
              k1, k2, dopls, domin, dcpls, dcmin, dbpls, dbmin, resp_c, cal_c, ttrorg, ttrcal, ttrtc,  &
              ttral, diftc, difal, self%nzmax, counter)
          
          call bury (zsed, delz, pore, self%kmax, calgg, orggg, sed_ml_mass, rain_cal_p, ttrcal, &
              rain_clay_p, buried_mass, buried_calfrac, self%nzmax, self%ibmax, counter)
          
         ! save internal variables
         _SET_BOTTOM_DIAGNOSTIC_(self%id_zrct,zrct)
         do i=1,self%nzmax              
            _SET_BOTTOM_DIAGNOSTIC_(self%id_co2(i),  carb(i,1))
            _SET_BOTTOM_DIAGNOSTIC_(self%id_hco3(i), carb(i,2))
            _SET_BOTTOM_DIAGNOSTIC_(self%id_co3(i),  carb(i,3))
            _SET_BOTTOM_DIAGNOSTIC_(self%id_o2s(i),  o2(i))
            _SET_BOTTOM_DIAGNOSTIC_(self%id_orgml(i),orgml(i))
            _SET_BOTTOM_DIAGNOSTIC_(self%id_calgg(i),calgg(i))
            _SET_BOTTOM_DIAGNOSTIC_(self%id_orggg(i),orggg(i))
            _SET_BOTTOM_DIAGNOSTIC_(self%id_cal_c(i),cal_c(i))
         end do
         do i=1,self%ibmax              
            _SET_BOTTOM_DIAGNOSTIC_(self%id_bmass(i),   buried_mass(i))
            _SET_BOTTOM_DIAGNOSTIC_(self%id_bcalfrac(i),buried_calfrac(i))
         end do

         ! update boundary fluxes
         tmp = (ttrcal - rain_cal_p)/(self%dtsed*86400._rk)
         weathflx = -tmp
         tmp = tmp/(dzt*100._rk) ! convert layer height to cm and convert tmp from umol cm-2 s-1 to umol cm-3 s-1
         
         _ADD_BOTTOM_FLUX_(self%id_dic,  tmp) ! flux from bottom to interior
         _ADD_BOTTOM_FLUX_(self%id_alk,  tmp*2._rk) ! flux from bottom to interior

         _SET_BOTTOM_DIAGNOSTIC_(self%id_sed_ml_mass, sed_ml_mass)
         _SET_BOTTOM_DIAGNOSTIC_(self%id_weathflx, weathflx)
      _BOTTOM_LOOP_END_
   end subroutine do_bottom
   
   subroutine calc_k (temp, sal, z, k1, k2, kb, csat)

      implicit none

      real(rk) kprime, delv, rr, dk

      parameter (kprime = 4.75e-7_rk) ! mol2 / kg2
      parameter (delv = -44._rk)      ! cm3 / mol
      parameter (rr = 83.14_rk)       ! cm3 bar / k mol
      parameter (dk = -.0133_rk)      ! cm3 / bar mol

!     arguments
      real(rk) temp, sal, z, k1, k2
      real(rk) kb, csat

!     local variables
      real(rk) tk, cp, prat, pres, kpres

      tk = temp + 273.15_rk

!       k1 and k2 (apparent), from mehrbach

      k1 = 13.7201_rk - 0.031334_rk*tk - 3235.76_rk/tk - 1.3e-5_rk*sal*tk + 0.1032_rk*sal**(0.5_rk)

      k1 = 10._rk**(k1)
      cp = (z/10._rk)/83.143_rk/tk
      prat = (24.2_rk - 0.085_rk*temp)*cp
      prat = exp(prat)
      k1 = k1*prat

      k2 = - 5371.9645_rk - 1.671221_rk*tk + 128375.28_rk/tk &
          + 2194.3055_rk*log(tk)/2.30259_rk - 0.22913_rk*sal &
          - 18.3802_rk*log(sal)/2.30259_rk + 8.0944e-4_rk*sal*tk &
          + 5617.11_rk*log(sal)/tk/2.30259_rk - 2.136_rk*sal/tk

      k2 = 10._rk**(k2)
      prat = (16.4_rk - 0.04_rk*temp)*cp
      prat = exp(prat)
      k2 = k2*prat

!    lyman's kb
     kb = 2291.9_rk/(temp + 273._rk) + 0.01756_rk*(temp + 273._rk) &
         - 3.385_rk - .32051_rk*(sal/1.80655_rk)**(1._rk/3._rk)
     kb = 10._rk**(-kb)
     prat = (27.5_rk - 0.095_rk*temp)*cp
     prat = exp(prat)
     kb = kb*prat

!    calculate sayles calcite saturation state
     pres = z/10._rk ! bar
     kpres = LOG(4.75e-7_rk) - delv/(rr*(temp+273._rk))*(pres) &
         + 0.5_rk*dk/(rr*(temp + 273._rk))*(pres)**2_rk
     kpres = eXP(kpres)
     csat = kpres/0.01_rk

     return
   end subroutine calc_k
   
   subroutine calc_buff (alk, tco2, sal, k1, k2, kb, co2, hco3, co3)

      implicit none

      integer icnt

      real(rk) alk, tco2, sal, k1, k2
      real(rk) kb, co2, hco3, co3

      real(rk) tbor, tkt, tk, c1, c2
      real(rk) c4, a, x, ah1, aht

!     all units of moles/l

      tbor = 4.106e-4_rk*sal/35._rk
      c1 = k1/2.0_rk
      c2 = 1.0_rk - 4.0_rk*k2/k1
      c4 = tbor*kb
      aht = 0.74e-8_rk

      do icnt=1,100
          a = alk - c4/(kb + aht)
          x = a/tco2
          ah1 = c1/x*(1._rk - x + sqrt(1._rk + c2*x*(-2._rk + x)))
          aht=ah1
      enddo

      co3 = (a - tco2)/(1.0_rk - (ah1*ah1)/(k1*k2))
      hco3 = tco2/(1._rk + ah1/k1 + k2/ah1)
      co2 = tco2/(1._rk + k1/ah1 + k1*k2/(ah1*ah1))
      return
   end subroutine calc_buff
   
   subroutine sed_const_cal (zsed, delz, form, pore, kmax, o2, zrct, carb, orgml, orggg, calgg, &
              raincal, rainorg, rainclay, rc, dissc, dissn, csat, &
              u1, u2, dopls, domin, dcpls, dcmin, dbpls, dbmin, resp_c, cal_c, ttrorg, ttrcal, ttrtc,  &
              ttral, diftc, difal, nzmax, counter)
!-----------------------------------------------------------------------
!  finds the steady state pore water chemistry with const calcite
!-----------------------------------------------------------------------

      implicit none

      integer inner_loop_limit, iter, counter
      integer k, kmax, nzmax

      !   arguments
      real(rk) zsed(nzmax),delz(nzmax), form(nzmax), pore(nzmax)
      real(rk) o2(nzmax), zrct, carb(nzmax,3)
      real(rk) orgml(nzmax), orggg(nzmax)
      real(rk) calgg(nzmax), raincal, rainorg
      real(rk) rainclay, rc
      real(rk) dissc, dissn, csat, u1, u2
      real(rk) dopls(nzmax), domin(nzmax), dcpls(nzmax,3)
      real(rk) dcmin(nzmax,3), dbpls(nzmax), dbmin(nzmax)
      real(rk) resp_c(nzmax,3), cal_c(nzmax)
      real(rk) raincal_cutoff

      ! results
      real(rk) ttrorg, ttrcal, ttrtc, ttral
      real(rk) difal, diftc

      ! internal variables
      logical wrong
      real(rk) expb, db, difo2, difc(3), pore_max, exp_pore

      expb = 3.0_rk
      difo2 = 12.1e-6_rk
      difc(1) = 10.5e-6_rk
      difc(2) = 6.4e-6_rk
      difc(3) = 5.2e-6_rk
      db = 0.15_rk

      pore_max = 1._rk - (0.483_rk + 0.45_rk*calgg(kmax))/2.5_rk
      exp_pore = 0.25_rk*calgg(kmax) + 3._rk*(1._rk - calgg(kmax))
      do k=2,kmax
        pore(k) = exp(-zsed(k)/exp_pore)*(1._rk-pore_max) + pore_max
      enddo
      
      do k=1,kmax
          form(k) = pore(k)**expb
      enddo
      
      if (counter==1) then
          do k=2,kmax
              orgml(k) = orggg(k)*2.5_rk*(1._rk-pore(k))*1000._rk/12._rk
          enddo
      endif
      
      call calc_do2 (difo2, form, pore, delz, kmax, dopls, domin, nzmax)
     
      call calc_dc (difc, form, pore, delz, kmax, dcpls, dcmin, nzmax)
      
      call calc_db (db, pore, zsed, delz, kmax, dbpls, dbmin, nzmax)
      
      call o2org (rainorg, rc, kmax, zsed, delz, form, pore, dopls, domin, &
          dbpls, dbmin, o2, zrct, orgml, orggg, resp_c, nzmax)

      raincal_cutoff = 0.1e-6_rk

      wrong = .false.
      if (raincal .gt. raincal_cutoff) then
        wrong = .true.
      else
        ttrcal = raincal
        do k=2,kmax
          calgg(k) = 0.
        enddo
      endif

      inner_loop_limit = 200
      
      call co3ss (resp_c, cal_c, dissc, dissn, csat, u1, u2, zsed, delz, &
          form, pore, kmax, dcpls, dcmin, calgg, carb, ttrorg, &
          ttrcal, ttral, ttrtc, difal, diftc, nzmax, inner_loop_limit, &
          wrong, iter)

      if ((ttrcal.gt.raincal) .and. (calgg(kmax).lt.0.001_rk)) ttrcal = raincal

      return
   end subroutine sed_const_cal
              
   subroutine calc_do2 (difo2, form, pore, delz, kmax, dopls, domin, nzmax)
      implicit none

      integer i, kmax, nzmax

      real(rk) difo2, form(nzmax), pore(nzmax)
      real(rk) dopls(nzmax), domin(nzmax), delz(kmax)

      difo2 = 12.e-6_rk

      do i=3,kmax-1
        dopls(i) = difo2*((form(i+1) + form(i))*0.5_rk)*1._rk/pore(i)*2._rk/((delz(i+1) + delz(i))*delz(i))
        domin(i) = difo2*((form(i-1) + form(i))*0.5_rk)*1._rk/pore(i)*2._rk/((delz(i-1) + delz(i))*delz(i))
      enddo

      i=kmax
      dopls(i) = 0.
      domin(i) = difo2*((form(i-1) + form(i))*0.5_rk)*1._rk/pore(i)*2._rk/((delz(i-1) + delz(i))*delz(i))
      i=2
      dopls(i) = difo2*((form(i+1) + form(i))*0.5_rk)*1._rk/pore(i)*2._rk/((delz(i+1) + delz(i))*delz(i))
      domin(i) = difo2*(form(i)+1)*0.5_rk*1._rk/pore(i)*1._rk/delz(i)**2

      return
   end subroutine calc_do2
    
   subroutine calc_dc (difc, form, pore, delz, kmax, dcpls, dcmin, nzmax)
      implicit none

      integer i, j, kmax, nzmax

      real(rk) difc(3), form(nzmax), pore(nzmax), delz(kmax)
      real(rk) dcpls(nzmax,3), dcmin(nzmax,3)

      difc(1) = 10.5e-6_rk
      difc(2) = 6.4e-6_rk
      difc(3) = 5.2e-6_rk

      do i=3,kmax-1
        do j=1,3
          dcpls(i,j) = difc(j)*(delz(i)*form(i+1) + delz(i+1)*form(i))/(delz(i)+delz(i+1))*1._rk/pore(i)*(2._rk/((delz(i+1) + delz(i))*delz(i)))
          dcmin(i,j) = difc(j)*(delz(i)*form(i-1) + delz(i-1)*form(i))/(delz(i)+delz(i-1))*1._rk/pore(i)*(2._rk/((delz(i-1) + delz(i))*delz(i)))
        enddo
      enddo

      do j=1,3
        i = kmax
        dcpls(i,j) = 0.
        dcmin(i,j) = difc(j)*(delz(i)*form(i-1) + delz(i-1)*form(i))/(delz(i) + delz(i-1))*1._rk/pore(i)*(2._rk/((delz(i-1) + delz(i))*delz(i)))
        i = 2
        dcpls(i,j) = difc(j)*(delz(i)*form(i+1) + delz(i+1)*form(i))/(delz(i) + delz(i+1))*1._rk/pore(i)*(2._rk/((delz(i+1) + delz(i))*delz(i)))
        dcmin(i,j) = difc(j)*(form(i)+1)*0.5_rk*1._rk/pore(i)*(1._rk/(delz(i)**2._rk))
      enddo

      return
   end subroutine calc_dc
   
   subroutine calc_db (db, pore, zsed, delz, kmax, dbpls, dbmin, nzmax)
      implicit none

      integer k, kmax, nzmax

      !     external input variables
      real(rk) pore(nzmax), zsed(kmax+1), delz(kmax)

      !     results; external variables
      real(rk) dbpls(nzmax), dbmin(nzmax)

      real(rk) db

      db = 0.15_rk

      zsed(kmax+1) = zsed(kmax) + 1._rk

      do k=3,kmax-1
        dbpls(k) = db*2._rk/((delz(k) + delz(k+1))*delz(k))*(1._rk-pore(k)+1._rk-pore(k+1))/(1._rk-pore(k))
        dbmin(k) = db*2._rk/((delz(k) + delz(k-1))*delz(k))*(1._rk-pore(k)+1._rk-pore(k-1))/(1._rk-pore(k))
      enddo

      k = 2
      dbpls(k) = db*2._rk/((delz(k) + delz(k+1))*delz(k))*(2._rk - pore(k) - pore(k+1))/(1._rk - pore(k))
      dbmin(k) = 0._rk
      k = kmax
      dbpls(k) = 0._rk
      dbmin(k) = db*2._rk/((delz(k) + delz(k-1))*delz(k))*(2._rk - pore(k) - pore(k-1))/(1._rk - pore(k))

      return
   end subroutine calc_db
   
   subroutine o2org (rainorg, rc, kmax, zsed, delz, form, pore, dopls, domin, dbpls, dbmin, o2, &
       zrct, orgml, orggg, resp_c, nzmax)
      implicit none

      integer i, kmax, loop_limit, m, nzmax

!     arguments
      real(rk) rainorg, rc, zsed(kmax), delz(kmax)
      real(rk) form(nzmax), pore(nzmax), dopls(nzmax)
      real(rk) domin(nzmax), dbpls(nzmax), dbmin(nzmax)
      real(rk) o2(nzmax),zrct, orgml(nzmax)
      real(rk) orggg(nzmax)

!     results
      real(rk) resp_c(nzmax,3)

!     local variables -- for diagnostics
      real(rk) rct_c, rms_c, rct_o2, rms_o2

      loop_limit = 3

      do m=1,loop_limit
        call orgc (rainorg, rc, zsed, dbpls, dbmin, pore, zrct, delz, &
            orggg, orgml, rct_c, rms_c, kmax, nzmax)

        call o2ss (zsed, delz, orgml, pore, dopls, domin, rc, kmax, &
            zrct, o2, rct_o2, rms_o2, nzmax)

        zrct = min(zsed(kmax), zrct*o2(1)/(o2(1) - o2(kmax) + 1.e-20_rk))

        if (zrct .lt. 0.1_rk) zrct = 0.1_rk
      enddo

      do i=1,kmax
        if (zsed(i) .le. zrct) then
          resp_c(i,1) = rc*orgml(i)
        elseif (zsed(i-1) .le. zrct) then
          resp_c(i,1) = rc*orgml(i)*(zrct - zsed(i-1))/delz(i)
        elseif (zsed(i-1) .gt. zrct) then
          resp_c(i,1) = 0.
        endif
        resp_c(i,2) = 0.
        resp_c(i,3) = 0.
      enddo

      return
   end subroutine o2org
    
   subroutine orgc (rainorg, rc, zsed, dbpls, dbmin, pore, zrct, delz, &
      orggg, orgml, smrct, rmserr, kmax, nzmax)
      implicit none

      integer i, j, k, kmax, nzmax

      ! arguments
      real(rk) rainorg, rc, zrct, dbpls(nzmax)
      real(rk) dbmin(nzmax), pore(nzmax), zsed(kmax), delz(kmax)

      ! results
      real(rk) orggg(nzmax), orgml(nzmax), smrct
      real(rk) rmserr

      ! local variables
      real(rk) res(kmax), dres(kmax,3), react(kmax)
      real(rk) dreac(kmax)

      ! arrays for tridiag
      real(rk) a(kmax), b(kmax), c(kmax), r(kmax)
      real(rk) u(kmax), bet, gam(kmax-1)

      smrct = 0.

      do i=2,kmax
        if (zsed(i) .le. zrct .and. orggg(i) .gt. 0) then
          react(i) = - rc*orggg(i)* 3.15e7_rk
          dreac(i) = - rc*3.15e7_rk
        elseif (zsed(i-1) .le. zrct) then
          react(i) = - rc*orggg(i)*3.15e7_rk*(zrct - zsed(i-1))/delz(i)
          dreac(i) = - rc*3.15e7_rk*(zrct - zsed(i-1))/delz(i)
        elseif (zsed(i-1) .gt. zrct) then
          react(i) = 0._rk
          dreac(i) = 0._rk
        endif
      enddo

      do i=3,kmax-1
        res(i) = dbpls(i)*(orggg(i+1) - orggg(i)) - dbmin(i)*(orggg(i) - orggg(i-1)) + react(i)
        dres(i,1) = dbpls(i)
        dres(i,2) = -dbpls(i) - dbmin(i) + dreac(i)
        dres(i,3) = dbmin(i)
      enddo

      i=2
      res(i) = dbpls(i)*(orggg(i+1) - orggg(i)) + react(i) + rainorg*12._rk/delz(i)/(1._rk - pore(i))/2.5_rk

      dres(i,1) = dbpls(i)
      dres(i,2) = -dbpls(i) + dreac(i)

      i = kmax
      res(i) = -dbmin(i)*(orggg(i) - orggg(i-1)) + react(i)
      dres(i,2) = -dbmin(i) + dreac(i)
      dres(i,3) = dbmin(i)

      !     set up the residual array
      do i=1,kmax-1
        r(i) = -res(i+1)
      enddo

      !     lower off-diagonal, dri/dxi-1
      do i=1,kmax-2
        a(i+1) = dres(i+2,3)
      enddo

      !     diagonal, dri/dxi
      do i=1,kmax-1
        b(i) = dres(i+1,2)
      enddo

      !     upper off-diagonal, dri/dxi+1
      do i=1,kmax-2
        c(i) = dres(i+1,1)
      enddo

      bet = b(1)
      u(1) = r(1)/bet
      do j=2,kmax-1
        gam(j) = c(j-1)/bet
        bet = b(j) - a(j)*gam(j)
        u(j) =(r(j) - a(j)*u(j-1))/bet
      enddo
      do j=kmax-2,1,-1
        u(j) = u(j) - gam(j+1)*u(j+1)
      enddo

      !     update the concentration array
      rmserr = 0._rk
      do i=1,kmax-1
        orggg(i+1) = orggg(i+1) + u(i)
        rmserr = rmserr + res(i+1)**2._rk
      enddo

      do k=2,kmax
        orgml(k) = orggg(k)*2.5_rk*(1._rk-pore(k))*1000._rk/12._rk
      enddo
      
      do i=2,kmax
        if (zsed(i) .le. zrct) then
          react(i) = -rc*orggg(i)*3.15e7_rk
        elseif (zsed(i-1) .le. zrct) then
          react(i) = -rc*orggg(i)*3.15e7_rk/pore(i)*(zrct - zsed(i-1))/delz(i)
        elseif (zsed(i-1) .gt. zrct) then
          react(i) = 0._rk
        endif
        smrct = smrct + react(i)*2.5_rk*(1._rk - pore(i))*delz(i)
      enddo

      !     into units of moles / cm2 yr
      smrct = -smrct/12._rk

      do i=3,kmax-1
        res(i) = dbpls(i)*(orggg(i+1) - orggg(i)) - dbmin(i)*(orggg(i) - orggg(i-1)) + react(i)
      enddo

      i = 2
      res(i) = dbpls(i)*(orggg(i+1) - orggg(i)) + react(i) + rainorg*12._rk/delz(i)/(1._rk - pore(i))/2.5_rk
      
      i = kmax
      res(i) = -dbmin(i)*(orggg(i) - orggg(i-1)) + react(i)

      rmserr = 0._rk
      do i=2,kmax
        rmserr = rmserr + res(i)**2._rk
        if (orggg(i) .gt. 1_rk) orggg(i) = 1._rk
        if (orggg(i) .lt. 0_rk) orggg(i) = 0._rk
      enddo
      rmserr = rmserr**(0.5_rk)

      return
   end subroutine orgc
   
   subroutine o2ss (zsed, delz, orgml, pore, dopls, domin, rate, kmax, zrct, o2, smrct, rmserr, nzmax)
      implicit none

      integer j, kmax, nzmax

      !     arguments
      real(rk) zsed(kmax), delz(kmax), pore(nzmax)
      real(rk) orgml(nzmax), dopls(nzmax)
      real(rk) domin(nzmax), rate, zrct

      !     results
      real(rk) o2(nzmax), smrct,rmserr

      !    internal arrays
      real(rk) dfplus(kmax), dfzero(kmax), dfmins(kmax)
      real(rk) res(kmax)

      !    arrays for tridiag
      real(rk) a(kmax), b(kmax), c(kmax), r(kmax)
      real(rk) u(kmax), bet, gam(kmax-1)

      do j=2,kmax-1
        res(j) = (dopls(j)*(o2(j+1) - o2(j)) - domin(j)*(o2(j) - o2(j-1)))
        !       units of m / cm2 (total) s
        if (zsed(j) .le. zrct) then
          res(j) = res(j) - 1.3_rk*rate*orgml(j)/pore(j)
        elseif (zsed(j-1) .le. zrct) then
          res(j) = res(j) - 1.3_rk*rate*orgml(j)/pore(j)*(zrct - zsed(j-1))/delz(j)
        endif
        !       units of m / cm2 (total) s
        dfplus(j) = dopls(j)
        dfzero(j) = -dopls(j) - domin(j)
        dfmins(j) = domin(j)
        !       appropriate to units of m / cm2 s
      enddo

      j=kmax
      res(j) = (-domin(j)*(o2(j) - o2(j-1)))
      if (zsed(j) .lt. zrct) then
        res(j) = res(j) - 1.3_rk*rate*orgml(j)/pore(j)
      elseif (zsed(j-1) .lt. zrct) then
        res(j) = res(j) - 1.3_rk*rate*orgml(j)*(zrct - zsed(j-1))/delz(j)
       endif
       dfplus(j) = 0._rk
       dfzero(j) = -domin(j)
       dfmins(j) = domin(j)

      do j=1,kmax-2
        a(j+1) = dfmins(j+2)
        b(j) = dfzero(j+1)
        c(j) = dfplus(j+1)
      enddo

      b(kmax-1) = dfzero(kmax) + dfplus(kmax)

      do j=1,kmax-1
        r(j) = -res(j+1)
      enddo   
      
      bet = b(1)
      u(1) = r(1)/bet
      do j=2,kmax-1
        gam(j) = c(j-1)/bet
        bet = b(j) - a(j)*gam(j)
        u(j) =(r(j) - a(j)*u(j-1))/bet
      enddo
      do j=kmax-2,1,-1
        u(j) = u(j) - gam(j+1)*u(j+1)
      enddo

      smrct = 0._rk
      do j=2,kmax
        o2(j) = o2(j) + u(j-1)
        if (zsed(j) .le. zrct) then
          !         convert from [moles O2/l total second] to [moles O2 / cm2 yr]
          smrct = smrct - 1.3_rk*rate*orgml(j)*3.15e7_rk/1000._rk*delz(j)
        elseif (zsed(j-1) .le. zrct) then
          smrct = smrct - 1.3_rk*rate*orgml(j)*(zrct-zsed(j-1))/delz(j)*3.15e7_rk/1000._rk*delz(j)
         endif
      enddo

      return
    end subroutine o2ss
   
    subroutine co3ss (resp_c, cal_c, dissc, dissn, csat, u1, u2, zsed, delz, form, &
        pore, kmax, dcpls, dcmin, calgg, carb, ttrorg, ttrcal, ttral, ttrtc, &
        difal, diftc, nzmax, loop_limit, wrong, l)

      implicit none

      integer i, i_l, i_loop, j, k, kmax, l, loop_limit, nzmax

      !     arguments
      real(rk) resp_c(nzmax,3), csat,u1, u2
      real(rk) zsed(kmax),delz(kmax), form(nzmax)
      real(rk) pore(nzmax), dcpls(nzmax,3)
      real(rk) dcmin(nzmax,3), dissc, dissn

      logical wrong, wrong_new

      !     results
      real(rk) cal_c(nzmax), calgg(nzmax)
      real(rk) carb(nzmax,3), ttrorg, ttrcal
      real(rk) ttral,ttrtc, difal, diftc

      real(rk) alkbal, tcbal

      !     internals
      logical done

      l = 1
      done = .false.
      do while (l .le. loop_limit .and. .not. done)
        l = l + 1
        
        call co3 (resp_c, dissc, dissn, csat, u1, u2, dcpls, dcmin, form, pore, kmax, &
            calgg, carb, cal_c, nzmax, wrong)
     
        call sed_diag (ttral, difal, ttrtc, diftc, ttrcal, ttrorg, resp_c, cal_c, &
            carb, pore, delz, dcmin, kmax, nzmax)

        wrong = .false.

        tcbal = abs(1._rk - abs(ttrtc/(diftc + 1.e-20_rk)))
        alkbal = abs(1._rk - abs(ttral/(difal + 1.e-20_rk)))

        if (((alkbal .gt. 0.01_rk) .and. (ttral .gt. 1.e-12_rk)) .or. (tcbal .gt. 0.01_rk)) then
          wrong = .true.
        endif

        if ((ttral .lt. 0._rk) .or. (carb(kmax,3) .lt. 0._rk) .or. (carb(kmax,3) .gt. 500.e-6_rk)) then
            carb(2,3) = (2._rk*carb(1,3) + (csat - 15.e-6_rk))/3._rk
            carb(3,3) = (carb(1,3) + 2._rk*(csat - 15.e-6_rk))/3._rk
            
            do k=4,kmax
              carb(k,3) = csat - 10.e-6_rk
            enddo
            
            do k=2,kmax
              carb(k,1) = carb(1,1)
              carb(k,2) = carb(1,2)
            enddo
        endif

        if (.not. wrong) done = .true.
      enddo

      if (wrong) then
          ttrcal = 0._rk
      endif

      return
    end subroutine co3ss
    
    subroutine co3 (resp_c, dissc, dissn, csat, disoc1, disoc2, dplus, dminus, &
        form, pore, kmax, calgg, carb, cal_c, nzmax, wrong)
      !-----------------------------------------------------------------------
      !     routine co3, which calculates a single iteration
      !     of the carbonate system chemistry.  must be run
      !     several times because of the non-linearity of
      !     calcite dissolution kinetics.
      !-----------------------------------------------------------------------
      implicit none

      integer i, ib, i_loop, m, info, i_row
      integer i_col, im, k, kmax, l, lda, nzmax

      parameter (lda = 16)

!     arguments
      real(rk) resp_c(nzmax,3), csat, disoc1
      real(rk) disoc2, dplus(nzmax,3), dminus(nzmax,3)
      real(rk) form(nzmax), pore(nzmax)
      real(rk) calgg(nzmax), dissc, dissn

      logical wrong

!     results
      real(rk) cal_c(nzmax), carb(nzmax,3)

!     local variables
      real(rk) r(kmax,3), dr(3,3,3,kmax)
      real(rk) rmstc, rmsal, rmsph, weight
      real(rk) trialw

      integer weight_diag

!     linpack variables
      integer ipvt(3*kmax)
      real(rk) bbd(3*kmax), abd(lda,3*kmax)
      
      if (wrong) then
        do i=1,3
!         for the bottom boundary condition, no flux
          carb(kmax+1,i) = carb(kmax,i)
        enddo

!     the residual terms: array (depth; tc, alk, ph)

        do k=2,kmax
!         total co2 equation
          r(k-1,1) = (dplus(k,1)*(carb(k+1,1) - carb(k,1)) - dminus(k,1)*(carb(k,1) - carb(k-1,1))) &
              + (dplus(k,2)*(carb(k+1,2) - carb(k,2)) - dminus(k,2)*(carb(k,2) - carb(k-1,2))) &
              + (dplus(k,3)*(carb(k+1,3) - carb(k,3)) - dminus(k,3)*(carb(k,3) - carb(k-1,3)))
!                     units of moles / l *porewater* sec
          r(k-1,1) = r(k-1,1) + resp_c(k,1)/pore(k) + resp_c(k,2)/pore(k) + resp_c(k,3)/pore(k)
!                     units of moles / l *porewater* sec
          if (carb(k,3) .lt. csat) then
            r(k-1,1) = r(k-1,1) + dissc *((1._rk - (carb(k,3)/csat))**dissn) &
                *(1._rk - pore(k))/pore(k)*calgg(k)*(2.5_rk*1000._rk)/(100._rk)
          endif

!       alkalinity equation
          r(k-1,2) = dplus(k,3)*(carb(k+1,3)-carb(k,3)) - dminus(k,3)*(carb(k,3)-carb(k-1,3)) &
              + 0.5_rk*dplus(k,2)*(carb(k+1,2)-carb(k,2)) - 0.5_rk*dminus(k,2)*(carb(k,2)-carb(k-1,2))
          r(k-1,2) = r(k-1,2) + resp_c(k,3)/pore(k) + 0.5_rk*resp_c(k,2)/pore(k)
          if (carb(k,3) .lt. csat) then
            r(k-1,2) = r(k-1,2) + dissc*((1._rk - (carb(k,3)/csat))**dissn) &
                *(1._rk - pore(k))/pore(k)*calgg(k)*(2.5_rk*1000_rk)/(100_rk)
          endif

          r(k-1,3) = carb(k,1)*carb(k,3)/carb(k,2)**2._rk - disoc2/disoc1
        enddo

!     the derivitive terms: array (function, variable,
!     'k+'= 3 to 'k-' = 1, and depth level k)

        do k=2,kmax-1
          dr(1,1,3,k-1) = dplus(k,1)
          dr(1,1,2,k-1) = -dplus(k,1) - dminus(k,1)
          dr(1,1,1,k-1) = dminus(k,1)

          dr(1,2,3,k-1) = dplus(k,2)
          dr(1,2,2,k-1) = -dplus(k,2) - dminus(k,2)
          dr(1,2,1,k-1) = dminus(k,2)

          dr(1,3,3,k-1) = dplus(k,3)
          dr(1,3,2,k-1) = -dplus(k,3) - dminus(k,3)

          if (carb(k,3) .lt. csat) then
              dr(1,3,2,k-1) = dr(1,3,2,k-1) - dissc*dissn &
              *(1-pore(k))/pore(k)*calgg(k)*(2.5_rk*1000._rk)/100._rk/csat &
              *((1._rk-(carb(k,3)/csat))**(dissn-1._rk))
          endif
          dr(1,3,1,k-1) = dminus(k,3)

          dr(2,1,3,k-1) = 0._rk
          dr(2,1,2,k-1) = 0._rk
          dr(2,1,1,k-1) = 0._rk

          dr(2,2,3,k-1) = 0.5_rk*dplus(k,2)
          dr(2,2,2,k-1) = -0.5_rk*dplus(k,2) - 0.5_rk*dminus(k,2)
          dr(2,2,1,k-1) = 0.5_rk*dminus(k,2)

          dr(2,3,3,k-1) = dplus(k,3)
          dr(2,3,2,k-1) = -dplus(k,3) - dminus(k,3)
          if (carb(k,3) .lt. csat) then
            dr(2,3,2,k-1) = dr(2,3,2,k-1) - dissc*dissn*(1._rk-pore(k))/pore(k)*calgg(k) &
                *(2.5_rk*1000._rk)/100._rk/csat*((1._rk-(carb(k,3)/csat))**(dissn-1._rk))
          endif
          dr(2,3,1,k-1) = dminus(k,3)

          dr(3,1,3,k-1) = 0._rk
          dr(3,1,2,k-1) = carb(k,3)/carb(k,2)**2._rk
          dr(3,1,1,k-1) = 0._rk

          dr(3,2,3,k-1) = 0._rk
          dr(3,2,2,k-1) = -0.5_rk*carb(k,1)*carb(k,3)/carb(k,2)**3._rk
          dr(3,2,1,k-1) = 0._rk

          dr(3,3,3,k-1) = 0._rk
          dr(3,3,2,k-1) = carb(k,1)/carb(k,2)**2._rk
          dr(3,3,1,k-1) = 0._rk
        enddo
      endif

!     bottom special conditions

      do l=1,3 ! function
        do i=1,3 ! variable
          do m=1,3 ! above, below
            if (wrong) then
              dr(l,i,m,kmax-1) = dr(l,i,m,kmax-2)
            endif
          enddo
        enddo
      enddo

!     load the big array

      do k=1,3*kmax
        do l=1,lda
          abd(l,k) = 0._rk
        enddo
      enddo

      if (wrong) then
        do k=2,kmax ! depth level
          do l=1,3 ! function
            do m=1,3 ! up, down
              do i=1,3 !variable
                i_row = (k-2)*3 + l ! row number
                i_col = (m+k-4)*3 + i ! column number
                im = 11
                ib = i_row - i_col + im
                if ((i_col .gt. 0._rk) .and. (i_col .le. (kmax-1)*3)) then
                  abd(ib,i_col) = dr(l,i,m,k-1)
                  if (k .eq. kmax) then
                    if (m .eq. 2) then
                      abd(ib,i_col) = abd(ib,i_col) + dr(l,i,m+1,k-1)
                    endif
                  endif
                endif
              enddo
            enddo
          enddo
        enddo
!     load the residual array

        do l=1,3
          do k=2,kmax
            i_row = (k-2)*3 + l
            bbd(i_row) = -r(k-1,l)
          enddo
        enddo
      
        call my_sgbfa (abd, lda, (kmax-1)*3, 5, 5, ipvt, info, 3*kmax)
        
        call my_sgbsl (abd, lda, (kmax-1)*3, 5, 5, ipvt, bbd)

        weight = 0.5_rk !1.0_rk
        weight_diag = 0

        do k=2,kmax

          i_row = (k-2)*3 + 3
!         CO3= can't go up more than 75%
          trialw = -0.75_rk*carb(k,3)/(bbd(i_row) + 1.e-20_rk)
          if ((trialw .gt. 0._rk) .and. (trialw .lt. weight)) then
            weight = trialw
            weight_diag = 1
          endif

          i_row = (k-2)*3 + 1
!         CO3= can't go up more than 75%
          trialw = -0.75_rk*carb(k,1)/(bbd(i_row) + 1.e-20_rk)
          if ((trialw .gt. 0._rk) .and. (trialw .lt. weight)) then
            weight = trialw
            weight_diag = 2
          endif
        enddo

        do k=2,kmax
          do i=1,3
            i_row = (k-2)*3 + i
            carb(k,i) = carb(k,i) + bbd(i_row)*weight
          enddo
        enddo

        do k=1,kmax
          if (carb(k,3) .lt. csat) then
            cal_c(k) = dissc*((1._rk - (carb(k,3)/csat))**dissn)*(1._rk - pore(k))*calgg(k)*(2.5_rk*1000._rk)/(100._rk)
          else
            cal_c(k) = 0._rk
          endif
        enddo
      endif 

      return
    end subroutine co3 
    
    subroutine my_sgbfa (abd, lda, n, ml, mu, ipvt, info, ncolmax)
!-----------------------------------------------------------------------
!     modified by archer, Oct 92, for vectorized co3 solution
!     no pivoting.   operates in vector form on entire array of
!     similar matrices.

!     expects abd to be filled continuously to first index = n_matrices

!     sgbfa factors a real band matrix by elimination.

!     sgbfa is usually called by sgbco, but it can be called
!     directly with a saving in time if  rcond  is not needed.

!     on entry

!        abd     real(lda, n)
!                contains the matrix in band storage.  the columns
!                of the matrix are stored in the columns of  abd  and
!                the diagonals of the matrix are stored in rows
!                ml+1 through 2*ml+mu+1 of  abd .
!                see the comments below for details.

!        lda     integer
!                the leading dimension of the array  abd .
!                lda must be .ge. 2*ml + mu + 1 .

!        n       integer
!                the order of the original matrix.

!        ml      integer
!                number of diagonals below the main diagonal.
!                0 .le. ml .lt. n .

!        mu      integer
!                number of diagonals above the main diagonal.
!                0 .le. mu .lt. n .
!                more efficient if  ml .le. mu .
!     on return

!        abd     an upper triangular matrix in band storage and
!                the multipliers which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.

!        ipvt    integer(n)
!                an integer vector of pivot indices.

!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that sgbsl will divide by zero if
!                     called.  use  rcond  in sgbco for a reliable
!                     indication of singularity.
!-----------------------------------------------------------------------

      implicit none

      integer lda, n, ml, mu, ipvt(*), info, nmat
      integer ncolmax

      real(rk) abd(lda,ncolmax), t

      integer i, isamax, i0, j, ju, jz, j0, j1, k, kp1, l, lm, m, mm
      integer nm1

      m = ml + mu + 1
      info = 0

!     zero initial fill-in columns
      j0 = mu + 2
      j1 = min0(n,m) - 1
      if (j1 .gt. j0) then
        do jz=j0,j1
          i0 = m + 1 - jz
          do i=i0,ml
            abd(i,jz) = 0._rk
          enddo
        enddo
      endif

      jz = j1
      ju = 0

!     gaussian elimination with partial pivoting
      nm1 = n - 1
      if (nm1 .ge. 1) then

        do k=1,nm1
          kp1 = k + 1

!         zero next fill-in column
          jz = jz + 1
          if ((jz .le. n) .and. (ml .ge. 1)) then
            do i=1,ml
              abd(i,jz) = 0._rk
            enddo
          endif

!         hardwire l = pivot index to no pivoting
          lm = min0(ml,n-k)
          l = m
          ipvt(k) = k

!         zero pivot implies this column already triangularized
!         do test on just one element
          if (abd(l,k) .ne. 0._rk) then
!           compute multipliers
            if (abd(m,k) .eq. 0._rk) write(6,*) 'div by 0', m, k
            t = -1./abd(m,k)

            call my_sscal (lm, t, abd(m+1,k))
            
!           row elimination with column indexing
            ju = min0(max0(ju,mu+ipvt(k)),n)
            mm = m
            if (ju .ge. kp1) then
              do j=kp1,ju
                l = l - 1
                mm = mm - 1
                t = abd(l,j)
                if (l .ne. mm) then
                  abd(l,j) = abd(mm,j)
                  abd(mm,j) = t
                endif
                call my_saxpy (lm, t, abd(m+1,k), abd(mm+1,j))
              enddo
            endif
          else
            info = k
          endif
        enddo
      endif

      ipvt(n) = n
      if (abd(m,n) .eq. 0.) info = n

      return
    end subroutine my_sgbfa
    
    subroutine my_sgbsl (abd, lda, n, ml, mu, ipvt, b)

!-----------------------------------------------------------------------

!     sgbsl solves the real band system
!     a * x = b  or  trans(a) * x = b
!     using the factors computed by sgbco or sgbfa.

!     modified by archer for carbonate solver, Oct 92
!     only job = 0 (solve a*x = b) has been kept

!     on entry

!        abd     real(lda, n)
!                the output from sgbco or sgbfa.

!        lda     integer
!                the leading dimension of the array  abd .

!        n       integer
!                the order of the original matrix.

!        ml      integer
!                number of diagonals below the main diagonal.

!        mu      integer
!                number of diagonals above the main diagonal.

!        ipvt    integer(n)
!                the pivot vector from sgbco or sgbfa.

!        b       real(n)
!                the right hand side vector.

!        job     integer
!                = 0         to solve  a*x = b ,
!                = nonzero   to solve  trans(a)*x = b , where
!                            trans(a)  is the transpose.

!     on return

!        b       the solution vector  x .
!-----------------------------------------------------------------------

      implicit none

      integer lda, n, ml, mu, ipvt(*), job

      real(rk) abd(lda,*), b(*), t

      integer k, kb, l, la, lb, lm, m, nm1

      m = mu + ml + 1
      nm1 = n - 1

!     job = 0 , solve  a * x = b
!     first solve l*y = b
      if ((ml .ne. 0) .and. (nm1 .ge. 1)) then
        do k=1,nm1
          lm = min0(ml,n-k)
          l = ipvt(k)
          t = b(l)
          if (l .ne. k) then
            b(l) = b(k)
            b(k) = t
          endif
          call my_saxpy (lm, t, abd(m+1,k), b(k+1))
        enddo
      endif

!     now solve  u*x = y
      do kb=1,n
        k = n + 1 - kb
        b(k) = b(k)/abd(m,k)
        lm = min0(k,m) - 1
        la = m - lm
        lb = k - lm
        t = -b(k)
        call my_saxpy (lm, t, abd(la,k), b(lb))
      enddo

      return
    end subroutine my_sgbsl
    
    subroutine my_sscal (n, sa, sx)

!-----------------------------------------------------------------------
!     scales a vector by a constant.
!-----------------------------------------------------------------------

      implicit none

      integer i, n

      real(rk) sa, sx(*)

      if (n .le. 0) return
      do i=1,n
        sx(i) = sa*sx(i)
      enddo

      return
    end subroutine my_sscal
    
    subroutine my_saxpy (n, sa, sx, sy)

!-----------------------------------------------------------------------
!     constant times a vector plus a vector.
!-----------------------------------------------------------------------

      implicit none

      integer i, n

      real(rk) sx(*), sy(*), sa

      if (n .le. 0) return
      do i=1,n
        sy(i) = sy(i) + sa*sx(i)
      enddo

      return
    end subroutine my_saxpy
    
    subroutine sed_diag (ttral, difal, ttrtc, diftc, ttrcal, ttrorg, resp_c, &
        cal_c, carb, pore, delz, dcmin, kmax, nzmax)

!-----------------------------------------------------------------------
!     file 'diag.for', which calculates the diffusive fluxes of
!     o2, total co2, and alkalinity at the sediment-water
!     interface, and also the integrated reaction rates of
!     those quantities.  used by co3main to determine when to
!     stop repeating the co3 subroutine.
!-----------------------------------------------------------------------

      implicit none

      integer j, k, kmax, nzmax

!     arguments
      real(rk) resp_c(nzmax,3), cal_c(nzmax)
      real(rk) carb(nzmax,3), pore(nzmax), delz(kmax)
      real(rk) dcmin(nzmax,3)

!     results
      real(rk) ttral, difal, ttrtc, diftc
      real(rk) ttrcal,ttrorg

!     internal arrays
      real(rk) ttreac(3), difflx(3)

!     zero the diagnostics variables, ttreac and flux
      do j=1,3
        ttreac(j) = 0._rk
        difflx(j) = 0._rk
        ttrcal = 0._rk
        ttrorg = 0._rk
      enddo

!     reaction rates are in units of mol species/cm2 (total) y
      do k=1,kmax
!       replaced 3.15e7/1e3 with 3.15e4
        ttreac(1) = ttreac(1) + resp_c(k,1)*delz(k)*3.15e4_rk
        ttreac(2) = ttreac(2) + resp_c(k,2)*delz(k)*3.15e4_rk
        ttreac(3) = ttreac(3) + (resp_c(k,3) + cal_c(k))*delz(k)*3.15e4_rk
        ttrcal = ttrcal + cal_c(k)*delz(k)*3.15e4_rk
        ttrorg = ttrorg + (resp_c(k,1) + resp_c(k,2) + resp_c(k,3))*delz(k)*3.15e4_rk
      enddo

!     the diffusive fluxes

!     replaced 3.15e7/1e3 with 3.15e4
      difflx(1) = dcmin(2,1)*(carb(1,1) - carb(2,1))*pore(2)*delz(2)*3.15e4_rk
      difflx(2) = dcmin(2,2)*(carb(1,2) - carb(2,2))*pore(2)*delz(2)*3.15e4_rk
      difflx(3) = dcmin(2,3)*(carb(1,3) - carb(2,3))*pore(2)*delz(2)*3.15e4_rk

      ttrtc = ttreac(1) + ttreac(3)
      ttral = ttreac(3)*2._rk
      diftc = difflx(1) + difflx(2) + difflx(3)
      difal = difflx(2) + difflx(3)*2._rk

      return
    end subroutine sed_diag
    
    subroutine bury (zsed, delz, pore, kmax, calgg, orggg, sed_ml_mass, &
        rain_cal_p, ttrcal, rain_clay_p, buried_mass, buried_calfrac, &
        nzmax, ibmax, counter)

      implicit none

      integer ibmax, k, counter
      integer kmax, n_coretop, n_depth, nzmax, n_accum!, nyear

!     arguments
      real(rk) zsed(nzmax), delz(nzmax), pore(nzmax)
      real(rk) calgg(nzmax), orggg(nzmax)
      real(rk) sed_ml_mass, rain_cal_p, ttrcal
      real(rk) rain_clay_p, exp_pore, pore_max

!     results
      real(rk) buried_mass(ibmax), buried_calfrac(ibmax)

!     internal variables (temporary)
      real(rk) buried_mass_step, calgg_est
      real(rk) sed_ml_mass_new, up_cal_frac, cal_total_new
            
!     estimate the change in calgg, then calculate new porosities
      calgg_est = (calgg(kmax)*sed_ml_mass + rain_cal_p*100._rk - ttrcal*100._rk)/sed_ml_mass
     
      pore_max = 1._rk - (0.483_rk + 0.45_rk*calgg_est)/2.5_rk
      exp_pore = 0.25_rk*calgg_est + 3._rk*(1._rk-calgg_est)
      do k=2,kmax
        pore(k) = exp(-zsed(k)/exp_pore)*(1._rk-pore_max) + pore_max
      enddo
     
      sed_ml_mass_new = 0._rk
      do k=1,kmax
        sed_ml_mass_new = sed_ml_mass_new+delz(k)*(1._rk-pore(k))*2.5_rk
      enddo

      n_depth = 2
      !do k=ibmax-2,1,-1                                                      ! NIC: seems not to be used in standard UVic setups. 
      !  if (float(nyear) - depth_age(k) .gt. -1.e-6_rk) n_depth = k + 1      ! NIC: seems not to be used in standard UVic setups. 
      !enddo                                                                  ! NIC: seems not to be used in standard UVic setups. 

!     calculate mass buried this time step
!     units of g/cm2 per this time step
      buried_mass_step = rain_cal_p*100._rk - ttrcal*100._rk + rain_clay_p + sed_ml_mass - sed_ml_mass_new

!     now do some burial
      do k=1,n_depth
        if (buried_mass(k) .gt. 0._rk) n_coretop = k
      enddo

!     find the "latest" box with mass in it, one ip at a time
      if (buried_mass_step .gt. 0.) then
!       then we're accumulating -- that's easy

!       add buried calcite to the sediment record
        buried_calfrac(n_depth) = (buried_mass(n_depth)*buried_calfrac(n_depth) + buried_mass_step &
            *calgg(kmax))/(buried_mass(n_depth) + buried_mass_step)

!       add buried mass to sediment record
        buried_mass(n_depth) = buried_mass(n_depth) + buried_mass_step

!       adjust calgg
        do k=1,kmax
          calgg(k) = (calgg(k)*sed_ml_mass + rain_cal_p*100._rk - &
              ttrcal*100._rk - buried_mass_step*calgg(kmax))/sed_ml_mass_new
        enddo

!       all business is complete at this point

      else
!       or else there's chemical erosion -- bms < 0

        if (-buried_mass_step .gt. buried_mass(n_coretop)) then
!         then we're using up the entire coretop box by erosion
!         assume that this can only happen once per year timestep.

!         calculate the fraction of calcite in the upwelling
!         material by combining the fractions of boxes n_coretop
!         (until it's used up) and n_coretop - 1 (below top)

!         write(6,*) "into box depletion in bury"

          up_cal_frac = (buried_mass(n_coretop)*buried_calfrac(n_coretop) + ((-buried_mass_step - &
              buried_mass(n_coretop))*buried_calfrac(n_coretop-1)))/(-buried_mass_step)

!         adjust the sediment record to reflect erosion
          buried_mass(n_coretop-1) = buried_mass(n_coretop-1) + (buried_mass_step + buried_mass(n_coretop))
          buried_mass(n_coretop) = 0._rk
          n_coretop = n_coretop - 1

        else
!         or else erosion can come entirely out of the box
!         n_coretop, and the upwelling calcite fraction is
!         simply  buried_calfrac(n_coretop,ip)

          up_cal_frac = buried_calfrac(n_coretop)

!         adjust the sediment record

          buried_mass(n_coretop) = buried_mass(n_coretop) + buried_mass_step

        endif

!       now gotta adjust calgg
        cal_total_new = sed_ml_mass*calgg(kmax) + rain_cal_p*100._rk - ttrcal*100._rk - buried_mass_step*up_cal_frac

        do k=1,kmax
          calgg(k) = cal_total_new/sed_ml_mass_new
        enddo
      endif
!     end of the "erosion if"

      do k=1,kmax
        if (calgg(k) .lt. 0._rk) then
          calgg(k) = 1.e-6_rk
        endif
      enddo

      sed_ml_mass = sed_ml_mass_new

      return
    end subroutine bury
end module