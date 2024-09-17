#include "fabm_driver.h"

module uvic_nut_chem
   use fabm_types
   use fabm_particle
   use fabm_builtin_models!, only: type_depth_integral
   
   use uvic_shared

   implicit none

   private

   type, extends(type_particle_model), public :: type_uvic_nut_chem
      type (type_state_variable_id)                            :: id_det, id_alk, id_o2, id_phosphorus, id_no3, id_dic, id_c14!id_det, id_o2, id_phosphorus, id_no3
      type (type_dependency_id)                                :: id_o2_sms, id_phosphorus_sms, id_depth, id_dzt, id_calc_sms, id_sss, id_sst!, id_temp, id_depth, id_wd_in
      type (type_surface_dependency_id)                        :: id_intcalc, id_ws, id_aice, id_dc14ccn, id_atco2
      type (type_bottom_dependency_id)                         :: id_bdepth
      type (type_diagnostic_variable_id)                       :: id_o2_dummy, id_phosphorus_dummy, id_calc_dummy, id_deni
      type (type_diagnostic_variable_id)                       :: id_rain_cal
      
      real(rk)                            :: dcaco3, rstd
      logical                             :: nitrogen
   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_column
      procedure :: do_surface
   end type

contains

   subroutine initialize(self, configunit)
      class (type_uvic_nut_chem), intent(inout), target :: self
      integer,                    intent(in)            :: configunit
      class (type_depth_integral), pointer :: calc_sms_integrator
      allocate(calc_sms_integrator)
      
      ! parameters
      call self%get_parameter(self%nitrogen,     'nitrogen',     '',            'turn on nitrogen dependency',                    default=.true.)
      call self%get_parameter(self%dcaco3,       'dcaco3',       'cm',          'calcite remineralization depth',                 default=650000.0_rk, scale_factor=1.e-2_rk)
      call self%get_parameter(self%rstd,         'rstd',         '',            'standard c14/c12 ratio',                         default=1.176e-12_rk)

      call self%register_state_variable(self%id_o2,         'o2',         'umol O cm-3',   'oxygen concentration')
      call self%register_state_variable(self%id_phosphorus, 'p',          'mmol P m-3',    'dissolved inorganic phosphorous concentration')
      call self%register_state_variable(self%id_dic,        'dic',        'umol C cm-3',    'dissolved inorganic carbon concentration')
      call self%register_state_variable(self%id_alk,        'alkalinity', 'umol cm-3',     'ocean alkalinity')
      call self%register_state_variable(self%id_c14,        'c14',        'umol c14 cm-3', 'carbon 14 concentration')
      if (self%nitrogen) then
          call self%register_state_variable(self%id_no3,    'no3',      'mmol n m-3',  'dissolved inorganic nitrogen concentration')
          call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_no3)
          call self%register_diagnostic_variable(self%id_deni, 'denitrification',  'mmol m-3 s-1', 'denitrification rate')
      endif
      
      call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_phosphorus)
      
      ! set up dummy state vars to pick up sources from the npzd model
      call self%register_dependency(self%id_o2_sms,  'oxygen_sms',   'umol cm-3 s-1',  'oxygen sms')
      call self%register_diagnostic_variable(self%id_o2_dummy, 'o2_dummy','umol cm-3', 'target pool for net respiration',source=source_constant, output=output_none, act_as_state_variable=.true.)
      call self%request_coupling(self%id_o2_sms,'./o2_dummy_sms_tot')
      
      call self%register_dependency(self%id_calc_sms,  'calcite_sms',   'umol cm-3 s-1',  'calc sms')
      call self%register_diagnostic_variable(self%id_calc_dummy, 'calc_dummy','umol cm-3', 'target pool for calcite production',source=source_constant, output=output_none, act_as_state_variable=.true.)
      call self%request_coupling(self%id_calc_sms,'./calc_dummy_sms_tot')

      call self%register_dependency(self%id_phosphorus_sms,  'phosphorus_sms',   'mmol P m-3 s-1',  'phosphorus sms')
!      call self%register_diagnostic_variable(self%id_phosphorus_dummy, 'phosphorus_dummy','mmol P m-3', 'target pool for p release',source=source_constant, output=output_none, act_as_state_variable=.true.)
      call self%request_coupling(self%id_phosphorus_sms,'./p_sms_tot')
      self%id_phosphorus%sms%link%target%source = source_constant
      
      ! register calcite bottom flux diagnostic
      call self%register_diagnostic_variable(self%id_rain_cal,  'rain_cal', 'umol cm-2 s-1', 'calcite sedimentation flux', source=source_do_column)
      
      ! set up depth integral of calcite production
      call self%add_child(calc_sms_integrator, 'calc_sms_integrator')
      call calc_sms_integrator%request_coupling('source','../calc_dummy_sms_tot')
      
      ! couple depth-integrated calcite production dependency to calc_sms_integrator 
      call self%register_dependency(self%id_intcalc, 'intcalc', 'umol cm-2', 'depth-integrated calcite production')
      call self%request_coupling(self%id_intcalc, './calc_sms_integrator/result')
      
      ! environmental dependencies
      call self%register_dependency(self%id_depth,   standard_variables%depth)
      call self%register_dependency(self%id_dzt,     standard_variables%cell_thickness)
      call self%register_dependency(self%id_bdepth,  standard_variables%bottom_depth)
      call self%register_dependency(self%id_atco2,   standard_variables%mole_fraction_of_carbon_dioxide_in_air) ! 'atco2', '1e-6', 'atmospheric carbon concentration')
      call self%register_dependency(self%id_sst,     standard_variables%temperature)
      call self%register_dependency(self%id_sss,     standard_variables%practical_salinity)
      call self%register_dependency(self%id_ws,      standard_variables%wind_speed)
      call self%register_dependency(self%id_aice,    'aice', '-', 'ice area fraction')
      call self%register_dependency(self%id_dc14ccn, 'dc14ccn',  '1e-6', 'atmospheric C14 concentration')
      

! declare p, calcite, optional n, alk and dic
! get roc of phosphorous and apply to P, alk, and dic
      
   end subroutine initialize
   
   subroutine do(self, _ARGUMENTS_DO_)
      class (type_uvic_nut_chem), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_
      
      real(rk) :: p, p_release, o2, no3, nud, temp, bct, det, o2_demand, fo2, o2_roc, no3flag
      real(rk) :: deni, calc_prod, dic_roc, alk_roc, so2, Paulmier_a, Paulmier_z, Paulmier_R0
      real(rk) :: c14, c14_roc
      
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
      
      real(rk) :: dzt, depth, rcak, rcab, zw, zw_prev, intcalc, bdepth, dic_roc, alk_roc, rain_cal
      _GET_SURFACE_(self%id_intcalc,intcalc)
      intcalc = intcalc*100._rk              ! we multiply with [cm m-1], to convert back to cm-2 (because of unit mismatch, type_depth_integral returns [umol cm-3 m])
      _GET_BOTTOM_(self%id_bdepth,bdepth)

      rcak = 0.0_rk 
      rcab = 0.0_rk
      _DOWNWARD_LOOP_BEGIN_
         _GET_(self%id_dzt,dzt)
         _GET_(self%id_depth,depth)
         dzt = dzt*100._rk ! convert from m to cm
         depth = depth*100._rk ! convert from m to cm
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
             rain_cal = 0.0_rk
         else
             dic_roc = intcalc*rcab
             alk_roc = intcalc*rcab*2._rk
             rain_cal = intcalc*rcab - intcalc*rcak*dzt
         endif
         
         _ADD_SOURCE_(self%id_dic, dic_roc) !NIC: NOT FINISHED we still need to interact with the sediment (line 742 of tracer.f90) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         _ADD_SOURCE_(self%id_alk, alk_roc)
         _SET_DIAGNOSTIC_(self%id_rain_cal,rain_cal)
      _DOWNWARD_LOOP_END_
      
   end subroutine do_column
   
   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_uvic_nut_chem), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: dic_in, ta_in, o2_in, co2_in, t_in, s_in, pt_in, sit_in, atmpres, pHlo, pHhi, co2star
      real(rk) :: dco2star, pCO2, hplus, CO3, Omega_c, Omega_a, ws, dicflx, c14, aice, ao, dc14ccn, sspH
      real(rk) :: ssCO3, ssOc, ssOa, sspCO2, scco2, piston_vel, c14flx, f1, f2, f3, f4, f5, o2sat, o2flx
      real(rk) :: sco2, piston_o2
      
      
      _SURFACE_LOOP_BEGIN_
         _GET_(self%id_dic, dic_in)
         _GET_(self%id_o2, o2_in)
         _GET_(self%id_alk, ta_in)
         _GET_SURFACE_(self%id_atco2, co2_in)
         _GET_SURFACE_(self%id_dc14ccn, dc14ccn)
         _GET_(self%id_sst, t_in)
         _GET_(self%id_sss, s_in)
         _GET_SURFACE_(self%id_ws, ws)
         _GET_SURFACE_(self%id_aice, aice)
         ao = 1._rk-aice
         ws = ws*100._rk ! convert from m s-1 to cm s-1
         ! put reasonable limits on sst and sss for chemistry flux calculations
         t_in = min(35._rk,max(t_in,-2._rk))
         s_in = min(45._rk,max(s_in,0._rk))
         
         pt_in = 0.5125e-3_rk     ! mol/m^3
         sit_in = 7.6875e-03_rk   ! mol/m^3
         atmpres = 1.0_rk         ! atm
         pHlo = 6._rk
         pHhi = 10._rk
         
         call co2calc_SWS (t_in, s_in, dic_in, ta_in, co2_in, pt_in, &
             sit_in, atmpres, pHlo, pHhi, &
             co2star, dco2star, pCO2, hplus, CO3, &
             Omega_c, Omega_a)
         
         ! sspH is really H+ until the final averaging
         sspH = hplus
         ssCO3 = CO3
         ssOc = Omega_c
         ssOa = Omega_a
         sspCO2 = pCO2
                  
         ! Schmidt number for CO2
         scco2 = 2073.1_rk - 125.62_rk*t_in + 3.6276_rk*t_in**2._rk &
             - 0.043219_rk*t_in**3._rk
         piston_vel = ao*xconv*((ws*0.01_rk)**2._rk) &
             *((scco2/660._rk)**(-0.5_rk))
         ! dic in umol cm-3 or (mol m-3) => flux in umol cm-2 s-1
         dicflx = piston_vel*dco2star
         
         c14flx = piston_vel*((dco2star + co2star) &
             *(1._rk + dc14ccn*0.001_rk)*self%rstd &
             - co2star*c14/dic_in)
         
         !-----------------------------------------------------------------------
         !           calculate ocean oxygen fluxes
         !-----------------------------------------------------------------------
         ! Schmidt number for O2
         sco2 = 1638.0_rk - 81.83_rk*t_in + 1.483_rk*t_in**2._rk - 0.008004_rk*t_in**3._rk
         ! piston velocity for O2
         piston_o2 = ao*xconv*((ws*0.01_rk)**2._rk) &
             *(sco2/660.0_rk)**(-0.5_rk)
         ! oxygen saturation concentration [mol/m^3]
         f1 = log((298.15_rk - t_in)/(C2K + t_in))
         f2 = f1*f1
         f3 = f2*f1
         f4 = f3*f1
         f5 = f4*f1
         o2sat = exp (2.00907_rk + 3.22014_rk*f1 + 4.05010_rk*f2 &
             + 4.94457_rk*f3 - 2.56847E-1_rk*f4 + 3.88767_rk*f5  &
             + s_in*(-6.24523e-3_rk - 7.37614e-3_rk*f1 - 1.03410e-2_rk*f2 &
             - 8.17083E-3_rk*f3) - 4.88682E-7_rk*s_in*s_in)
         ! Convert from ml/l to mol/m^3
         o2sat = o2sat/22391.6_rk*1000.0_rk
         
         o2flx = piston_o2*(o2sat - o2_in)

         _ADD_SURFACE_FLUX_(self%id_o2, o2flx)  
         _ADD_SURFACE_FLUX_(self%id_dic, dicflx)  
         _ADD_SURFACE_FLUX_(self%id_c14, c14flx)  
      _SURFACE_LOOP_END_
   end subroutine do_surface
   
   subroutine co2calc_SWS (t, s, dic_in, ta_in, xco2_in, pt_in, &
       sit_in, atmpres, phlo, phhi, &
       co2star, dco2star, pco2, hSWS, &
       CO3, Omega_c, Omega_a)
      !-------------------------------------------------------------------------
      ! Modified from co2calc.f (RCS version 1.8, OCMIP-2)
      ! Constants are given on "seawater" H scale (hSWS) except for the "free" 
      ! H scale dissociation constants for S and F (necessary to work with hSWS).
      
      ! PURPOSE: Calculate delta co2* (dco2star)
      
      ! INPUT:
      !        t       = temperature (degrees C)
      !        s       = salinity (PSU)
      !        dic_in  = total inorganic carbon (mol/m^3)
      !        ta_in   = total alkalinity (eq/m^3)
      !        xco2_in = atmospheric CO2 concentration (ppmv)
      !        pt_in   = inorganic phosphate (mol/m^3)
      !        sit_in  = inorganic silicate (mol/m^3)
      !        atmpres = atmospheric pressure in atmospheres (1 atm==1013.25mbar)
      !        depth   = ocean depth (m)
      !        phlo    = lower limit of pH range
      !        phhi    = upper limit of pH range
      
      ! OUTPUT:
      !        co2star  = ocean CO2* (mol/m^3)
      !        dco2star = delta (atm-ocn) CO2* (mol/m^3)
      !        pco2     = oceanic pCO2 (uatm)
      !        hsws     = H+ on seawater H scale
      !--------------------------------------------------------------------------

      implicit none

      real(rk) permil, pt, pt_in, sit, sit_in, ta, ta_in, dic, dic_in
      real(rk) permeg, xco2, xco2_in, tk, t, tk100, tk1002, dlogtk, s
      real(rk) sqrtis, s2, sqrts, s15, scl, bt, st, ft, ff, x1, phhi
      real(rk) x2, phlo, xacc, hSWS, hSWS2, co2star, co2starair
      real(rk) atmpres, dco2star, ph, pco2, invtk, is, is2, k0, k1
      real(rk) k12, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi
      real(rk) Kspc, Kspa, ca, co3, omega_c, omega_a, depth, pres
      real(rk) t2, dvc, dva, dk, pitkR, p2itkR, Rgas, fugf

      !--------------------------------------------------------------------------
      ! Models carry tracers in mol/m^3 but this routine uses umol/kg
      ! Change units from the input of mol/m^3 -> mol/kg:
      ! (1 mol/m^3) x (1 m^3/1024.5 kg)
      ! where the ocean's mean surface density is 1024.5 kg/m^3

      permil = 1._rk/1024.5_rk ! to convert mol/m^3 to mol/kg
      permeg = 1.e-6_rk        ! to convert uatm to atm.
      Rgas = 83.1451_rk        ! ideal gas constant R in cm3 bar k-1

      ! Convert input
      dic = dic_in*permil   ! change mol/m^3 to mol/kg
      ta = ta_in*permil     ! change mol/m^3 to mol/kg
      xco2 = xco2_in*permeg ! change units from ppmv to parts by volume.
      pt = pt_in*permil     ! change mol/m^3 to mol/kg
      sit = sit_in*permil   ! change mol/m^3 to mol/kg
      pres = depth*0.1_rk   ! change from meters to decibars

      !--------------------------------------------------------------------------
      ! Calculate all constants needed to convert between various measured
      ! carbon species. References for each equation are noted in the code.
      ! Version 2 of "Handbook of Methods for the Analysis of the Various
      ! Parameters of the Carbon Dioxide System in Seawater", DOE, 1994
      ! (SOP No. 3, p25-26).
      
      ! Derive simple terms used more than once
      tk = C2K + t
      tk100 = tk/100._rk
      tk1002 = tk100*tk100
      invtk = 1._rk/tk
      dlogtk = log(tk)
      is = 19.924_rk*s/(1000._rk - 1.005_rk*s)
      is2 = is*is
      sqrtis = sqrt(is)
      s2 = s*s
      t2 = t*t
      sqrts = sqrt(s)
      s15 = s**1.5_rk
      scl = s/1.80655_rk
      pitkR = pres/tk/Rgas
      p2itkR = pres*pitkR

      !------------------------------------------------------------------------
      ! Calculate concentrations for borate, sulfate, and fluoride
      ! Uppstrom (1974)
      bt = 0.000232_rk*scl/10.811_rk
      ! Morris & Riley (1966)
      st = 0.14*scl/96.062_rk
      ! Riley (1965)
      ft = 0.000067_rk*scl/18.9984_rk

      !------------------------------------------------------------------------
      ! f = k0(1-pH2O)*correction term for non-ideality
      ! Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)
      ff = exp(-162.8301_rk + 218.2968_rk/tk100  + 90.9241_rk*log(tk100) &
          - 1.47696_rk*tk1002 + s*(.025695_rk - .025225_rk*tk100 &
          + 0.0049867_rk*tk1002))

      !------------------------------------------------------------------------
      ! k1 = [H][HCO3]/[H2CO3]
      ! k2 = [H][CO3]/[HCO3] on hSWS
      ! Millero p.664 (1995) using Mehrbach et al. data on SEAWATER scale
      ! (Original reference: Dickson and Millero, DSR, 1987)
      k1 = 10._rk**(-1._rk*(3670.7_rk*invtk - 62.008_rk + 9.7944_rk*dlogtk &
          - 0.0118_rk*s + 0.000116_rk*s2))
      k2 = 10_rk**(-1._rk*(1394.7_rk*invtk + 4.777_rk - 0.0184_rk*s + 0.000118_rk*s2))

      !------------------------------------------------------------------------
      ! k1p = [H][H2PO4]/[H3PO4] on hSWS
      ! Millero p.670 (1995)
      k1p = exp(-4576.752_rk*invtk + 115.540_rk - 18.453_rk*dlogtk &
          + (-106.736_rk*invtk + 0.69171_rk)*sqrts & 
          + (-0.65643_rk*invtk - 0.01844_rk)*s)

      !------------------------------------------------------------------------
      ! k2p = [H][HPO4]/[H2PO4] on hSWS
      ! Millero p.670 (1995)
      k2p = exp(-8814.715_rk*invtk + 172.1033_rk - 27.927_rk*dlogtk &
          + (-160.340_rk*invtk + 1.3566_rk)*sqrts &
          + (0.37335_rk*invtk - 0.05778_rk)*s)

      !------------------------------------------------------------------------
      ! k3p = [H][PO4]/[HPO4] on hSWS
      ! Millero p.670 (1995)
      k3p = exp(-3070.75_rk*invtk - 18.126_rk &
          + (17.27039_rk*invtk + 2.81197_rk)*sqrts &
          + (-44.99486_rk*invtk - 0.09984_rk)*s)

      !------------------------------------------------------------------------
      ! ksi = [H][SiO(OH)3]/[Si(OH)4] on hSWS
      ! Millero p.671 (1995) using data from Yao and Millero (1995)
      ! change to (mol/ kg soln)
      ! depth dependancy assumed to be the same as boric acid
      ! typo in Millero 1994 corrected in sign of 0.1622

      ksi = exp(-8904.2_rk*invtk + 117.400_rk - 19.334_rk*dlogtk &
          + (-458.79_rk*invtk + 3.5913_rk)*sqrtis &
          + (188.74_rk*invtk - 1.5998_rk)*is &
          + (-12.1652_rk*invtk + 0.07871_rk)*is2 &
          + log(1._rk - 0.001005_rk*s))

      !------------------------------------------------------------------------
      ! kw = [H][OH] on hSWS
      ! Millero p.670 (1995) using composite data
      ! pressure dependancy in Millero 1994 corrected for sea water from
      ! Millero 1983
      kw = exp(-13847.26_rk*invtk + 148.9802_rk - 23.6521_rk*dlogtk &
          + (118.67_rk*invtk - 5.977_rk + 1.0495_rk*dlogtk)*sqrts - 0.01615_rk*s)

      !------------------------------------------------------------------------
      ! ks = [H][SO4]/[HSO4] on free H scale
      ! Dickson (1990, J. chem. Thermodynamics 22, 113)
      ! change to (mol/kg soln)
      ks = exp(-4276.1_rk*invtk + 141.328_rk - 23.093_rk*dlogtk &
          + (-13856._rk*invtk + 324.57_rk - 47.986_rk*dlogtk)*sqrtis &
          + (35474._rk*invtk - 771.54_rk + 114.723_rk*dlogtk)*is &
          - 2698._rk*invtk*is**1.5_rk + 1776._rk*invtk*is2 &
          + log(1._rk - 0.001005_rk*s))

      !------------------------------------------------------------------------
      ! kf = [H][F]/[HF] on free H scale
      ! Dickson and Riley (1979)
      ! change to (mol/ kg soln)
      kf = exp(1590.2_rk*invtk - 12.641_rk + 1.525_rk*sqrtis &
          + log(1._rk - 0.001005_rk*s))

      !------------------------------------------------------------------------
      ! kb = [H][BO2]/[HBO2] on hSWS
      ! Dickson p.673 (1990)
      ! change from htotal to hSWS
      ! typo in Millero 1994 corrected in sign of 0.1622
      kb = exp((-8966.90_rk - 2890.53_rk*sqrts - 77.942_rk*s &
          + 1.728_rk*s15 - 0.0996_rk*s2)*invtk &
          + (148.0248_rk + 137.1942_rk*sqrts + 1.62142_rk*s) &
          + (-24.4344_rk - 25.085_rk*sqrts - 0.2474_rk*s)*dlogtk &
          + 0.053105_rk*sqrts*tk &
          + log((1._rk+(st/ks)+(ft/kf))/(1._rk+(st/ks))))

      !-------------------------------------------------------------------------
      ! Calculate [H+] SWS when DIC and TA are known at T, S and 1 atm.
      
      ! The solution converges to err of xacc. The solution must be within
      ! the range x1 to x2.
      
      ! If DIC and TA are known then either a root finding or iterative method
      ! must be used to calculate hSWS. In this case we use the Newton-Raphson
      ! "safe" method taken from "Numerical Recipes" (function "rtsafe.f" with
      ! error trapping removed).
      
      ! As currently set, this procedure iterates about 12 times. The x1 and x2
      ! values set below will accomodate ANY oceanographic values. If an initial
      ! guess of the pH is known, then the number of iterations can be reduced to
      ! about 5 by narrowing the gap between x1 and x2. It is recommended that
      ! the first few time steps be run with x1 and x2 set as below. After that,
      ! set x1 and x2 to the previous value of the pH +/- ~0.5. The current
      ! setting of xacc will result in co2star accurate to 3 significant figures
      ! (xx.y). Making xacc bigger will result in faster convergence also, but this
      ! is not recommended (xacc of 10**-9 drops precision to 2 significant figures).

      x1 = 10._rk**(-phhi)
      x2 = 10._rk**(-phlo)
      xacc = 1.e-10_rk
      hSWS = drtsafe (k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi, bt, st, &
          ft, sit, pt, dic, ta, x1, x2, xacc)

      !-------------------------------------------------------------------------
      ! Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2,
      ! ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
      hSWS2 = hSWS*hSWS
      co2star = (dic*hSWS2/(hSWS2 + k1*hSWS + k1*k2))
      co2starair = xco2*ff*atmpres
      dco2star = co2starair - co2star

      ph = -log10(hSWS)
      ! Calculate [CO3] in mol/kg from co2star and hSWS
      CO3 = k1*k2*co2star/hSWS2
      ! K0 and fugf from Weiss 1974 Marine Chemistry
      k0 = exp(93.4517_rk/tk100 - 60.2409_rk + 23.3585_rk*log(tk100) &
          + s*(.023517_rk - 0.023656_rk*tk100 + 0.0047036_rk*tk1002))
      ! Fugacity factor (assumes pressure of 1 atm or 1.01325 bars)
      fugf = exp((-1636.75_rk + 12.0408_rk*tk - 0.0327957_rk*tk**2._rk &
          + 3.16528e-5_rk*tk**3._rk + 2._rk*(57.7_rk - 0.118_rk*tk))*1.01325_rk/(Rgas*tk))
      ! Calculate the partial pressures
      pCO2 = co2star/(k0*fugf)/permeg ! convert atm to uatm

      !-------------------------------------------------------------------------
      ! Solubility products of calcite and aragonite at sea level
      ! CalciteSolubility:
      ! Mucci, Alphonso, Amer. J. of Science 283:781-799, 1983.
      Kspc = 10._rk**(-171.9065_rk - 0.077993_rk*tk + 2839.319_rk/tk &
          + 71.595_rk*dlogtk/log(10._rk) + (-0.77712_rk + 0.0028426_rk*tk &
          + 178.34_rk/tk)*sqrts - 0.07711_rk*s + 0.0041249_rk*sqrts*s)

      ! Aragonite Solubility
      ! Mucci, Alphonso, Amer. J. of Science 283:781-799, 1983
      Kspa = 10._rk**(-171.945_rk - 0.077993_rk*tk + 2903.293_rk/tk &
          + 71.595_rk*dlogtk/log(10._rk) &
          + (-0.068393_rk + 0.0017276_rk*tk + 88.135_rk/tk)*sqrts &
          - 0.10018_rk*s + 0.0059415_rk*sqrts*s)


      ! Calculate saturation state (Omega)
      ! Riley, J. P. and Tongudai, M., Chemical Geology 2:263-269, 1967
      Ca = (0.02128_rk/40.087_rk)*(s/1.80655_rk)

      ! calcite (Omgc)
      Omega_c = Ca*CO3/Kspc

      ! aragonite (Omga)
      Omega_a = Ca*CO3/Kspa

      CO3 = CO3/permil ! convert from mol/kg -> mol/m^3


      !-------------------------------------------------------------------------
      co2star = co2star/permil ! convert from mol/kg -> mol/m^3
      dco2star = dco2star/permil ! convert from mol/kg -> mol/m^3
 
      return
   end subroutine co2calc_SWS


   real(rk) function drtsafe (k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi, &
       bt, st, ft, sit, pt, dic, ta, x1, x2, xacc)
      !-------------------------------------------------------------------------
      !     File taken from Numerical Recipes. Modified  R.M.Key 4/94

      implicit none

      integer maxit, j

      real(rk) k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi, bt, st, ft, sit
      real(rk) pt, dic, ta
      real(rk) x1, fl, df, x2, fh, xl, xh, swap, dxold, dx, f, temp, xacc

      maxit = 100
      call ta_iter_SWS (k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi, &
          bt, st, ft, sit, pt, dic, ta, x1, fl, df)
      call ta_iter_SWS (k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi, &
          bt, st, ft, sit, pt, dic, ta, x2, fh, df)
      if (fl .lt. 0._rk) then
        xl = x1
        xh = x2
      else
        xh = x1
        xl = x2
        swap = fl
        fl = fh
        fh = swap
      endif
      drtsafe = 0.5_rk*(x1 + x2)
      dxold = abs(x2 - x1)
      dx = dxold
      call ta_iter_SWS (k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi, &
          bt, st, ft, sit, pt, dic, ta, drtsafe, f, df)
      do j=1,maxit
        if (((drtsafe - xh)*df - f)*((drtsafe - xl)*df - f) .ge. 0._rk &
            .or. abs(2._rk*f) .gt. abs(dxold*df)) then
          dxold = dx
          dx = 0.5_rk*(xh - xl)
          drtsafe = xl + dx
          if (xl .eq. drtsafe) return
        else
          dxold = dx
          dx = f/df
          temp = drtsafe
          drtsafe = drtsafe - dx
          if (temp .eq. drtsafe) return
        endif
        if (abs(dx) .lt. xacc) return
        call ta_iter_SWS (k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi, &
            bt, st, ft, sit, pt, dic, ta, drtsafe, f, df)
        if (f .lt. 0._rk) then
          xl = drtsafe
          fl = f
        else
          xh = drtsafe
          fh = f
        endif
      enddo

      return
   end


   subroutine ta_iter_SWS (k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi, &
       bt, st, ft, sit, pt, dic, ta, x, fn, df)
      !-------------------------------------------------------------------------
      ! This routine expresses TA as a function of DIC, hSWS and constants.
      ! It also calculates the derivative of this function with respect to
      ! hSWS. It is used in the iterative solution for hSWS. In the call
      ! "x" is the input value for hSWS, "fn" is the calculated value for TA
      ! and "df" is the value for dTA/dhSWS
      ! Modified from ta_iter_1.f (RCS version 1.2, OCMIP-2)

      implicit none

      real(rk) k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi, bt, st, ft, sit
      real(rk) pt, dic, ta, x, fn, df
      real(rk) rx, x2, x3, k12, k12p, k123p, c, rc, a, ra, da, b, rb, db
      real(rk) rkb, rksi, kwrx, ptra, xrc, crx, crx2, t1, t2, t3, t4, t5, t6

      rx = 1._rk/x
      x2 = x*x
      x3 = x2*x
      k12 = k1*k2
      k12p = k1p*k2p
      k123p = k12p*k3p
      c = 1._rk + st/ks + ft/kf
      rc = 1._rk/c
      a = x3 + k1p*x2 + k12p*x + k123p
      ra = 1._rk/a
      da = 3._rk*x2 + 2._rk*k1p*x + k12p
      b = x2 + k1*x + k12
      rb = 1._rk/b
      db = 2._rk*x + k1
      rkb = 1._rk/kb
      rksi = 1._rk/ksi
      kwrx = kw*rx
      ptra = pt*ra
      xrc = x*rc
      crx = c*rx
      crx2 = crx*rx
      t1 = k1*dic*rb
      t2 = 2._rk*dic*k12*rb
      t3 = 1._rk + x*rkb
      t4 = 1._rk + x*rksi 
      t5 = 1._rk + ks*crx
      t6 = 1._rk + kf*crx

      !-------------------------------------------------------------------------
      ! fn = hco3+co3+borate+oh+hpo4+2*po4+silicate-hfree-hso4-hf-h3po4-ta
      ! df = dfn/dx

      fn = t1*x + t2 + bt/t3 + kwrx + ptra*k12p*x + 2._rk*ptra*k123p &
          + sit/t4 - xrc - st/t5 - ft/t6 - ptra*x3 - ta

      df = t1 - t1*x*db*rb - t2*db*rb - bt*rkb/t3**2._rk - kwrx*rx &
          + (ptra*k12p*(a - x*da))*ra - 2._rk*ptra*k123p*da*ra &
          - sit*rksi/t4**2._rk - rc - st*t5**(-2._rk)*(ks*crx2) &
          - ft*t6**(-2._rk)*(kf*crx2) - ptra*x2*(3._rk*a - x*da)*ra
      return
   end


   subroutine calc_hSWS_approx (dic, ta, pt, sit, bt, k1, k2, &
       k1p, k2p, k3p, kb, kw, ksi, H)
      !-------------------------------------------------------------------------
      ! solve carbonate system for H+
      ! based on code from M. Follows, T. Ito, S. Dutkiewicz
      ! see: Ocean Modelling 12 (2006), 290-301
      !-------------------------------------------------------------------------

      implicit none

      real(rk) dic, ta, pt, sit, bt, k1, k2, k1p, k2p, k3p, kb, kw, ksi
      real(rk) H, gamm, hg, cag, bohg, h3po4g, h2po4g, hpo4g, po4g
      real(rk) siooh3g, recip, dummy

      ! dic           = dissolved inorganic carbon
      ! ta            = total alkalinity
      ! pt            = dissolved inorganic phosphorus 
      ! sit           = dissolved inorganic silica
      ! bt            = dissolved inorganic boron
      ! ta            = total alkalinity
      ! k1, k2        = carbonate equilibrium coefficients
      ! kw            = dissociation of water
      ! klp, k2p, k3p = phosphate equilibrium coefficients
      ! ksi, kb       = silicate and borate equilibrium coefficients
      ! H             = [H+]
      
      ! First guess of [H+]: from last timestep
      hg = H
      ! estimate contributions to total alk from borate, silicate, phosphate
      bohg = bt*kb/(hg + kb)
      siooh3g = sit*ksi/(ksi + hg)
      recip = 1._rk/(hg*hg*hg + k1p*hg*hg + k1p*k2p*hg + k1p*k2p*k3p)
      h3po4g = (pt*hg*hg*hg)*recip
      h2po4g = (pt*k1p*hg*hg)*recip
      hpo4g = (pt*k1p*k2p*hg)*recip
      po4g = (pt*k1p*k2p*k3p)*recip
      ! estimate carbonate alkalinity
      cag = ta -bohg -(kw/hg) +hg -hpo4g -(2._rk*po4g) +h3po4g -siooh3g
      ! improved estimate of hydrogen ion concentration
      gamm = dic/cag
      dummy = (1._rk - gamm)*(1._rk - gamm)*k1*k1 - 4._rk*k1*k2*(1._rk - 2._rk*gamm)
      H = 0.5_rk*((gamm - 1._rk)*k1 + sqrt(dummy))
      return
   end   
end module
    