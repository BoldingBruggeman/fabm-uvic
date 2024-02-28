module ihamocc_model_library

   use fabm_types, only: type_base_model_factory, type_base_model

   use ihamocc_oxygen
   use ihamocc_carbon
   use ihamocc_bromo
   use ihamocc_cfc
   use ihamocc_dms
   use ihamocc_alkalinization
   use ihamocc_iron
   use ihamocc_nitrogen
   use ihamocc_preformed_tracer
   use ihamocc_phytoplankton
   use ihamocc_zooplankton
   use ihamocc_detritus
   use ihamocc_cisonew
   use ihamocc_natdic
   use ihamocc_mixed_layer
   use ihamocc_sediment_bypass
   use ihamocc_light
   use ihamocc_tracer
   implicit none

   private

   type, extends(type_base_model_factory) :: type_factory
   contains
      procedure :: create
   end type

   type (type_factory), save, target, public :: ihamocc_model_factory

contains

   subroutine create(self,name,model)
      class (type_factory), intent(in) :: self
      character(*),         intent(in) :: name
      class (type_base_model), pointer :: model

      select case (name)
         case ('oxygen');            allocate(type_ihamocc_oxygen::model)
         case ('carbon');            allocate(type_ihamocc_carbon::model)
         case ('bromo');             allocate(type_ihamocc_bromo::model)
         case ('cfc');               allocate(type_ihamocc_cfc::model)
         case ('dms');               allocate(type_ihamocc_dms::model)
         case ('alkalinization');    allocate(type_ihamocc_alkalinization::model)    
         case ('iron');              allocate(type_ihamocc_iron::model)               
         case ('nitrogen');          allocate(type_ihamocc_nitrogen::model) 
         case ('preformed_tracer');  allocate(type_ihamocc_preformed_tracer::model)
         case ('phytoplankton');     allocate(type_ihamocc_phytoplankton::model)
         case ('zooplankton');       allocate(type_ihamocc_zooplankton::model)
         case ('detritus');          allocate(type_ihamocc_detritus::model)    
         case ('cisonew');           allocate(type_ihamocc_cisonew::model)    
         case ('natdic');            allocate(type_ihamocc_natdic::model)    
         case ('mixed_layer');       allocate(type_ihamocc_mixed_layer::model)
         case ('sediment_bypass');   allocate(type_ihamocc_sediment_bypass::model)
         case ('light');             allocate(type_ihamocc_light::model)   
         case ('tracer');            allocate(type_ihamocc_tracer::model)
         case default
            call self%type_base_model_factory%create(name, model)
      end select
   end subroutine create
end module
