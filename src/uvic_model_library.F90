module uvic_model_library

   use fabm_types, only: type_base_model_factory, type_base_model
   
   use uvic_detritus
   use uvic_light
   use uvic_nut_chem
   use uvic_phytoplankton
   use uvic_solar
   use uvic_zooplankton
   use uvic_sediment

   implicit none

   private

   type, extends(type_base_model_factory) :: type_factory
   contains
      procedure :: create
   end type

   type (type_factory), save, target, public :: uvic_model_factory

contains

   subroutine create(self,name,model)
      class (type_factory), intent(in) :: self
      character(*),         intent(in) :: name
      class (type_base_model), pointer :: model

      select case (name)
         case ('detritus');          allocate(type_uvic_detritus::model)
         case ('light');             allocate(type_uvic_light::model)
         case ('nut_chem');          allocate(type_uvic_nut_chem::model)
         case ('phytoplankton');     allocate(type_uvic_phytoplankton::model)
         case ('solar');             allocate(type_uvic_solar::model)
         case ('zooplankton');       allocate(type_uvic_zooplankton::model) 
         case ('sediment');          allocate(type_uvic_sediment::model)
         case default
            call self%type_base_model_factory%create(name, model)
      end select
   end subroutine create
end module
