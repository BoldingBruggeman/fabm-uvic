module uvic_model_library

   use fabm_types, only: type_base_model_factory, type_base_model
   
   use uvic_tracer
   use uvic_phytoplankton
   use uvic_zooplankton

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
         case ('tracer');            allocate(type_uvic_tracer::model)
         case ('phytoplankton');     allocate(type_uvic_phytoplankton::model) 
         case ('zooplankton');       allocate(type_uvic_zooplankton::model) 
         case default
            call self%type_base_model_factory%create(name, model)
      end select
   end subroutine create
end module
