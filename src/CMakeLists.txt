add_library(fabm_models_uvic OBJECT
            uvic_model_library.F90
			detritus.F90
			light.F90
            nut_chem.F90
			phytoplankton.F90
			shared.F90
			solar.F90
            zooplankton.F90
           )

add_dependencies(fabm_models_uvic fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2024 - Bolding & Bruggeman ApS
#-----------------------------------------------------------------------
