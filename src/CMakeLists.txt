add_library(fabm_models_uvic OBJECT
            uvic_model_library.F90
			tracer.F90
			phytoplankton.F90
            zooplankton.F90
           )

add_dependencies(fabm_models_uvic fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2024 - Bolding & Bruggeman ApS
#-----------------------------------------------------------------------
