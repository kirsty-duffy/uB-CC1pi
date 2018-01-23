
art_make(	
	  LIB_LIBRARIES
			${ART_FRAMEWORK_CORE}
			art_Utilities canvas
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			cetlib cetlib_except
                   	larpandora_LArPandoraInterface	
                        ${Boost_SYSTEM_LIBRARY}
			${ROOT_BASIC_LIB_LIST}
          MODULE_LIBRARIES 
		   	larpandora_LArPandoraInterface	
			larcorealg_Geometry
			larcore_Geometry_Geometry_service
			lardataobj_RecoBase
			lardataobj_AnalysisBase	
			lardata_RecoObjects larreco_Deprecated  
			larreco_RecoAlg
			lardata_Utilities
			larevt_Filters
			nusimdata_SimulationBase
			${LARRECO_LIB}
			${LARDATA_LIB}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			art_Utilities canvas
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	                ${FHICLCPP}
			cetlib cetlib_except
			${CLHEP}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
			${ROOT_BASIC_LIB_LIST}
        )

install_headers()
install_fhicl()
install_source()

