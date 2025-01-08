########################################################
# cmake file for building MuonCVXDDigitiser Package
CMAKE_MINIMUM_REQUIRED(VERSION 3.11 FATAL_ERROR)
########################################################

# project name
PROJECT( MuonCVXDDigitiser )


# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 2 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

### BUILD OPTIONS ###########################################################

set(PACKAGE_INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(PACKAGE_INSTALL_INC_DIR include CACHE PATH "Installation directory for headers")
option(ZSEGMENTED "Sensitives are segmented along z" ON)
option(OPENMP_ENABLED "Activate OpenMP parallel compilation" ON)
if(OPENMP_FOUND AND OPENMP_ENABLED)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
endif()

# definitions to pass to the compiler
ADD_DEFINITIONS( "-Wno-long-long -Wno-effc++ -Wno-unused-parameter -Wno-narrowing" )
if( ZSEGMENTED )
    message ( STATUS "Enabled sensor segmentation" )
    add_compile_definitions( "ZSEGMENTED" )
endif()

### LIBRARY #################################################################

FILE ( GLOB _plugin_sources src/*.c* )

gaudi_add_module(MuonCVXDigitiserPlugins
                 SOURCES	${_plugin_sources}
                 LINK   	Gaudi::GaudiKernel
                        	EDM4HEP::edm4hep
                        	ROOT::Hist
                        	${ROOT_LIBRARIES}
                        	DD4hep::DDCore
                        	k4FWCore::k4FWCore
                        	${GSL_LIBRARIES}
				${CLHEP_LIBRARIES})

INSTALL( TARGETS MuonCVXDigitiserPlugins
         EXPORT ${CMAKE_PROJECT_NAME}TARGETS
         RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
         LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
         PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/@{CMAKE_PROJECT_NAME}"
         COMPONENT dev )


# include directories
TARGET_INCLUDE_DIRECTORIES(MuonCVXDigitiserPlugins 
                   PRIVATE ./include 
                           ${CMAKE_CURRENT_BINARY_DIR}
                           ${Gaudi_INCLUDE_DIRS}
                           ${ROOT_INCLUDE_DIRS}
			   ${GSL_INCLUDE_DIRS})

SET( GAUDI_GENCONF_DIR "genConfDir" )
if (${Gaudi_VERSION} VERSION_LESS 35.1)
	SET( GAUDI_GENCONF_DIR "genConf" )
endif()

