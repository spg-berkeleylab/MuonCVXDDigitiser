########################################################
# cmake file for building MuonCVXDDigitiser Package
CMAKE_MINIMUM_REQUIRED(VERSION 3.11 FATAL_ERROR)
########################################################

# project name
PROJECT( MuonCVXDDigitiser )


# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 1 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

### BUILD OPTIONS ###########################################################

set(PACKAGE_INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(PACKAGE_INSTALL_INC_DIR include CACHE PATH "Installation directory for headers")
option(ZSEGMENTED "Sensitives are segmented along z")

### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )

# load default settings from ILCSOFT_CMAKE_MODULES
#INCLUDE( ilcsoft_default_settings )

FIND_PACKAGE( Marlin 1.0 REQUIRED ) # minimum required Marlin version
INCLUDE_DIRECTORIES( ${Marlin_INCLUDE_DIRS} )
LINK_LIBRARIES( ${Marlin_LIBRARIES} )
ADD_DEFINITIONS( ${Marlin_DEFINITIONS} )

FIND_PACKAGE( ROOT REQUIRED )
INCLUDE_DIRECTORIES( ${ROOT_INCLUDE_DIRS} )
LINK_LIBRARIES( ${ROOT_LIBRARIES} )

FIND_PACKAGE( LCIO REQUIRED )
INCLUDE_DIRECTORIES( ${LCIO_INCLUDE_DIRS} )
LINK_LIBRARIES( ${LCIO_LIBRARIES} )

FIND_PACKAGE( GSL REQUIRED )
INCLUDE_DIRECTORIES( ${GSL_INCLUDE_DIRS} )
LINK_LIBRARIES( ${GSL_LIBRARIES} )


find_package( DD4hep REQUIRED COMPONENTS DDRec)
INCLUDE_DIRECTORIES( SYSTEM ${DD4hep_INCLUDE_DIRS} )
LINK_LIBRARIES( ${DD4hep_LIBRARIES} ${DD4hep_COMPONENT_LIBRARIES}  )

### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ENDIF()

### LIBRARY #################################################################

# definitions to pass to the compiler
ADD_DEFINITIONS( "-Wall -ansi -pedantic -std=gnu++11" )
ADD_DEFINITIONS( "-Wno-long-long -Wno-effc++ -Wno-unused-parameter" )
if( ZSEGMENTED )
    message ( STATUS "Enabled sensor segmentation" )
    add_compile_definitions( "ZSEGMENTED" )
endif()

# include directories
INCLUDE_DIRECTORIES( ./include )
#INSTALL_DIRECTORY( ./include
#                   DESTINATION ${PACKAGE_INSTALL_INC_DIR}
#                   FILES_MATCHING PATTERN "*.h" )

# add library
FILE(GLOB_RECURSE library_sources "src/*.cc")
ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources} )
INSTALL_SHARED_LIBRARY( ${PROJECT_NAME} DESTINATION ${PACKAGE_INSTALL_LIB_DIR} )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

