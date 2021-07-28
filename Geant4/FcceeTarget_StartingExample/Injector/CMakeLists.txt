# Created by hand by Yongke ZHAO on 17/10/2018 #

# Setup the project
cmake_minimum_required(VERSION 2.6)
project(Injector)
set(CMAKE_CXX_STANDARD 14) # Enable c++14 standard

if(MSVC)
  #add_compile_options(/W4 /WX)
else()
  #add_compile_options(-Wall -Wextra -pedantic -Werror)
  add_compile_options(-Wno-unused)
endif()

# Find Geant4 package, activating all available UI and Vis drivers by default
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

# Find ROOT (required package)
find_package(ROOT REQUIRED)

# Setup Geant4 include directories and compile definitions
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/Fot ${ROOT_INCLUDE_DIRS})

# Locate sources and headers for this project
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc ${PROJECT_SOURCE_DIR}/Fot/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/Fot/include/*.hh)

# Add the executable, and link it to the Geant4 libraries
add_executable(injector injector.cpp ${sources} ${headers})
target_link_libraries(injector ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

# For internal Geant4 use - but has no effect if you build it standalone
add_custom_target(Injector DEPENDS injector)

# Copy files to the build directory
#set(injector_copyfiles
#  mac
#  )
#foreach(_copyfile ${injector_copyfiles})
#  configure_file(
#    ${PROJECT_SOURCE_DIR}/${_copyfile}
#    ${PROJECT_BINARY_DIR}/${_copyfile}
#    COPYONLY
#    )
#endforeach()

# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#install(TARGETS injector DESTINATION bin)
