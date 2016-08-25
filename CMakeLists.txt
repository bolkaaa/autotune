cmake_minimum_required (VERSION 2.8)

project (AutoTune)

include_directories(/home/matus/Downloads/supercollider/include/plugin_interface)
include_directories(/home/matus/Downloads/supercollider/include/common)
include_directories(/home/matus/Downloads/supercollider/external_libraries/libsndfile/)

set(CMAKE_SHARED_MODULE_PREFIX "") 
set(CMAKE_CXX_FLAGS	"-std=c++11 -g -Wall -Wextra -Wconversion")

if(APPLE OR WIN32)
	set(CMAKE_SHARED_MODULE_SUFFIX ".scx")
endif()

link_directories(/usr/lib/x86_64-linux-gnu/)
add_library(AutoTune MODULE AutoTune.cpp mayer_fft.c mayer_fft.h)
