cmake_minimum_required(VERSION 3.16)

if (NOT TARGET cage-core)
	message(FATAL_ERROR "include cage first")
endif()

add_subdirectory(sources)
