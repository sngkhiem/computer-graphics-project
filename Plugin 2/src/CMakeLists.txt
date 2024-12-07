# Copyright 2019-2020, Collabora, Ltd.
# SPDX-License-Identifier: BSL-1.0


# Only build if we have qhull
if(TARGET external-qhull)

	set(SOURCES filter_idt.cpp)

	set(HEADERS filter_idt.h)

	add_meshlab_plugin(filter_idt ${SOURCES} ${HEADERS})

	target_link_libraries(filter_idt PRIVATE external-qhull)

else()
	message(STATUS "Skipping filter_qhull - missing qhull")
endif()