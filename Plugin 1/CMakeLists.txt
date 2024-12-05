# Copyright 2019-2020, Collabora, Ltd.
# SPDX-License-Identifier: BSL-1.0


# Only build if we have qhull
if(TARGET external-qhull)

	set(SOURCES filter_normal.cpp)

	set(HEADERS filter_normal.h)

	add_meshlab_plugin(filter_normal ${SOURCES} ${HEADERS})

	target_link_libraries(filter_normal PRIVATE external-qhull)

else()
	message(STATUS "Skipping filter_qhull - missing qhull")
endif()