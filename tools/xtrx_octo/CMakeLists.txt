# Copyright (C) 2021 J. Arribas and C. Fernandez-Prades. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in
#    the documentation and/or other materials provided with the
#    distribution.
# 3. Neither the name GNSS-SDR-SIM nor the names of its contributors may be
#    used to endorse or promote products derived from this software
#    without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
# AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# This file is part of DONGLE-CAPTURE

find_package(libxtrxll REQUIRED)
find_package(LibXTRX REQUIRED)


set(test_xtrx_octo_src
	main.cc
    adf4355.c
)
	
add_executable(test_xtrx_octo ${test_xtrx_octo_src})

target_include_directories(test_xtrx_octo 
	PUBLIC 
		${LIBXTRX_INCLUDE_DIRS}
		${LIBXTRXLL_INCLUDE_DIRS})

target_link_libraries(test_xtrx_octo
    PUBLIC
        ${LIBXTRX_LIBRARIES}
        ${LIBXTRXLL_LIBRARIES}
        Gflags::gflags
        Glog::glog
        Threads::Threads
        ${Boost_LIBRARIES}
)

target_compile_definitions(test_xtrx_octo PRIVATE -DDONGLE_CAPTURE_VERSION="${VERSION}")

if(ENABLE_CLANG_TIDY)
    if(CLANG_TIDY_EXE)
        set_target_properties(test_xtrx_octo
            PROPERTIES
                CXX_CLANG_TIDY "${DO_CLANG_TIDY}"
        )
    endif()
endif()

add_custom_command(TARGET test_xtrx_octo
    POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:test_xtrx_octo>
        ${CMAKE_SOURCE_DIR}/install/$<TARGET_FILE_NAME:test_xtrx_octo>
)

install(TARGETS test_xtrx_octo
    RUNTIME DESTINATION bin
    COMPONENT "test_xtrx_octo"
)
