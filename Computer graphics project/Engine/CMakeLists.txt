set(SOURCES
        ../lib/Exceptions/math_exceptions.cpp)

set(HEADERS
        ../lib/Math/base.h
        ../lib/Exceptions/engine_exceptions.h ../lib/Engine/Entities.h)

add_library(${CMAKE_PROJECT_NAME}_lib STATIC ${SOURCES} ${HEADERS})