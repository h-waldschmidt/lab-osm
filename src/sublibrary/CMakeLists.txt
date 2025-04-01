set(HEADER_FILES sublib_component.h)
set(SOURCE_FILES sublib_component.cpp)

add_library(sublibrary SHARED ${HEADER_FILES} ${SOURCE_FILES})
target_link_libraries(sublibrary PUBLIC quill::quill)