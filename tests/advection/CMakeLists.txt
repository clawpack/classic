include_directories(${PROJECT_BINARY_DIR}/src)

macro(ADDTEST name)
    add_executable(${name} ${name}.f90 rp1_advection.f90 rp1_acoustics.f90)
    target_link_libraries(${name} classic ${LIBS})
    add_test(${name} ${PROJECT_BINARY_DIR}/${name})
endmacro(ADDTEST)

project(tests)

ADDTEST(main)
