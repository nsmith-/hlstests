## Set the top level module
set_top hls_demo
#
### Add source code
add_files ${PROJ_DIR}/src/hls_demo.cpp

## Add testbed files
add_files -tb ${PROJ_DIR}/src/hls_demo_tb.cpp -cflags ${CFLAGS}
