set TOP tau_mass

# Set the top level module
set_top ${TOP}

# Add source code
add_files ${PROJ_DIR}/src/${TOP}.cpp

# Add testbed files
add_files -tb ${PROJ_DIR}/src/${TOP}_tb.cpp -cflags ${CFLAGS}

add_files -tb ${PROJ_DIR}/src/3prong.csv
