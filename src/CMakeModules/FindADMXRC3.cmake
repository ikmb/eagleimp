
set(ADMXRC3_INC_PATHS /usr/include /opt/alphadata/driver/latest/include)
set(ADMXRC3_LIB_PATHS /usr/lib)

find_library(ADMXRC3_LIBRARY NAMES admxrc3)
find_path(ADMXRC3_INCLUDE_DIR admxrc3.h PATHS ${ADMXRC3_INC_PATHS})

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(ADMXRC3 DEFAULT_MSG ADMXRC3_LIBRARY ADMXRC3_INCLUDE_DIR)

