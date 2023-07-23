# # cmake_minimum_required(VERSION 3.0)

# cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
# project(matplotlib_cpp LANGUAGES CXX)
# include_directories(wamv14/inc)
# include_directories(modeltest)
# include(GNUInstallDirs)
# set(PACKAGE_NAME matplotlib_cpp)
# set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/${PACKAGE_NAME}/cmake)


# # Library target
# add_library(matplotlib_cpp INTERFACE)
# target_include_directories(matplotlib_cpp
#   INTERFACE
#     $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/examples>
#     $<INSTALL_INTERFACE:include>
# )
# target_compile_features(matplotlib_cpp INTERFACE
#   cxx_std_11
# )
# # TODO: Use `Development.Embed` component when requiring cmake >= 3.18
# find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
# target_link_libraries(matplotlib_cpp INTERFACE
#   Python3::Python
#   Python3::Module
# )
# find_package(Python3 COMPONENTS NumPy)
# if(Python3_NumPy_FOUND)
#   target_link_libraries(matplotlib_cpp INTERFACE
#     Python3::NumPy
#   )
# else()
#   target_compile_definitions(matplotlib_cpp INTERFACE WITHOUT_NUMPY)
# endif()
# install(
#   TARGETS matplotlib_cpp
#   EXPORT install_targets
# )


# # Examples
# add_executable(test modeltest/main1.cpp wamv14/src/asv.cpp)
# target_link_libraries(test PRIVATE matplotlib_cpp)
# # set_target_properties(test PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")


# 声明要求的 cmake 最低版本
cmake_minimum_required(VERSION 3.0)
set(CMAKE_CXX_STANDARD 11)

project(TEST)

# 调试相关设置
set(CMAKE_BUILD_TYPE DEBUG)

# 包含项目需要的所有头文件目录
include_directories(commonfun/inc)
include_directories(wamv14/inc)
include_directories(nmpc/inc)
include_directories(ttg/inc)

# 查找python库版本
find_package(PythonLibs 3.8)


# 指定头文件路径
set(PYTHON3.8_INLCUDE_DIRS "/usr/include/python3.8")

# 添加头文件到工程
include_directories(
    ${PYTHON3.8_INLCUDE_DIRS}
    )
#这样配置可以使得matplotlib.h中找到python.h文件，通常这样设置就可以。

# 添加一个可执行程序
add_executable(test modeltest/main1.cpp wamv14/src/asv.cpp commonfun/src/comfun.cpp 
               nmpc/src/nmpc.cpp nmpc/src/nmpc1.cpp ttg/src/single_trajectory.cpp)

# 添加相关库文件链接到工程
target_include_directories(test PRIVATE ${PYTHON2_INCLUDE_DIRS})
target_link_libraries(test ${PYTHON_LIBRARIES})
target_link_libraries(test /usr/local/lib/libcasadi.so.3.7) 
# 设置输出的可执行文件存放目录
set_target_properties(test PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
