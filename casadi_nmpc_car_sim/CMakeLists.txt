# 声明要求的 cmake 最低版本
cmake_minimum_required(VERSION 3.0)
set(CMAKE_CXX_STANDARD 11)

project(TEST)

# 调试相关设置
set(CMAKE_BUILD_TYPE DEBUG)

# 包含项目需要的所有头文件目录

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
add_executable(test test.cpp comfun.cpp car_model.cpp nmpc1.cpp)

# 添加相关库文件链接到工程
target_include_directories(test PRIVATE ${PYTHON2_INCLUDE_DIRS})
target_link_libraries(test ${PYTHON_LIBRARIES})
target_link_libraries(test /usr/local/lib/libcasadi.so.3.7) 
# 设置输出的可执行文件存放目录
set_target_properties(test PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
