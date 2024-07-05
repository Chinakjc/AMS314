#!/bin/bash

# 默认的mesh文件名
MESH_NAME="carre_4h.mesh"

# 检查是否有命令行参数
if [ "$#" -eq 2 ]; then
    if [ "$1" == "--msh" ]; then
        MESH_NAME=$2
    fi
fi

# 编译和运行程序
cd src/
make -f makemesh.make
cd ..
./src/mesh data/${MESH_NAME}
