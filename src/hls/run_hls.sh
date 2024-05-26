
#!/bin/bash

PLATFORM=xilinx_u250_gen3x16_xdma_4_1_202210_1


# Execute commands for hostname not equal to str1
echo "you are on nsc"
. /ceph/grid/software/xilinx/tools/Vitis/2023.1/settings64.sh
. /opt/xilinx/xrt/setup.sh
export PLATFORM_REPO_PATHS=/opt/xilinx/platforms/${PLATFORM}/
module load foss/2022a



HOME=$PWD
DIR=../src
BUILD_DIR=./build
mkdir -p $BUILD_DIR

KERNEL=dwt

cp ./src/xrt.ini $BUILD_DIR/xrt.ini
cp ./src/run_hls_dwt.tcl $BUILD_DIR/run_hls_dwt.tcl
cp ./src/run_hls_idwt.tcl $BUILD_DIR/run_hls_idwt.tcl




cd $BUILD_DIR


# uncomment to run HLS for dwt
#vitis_hls -f run_hls_dwt.tcl


# uncomment to run HLS for idwt
vitis_hls -f run_hls_idwt.tcl

cd $HOME