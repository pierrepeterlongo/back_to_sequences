#!/bin/bash

# get the full path of the current directory
install_bin_path=`pwd`

cd ./thirdparty/kmindex/ && sh install.sh -p ${install_bin_path}