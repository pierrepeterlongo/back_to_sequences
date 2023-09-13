#!/bin/bash

# get the full path of the current directory
install_bin_path=`pwd`

# install kmindex
cmd="cd ./thirdparty/kmindex/ && sh install.sh -p ${install_bin_path}"
echo $cmd
eval $cmd