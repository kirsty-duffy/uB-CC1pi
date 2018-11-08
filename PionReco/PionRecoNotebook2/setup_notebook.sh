#!/bin/bash

setup_uboone
setup uboonecode v07_07_00 -qe17:prof
export PYTHONUSERBASE=/uboone/app/users/kduffy/python_libs
export PYTHONPATH=$PYTHONUSERBASE/bin:$PYTHONPATH
export PATH=$PYTHONUSERBASE/bin:$PATH

root --notebook