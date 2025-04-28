#! /bin/bash

#apptainer exec --cleanenv --env USER=${USER} --env DISPLAY=${DISPLAY} --env SSH_AUTH_SOCK=${SSH_AUTH_SOCK} /mnt/rclsfserv005/local/dipole_fitting.sif /bin/bash --norc --noprofile
apptainer exec --cleanenv --env USER=${USER} --env DISPLAY=${DISPLAY} --env SSH_AUTH_SOCK=${SSH_AUTH_SOCK} --env AGREE_TO_MATLAB_RUNTIME_LICENSE=yes /mnt/rclsfserv005/local/dipole_fitting.sif /bin/bash --norc --noprofile
