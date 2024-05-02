#!/bin/bash
echo "Hello, World!"

module load comsol matlab
comsol server &
matlab -batch "addpath('/apps/generic/comsol/6.2/mli/');pause(10);disp('OK1');mphstart(2036);disp('OK2');import com.comsol.model.util.*;modelname='MPH/Scintillator3D_1DStudy_2Dgeomv2 - Copy';mphopen(modelname);exit;" 

