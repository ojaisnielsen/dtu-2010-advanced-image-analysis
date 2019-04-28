function [] = initvlfeat()

p = pwd;
cd 'vlfeat/toolbox';
vl_setup;
cd(p);
clear p

end