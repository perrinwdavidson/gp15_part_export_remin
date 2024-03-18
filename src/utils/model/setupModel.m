%%  setup environment
% cd('/net/ocean-dev.mit.edu/export/home/users/perrinwd/Research/whoi/gp15');
cd('/Users/perrindavidson/research/whoi/current/gp15'); 
addpath(genpath(pwd));
clear; 
close('all');
warning('off')

%%  set in and out paths
input_basepath = [pwd '/data/data_raw/'];
pro_output_basepath = [pwd '/data/data_pro/'];
sim_output_basepath = [pwd '/data/sim/'];
plot_output_basepath = [pwd '/plots/'];

%%  end subroutine
