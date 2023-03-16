clc;
clear all;

% load data
% load 50g_1rec_scale-free_x_1_1000measurements.mat;
% data=x;
load Y200_new.mat;
data=Y;
 

   n_time_points = size(Y,1); 
   kfold=5; %must greater than one;
   maxCenter=5;
  
   
   [z, C]= estimate_numC_NCGC( data, kfold, maxCenter ,n_time_points);