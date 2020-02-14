
%% Stable MPC short preview
clear;clc;close;figure;

foot_distance_x = 0.20;
foot_distance_y = 0.18;
S = 30;
D = 20;
N = 2*(S+D);
omega = sqrt(9.8/0.8);
fs_matrix = [0,-foot_distance_y;%middle of rhe foot point
             foot_distance_x,foot_distance_y;
             2*foot_distance_x,-foot_distance_y;
             3*foot_distance_x,foot_distance_y;
             4*foot_distance_x,-foot_distance_y;
             5*foot_distance_x,foot_distance_y;
             6*foot_distance_x,-foot_distance_y;
             7*foot_distance_x,foot_distance_y;
             8*foot_distance_x,-foot_distance_y;
             9*foot_distance_x,+foot_distance_y;
             10*foot_distance_x,-foot_distance_y;
             11*foot_distance_x,+foot_distance_y;
             12*foot_distance_x,-foot_distance_y];
run('general_graphs.m')
