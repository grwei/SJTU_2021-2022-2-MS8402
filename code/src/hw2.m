%% hw2.m
% Description: MATLAB code for homework 2 (MS8402, 2022 Spring)
% Author: Guorui Wei (危国锐) (313017602@qq.com; weiguorui@sjtu.edu.cn)
% Student ID: 120034910021
% Created: 2022-03-09
% Last modified: 2022-03-

%% Initialize project

clc; clear; close all
init_env();

%% Question 2

U_0 = 2;    % [m/s]
L_0 = 300;  % [m]
H_0 = 200;  % [m]
H_M = 300;  % [m]
H_S = 100;  % [m]
H_N = 500;  % [m]
f = 1e-4;   % [1/s]

V_M = sqrt(2*U_0*L_0*(1-H_0/H_M)*f);
U_S = sqrt(2*U_0*L_0*(H_0/H_S-1)*f);
V_N = sqrt(2*U_0*L_0*(1-H_0/H_N)*f);
L_M = sqrt(2*U_0*L_0*H_0^2/(H_M*(H_M-H_0)*f));
L_S = sqrt(2*U_0*L_0*H_0^2/(H_S*(H_0-H_S)*f));
L_N = sqrt(2*U_0*L_0*H_0^2/(H_N*(H_N-H_0)*f));

fprintf("V_M = %.2f (m/s),\t L_M = %.0f (m),\n" + ...
    "U_S = %.2f (m/s),\t L_S = %.0f (m),\n" + ...
    "V_N = %.2f (m/s),\t L_N = %.0f (m).\n", ...
    V_M,L_M,U_S,L_S,V_N,L_N);

%% local functions

%% Initialize environment
function [] = init_env()
    % set up project directory
    if ~isfolder("../doc/fig/")
        mkdir ../doc/fig/
    end
    % configure searching path
    mfile_fullpath = mfilename('fullpath'); % the full path and name of the file in which the call occurs, not including the filename extension.
    mfile_fullpath_without_fname = mfile_fullpath(1:end-strlength(mfilename));
    addpath(genpath(mfile_fullpath_without_fname + "../data"), ...
            genpath(mfile_fullpath_without_fname + "../inc")); % adds the specified folders to the top of the search path for the current MATLAB® session.
end
