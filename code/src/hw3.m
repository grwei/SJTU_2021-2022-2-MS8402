%% hw3.m
% Description: MATLAB code for homework 3 (MS8402, 2022 Spring)
% Author: Guorui Wei (危国锐) (313017602@qq.com; weiguorui@sjtu.edu.cn)
% Student ID: 120034910021
% Created: 2022-04-01
% Last modified: 2022-04-

%% Initialize project

clc; clear; close all
init_env();

%% Question 1

epsilon = .1;
tau_0 = 1;
x = linspace(0,1,127);
y = x;
[X,Y] = meshgrid(x,y);
psi = pi*tau_0*(1-X-exp(-X/epsilon)).*sin(pi*Y);

% fig.
figure("Name","Question 1")
t_TCL = tiledlayout(1,1,"TileSpacing","tight","Padding","tight");
t_Axes = nexttile(t_TCL,1);
[~,c_contour] = contour(t_Axes,X,Y,psi,"Fill","off");
clabel([],c_contour,c_contour.TextList([1,idivide(end+1,uint8(2)),end]),"Interpreter",'latex')
c = colorbar(t_Axes,"Location","eastoutside","TickLabelInterpreter","latex");
c.Label.String = "\psi";
set(t_Axes,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off',"TickDir","out");
ylabel(t_Axes,"$y$","Interpreter",'latex');
xlabel(t_Axes,"$x$","Interpreter",'latex');
title(t_Axes,sprintf("$\\psi = (1-x-\\mathrm{e}^{-x/\\varepsilon}) \\pi \\tau_0 \\sin{\\pi y},\\quad \\tau_0 = %.3g$, $\\varepsilon = %.3g.$",tau_0,epsilon),"Interpreter",'latex')
[t_title_t,t_title_s] = title(t_TCL,"\bf 2022 Spring MS8402 Homework 3 Q1","Guorui Wei 120034910021","Interpreter",'latex');
set(t_title_s,'FontSize',8)
exportgraphics(t_TCL,"..\\doc\\fig\\hw3_Q1.emf",'Resolution',600,'ContentType','auto','BackgroundColor','none','Colorspace','rgb')
exportgraphics(t_TCL,"..\\doc\\fig\\hw3_Q1.png",'Resolution',600,'ContentType','auto','BackgroundColor','none','Colorspace','rgb')

%% Question 2

clear; clc; close all;
path_votemper = "..\data\votemper_control_monthly_highres_3D_202201_OPER_v0.1.nc";
finfo = ncinfo(path_votemper);
nav_lat = ncread(path_votemper,'nav_lat');
nav_lon = ncread(path_votemper,'nav_lon');
deptht = ncread(path_votemper,'deptht');
lat = nav_lat(1,:);
lon = nav_lon(:,1);

% zonal average at sea surface
sosstsst = ncread("..\data\sosstsst_control_monthly_highres_2D_202201_OPER_v0.1.nc",'sosstsst'); % [deg C]
sosaline = ncread("..\data\sosaline_control_monthly_highres_2D_202201_OPER_v0.1.nc",'sosaline'); % [psu]
sosstsst_zonal_avg = mean(sosstsst,1,"omitnan"); % [deg C]
sosaline_zonal_avg = mean(sosstsst,1,"omitnan"); % [psu]

%
figure("Name","Q2_sst_zonal_avg")
t_TCL_1 = tiledlayout(1,1,"TileSpacing","tight","Padding","tight");
t_Axes_1 = nexttile(t_TCL_1,1);
t_plot_1 = plot(t_Axes_1,lat,sosstsst_zonal_avg,'-',"DisplayName",'SST');
set(t_Axes_1,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off');
xlabel(t_Axes_1,"latitude $(^{\circ}{\rm{N}})$","Interpreter",'latex');
ylabel(t_Axes_1,"sea surface temperature","Interpreter",'latex');
title(t_Axes_1,"sea surface temperature zonal average",'Interpreter','latex')
%
legend(t_plot_1,"Location",'east','Interpreter','latex',"Box","off",'FontSize',10);
[t_title_t,t_title_s] = title(t_TCL_1,"\bf 2022 Spring MS8402 Hw1 Q3","Guorui Wei 120034910021","Interpreter",'latex');
set(t_title_s,'FontSize',8)
%
exportgraphics(t_TCL_1,"..\\doc\\fig\\hw3_Q2_sst_zonal_avg.emf",'Resolution',600,'ContentType','auto','BackgroundColor','none','Colorspace','rgb')
exportgraphics(t_TCL_1,"..\\doc\\fig\\hw3_Q2_sst_zonal_avg.png",'Resolution',600,'ContentType','auto','BackgroundColor','none','Colorspace','rgb')

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
