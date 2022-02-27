%% hw1.m
% Description: MATLAB code used in homework 1 (MS8402, 2022 Spring)
% Author: Guorui Wei (危国锐) (313017602@qq.com; weiguorui@sjtu.edu.cn)
% Student ID: 120034910021
% Created: 2022-02-27
% Last modified: 2022-02-28

%% Initialize project

clc; clear; close all
init_env();

%% Question 2

%
u = [-.25,-.25,.05,.03]; % [m/s]
v = [0,-.01,.25,.30]; % [m/s]
delta_x = 5e5; % [m]
delta_y = 5.4e5; % [m]
%
u_wrt_x = (u(2)-u(1)+u(3)-u(4))/(2*delta_x); % [1/s]
v_wrt_y = (v(1)-v(4)+v(2)-v(3))/(2*delta_y); % [1/s]
w = @(z) -(u_wrt_x+v_wrt_y)*z; % [m/s]
fprintf("w(-50) = %.2d cm/s\n",w(-50)*100);

%% Question 3

%
site_A.p = [0,20,40,60,80,100,200,250,500,1000].'; % [dbar]
site_A.SA = [35.10,34.99,34.88,34.78,34.68,34.60,34.45,34.35,34.25,34.53].'; % [g/kg]
site_A.CT = [28.50,28.45,28.35,24.55,22.75,20.55,15.50,13.00,6.58,4.20].'; % [deg C]
site_B.p = site_A.p; % [dbar]
site_B.SA = [33.50,33.50,34.25,34.55,34.65,34.74,34.90,35.10,35.23,35.40].'; % [g/kg]
site_B.CT = [2.50,3.74,4.02,4.10,4.15,4.20,4.30,4.35,4.25,3.75].'; % [deg C]
p_i = (site_A.p(1):1:site_A.p(end)).'; % specific query points at which the interpolated SA_i and CT_i are required [ dbar ]

%% Figure.

t_TCL = tiledlayout(1,2,"TileSpacing","tight","Padding","tight");
%
site_A = hw1_3(site_A,p_i,t_TCL,1,"\bf Site A");
site_B = hw1_3(site_B,p_i,t_TCL,2,"\bf Site B");
%
xlabel(t_TCL,"density (kg/$\rm{m}^3$)","Interpreter",'latex');
ylabel(t_TCL,"pressure (dbar)","Interpreter",'latex');
[t_title_t,t_title_s] = title(t_TCL,"\bf 2022 Spring MS8402 Hw1 Q3","Guorui Wei 120034910021","Interpreter",'latex');
set(t_title_s,'FontSize',8)
%
exportgraphics(t_TCL,"..\\doc\\fig\\hw1_Q3.emf",'Resolution',800,'ContentType','auto','BackgroundColor','none','Colorspace','rgb')
exportgraphics(t_TCL,"..\\doc\\fig\\hw1_Q3.png",'Resolution',800,'ContentType','auto','BackgroundColor','none','Colorspace','rgb')

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

%% Question 3

function [site_struct] = hw1_3(site_struct,p_i,t_TCL,num_Tile,textbox_string)
%% hw1_3
% Description.
arguments
    site_struct
    p_i
    t_TCL
    num_Tile = 1
    textbox_string = "\bf Site"
end

site_struct.p_i = p_i; % specific query points at which the interpolated SA_i and CT_i are required [ dbar ]
[site_struct.SA_i,site_struct.CT_i] = gsw_SA_CT_interp(site_struct.SA,site_struct.CT,site_struct.p,p_i); % SA and CT interpolation to p_i on a cast

%%% 1. Calculate density with the linear Equations of State (EOS).

rho_0 = 1027; % [kg/m^3]
beta_T = 0.15; % [kg/m^3/(deg C)]
beta_S = 0.78; % [kg/m^3/(g/kg)]
beta_p = 4.5;  % [kg/m^3/dbar]

func_rho_linear = @(T,S,p) rho_0 * (1 - beta_T*(T-T(1)) + beta_S*(S-S(1)) + beta_p*(p-p(1)));
site_struct.rho_linear = func_rho_linear(site_struct.CT, site_struct.SA, site_struct.p);
site_struct.rho_i_linear = func_rho_linear(site_struct.CT_i,site_struct.SA_i,site_struct.p_i);

%%% 2. Use TEOS-10 to compute the density.
% NOTE: Since the location of site A & B is not provided, I have to let SA, CT
% (TEOS-10) be SP, t, and therefore the results here should not be correct.

site_struct.rho_TEOS_10 = gsw_rho_CT_exact(site_struct.SA,site_struct.CT,site_struct.p); % Calculates in-situ density from Absolute Salinity and Conservative Temperature.
site_struct.rho_i_TEOS_10 = gsw_rho_CT_exact(site_struct.SA_i,site_struct.CT_i,site_struct.p_i);

%%% plot

%
t_Axes_TEOS = nexttile(t_TCL,num_Tile);
t_plot_1 = plot(t_Axes_TEOS,site_struct.rho_i_TEOS_10,site_struct.p_i,'-',"color",'#0072BD',"DisplayName",'TEOS-10');
hold on
t_plot_2 = plot(site_struct.rho_TEOS_10,site_struct.p,'.',"color",'#0072BD',"MarkerSize",10,"DisplayName",'');
set(t_Axes_TEOS,"YDir",'reverse',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off',"XColor",'#0072BD',"YColor",'#0072BD');
%
t_Axes_linear = axes(t_TCL);
t_Axes_linear.Layout.Tile = num_Tile;
t_plot_3 = plot(t_Axes_linear,site_struct.rho_i_linear,site_struct.p_i,'--','Color','#D95319',"DisplayName",'LEOS');
hold on
t_plot_4 = plot(site_struct.rho_linear,site_struct.p,'o',"MarkerSize",4,'Color','#D95319',"DisplayName",'');
set(t_Axes_linear,'YDir','reverse','FontSize',10,'TickLabelInterpreter','latex','XAxisLocation','top','YAxisLocation','right','YTickLabel',{},'Box','off','Color','none','XColor','#D95319','YColor','#D95319','YLimitMethod','tight')
%
linkaxes([t_Axes_TEOS,t_Axes_linear],'y');
legend([t_plot_2,t_plot_4],["TEOS-10","Linear EOS"],"Location",'southwest','Interpreter','latex',"Box","off");
annotation('textbox',[0.36+0.5*(num_Tile-1) 0.68 0.10 0.06],'String',textbox_string,'LineStyle','none','FontWeight','bold','Interpreter','latex');

end
