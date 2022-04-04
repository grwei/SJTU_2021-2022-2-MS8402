%% hw3.m
% Description: MATLAB code for homework 3 (MS8402, 2022 Spring)
% Author: Guorui Wei (危国锐) (313017602@qq.com; weiguorui@sjtu.edu.cn)
% Student ID: 120034910021
% Created: 2022-04-01
% Last modified: 2022-04-
% Data: [1] [ORAS5 global ocean reanalysis monthly data from 1958 to present](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-oras5)

%% Initialize project

clc; clear; close all
init_env();

%% Question 1

epsilon = .1;
tau_0 = 1;
N = 127;
x = linspace(0,1,N);
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
[~,t_title_s] = title(t_TCL,"\bf 2022 Spring MS8402 Homework 3 Q1","Guorui Wei 120034910021","Interpreter",'latex');
set(t_title_s,'FontSize',8)
exportgraphics(t_TCL,"..\\doc\\fig\\hw3_Q1.emf",'Resolution',600,'ContentType','auto','BackgroundColor','none','Colorspace','rgb')
exportgraphics(t_TCL,"..\\doc\\fig\\hw3_Q1.png",'Resolution',600,'ContentType','auto','BackgroundColor','none','Colorspace','rgb')

%% Question 2

% import data
clear; clc; close all;
path_votemper = "..\data\votemper_control_monthly_highres_3D_202201_OPER_v0.1.nc"; % The temperature of a parcel of sea water would have if moved adiabatically to sea level pressure. This variable is a 3D field.
path_vosaline = "..\data\vosaline_control_monthly_highres_3D_202201_OPER_v0.1.nc"; % The salt content of sea water as measured on the practical salinity units (PSU) scale. This variable is a 3D field.
finfo = ncinfo(path_votemper);
nav_lat = ncread(path_votemper,'nav_lat');
nav_lon = ncread(path_votemper,'nav_lon');
deptht = ncread(path_votemper,'deptht');
%
n_lat_bin = 255;
[lat_bin_num,lat_edges] = discretize(nav_lat,n_lat_bin);
Q2_zonal_average(nav_lat,lat_bin_num,lat_edges,'off');
%
lon = -75;
lat = 36;
[x,y,arclen] = Q2_location_query(lon,lat,nav_lon,nav_lat);

votemper = ncread(path_votemper,'votemper');
vosaline = ncread(path_vosaline,'vosaline');
for y = 495:1021
    for x = 1:1442
        [SP,pt,p,SA,CT,N2,p_mid] = Q2_gsw_N2(x,y,nav_lon,nav_lat,deptht,votemper,vosaline,path_votemper,path_vosaline,'on');
        fprintf("min(N2) = %.3g, x = %i, y = %i.\n",min(N2,[],'all'),x,y);
    end
end

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

%% surface zonal average 
function [] = Q2_zonal_average(nav_lat,lat_bin_num,lat_edges,fig_EN)
% hw3 Q2
%
    arguments
        nav_lat
        lat_bin_num
        lat_edges
        fig_EN string = 'on'
    end

    sosstsst = ncread("..\data\sosstsst_control_monthly_highres_2D_202201_OPER_v0.1.nc",'sosstsst'); % [deg C] Water temperature close to the ocean surface. This variable is a 2D field.
    sosaline = ncread("..\data\sosaline_control_monthly_highres_2D_202201_OPER_v0.1.nc",'sosaline'); % [psu] Salt concentration close to the ocean surface. This variable is a 2D field.
    %
    TF_so_avail = ~isnan(sosstsst);
    n_lat_bin = length(lat_edges)-1;
    lat_bin_avg = zeros(n_lat_bin,1);
    N_pts_bin = lat_bin_avg;
    sosstsst_zonal_avg = lat_bin_avg;
    sosaline_zonal_avg = lat_bin_avg;
    for bin_num = 1:n_lat_bin
        TF_lat_avail = (lat_bin_num == bin_num) & TF_so_avail;
        N_pts_bin(bin_num) = nnz(TF_lat_avail);
        lat_bin_avg(bin_num) = mean(nav_lat(TF_lat_avail));
        sosstsst_zonal_avg(bin_num) = mean(sosstsst(TF_lat_avail));
        sosaline_zonal_avg(bin_num) = mean(sosaline(TF_lat_avail));
    end
    if (strcmpi(fig_EN,"off"))
        return;
    end

    %
    figure("Name","Q2_sst_zonal_avg")
    t_TCL_1 = tiledlayout(1,1,"TileSpacing","tight","Padding","tight");
    t_Axes_1 = nexttile(t_TCL_1,1);
    yyaxis(t_Axes_1,"left")
    t_plot_sst = plot(t_Axes_1,lat_bin_avg,sosstsst_zonal_avg,'-',"DisplayName",'temperature');
    set(t_Axes_1,'YColor','#0072BD','XLim',[-90 90],"YDir",'normal','XTick',linspace(-90,90,7),"TickLabelInterpreter",'latex',"FontSize",10,'Box','off');
    xlabel(t_Axes_1,"latitude $(^{\circ}{\rm{N}})$","Interpreter",'latex');
    ylabel(t_Axes_1,"Sea Surface temperature $(^{\circ}{\rm{C}})$","Interpreter",'latex');
    yyaxis(t_Axes_1,"right")
    t_plot_sal = plot(t_Axes_1,lat_bin_avg,sosaline_zonal_avg,'-',"DisplayName",'Salinity');
    set(t_Axes_1,"YDir",'normal');
    ylabel(t_Axes_1,"Sea Surface Salinity (psu)","Interpreter",'latex');
    title(t_Axes_1,"\bf Zonal Average (Jan, 2022)",'Interpreter','latex')
    grid on
    %
    legend([t_plot_sst,t_plot_sal],"Location",'south','Interpreter','latex',"Box","off",'FontSize',10);
    [~,t_title_s] = title(t_TCL_1,"\bf 2022 Spring MS8402 Homework 3 Q2","Guorui Wei 120034910021","Interpreter",'latex');
    set(t_title_s,'FontSize',8)
    %
    exportgraphics(t_TCL_1,"..\\doc\\fig\\hw3_Q2_sea_surface_zonal_avg.emf",'Resolution',600,'ContentType','auto','BackgroundColor','none','Colorspace','rgb')
    exportgraphics(t_TCL_1,"..\\doc\\fig\\hw3_Q2_sea_surface_zonal_avg.png",'Resolution',600,'ContentType','auto','BackgroundColor','none','Colorspace','rgb')

    return;
end

%% find nearest data point
function [x,y,arclen] = Q2_location_query(lon,lat,nav_lon,nav_lat)
%%% find nearest data point
% OUTPUT:
% x: lon index
% y: lat index
% arclen: the lengths, arclen, of the great circle arcs connecting pairs 
% of points on the surface of a sphere.
    arguments
        lon
        lat
        nav_lon
        nav_lat
    end
    arclen = nan(size(nav_lon));
    for j = 1:size(nav_lon,2)
        arclen(:,j) = distance(lat,lon,nav_lat(:,j),nav_lon(:,j));
    end
    [arclen,I] = min(arclen,[],"all","omitnan","linear");
    [x,y] = ind2sub(size(nav_lon),I);

    return;
end

%%
function [SP,pt,p,SA,CT,N2,p_mid] = Q2_gsw_N2(x,y,nav_lon,nav_lat,deptht,votemper,vosaline,path_votemper,path_vosaline,mem_EN)
% Q2
%
    arguments
        x
        y
        nav_lon
        nav_lat
        deptht
        votemper
        vosaline
        path_votemper = "..\data\votemper_control_monthly_highres_3D_202201_OPER_v0.1.nc";
        path_vosaline = "..\data\vosaline_control_monthly_highres_3D_202201_OPER_v0.1.nc";
        mem_EN string = 'off' % accelerate by storing large data in memory
    end
    
    if (strcmpi(mem_EN,"on"))
        SP = squeeze(vosaline(x,y,:));
        pt = squeeze(votemper(x,y,:));
    else
        SP = squeeze(ncread(path_vosaline,'vosaline',[x,y,1,1],[1,1,Inf,1]));
        pt = squeeze(ncread(path_votemper,'votemper',[x,y,1,1],[1,1,Inf,1]));
    end
    if (nargout < 3)
        p = uint8(0); SA = p; CT = p; N2 = p; p_mid = p;
        return;
    end

    z = -deptht; % [m] Height (z) is NEGATIVE in the ocean.
    lat = nav_lat(x,y); % [deg N]
    p = gsw_p_from_z(z,lat);
    if (nargout < 4)
        SA = uint8(0); CT = SA; N2 = SA; p_mid = SA;
        return;
    end

    if (min(isnan(SP),[],"all"))
        SA = uint8(0); CT = SA; N2 = SA; p_mid = SA;
        warning("warning: no data! (x,y) = (%i,%i)\n",x,y);
        return;
    end
    lon = nav_lon(x,y); % [deg E]
    [SA,in_ocean] = gsw_SA_from_SP(SP,p,lon,lat);
    CT = gsw_CT_from_pt(SA,pt);
    if(~min(gsw_infunnel(SA,CT,p)))
%         warning("warning: not in funnel! (x,y) = (%i,%i)\n",x,y);
    end
    if (nargout < 6)
        N2 = uint8(0); p_mid = N2;
        return;
    end

    [N2,p_mid] = gsw_Nsquared(SA,CT,p,lat);

    return;
end
