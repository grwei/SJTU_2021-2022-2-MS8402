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
t_TCL_pt_SP = tiledlayout(1,1,"TileSpacing","tight","Padding","tight");
t_Axes = nexttile(t_TCL_pt_SP,1);
[~,c_contour] = contour(t_Axes,X,Y,psi,"Fill","off");
clabel([],c_contour,c_contour.TextList([1,idivide(end+1,uint8(2)),end]),"Interpreter",'latex')
c = colorbar(t_Axes,"Location","eastoutside","TickLabelInterpreter","latex");
c.Label.String = "\psi";
set(t_Axes,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off',"TickDir","out");
ylabel(t_Axes,"$y$","Interpreter",'latex');
xlabel(t_Axes,"$x$","Interpreter",'latex');
title(t_Axes,sprintf("$\\psi = (1-x-\\mathrm{e}^{-x/\\varepsilon}) \\pi \\tau_0 \\sin{\\pi y},\\quad \\tau_0 = %.3g$, $\\varepsilon = %.3g.$",tau_0,epsilon),"Interpreter",'latex')
[~,t_title_s] = title(t_TCL_pt_SP,"\bf 2022 Spring MS8402 Homework 3 Q1","Guorui Wei 120034910021","Interpreter",'latex');
set(t_title_s,'FontSize',8)
exportgraphics(t_TCL_pt_SP,"..\\doc\\fig\\hw3_Q1.emf",'Resolution',600,'ContentType','auto','BackgroundColor','none','Colorspace','rgb')
exportgraphics(t_TCL_pt_SP,"..\\doc\\fig\\hw3_Q1.png",'Resolution',600,'ContentType','auto','BackgroundColor','none','Colorspace','rgb')

%% Question 2: init

% import data
clear; clc; close all;
path_votemper = "global-reanalysis-phy-001-031-grepv2-mnstd-monthly_1649316520167.nc"; % The temperature of a parcel of sea water would have if moved adiabatically to sea level pressure. This variable is a 3D field.
path_vosaline = path_votemper; % The salt content of sea water as measured on the practical salinity units (PSU) scale. This variable is a 3D field.
finfo = ncinfo(path_votemper);
vec_lat = ncread(path_votemper,'latitude');
vec_lon = ncread(path_votemper,'longitude');
[nav_lat,nav_lon] = meshgrid(vec_lat,vec_lon); 
deptht = ncread(path_votemper,'depth');
votemper = ncread(path_votemper,'thetao_mean',[1,1,1,1],[Inf,Inf,Inf,1]);
vosaline = ncread(path_vosaline,'so_mean',[1,1,1,1],[Inf,Inf,Inf,1]);
%
n_lat_bin = 255;

Q2_zonal_average(n_lat_bin,'on');
%
lat_center = [6,45,74,82];
lat_tol = 0.5;
%
ind_depth_thrs = find(deptht > 750,1);
TF_is_ocean = ~isnan(vosaline(:,:,ind_depth_thrs));
ind_N2_min = nan(size(lat_center));
ind_N2_max = ind_N2_min;
ind_N2_avg_min = ind_N2_min;
ind_N2_avg_max = ind_N2_min;
val_N2_min = ind_N2_min;
val_N2_max = ind_N2_min;
val_N2_avg_min = ind_N2_min;
val_N2_avg_max = ind_N2_min;
for i = 1:length(lat_center)
    [ind_N2_min(i),ind_N2_max(i),ind_N2_avg_min(i),ind_N2_avg_max(i),val_N2_min(i),val_N2_max(i),val_N2_avg_min(i),val_N2_avg_max(i)] = Q2_gsw_N2_min_max(lat_center(i),lat_tol,TF_is_ocean,nav_lon,nav_lat,deptht,votemper,vosaline);
end

%% Question 2

ind_tropi = ind_N2_max(1);
ind_mid = ind_N2_max(2);
ind_polar = ind_N2_max(4);
ind_vec = [ind_tropi,ind_mid,ind_polar];
%
SP = nan(length(deptht),length(ind_vec));
pt = SP;
p = SP;
SA = SP;
CT = SP;
pot_rho = SP;
N2 = nan(length(deptht)-1,length(ind_vec));
p_mid = N2;
for i = 1:length(ind_vec)
    [x,y] = ind2sub(size(nav_lon),ind_vec(i));
    [SP(:,i),pt(:,i),p(:,i),SA(:,i),CT(:,i),N2(:,i),p_mid(:,i),pot_rho(:,i)] = Q2_gsw_N2(x,y,nav_lon,nav_lat,deptht,votemper,vosaline);
end
%
figure("Name","Q2_pt_SP_vertical_profile")
t_TCL_pt_SP = tiledlayout(1,2,"TileSpacing","tight","Padding","tight");
t_Axes_pt = nexttile(t_TCL_pt_SP);
t_Axes_pt = Q2_plot_vertical(t_Axes_pt,pt,p,ind_vec,nav_lat,nav_lon,"potential temperature $(^{\circ}{\rm{C}})$","pressure (dbar)","temperature");
t_Axes_SP = nexttile(t_TCL_pt_SP);
t_Axes_SP = Q2_plot_vertical(t_Axes_SP,SP,p,ind_vec,nav_lat,nav_lon,"Salinity (psu)","pressure (dbar)","Salinity");
%
[~,t_title_s] = title(t_TCL_pt_SP,"\bf 2022 Spring MS8402 Homework 3 Q2","Guorui Wei 120034910021","Interpreter",'latex');
set(t_title_s,'FontSize',8)
%
exportgraphics(t_TCL_pt_SP,"..\\doc\\fig\\hw3_Q2_pt_SP_vertical_profile.emf",'Resolution',600,'ContentType','auto','BackgroundColor','none','Colorspace','rgb')
exportgraphics(t_TCL_pt_SP,"..\\doc\\fig\\hw3_Q2_pt_SP_vertical_profile.png",'Resolution',600,'ContentType','auto','BackgroundColor','none','Colorspace','rgb')
%
figure("Name","Q2_rho_N2_vertical_profile")
t_TCL_rho_N2 = tiledlayout(1,2,"TileSpacing","tight","Padding","tight");
t_Axes_rho = nexttile(t_TCL_rho_N2);
t_Axes_rho = Q2_plot_vertical(t_Axes_rho,pot_rho,p,ind_vec,nav_lat,nav_lon,"potential density $(\rm{kg}/\rm{m}^3)$","pressure (dbar)","density");
t_Axes_N2 = nexttile(t_TCL_rho_N2);
t_Axes_N2 = Q2_plot_vertical(t_Axes_N2,N2,p_mid,ind_vec,nav_lat,nav_lon,"$N^2$ $(\rm{rad}^2 / s^2)$","pressure (dbar)","Stratification frequency");
%
[~,t_title_s] = title(t_TCL_rho_N2,"\bf 2022 Spring MS8402 Homework 3 Q2","Guorui Wei 120034910021","Interpreter",'latex');
set(t_title_s,'FontSize',8)
%
exportgraphics(t_TCL_rho_N2,"..\\doc\\fig\\hw3_Q2_rho_N2_vertical_profile.emf",'Resolution',600,'ContentType','auto','BackgroundColor','none','Colorspace','rgb')
exportgraphics(t_TCL_rho_N2,"..\\doc\\fig\\hw3_Q2_rho_N2_vertical_profile.png",'Resolution',600,'ContentType','auto','BackgroundColor','none','Colorspace','rgb')

%% local functions

%% Initialize environment
function [] = init_env()
% Initialize environment
%
    % set up project directory
    if ~isfolder("../doc/fig/")
        mkdir ../doc/fig/
    end
    % configure searching path
    mfile_fullpath = mfilename('fullpath'); % the full path and name of the file in which the call occurs, not including the filename extension.
    mfile_fullpath_without_fname = mfile_fullpath(1:end-strlength(mfilename));
    addpath(genpath(mfile_fullpath_without_fname + "../data"), ...
            genpath(mfile_fullpath_without_fname + "../inc")); % adds the specified folders to the top of the search path for the current MATLAB® session.

    return;
end

%% surface zonal average 
function [] = Q2_zonal_average(n_lat_bin,fig_EN)
% hw3 Q2
%
    arguments
        n_lat_bin
        fig_EN string = 'on'
    end
    
    nav_lat = ncread("..\data\sosstsst_control_monthly_highres_2D_202201_OPER_v0.1.nc",'nav_lat');
    [lat_bin_num,lat_edges] = discretize(nav_lat,n_lat_bin);
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
% find nearest data point
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
function [SP,pt,p,SA,CT,N2,p_mid,pot_rho] = Q2_gsw_N2(x,y,nav_lon,nav_lat,deptht,votemper,vosaline,p_ref,mem_EN,path_votemper,path_vosaline)
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
        p_ref = 0; % reference pressure of potential density
        mem_EN string = 'on' % accelerate by storing large data in memory
        path_votemper = "..\data\votemper_control_monthly_highres_3D_202201_OPER_v0.1.nc";
        path_vosaline = "..\data\vosaline_control_monthly_highres_3D_202201_OPER_v0.1.nc";
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
%     if(~min(gsw_infunnel(SA,CT,p)))
%         warning("warning: not in funnel! (x,y) = (%i,%i)\n",x,y);
%     end
    if (nargout < 6)
        N2 = uint8(0); p_mid = N2;
        return;
    end
    
    [N2,p_mid] = gsw_Nsquared(SA,CT,p,lat);
    if (nargout < 8)
        return;
    end

    pot_rho = gsw_rho(SA,CT,p_ref);

    return;
end

%%
function [ind_N2_min,ind_N2_max,ind_N2_avg_min,ind_N2_avg_max,val_N2_min,val_N2_max,val_N2_avg_min,val_N2_avg_max] = Q2_gsw_N2_min_max(lat_center,lat_tol,TF_is_ocean,nav_lon,nav_lat,deptht,votemper,vosaline,p_ref,mem_EN,path_votemper,path_vosaline)
% Q2
%
    arguments
        lat_center
        lat_tol
        TF_is_ocean
        nav_lon
        nav_lat
        deptht
        votemper
        vosaline
        p_ref = 0; % reference pressure of potential density
        mem_EN string = 'on' % accelerate by storing large data in memory
        path_votemper = "..\data\votemper_control_monthly_highres_3D_202201_OPER_v0.1.nc";
        path_vosaline = "..\data\vosaline_control_monthly_highres_3D_202201_OPER_v0.1.nc";
    end

    % params
    depth_max = 750;
    %
    tStart = tic;
    loc_ind_linear = find(abs(nav_lat-lat_center) < lat_tol & TF_is_ocean);
    N2_min = nan(size(loc_ind_linear));
    N2_max = N2_min;
    N2_avg = N2_min;
    for i = 1:length(loc_ind_linear)
        [x,y] = ind2sub(size(nav_lon),loc_ind_linear(i));
        [~,~,~,~,~,N2,~,~] = Q2_gsw_N2(x,y,nav_lon,nav_lat,deptht,votemper,vosaline,p_ref,mem_EN,path_votemper,path_vosaline);
        N2_min(i) = min(N2(deptht<depth_max),[],"omitnan");
        N2_max(i) = max(N2(deptht<depth_max),[],"omitnan");
        N2_avg(i) = mean(N2(deptht<depth_max),"omitnan");
%         fprintf("\rDone: %i/%i, (%.2f N, %.2f E), N2_min = %.2g, N2_max = %.2g.", ...
%             i,length(loc_ind_linear),nav_lat(loc_ind_linear(i)),nav_lon(loc_ind_linear(i)),N2_min(i),N2_max(i));
    end
    [val_N2_min,ind_min] = min(N2_min,[],"omitnan",'linear'); 
    [val_N2_max,ind_max] = max(N2_max,[],"omitnan",'linear'); 
    [val_N2_avg_min,ind_avg_min] = min(N2_avg,[],"omitnan",'linear'); 
    [val_N2_avg_max,ind_avg_max] = max(N2_avg,[],"omitnan",'linear'); 
    ind_N2_min = loc_ind_linear(ind_min);
    ind_N2_max = loc_ind_linear(ind_max);
    ind_N2_avg_min = loc_ind_linear(ind_avg_min);
    ind_N2_avg_max = loc_ind_linear(ind_avg_max);
    fprintf("\nSummary: %.1f secs used, %i points processed.\n" + ...
        "lat = %.2f N, lat_tol = %.2f deg, depth_max = %.2f m.\n" + ...
        "N2_min = %.2e (%.2f N, %.2f E),\n" + ...
        "N2_max = %.2e (%.2f N, %.2f E),\n" + ...
        "N2_avg_min = %.2e (%.2f N, %.2f E),\n" + ...
        "N2_avg_max = %.2e (%.2f N, %.2f E).\n", ...
        toc(tStart),length(loc_ind_linear), ...
        lat_center,lat_tol,depth_max, ...
        val_N2_min,nav_lat(ind_N2_min),nav_lon(ind_N2_min), ...
        val_N2_max,nav_lat(ind_N2_max),nav_lon(ind_N2_max), ...
        val_N2_avg_min,nav_lat(ind_N2_avg_min),nav_lon(ind_N2_avg_min), ...
        val_N2_avg_max,nav_lat(ind_N2_avg_max),nav_lon(ind_N2_avg_max));

    return;
end

%%
function [t_Axes] = Q2_plot_vertical(t_Axes,x_data,y_data,ind_vec,nav_lat,nav_lon,xlabel_str,ylabel_str,axes_title_str)
% Q2
%
    arguments
        t_Axes
        x_data
        y_data
        ind_vec
        nav_lat
        nav_lon
        xlabel_str
        ylabel_str
        axes_title_str
    end

    hold on
    for i = 1:length(ind_vec)
        t_plot_pt = plot(t_Axes,x_data(:,i),y_data(:,i),'-',"DisplayName",sprintf("%.2f$^{\\circ}{\\rm{N}}$, %.2f$^{\\circ}{\\rm{E}}$",nav_lat(ind_vec(i)),nav_lon(ind_vec(i))));
    end
    hold off
    grid on
    set(t_Axes,"YDir",'reverse',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off');
    xlabel(t_Axes,xlabel_str,"Interpreter",'latex');
    ylabel(t_Axes,ylabel_str,"Interpreter",'latex');
    legend(t_Axes,"Location",'best','Interpreter','latex',"Box","off",'FontSize',10);
%     title(t_Axes,axes_title_str,"Interpreter",'latex')
    
    return;
end
