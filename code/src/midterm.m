%% midterm.m
% Description: MATLAB code for mid-term exam (MS8402, 2022 Spring)
% Author: Guorui Wei (危国锐) (313017602@qq.com; weiguorui@sjtu.edu.cn)
% Student ID: 120034910021
% Created: 2022-03-24
% Last modified: 2022-03-

%% Initialize project

clc; clear; close all
init_env();

%% Problem 1

%
f = 1e-4;           % [s^{-1}]
h_1 = 200;          % [m]
U_1 = 0.5;          % [m/s]
L_1 = 10e3;         % [m]
h_2 = [160,100];    % [m]
%
h_c = L_1.*h_1.*f./(U_1+L_1.*f);
U_2_sup = U_1.*sqrt(1+L_1.*f./U_1);
U_2 = @(h_2) U_1.*sqrt(1-(L_1.*(h_1./h_2 -1).*f)./U_1);
L_2 = @(h_2) L_1.*sqrt((U_1.*h_1.^2)./(h_2.*(U_1.*h_2-(h_1-h_2).*L_1.*f)));
%
fprintf("h_2 = %g:\t U_2 = %.4g, L_2 = %.4g,\n" + ...
    "h_2 = %g:\t U_2 = %.3g, L_2 = %.3g.\n" + ...
    "h_2 = %g:\t U_2 = %.3g, L_2 = %.3g.\n", ...
    h_2(1),U_2(h_2(1)),L_2(h_2(1)), ...
    h_c+eps,U_2(h_c+eps),L_2(h_c+eps), ...
    h_2(2),U_2(h_2(2)),L_2(h_2(2)));
%
figure("Name","Problem 1")
t_TCL_1 = tiledlayout(1,1,"TileSpacing","tight","Padding","tight");
t_Axes_1 = nexttile(t_TCL_1,1);
h_2 = linspace(7*h_c,h_c,6001);
t_plot_1_L2 = plot(t_Axes_1,h_2,L_2(h_2),'-',"DisplayName",'$L_2$');
set(t_Axes_1,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off');
ylabel(t_Axes_1,"$L_2$ (m)","Interpreter",'latex');
xlabel(t_Axes_1,"$h_2$ (m)","Interpreter",'latex');
yyaxis(t_Axes_1,'right');
t_plot_1_U2 = plot(t_Axes_1,h_2,U_2(h_2),'-',"DisplayName",'$U_2$');
ylabel(t_Axes_1,"$U_2$ (m/s)","Interpreter",'latex');
set(t_Axes_1,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off',"XLimitMethod",'tickaligned');
hold on
t_plot_1_hc = plot(t_Axes_1,[h_c,h_c],[0,0],'o',"DisplayName","$h_{\rm{c}}$");
t_plot_1_U2_sup = plot(t_Axes_1,[h_2(ceil(1*length(h_2)/4)),h_2(1)],[U_2_sup,U_2_sup],'-.',"DisplayName","$\sup\limits_{h_2 \ge h_{\rm{c}}}{U_2}$");
hold off
%
legend([t_plot_1_hc,t_plot_1_L2,t_plot_1_U2,t_plot_1_U2_sup],"Location",'east','Interpreter','latex',"Box","off",'FontSize',10);
[t_title_t,t_title_s] = title(t_TCL_1,"\bf 2022 Spring MS8402 Mid-term Exam Pblm 1","Guorui Wei 120034910021","Interpreter",'latex');
set(t_title_s,'FontSize',8)
%
exportgraphics(t_TCL_1,"..\\doc\\fig\\midterm_P1.emf",'Resolution',1200,'ContentType','auto','BackgroundColor','none','Colorspace','rgb')
exportgraphics(t_TCL_1,"..\\doc\\fig\\midterm_P1.png",'Resolution',1200,'ContentType','auto','BackgroundColor','none','Colorspace','rgb')


%% Problem 2

% params
U = 40;             % [m/s]
L = 570e3;          % [m]
beta_0 = 1.61e-11;  % [m^{-1} s^{-1}]
%
beta_plus = beta_0 ./ (U./L.^2);
y_plus = linspace(-5,5,5001);
u_plus = @(y_plus) exp(-y_plus.^2 / 2);
u_plus_second_derivatives_wrt_y_plus = @(y_plus) -(1-y_plus.^2).*exp(-y_plus.^2 / 2);
[y_1,u_1,exitflag_1,output_1] = fzero(@(x) beta_plus-u_plus_second_derivatives_wrt_y_plus(x),[1 1.5]);
[y_2,u_2,exitflag_2,output_2] = fzero(@(x) beta_plus-u_plus_second_derivatives_wrt_y_plus(x),[2.5 3]);
%
t_fig_2 = figure("Name","Problem 2");
t_TCL = tiledlayout(1,2,"TileSpacing","tight","Padding","tight");
t_Axes_a = nexttile(t_TCL,1);
t_plot_a = plot(t_Axes_a,u_plus(y_plus),y_plus,'-',"color",'#0072BD',"DisplayName",'$u_+$');
set(t_Axes_a,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off',"XColor",'black',"YColor",'black');
xlabel(t_Axes_a,"$u_+$","Interpreter",'latex');
%
t_Axes_b = nexttile(t_TCL,2);
t_plot_b = plot(t_Axes_b,beta_plus-u_plus_second_derivatives_wrt_y_plus(y_plus),y_plus,'-',"DisplayName",'$u_+$');
hold on
plot(t_Axes_b,[0,0],t_Axes_b.YLim,'--');
t_plot_b_yA = plot(t_Axes_b,[u_1,-u_1],[y_1,-y_1],'o',"DisplayName",sprintf("$y^{\\rm{A}}_{1,2} = \\pm %.3f$",y_1));
t_plot_b_yB = plot(t_Axes_b,[u_2,-u_2],[y_2,-y_2],'^',"DisplayName",sprintf("$y^{\\rm{B}}_{1,2} = \\pm %.3f$",y_2));
hold off
t_Axes_b.YAxis.Visible = 'on';
legend([t_plot_b_yA,t_plot_b_yB],"Location",'best','Interpreter','latex',"Box","off",'FontSize',10);
set(t_Axes_b,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off',"XColor",'black',"YColor",'black');
xlabel(t_Axes_b,"$\beta_{0+} - \frac{{\rm d}^2 u_+}{({\rm d} y_+)^2}$","Interpreter",'latex');
annotation(t_fig_2,'textbox',[t_Axes_a.Position(1)+t_Axes_a.Position(3)-.07 .13 .07 .06],'String','\bf{(a)}','LineStyle','none','Interpreter','latex','FontWeight','bold','FitBoxToText','on');
annotation(t_fig_2,'textbox',[t_Axes_b.Position(1)+t_Axes_b.Position(3)-.07 .13 .07 .06],'String','\bf{(b)}','LineStyle','none','Interpreter','latex','FontWeight','bold','FitBoxToText','on');
%
ylabel(t_TCL,"$y_+$","Interpreter",'latex');
[t_title_t,t_title_s] = title(t_TCL,"\bf 2022 Spring MS8402 Mid-term Exam Pblm 2","Guorui Wei 120034910021","Interpreter",'latex');
set(t_title_s,'FontSize',8)
colors_default = colororder('default');
colororder(t_Axes_b,colors_default(2:end,:));
%
exportgraphics(t_TCL,"..\\doc\\fig\\midterm_P2.emf",'Resolution',800,'ContentType','auto','BackgroundColor','none','Colorspace','rgb')
exportgraphics(t_TCL,"..\\doc\\fig\\midterm_P2.png",'Resolution',800,'ContentType','auto','BackgroundColor','none','Colorspace','rgb')

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
