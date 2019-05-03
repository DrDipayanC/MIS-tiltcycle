% Created by Dipayan Choudhury on 11/04/2019
% Code to get a stable equilibriated ice sheet profile following Robel et
% al., 2018 paper 
% REPLICATES ROBEL'S 2 STAGE MODEL - even with a noise in the accumulation
% CURRENTLY WORKS GREAT FOR PROGRADE
% CHANGELOG (compared to previous 2 iterations):
% - Accumulation updated to have per year values
% - Bed elevation at Icedivide is not zero

clc
clear

%% The initial conditions
% Initial ice sheet configurations
H_init = 2000 ; %2000 %initial mean thickness of glacier % -1*bg_init_icedivide + 
L_init = 600e3 ; % 450000 %initial length of glacier (terminus position)
bedslope = -1e-3 ; % Negative is prograde
bg_init_icedivide = -100 ; % Bed height at the center of the ice sheet
theta_buttr = 0.6 ; % Currently unused **

rho_i = 917.0 ;
rho_w = 1028.0 ;
rho_b = 2650.0 ; % 2650.0    # average density of rock 
lambda_rho=rho_w/rho_i ;
bg_init = bg_init_icedivide + (bedslope * L_init) ; % Depth of the bed at grounding line % (bedslope * L_init)
hg_init = -1*(lambda_rho * bg_init) ; % Ice Height at Grounding line % Just based on buoyancy
% Accumulation forcing 
accum_mean = 0.5/year ; %0.3 %/year
accum_var = 0.1/year ;
% Noise in Accum 
P_accum = accum_var*randn(num_time,1) + accum_mean ;

num_years = 10e3 ; % Number of steps to run for - Robel ran for 3000,000 years
year = 3600*24*365 ; % Convert year to seconds
% tf = 10e3.*year;
% nt = 10e3;
% dt = tf/nt;
num_time=num_years ; % 
dt=year ;

% Other constants used here
g = 9.81 ;
n = 3.0 ; % Nye-Glen law exponenet
m = 1/3.0 ; % Weertman frictin law exponent
A_glen = 4.227e-25 ; % Nye-Glen law coeff - Converting to per year
C = 7.624e6 ; % Basal friction coeff
theta_buttr = 0.9 ; % 
% accum_mean = 0.5 ; % 0.3 #/year
% accum_var = 0.1; % 

% The grounding line flux is Qg=omega*hg^beta
omega=((A_glen*(rho_i*g)^(n+1) * (theta_buttr*(1-1/lambda_rho))^ n)/(4^n*C) )^(1/(m+1)) ;
beta = (m+n+3)/(m+1) ;
% The interior ice flux
nu= (rho_i*g/C)^n ;
% Q_init = nu*H_init^(2*n+1)/(L_init^n)

%% The parabolic ice sheet profile
x_init=(1:1:L_init)' ;
hx_init = zeros(max(size(x_init)),1) ; % Gets the ice parabolic profile

for i_x = 1 : max(size(x_init))
    hx_init(i_x) = (3*(H_init-hg_init)*sqrt(L_init-i_x))/(2*sqrt(L_init)) + hg_init ;
end

%%
close all ;
set(gcf,'color','w') ;
plot(x_init/1000, hx_init, 'color',[0 0.4470 0.7410],'LineWidth',1.5) ;
hold on
plot([x_init(1)/1000 x_init(end)/1000], [0 0], '-.k')
plot([x_init(1)/1000 x_init(end)/1000], [bg_init_icedivide bg_init],'color' ...
    , [0.8500 0.3250 0.0980], 'LineWidth',1.5) ;
plot([x_init(end)/1000 x_init(end)/1000], [bg_init hg_init], 'color', ...
    [0 0.4470 0.7410],'LineWidth',1.5) ;
set(gca, 'YGrid', 'on', 'XGrid', 'off') ;
xlabel('Grounding Line Length (km)','fontsize',12) ;
ylabel('Ice Height (m)','fontsize',12) ;

%% Running it forward in time

H_t = nan(num_time,1) ;
L_t = nan(num_time,1) ;
hg_t = nan(num_time,1) ;
bg_t = nan(num_time,1) ;
Q_t = nan(num_time,1) ;
Qg_t = nan(num_time,1) ;

H_t(1)=H_init ;
L_t(1)=L_init ;

for i_time = 1 : num_time-1
    
    P = P_accum(i_time) ;
%     P = accum_mean ;
    H = H_t(i_time) ;
    L = L_t(i_time) ;
    bg = (bg_init_icedivide + bedslope * L) ; % abs(bg_init_icedivide + (lambda_rho * bedslope * L)) ;
    hg = -1*((lambda_rho * (bg_init_icedivide + bedslope * L))) ; % abs(bg_init_icedivide + (lambda_rho * bedslope * L)) ;
    Q = nu*H^(2*n+1)/(L^n) ;
    %Q = (rho_i*g/(C*xg))^n * (h^(2*n + 1));
    Qg = omega * hg^beta ;
    
    bg_t(i_time) = bg ;
    hg_t(i_time) = hg ;
    Q_t(i_time) = Q ;
    Qg_t(i_time) = Qg ;
    
    dH = P - Qg/L - H*(Q-Qg)/(hg*L) ; % Sensitive to using P instead of P/year
    dL = (Q-Qg)/hg ;
    
    H_t(i_time+1) = H + dH*dt ;
    L_t(i_time+1) = L + dL*dt ;
end

hg_t(end) = hg ;
bg_t(end) = bg ;

%%
close all;
subplot(2,2,1)
plot(H_t(2:end) - H_t(1:end-1)) ;
ylabel('del H (m)') ;
subplot(2,2,2)
plot(L_t(2:end) - L_t(1:end-1))
ylabel('del L (m)') ;
subplot(2,2,3)
plot(H_t) ;
ylabel('H (m)') ;
subplot(2,2,4)
plot(L_t)
ylabel('L (m)') ;

%% Will make a movie of the ice profile

tic
vidfile = VideoWriter('IceProfile_Prograde.mp4','MPEG-4');
open(vidfile);
for  i_time = 1 : 1 : 1000 % num_time
    
    i_time % Just to print
    H = H_t(i_time) ;
    L = L_t(i_time) ;
    bg = bg_init_icedivide + (bedslope * L) ; % Depth of the bed at grounding line 
    hg = -1*(lambda_rho * bg) ; % Ice Height at the grounding line

    x=(1:1:L)' ;
    hx = zeros(max(size(x)),1) ; % Gets the ice parabolic profile

    for i_x = 1 : max(size(x))
        hx(i_x) = (3*(H-hg)*sqrt(L-i_x))/(2*sqrt(L)) + hg ;
    end
    
    close all ;
    h1 = figure();set(gcf,'Visible', 'off');
    set(gcf,'color','w') ;
    plot(x_init/1000, hx_init, '--', 'color',[0 0.4470 0.7410],'LineWidth',1.5) ;
    hold on
    plot(x/1000, hx, 'color',[0 0.4470 0.7410],'LineWidth',1.5) ;
    plot([x_init(1)/1000 x_init(end)/1000], [0 0], '-.k')
    plot([x_init(1)/1000 x_init(end)/1000], [bg_init_icedivide bg_init],'--', 'color' ...
        , [0.8500 0.3250 0.0980], 'LineWidth',1.5) ;
    plot([x(1)/1000 x(end)/1000], [bg_init_icedivide bg],'color' ...
        , [0.8500 0.3250 0.0980], 'LineWidth',1.5) ;
    plot([x_init(end)/1000 x_init(end)/1000], [bg_init hg_init],'--', 'color', ...
        [0 0.4470 0.7410],'LineWidth',1.5) ;
    plot([x(end)/1000 x(end)/1000], [bg hg], 'color', ...
        [0 0.4470 0.7410],'LineWidth',1.5) ;
    
    set(gca, 'YGrid', 'on', 'XGrid', 'off') ;
    ylim([-1000, 3500]) ;
    xlabel('Grounding Line Length (km)','fontsize',12) ;
    ylabel('Ice Height (m)','fontsize',12) ;
    title(sprintf('Ice Profile at t = %i', i_time)) ;
    
    f=getframe(gcf) ;
    writeVideo(vidfile, f);
end
close(vidfile)
toc





