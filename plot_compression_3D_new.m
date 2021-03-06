clear all
close all
format long

N = 100; % number of spheres
R = 1; % radius of each sphere

Nphi = 1000; % number of packing fractions to plot
phi_vec = linspace(0.01,0.72,Nphi); % range of packing fractions
vor_vec = (4/3*pi*R^3)./phi_vec; % corresponding range of Voronoi volumes

%% set constants

% critical lattice spacings
a1_fcc = 2*R; % lattice spacing at which exclusion spheres just intersect
a2_fcc = 2*sqrt(2)*R; % lattice spacing at percolation
a3_fcc = 2*sqrt(3)*R; % lattice spacing at leaky packing fraction
a4_fcc = 4*R; % lattice spacing at close packing fraction

% corresponding Voronoi volumes
v1_fcc = 1/sqrt(2)*a1_fcc^3;
v2_fcc = 1/sqrt(2)*a2_fcc^3;
v3_fcc = 1/sqrt(2)*a3_fcc^3;
v4_fcc = 1/sqrt(2)*a4_fcc^3;

% corresponding packing fractions
phi1_fcc = 4/3*pi*R^3/v1_fcc;
phi2_fcc = 4/3*pi*R^3/v2_fcc;
phi3_fcc = 4/3*pi*R^3/v3_fcc;
phi4_fcc = 4/3*pi*R^3/v4_fcc;

%% compute quantities of interest from cell theory

% free volumes
F_vec_cubic = my_F_3D_cubic(vor_vec,R); %free volume in SC lattice
dF_vec_cubic = my_F_deriv_3D_cubic(vor_vec,R); %derivative of free volume in SC lattice

F_vec_fcc = my_F_3D_fcc(vor_vec,R); %free volume in FCC lattice
dF_vec_fcc = my_F_deriv_3D_fcc(vor_vec,R); %derivative of free volume in FCC lattice

F_vec_bcc = my_F_3D_bcc(vor_vec,R); %free volume in BCC lattice
dF_vec_bcc = my_F_deriv_3D_bcc(vor_vec,R); %derivative of free volume in BCC lattice

% free energy densities
f_fcc = -phi_vec/(4/3*pi*R^3).*(-3*log(2*R)+log(F_vec_fcc));
f_bcc = -phi_vec/(4/3*pi*R^3).*(-3*log(2*R)+log(F_vec_bcc));
f_cubic = -phi_vec/(4/3*pi*R^3).*(-3*log(2*R)+log(F_vec_cubic));

% alternative free energy density with calibration to low-density
f_cubic_w1 = -phi_vec/(4/3*pi*R^3).*(1-3*log(2*R)+log(F_vec_cubic));
f_fcc_w1 = -phi_vec/(4/3*pi*R^3).*(1-3*log(2*R)+log(F_vec_fcc));
f_bcc_w1 = -phi_vec/(4/3*pi*R^3).*(1-3*log(2*R)+log(F_vec_bcc));

% derivatives of free energy densities
dfdrho_fcc = -(-3*log(2*R)+log(F_vec_fcc))+...
    +vor_vec./(F_vec_fcc).*dF_vec_fcc;
dfdrho_bcc = -(-3*log(2*R)+log(F_vec_bcc))+...
    +vor_vec./(F_vec_bcc).*dF_vec_bcc;
dfdrho_cubic = -(-3*log(2*R)+log(F_vec_cubic))+...
    +vor_vec./F_vec_cubic.*dF_vec_cubic;

% corresponding function handles
df_fcc_fun = @(x) interp1(phi_vec/(4/3*pi*R^3),dfdrho_fcc,x);
df_bcc_fun = @(x) interp1(phi_vec/(4/3*pi*R^3),dfdrho_bcc,x);
df_cubic_fun = @(x) interp1(phi_vec/(4/3*pi*R^3),dfdrho_cubic,x);
f_fcc_fun = @(x) interp1(phi_vec/(4/3*pi*R^3),f_fcc,x);
f_bcc_fun = @(x) interp1(phi_vec/(4/3*pi*R^3),f_bcc,x);
f_cubic_fun = @(x) interp1(phi_vec/(4/3*pi*R^3),f_cubic,x);

% Percus-Yevick equation of state
f_py = phi_vec/(4/3*pi*R^3).*(-1+log(6*phi_vec/pi)-log(1-phi_vec)+1.5*(1./(1-phi_vec).^2-1));
dfdrho_py = (-1+log(6*phi_vec/pi)-log(1-phi_vec)+1.5*(1./(1-phi_vec).^2-1))...
    +phi_vec.*(1./phi_vec+1./(1-phi_vec)+3./(1-phi_vec).^3);
f_py_fun = @(x) interp1(phi_vec/(4/3*pi*R^3),f_py,x);
df_py_fun = @(x) interp1(phi_vec/(4/3*pi*R^3),dfdrho_py,x);

% Carnahan-Starling equation of state
f_cs = phi_vec/(4/3*pi*R^3).*(-1+log(6*phi_vec/pi)+(1./(1-phi_vec).^2-1)+2*(1./(1-phi_vec)-1)); %Carnahan-Starling

%% make plots

% free volumes
figure(1);clf;
semilogy(phi_vec,F_vec_fcc,'k','LineWidth',2);
hold on
semilogy(phi_vec,F_vec_bcc,'b','LineWidth',2);
semilogy(phi_vec,F_vec_cubic,'r','LineWidth',2);
legend('FCC','BCC','SC','Location','Best')
title('Free Volume');
set(gca, 'fontsize', 18);
xlim([0 0.72])
ylim([10^-6 10^4])

% compressibility factors
figure(3);clf;
C1pf=plot(phi_vec,vor_vec./F_vec_fcc.*dF_vec_fcc,'k','LineWidth',2);
hold on
C1pb=plot(phi_vec,vor_vec./F_vec_bcc.*dF_vec_bcc,'b','LineWidth',2);
C1pc=plot(phi_vec,vor_vec./F_vec_cubic.*dF_vec_cubic,'r','LineWidth',2);
ylim([0 20])
CSpy = (1+phi_vec+phi_vec.^2)./(1-phi_vec).^3;
CSpy=plot(phi_vec,CSpy,'--','Color',[0.4940 0.1840 0.5560],'LineWidth',2);
CS = (1+phi_vec+phi_vec.^2-phi_vec.^3)./(1-phi_vec).^3;
CSp=plot(phi_vec,CS,'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2);
legend([C1pf C1pb C1pc CSpy CSp],'FCC','BCC','SC','PY','CS','Location','NW');
title('Compressibility');
set(gca, 'fontsize', 18);

% free energy densities
figure(11);clf;
plot(phi_vec,f_fcc,'k','LineWidth',2)
hold on
plot(phi_vec,f_bcc,'b','LineWidth',2)
plot(phi_vec,f_cubic,'r','LineWidth',2)
plot(phi_vec,f_py,'--','Color',[0.4940 0.1840 0.5560],'LineWidth',2)
plot(phi_vec,f_cs,'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2);
legend('FCC','BCC','SC','PY','CS','Location','Best');
title('Free Energy Density');
set(gca, 'fontsize', 18);
xlim([phi_vec(1) phi_vec(end)])
ylim([-0.1 3])

% concatenate free energy densities obtained using both
% high-density and low-density calibrations
Lsave_cubic = [f_cubic; f_cubic_w1];
Lsave_fcc = [f_fcc; f_fcc_w1];
Lsave_bcc = [f_bcc; f_bcc_w1];

Lsave_cubic_mean = mean([f_cubic; f_cubic_w1],1);
Lsave_fcc_mean = mean([f_fcc; f_fcc_w1],1);
Lsave_bcc_mean = mean([f_bcc; f_bcc_w1],1);

% Workaround to show right range (typically, shadedErrorBar shows the
% standard deviation, but here we want to show the actual range between two limits)
% The shadedErrorBar package is from https://github.com/raacampbell/shadedErrorBar
Lsave_cubic = Lsave_cubic_mean+([f_cubic; f_cubic_w1]-Lsave_cubic_mean)/sqrt(2);
Lsave_fcc = Lsave_fcc_mean+([f_fcc; f_fcc_w1]-Lsave_fcc_mean)/sqrt(2);
Lsave_bcc = Lsave_bcc_mean+([f_bcc; f_bcc_w1]-Lsave_bcc_mean)/sqrt(2);

% plot comparison between free energy densities
figure(25);clf;
plot(phi_vec,Lsave_cubic_mean,'r','LineWidth',1.5)
leg_cubic = findobj(gca, 'Type', 'Line', 'Color', 'r');
trans = 0.3;
leg_cubic.Color(4) = trans;
hold on
plot(phi_vec,Lsave_fcc_mean,'k','LineWidth',1.5)
leg_fcc = findobj(gca, 'Type', 'Line', 'Color', 'k');
trans = 0.3;
leg_fcc.Color(4) = trans;
plot(phi_vec,Lsave_bcc_mean,'b','LineWidth',1.5)
leg_bcc = findobj(gca, 'Type', 'Line', 'Color', 'b');
trans = 0.3;
leg_bcc.Color(4) = trans;
sh_cubic=shadedErrorBar(phi_vec,Lsave_cubic,{@mean,@std},'lineprops','-r','patchSaturation',0.06);
g_cubic=hggroup;
set(sh_cubic.patch,'Parent',g_cubic);
set(sh_cubic.mainLine,'Parent',g_cubic);
delete(sh_cubic.mainLine);
sh_cubic.edge = Lsave_cubic;
hold on
sh_fcc=shadedErrorBar(phi_vec,Lsave_fcc,{@mean,@std},'lineprops','-k','patchSaturation',0.06);
g_fcc=hggroup;
set(sh_fcc.patch,'Parent',g_fcc);
set(sh_fcc.mainLine,'Parent',g_fcc);
sh_fcc.edge = Lsave_fcc;
delete(sh_fcc.mainLine);
sh_bcc=shadedErrorBar(phi_vec,Lsave_bcc,{@mean,@std},'lineprops','-b','patchSaturation',0.06);
g_bcc=hggroup;
set(sh_bcc.patch,'Parent',g_bcc);
set(sh_bcc.mainLine,'Parent',g_bcc);
sh_bcc.edge = Lsave_bcc;
delete(sh_bcc.mainLine);
leg_py=plot(phi_vec,f_py,'--','Color',[0.4940 0.1840 0.5560],'LineWidth',2);
leg_cs=plot(phi_vec,f_cs,'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2);
set(gca, 'fontsize', 18);
[leg,hl]=legend([leg_cubic leg_fcc leg_bcc leg_py leg_cs],'SC','FCC','BCC','PY','CS','location', 'Best');
title('Effect of Calibration');
 LineInLegend = findobj(hl, 'type', 'Line','Color','k');
LineInLegend(1).Color(4) = trans;
 LineInLegend = findobj(hl, 'type', 'Line','Color','b');
LineInLegend(1).Color(4) = trans;
 LineInLegend = findobj(hl, 'type', 'Line','Color','r');
LineInLegend(1).Color(4) = trans;
xlim([phi_vec(1) phi_vec(end)])
ylim([-0.1 4])


%% PY-BCC intersection 
% here we use the common tangent construction to identify the coexistence
% between liquid and solid phases, given respectively by the
% Percus-Yevick and cell theory equations of state
% (The existence of a common tangent depends on whether the
% high-density or low-density calibration is used
myfun = @(x) fsolve_fun(x,f_py_fun,f_bcc_fun,df_py_fun,df_bcc_fun);
x0 = [.45 .55]/(4/3*pi*R^3); %initial guess of common tangent
x = fsolve(myfun,x0); %solve for common tangent

% output the phase transition identified
'freezing is at liquid and solid packing fractions of'
x*(4/3*pi*R^3)

% plot the common tangent
figure(21);clf;
plot(phi_vec,f_bcc,'k','LineWidth',2)
hold on
plot(phi_vec,f_py,'r','LineWidth',2)
scatter(x(1)*(4/3*pi*R^3),f_py_fun(x(1)),100,'r');
scatter(x(2)*(4/3*pi*R^3),f_bcc_fun(x(2)),100,'k');
%draw a line to illustrate the common tangent
Nparams = 100; %number of points to plot in line
params = linspace(-0.5,1.5,Nparams);
param_fn = f_py_fun(x(1))+params*(x(2)-x(1))*df_py_fun(x(1));
plot((x(1)+params*(x(2)-x(1)))*(4/3*pi*R^3),param_fn,'b','LineWidth',1);
xlim([0 max(phi_vec)])
legend('BCC','Percus-Yevick','Location','Best');
title('Common Tangent Construction');
set(gca, 'fontsize', 18);

%function used to identify the common tangent
function F = fsolve_fun(x,f1,f2,df1,df2)
  x1 = x(1); x2 = x(2);
  out1 = df1(x1)-df2(x2);
  out2 = f1(x1)-f2(x2)+(x2-x1)*df1(x1);
  F = [out1 out2];
end