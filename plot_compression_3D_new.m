clear all
close all
format long

N = 100;
R = 1;

Nphi = 1000;

phi_vec = linspace(0.01,0.72,Nphi); %PYfreezing range

Nphi_ds = 20;
phi_vec_ds = downsample(phi_vec,Nphi_ds);

a1_cubic = 2*R;
a2_cubic = sqrt(6)*R;
a3_cubic = 2*sqrt(2)*R;
a4_cubic = 4*R;

phi1_cubic = 4/3*pi*R^3/(a1_cubic^3);
phi2_cubic = 4/3*pi*R^3/(a2_cubic^3);
phi3_cubic = 4/3*pi*R^3/(a3_cubic^3);
phi4_cubic = 4/3*pi*R^3/(a4_cubic^3);

a1_fcc = 2*R;
a2_fcc = 2*sqrt(2)*R;
a3_fcc = 2*sqrt(3)*R;
a4_fcc = 4*R;

v1_fcc = 1/sqrt(2)*a1_fcc^3;
v2_fcc = 1/sqrt(2)*a2_fcc^3;
v3_fcc = 1/sqrt(2)*a3_fcc^3;
v4_fcc = 1/sqrt(2)*a4_fcc^3;

phi1_fcc = 4/3*pi*R^3/v1_fcc;
phi2_fcc = 4/3*pi*R^3/v2_fcc;
phi3_fcc = 4/3*pi*R^3/v3_fcc;
phi4_fcc = 4/3*pi*R^3/v4_fcc;

vor_vec = (4/3*pi*R^3)./phi_vec;

F_vec_cubic = my_F_3D_cubic(vor_vec,R,N);
dF_vec_cubic = my_F_deriv_3D_cubic(vor_vec,R,N);

F_vec_fcc = my_F_3D_fcc(vor_vec,R,N);
dF_vec_fcc = my_F_deriv_3D_fcc(vor_vec,R,N);

F_vec_bcc = my_F_3D_bcc(vor_vec,R,N);
dF_vec_bcc = my_F_deriv_3D_bcc(vor_vec,R,N);

f_py = phi_vec/(4/3*pi*R^3).*(-1+log(6*phi_vec/pi)-log(1-phi_vec)+1.5*(1./(1-phi_vec).^2-1));

dim=3;
crit_dens = (4/3*pi*R^3)/v1_fcc;
F_vec2 = log(phi_vec)-1-log(pi/6)-log(N*vor_vec)+...
    4/3*pi*R^3*2^dim*((crit_dens./((4/3*pi*R^3)./vor_vec)).^(1/dim)-1).^dim./crit_dens;
dF_vec2 = -2./vor_vec+2^dim*((crit_dens./((4/3*pi*R^3)./vor_vec)).^(1/dim)-1).^dim...
    +2^dim*((crit_dens./((4/3*pi*R^3)./vor_vec)).^(1/dim)-1).^(dim-1).*(crit_dens./((4/3*pi*R^3)./vor_vec)).^(1/dim-1);

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


dfdrho_fcc = -(-3*log(2*R)+log(F_vec_fcc))+...
    +vor_vec./(F_vec_fcc).*dF_vec_fcc;
dfdrho_bcc = -(-3*log(2*R)+log(F_vec_bcc))+...
    +vor_vec./(F_vec_bcc).*dF_vec_bcc;
dfdrho_cubic = -(-3*log(2*R)+log(F_vec_cubic))+...
    +vor_vec./F_vec_cubic.*dF_vec_cubic;
f_fcc = -phi_vec/(4/3*pi*R^3).*(-3*log(2*R)+log(F_vec_fcc));
f_bcc = -phi_vec/(4/3*pi*R^3).*(-3*log(2*R)+log(F_vec_bcc));
f_cubic = -phi_vec/(4/3*pi*R^3).*(-3*log(2*R)+log(F_vec_cubic));

df_fcc_fun = @(x) interp1(phi_vec/(4/3*pi*R^3),dfdrho_fcc,x);
df_bcc_fun = @(x) interp1(phi_vec/(4/3*pi*R^3),dfdrho_bcc,x);
df_cubic_fun = @(x) interp1(phi_vec/(4/3*pi*R^3),dfdrho_cubic,x);
f_fcc_fun = @(x) interp1(phi_vec/(4/3*pi*R^3),f_fcc,x);
f_bcc_fun = @(x) interp1(phi_vec/(4/3*pi*R^3),f_bcc,x);
f_cubic_fun = @(x) interp1(phi_vec/(4/3*pi*R^3),f_cubic,x);

f_ct = phi_vec/(4/3*pi*R^3).*(log(sqrt(2))-3*log((pi*sqrt(2)./(6*phi_vec)).^(1/3)-1));

dfdrho_ct = (log(sqrt(2))-3*log((pi*sqrt(2)./(6*phi_vec)).^(1/3)-1))...
    +phi_vec.*(1./((pi*sqrt(2)./(6*phi_vec)).^(1/3)-1).*(pi*sqrt(2)./(6*phi_vec)).^(4/3)*(6/(pi*sqrt(2))));

f_ct_fun = @(x) interp1(phi_vec/(4/3*pi*R^3),f_ct,x);
df_ct_fun = @(x) interp1(phi_vec/(4/3*pi*R^3),dfdrho_ct,x);

f_py = phi_vec/(4/3*pi*R^3).*(-1+log(6*phi_vec/pi)-log(1-phi_vec)+1.5*(1./(1-phi_vec).^2-1));
dfdrho_py = (-1+log(6*phi_vec/pi)-log(1-phi_vec)+1.5*(1./(1-phi_vec).^2-1))...
    +phi_vec.*(1./phi_vec+1./(1-phi_vec)+3./(1-phi_vec).^3);
f_py_fun = @(x) interp1(phi_vec/(4/3*pi*R^3),f_py,x);
df_py_fun = @(x) interp1(phi_vec/(4/3*pi*R^3),dfdrho_py,x);

f_cs = phi_vec/(4/3*pi*R^3).*(-1+log(6*phi_vec/pi)+(1./(1-phi_vec).^2-1)+2*(1./(1-phi_vec)-1)); %Carnahan-Starling

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
ylim([-0.1 4])
ylim([-0.1 3])

f_cubic_w1 = -phi_vec/(4/3*pi*R^3).*(1-log(4/3*pi*R^3)-3*log(2)+log(pi/6)+log(F_vec_cubic));
f_fcc_w1 = -phi_vec/(4/3*pi*R^3).*(1-log(4/3*pi*R^3)-3*log(2)+log(pi/6)+log(F_vec_fcc));
f_bcc_w1 = -phi_vec/(4/3*pi*R^3).*(1-log(4/3*pi*R^3)-3*log(2)+log(pi/6)+log(F_vec_bcc));
f_cubic_w1 = -phi_vec/(4/3*pi*R^3).*(1-3*log(2*R)+log(F_vec_cubic));
f_fcc_w1 = -phi_vec/(4/3*pi*R^3).*(1-3*log(2*R)+log(F_vec_fcc));
f_bcc_w1 = -phi_vec/(4/3*pi*R^3).*(1-3*log(2*R)+log(F_vec_bcc));

Lsave_cubic = [f_cubic; f_cubic_w1];
Lsave_fcc = [f_fcc; f_fcc_w1];
Lsave_bcc = [f_bcc; f_bcc_w1];

Lsave_cubic_mean = mean([f_cubic; f_cubic_w1],1);
Lsave_fcc_mean = mean([f_fcc; f_fcc_w1],1);
Lsave_bcc_mean = mean([f_bcc; f_bcc_w1],1);

%HACK to show right range
Lsave_cubic = Lsave_cubic_mean+([f_cubic; f_cubic_w1]-Lsave_cubic_mean)/sqrt(2);
Lsave_fcc = Lsave_fcc_mean+([f_fcc; f_fcc_w1]-Lsave_fcc_mean)/sqrt(2);
Lsave_bcc = Lsave_bcc_mean+([f_bcc; f_bcc_w1]-Lsave_bcc_mean)/sqrt(2);

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


%% PY-BCC intersection (no solution depending on factor of e)
myfun = @(x) fsolve_fun(x,f_py_fun,f_bcc_fun,df_py_fun,df_bcc_fun);
x0 = [.55 .7]/(4/3*pi*R^3);
x0 = [.6 .67]/(4/3*pi*R^3);
x0 = [.45 .55]/(4/3*pi*R^3);
x = fsolve(myfun,x0);

'freezing is'
x*(4/3*pi*R^3)

figure(21);clf;
plot(phi_vec,f_bcc,'k','LineWidth',2)
hold on
plot(phi_vec,f_py,'r','LineWidth',2)
% plot(phi_vec,f_ct,'g','LineWidth',2)
% plot(phi_vec,f_cubic,'b','LineWidth',2)

scatter(x(1)*(4/3*pi*R^3),f_py_fun(x(1)),100,'r');
scatter(x(2)*(4/3*pi*R^3),f_bcc_fun(x(2)),100,'k');
Nparams = 100;
params = linspace(-0.5,1.5,Nparams);
param_fn = f_py_fun(x(1))+params*(x(2)-x(1))*df_py_fun(x(1));
plot((x(1)+params*(x(2)-x(1)))*(4/3*pi*R^3),param_fn,'b','LineWidth',1);
% plot(phi_vec,dfdrho_bcc,'k','LineWidth',2)
% hold on
% plot(phi_vec,dfdrho_cubic,'r','LineWidth',2)
xlim([0 max(phi_vec)])
legend('BCC','Percus-Yevick','Location','Best');
title('Common Tangent Construction');

set(gca, 'fontsize', 18);


function F = fsolve_fun(x,f1,f2,df1,df2)
  x1 = x(1); x2 = x(2);
  out1 = df1(x1)-df2(x2);
  out2 = f1(x1)-f2(x2)+(x2-x1)*df1(x1);
  F = [out1 out2];
end

function F = trans_fun(xvec,x_crit,f1_fun,f2_fun)
  F = zeros(size(xvec));
  for i=1:numel(xvec)
      x=xvec(i);
  if (x < x_crit)
      F(i) = f1_fun(x);
  else
      F(i) = f2_fun(x);
  end
  end
end