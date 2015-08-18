%   Second Year Paper: Running Programs of Correia et. al (2013) with
%   decentralized U.S. fiscal federalism

%   Andrew Fieldhouse
%   Cornell University
%   June 30, 2015

%   Analysis of Unconventional Fiscal Policy

%   Section IV: Unconventional Fiscal Policy

clear all
clear cfc

addpath(genpath('/Applications/Dynare/4.4.3/matlab'))
addpath('/Users/afieldhouse/Documents/MATLAB/Correia_et_al_(2013)/Fieldhouse_adapted')
cd '/Users/afieldhouse/Documents/MATLAB/Correia_et_al_(2013)/Fieldhouse_adapted'


%% Run Liquidity Trap, Allowing for negative nominal interest rates
dynare Fieldhouse_sticPflexW_adj0.mod noclearall

M0_BA = [Z, dI, dn, PI-1, i, dy, dc, deta - beta, i - PI+1];
M0_BP = [dg g dy y];
M0_W = [c, n, deta];

%% Produce Figure 2 (PosNom_consumerPI.m)

dynare Fieldhouse_sticPflexW_adj0.mod noclearall

M0_BA_UFP = [Z, dI, dn, PI-1, i, dy, dc, deta - beta, i - PI+1];
M0_BP_UFP = [dg g dy y];
M0_W_UFP = [c, n, deta];

% Save Figure 9: Allocation Allowing Negative Nominal Interest Rates
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Figure 9: Allocation IRFs, Centralized Fiscal Federalism Without Zero Bound Constraint', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
saveas(gcf,'SYP_F9','psc2')

% Compare allocations

BZLB_Less_UFP = M0_BA - M0_BA_UFP;
BZLB_Rel_UFP = M0_BA ./ M0_BA_UFP-1;
UFP_Less_BZLB = M0_BA_UFP - M0_BA;
UFP_Rel_BZLB = M0_BA_UFP ./ M0_BA -1 ;

% Compare welfare

[W0,CW0, LTW0] = CalcWelfare(M0_W,gamma,sigma);
[W0_UFP,CW0_UFP, LTW0_UFP] = CalcWelfare(M0_W_UFP,gamma,sigma);

% Make sure i = Z in Fieldhouse_sticPflexW_adj0.mod

%function [Nomi] = PosNom_consumerPI(Z,dI,dn,PI,i,dy,dc,beta,deta,y,u_c,u_n, I, R, k,sigma_I,delta,tau_k,c,n,mu);
Par_tau_c = 0.096;
Par_tau_n = 0.36;
Par_tau_k = 0.24;

x = length(i);
Nomi = zeros(x,1);
nrep    = 20;
smpl    = 2+(1:nrep);
for o = 1:x
    if i(o)<=0
        Nomi(o)=0;
    else
        Nomi(o)=i(o);
    end
end


tau_c=ones(x,1)*Par_tau_c;
tau_n=ones(x,1)*Par_tau_n;
tau_k=ones(x,1)*Par_tau_k;
u=ones(x,1)*0;      %u() is sI_tilde

tau_c_alt=ones(x,1)*Par_tau_c;
tau_n_alt=ones(x,1)*Par_tau_n;

tau_k_tilde=ones(x,1)*0; 
tau_k_tilde_CBO=ones(x,1)*0; 
Pt = ones(x,1);     %Initialize Pt, the price level 
Pt(1) = u_c(1) / u_n(1) * (1-tau_n(1))/(1+tau_c(1))*w(1);   % Set date-zero price level based on labor supply equation
taxSS = zeros(x,1);
taxGDP = zeros(x,1);
tau_c_growth = ones(x,1);
tau_c_alt_growth = ones(x,1);
PI_tilde = ones(x,1);
PI_tilde_alt = ones(x,1);


for o = 1:x-1
    % Correia et al code version
    if i(o+1)<=0
        tau_c(o+1) = (1+tau_c(o))/(1+i(o+1))-1;
        tau_n(o+1) = 1-(1-Par_tau_n)/((1+Par_tau_c))*(1+tau_c(o+1));
    else
       tau_c(o+1)=tau_c(o);
        tau_n(o+1)=tau_n(o);
    end
    
    % Alternative version backed straight out of Euler equationand labor supply
    if i(o+1)<=0
        tau_c_alt(o+1) = (1+tau_c_alt(o))*(1+Nomi(o))*deta(o)*u_c(o+1)/u_c(o)*(1/PI(o))-1;
        tau_n_alt(o+1) = 1-(1-Par_tau_n)/((1+Par_tau_c))*(1+tau_c_alt(o+1));
    else
       tau_c_alt(o+1)=tau_c_alt(o);
        tau_n_alt(o+1)=tau_n_alt(o);
    end
    
    % Back out price level from inflation
    Pt(o+1) = Pt(o)*PI(o+1); 
     
    %Calculate investment subsidy 
    for q = x-1:-1:1
          u(q) = 1 - (deta(q)*(1+tau_c(q))*(1-sigma_I*(I(q)/k(q)-delta))/(u_c(q)))*(u_c(q+1)/(1+tau_c(q+1))*(R(q+1)-Par_tau_k*(R(q+1)-delta))+u_c(q+1)*((1-u(q+1))/((1+tau_c(q+1))*(1-sigma_I*(I(q+1)/k(q+1)-delta))))*(1-delta-sigma_I/2*(I(q+1)/k(q+1)-delta)^2+sigma_I*(I(q+1)/k(q+1)-delta)*(I(q+1)/k(q+1))));  
    end
    
    %Calculate alternative capital income tax schedule
    tau_k_tilde(o+1) = ((u_c(o)/u_c(o+1))*(1+tau_c(o))/(1+tau_c(o))*(1/deta(o)) - (1-delta) - R(o+1)/Pt(o+1))/(delta-R(o+1)/Pt(o+1));
    %Back out alternative capital income tax schedule from Correia code, setting u(q), u(q+1) = 0
    tau_k_tilde_CBO(o+1) = R(q+1)/(R(q+1)-delta)*(deta(q)*(1+tau_c(q))*(1-sigma_I*(I(q)/k(q)-delta))/(u_c(q)))*(u_c(q+1)/(1+tau_c(q+1))/(1 + u_c(q+1)*(1/((1+tau_c(q+1))*(1-sigma_I*(I(q+1)/k(q+1)-delta))))*(1-delta-sigma_I/2*(I(q+1)/k(q+1)-delta)^2+sigma_I*(I(q+1)/k(q+1)-delta)*(I(q+1)/k(q+1)))));
    
end

% Inflation growth with consumption tax path
for o = 2:x
    tau_c_growth(o) = (1+tau_c(o))/(1+tau_c(o-1))-1;
    tau_c_alt_growth(o) = (1+tau_c_alt(o))/(1+tau_c_alt(o-1))-1;
    PI_tilde(o) = PI(o)*tau_c_growth(o);
    PI_tilde_alt(o) = PI(o)*tau_c_alt_growth(o);
end


% More failed attempts to back out Pt(+1) for a more sensible tau_k_tilde

% for o = 1:x-1
%     
%     % Back out price level from inflation
%     %Pt(o+1) = Pt(o)*(tau_c_growth(o+1)+PI(o+1)-1); 
%     Pt(o+1) = Pt(o)*PI(o+1);
%     
%     %Calculate alternative capital income tax schedule
%     tau_k_tilde(o+1) = ((u_c(o)/u_c(o+1))*(1+tau_c(o))/(1+tau_c(o))*(1/beta) - (1-delta) - R(o+1)/Pt(o+1))/(delta-R(o+1)/Pt(o+1));
%     tau_k_tilde_alt(o+1) = ((u_c(o)/u_c(o+1))*(1/beta) - (1-delta) - R(o+1)/Pt(o+1))/(delta-R(o+1)/Pt(o+1));
%     
% end
% 
% Pt = ones(x,1);
% 
% syms BackedPt
% for i=1:length(tau_n)-1
%     eqn3 = u_c(i) /(1+tau_c(i))  == beta*u_c(i+1) /(1+tau_c(i+1))*(1 - delta + (1-tau_k(i))*R(i)/BackedPt + tau_k(i)*delta);
%     solx3 = solve(eqn3,BackedPt);
%     Pt(i) = solx3;
% end
% 
% PI_alt = ones(x,1);
% for i = 2:length(tau_n)
%     PI_alt(i) = Pt(i)/Pt(i-1)-1;
% end


% SS = Par_tau_c * c(1) + Par_tau_k * (R(1)-delta)* k(1) + Par_tau_n * n(1);
% taxSS = (tau_c(1:402).*c + tau_n(1:402).*n-u(1:402).*I + tau_k * (R-delta).* k)/SS;
% taxGDP = (tau_c(1:402).*c + tau_n(1:402).*n-u(1:402).*I + tau_k * (R-delta).* k)./y;

% Correia version of plot

% figure
% subplot(331);plot(4*Z(smpl),'LineWidth',2);title('Z')
% subplot(332);plot(dI(smpl),'LineWidth',2);title('Investment')
% subplot(333);plot(dn(smpl),'LineWidth',2);title('Hours')
% subplot(334);plot(4*(tau_c_growth(smpl)+(PI(smpl)-1)),'LineWidth',2);title('Price Inflation (w/ Consumption Taxes)')
% subplot(335);plot(4*(Nomi(smpl)),'LineWidth',2);title('Nominal Interest Rate')
% subplot(336);plot(dy(smpl),'LineWidth',2);title('Output')
% subplot(337);plot(dc(smpl),'LineWidth',2);title('Consumption')
% subplot(338);plot(4*(deta(smpl)-beta),'LineWidth',2);title('Discount Factor Shock')
% subplot(339);plot(4*(i(smpl)-(PI(smpl)-1)),'LineWidth',2);title('Real Interest (w/ Consumption Taxes)')

figure
subplot(331);plot(4*Z(smpl),'LineWidth',2);title('Z')
subplot(332);plot(dI(smpl),'LineWidth',2);title('Investment')
subplot(333);plot(dn(smpl),'LineWidth',2);title('Hours')
subplot(334);plot(4*((PI_tilde(smpl))),'LineWidth',2);title('Effective Price Inflation')
subplot(335);plot(4*(Nomi(smpl)),'LineWidth',2);title('Nominal Interest Rate')
subplot(336);plot(dy(smpl),'LineWidth',2);title('Output')
subplot(337);plot(dc(smpl),'LineWidth',2);title('Consumption')
subplot(338);plot(4*(deta(smpl)-beta),'LineWidth',2);title('Discount Factor Shock')
subplot(339);plot(4*(Nomi(smpl)-PI_tilde(smpl)),'LineWidth',2);title('Effective Real Interest')

% Save Figure 10: Allocation With Unconventional Fiscal Policy
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Figure 10: Allocation IRFs, Centralized Fiscal Federalism With Unconventional Fiscal Policy', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
saveas(gcf,'SYP_F10','psc2')

%% Produce Figure 3 (sticPflexW_ITC_adj1_taxes.m)

%Plot Unconventional Fiscal Policy Instruments

close all
figure
subplot(311);plot(tau_c(2:21),'LineWidth',2);title('Consumption Tax')
ylabel('Percent')
subplot(312);plot(u(2:21),'LineWidth',2);title('Investment Tax Credit')
ylabel('Percent')
subplot(313);plot(tau_n(2:21),'LineWidth',2);title('Labor Income Tax')
ylabel('Percent')
xlabel('Quarters')

% Save Figure 11: Allocation Allowing Negative Nominal Interest Rates
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Figure 11: Policy Instruments, Centralized Fiscal Federalism With Unconventional Fiscal Policy', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
saveas(gcf,'SYP_F11','psc2')

%% Back out implied Fed. vs S&L income tax
tau_sn = zeros(402,1);
tau_fn = zeros(402,1);
tau_fn_mid = zeros(402,1);

syms x
for i=1:length(tau_n)
    tn = tau_n(i);
    eqn1 = tn == (1-x)*0.33 + x;
    solx1 = solve(eqn1,x);
    tau_sn(i,:) = solx1;
    
    eqn2 = tn == (1-0.044)*x + 0.044;
    solx2 = solve(eqn2,x);
    tau_fn(i,:) = solx2;
    
    eqn3 = tn == (1-max(tau_sn(i),0))*x + max(tau_sn(i),0);
    solx3 = solve(eqn3,x);
    tau_fn_mid(i,:) = solx3;
end

tau_sn_bound = max(tau_sn,0);

% Save unconventional policy

UnconvPolicy = [tau_c, tau_n, tau_sn, tau_fn, tau_sn_bound];

% Plot comparison of conventional and unconventional fiscal policy

kill
clf

figure
subplot(411);plot(UnconvPolicy(2:21,1),'LineWidth',2);title('Unconventional Consumption Tax')
ylabel('Percent')
subplot(412);plot(M2a_BP(2:21,11),'Color','m','LineWidth',1.2);title('Conventional Consumption Tax')
hold on
ylabel('Percent')
ylim([0.05,0.15])
subplot(412);plot(M5a_BP(2:21,11),'Color','k','LineWidth',1.2);
subplot(412);plot(M1a_BP(2:21,11),'LineWidth',1.5);
legend('Endogenous TauC','Great Recession','Other','Location','northeast','Orientation','vertical')
hold off
subplot(413);plot(UnconvPolicy(2:21,3),'LineWidth',2);title('Unconventional Labor Income Tax')
hold on
ylabel('Percent')
ylim([-0.1,0.4])
subplot(413);plot(UnconvPolicy(2:21,4),'Color','r','LineWidth',2);
legend('Just Sub-federal','Just Federal','Location','southeast','Orientation','vertical')
hline = refline(0,0);
hline.Color = 'k';
hold off
subplot(414);plot(M3a_BP(2:21,10),'Color','m','LineWidth',1.2);title('Conventional Sub-federal Labor Income Tax')
ylabel('Percent')
ylim([0.03,0.09])
hold on
subplot(414);plot(M5a_BP(2:21,10)+tau_nsDisc(2:21),'Color','k','LineWidth',1.2);
subplot(414);plot(M1a_BP(2:21,10),'LineWidth',1.5);
legend('Endogenous TauNS','Great Recession','Other','Location','southeast','Orientation','vertical')
hold off

ylabel('Percent')
xlabel('Quarters')

% Save Figure 12: Allocation Allowing Negative Nominal Interest Rates
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Figure 12: Decentralzied Policy Instruments: Conventional vs. Unconventional Fiscal Policy', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
saveas(gcf,'SYP_F12','psc2')

% Test SLG budget 
NeededTT_a = zeros(402,1);
NeededTT_b = zeros(402,1);
NeededTT_c = zeros(402,1);
NeededGS = zeros(402,1);
for i = 1:length(tau_n)
    NeededTT_a(i) = g(i)/2 - (tau_sn(i))*w(i)*n(i) - (tau_c(i))*c(i);
    NeededTT_b(i) = g(i)/2 - (0.044)*w(i)*n(i) - (tau_c(i))*c(i);
    NeededGS(i) = max(NeededTT_b(i),0) + (0.044)*w(i)*n(i) + (tau_c(i))*c(i);
    NeededTT_c(i) = g(i)/2 - (tau_sn_bound(i))*w(i)*n(i) - (tau_c(i))*c(i);
end

feasTT_b = max(NeededTT_b, 0);
reqTau_sn = (g/2 - tau_c.*c - feasTT_b)./(w.*n);
reqTau_fn = (tau_n - reqTau_sn) ./ (1 - reqTau_sn);

dTf_unconv = feasTT_b./M1a_BP(:,8)-1;

% Plot comparison of conventional and unconventional fiscal policy

kill
clf

figure
subplot(411);plot(100*tau_c(2:21),'LineWidth',2);title('Unconventional Sub-federal Consumption Tax')
ylabel('Percent')
%ylim([0.05,0.2])
subplot(412);plot(100*reqTau_sn(2:21),'LineWidth',2);title('Unconventional Sub-federal Labor Income Tax')
ylabel('Percent')
hline = refline(0,0);
hline.Color = 'k';
%ylim([-0.1,0.4])
subplot(413);plot(100*reqTau_fn(2:21),'Color','r','LineWidth',2);title('Unconventional Federal Labor Income Tax')
ylabel('Percent')
subplot(414);plot(dTf_unconv(2:21),'Color','r','LineWidth',2);title('Unconvention Federal Transfers to States')
ylabel('Percent')
hline = refline(0,0);
hline.Color = 'k';
ylim([-1,0.5])
ylabel('Deviation from Steady State')
xlabel('Quarters')

% Save Figure 13: Coordinated Unconventional Fiscal Policy
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Figure 13: Coordinated Unconventional Fiscal Policy: Decentralized Fiscal Federalism', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
saveas(gcf,'SYP_F13','psc2')

save parameterfile tau_c tau_n tau_sn tau_fn


%% Test three-way split
%tau_l
(0.035-0.03)*M0_BP(1,4)/M0_W(1,1); 


%% Section 4.2.1 

% Rerun unconventional fiscal policy model with Great Recession Mix (Model 3.2.5)
% FSL for Federal Stackelberg Leader

% Run Liquidity Trap

dynare Fieldhouse_sticPflexW_adj5.mod noclearall

if length(tau_nf) == 1
    tau_nf = repmat(tau_nf,length(M5a_BA),1);
end
if length(tau_ns) == 1
    tau_ns = repmat(tau_ns,length(M5a_BA),1);
end
if length(tau_c) == 1
    tau_c = repmat(tau_c,length(M5a_BA),1);
end
if length(gs) == 1
    gs = repmat(gs,length(M45a_BA),1);
end

M5_BA = [Z, dI, dn, PI-1, i, dy, dc, deta - beta, i - PI+1];
M5_BP  = [dgf, gf, dgs, gs, dy, y, dTf, Tf, tau_nf, tau_ns, tau_c, TT];
M5_W = [c, n, deta];

% ReRun Liquidity Trap, Allowing for negative nominal interest rates

dynare Fieldhouse_sticPflexW_adj5.mod noclearall

if length(tau_nf) == 1
    tau_nf = repmat(tau_nf,length(M5a_BA),1);
end
if length(tau_ns) == 1
    tau_ns = repmat(tau_ns,length(M5a_BA),1);
end
if length(tau_c) == 1
    tau_c = repmat(tau_c,length(M5a_BA),1);
end
if length(gs) == 1
    gs = repmat(gs,length(M5a_BA),1);
end

M5_BA_UFP = [Z, dI, dn, PI-1, i, dy, dc, deta - beta, i - PI+1];
M5_BP_UFP = [dgf, gf, dgs + gsDisc./gs(1,1), gs + gsDisc, dy, y, dTf, Tf, tau_nf, tau_ns + tau_nsDisc, tau_c, TT];
M5_W_UFP = [c, n, deta];

% Save Figure 9: Allocation Allowing Negative Nominal Interest Rates
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Figure 9: Allocation IRFs, Centralized Fiscal Federalism Without Zero Bound Constraint', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
%saveas(gcf,'SYP_F9','psc2')

% Compare allocations

FSL_BZLB_Less_UFP = M5_BA - M5_BA_UFP;
FSL_BZLB_Rel_UFP = M5_BA ./ M5_BA_UFP-1;
FSL_UFP_Less_BZLB = M5_BA_UFP - M5_BA;
FSL_UFP_Rel_BZLB = M5_BA_UFP ./ M5_BA -1;

% Compare welfare

[W5,CW5, LTW5] = CalcWelfare(M5_W,gamma,sigma);
[W5_UFP,CW5_UFP, LTW5_UFP] = CalcWelfare(M5_W_UFP,gamma,sigma);

% Make sure i = Z in Fieldhouse_sticPflexW_adj0.mod

%function [Nomi] = PosNom_consumerPI(Z,dI,dn,PI,i,dy,dc,beta,deta,y,u_c,u_n, I, R, k,sigma_I,delta,tau_k,c,n,mu);
M5_tau_c = M5_BP_UFP(:,11);
M5_tau_n = M5_BP_UFP(:,9).*(1-M5_BP_UFP(:,10)) + M5_BP_UFP(:,10);
M5_tau_k = 0.24;

x = length(i);
FSL_Nomi = zeros(x,1);

for o = 1:x
    if i(o)<=0
        FSL_Nomi(o)=0;
    else
        FSL_Nomi(o)=i(o);
    end
end


FSL_tau_c=ones(x,1)*Par_tau_c;
FSL_tau_n=ones(x,1)*Par_tau_n;
FSL_tau_k=ones(x,1)*Par_tau_k;
FSL_u=ones(x,1)*0;     

FSL_tau_c_alt=ones(x,1)*Par_tau_c;
FSL_tau_n_alt=ones(x,1)*Par_tau_n;

FSL_tau_k_tilde=ones(x,1)*0; 
FSL_tau_k_tilde_CBO=ones(x,1)*0; 
FSL_Pt = ones(x,1);     %Initialize Pt, the price level 
FSL_Pt(1) = u_c(1) / u_n(1) * (1-tau_n(1))/(1+tau_c(1))*w(1);   % Set date-zero price level based on labor supply equation
FSL_taxSS = zeros(x,1);
FSL_taxGDP = zeros(x,1);
FSL_tau_c_growth = ones(x,1);
FSL_tau_c_alt_growth = ones(x,1);
FSL_PI_tilde = ones(x,1);
FSL_PI_tilde_alt = ones(x,1);


for o = 1:x-1
    % Correia et al code version
    if i(o+1)<=0
        FSL_tau_c(o+1) = (1+FSL_tau_c(o))/(1+i(o+1))-1;
        FSL_tau_n(o+1) = 1-(1-M5_tau_n(o))/((1+M5_tau_c(o)))*(1+FSL_tau_c(o+1));
    else
        FSL_tau_c(o+1)=FSL_tau_c(o);
        FSL_tau_n(o+1)=FSL_tau_n(o);
    end
    
    % Alternative version backed straight out of Euler equationand labor supply
    if i(o+1)<=0
        FSL_tau_c_alt(o+1) = (1+FSL_tau_c_alt(o))*(1+FSL_Nomi(o))*deta(o)*u_c(o+1)/u_c(o)*(1/PI(o))-1;
        FSL_tau_n_alt(o+1) = 1-(1-M5_tau_n(o))/((1+M5_tau_c(o)))*(1+FSL_tau_c_alt(o+1));
    else
        FSL_tau_c_alt(o+1)=FSL_tau_c_alt(o);
        FSL_tau_n_alt(o+1)=FSL_tau_n_alt(o);
    end
    
    % Back out price level from inflation
    FSL_Pt(o+1) = FSL_Pt(o)*PI(o+1); 
     
    %Calculate investment subsidy 
    for q = x-1:-1:1
          FSL_u(q) = 1 - (deta(q)*(1+tau_c(q))*(1-sigma_I*(I(q)/k(q)-delta))/(u_c(q)))*(u_c(q+1)/(1+tau_c(q+1))*(R(q+1)-Par_tau_k*(R(q+1)-delta))+u_c(q+1)*((1-u(q+1))/((1+tau_c(q+1))*(1-sigma_I*(I(q+1)/k(q+1)-delta))))*(1-delta-sigma_I/2*(I(q+1)/k(q+1)-delta)^2+sigma_I*(I(q+1)/k(q+1)-delta)*(I(q+1)/k(q+1))));  
    end
    
    %Calculate alternative capital income tax schedule
    FSL_tau_k_tilde(o+1) = ((u_c(o)/u_c(o+1))*(1+FSL_tau_c(o))/(1+FSL_tau_c(o))*(1/deta(o)) - (1-delta) - R(o+1)/Pt(o+1))/(delta-R(o+1)/Pt(o+1));
    %Back out alternative capital income tax schedule from Correia code, setting u(q), u(q+1) = 0
    FSL_tau_k_tilde_CBO(o+1) = R(q+1)/(R(q+1)-delta)*(deta(q)*(1+FSL_tau_c(q))*(1-sigma_I*(I(q)/k(q)-delta))/(u_c(q)))*(u_c(q+1)/(1+FSL_tau_c(q+1))/(1 + u_c(q+1)*(1/((1+FSL_tau_c(q+1))*(1-sigma_I*(I(q+1)/k(q+1)-delta))))*(1-delta-sigma_I/2*(I(q+1)/k(q+1)-delta)^2+sigma_I*(I(q+1)/k(q+1)-delta)*(I(q+1)/k(q+1)))));
    
end

% Inflation growth with consumption tax path
for o = 2:x
    FSL_tau_c_growth(o) = (1+FSL_tau_c(o))/(1+FSL_tau_c(o-1))-1;
    FSL_tau_c_alt_growth(o) = (1+FSL_tau_c_alt(o))/(1+FSL_tau_c_alt(o-1))-1;
    FSL_PI_tilde(o) = PI(o)*FSL_tau_c_growth(o);
    FSL_PI_tilde_alt(o) = PI(o)*FSL_tau_c_alt_growth(o);
end

%% Feasability? 

% Figure 14: Desired Tax Rates for UFP with Federal Stackelberg Leader
clf
plot(FSL_tau_c(1:21),'--b')
hold on
title('Figure 14: Consumption Taxes')
xlabel('Quarters')
ylabel('Percent')
ylim([0.05,0.2])
plot(M5_BP(1:21,11),'g')
plot(M5_tau_c(1:21),'r')
legend('Unconventional Fiscal Policy Target','Conventional Policy','Conventional Policy with Unconventional Policy Allocation','Location','southeast','Orientation','vertical')
hold off
saveas(gcf,'SYP_F14','psc2')

%Feasible Fiscal Coersion: Bi-level
FSL_SLGdef= -Tf;
FSL_gsDiscA = psi*FSL_SLGdef;
FSL_tau_nsDiscA  = -nu*(1-psi)*FSL_SLGdef./(w.*n);
FSL_tau_c_FeasA  = -(1-nu)*(1-psi)*FSL_SLGdef./(c);

%Feasible Fiscal Coersion: Bi-level
psiB = 0.13;

FSL_gsDiscB = (1-psiB)*psi*FSL_SLGdef;
FSL_tau_nsDiscB  = -nu*(1-psi)*FSL_SLGdef./(w.*n);
FSL_tau_c_FeasB  = -((1-nu)*(1-psi)*FSL_SLGdef+psiB*psi*FSL_SLGdef)./(c); 


%Feasible Fiscal Coersion: Just Consumption Taxes
FSL_tau_c_FeasC  = -FSL_SLGdef./(c);
shortfall = FSL_tau_c - (tau_c + FSL_tau_c_FeasC);





