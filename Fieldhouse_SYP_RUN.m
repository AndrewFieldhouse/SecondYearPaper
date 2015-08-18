%   Fiscal Policy Under Decentralized Fiscal Federalism: 
%   Complications for Liquidity Trap Management
%   
%   Second Year Paper: Extension of Correia et. al (2013) with decentralized U.S. fiscal federalism

%   Andrew Fieldhouse
%   Cornell University
%   June 30, 2015

%%   Code for Section III: Centralized vs. Decentralized Federalism: Multipliers Above and Below ZLB
%
%    Note:  Code for Section IV is contained in UnconventionalFP.m
%           Code for Section V (Welfare Analysis) is included in this file
%
%    Other required files: 
%           ComputeWelfare.m
%           Dynare 

clear all
clear cfc

addpath(genpath('/Applications/Dynare/4.4.3/matlab'))
addpath('/Users/afieldhouse/Documents/MATLAB/Correia_et_al_(2013)/Fieldhouse_adapted')
cd '/Users/afieldhouse/Documents/MATLAB/Correia_et_al_(2013)/Fieldhouse_adapted'

%  Naming conventions: 'M[X]' for model, 'a' for baseline liquidity trap, 'b' for additional government spending shock, 'c' for F-to-SLG transfer shock

%% M0: Baseline Model of Decentralized Federalism

% M0a) Run baseline liquidity trap with preference shock, centralized government
dynare Fieldhouse_sticPflexW_adj0.mod noclearall

% annotation('textbox', [0 0.9 1 0.1], ...
%     'String', 'Figure 1: Allocation IRFs, Benchmark Model of Centralized Fiscal Federalism', ...
%     'EdgeColor', 'none', ...
%     'HorizontalAlignment', 'center')
% saveas(gcf,'SYP_F1','psc2')

% Save M0a benchmark path w/ full centralization, liquidity trap
M0a_BA = [Z, dI, dn, PI-1, i, dy, dc, deta - beta, i - PI+1];
M0a_BP = [dg g dy y];
M0a_W = [c, n, deta];

% Test GBC

imp_trans_M0a = tau_n*w(1)*n(1) + tau_c*c(1) + tau_k*k(1)*(R(1) - delta) - g(1);
profits = y(1) - R(1)*k(1) - w(1)*n(1);
ITpercY = imp_trans_M0a/y(1);
tau_d = 0.272;
imp_trans_M0a_alt = tau_n*w(1)*n(1) + tau_c*c(1) + tau_k*k(1)*(R(1) - delta) + profits*tau_d - g(1);
ITpercY_alt = imp_trans_M0a_alt/y(1);

% Figure 1: Plot Benchmark

x = 0:1:20;
clf;
subplot(3,3,1)
plot(x,400*M0a_BA(1:21,1),'b','linewidth',2)
ylabel('Percent (Annualized)')
hold on
title('Z')
subplot(3,3,2)
plot(x,100*M0a_BA(1:21,2),'b','linewidth',2)
title('Investment')
ylabel('Percentage Deviation')
subplot(3,3,3)
plot(x,100*M0a_BA(1:21,3),'b','linewidth',2)
title('Hours')
ylabel('Percentage Deviation')
subplot(3,3,4)
plot(x,400*M0a_BA(1:21,4),'b','linewidth',2)
title('Inflation')
ylabel('Percent (Annualized)')
subplot(3,3,5)
plot(x,400*M0a_BA(1:21,5),'b','linewidth',2)
title('Nominal Interest Rate')
ylabel('Percent (Annualized)')
subplot(3,3,6)
plot(x,100*M0a_BA(1:21,6),'b','linewidth',2)
title('Output')
ylabel('Percentage Deviation')
subplot(3,3,7)
plot(x,100*M0a_BA(1:21,7),'b','linewidth',2)
title('Consumption')
xlabel('Quarters')
ylabel('Percentage Deviation')
subplot(3,3,8)
plot(x,400*M0a_BA(1:21,8),'b','linewidth',2)
xlabel('Quarters')
ylabel('Percent (Annualized)')
title('Discount Factor Shock')
subplot(3,3,9)
plot(x,400*M0a_BA(1:21,9),'b','linewidth',2)
title('Real Interest Rate')
xlabel('Quarters')
ylabel('Percent (Annualized)')
hold off
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Figure 1: Liquidity Trap Allocation IRFs, Benchmark Model of Centralized Fiscal Federalism', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
saveas(gcf,'SYP_F1','psc2')


% M0b) Run liquidity simulation trap with additional government spending shock shocks

% Make sure exG is toggled on in Fieldhouse_sticPflexW_adj0.mod (Line 145-147)

dynare Fieldhouse_sticPflexW_adj0.mod noclearall

% Save M0b benchmark path w/ full centralization + govt spending shock

M0b_BA = [Z, dI, dn, PI-1, i, dy, dc, deta - beta, i - PI+1];
M0b_BP = [dg g dy y];
M0b_W = [c, n, deta];
M0a_W_AZLB = [repmat(c(1),402,1), repmat(n(1),402,1), repmat(deta(1),402,1)];
M0a_W_AZLB_altA = [repmat(c(1),402,1), repmat(n(1),402,1), repmat(deta(1),402,1), repmat(M0b_BP(1,2),402,1)];

% Check response of output when G is shocked
M0b_Less_M0a = M0b_BA-M0a_BA;
dM0b_BP = M0b_BP(:,4)./M0a_BP(:,4)-1;

%Calculate baseline government spending multiplier in liquidity trap ($1.28)
mG_LT_0b = (M0b_BP(:,4)-M0a_BP(:,4))./(M0b_BP(:,2)-M0a_BP(:,2));

%Calculate baseline government spending multiplier above ZLB ($0.96)
mG_0b = (M0b_BP(:,4)-M0b_BP(1,4))./(M0b_BP(:,2)-M0b_BP(1,2));



%% 1 Baseline Model of Decentralized Federalism: Endogenous G^S

% 1a) Run baseline liquidity trap with preference shock, decentralized government
dynare Fieldhouse_sticPflexW_adj1.mod noclearall

% annotation('textbox', [0 0.9 1 0.1], ...
%     'String', 'Figure 1c: Allocation IRFs with Centralized Fiscal Federalism, Endogenous GS', ...
%     'EdgeColor', 'none', ...
%     'HorizontalAlignment', 'center')
% saveas(gcf,'SYP_M1','psc2')

% Save M1a path w/ decentralization, endogenous GS

M1a_BA = [Z, dI, dn, PI-1, i, dy, dc, deta - beta, i - PI+1];
M1a_W = [c, n, deta];

% Make static policy instruments conform in length to endogenous variables

if length(tau_nf) == 1
    tau_nf = repmat(tau_nf,length(M1a_BA),1);
end
if length(tau_ns) == 1
    tau_ns = repmat(tau_ns,length(M1a_BA),1);
end
if length(tau_c) == 1
    tau_c = repmat(tau_c,length(M1a_BA),1);
end

M1a_BP  = [dgf, gf, dgs, gs, dy, y, dTf, Tf, tau_nf, tau_ns, tau_c, TT];

% Check relative allocation and policy responses 
M1a_Rel_M0a = M1a_BA./M0a_BA-1;
M1a_Less_M0a = M1a_BA-M0a_BA;

% Endogenous Policy Response Rel. Benchmark
dPR1a = [M1a_BP(:,1) + M1a_BP(:,3) - M0a_BP(:,1); M1a_BP(:,2) + M1a_BP(:,4) - M0a_BP(:,2)];

% Test GBC 
% Fed

F_imp_trans_M1a = tau_nf*(1-tau_ns(1))*w(1)*n(1)  + tau_k*k(1)*(R(1) - delta) - gf(1) - TT(1);
profits_M1a = y(1) - R(1)*k(1) - w(1)*n(1);
F_ITpercY_M1a = F_imp_trans_M1a/y(1);
tau_d = 0.272;
F_imp_trans_M1a_alt = tau_nf*(1-tau_ns(1))*w(1)*n(1) + tau_k*k(1)*(R(1) - delta) + profits*tau_d - gf(1) - TT(1);
F_ITpercY_alt = F_imp_trans_M1a_alt/y(1);

%subfed

SF_imp_trans_M0a = tau_ns(1)*w(1)*n(1) + tau_c*c(1) + TT(1) - gs(1);
SF_ITpercY = SF_imp_trans_M0a/y(1);


xaxbar = zeros(21,1);

% Figure 2: Plot comparison of M0a, M1a

clf;
subplot(3,3,1)
plot(x,400*M1a_Less_M0a(1:21,1),'r','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
ylabel('Percentage Points')
title('Monetary Policy Without ZLB')
hold off
subplot(3,3,2)
plot(x,100*M1a_Less_M0a(1:21,2),'r','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
title('Investment')
hold off
subplot(3,3,3)
plot(x,100*M1a_Less_M0a(1:21,3),'r','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
title('Hours')
hold off
subplot(3,3,4)
plot(x,400*M1a_Less_M0a(1:21,4),'r','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
title('Inflation')
ylabel('Percentage Points')
hold off
subplot(3,3,5)
plot(x,400*M1a_Less_M0a(1:21,5),'r','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
title('Nominal Interest Rate')
hold off
subplot(3,3,6)
plot(x,100*M1a_Less_M0a(1:21,6),'r','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
title('Output')
hold off
subplot(3,3,7)
plot(x,100*M1a_Less_M0a(1:21,7),'r','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
title('Consumption')
xlabel('Quarters')
ylabel('Percentage Points')
hold off
subplot(3,3,8)
plot(x,400*M1a_Less_M0a(1:21,9),'r','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
title('Real Interest Rate')
xlabel('Quarters')
hold off
subplot(3,3,9)
plot(x,100*dPR1a(1:21),'r','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
title('Sub-federal Spending')
xlabel('Quarters')
ylabel('Percentage Deviation')
hold off
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Figure 2: Decentralized Fiscal Federalism IRFs Less Centralized (Model 3.2.1: Endogenous GS)', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
saveas(gcf,'SYP_F2','psc2')


% 1b) Run decentralized liquidity simulation trap with additional federal spending shocks

% Make sure exGf is toggled on in Fieldhouse_sticPflexW_adj1.mod
% (lines 184-186)

dynare Fieldhouse_sticPflexW_adj1.mod noclearall

% Save M1b path w/ decentralization, endogenous GS + GF shock
M1b_BA = [Z, dI, dn, PI-1, i, dy, dc, deta - beta, i - PI+1];
M1b_W = [c, n, deta];
M1a_W_AZLB = [repmat(c(1),402,1), repmat(n(1),402,1), repmat(deta(1),402,1)];

if length(tau_nf) == 1
    tau_nf = repmat(tau_nf,length(M1b_BA),1);
end
if length(tau_ns) == 1
    tau_ns = repmat(tau_ns,length(M1b_BA),1);
end
if length(tau_c) == 1
    tau_c = repmat(tau_c,length(M1b_BA),1);
end

M1b_BP  = [dgf, gf, dgs, gs, dy, y, dTf, Tf, tau_nf, tau_ns, tau_c, TT];
M1a_W_AZLB_altA = [M1a_W_AZLB, repmat(M1b_BP(1,2),402,1) + repmat(M1b_BP(1,4),402,1)];
M1a_W_AZLB_altB = [M1a_W_AZLB, repmat(M1b_BP(1,2),402,1) , repmat(M1b_BP(1,4),402,1)];
dM1b_BP = M1b_BP-M1a_BP;

% Check relative output and allocation responses 
M1b_Rel_M1a = M1b_BA./M1a_BA-1;
M1b_Less_M1a = M1b_BA-M1a_BA;

% Check level response of output
y_level_1b = M1b_BP(:,6)./M1a_BP(:,4)-1;

%Calculate baseline government spending multiplier in liquidity trap ($1.348)

mG_LT_1b = (M1b_BP(:,6)-M1a_BP(:,6))./(M1b_BP(:,2)-M1a_BP(:,2));
mGdelta1b = mG_LT_1b - mG_LT_0b;
mGdeltaSR1b = mean(mGdelta1b(2:21));

%Calculate baseline government spending multiplier above ZLB ($0.999)
mG_1b = (M1b_BP(:,6)-M1b_BP(1,6))./(M1b_BP(:,2)-M1b_BP(1,2));

% 1c) Run decentralized liquidity simulation trap with additional federal transfers shocks

% Make sure exTf is toggled on in Fieldhouse_sticPflexW_adj1.mod
% (lines 184-186)

dynare Fieldhouse_sticPflexW_adj1.mod noclearall

% Save M1c path w/ decentralization, endogenous GS + Tf shock
M1c_BA = [Z, dI, dn, PI-1, i, dy, dc, deta - beta, i - PI+1];
M1c_W = [c, n, deta];

if length(tau_nf) == 1
    tau_nf = repmat(tau_nf,length(M1c_BA),1);
end
if length(tau_ns) == 1
    tau_ns = repmat(tau_ns,length(M1c_BA),1);
end
if length(tau_c) == 1
    tau_c = repmat(tau_c,length(M1c_BA),1);
end

M1c_BP  = [dgf, gf, dgs, gs, dy, y, dTf, Tf, tau_nf, tau_ns, tau_c, TT];
dM1c_BP = M1c_BP-M1a_BP;

% Check relative output and allocation responses 
M1c_Rel_M1a = M1c_BA./M1a_BA-1;
M1c_Less_M1a = M1c_BA-M1a_BA;

% Check level response of output
y_level_1c = M1c_BP(:,6)./M1a_BP(:,6)-1;

%Calculate traditional government spending multiplier in liquidity trap ($1.354)

mG_LT_1c = (M1c_BP(:,6)-M1a_BP(:,6))./(M1c_BP(:,8)-M1a_BP(:,8));
mGdelta1c = mG_LT_1c - mG_LT_0b;
mGdeltaSR1c = mean(mGdelta1c(2:21));

%Calculate baseline government spending multiplier above ZLB ($0.998)
mG_1c = (M1c_BP(:,6)-M1c_BP(1,6))./(M1c_BP(:,8)-M1c_BP(1,8));



% 1d) Run decentralized liquidity simulation trap with additional subfederal balanced budget spending shock

% Make sure exGs is toggled on in Fieldhouse_sticPflexW_adj1.mod
% (lines 195-197)

dynare Fieldhouse_sticPflexW_adj1.mod noclearall

% Save M1c path w/ decentralization, endogenous GS + Tf shock
M1d_BA = [Z, dI, dn, PI-1, i, dy, dc, deta - beta, i - PI+1];
M1d_W = [c, n, deta];

if length(tau_nf) == 1
    tau_nf = repmat(tau_nf,length(M1d_BA),1);
end
if length(tau_ns) == 1
    tau_ns = repmat(tau_ns,length(M1d_BA),1);
end
if length(tau_c) == 1
    tau_c = repmat(tau_c,length(M1d_BA),1);
end

M1d_BP  = [dgf, gf, dgs, gs, dy, y, dTf, Tf, tau_nf, tau_ns, tau_c, TT];
dM1d_BP = M1d_BP-M1a_BP;

% Check relative output and allocation responses 
M1d_Rel_M1a = M1d_BA./M1a_BA-1;
M1d_Less_M1a = M1d_BA-M1a_BA;

% Calculate sub-federal balanced budget govt spending multipler
% in liquidity trap ($1.354):

mG_LT_1d = (M1d_BP(:,6)-M1a_BP(:,6))./(M1d_BP(:,4)-M1a_BP(:,4));   % Trivially zero / undefined save first period of $1.05
mGdelta1d = mG_LT_1d - mG_LT_0b;
mGdeltaSR1d = mean(mGdelta1d(2:21));

%Calculate sub-federal balanced budget govt spending multipler above ZLB ($0.958)
mG_1d = (M1d_BP(:,6)-M1d_BP(1,6))./(M1d_BP(:,4)-M1d_BP(1,4));

%% 2) Baseline Model of Decentralized Federalism: Endogenous Tau_c

% 2a) Run baseline liquidity trap with preference shock, decentralized government
dynare Fieldhouse_sticPflexW_adj2.mod noclearall

M2a_BA = [Z, dI, dn, PI-1, i, dy, dc, deta - beta, i - PI+1];
M2a_W = [c, n, deta];

% Save benchmark path w/ decentralization
if length(tau_nf) == 1
    tau_nf = repmat(tau_nf,length(M2a_BA),1);
end
if length(tau_ns) == 1
    tau_ns = repmat(tau_ns,length(M2a_BA),1);
end
if length(tau_c) == 1
    tau_c = repmat(tau_c,length(M2a_BA),1);
end
if length(gs) == 1
    gs = repmat(gs,length(M2a_BA),1);
end

M2a_BP  = [dgf, gf, dgs, gs, dy, y, dTf, Tf, tau_nf, tau_ns, tau_c, TT];

% Check relative output and allocation responses rel. benchmark
M2a_Rel_M0a = M2a_BA./M0a_BA-1;
M2a_Less_M0a = M2a_BA-M0a_BA;

% Check level response of output rel. benchmark
y_level_2a_0a = M2a_BP(:,6)./M0a_BP(:,4)-1;
y_level_2a_1a = M2a_BP(:,6)./M1a_BP(:,6)-1;

% Check relative output and allocation responses rel. benchmark
M2a_Rel_M1a = M2a_BA./M1a_BA-1;
M2a_Less_M1a = M2a_BA-M1a_BA;

% Endogenous Policy Response Rel. Benchmark
dPR2a = [M2a_BP(:,11) - 0.096];

% Figure 3: Plot comparison of M0a, M2a
clf;
subplot(3,3,1)
plot(x,400*M2a_Less_M0a(1:21,1),'g','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
ylabel('Percentage Points')
title('Monetary Policy Without ZLB')
hold off
subplot(3,3,2)
plot(x,100*M2a_Less_M0a(1:21,2),'g','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
title('Investment')
hold off
subplot(3,3,3)
plot(x,100*M2a_Less_M0a(1:21,3),'g','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
title('Hours')
hold off
subplot(3,3,4)
plot(x,400*M2a_Less_M0a(1:21,4),'g','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
title('Inflation')
ylabel('Percentage Points')
hold off
subplot(3,3,5)
plot(x,400*M2a_Less_M0a(1:21,5),'g','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
title('Nominal Interest Rate')
hold off
subplot(3,3,6)
plot(x,100*M2a_Less_M0a(1:21,6),'g','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
title('Output')
hold off
subplot(3,3,7)
plot(x,100*M2a_Less_M0a(1:21,7),'g','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
title('Consumption')
xlabel('Quarters')
ylabel('Percentage Points')
hold off
subplot(3,3,8)
plot(x,400*M2a_Less_M0a(1:21,9),'g','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
title('Real Interest Rate')
xlabel('Quarters')
hold off
subplot(3,3,9)
plot(x,100*dPR2a(1:21),'g','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
title('Consumption Taxes')
xlabel('Quarters')
hold off
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Figure 3: Decentralized Fiscal Federalism IRFs Less Centralized (Model 3.2.2: Endogenous TauC)', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
saveas(gcf,'SYP_F3','psc2')


% 2b) Run decentralized liquidity simulation trap with additional federal spending shocks

% Make sure exGf is toggled on in Fieldhouse_sticPflexW_adj2.mod
% (lines 184-186)

dynare Fieldhouse_sticPflexW_adj2.mod noclearall

% Save benchmark path w/ decentralization
M2b_BA = [Z, dI, dn, PI-1, i, dy, dc, deta - beta, i - PI+1];
M2b_W = [c, n, deta];
M2a_W_AZLB = [repmat(c(1),402,1), repmat(n(1),402,1), repmat(deta(1),402,1)];

if length(tau_nf) == 1
    tau_nf = repmat(tau_nf,length(M2b_BA),1);
end
if length(tau_ns) == 1
    tau_ns = repmat(tau_ns,length(M2b_BA),1);
end
if length(tau_c) == 1
    tau_c = repmat(tau_c,length(M2b_BA),1);
end
if length(gs) == 1
    gs = repmat(gs,length(M2a_BA),1);
end

M2b_BP  = [dgf, gf, dgs, gs, dy, y, dTf, Tf, tau_nf, tau_ns, tau_c, TT];
dM2b_BP = M2b_BP-M2a_BP;
M2a_W_AZLB_altA = [M2a_W_AZLB, repmat(M2b_BP(1,2),402,1) + repmat(M2b_BP(1,4),402,1)];
M2a_W_AZLB_altB = [M2a_W_AZLB, repmat(M2b_BP(1,2),402,1) , repmat(M2b_BP(1,4),402,1)];

% Check relative output and allocation responses 
M2b_Relative_M0b = M2b_BA./M0b_BA-1;
M2b_Less_M0b = M2b_BA-M0b_BA;

M2b_Less_M2a = M2b_BA-M2a_BA;

% Calculate traditional federal govt spending multipler in liquidity trap ($1.326):

mG_LT_2b = (M2b_BP(:,6)-M2a_BP(:,6))./(M2b_BP(:,2)-M2a_BP(:,2));
mGdelta2b = mG_LT_2b - mG_LT_0b;
mGdeltaSR2b = mean(mGdelta2b(2:21));

%Calculate sub-federal balanced budget govt spending multipler above ZLB ($0.992)
mG_2b = (M2b_BP(:,6)-M2b_BP(1,6))./(M2b_BP(:,2)-M2b_BP(1,2));

% 2c) Run decentralized liquidity simulation trap with additional federal transfers shocks

% Make sure exTf is toggled on in Fieldhouse_sticPflexW_adj2.mod
% (lines 189-191)

dynare Fieldhouse_sticPflexW_adj2.mod noclearall

% Save M2c path w/ decentralization, endogenous GS + Tf shock
M2c_BA = [Z, dI, dn, PI-1, i, dy, dc, deta - beta, i - PI+1];
M2c_W = [c, n, deta];

if length(tau_nf) == 1
    tau_nf = repmat(tau_nf,length(M2c_BA),1);
end
if length(tau_ns) == 1
    tau_ns = repmat(tau_ns,length(M2c_BA),1);
end
if length(tau_c) == 1
    tau_c = repmat(tau_c,length(M2c_BA),1);
end

M2c_BP  = [dgf, gf, dgs, gs, dy, y, dTf, Tf, tau_nf, tau_ns, tau_c, TT];
dM2c_BP = M2c_BP-M2a_BP;

% Check relative output and allocation responses 
M2c_Rel_M2a = M2c_BA./M2a_BA-1;
M2c_Less_M2a = M2c_BA-M2a_BA;

% Calculate federal transfer spending multipler below ZLB ($0.916):

mG_LT_2c = (M2c_BP(:,6)-M2a_BP(:,6))./(M2c_BP(:,8)-M2a_BP(:,8));
mGdelta2c = mG_LT_2c - mG_LT_0b;
mGdeltaSR2c = mean(mGdelta2c(2:21));

%Calculate federal transfer spending multipler above ZLB ($0.824)
mG_2c = (M2c_BP(:,6)-M2c_BP(1,6))./(M2c_BP(:,8)-M2c_BP(1,8));

% 2d) Run decentralized liquidity simulation trap with additional subfederal balanced budget spending shock

% Make sure exGs is toggled on in Fieldhouse_sticPflexW_adj2.mod
% (lines 195-197)

dynare Fieldhouse_sticPflexW_adj2.mod noclearall

% Save M2d path w/ decentralization, endogenous GS + GS shock
M2d_BA = [Z, dI, dn, PI-1, i, dy, dc, deta - beta, i - PI+1];
M2d_W = [c, n, deta];

if length(tau_nf) == 1
    tau_nf = repmat(tau_nf,length(M2d_BA),1);
end
if length(tau_ns) == 1
    tau_ns = repmat(tau_ns,length(M2d_BA),1);
end
if length(tau_c) == 1
    tau_c = repmat(tau_c,length(M2d_BA),1);
end

M2d_BP  = [dgf, gf, dgs, gs, dy, y, dTf, Tf, tau_nf, tau_ns, tau_c, TT];
dM2d_BP = M2d_BP-M2a_BP;

% Check relative output and allocation responses 
M2d_Rel_M2a = M2d_BA./M2a_BA-1;
M2d_Less_M2a = M2d_BA-M2a_BA;

% Calculate sub-federal balanced budget govt spending multipler below ZLB ($0.916):

mG_LT_2d = (M2d_BP(:,6)-M2a_BP(:,6))./(M2d_BP(:,4)-M2a_BP(:,4));
mGdelta2d = mG_LT_2d - mG_LT_0b;
mGdeltaSR2d = mean(mGdelta2d(2:21));

%Calculate sub-federal balanced budget govt spending multipler above ZLB (N/A)
mG_2d = (M2d_BP(:,6)-M2d_BP(1,6))./(M2d_BP(:,4)-M2d_BP(1,4));


%% 3) Baseline Model of Decentralized Federalism: Endogenous Tau_n^s

% 3a) Run baseline liquidity trap with preference shock, decentralized
% government, endogenous SLG labor income tax response

dynare Fieldhouse_sticPflexW_adj3.mod noclearall

M3a_BA = [Z, dI, dn, PI-1, i, dy, dc, deta - beta, i - PI+1];
M3a_W = [c, n, deta];

% Save benchmark path w/ decentralization
if length(tau_nf) == 1
    tau_nf = repmat(tau_nf,length(M3a_BA),1);
end
if length(tau_ns) == 1
    tau_ns = repmat(tau_ns,length(M3a_BA),1);
end
if length(tau_c) == 1
    tau_c = repmat(tau_c,length(M3a_BA),1);
end
if length(gs) == 1
    gs = repmat(gs,length(M3a_BA),1);
end

M3a_BP  = [dgf, gf, dgs, gs, dy, y, dTf, Tf, tau_nf, tau_ns, tau_c, TT];

% Check relative output and allocation responses rel. M2, M1, M0
M3a_Rel_M0a = M3a_BA./M0a_BA-1;
M3a_Less_M0a = M3a_BA-M0a_BA;
M3a_Rel_M1a = M3a_BA./M1a_BA-1;
M3a_Less_M1a = M3a_BA-M1a_BA;
M3a_Rel_M2a = M3a_BA./M2a_BA-1;
M3a_Less_M2a = M3a_BA-M2a_BA;

% Check hours

nethours = sum(M3a_Less_M0a(2:7,3));

% Check level response of output rel. benchmark, M1, M2
y_level_3a_0 = M3a_BP(:,6)./M0a_BP(:,4)-1;
y_level_3a_1 = M3a_BP(:,6)./M1a_BP(:,6)-1;
y_level_3a_2 = M3a_BP(:,6)./M2a_BP(:,6)-1;

% Endogenous Policy Response Rel. Benchmark
dPR3a = [tau_nf .*(1-tau_ns) + tau_ns - 0.36];

% Figure 4: Plot comparison of M0a, M3a
clf;
subplot(3,3,1)
plot(x,400*M3a_Less_M0a(1:21,1),'m','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
ylabel('Percentage Points')
title('Monetary Policy Without ZLB')
hold off
subplot(3,3,2)
plot(x,100*M3a_Less_M0a(1:21,2),'m','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
title('Investment')
hold off
subplot(3,3,3)
plot(x,100*M3a_Less_M0a(1:21,3),'m','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
title('Hours')
hold off
subplot(3,3,4)
plot(x,400*M3a_Less_M0a(1:21,4),'m','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
title('Inflation')
ylabel('Percentage Points')
hold off
subplot(3,3,5)
plot(x,400*M3a_Less_M0a(1:21,5),'m','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
title('Nominal Interest Rate')
hold off
subplot(3,3,6)
plot(x,100*M3a_Less_M0a(1:21,6),'m','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
title('Output')
hold off
subplot(3,3,7)
plot(x,100*M3a_Less_M0a(1:21,7),'m','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
title('Consumption')
xlabel('Quarters')
ylabel('Percentage Points')
hold off
subplot(3,3,8)
plot(x,400*M3a_Less_M0a(1:21,9),'m','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
title('Real Interest Rate')
xlabel('Quarters')
hold off
subplot(3,3,9)
plot(x,100*dPR3a(1:21),'m','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
title('Sub-federal Labor Taxes')
xlabel('Quarters')
%ylabel('Percentage Deviation')
hold off
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Figure 4: Decentralized Fiscal Federalism IRFs Less Centralized (Model 3.2.3: Endogenous TauSFN)', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
saveas(gcf,'SYP_F4','psc2')

% 3b) Run decentralized liquidity simulation trap with additional federal spending shocks

% Make sure exGf is toggled on in Fieldhouse_sticPflexW_adj3.mod
% (lines 184-186)

dynare Fieldhouse_sticPflexW_adj3.mod noclearall

M3b_BA = [Z, dI, dn, PI-1, i, dy, dc, deta - beta, i - PI+1];
M3b_W = [c, n, deta];
M3a_W_AZLB = [repmat(c(1),402,1), repmat(n(1),402,1), repmat(deta(1),402,1)];

if length(tau_nf) == 1
    tau_nf = repmat(tau_nf,length(M3b_BA),1);
end
if length(tau_ns) == 1
    tau_ns = repmat(tau_ns,length(M3b_BA),1);
end
if length(tau_c) == 1
    tau_c = repmat(tau_c,length(M3b_BA),1);
end
if length(gs) == 1
    gs = repmat(gs,length(M3a_BA),1);
end

M3b_BP  = [dgf, gf, dgs, gs, dy, y, dTf, Tf, tau_nf, tau_ns, tau_c, TT];
M3a_W_AZLB_altA = [M3a_W_AZLB, repmat(M3b_BP(1,2),402,1) + repmat(M3b_BP(1,4),402,1)];
M3a_W_AZLB_altB = [M3a_W_AZLB, repmat(M3b_BP(1,2),402,1) , repmat(M3b_BP(1,4),402,1)];
dM3b_BP = M3b_BP-M3a_BP;

% Check relative output and allocation responses 
M3b_Relative_M3a = M3b_BA./M3a_BA-1;
M3b_Less_M3a = M3b_BA-M3a_BA;

% Calculate traditional federal govt spending multipler below ZLB ($0.916):

mG_LT_3b = (M3b_BP(:,6)-M3a_BP(:,6))./(M3b_BP(:,2)-M3a_BP(:,2));
mGdelta3b = mG_LT_3b - mG_LT_0b;
mGdeltaSR3b = mean(mGdelta3b(2:21));

%Calculate traditional federal govt spending multipler above ZLB (N/A)
mG_3b = (M3b_BP(:,6)-M3b_BP(1,6))./(M3b_BP(:,2)-M3b_BP(1,2));

% 3c) Run decentralized liquidity simulation trap with additional federal transfer shock

% Make sure exGf is toggled on in Fieldhouse_sticPflexW_adj3b.mod
% (lines 184-186)

dynare Fieldhouse_sticPflexW_adj3.mod noclearall

M3c_BA = [Z, dI, dn, PI-1, i, dy, dc, deta - beta, i - PI+1];
M3c_W = [c, n, deta];

% Save benchmark path w/ decentralization
if length(tau_nf) == 1
    tau_nf = repmat(tau_nf,length(M3c_BA),1);
end
if length(tau_ns) == 1
    tau_ns = repmat(tau_ns,length(M3c_BA),1);
end
if length(tau_c) == 1
    tau_c = repmat(tau_c,length(M3c_BA),1);
end
if length(gs) == 1
    gs = repmat(gs,length(M3a_BA),1);
end

M3c_BP  = [dgf, gf, dgs, gs, dy, y, dTf, Tf, tau_nf, tau_ns, tau_c, TT];
dM3c_BP = M3c_BP-M3a_BP;

% Check relative output and allocation responses 
M3c_Rel_M3a = M3c_BA./M3a_BA-1;
M3c_Less_M3a = M3c_BA-M3a_BA;

% Calculate federal transfer multipler below ZLB ($0.916):

mG_LT_3c = (M3c_BP(:,6)-M3a_BP(:,6))./(M3c_BP(:,8)-M3a_BP(:,8));
mGdelta3c = mG_LT_3c - mG_LT_0b;
mGdeltaSR3c = mean(mGdelta3c(2:21));

%Calculate federal transfer multipler above ZLB ($0.3567)
mG_3c = (M3c_BP(:,6)-M3c_BP(1,6))./(M3c_BP(:,8)-M3c_BP(1,8));

% 3d) Run decentralized liquidity simulation trap with additional subfederal balanced budget spending shock

% Make sure exGs is toggled on in Fieldhouse_sticPflexW_adj3.mod
% (lines 195-197)

dynare Fieldhouse_sticPflexW_adj3.mod noclearall

% Save M3d path w/ decentralization, endogenous GS + GS shock
M3d_BA = [Z, dI, dn, PI-1, i, dy, dc, deta - beta, i - PI+1];
M3d_W = [c, n, deta];

if length(tau_nf) == 1
    tau_nf = repmat(tau_nf,length(M3d_BA),1);
end
if length(tau_ns) == 1
    tau_ns = repmat(tau_ns,length(M3d_BA),1);
end
if length(tau_c) == 1
    tau_c = repmat(tau_c,length(M3d_BA),1);
end

M3d_BP  = [dgf, gf, dgs, gs, dy, y, dTf, Tf, tau_nf, tau_ns, tau_c, TT];
dM3d_BP = M3d_BP-M3a_BP;
M4a_W_AZLB_altA = [M4a_W_AZLB, repmat(M4b_BP(1,2),402,1) + repmat(M4b_BP(1,4),402,1)];
M4a_W_AZLB_altB = [M4a_W_AZLB, repmat(M4b_BP(1,2),402,1) , repmat(M4b_BP(1,4),402,1)];

% Check relative output and allocation responses 
M3d_Rel_M3a = M3d_BA./M3a_BA-1;
M3d_Less_M3a = M3d_BA-M3a_BA;

% Calculate sub-federal balanced budget govt spending multipler below ZLB ($0.916):

mG_LT_3d = (M3d_BP(:,6)-M3a_BP(:,6))./(M3d_BP(:,4)-M3a_BP(:,4));
mGdelta3d = mG_LT_3d - mG_LT_0b;
mGdeltaSR3d = mean(mGdelta3d(2:21));

%Calculate sub-federal balanced budget govt spending multipler above ZLB ($0.626)
mG_3d = (M3d_BP(:,6)-M3d_BP(1,6))./(M3d_BP(:,4)-M3d_BP(1,4));

%% 4) Model 4: Decentralized Federalism: Endogenous Federal Transfer Response

% 4a) Run baseline liquidity trap with preference shock, decentralized
% government, endogenous SLG labor income tax response

dynare Fieldhouse_sticPflexW_adj4.mod noclearall

M4a_BA = [Z, dI, dn, PI-1, i, dy, dc, deta - beta, i - PI+1];
M4a_W = [c, n, deta];

% Save benchmark path w/ decentralization
if length(tau_nf) == 1
    tau_nf = repmat(tau_nf,length(M4a_BA),1);
end
if length(tau_ns) == 1
    tau_ns = repmat(tau_ns,length(M4a_BA),1);
end
if length(tau_c) == 1
    tau_c = repmat(tau_c,length(M4a_BA),1);
end
if length(gs) == 1
    gs = repmat(gs,length(M4a_BA),1);
end

M4a_BP  = [dgf, gf, dgs, gs, dy, y, dTf, Tf, tau_nf, tau_ns, tau_c, TT];

dPR4a = M4a_BP(:,7) - M1a_BP(:,7);

% Check relative output and allocation responses rel. M2, M1, M0
M4a_Rel_M2a = M4a_BA./M2a_BA-1;
M4a_Less_M2a = M4a_BA-M2a_BA;
M4a_Rel_M1a = M4a_BA./M1a_BA-1;
M4a_Less_M1a = M4a_BA-M1a_BA;
M4a_Rel_M0a = M4a_BA./M0a_BA-1;
M4a_Less_M0a = M4a_BA-M0a_BA;

% Check level response of output rel. benchmark, M1, M2
y_level_4a_0 = M4a_BP(:,6)./M0a_BP(:,4)-1;
y_level_4a_1 = M4a_BP(:,6)./M1a_BP(:,6)-1;
y_level_4a_2 = M4a_BP(:,6)./M2a_BP(:,6)-1;

% % Figure 5: Plot comparison of M0a, M4a
% clf;
% subplot(3,3,1)
% plot(x,400*M4a_Less_M0a(1:21,1),'c','linewidth',1.5)
% hold on
% plot(x,xaxbar,'k--')
% ylabel('Percentage Points')
% title('Monetary Policy Without ZLB')
% hold off
% subplot(3,3,2)
% plot(x,100*M4a_Less_M0a(1:21,2),'c','linewidth',1.5)
% hold on
% plot(x,xaxbar,'k--')
% title('Investment')
% hold off
% subplot(3,3,3)
% plot(x,100*M4a_Less_M0a(1:21,3),'c','linewidth',1.5)
% hold on
% plot(x,xaxbar,'k--')
% title('Hours')
% hold off
% subplot(3,3,4)
% plot(x,400*M4a_Less_M0a(1:21,4),'c','linewidth',1.5)
% hold on
% plot(x,xaxbar,'k--')
% title('Inflation')
% ylabel('Percentage Points')
% hold off
% subplot(3,3,5)
% plot(x,400*M4a_Less_M0a(1:21,5),'c','linewidth',1.5)
% hold on
% plot(x,xaxbar,'k--')
% title('Nominal Interest Rate')
% hold off
% subplot(3,3,6)
% plot(x,100*M4a_Less_M0a(1:21,6),'c','linewidth',1.5)
% hold on
% plot(x,xaxbar,'k--')
% title('Output')
% hold off
% subplot(3,3,7)
% plot(x,100*M4a_Less_M0a(1:21,7),'c','linewidth',1.5)
% hold on
% plot(x,xaxbar,'k--')
% title('Consumption')
% xlabel('Quarters')
% ylabel('Percentage Points')
% hold off
% subplot(3,3,8)
% plot(x,400*M4a_Less_M0a(1:21,9),'c','linewidth',1.5)
% hold on
% plot(x,xaxbar,'k--')
% title('Real Interest Rate')
% xlabel('Quarters')
% hold off
% subplot(3,3,9)
% plot(x,100*dPR4a(1:21),'c','linewidth',1.5)
% hold on
% plot(x,xaxbar,'k--')
% title('Federal Transfers to States')
% xlabel('Quarters')
% ylabel('Percentage Deviation')
% hold off
% annotation('textbox', [0 0.9 1 0.1], ...
%     'String', 'Figure 5: Decentralized Fiscal Federalism IRFs Less Centralized (Model 3.2.4: Endogenous Federal Transfers)', ...
%     'EdgeColor', 'none', ...
%     'HorizontalAlignment', 'center')
% saveas(gcf,'SYP_F5','psc2')

% 4b) Run decentralized liquidity simulation trap with additional federal spending shocks

% Make sure exGf is toggled on in Fieldhouse_sticPflexW_adj4.mod
% (lines 184-186)

dynare Fieldhouse_sticPflexW_adj4.mod noclearall

M4b_BA = [Z, dI, dn, PI-1, i, dy, dc, deta - beta, i - PI+1];
M4b_W = [c, n, deta];
M4a_W_AZLB = [repmat(c(1),402,1), repmat(n(1),402,1), repmat(deta(1),402,1)];

if length(tau_nf) == 1
    tau_nf = repmat(tau_nf,length(M4b_BA),1);
end
if length(tau_ns) == 1
    tau_ns = repmat(tau_ns,length(M4b_BA),1);
end
if length(tau_c) == 1
    tau_c = repmat(tau_c,length(M4b_BA),1);
end
if length(gs) == 1
    gs = repmat(gs,length(M4a_BA),1);
end

M4b_BP  = [dgf, gf, dgs, gs, dy, y, dTf, Tf, tau_nf, tau_ns, tau_c, TT];
dM4b_BP = M4b_BP-M4a_BP;

% Check relative output and allocation responses 
M4b_Relative_M4a = M4b_BA./M4a_BA-1;
M4b_Less_M4a = M4b_BA-M4a_BA;

% Calculate traditional federal govt spending multipler below ZLB ($0.916):

mG_LT_4b = (M4b_BP(:,6)-M4a_BP(:,6))./(M4b_BP(:,2)-M4a_BP(:,2));
mGdelta4b = mG_LT_4b - mG_LT_0b;
mGdeltaSR4b = mean(mGdelta4b(2:21));

%Calculate traditional federal govt spending multipler above ZLB (N/A)
mG_4b = (M4b_BP(:,6)-M4b_BP(1,6))./(M4b_BP(:,2)-M4b_BP(1,2));


% 4c) Run decentralized liquidity simulation trap with additional federal transfer shock

% Make sure exGf is toggled on in Fieldhouse_sticPflexW_adj4.mod
% (lines 184-186)

dynare Fieldhouse_sticPflexW_adj4.mod noclearall

M4c_BA = [Z, dI, dn, PI-1, i, dy, dc, deta - beta, i - PI+1];
M4c_W = [c, n, deta];

% Save benchmark path w/ decentralization
if length(tau_nf) == 1
    tau_nf = repmat(tau_nf,length(M4c_BA),1);
end
if length(tau_ns) == 1
    tau_ns = repmat(tau_ns,length(M4c_BA),1);
end
if length(tau_c) == 1
    tau_c = repmat(tau_c,length(M4c_BA),1);
end
if length(gs) == 1
    gs = repmat(gs,length(M4a_BA),1);
end

M4c_BP  = [dgf, gf, dgs, gs, dy, y, dTf, Tf, tau_nf, tau_ns, tau_c, TT];
dM4c_BP = M4c_BP-M4a_BP;

% Check relative output and allocation responses 
M4c_Rel_M4a = M4c_BA./M4a_BA-1;
M4c_Less_M4a = M4c_BA-M4a_BA;

% Calculate federal transfer multipler below ZLB ($0.916):

mG_LT_4c = (M4c_BP(:,6)-M4a_BP(:,6))./(M4c_BP(:,8)-M4a_BP(:,8));
mGdelta4c = mG_LT_4c - mG_LT_0b;
mGdeltaSR4c = mean(mGdelta4c(2:21));

%Calculate federal transfer multipler above ZLB ($0.3567)
mG_4c = (M4c_BP(:,6)-M4c_BP(1,6))./(M4c_BP(:,8)-M4c_BP(1,8));

% 4d) Run decentralized liquidity simulation trap with additional subfederal balanced budget spending shock

% Make sure exGs is toggled on in Fieldhouse_sticPflexW_adj3.mod
% (lines 195-197)

dynare Fieldhouse_sticPflexW_adj4.mod noclearall

% Save M4d path w/ decentralization, endogenous GS + GS shock
M4d_BA = [Z, dI, dn, PI-1, i, dy, dc, deta - beta, i - PI+1];
M4d_W = [c, n, deta];

if length(tau_nf) == 1
    tau_nf = repmat(tau_nf,length(M4d_BA),1);
end
if length(tau_ns) == 1
    tau_ns = repmat(tau_ns,length(M4d_BA),1);
end
if length(tau_c) == 1
    tau_c = repmat(tau_c,length(M4d_BA),1);
end

M4d_BP  = [dgf, gf, dgs, gs, dy, y, dTf, Tf, tau_nf, tau_ns, tau_c, TT];
dM4d_BP = M4d_BP-M4a_BP;

% Check relative output and allocation responses 
M4d_Rel_M4a = M4d_BA./M4a_BA-1;
M4d_Less_M4a = M4d_BA-M4a_BA;

% Calculate sub-federal balanced budget govt spending multipler below ZLB ($0.916):

mG_LT_4d = (M4d_BP(:,6)-M4a_BP(:,6))./(M4d_BP(:,4)-M4a_BP(:,4));
mGdelta4d = mG_LT_4d - mG_LT_0b;
mGdeltaSR4d = mean(mGdelta4d(2:21));

%Calculate sub-federal balanced budget govt spending multipler above ZLB ($0.626)
mG_4d = (M4d_BP(:,6)-M4d_BP(1,6))./(M4d_BP(:,4)-M4d_BP(1,4));

%% 5) Model 3.2.5: Decentralized Federalism: Endogenous Federal Transfer Response

% 5a) Run baseline liquidity trap with preference shock, decentralized
% government, endogenous SLG labor income tax response

dynare Fieldhouse_sticPflexW_adj5.mod noclearall

M5a_BA = [Z, dI, dn, PI-1, i, dy, dc, deta - beta, i - PI+1];
M5a_W = [c, n, deta];

% Save benchmark path w/ decentralization
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

M5a_BP  = [dgf, gf, dgs - gsDisc./gs(1,1), gs - gsDisc, dy, y, dTf, Tf, tau_nf, tau_ns + tau_nsDisc, tau_c, TT];

% Check discretionary sub-federal policy response ratios
dLTax = n.*w.*(tau_nsDisc);
dCTax = c.*(tau_c-0.0968798);
dGs = gsDisc./(dLTax+ dCTax);
CLrevratio = dCTax./dLTax;

% Check policy instruments
ntauns = tau_ns + tau_nsDisc;
nettau = tau_nf.*(1 - ntauns) + ntauns;
nGS = gs - gsDisc;
dnGS = nGS./gs;

% Check relative output and allocation responses rel. benchmark
M5a_Rel_M0a = M5a_BA./M0a_BA-1;
M5a_Less_M0a = M5a_BA-M0a_BA;

% Check level response of output rel. benchmark
y_level_5a_0a = M5a_BP(:,6)./M0a_BP(:,4)-1;
y_level_5a_1a = M5a_BP(:,6)./M1a_BP(:,6)-1;

% Check relative output and allocation responses rel. benchmark
M5a_Rel_M1a = M5a_BA./M1a_BA-1;
M5a_Less_M1a = M5a_BA-M1a_BA;


dPR5a = [-gsDisc./gs(1,1), M5a_BP(:,11) - M5a_BP(1,11), M5a_BP(:,10) - M1a_BP(:,10)];


% Figure 5: Plot comparison of M0a, M5a
clf;
subplot(3,3,1)
plot(x,400*M5a_Less_M0a(1:21,1),'c','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
ylabel('Percentage Points')
title('Monetary Policy Without ZLB')
hold off
subplot(3,3,2)
plot(x,100*M5a_Less_M0a(1:21,2),'c','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
title('Investment')
hold off
subplot(3,3,3)
plot(x,100*M5a_Less_M0a(1:21,3),'c','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
title('Hours')
hold off
subplot(3,3,4)
plot(x,400*M5a_Less_M0a(1:21,4),'c','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
title('Inflation')
ylabel('Percentage Points')
hold off
subplot(3,3,5)
plot(x,400*M5a_Less_M0a(1:21,5),'c','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
title('Nominal Interest Rate')
hold off
subplot(3,3,6)
plot(x,100*M5a_Less_M0a(1:21,6),'c','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
title('Output')
hold off
subplot(3,3,7)
plot(x,100*M5a_Less_M0a(1:21,7),'c','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
title('Consumption')
xlabel('Quarters')
ylabel('Percentage Points')
hold off
subplot(3,3,8)
plot(x,400*M5a_Less_M0a(1:21,9),'c','linewidth',1.5)
hold on
plot(x,xaxbar,'k--')
title('Real Interest Rate')
xlabel('Quarters')
hold off
subplot(3,3,9)
plot(x,100*dPR5a(1:21,1),'r','linewidth',1.2)
ylim([-14,2])
hold on
plot(x,100*dPR5a(1:21,2),'g','linewidth',1.2)
plot(x,100*dPR5a(1:21,3),'m','linewidth',1.2)
plot(x,xaxbar,'k--')
title('Mixed Policy Responses')
xlabel('Quarters')
legend('Gsf','TauC','TauNsf','Location','southeast','Orientation','vertical')
%ylabel('Percentage Deviation')
hold off
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Figure 5: Decentralized Fiscal Federalism IRFs Less Centralized (Model 3.2.5: Great Recession Mix)', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
saveas(gcf,'SYP_F5','psc2')


% 5b) Run decentralized liquidity simulation trap with additional federal spending shocks

% Make sure exGf is toggled on in Fieldhouse_sticPflexW_adj5.mod
% (lines 184-186)

dynare Fieldhouse_sticPflexW_adj5.mod noclearall

M5b_BA = [Z, dI, dn, PI-1, i, dy, dc, deta - beta, i - PI+1];
M5b_W = [c, n, deta];
M5a_W_AZLB = [repmat(c(1),402,1), repmat(n(1),402,1), repmat(deta(1),402,1)];

if length(tau_nf) == 1
    tau_nf = repmat(tau_nf,length(M5b_BA),1);
end
if length(tau_ns) == 1
    tau_ns = repmat(tau_ns,length(M5b_BA),1);
end
if length(tau_c) == 1
    tau_c = repmat(tau_c,length(M5b_BA),1);
end
if length(gs) == 1
    gs = repmat(gs,length(M5a_BA),1);
end

M5b_BP  = [dgf, gf, dgs - gsDisc./gs(1,1), gs - gsDisc, dy, y, dTf, Tf, tau_nf, tau_ns + tau_nsDisc, tau_c, TT];
dM5b_BP = M5b_BP-M5a_BP;

% Check relative output and allocation responses 
M5b_Relative_M5a = M5b_BA./M5a_BA-1;
M5b_Less_M5a = M5b_BA-M5a_BA;

% Calculate traditional federal govt spending multipler below ZLB ($0.916):

mG_LT_5b = (M5b_BP(:,6)-M5a_BP(:,6))./(M5b_BP(:,2)-M5a_BP(:,2));
mGdelta5b = mG_LT_5b - mG_LT_0b;
mGdeltaSR5b = mean(mGdelta5b(2:21));

%Calculate traditional federal govt spending multipler above ZLB (0.964)
mG_5b = (M5b_BP(:,6)-M5b_BP(1,6))./(M5b_BP(:,2)-M5b_BP(1,2));


% 5c) Run decentralized liquidity simulation trap with additional federal transfer shock

% Make sure exGf is toggled on in Fieldhouse_sticPflexW_adj5.mod
% (lines 184-186)

dynare Fieldhouse_sticPflexW_adj5.mod noclearall

M5c_BA = [Z, dI, dn, PI-1, i, dy, dc, deta - beta, i - PI+1];
M5c_W = [c, n, deta];

% Save benchmark path w/ decentralization
if length(tau_nf) == 1
    tau_nf = repmat(tau_nf,length(M5c_BA),1);
end
if length(tau_ns) == 1
    tau_ns = repmat(tau_ns,length(M5c_BA),1);
end
if length(tau_c) == 1
    tau_c = repmat(tau_c,length(M5c_BA),1);
end
if length(gs) == 1
    gs = repmat(gs,length(M5a_BA),1);
end

M5c_BP  = [dgf, gf, dgs - gsDisc./gs(1,1), gs - gsDisc, dy, y, dTf, Tf, tau_nf, tau_ns + tau_nsDisc, tau_c, TT];
dM5c_BP = M5c_BP-M5a_BP;

% Check relative output and allocation responses 
M5c_Rel_M5a = M5c_BA./M5a_BA-1;
M5c_Less_M5a = M5c_BA-M5a_BA;

% Calculate federal transfer multipler below ZLB ($0.916):

mG_LT_5c = (M5c_BP(:,6)-M5a_BP(:,6))./(M5c_BP(:,8)-M5a_BP(:,8));
mGdelta5c = mG_LT_5c - mG_LT_0b;
mGdeltaSR5c = mean(mGdelta5c(2:21));

%Calculate federal transfer multipler above ZLB ($0.3567)
mG_5c = (M5c_BP(:,6)-M5c_BP(1,6))./(M5c_BP(:,8)-M5c_BP(1,8));

% 5d) Run decentralized liquidity simulation trap with additional subfederal balanced budget spending shock

% Make sure exGs is toggled on in Fieldhouse_sticPflexW_adj3.mod
% (lines 195-197)

dynare Fieldhouse_sticPflexW_adj5.mod noclearall

% Save M5d path w/ decentralization, endogenous GS + GS shock
M5d_BA = [Z, dI, dn, PI-1, i, dy, dc, deta - beta, i - PI+1];
M5d_W = [c, n, deta];

if length(tau_nf) == 1
    tau_nf = repmat(tau_nf,length(M5d_BA),1);
end
if length(tau_ns) == 1
    tau_ns = repmat(tau_ns,length(M5d_BA),1);
end
if length(tau_c) == 1
    tau_c = repmat(tau_c,length(M5d_BA),1);
end

M5d_BP  = [dgf, gf, dgs - gsDisc./gs(1,1), gs - gsDisc, dy, y, dTf, Tf, tau_nf, tau_ns + tau_nsDisc, tau_c, TT];
dM5d_BP = M5d_BP-M5a_BP;

% Check relative output and allocation responses 
M5d_Rel_M5a = M5d_BA./M5a_BA-1;
M5d_Less_M5a = M5d_BA-M5a_BA;

% Calculate sub-federal balanced budget govt spending multipler below ZLB ($):

mG_LT_5d = (M5d_BP(:,6)-M5a_BP(:,6))./(M5d_BP(:,4)-M5a_BP(:,4));
mGdelta5d = mG_LT_5d - mG_LT_0b;
mGdeltaSR5d = mean(mGdelta5d(2:21));

%Calculate sub-federal balanced budget govt spending multipler above ZLB ($)
mG_5d = (M5d_BP(:,6)-M5d_BP(1,6))./(M5d_BP(:,4)-M5d_BP(1,4));

%% 6) Store Multipliers and Compare

mG_0fill = zeros(length(mG_0b),1);

AboveZLBb = [mG_0b, mG_1b, mG_2b, mG_3b, mG_4b, mG_5b];
AboveZLBc = [mG_0fill, mG_1c, mG_2c, mG_3c, mG_4c, mG_5c];
AboveZLBd = [mG_0fill, mG_1d, mG_2d, mG_3d, mG_4d, mG_5d];

T2a = [AboveZLBb(2,:)', AboveZLBc(2,:)', AboveZLBd(2,:)'];
T2a(2,3) = 0;       

BelowZLBb = [mG_LT_0b, mG_LT_1b, mG_LT_2b, mG_LT_3b, mG_LT_4b, mG_LT_5b];
BelowZLBc = [mG_0fill, mG_LT_1c, mG_LT_2c, mG_LT_3c, mG_LT_4c, mG_LT_5c];
BelowZLBd = [mG_0fill, mG_LT_1d, mG_LT_2d, mG_LT_3d, mG_LT_4d, mG_LT_5d];

T2b = [BelowZLBb(2,:)', BelowZLBc(2,:)', BelowZLBd(2,:)'];
T2b(2,3) = 0;   

MultDelt = 100*[(BelowZLBb(2,:)./AboveZLBb(2,:)-1)', (BelowZLBc(2,:)./AboveZLBc(2,:)-1)', (BelowZLBd(2,:)./AboveZLBd(2,:)-1)'];

GRmixD = mG_LT_5b(2)/mG_LT_0b(2);




%% 7) Comparison of multipliers and output IRFs accross Benchmark and four models

xs = 1:1:8;
xsaxis = zeros(8,1);

% Figure 6: Fiscal Multipliers
clf
subplot(3,1,1)
plot(xs,mG_LT_0b(2:9),'--b','linewidth',1.3)
hold on
title('Federal Spending')
ylabel('dY/dGf')
ylim([1, 1.4])
plot(xs,mG_LT_1b(2:9),'r','linewidth',1.3)
plot(xs,mG_LT_2b(2:9),'g','linewidth',1.3)
plot(xs,mG_LT_3b(2:9),'m','linewidth',1.3)
plot(xs,mG_LT_4b(2:9),'y','linewidth',1.3)
plot(xs,mG_LT_5b(2:9),'--c','linewidth',1.3)
%legend('Centralized Model','Decentralized: Endogenous GS ','Decentralized: Endogenous TauC','Decentralized: Endogenous TauN','Location','northeast','Orientation','vertical')
subplot(3,1,2)
plot(xs,mG_LT_1c(2:9),'r','linewidth',1.3)
title('Federal Transfers')
hold on
ylabel('dY/dTF')
plot(xs,mG_LT_2c(2:9),'g','linewidth',1.3)
plot(xs,mG_LT_3c(2:9),'m','linewidth',1.3)
plot(xs,mG_LT_5c(2:9),'--c','linewidth',1.3)
plot(xs, xsaxis,'k')
hold off
subplot(3,1,3)
plot(xs,mG_LT_2d(2:9),'g','linewidth',1.3)
hold on
title('Sub-Federal Spending')
xlabel('Quarters')
ylabel('dY/dGsg')
ylim([0, 1.75])
plot(xs,mG_LT_3d(2:9),'m','linewidth',1.3)
plot(xs,mG_LT_4d(2:9),'y','linewidth',1.3)
plot(xs,mG_LT_5d(2:9),'--c','linewidth',1.3)
plot(xs, xsaxis,'k')
hold off
xlabel('Quarters')
saveas(gcf,'SYP_F6','psc2')
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Figure 6: Fiscal Multipliers in a Liquidity Trap with Decentralized Federalism', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
saveas(gcf,'SYP_F6','psc2')

% Calculate Impulse Response Functions

IRFy_Gb0 = 100*(M0b_BP(:,3)-M0a_BP(:,3));
IRFy_Gf1b = 100*(M1b_BP(:,5)-M1a_BP(:,5));
IRFy_Gf2b = 100*(M2b_BP(:,5)-M2a_BP(:,5));
IRFy_Gf3b = 100*(M3b_BP(:,5)-M3a_BP(:,5));
IRFy_Gf4b = 100*(M4b_BP(:,5)-M4a_BP(:,5));
IRFy_Gf5b = 100*(M5b_BP(:,5)-M5a_BP(:,5));

IRFy_GT1c = 100*(M1c_BP(:,5)-M1a_BP(:,5));
IRFy_GT2c = 100*(M2c_BP(:,5)-M2a_BP(:,5));
IRFy_GT3c = 100*(M3c_BP(:,5)-M3a_BP(:,5));
IRFy_GT5c = 100*(M5c_BP(:,5)-M5a_BP(:,5));

IRFy_Gs2d = 100*(M2d_BP(:,5)-M2a_BP(:,5));
IRFy_Gs3d = 100*(M3d_BP(:,5)-M3a_BP(:,5));
IRFy_Gs4d = 100*(M4d_BP(:,5)-M4a_BP(:,5));
IRFy_Gs5d = 100*(M5d_BP(:,5)-M5a_BP(:,5));

% Figure 7: Output Impulse Response Functions to Federal Spending Shocks
clf
subplot(3,1,1)
plot(xs,IRFy_Gb0(2:9),'--b','linewidth',1)
hold on
title('Federal Spending')
ylabel('dY/dGf')
ylim([0, 0.75])
plot(xs,IRFy_Gf1b(2:9),'r','linewidth',1)
plot(xs,IRFy_Gf2b(2:9),'g','linewidth',1)
plot(xs,IRFy_Gf3b(2:9),'m','linewidth',1)
plot(xs,IRFy_Gf4b(2:9),'y','linewidth',1)
plot(xs,IRFy_Gf5b(2:9),'--c','linewidth',1)
subplot(3,1,2)
plot(xs,IRFy_GT1c(2:9),'r','linewidth',1)
title('Federal Transfers')
hold on
ylabel('dY/dTF')
plot(xs,IRFy_GT2c(2:9),'g','linewidth',1)
plot(xs,IRFy_GT3c(2:9),'m','linewidth',1)
plot(xs,IRFy_GT5c(2:9),'--c','linewidth',1)
plot(xs, xsaxis,'k')
hold off
subplot(3,1,3)
plot(xs,IRFy_Gs2d(2:9),'g','linewidth',1)
hold on
title('Sub-Federal Spending')
xlabel('Quarters')
ylabel('dY/dGsg')
%ylim([-0.25, 1.75])
plot(xs,IRFy_Gs3d(2:9),'m','linewidth',1)
plot(xs,IRFy_Gs4d(2:9),'y','linewidth',1)
plot(xs,IRFy_Gs5d(2:9),'--c','linewidth',1)
plot(xs, xsaxis,'k')
hold off
xlabel('Quarters')
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Figure 7: Output Impulse Response Functions for Fiscal Shocks a Liquidity Trap', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
saveas(gcf,'SYP_F7','psc2')

% Figure 8: Fiscal Multipliers Above ZLB
clf
subplot(3,1,1)
plot(xs,mG_0b(2:9),'--b','linewidth',1.3)
hold on
title('Federal Spending')
ylabel('dY/dGf')
%ylim([1, 1.4])
plot(xs,mG_1b(2:9),'r','linewidth',1.3)
plot(xs,mG_2b(2:9),'g','linewidth',1.3)
plot(xs,mG_3b(2:9),'m','linewidth',1.3)
plot(xs,mG_4b(2:9),'y','linewidth',1.3)
plot(xs,mG_5b(2:9),'--c','linewidth',1.3)
subplot(3,1,2)
plot(xs,mG_1c(2:9),'r','linewidth',1.3)
title('Federal Transfers')
hold on
ylabel('dY/dTF')
plot(xs,mG_2c(2:9),'g','linewidth',1.3)
plot(xs,mG_3c(2:9),'m','linewidth',1.3)
plot(xs,mG_5c(2:9),'--c','linewidth',1.3)
plot(xs, xsaxis,'k')
hold off
subplot(3,1,3)
plot(xs,mG_2d(2:9),'g','linewidth',1.3)
hold on
title('Sub-Federal Spending')
xlabel('Quarters')
ylabel('dY/dGsg')
%ylim([0, 1.75])
plot(xs,mG_3d(2:9),'m','linewidth',1.3)
plot(xs,mG_4d(2:9),'y','linewidth',1.3)
plot(xs,mG_5d(2:9),'--c','linewidth',1.3)
plot(xs, xsaxis,'k')
hold off
xlabel('Quarters')
saveas(gcf,'SYP_F6','psc2')
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Figure 8: Fiscal Multipliers Above the Zero Lower Bound', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
saveas(gcf,'SYP_F8','psc2')

%% 8) Calculate impact on welfare with government spending enters utility function

% 8a) Perfect substitues:

M0a_W_altA = [M0a_W,M0a_BP(:,2)];
M0b_W_altA = [M0b_W,M0b_BP(:,2)];

M1a_W_altA = [M1a_W,M1a_BP(:,2) + M1a_BP(:,4)];
M1b_W_altA = [M1b_W,M1b_BP(:,2) + M1b_BP(:,4)];
M1c_W_altA = [M1c_W,M1c_BP(:,2) + M1c_BP(:,4)];
M1d_W_altA = [M1d_W,M1d_BP(:,2) + M1d_BP(:,4)];

M2a_W_altA = [M2a_W,M2a_BP(:,2) + M2a_BP(:,4)];
M2b_W_altA = [M2b_W,M2b_BP(:,2) + M2b_BP(:,4)];
M2c_W_altA = [M2c_W,M2c_BP(:,2) + M2c_BP(:,4)];
M2d_W_altA = [M2d_W,M2d_BP(:,2) + M2d_BP(:,4)];

M3a_W_altA = [M3a_W,M3a_BP(:,2) + M3a_BP(:,4)];
M3b_W_altA = [M3b_W,M3b_BP(:,2) + M3b_BP(:,4)];
M3c_W_altA = [M3c_W,M3c_BP(:,2) + M3c_BP(:,4)];
M3d_W_altA = [M3d_W,M3d_BP(:,2) + M3d_BP(:,4)];

M4a_W_altA = [M4a_W,M4a_BP(:,2) + M4a_BP(:,4)];
M4b_W_altA = [M4b_W,M4b_BP(:,2) + M4b_BP(:,4)];
M4c_W_altA = [M4c_W,M4c_BP(:,2) + M4c_BP(:,4)];
M4d_W_altA = [M4d_W,M4d_BP(:,2) + M4d_BP(:,4)];

M5a_W_altA = [M5a_W,M5a_BP(:,2) + M5a_BP(:,4)];
M5b_W_altA = [M5b_W,M5b_BP(:,2) + M5b_BP(:,4)];
M5c_W_altA = [M5c_W,M5c_BP(:,2) + M5c_BP(:,4)];
M5d_W_altA = [M5d_W,M5d_BP(:,2) + M5d_BP(:,4)];

[W0a_AZLB_altA,CW0a_AZLB_altA, LTW0a_AZLB_altA] = CalcWelfare(M0a_W_AZLB_altA,gamma,sigma);
[W1a_AZLB_altA,CW1a_AZLB_altA, LTW1a_AZLB_altA] = CalcWelfare(M1a_W_AZLB_altA,gamma,sigma);
[W2a_AZLB_altA,CW2a_AZLB_altA, LTW2a_AZLB_altA] = CalcWelfare(M2a_W_AZLB_altA,gamma,sigma);
[W3a_AZLB_altA,CW3a_AZLB_altA, LTW3a_AZLB_altA] = CalcWelfare(M3a_W_AZLB_altA,gamma,sigma);
[W4a_AZLB_altA,CW4a_AZLB_altA, LTW4a_AZLB_altA] = CalcWelfare(M4a_W_AZLB_altA,gamma,sigma);
[W5a_AZLB_altA,CW5a_AZLB_altA, LTW5a_AZLB_altA] = CalcWelfare(M5a_W_AZLB_altA,gamma,sigma);

[W0a_altA,CW0a_altA, LTW0a_altA] = CalcWelfare(M0a_W_altA,gamma,sigma);
[W0b_altA,CW0b_altA, LTW0b_altA] = CalcWelfare(M0b_W_altA,gamma,sigma);

[W1a_altA,CW1a_altA, LTW1a_altA] = CalcWelfare(M1a_W_altA,gamma,sigma);
[W1b_altA,CW1b_altA, LTW1b_altA] = CalcWelfare(M1b_W_altA,gamma,sigma);
[W1c_altA,CW1c_altA, LTW1c_altA] = CalcWelfare(M1c_W_altA,gamma,sigma);
[W1d_altA,CW1d_altA, LTW1d_altA] = CalcWelfare(M1d_W_altA,gamma,sigma);

[W2a_altA,CW2a_altA, LTW2a_altA] = CalcWelfare(M2a_W_altA,gamma,sigma);
[W2b_altA,CW2b_altA, LTW2b_altA] = CalcWelfare(M2b_W_altA,gamma,sigma);
[W2c_altA,CW2c_altA, LTW2c_altA] = CalcWelfare(M2c_W_altA,gamma,sigma);
[W2d_altA,CW2d_altA, LTW2d_altA] = CalcWelfare(M2d_W_altA,gamma,sigma);

[W3a_altA,CW3a_altA, LTW3a_altA] = CalcWelfare(M3a_W_altA,gamma,sigma);
[W3b_altA,CW3b_altA, LTW3b_altA] = CalcWelfare(M3b_W_altA,gamma,sigma);
[W3c_altA,CW3c_altA, LTW3c_altA] = CalcWelfare(M3c_W_altA,gamma,sigma);
[W3d_altA,CW3d_altA, LTW3d_altA] = CalcWelfare(M3d_W_altA,gamma,sigma);

[W4a_altA,CW4a_altA, LTW4a_altA] = CalcWelfare(M4a_W_altA,gamma,sigma);
[W4b_altA,CW4b_altA, LTW4b_altA] = CalcWelfare(M4b_W_altA,gamma,sigma);
[W4c_altA,CW4c_altA, LTW4c_altA] = CalcWelfare(M4c_W_altA,gamma,sigma);
[W4d_altA,CW4d_altA, LTW4d_altA] = CalcWelfare(M4d_W_altA,gamma,sigma);

[W5a_altA,CW5a_altA, LTW5a_altA] = CalcWelfare(M5a_W_altA,gamma,sigma);
[W5b_altA,CW5b_altA, LTW5b_altA] = CalcWelfare(M5b_W_altA,gamma,sigma);
[W5c_altA,CW5c_altA, LTW5c_altA] = CalcWelfare(M5c_W_altA,gamma,sigma);
[W5d_altA,CW5d_altA, LTW5d_altA] = CalcWelfare(M5d_W_altA,gamma,sigma);

dW0b_altA = -100*(W0b_altA(1:21)./W0a_altA(1:21)-1);

dW1b_altA = -100*(W1b_altA(1:21)./W1a_altA(1:21)-1);
dW1c_altA = -100*(W1c_altA(1:21)./W1a_altA(1:21)-1);
dW1d_altA = -100*(W1d_altA(1:21)./W1a_altA(1:21)-1);

dW2b_altA = -100*(W2b_altA(1:21)./W2a_altA(1:21)-1);
dW2c_altA = -100*(W2c_altA(1:21)./W2a_altA(1:21)-1);
dW2d_altA = -100*(W2d_altA(1:21)./W2a_altA(1:21)-1);

dW3b_altA = -100*(W3b_altA(1:21)./W3a_altA(1:21)-1);
dW3c_altA = -100*(W3c_altA(1:21)./W3a_altA(1:21)-1);
dW3d_altA = -100*(W3d_altA(1:21)./W3a_altA(1:21)-1);

dW4b_altA = -100*(W4b_altA(1:21)./W4a_altA(1:21)-1);
dW4c_altA = -100*(W4c_altA(1:21)./W4a_altA(1:21)-1);
dW4d_altA = -100*(W4d_altA(1:21)./W4a_altA(1:21)-1);

dW5b_altA = -100*(W5b_altA(1:21)./W5a_altA(1:21)-1);
dW5c_altA = -100*(W5c_altA(1:21)./W5a_altA(1:21)-1);
dW5d_altA = -100*(W5d_altA(1:21)./W5a_altA(1:21)-1);

% Table 4a (Welfare, perfect substitutes)

T4a = [dW0b_altA(2), 0, 0;
       dW1b_altA(2), dW1c_altA(2), 0;
       dW2b_altA(2), dW2c_altA(2), dW2d_altA(2);
       dW3b_altA(2), dW3c_altA(2), dW3d_altA(2);
       dW4b_altA(2), 0, dW4d_altA(2);
       dW5b_altA(2), dW5c_altA(2), dW5d_altA(2)];

% 8b) Perfect compliments:   

% Make sure sw = 0 is toggled on in CalcWelfare.m

M1a_W_altB = [M1a_W,M1a_BP(:,2) , M1a_BP(:,4)];
M1b_W_altB = [M1b_W,M1b_BP(:,2) , M1b_BP(:,4)];
M1c_W_altB = [M1c_W,M1c_BP(:,2) , M1c_BP(:,4)];
M1d_W_altB = [M1d_W,M1d_BP(:,2) , M1d_BP(:,4)];

M2a_W_altB = [M2a_W,M2a_BP(:,2) , M2a_BP(:,4)];
M2b_W_altB = [M2b_W,M2b_BP(:,2) , M2b_BP(:,4)];
M2c_W_altB = [M2c_W,M2c_BP(:,2) , M2c_BP(:,4)];
M2d_W_altB = [M2d_W,M2d_BP(:,2) , M2d_BP(:,4)];

M3a_W_altB = [M3a_W,M3a_BP(:,2) , M3a_BP(:,4)];
M3b_W_altB = [M3b_W,M3b_BP(:,2) , M3b_BP(:,4)];
M3c_W_altB = [M3c_W,M3c_BP(:,2) , M3c_BP(:,4)];
M3d_W_altB = [M3d_W,M3d_BP(:,2) , M3d_BP(:,4)];

M4a_W_altB = [M4a_W,M4a_BP(:,2) , M4a_BP(:,4)];
M4b_W_altB = [M4b_W,M4b_BP(:,2) , M4b_BP(:,4)];
M4c_W_altB = [M4c_W,M4c_BP(:,2) , M4c_BP(:,4)];
M4d_W_altB = [M4d_W,M4d_BP(:,2) , M4d_BP(:,4)];

M5a_W_altB = [M5a_W,M5a_BP(:,2) , M5a_BP(:,4)];
M5b_W_altB = [M5b_W,M5b_BP(:,2) , M5b_BP(:,4)];
M5c_W_altB = [M5c_W,M5c_BP(:,2) , M5c_BP(:,4)];
M5d_W_altB = [M5d_W,M5d_BP(:,2) , M5d_BP(:,4)];

[W1a_AZLB_altB,CW1a_AZLB_altB, LTW1a_AZLB_altB] = CalcWelfare(M1a_W_AZLB_altB,gamma,sigma);
[W2a_AZLB_altB,CW2a_AZLB_altB, LTW2a_AZLB_altB] = CalcWelfare(M2a_W_AZLB_altB,gamma,sigma);
[W3a_AZLB_altB,CW3a_AZLB_altB, LTW3a_AZLB_altB] = CalcWelfare(M3a_W_AZLB_altB,gamma,sigma);
[W4a_AZLB_altB,CW4a_AZLB_altB, LTW4a_AZLB_altB] = CalcWelfare(M4a_W_AZLB_altB,gamma,sigma);
[W5a_AZLB_altB,CW5a_AZLB_altB, LTW5a_AZLB_altB] = CalcWelfare(M5a_W_AZLB_altB,gamma,sigma);

[W1a_altB,CW1a_altB, LTW1a_altB] = CalcWelfare(M1a_W_altB,gamma,sigma);
[W1b_altB,CW1b_altB, LTW1b_altB] = CalcWelfare(M1b_W_altB,gamma,sigma);
[W1c_altB,CW1c_altB, LTW1c_altB] = CalcWelfare(M1c_W_altB,gamma,sigma);
[W1d_altB,CW1d_altB, LTW1d_altB] = CalcWelfare(M1d_W_altB,gamma,sigma);

[W2a_altB,CW2a_altB, LTW2a_altB] = CalcWelfare(M2a_W_altB,gamma,sigma);
[W2b_altB,CW2b_altB, LTW2b_altB] = CalcWelfare(M2b_W_altB,gamma,sigma);
[W2c_altB,CW2c_altB, LTW2c_altB] = CalcWelfare(M2c_W_altB,gamma,sigma);
[W2d_altB,CW2d_altB, LTW2d_altB] = CalcWelfare(M2d_W_altB,gamma,sigma);

[W3a_altB,CW3a_altB, LTW3a_altB] = CalcWelfare(M3a_W_altB,gamma,sigma);
[W3b_altB,CW3b_altB, LTW3b_altB] = CalcWelfare(M3b_W_altB,gamma,sigma);
[W3c_altB,CW3c_altB, LTW3c_altB] = CalcWelfare(M3c_W_altB,gamma,sigma);
[W3d_altB,CW3d_altB, LTW3d_altB] = CalcWelfare(M3d_W_altB,gamma,sigma);

[W4a_altB,CW4a_altB, LTW4a_altB] = CalcWelfare(M4a_W_altB,gamma,sigma);
[W4b_altB,CW4b_altB, LTW4b_altB] = CalcWelfare(M4b_W_altB,gamma,sigma);
[W4c_altB,CW4c_altB, LTW4c_altB] = CalcWelfare(M4c_W_altB,gamma,sigma);
[W4d_altB,CW4d_altB, LTW4d_altB] = CalcWelfare(M4d_W_altB,gamma,sigma);

[W5a_altB,CW5a_altB, LTW5a_altB] = CalcWelfare(M5a_W_altB,gamma,sigma);
[W5b_altB,CW5b_altB, LTW5b_altB] = CalcWelfare(M5b_W_altB,gamma,sigma);
[W5c_altB,CW5c_altB, LTW5c_altB] = CalcWelfare(M5c_W_altB,gamma,sigma);
[W5d_altB,CW5d_altB, LTW5d_altB] = CalcWelfare(M5d_W_altB,gamma,sigma);

dW1b_altB = -100*(W1b_altB(1:21)./W1a_altB(1:21)-1);
dW1c_altB = -100*(W1c_altB(1:21)./W1a_altB(1:21)-1);
dW1d_altB = -100*(W1d_altB(1:21)./W1a_altB(1:21)-1);

dW2b_altB = -100*(W2b_altB(1:21)./W2a_altB(1:21)-1);
dW2c_altB = -100*(W2c_altB(1:21)./W2a_altB(1:21)-1);
dW2d_altB = -100*(W2d_altB(1:21)./W2a_altB(1:21)-1);

dW3b_altB = -100*(W3b_altB(1:21)./W3a_altB(1:21)-1);
dW3c_altB = -100*(W3c_altB(1:21)./W3a_altB(1:21)-1);
dW3d_altB = -100*(W3d_altB(1:21)./W3a_altB(1:21)-1);

dW4b_altB = -100*(W4b_altB(1:21)./W4a_altB(1:21)-1);
dW4c_altB = -100*(W4c_altB(1:21)./W4a_altB(1:21)-1);
dW4d_altB = -100*(W4d_altB(1:21)./W4a_altB(1:21)-1);

dW5b_altB = -100*(W5b_altB(1:21)./W5a_altB(1:21)-1);
dW5c_altB = -100*(W5c_altB(1:21)./W5a_altB(1:21)-1);
dW5d_altB = -100*(W5d_altB(1:21)./W5a_altB(1:21)-1);

% Table 4b (Welfare, perfect compliments)

T4b = [dW0b_altA(2), 0, 0;
       dW1b_altB(2), dW1c_altB(2), 0;
       dW2b_altB(2), dW2c_altB(2), dW2d_altB(2);
       dW3b_altB(2), dW3c_altB(2), dW3d_altB(2);
       dW4b_altB(2), 0, dW4d_altB(2);
       dW5b_altB(2), dW5c_altB(2), dW5d_altB(2)];
   
% 8c) Seperate goods

% Make sure sw = 1 is toggled on in CalcWelfare.m

[W1a_altC,CW1a_altC, LTW1a_altC] = CalcWelfare(M1a_W_altB,gamma,sigma);
[W1b_altC,CW1b_altC, LTW1b_altC] = CalcWelfare(M1b_W_altB,gamma,sigma);
[W1c_altC,CW1c_altC, LTW1c_altC] = CalcWelfare(M1c_W_altB,gamma,sigma);
[W1d_altC,CW1d_altC, LTW1d_altC] = CalcWelfare(M1d_W_altB,gamma,sigma);

[W2a_altC,CW2a_altC, LTW2a_altC] = CalcWelfare(M2a_W_altB,gamma,sigma);
[W2b_altC,CW2b_altC, LTW2b_altC] = CalcWelfare(M2b_W_altB,gamma,sigma);
[W2c_altC,CW2c_altC, LTW2c_altC] = CalcWelfare(M2c_W_altB,gamma,sigma);
[W2d_altC,CW2d_altC, LTW2d_altC] = CalcWelfare(M2d_W_altB,gamma,sigma);

[W3a_altC,CW3a_altC, LTW3a_altC] = CalcWelfare(M3a_W_altB,gamma,sigma);
[W3b_altC,CW3b_altC, LTW3b_altC] = CalcWelfare(M3b_W_altB,gamma,sigma);
[W3c_altC,CW3c_altC, LTW3c_altC] = CalcWelfare(M3c_W_altB,gamma,sigma);
[W3d_altC,CW3d_altC, LTW3d_altC] = CalcWelfare(M3d_W_altB,gamma,sigma);

[W4a_altC,CW4a_altC, LTW4a_altC] = CalcWelfare(M4a_W_altB,gamma,sigma);
[W4b_altC,CW4b_altC, LTW4b_altC] = CalcWelfare(M4b_W_altB,gamma,sigma);
[W4c_altC,CW4c_altC, LTW4c_altC] = CalcWelfare(M4c_W_altB,gamma,sigma);
[W4d_altC,CW4d_altC, LTW4d_altC] = CalcWelfare(M4d_W_altB,gamma,sigma);

[W5a_altC,CW5a_altC, LTW5a_altC] = CalcWelfare(M5a_W_altB,gamma,sigma);
[W5b_altC,CW5b_altC, LTW5b_altC] = CalcWelfare(M5b_W_altB,gamma,sigma);
[W5c_altC,CW5c_altC, LTW5c_altC] = CalcWelfare(M5c_W_altB,gamma,sigma);
[W5d_altC,CW5d_altC, LTW5d_altC] = CalcWelfare(M5d_W_altB,gamma,sigma);

dW1b_altC = -100*(W1b_altC(1:21)./W1a_altC(1:21)-1);
dW1c_altC = -100*(W1c_altC(1:21)./W1a_altC(1:21)-1);
dW1d_altC = -100*(W1d_altC(1:21)./W1a_altC(1:21)-1);

dW2b_altC = -100*(W2b_altC(1:21)./W2a_altC(1:21)-1);
dW2c_altC = -100*(W2c_altC(1:21)./W2a_altC(1:21)-1);
dW2d_altC = -100*(W2d_altC(1:21)./W2a_altC(1:21)-1);

dW3b_altC = -100*(W3b_altC(1:21)./W3a_altC(1:21)-1);
dW3c_altC = -100*(W3c_altC(1:21)./W3a_altC(1:21)-1);
dW3d_altC = -100*(W3d_altC(1:21)./W3a_altC(1:21)-1);

dW4b_altC = -100*(W4b_altC(1:21)./W4a_altC(1:21)-1);
dW4c_altC = -100*(W4c_altC(1:21)./W4a_altC(1:21)-1);
dW4d_altC = -100*(W4d_altC(1:21)./W4a_altC(1:21)-1);

dW5b_altC = -100*(W5b_altC(1:21)./W5a_altC(1:21)-1);
dW5c_altC = -100*(W5c_altC(1:21)./W5a_altC(1:21)-1);
dW5d_altC = -100*(W5d_altC(1:21)./W5a_altC(1:21)-1);

% Table 4c (Distinct public goods)

T4c = [dW0b_altA(2), 0, 0;
       dW1b_altC(2), dW1c_altC(2), 0;
       dW2b_altC(2), dW2c_altC(2), dW2d_altC(2);
       dW3b_altC(2), dW3c_altC(2), dW3d_altC(2);
       dW4b_altC(2), 0, dW4d_altC(2);
       dW5b_altC(2), dW5c_altC(2), dW5d_altC(2)];
   
   
[W5a_altD,CW5a_altD, LTW5a_altD] = CalcWelfare(M5a_W_altB,gamma,sigma);
[W5b_altD,CW5b_altD, LTW5b_altD] = CalcWelfare(M5b_W_altB,gamma,sigma);
[W5c_altD,CW5c_altD, LTW5c_altD] = CalcWelfare(M5c_W_altB,gamma,sigma);
[W5d_altD,CW5d_altD, LTW5d_altD] = CalcWelfare(M5d_W_altB,gamma,sigma);

 CW5d_altD - CW5a_altD

%% 9) Welfare Analysis and Figures

% Figure 14: Per Period Utility Impact of Fiscal Shocks
xc = 1:1:20;
xcbar = zeros(length(xc),1);

clf
subplot(3,1,1)
plot(xc,dW0b_altA(2:21),'--b','linewidth',1)
hold on
title('Federal Spending')
ylabel('%dU')
ylim([0, 0.3])
plot(xc,dW1b_altC(2:21),'r','linewidth',1)
plot(xc,dW5b_altC(2:21),'c','linewidth',1)
legend('Centralized Model','Decentralized: Endogenous GS ','Decentralized: GR','Location','northeast','Orientation','vertical')
subplot(3,1,2)
plot(xc,dW1c_altC(2:21),'r','linewidth',1)
title('Federal Transfers')
hold on
ylim([0, 0.3])
ylabel('%dU')
plot(xc,dW5c_altC(2:21),'c','linewidth',1)
plot(xc, xcbar,'k')
hold off
subplot(3,1,3)
plot(xc,dW5d_altC(2:21),'c','linewidth',1)
title('Sub-Federal Spending')
hold on
xlabel('Quarters')
ylabel('%dU')
%ylim([-0.25, 1.75])
plot(xc, xcbar,'k')
xlabel('Quarters')
hold off
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Figure 14: Percent Change in Utility from Fiscal Shocks in a Liquidity Trap', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
saveas(gcf,'SYP_F14','psc2')

% Figure 15: Per Period Utility Impact of Fiscal Shocks

clf
subplot(3,1,1)
plot(xc,CW0b_altA(2:21) - CW0a_altA(2:21),'--b','linewidth',1)
hold on
title('Federal Spending')
ylabel('Utils')
%ylim([0, 0.3])
plot(xc,CW1b_altC(2:21) - CW1a_altC(2:21),'r','linewidth',1)
plot(xc,CW5b_altC(2:21) - CW5a_altC(2:21),'c','linewidth',1)
legend('Centralized Model','Decentralized: Endogenous GS ','Decentralized: GR','Location','southeast','Orientation','vertical')
subplot(3,1,2)
plot(xc,CW1c_altC(2:21) - CW1a_altC(2:21),'r','linewidth',1)
title('Federal Transfers')
hold on
%ylim([0, 0.3])
ylabel('Utils')
plot(xc,CW5c_altC(2:21) - CW5a_altC(2:21),'c','linewidth',1)
plot(xc, xcbar,'k')
hold off
subplot(3,1,3)
plot(xc,CW5d_altC(2:21) - CW5a_altC(2:21),'c','linewidth',1)
title('Sub-Federal Spending')
hold on
xlabel('Quarters')
ylabel('Utils')
%ylim([-0.25, 1.75])
plot(xc, xcbar,'k')
xlabel('Quarters')
hold off
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Figure 15: Cumulative Change in Utility from Fiscal Shocks in a Liquidity Trap', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')
saveas(gcf,'SYP_F15','psc2')

% 
% CRFy_Gb0 = zeros(21);
% CRFy_Gf1b = zeros(21);
% CRFy_Gf2b = zeros(21);
% CRFy_Gf3b = zeros(21);
% 
% for i =2:21
%     CRFy_Gb0(i) = CRFy_Gb0(i-1)+ IRFy_Gb0(i);
%     CRFy_Gf1b(i) = CRFy_Gf1b(i-1) + IRFy_Gf1b(i);
%     CRFy_Gf2b(i) = CRFy_Gf2b(i-1) + IRFy_Gf2b(i);
%     CRFy_Gf3b(i) = CRFy_Gf3b(i-1) + IRFy_Gf3b(i);
% end
% 
% % Figure 14: Cumulative Response Functions to Fiscal Shocks
% clf
% plot(x,CRFy_Gb0(1:21),'--b')
% hold on
% title('Figure 14: Output Cumulative Response Functions to Fiscal Shocks')
% xlabel('Quarters')
% ylabel('y*-y')
% plot(x,CRFy_Gf1b(1:21),'r')
% plot(x,CRFy_Gf2b(1:21),'g')
% plot(x,CRFy_Gf3b(1:21),'m')
% hold off
% saveas(gcf,'SYP_F14','psc2')
% 
% % Figure 15: Federal Tranfser Multipliers
% clf
% plot(x,mG_LT_1c(1:21),'r')
% hold on
% title('Figure 12: Federal Transfer Multipliers in a Liquidity Trap')
% xlabel('Quarters')
% ylabel('dy/dg')
% plot(x,mG_LT_2c(1:21),'g')
% plot(x,mG_LT_3c(1:21),'m')
% legend('Endogenous SLG Spending ','Endogenous SLG Tau C','Endogenous SLG Tau N','Location','southeast','Orientation','vertical')
% hold off
% saveas(gcf,'SYP_F15','psc2')
% 
% 
% IRFy_Gf1c = 100*(M1c_BP(:,5)-M1a_BP(:,5));
% IRFy_Gf2c = 100*(M2c_BP(:,5)-M2a_BP(:,5));
% IRFy_Gf3c = 100*(M3c_BP(:,5)-M3a_BP(:,5));
% 
% % Figure 16: Output Impulse Response Functions to Federal Transfer Shocks
% clf
% plot(x,IRFy_Gf1c(1:21),'r')
% hold on
% title('Figure 13: Output Impulse Response Functions to Federal Transfer Shocks')
% xlabel('Quarters')
% ylabel('Percentage Deviation, With Stimulus Less Baseline Liquidity Trap ')
% plot(x,IRFy_Gf2c(1:21),'g')
% plot(x,IRFy_Gf3c(1:21),'m')
% legend('Endogenous SLG Spending ','Endogenous SLG Tau C','Endogenous SLG Tau N','Location','northeast','Orientation','vertical')
% hold off
% saveas(gcf,'SYP_F16','psc2')
% 
% %% Plot Endogenous Policy Instruments, Inflation, and Labor Supply  
% 
% % Figure 17: Model 1: Compare Transfer and Spending Shocks
% 
% % Plot comparison of M0b, M4b
% x = 0:1:20;
% clf;
% subplot(3,3,1)
% plot(x,dM1b_BP(1:21,3),'b',x,dM1c_BP(1:21,3),'r')
% title('SLG Spending')
% ylabel('Percent Deviation')
% hold on
% title('Z')
% subplot(3,3,2)
% plot(x,4*M1b_Less_M1a(1:21,4),'b',x,4*M1c_Less_M1a(1:21,4),'r')
% title('Inflation')
% ylabel('Percent (Annualized)')
% subplot(3,3,3)
% plot(x,M1b_Less_M1a(1:21,3),'b',x,M1c_Less_M1a(1:21,3),'r')
% title('Hours')
% ylabel('Percent Deviation')
% legend('Federal Spending','Federal Transfers','Location','northeast','Orientation','vertical')
% subplot(3,3,4)
% plot(x,100*dM2b_BP(1:21,11),'b',x,100*dM2c_BP(1:21,11),'r')
% title('Consumption Tax Rate')
% ylabel('Percentage Points')
% subplot(3,3,5)
% plot(x,4*M2b_Less_M2a(1:21,4),'b',x,4*M2c_Less_M2a(1:21,4),'r')
% title('Inflation')
% ylabel('Percent (Annualized)')
% subplot(3,3,6)
% plot(x,M2b_Less_M2a(1:21,3),'b',x,M2c_Less_M2a(1:21,3),'r')
% title('Hours')
% ylabel('Percent Deviation')
% subplot(3,3,7)
% plot(x,100*dM3b_BP(1:21,10),'b',x,100*dM3c_BP(1:21,10),'r')
% title('Labor Income Tax Rate')
% xlabel('Quarters')
% ylabel('Percentage Points')
% subplot(3,3,8)
% plot(x,4*M3b_Less_M3a(1:21,4),'b',x,4*M3c_Less_M3a(1:21,4),'r')
% xlabel('Quarters')
% title('Inflation')
% ylabel('Percent (Annualized)')
% subplot(3,3,9)
% plot(x,M3b_Less_M3a(1:21,3),'b',x,M3c_Less_M3a(1:21,3),'r')
% title('Hours')
% xlabel('Quarters')
% ylabel('Percent Deviation')
% annotation('textbox', [0 0.9 1 0.1], ...
%     'String', 'Figure 17: Endogenous Responses to Transfer vs Federal Spending Shocks', ...
%     'EdgeColor', 'none', ...
%     'HorizontalAlignment', 'center')
% annotation('textbox', [0.9 0.8 0.5 0.5], ...
%     'String', 'All Responses Are Measured as the Difference Relative to a Liquidity Trap without Any Fiscal Stabilization Policy ', ...
%     'EdgeColor', 'none', ...
%     'HorizontalAlignment', 'center')
% hold off
% saveas(gcf,'SYP_F17','psc2')
% 
% % Figure 18: Output Impulse Response Functions to Federal Transfer Shocks
% 
% IRFy_Gf1d = 100*(M1d_BP(:,5)-M1a_BP(:,5));
% IRFy_Gf2d = 100*(M2d_BP(:,5)-M2a_BP(:,5));
% IRFy_Gf3d = 100*(M3d_BP(:,5)-M3a_BP(:,5));
% 
% clf
% plot(x,IRFy_Gb0(1:21),'--b')
% hold on
% title('Figure 18: Output Impulse Response Functions to Federal Transfer Shocks')
% xlabel('Quarters')
% ylabel('Percentage Deviation, With Stimulus Less Baseline Liquidity Trap ')
% plot(x,IRFy_Gf1d(1:21),'r')
% plot(x,IRFy_Gf2d(1:21),'g')
% plot(x,IRFy_Gf3d(1:21),'m')
% legend('Centralized Model','Decentralized Model: Endog. GS ','Decentralized Model: Endog. TauC','Decentralized Model: Endog. TauN','Location','northeast','Orientation','vertical')
% hold off
% saveas(gcf,'SYP_F18','psc2')

%% 9) Generate Appendix Figures on F, SLG Net Lending

% Import data
NetLendingData = xlsread('/Users/afieldhouse/Documents/MATLAB/Correia_et_al_(2013)/Fieldhouse_adapted/Fig20');
NetLendingData = NetLendingData(2:53,:);

yr = NetLendingData(:,1);
xaxis = zeros(length(yr),1);
F_NL_PTR = 100*NetLendingData(:,2);
SL_NL_PTR = 100*NetLendingData(:,3);
S_NL_PTR = 100*NetLendingData(:,4);
L_NL_PTR = 100*NetLendingData(:,5);

F_NL_PGDP = 100*NetLendingData(:,6)./ NetLendingData(:,10);
SL_NL_PGDP = 100*NetLendingData(:,7)./ NetLendingData(:,10);
S_NL_PGDP = 100*NetLendingData(:,8)./ NetLendingData(:,10);
L_NL_PGDP = 100*NetLendingData(:,9)./ NetLendingData(:,10);

% Figure A1 - Net Lending % PotGDP

clf
plot(yr,F_NL_PGDP,'k')
hold on
title('Figure A1: Net Governemnt Lending/Borrowing as a Percent of Potential GDP')
xlabel({' ','Source: Bureau of Economic Analysis, Congressional Budget Office','Net Lending (+) / Borrowing (-) = Total Receipts - Total Expenditures'})
ylabel('Percent')
xlim([1961,2012])
plot(yr,SL_NL_PGDP,'r')
plot(yr,S_NL_PGDP,'r--')
plot(yr,L_NL_PGDP,'m--')
plot(yr,xaxis,'k')
legend('Federal','State and Local','State','Local','Location','southwest','Orientation','vertical')
hold off
saveas(gcf,'SYP_FA1','psc2')

% Figure A2 - Net Lending % PotReciepts

clf
plot(yr,F_NL_PTR,'k')
hold on
title('Figure A2: Net Governemnt Lending/Borrowing as a Percent of Total Receipts')
xlabel({'Source: Bureau of Economic Analysis, Congressional Budget Office','Net Lending (+) / Borrowing (-) = Total Receipts - Total Expenditures'})
ylabel('Percent')
xlim([1961,2012])
plot(yr,SL_NL_PTR,'r')
plot(yr,S_NL_PTR,'r--')
plot(yr,L_NL_PTR,'m--')
legend('Federal','State and Local','State','Local','Location','southwest','Orientation','vertical')
plot(yr,xaxis,'k')
hold off
saveas(gcf,'SYP_FA2','psc2')

% Import data
TaylorSLG = xlsread('/Users/afieldhouse/Documents/MATLAB/Correia_et_al_(2013)/Fieldhouse_adapted/FigA3');

xaxisQ= 2007:0.25:2015;

% Figure A3: State and Local Budget Composition, 2007Q1-2015Q1

clf
plot(xaxisQ,TaylorSLG(:,2),'r')
hold on
title('Figure A3: Change in State and Local Budget Composition from 2008Q4')
xlabel({' ','Source: Bureau of Economic Analysis','Net Lending (+) / Borrowing (-) = Total Receipts - Total Expenditures'})
ylabel('Billions of $2009 Dollars (Annualized)')
plot(xaxisQ,TaylorSLG(:,3),'b--')
plot(xaxisQ,TaylorSLG(:,4),'b')
plot(xaxisQ,TaylorSLG(:,6),'k--')
plot(xaxisQ,TaylorSLG(:,5),'k')
legend('Net Borrowing','Government Purchases','Other Expenditure','Fed. Transfers','Other Receipts','Location','northwest','Orientation','vertical')
hold off
saveas(gcf,'SYP_FA3','psc2')



%%  OLD CODE 


% Cold Model 4 (mixed policy instead of endogenous fed. transfers)
%% 4) Model 4: Decentralized Federalism: Quadratic Loss Function for Endogenous Policy Response

% % 4a) Run baseline liquidity trap with preference shock, decentralized government
% dynare Fieldhouse_sticPflexW_adj4.mod noclearall
% 
% 
% % Testing likeness to Great Recession policy responses
% PRatios = [gs/gsss, tau_c/tau_css, tau_n/tau_nss]
% 
% 
% 
% M4a_BA = [Z, dI, dn, PI-1, i, dy, dc, deta - beta, i - PI+1];
% 
% % Save benchmark path w/ decentralization
% if length(tau_nf) == 1
%     tau_nf = repmat(tau_nf,length(M4a_BA),1);
% end
% if length(tau_ns) == 1
%     tau_ns = repmat(tau_ns,length(M4a_BA),1);
% end
% if length(tau_c) == 1
%     tau_c = repmat(tau_c,length(M4a_BA),1);
% end
% if length(gs) == 1
%     gs = repmat(gs,length(M4a_BA),1);
% end
% 
% M4a_BP  = [dgf, gf, dgs, gs, dy, y, dTf, tau_nf, tau_ns, tau_c];
% 
% annotation('textbox', [0 0.9 1 0.1], ...
%     'String', 'Figure 9: Allocation IRFs with Centralized Fiscal Federalism, Endogenous TauNS', ...
%     'EdgeColor', 'none', ...
%     'HorizontalAlignment', 'center')
% saveas(gcf,'SYP_F9','psc2')
% 
% % Plot comparison of M0a, M4a
% x = 0:1:20;
% clf;
% subplot(3,3,1)
% plot(x,4*M0a_BA(1:21,1),'b',x,4*M4a_BA(1:21,1),'r')
% ylabel('Percent (Annualized)')
% hold on
% title('Z')
% subplot(3,3,2)
% plot(x,M0a_BA(1:21,2),'b',x,M4a_BA(1:21,2),'r')
% title('Investment')
% ylabel('Percentage Deviation')
% subplot(3,3,3)
% plot(x,M0a_BA(1:21,3),'b',x,M4a_BA(1:21,3),'r')
% title('Hours')
% ylabel('Percentage Deviation')
% subplot(3,3,4)
% plot(x,4*M0a_BA(1:21,4),'b',x,4*M4a_BA(1:21,4),'r')
% title('Inflation')
% ylabel('Percent (Annualized)')
% subplot(3,3,5)
% plot(x,4*M0a_BA(1:21,5),'b',x,4*M4a_BA(1:21,5),'r')
% title('Nominal Interest Rate')
% ylabel('Percent (Annualized)')
% subplot(3,3,6)
% plot(x,M0a_BA(1:21,6),'b',x,M4a_BA(1:21,6),'r')
% title('Output')
% ylabel('Percentage Deviation')
% subplot(3,3,7)
% plot(x,M0a_BA(1:21,7),'b',x,M4a_BA(1:21,7),'r')
% title('Consumption')
% xlabel('Quarters')
% ylabel('Percentage Deviation')
% subplot(3,3,8)
% plot(x,4*M0a_BA(1:21,8),'b',x,4*M4a_BA(1:21,8),'r')
% xlabel('Quarters')
% ylabel('Percent (Annualized)')
% title('Discount Factor Shock')
% subplot(3,3,9)
% plot(x,4*M0a_BA(1:21,9),'b',x,4*M4a_BA(1:21,9),'r')
% title('Real Interest Rate')
% xlabel('Quarters')
% ylabel('Percent (Annualized)')
% hold off
% annotation('textbox', [0 0.9 1 0.1], ...
%     'String', 'Figure 10: Liquidity Trap IRFs with Centralized Fiscal Federalism, Endogenous TauN', ...
%     'EdgeColor', 'none', ...
%     'HorizontalAlignment', 'center')
% saveas(gcf,'SYP_F10','psc2')
% 
% % Check relative output and allocation responses rel. benchmark
% M4a_Relative_M0a = M4a_BA./M0a_BA-1;
% M4a_Less_M0a = M4a_BA-M0a_BA;
% 
% % Check level response of output rel. benchmark
% y_level_3_0a = M4a_BP(:,6)./M0a_BP(:,4)-1;
% 
% % Check relative output and allocation responses rel. M1
% M4a_Relative_M1a = M4a_BA./M1a_BA-1;
% M4a_Less_M1a = M4a_BA-M1a_BA;
% 
% % Check level response of output rel. M1
% y_level_3_1 = M4a_BP(:,6)./M1a_BP(:,6)-1;
% 
% % Check relative output and allocation responses rel. M2
% M4a_Relative_M2a = M4a_BA./M2a_BA-1;
% M4a_Less_M2a = M4a_BA-M2a_BA;
% 
% % Check level response of output rel. M2
% y_level_3_2 = M4a_BP(:,6)./M2a_BP(:,6)-1;
% 
% 
% 
% % % 4b) Run decentralized liquidity simulation trap with additional federal spending shocks
% % 
% % % Make sure exGf is toggled on in Fieldhouse_sticPflexW_adj1.mod
% % % (lines 184-186)
% % 
% % dynare Fieldhouse_sticPflexW_adj4.mod noclearall
% % 
% % M4b_BA = [Z, dI, dn, PI-1, i, dy, dc, deta - beta, i - PI+1];
% % 
% % % Save benchmark path w/ decentralization
% % if length(tau_nf) == 1
% %     tau_nf = repmat(tau_nf,length(M4b_BA),1);
% % end
% % if length(tau_ns) == 1
% %     tau_ns = repmat(tau_ns,length(M4b_BA),1);
% % end
% % if length(tau_c) == 1
% %     tau_c = repmat(tau_c,length(M4b_BA),1);
% % end
% % if length(gs) == 1
% %     gs = repmat(gs,length(M4a_BA),1);
% % end
% % 
% % M4b_BP  = [dgf, gf, dgs, gs, dy, y, dTf, tau_nf, tau_ns, tau_c];
% % 
% % % Plot comparison of M0b, M4b
% % x = 0:1:20;
% % clf;
% % subplot(3,3,1)
% % plot(x,4*M0b_BA(1:21,1),'b',x,4*M4b_BA(1:21,1),'r')
% % ylabel('Percent (Annualized)')
% % hold on
% % title('Z')
% % subplot(3,3,2)
% % plot(x,M0b_BA(1:21,2),'b',x,M4b_BA(1:21,2),'r')
% % title('Investment')
% % ylabel('Percentage Deviation')
% % subplot(3,3,3)
% % plot(x,M0b_BA(1:21,3),'b',x,M4b_BA(1:21,3),'r')
% % title('Hours')
% % ylabel('Percentage Deviation')
% % subplot(3,3,4)
% % plot(x,4*M0b_BA(1:21,4),'b',x,4*M4b_BA(1:21,4),'r')
% % title('Inflation')
% % ylabel('Percent (Annualized)')
% % subplot(3,3,5)
% % plot(x,4*M0b_BA(1:21,5),'b',x,4*M4b_BA(1:21,5),'r')
% % title('Nominal Interest Rate')
% % ylabel('Percent (Annualized)')
% % subplot(3,3,6)
% % plot(x,M0b_BA(1:21,6),'b',x,M4b_BA(1:21,6),'r')
% % title('Output')
% % ylabel('Percentage Deviation')
% % subplot(3,3,7)
% % plot(x,M0b_BA(1:21,7),'b',x,M4b_BA(1:21,7),'r')
% % title('Consumption')
% % xlabel('Quarters')
% % ylabel('Percentage Deviation')
% % subplot(3,3,8)
% % plot(x,4*M0b_BA(1:21,8),'b',x,4*M4b_BA(1:21,8),'r')
% % xlabel('Quarters')
% % ylabel('Percent (Annualized)')
% % title('Discount Factor Shock')
% % subplot(3,3,9)
% % plot(x,4*M0b_BA(1:21,9),'b',x,4*M4b_BA(1:21,9),'r')
% % title('Real Interest Rate')
% % xlabel('Quarters')
% % ylabel('Percent (Annualized)')
% % hold off
% % annotation('textbox', [0 0.9 1 0.1], ...
% %     'String', 'Figure 11: Federal Spending Stimulus in a Liquidity Trap: Benchmark v. Model 2', ...
% %     'EdgeColor', 'none', ...
% %     'HorizontalAlignment', 'center')
% % saveas(gcf,'SYP_F11','psc2')
% % 
% % % Check relative output and allocation responses 
% % M4b_Relative_M0b = M4b_BA./M0b_BA-1;
% % M4b_Less_M0b = M4b_BA-M0b_BA;
% % 
% % M4b_Less_M1b = M4b_BA-M1b_BA;
% % 
% % % Check level response of output
% % y_level_3_0b = M4b_BP(:,6)./M0b_BP(:,4)-1;
% % 
% % % Calculate traditional federal govt spending multipler:
% % 
% % mGf3 = (M4b_BP(:,6)-M4a_BP(:,6))./(M4b_BP(:,2)-M4a_BP(:,2));
% % mGdelta3 = mGf3 - mGb0;
% % mGdeltaSR3 = mean(mGdelta3(2:21));

% % % clf;
% % % subplot(6,1,1)
% % % plot(IRF_allocation3(:,6),'g')
% % % title('Output Response')
% % % ylabel('Percent Deviation')
% % % %xlabel('Quarters')
% % % xlim([0 400])
% % % hold on
% % % subplot(6,1,2)
% % % plot(IRF_policy3(:,1),'b')
% % % title('Federal Governemnt Spending')
% % % ylabel('Percent Deviation')
% % % %xlabel('Quarters')
% % % xlim([0 400])
% % % subplot(6,1,3)
% % % plot(IRF_policy3(:,2),'b')
% % % title('Net Sub-Federal Governemnt Spending')
% % % ylabel('Percent Deviation')
% % % %xlabel('Quarters')
% % % xlim([0 400])
% % % subplot(6,1,4)
% % % plot(IRF_policy3(:,3),'b')
% % % title('Sub-Federal Discretionary Austerity')
% % % ylabel('Percent Deviation')
% % % %xlabel('Quarters')
% % % xlim([0 400])
% % % subplot(6,1,5)
% % % plot(IRF_policy3(:,4),'r')
% % % title('Federal Governemnt Transfer Spending')
% % % ylabel('Percent Deviation')
% % % %xlabel('Quarters')
% % % xlim([0 400])
% % % subplot(6,1,6)
% % % plot(IRF_policy3(:,5),'r')
% % % title('Labor income tax')
% % % ylabel('Percent')
% % % %xlabel('Quarters')
% % % xlim([0 400])
% % % hold off
% % % annotation('textbox', [0 0.9 1 0.1], ...
% % %     'String', 'Figure 14: Output and Policy Response to GF Shock, Decentralized Govt. w/ Static Subfederal Taxes', ...
% % %     'EdgeColor', 'none', ...
% % %     'HorizontalAlignment', 'center')
% % % saveas(gcf,'SYP_F14','psc2')
% % % 
% % % %% Old Code
% % % 
% % % % Pt = [];
% % % % Pt(1) = 1;
% % % % for i = 2:length(PI)
% % % %     Pt(i) = Pt(i-1)*PI(i);
% % % % end
% % % 
% % % % Check for budget balance
% % % deficit_f = n.*w.*tau_n + tau_k.* k .* (R - delta) - gf - TT;
% % % deficit_s = c.*tau_c + TT - gs;
% % % %deficit = g - c.*tau_c - n.*w.*tau_n - tau_k.* k .* (R - Pt' * delta);



%% 6) Comparison of allocations accross Benchmark and four models

% % Figure 5: Plot comparison of M1a, M2a, M3a all rel. M0a
% x = 0:1:20;
% clf;
% subplot(3,3,1)
% plot(x,4*M1a_Less_M0a(1:21,1),'r',x,4*M2a_Less_M0a(1:21,1),'g',x,4*M3a_Less_M0a(1:21,1),'m')
% ylabel('Percent (Annualized)')
% ylim([-0.02, 0.12])
% hold on
% title('Monetary Policy Without ZLB')
% legend('Model 1: GS ','Model 2: TauC','Model 2: TauN','Location','northeast','Orientation','vertical')
% subplot(3,3,2)
% plot(x,M1a_Less_M0a(1:21,2),'r',x,M2a_Less_M0a(1:21,2),'g',x,M3a_Less_M0a(1:21,2),'m')
% title('Investment')
% ylabel('Percentage Deviation')
% subplot(3,3,3)
% plot(x,M1a_Less_M0a(1:21,3),'r',x,M2a_Less_M0a(1:21,3),'g',x,M3a_Less_M0a(1:21,3),'m')
% title('Hours')
% ylabel('Percentage Deviation')
% subplot(3,3,4)
% plot(x,4*M1a_Less_M0a(1:21,4),'r',x,4*M2a_Less_M0a(1:21,4),'g',x,4*M3a_Less_M0a(1:21,4),'m')
% title('Inflation')
% ylabel('Percent (Annualized)')
% subplot(3,3,5)
% plot(x,4*M1a_Less_M0a(1:21,5),'r',x,4*M2a_Less_M0a(1:21,5),'g',x,4*M3a_Less_M0a(1:21,5),'m')
% title('Nominal Interest Rate')
% ylabel('Percent (Annualized)')
% subplot(3,3,6)
% plot(x,M1a_Less_M0a(1:21,6),'r',x,M2a_Less_M0a(1:21,6),'g',x,M3a_Less_M0a(1:21,6),'m')
% title('Output')
% ylabel('Percentage Deviation')
% subplot(3,3,7)
% plot(x,M1a_Less_M0a(1:21,7),'r',x,M2a_Less_M0a(1:21,7),'g',x,M3a_Less_M0a(1:21,7),'m')
% title('Consumption')
% xlabel('Quarters')
% ylabel('Percentage Deviation')
% subplot(3,3,8)
% plot(x,4*M1a_Less_M0a(1:21,9),'r',x,4*M2a_Less_M0a(1:21,9),'g',x,4*M3a_Less_M0a(1:21,9),'m')
% title('Real Interest Rate')
% xlabel('Quarters')
% ylabel('Percent (Annualized)')
% subplot(3,3,9)
% plot(x,dPR1a(1:21),'r',x,dPR2a(1:21),'g',x,dPR3a(1:21),'m')
% title('Endogenous Fiscal Response')
% xlabel('Quarters')
% ylabel('Percent (Annualized)')
% hold off
% annotation('textbox', [0 0.9 1 0.1], ...
%     'String', 'Figure 5: Liquidity Trap IRFs, Decentralized Models of Fiscal Federalism Less Centralized', ...
%     'EdgeColor', 'none', ...
%     'HorizontalAlignment', 'center')
% saveas(gcf,'SYP_F5','psc2')

% % Figure 6: Plot comparison of M1b-M1a, M2b-M2a, M3b-M3a all rel. M0b-M0a
% x = 0:1:20;
% clf;
% subplot(3,3,1)
% plot(x,100*4*M1b_Less_M1a(1:21,1),'r',x,4*M2b_Less_M2a(1:21,1),'g',x,4*M3b_Less_M3a(1:21,1),'m',x,4*M0b_Less_M0a(1:21,1),'--k')
% ylabel('Percent (Annualized)')
% ylim([-02, 12])
% hold on
% title('Monetary Policy Without ZLB')
% legend('Model 1: GS ','Model 2: TauC','Model 3: TauN','Benchmark','Location','northeast','Orientation','vertical')
% subplot(3,3,2)
% plot(x,100*M1b_Less_M1a(1:21,2),'r',x,100*M2b_Less_M2a(1:21,2),'g',x,100*M3b_Less_M3a(1:21,2),'m',x,100*M0b_Less_M0a(1:21,2),'--k')
% title('Investment')
% ylabel('Percentage Deviation')
% subplot(3,3,3)
% plot(x,100*M1b_Less_M1a(1:21,3),'r',x,100*M2b_Less_M2a(1:21,3),'g',x,100*M3b_Less_M3a(1:21,3),'m',x,100*M0b_Less_M0a(1:21,3),'--k')
% title('Hours')
% ylabel('Percentage Deviation')
% subplot(3,3,4)
% plot(x,100*4*M1b_Less_M1a(1:21,4),'r',x,100*4*M2b_Less_M2a(1:21,4),'g',x,100*4*M3b_Less_M3a(1:21,4),'m',x,100*4*M0b_Less_M0a(1:21,4),'--k')
% title('Inflation')
% ylabel('Percent (Annualized)')
% subplot(3,3,5)
% plot(x,100*4*M1b_Less_M1a(1:21,5),'r',x,100*4*M2b_Less_M2a(1:21,5),'g',x,100*4*M3b_Less_M3a(1:21,5),'m',x,100*4*M0b_Less_M0a(1:21,5),'--k')
% title('Nominal Interest Rate')
% ylabel('Percent (Annualized)')
% subplot(3,3,6)
% plot(x,100*M1b_Less_M1a(1:21,6),'r',x,100*M2b_Less_M2a(1:21,6),'g',x,100*M3b_Less_M3a(1:21,6),'m',x,100*M0b_Less_M0a(1:21,6),'--k')
% title('Output')
% ylabel('Percentage Deviation')
% subplot(3,3,7)
% plot(x,100*M1b_Less_M1a(1:21,7),'r',x,100*M2b_Less_M2a(1:21,7),'g',x,100*M3b_Less_M3a(1:21,7),'m',x,100*M0b_Less_M0a(1:21,7),'--k')
% title('Consumption')
% xlabel('Quarters')
% ylabel('Percentage Deviation')
% subplot(3,3,8)
% plot(x,100*4*M1b_Less_M1a(1:21,9),'r',x,100*4*M2b_Less_M2a(1:21,9),'g',x,100*4*M3b_Less_M3a(1:21,9),'m',x,100*4*M0b_Less_M0a(1:21,9),'--k')
% title('Real Interest Rate')
% xlabel('Quarters')
% ylabel('Percent (Annualized)')
% subplot(3,3,9)
% plot(x,dM1b_BP(1:21,3),'r',x,dM2b_BP(1:21,11),'g',x,dM3b_BP(1:21,10),'m',x,dM0b_BP(1:21,1),'--k')
% title('Endogenous Fiscal Response')
% xlabel('Quarters')
% ylabel('Percent (Annualized)')
% hold off
% annotation('textbox', [0 0.9 1 0.1], ...
%     'String', 'Figure 6: Liquidity Trap IRFs, Government Spending Shock Less Liquidity Trap Baseline', ...
%     'EdgeColor', 'none', ...
%     'HorizontalAlignment', 'center')
% saveas(gcf,'SYP_F6','psc2')

