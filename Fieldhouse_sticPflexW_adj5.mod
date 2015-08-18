// Extension of Correia et al. (2013), based on their public code

// Andrew Fieldhouse 
// Cornell University

// Decentralized Model of Fiscal Federalism Version 5:

// The sub-federal fiscal adjustement mimics Great Recession policy response


warning off;

options_.maxit_=1000000;
var y, c, I, k, gf, gs, tau_c, tau_nsDisc, Tf, TT, n, v , k_hat, PI, PI_tilde, A, D, mc, w, R, lambda, mu, i, deta, dy, A_hat, D_hat, dI, dc, dR, Z, dn, u_c, u_n, deta_hat, U, L, dgs, dgf, dTf, tau_n, SLGdeficit, gsDisc;

// y is output, c is consumption, I is investment, k is capital, gf is fed government expenditure, gs is sub-federal government expenditure, 
// Tf is federal to subfederal transfers, TT = max(Tf,0), n is labor, tau_c is sub-federal consumption tax,
// v is the distortion of relative prices, k_hat is utilized capital (I've set this to equal capital)
// PI is inflation +1, PI_tilde is reset inflation +1, A and D are used for recursive sums in optimal
// price setting, mc is marginal cost, w is wage, R is return to capital, lambda is the lagrange multiplier
// on the budget constraint, mu is the lagrange multiplier on the capital accumulation constraint, 
// i is the interest rate, deta is the discount factor, xc is a persistent shock (unused)

varexo e_i exc exGf exGs exTf extau_c extau_nf extau_ns exi;

//exc is the discount factor shock, exGf is the federal expenditure shock
//exGf is the sub-federal expenditure shock, exTf is the federal transfer shock

parameters Lss, tau_k, tau_nf, rho_r,sigma, sigma_I,beta,gamma,alpha, theta, epsilon, DELTA, rho_g, phi_pi, phi_y, omega_f, omega_s, omega_t, psi, nu, delta, gsDiscss, tau_nsDiscss, tau_ns,
   nss,Rss,mcss, wss, k_hatss,  epsilon_w,  kss, css, yss, iss, lambdass, muss, Iss, uss, vss, PIss, PI_tildess, Ass, Dss, gsss, gfss, qss, rho_phi, rc, tau_nss, tau_css, SLGdeficitss;

beta = 0.99; //discount factor
alpha = 0.3; //capital share
omega_f = 0.11; //ratio of federal govt expenditure
omega_s = 0.11; //ratio of sub-federal govt expenditure
omega_t = 0.023; // ratio of federal govt transfers to sub-federal government
delta = 0.02; //depreciation 
gamma =0.29; //utility parameter
sigma = 2; //risk aversion
epsilon = 7; //intermediate good substitubility
sigma_I = 17; //adjustment parameter
rho_g = 0.8; //presistence in govt
phi_pi=1.5; //taylor for inflation
phi_y = 0.0; //taylor for output
theta = 0.85; //calvo
DELTA = (1/beta-1)/delta +1; //unused for now
rho_r = 0;
tau_nf =  0.33;
tau_ns =  0.044;
tau_k = 0.24;
//tau_c= 0.096;
epsilon_w = 3;
psi = 2.9/3.9;
nu = 1/2.2;

// tau_nf is federal labor income tax, tau_ns is sub-federal labor income tax

//shocks
rc =0.95; //persistence parameter on unused shock
// Benchmark calibration and steady state 
    nss = 0.33; //labor
    Rss = 1/beta-(1-delta);
    mcss = (epsilon-1)/epsilon;
    wss = mcss*(1-alpha)*(alpha*mcss/Rss)^(alpha/(1-alpha));
    k_hatss = ((alpha*mcss)/Rss)^(1/(1-alpha))*nss;
    kss = k_hatss;
    css = ((1-omega_f -omega_s)*(kss/nss)^alpha -delta*kss/nss)*nss;
    yss = (kss/nss)^alpha*nss;
    iss = 1/beta-1;
    lambdass = (css^gamma*(1-nss)^(1-gamma))^(-sigma)*gamma*((1-nss)/css)^(1-gamma);
    muss= beta*lambdass*Rss/(1+beta*(delta-1));
    Iss = delta*kss;
    uss = 1;
    vss =1;
    PIss = 1;
    PI_tildess = 1;
    Ass = lambdass*yss*mcss/(1-beta*theta);
    Dss = lambdass*yss/(1-beta*theta);
    gfss = omega_f*yss;
    gsss = omega_s*yss;
    Tfss = omega_t*yss;
    TTss = omega_t*yss;
    qss = 1;
    Lss = 1-nss;
    tau_nss = tau_nf*(1 - tau_nsss) + tau_nsss;
    tau_css = (gsss-Tfss-tau_ns*wss*nss)/css;
    gsDiscss = 0;
    SLGdeficitss = 0;
    tau_nsDiscss = 0;
    // tau_nfss = (gfss+Tfss-tau_k*kss*(Rss-delta))/(nss*wss);

model;
    L = 1-n;
    deta_hat = deta(-1);
    U = (c^gamma*L^(1-gamma))^(1-sigma)/(1-sigma) + deta_hat(+1)*U(+1);
    deta = beta + exc; //shock to the discount factor
    y = c+ I +gf + gs-gsDisc; //aggregate accounting identity
    y = (k_hat)^alpha*n^(1-alpha)/v; // aggregate production function
    k_hat = k; //capital utilization
    PI = ((1-theta)*(PI_tilde)^(1-epsilon)+theta)^(1/(1-epsilon)); //evolution of aggregate inflation
    v = (1-theta)*(PI/PI_tilde)^(epsilon)+theta*PI^(epsilon)*v(-1);//aggregate price dispresion  ??? - epsilon ???
    A = lambda*y*mc+theta*deta*((PI(+1))^epsilon)*A(+1);//auxiliary terms
    D = lambda*y + theta*deta*((PI(+1))^(epsilon-1))*D(+1);//auxiliary terms 
    A_hat = A(+1);
    D_hat = D(+1);
    PI_tilde = PI*(epsilon/(epsilon-1))*A_hat/D_hat;//reset inflation evolution
    w = mc*(1-alpha)*(k_hat/n)^alpha; //labor demand
    R = mc*alpha*(k_hat/n)^(alpha-1);//capital demand
    lambda*(1+tau_c+extau_c) = u_c;//marginal utility of income
    u_c = (c^gamma*(1-n)^(1-gamma))^(-sigma)*gamma*((1-n)/c)^(1-gamma);
    u_n = (c^gamma*(1-n)^(1-gamma))^(-sigma)*(1-gamma)*(c/(1-n))^gamma;
    tau_n = (tau_nf+extau_nf)*(1 - tau_ns - extau_ns) + (tau_ns+extau_ns);
    epsilon_w/(epsilon_w-1)*u_n = lambda*w*(1-tau_n);
    lambda = deta*lambda(+1)*(1+i(+1))*(PI(+1))^(-1);//Euler equation for bonds ??? i(+1) or just i ???
    mu = deta*(lambda(+1)*(R(+1)-tau_k*(R(+1)-delta))+mu(+1)*(1-delta-sigma_I/2*(I(+1)/k(+1)-delta)^2+sigma_I*(I(+1)/k(+1)-delta)*(I(+1)/k(+1)))); //FOC for capital
    lambda = mu*(1-sigma_I*(I/k-delta)); //FOC for investment
    k = I(-1)+(1-delta)*k(-1)-sigma_I/2*(I(-1)/k(-1)-delta)^2*k(-1); //capital accumulation
    log(gf) = (1-rho_g)*(log(omega_f)+log(STEADY_STATE(y)))+rho_g*log(gf(-1))+exGf; // federal government spending
    log(gs) = (1-rho_g)*(log(omega_s)+log(STEADY_STATE(y)))+rho_g*log(gs(-1))+exGs; //sub-federal government spending (when exogenous)
    dy = (y/STEADY_STATE(y))-1;
    dc = (c/STEADY_STATE(c))-1;
    dR= (R/STEADY_STATE(R))-1;
    dn = (n/STEADY_STATE(n))-1;
    dI = (I/STEADY_STATE(I))-1;
    dgf = (gf/STEADY_STATE(gf))-1;
    dgs = (gs/STEADY_STATE(gs))-1;
    dTf = (Tf/STEADY_STATE(TT))-1;
    gs-gsDisc = (tau_ns+tau_nsDisc+extau_ns)*w*n + (tau_c+extau_c)*c+Tf;         //Sub-federal budget constraint, with endogenous government spending response
    // gf = (tau_nf+extau_n)*w*n + tau_k*k(-1)*(R-delta)-Tf;       // Federal government budget constraint (redundant, missing transfers to HHs)
    Tf = max(TT+exTf,0);
    log(TT) = (1-rho_g)*(log(omega_t)+log(STEADY_STATE(y)))+rho_g*log(TT(-1)); // federal government transfer spending
    SLGdeficit = (tau_ns*STEADY_STATE(w)*STEADY_STATE(n) + STEADY_STATE(tau_c)*STEADY_STATE(c)+STEADY_STATE(TT))-((tau_ns+extau_ns)*w*n + (STEADY_STATE(tau_c)+extau_c)*c+Tf);
    gsDisc = psi*SLGdeficit;
    tau_nsDisc*w*n = nu*(1-psi)*SLGdeficit;
    
    // Policies

   i = max(Z,0); //Central bank policy
   // i = Z;    // Central bank policy without ZLB constraint
   Z = 1/beta*(PI(-1))^(phi_pi*(1-rho_r))*(y(-1)/STEADY_STATE(y))^(phi_y*(1-rho_r))*(beta*(1+i(-1)))^(rho_r)-1;
end;


initval;
    y   = yss;
    c = css;
    I = Iss;
    gf = gfss;
    gs = gsss;
    Tf = TTss;
    TT = TTss;
    tau_n = tau_nss;
    k_hat = kss;
    n = nss;
    v =vss;
    PI = 1;
    PI_tilde = 1;
    A = Ass;
    D = Dss;
    mc = mcss;
    w = wss;
    R = Rss;
    lambda = lambdass;
    mu = muss;
    i = iss;
    Z = iss;
    k = kss;
    deta = beta;
    A_hat = Ass;
    D_hat = Dss;
    deta_hat = deta;
    L = Lss;
    U = ((css^gamma*Lss^(1-gamma))^(1-sigma)/(1-sigma))/(1-beta);
    tau_c = tau_css;
    tau_nsDisc = tau_nsDiscss;
    gsDisc = gsDiscss;
    SLGdeficit = SLGdeficitss;

end;
options_.slowc=0.05; % you may lower this parameter if you want to stabilize the algorithm
resid
steady(solve_algo=2);
check;

// Hardcoding exTf and exGs shocks to mimick the persistence of the federal government shock:
xx = [0.00142678759906711; 0.00113684619417477; 0.000906557356008690; 0.000723384518670489; 0.000577520001118471; 0.000461257796381627; 0.000368521940412904; 0.000294508090019310; 0.000235408665374262; 0.000188200463495609; 0.000150479495947924; 0.000120331870203105; 9.62324082300486e-05; 7.69647590389069e-05; 6.15582644709761e-05; 4.92379464974921e-05; 3.93848127179866e-05; 3.15043023064801e-05; 2.52011715165118e-05; 2.01594843599190e-05];

// exogenous shocks: exc exGf exGs exTf extau_c extau_n
// shock the discount factor by 0.02 for 1 period, simulate and plot the results. 
shocks;
    var exc;
    periods 1:10;
    values 0.02;

    //var exGf;
    //periods 1;
    //values 0.04;

    //var exTf;
    //periods 1:20;
    //values(xx);

    //var exGs;
    //periods 1;
    //values 0.04;
    //periods 1:20;
    //values(xx);

end;
simul(periods=400);

nrep    = 20;
smpl    = 2+(1:nrep);

close all

subplot(331);plot(4*Z(smpl),'LineWidth',2);title('Z')
subplot(332);plot(dI(smpl),'LineWidth',2);title('Investment')
subplot(333);plot(dn(smpl),'LineWidth',2);title('Hours')
subplot(334);plot(4*(PI(smpl)-1),'LineWidth',2);title('Inflation')
subplot(335);plot(4*(i(smpl)),'LineWidth',2);title('Nominal Interest Rate')
subplot(336);plot(dy(smpl),'LineWidth',2);title('Output')
subplot(337);plot(dc(smpl),'LineWidth',2);title('Consumption')
subplot(338);plot(4*(deta(smpl)-beta),'LineWidth',2);title('Discount Factor Shock')
subplot(339);plot(4*(i(smpl)-(PI(smpl)-1)),'LineWidth',2);title('Real Interest Rate')