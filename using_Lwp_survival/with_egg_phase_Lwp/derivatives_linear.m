%% BYOM function derivatives.m (the model in ODEs)
%
%  Syntax: dX = derivatives(t,X,par,c,glo)
%
% This function calculates the derivatives for the model system. It is
% linked to the script files <byom_debkiss_with_egg.html
% byom_debkiss_with_egg.m>. As input, it gets:
%
% * _t_   is the time point, provided by the ODE solver
% * _X_   is a vector with the previous value of the states
% * _par_ is the parameter structure
% * _c_   is the external concentration (or scenario number)
% * _glo_ is the structure with information (normally global)
%
% Time _t_ and scenario name _c_ are handed over as single numbers by
% <call_deri.html call_deri.m> (you do not have to use them in this
% function). Output _dX_ (as vector) provides the differentials for each
% state at _t_.
%
% * Author: Tjalling Jager
% * Date: November 2020
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_debkiss.html>

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

%% Start

function dX = derivatives(t,X,par,c,glo)

% NOTE: glo is now no longer a global here, but passed on in the function
% call from call_deri. That saves 20% calculation time!

%% Unpack states
% The state variables enter this function in the vector _X_. Here, we give
% them a more handy name.

% DO = X(1); % state 1 is dissolved oxygen, do not use if oxygen is a zvd 
WV = X(1); % state 1 is the structural body mass
% cR = X(2); % state 2 is the cumulative reproduction (not used in this function)
WB = X(3); % state 3 is the egg buffer of assimilates
SS = X(4);  % state 4 is survival function



%% Unpack parameters
% The parameters enter this function in the structure _par_. The names in
% the structure are the same as those defined in the byom script file.
% The 1 between parentheses is needed as each parameter has 5 associated
% values.

dV    = glo.dV;       % dry weight density of structure
delM  = glo.delM;     % shape correction coefficient (if needed)
% delM is used to allow Lwp as a parameter instead of WVp

sJAm  = par.sJAm(1);  % specific assimilation rate 
sJM   = par.sJM(1);   % specific maintenance costs 
WB0   = par.WB0(1);   % initial weight of egg
Lwp   = par.Lwp(1);   % shell length at puberty
yAV   = par.yAV(1);   % yield of assimilates on volume (starvation)
yBA   = par.yBA(1);   % yield of egg buffer on assimilates
yVA   = par.yVA(1);   % yield of structure on assimilates (growth)
kap   = par.kap(1);   % allocation fraction to soma
f     = par.f(1);     % scaled food level
fB    = par.fB(1);    % scaled food level for the embryo
Lwf   = par.Lwf(1);   % half-saturation length for initial food limitation
mu_emb= par.mu_emb(1);% mortality rate for embryos
mu_lar= par.mu_lar(1);% mortality rate for larvae
A     = par.A(1);     % lower oxygen threshold
B     = par.B(1);     % upper oxygen threshold

if glo.len ~= 0 % only if the switch is set to 1 or 2!
    % Translate shell length at puberty to body weight
    WVp = dV*(Lwp * delM).^3; % initial length to initial dwt
    WVf = dV*(Lwf * delM).^3; % initial length to initial dwt
end

%% Calculate the derivatives
% This is the actual model, specified as a system of ODEs.
% Note: some unused fluxes are calculated because they may be needed for
% future modules (e.g., respiration and ageing).


% Define stress level s that ranges from 0 to 1 between a range of DO values

DO = read_scen(-1,c,t,glo);     % use read_scen to derive actual exposure (DO) concentration at time t

% A = 0.5;                      % set upper and lower oxygen thresholds for stress
% B = 4.8;                      % but don't use if fitting as parameters

if DO < A                       % if oxygen below lower threshold, stress = 1
    s = 1;
elseif (DO < B) && (DO > A)     % if oxygen between threshold, linear relationship to stress
    s = 1-(DO - A) / (B - A);   
else                            % if oxygen above upper threshold, stress = 0
    s = 0;
end

% Model equations with (1-s) applied to decrease the parameter(s) of interest

L = (WV/dV)^(1/3); % volumetric length

% Select what to do with maturity maintenance
if glo.mat == 1
    sJJ = sJM * (1-kap)/kap; % add specific maturity maintenance with the suggested value - Eq. 4.7
else
    sJJ = 0; % or ignore it
end

WVb = WB0 *yVA * kap; % body mass at birth

if WB > 0 % if we have an embryo
    f = fB; % assimilation at different rate
else
    if WVf > 0
        f = f / (1+WVf/WV); % hyperbolic relationship for f with body weight
    end
end

JA = f * sJAm * (1-s) * L^2;          % assimilation - Eq. 2.7
JM = sJM * L^3;               % somatic maintenance - Eq. 2.8
JV = yVA * (kap*JA-JM);       % growth - Eq. 2.8

if WV < WVp                   % below size at puberty
    JR = 0;                   % no reproduction
    JJ = sJJ * L^3;           % maturity maintenance flux - Eq. 4.1
    % JH = (1-kap) * JA - JJ; % maturation flux (not used!) 
else
    JJ = sJJ * (WVp/dV);      % maturity maintenance - Eq. 4.2
    JR = (1-kap) * JA - JJ;   % reproduction flux - Eq. 4.3
end

% Starvation rules may override these fluxes - Equations 2.13 to 2.16, 
% and Equations 4.4 to 4.6 (adding maturity maintenance)
if kap * JA < JM      % allocated flux to soma cannot pay maintenance
    if JA >= JM + JJ  % but still enough total assimilates to pay both maintenances
        JV = 0;       % stop growth
        if WV >= WVp  % for adults ...
            JR = JA - JM - JJ; % repro buffer gets what's left - Eq. 4.4
        else
            % JH = JA - JM - JJ; % maturation gets what's left (not used!)
        end
    elseif JA >= JM   % only enough to pay somatic maintenance
        JV = 0;       % stop growth
        JR = 0;       % stop reproduction
        % JH = 0;     % stop maturation for juveniles (not used)
        JJ = JA - JM; % maturity maintenance flux gets what's left - Eq. 4.5
    else              % we need to shrink
        JR = 0;       % stop reproduction
        JJ = 0;       % stop paying maturity maintenance
        % JH = 0;     % stop maturation for juveniles (not used)
        JV = (JA - JM) / yAV; % shrink; pay somatic maintenance from structure - Eq. 4.6
    end
end

% We do not work with a repro buffer here, so a negative JR does not make
% sense; therefore, minimise to zero
JR = max(0,JR);

% Calcululate the derivatives
dWV = JV;             % change in body mass - Eq. 2.2
dcR = yBA * JR / WB0; % continuous reproduction flux - Eq. 2.12
if WB > 0 
    dS = - mu_emb * SS;         % embryo mortality
else 
    dS = -mu_lar * SS;
end
if WB > 0             % for embryos ...
    dWB = -JA;        % decrease in egg buffer - Eq. 2.1
else
    dWB = 0;          % for juveniles/adults, no change
end
if WV < WVb/4 && dWB == 0 % dont shrink below quarter size at birth
    dWV = 0; % simply stop shrinking ...
end

% dDO = 0;            % oxygen doesn't change, do not use because oxygen is forcing variable

dX = [dWV;dcR;dWB;dS];       % collect all derivatives in one vector
