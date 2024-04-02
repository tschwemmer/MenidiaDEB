%% BYOM, byom_debkiss_with_egg.m
%
% * Author: Tjalling Jager
% * Date: November 2021
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_debkiss.html>
%
% BYOM is a General framework for simulating model systems in terms of
% ordinary differential equations (ODEs). The model itself needs to be
% specified in <derivatives.html derivatives.m>, and <call_deri.html
% call_deri.m> may need to be modified to the particular problem as well.
% The files in the engine directory are needed for fitting and plotting.
% Results are shown on screen but also saved to a log file (results.out).
%
% *The model:* standard DEBkiss model (see
% <http://www.debtox.info/book_debkiss.html>) following Jager et al (2013),
% <http://dx.doi.org/10.1016/j.jtbi.2013.03.011>. Parameters on length
% basis, and including the embryonic stage.
%
% *This script:* Example is for the pond snail at maximum food. This is a
% simulation, including the embryo phase. 
% 
%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

%% Initial things
% Make sure that this script is in a directory somewhere *below* the BYOM
% folder.

clear, clear global % clear the workspace and globals
global DATA W X0mat % make the data set and initial states global variables
global glo          % allow for global parameters in structure glo
diary off           % turn of the diary function (if it is accidentaly on)
% set(0,'DefaultFigureWindowStyle','docked'); % collect all figure into one window with tab controls
set(0,'DefaultFigureWindowStyle','normal'); % separate figure windows

pathdefine(0) % set path to the BYOM/engine directory (option 1 uses parallel toolbox)
glo.basenm  = mfilename; % remember the filename for THIS file for the plots
glo.saveplt = 0; % save all plots as (1) Matlab figures, (2) JPEG file or (3) PDF (see all_options.txt)

%% The data set
% Data are entered in matrix form, time in rows, scenarios (exposure
% concentrations) in columns. First column are the exposure times, first
% row are the concentrations or scenario numbers. The number in the top
% left of the matrix indicates how to calculate the likelihood:
%
% * -1 for multinomial likelihood (for survival data)
% * 0  for log-transform the data, then normal likelihood
% * 0.5 for square-root transform the data, then normal likelihood
% * 1  for no transformation of the data, then normal likelihood

% Total length in mm over time, at 24C. Time is days post fertilization.
% Sources: 16-110 dpf is Murray and Baumann 2020, 6dpf is Schwemmer unpub.
DATA{1} = [1        1
           6        5.1
           6        5.5
           16       8.9
           21       15.9
           21       15.7
           41       24.0
           56       30.0
           64       34.7
           89       48.7
           103      58.2
           110      55.6];    

% weight factors (number of replicates per observation)
W{1} = [10  
     60  
     60
     36  
     60
     60
     30  
     36 
     11  
     189  
     391];     

% Cumulative reproduction over time. (Number of eggs). Source: Concannon et
% al 2021. 
DATA{2} = [];

% the third state is for the egg buffer
DATA{3} = [1    1
           0   0.15
           6   0];

% weight factors (number of replicates per observation)
W{3} = [100
        100];

% add survival data, in proportion surviving at age in dpf. Sources: Cross
% et al 2019 and Murray et al 2017. 
DATA{4} = [1    1
           7    .88
           7    .95
           6    .74
           6    .62
           6    .666
           6    .902
           6    .576
           6    .562
           16   .32
           16   .72
           16   .41
           16   .33
           21   .5766
           21   .50668
           12   .30306
           21   .30274
           16   .83
           28   .2241
           47   .210654];

% weight factors (number of replicates per observation)
W{4} = [5
        5
        5
        3
        5
        5
        5
        5
        5
        5
        5
        3
        5
        5
        5
        5
        4
        4
        4];

% if weight factors are not specified, ones are assumed in start_calc.m

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat = [1       % the scenarios
         0.001   % initial shell length in mm
         0       % initial cumul. repro
         0.15    % initial weight of buffer in egg
         1.0];   % initial value for survival

% Put the position of the various states in globals, to make sure correct
% one is selected for extra things (e.g., to prevent shrinking in
% call_deri, and for using plot_tktd).
glo.locL = 1; % location of body size in the state variable list
glo.locR = 2; % location of cumulative reproduction in the state variable list

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 
% Global parameters as part of the structure glo
  
glo.delM  = 0.1066; % shape corrector (used in call_deri.m)
glo.dV    = 0.4;   % dry weight density (used in call_deri.m)
glo.len   = 2;     % switch to fit physical length (0=off, 1=on, 2=on and no shrinking) (used in call_deri.m)
glo.mat   = 1;     % include maturity maint. (0=off, 1=include)
% Note: if glo.len > 0 than the initial state for size in X0mat is length too!

% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.sJAm = [0.302    0 0 1e6]; % specific assimilation rate 
par.sJM  = [0.02143  0 0 1e6]; % specific maintenance costs - fix
par.WB0  = [0.15     0 0 1e6]; % initial weight of egg - fix
par.Lwp  = [115      0 0 1e6]; % total length at puberty - fix
par.yAV  = [0.8      0 0 1];   % yield of assimilates on volume (starvation)
par.yBA  = [0.95     0 0 1];   % yield of egg buffer on assimilates (repro)
par.yVA  = [0.35     1 0 1];   % yield of structure on assimilates (growth)
par.kap  = [0.7      0 0 1];   % allocation fraction to soma - fix
par.f    = [1        0 0 2];   % scaled food level - fix
par.fB   = [1        0 0 2];   % scaled food level, embryo - fix
par.Lwf  = [0        0 0 1e6]; % half-saturation total length - fix
par.mu_emb = [0.06624   0 0 1e6]; % mortality rate for embryos
par.mu_lar = [0.02809   0 0 1E6]; % mortality rate for larvae

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% used, based on the data set

glo.t   = linspace(0,400,100); % time vector for the model curves in days - MIDDLE ENTRY IS MAX TIME, LAST IS NO OF POINTS IN PLOT

% specify the y-axis labels for each state variable
glo.ylab{1} = 'total length (mm)';
glo.ylab{2} = 'cumulative reproduction (eggs)';
glo.ylab{3} = 'egg buffer (mg)';
glo.ylab{4} = 'survival';
% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = 'scenario '; % legend label before the 'scenario' number
glo.leglab2 = ''; % legend label after the 'scenario' number

prelim_checks % script to perform some preliminary checks and set things up
% Note: prelim_checks also fills all the options (opt_...) with defauls, so
% modify options after this call, if needed.

%% Calculations and plotting
% Here, the function is called that will do the calculation and the plotting.
% Options for the plotting can be set using opt_plot (see prelim_checks.m).
% Options for the optimsation routine can be set using opt_optim. Options
% for the ODE solver are part of the global glo. 
% 
% NOTE: for this package, the options useode and eventson in glo will not
% be functional: the ODE solver is always used, and the events function as
% well.

glo.stiff = [0 2]; % ODE solver 0) ode45 (standard), 1) ode113 (moderately stiff), 2) ode15s (stiff)
% Second argument is for normally tight (1), tighter (2), or very tight (3)
% tolerances. Use 1 for quick analyses, but check with 3 to see if there is
% a difference!

opt_optim.fit    = 1; % fit the parameters (1), or don't (0)
opt_optim.it     = 1; % show iterations of the optimisation (1, default) or not (0)
opt_plot.bw      = 0; % if set to 1, plots in black and white with different plot symbols
opt_plot.annot   = 0; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend
opt_plot.legsup  = 1; % if set to 1 omit legend, if 0 include legend (default)
opt_plot.limax   = 0; % if set to 1 limit axes to data, if 2 manually edit limit in calc_and_plot.m line 547

% optimise and plot (fitted parameters in par_out)
par_out = calc_optim(par,opt_optim); % start the optimisation
calc_and_plot(par_out,opt_plot); % calculate model lines and plot them

calc_proflik(par_out,'sJAm',opt_prof); % likelihood profiles
%calc_proflik(par_out,'kap',opt_prof);
%calc_proflik(par_out,'fB',opt_prof);
calc_ase(par_out,opt_ase); % calculate asymptotic standard errors and correlations between parameters

% % Code below demonstrates the use of the new plotting routine. However,
% % since we have only a single treatment here, the added benefit is not so
% % great.
% opt_tktd.repls = 0; % plot individual replicates (1) or means (0)
% opt_tktd.preds = 1; % set to 1 only plot predictions from X0mat without data
% plot_tktd(par_out,opt_tktd,[]);
% % Leave the options opt_conf empty to suppress all CIs for these plots.
% % Plotting CIs requires a sample as saved by the various methods available
% % in BYOM (see examples directory).

%% Other files: derivatives
% To archive analyses, publishing them with Matlab is convenient. To keep
% track of what was done, the file derivatives.m can be included in the
% published result.
% 
% <include>derivatives.m</include>

%% Other files: call_deri
% To archive analyses, publishing them with Matlab is convenient. To keep
% track of what was done, the file call_deri.m can be included in the
% published result.
%
% <include>call_deri.m</include>
  

