% Author: Luis Badesa - adapted by Samuel Tharp (ERCOT-style aggregation, RELAXED)

%%
clearvars
clc

%% =======================
%  Input data (ERCOT-like)
%  =======================

% --- Study controls (wind penetration sweep) ---
D = 24e3;                 % MW demand (single-interval snapshot)
alpha_wind = 0.10;        % wind availability as fraction of demand (e.g., 0.10/0.30/0.50/0.70)
InputData.P_Wind = alpha_wind * D;   % MW available wind
InputData.H_const_Wind = 0;          % s (as modeled here)

% --- ERCOT frequency + contingency assumptions ---
InputData.f_0       = 60;    % Hz (ERCOT)
InputData.PLossMax  = 2750;  % MW (ERCOT largest loss assumption)
InputData.nadir_req = 0.8;   % Hz (as used in Badesa formulation)
InputData.Rocof_max = 1;     % Hz/s (placeholder; adjust to your ERCOT assumption)

% --- Reindexed so Td is monotone (fastest first) ---
% Cluster order now: [Battery/Inverter, Peaker, Steam/Coal, CCGT]
InputData.Td = [2 7 10 10]'; % seconds

% =========================
%  Generator cluster design
% =========================
% Clusters (4) in the SAME order as Td above:
%  1) Battery / inverter fast response (zero inertia, fast FR capability)
%  2) SCGT/Peaker (low inertia, very flexible, expensive energy)
%  3) Steam/Coal-like synchronous (high inertia, inflexible)
%  4) CCGT (moderate inertia, flexible)

InputData.num_gen = [200 40 20 60]';   % number of units per cluster (aggregate count)

% No-load / commitment-like cost per online unit (in $)
NLHR = [0 50 800 200]';                % $  (battery ~0, peaker low, steam high, CCGT medium)

% Energy marginal costs (in $/MWh)
HRS  = [999 120 35 45]';               % $/MWh (battery high to discourage energy dispatch)

% Small bid/penalty for providing FR (in $/MW)
FR_bids = [0.001 0.001 0.001 0.001]';  % keep tiny unless you want to discourage over-procurement

% Generation limits for individual units in each cluster: [Pmin, Pmax] (MW)
InputData.Gen_limits = [0   10;        % 1) Battery (many small units)
                        0   200;       % 2) Peaker
                        300 800;       % 3) Steam/Coal-like
                        100 600];      % 4) CCGT

% Inertia constants (s) by cluster
InputData.H_const = [0.0 2.0 6.0 4.0]';  % battery 0, peaker low, steam high, CCGT medium

% FR capability parameters (dimension = num_Clusters x 1)
MaxFRpercentage = [1.00 0.90 0.50 0.70]'; % fraction of Pmax available as FR capacity per unit
TapperSlope     = [1.00 1.00 0.80 0.90]'; % tapering of FR vs headroom

% =========================
%  Largest-loss unit proxy
% =========================
HRS_Nuclear = 15;              % $/MWh (energy cost for the largest-loss unit proxy)
InputData.H_const_Nuclear = 5; % s
max_deloading = 0;             % MW

%% Create optimisation problem
num_Clusters = length(NLHR);

x  = sdpvar(num_Clusters,1);   % MW produced by each cluster
y  = sdpvar(num_Clusters,1);   % RELAXED commitment (continuous # units online)
FR = sdpvar(num_Clusters,1);   % MW FR provided by each cluster

x_Nuclear       = sdpvar(1);   % MW produced by largest-loss unit proxy
x_WindCurtailed = sdpvar(1);   % MW curtailed wind

% FR capacity per unit in each cluster (MW)
InputData.FR_capacity = MaxFRpercentage .* InputData.Gen_limits(:,2);

% Variable bounds
Bounds = [ 0 <= y  <= InputData.num_gen, ...
           0 <= x  <= InputData.num_gen .* InputData.Gen_limits(:,2), ...
           0 <= FR <= InputData.num_gen .* InputData.FR_capacity ];

Constraints_Gen_limits = [ ...
    y .* InputData.Gen_limits(:,1) <= x <= y .* InputData.Gen_limits(:,2), ...
    InputData.PLossMax - max_deloading <= x_Nuclear <= InputData.PLossMax, ...
    0 <= x_WindCurtailed <= InputData.P_Wind ];

% FR limits
FR_limits = [ ...
    0 <= FR <= y .* InputData.FR_capacity, ...
    FR <= TapperSlope .* (y .* InputData.Gen_limits(:,2) - x) ];

% Total inertia proxy
H_total = ((InputData.H_const .* y)' * InputData.Gen_limits(:,2)) / InputData.f_0 + ...
          (InputData.H_const_Wind * (InputData.P_Wind - x_WindCurtailed)) / InputData.f_0;

% Power balance
Power_balance = sum(x) + x_Nuclear + (InputData.P_Wind - x_WindCurtailed) == D;

% DV_PLoss equals the largest-loss unit output in this case study
DV_PLoss = sdpvar(1);
Bounds = [Bounds, InputData.PLossMax - max_deloading <= DV_PLoss <= InputData.PLossMax];
Deload_constraint = DV_PLoss == x_Nuclear;

% Nadir constraints
[DV_Total_FR_atTd, Inertia_term, PLoss_term, FR_term, Bounds, Nadir_constraints] = ...
    setNadir(InputData, FR, H_total, DV_PLoss, Bounds);

% Quasi-steady-state FR requirement and RoCoF
qss = sum(FR) >= DV_PLoss;
Rocof_constraint = H_total >= DV_PLoss / (2 * InputData.Rocof_max);

% Assemble all constraints
Constraints = [Bounds, ...
               Constraints_Gen_limits, ...
               FR_limits, ...
               Power_balance, ...
               Deload_constraint, ...
               Nadir_constraints, ...
               qss, ...
               Rocof_constraint];

%% Solve optimization
Cost = NLHR' * y + HRS' * x + HRS_Nuclear * x_Nuclear + FR_bids' * FR;
Objective = Cost;

options = sdpsettings('solver','mosek');
sol = optimize(Constraints, Objective, options);

ObjectiveFunction = value(Objective);

%% Results
Nadir_check = zeros(length(InputData.Td),1);
Nadir_difference = zeros(length(InputData.Td),1);
when_FR_delivered = zeros(length(InputData.Td),1);

for i = 1:length(InputData.Td)
    Nadir_check(i) = value(Inertia_term(i)) * value(FR_term(i)) ...
        >= value(PLoss_term(i))^2 / (4 * InputData.nadir_req);

    Nadir_difference(i) = value(Inertia_term(i)) * value(FR_term(i)) ...
        - value(PLoss_term(i))^2 / (4 * InputData.nadir_req);

    when_FR_delivered(i) = value(DV_Total_FR_atTd(i)) >= value(DV_PLoss);
end

Cost_k = value(Objective) * 1e-3
Generators_MW = value(x)
Commitment_units = value(y)
FR_MW = value(FR)
Wind_used_MW = InputData.P_Wind - value(x_WindCurtailed)
Wind_curtail_MW = value(x_WindCurtailed)
LargestLossUnit_MW = value(x_Nuclear)

when_FR_delivered
InputData.Td
