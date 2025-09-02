%% Robust and Nonlinear Control, EEN050:
% Template for assignment 1, 2, and 3
%-------------------------------------------------
%-------------------------------------------------
%Initialization
clear all;
close all;
clc
% READ THIS:
% - Run this section first and do not overwrite any of the variables here
% besides (if necessary) input names and output names for the ss-objects.

%Low fidelity F16 longitudinal model from Aircraft Control and Simulation (B.L.Stevens-F.L.Lewis) pp. 156
%The linearized dynamic of the airplane F16
A_n=[-0.127 -235 -32.2 -9.51 0.314;-7E-4 -0.969 0 0.908 -2E-4;0 0 0 1 0; 9E-4 -4.56 0 -1.58 0;0 0 0 0 -5];
%States
%[V(ft/s) speed, alpha(rad) angle of attack, theta(rad) pitch angle, q(rad/s) pitch rate, T(lb) engine_power]'
B_n=[0 -0.244;0 -0.00209; 0 0;10 -0.199; 1087 0];
%Control inputs
%[thrust (N); elevator_deflection(rad)]'
C_n=[0 57.3 0 0 0;0 0 0 1 0;0.0208 15.2 0 1.45 0];
%Measured outputs
%[alpha(deg); q(rad/s); normal_acceleration(ft/s^2) ]
D_n=[0 0;0 0;0 0.033];
%Note, elevator deflection has a direct effect on vertical/normal
%accelleration
Gn=ss(A_n,B_n,C_n,D_n);
Gn.InputName = 'utilde';
Gn.OutputName = 'y';
% Actuator dynamics (nominal):
GT = tf(1,[1/(2.5*10) (1/2.5+1/10) 1]); % Thrust
Ge = tf(1,[1/25 1]);                    % Elevator deflection
Ga  = ss([GT, 0;0, Ge]);                % Actuator dynamics
Ga.InputName = 'u';
Ga.OutputName = 'ytilde';

% Disturbance model:
dryden =  tf([0.9751 0.2491],[1 0.885 0.1958]);
Wd = ss([0;dryden]);
Wd.InputName = 'd';
Wd.OutputName = 'Wd';

% Noise filter Wn:
wn1 = rad2deg(0.001);
wn2 = 0.001;
ft2m = 0.3048; % feet/meter
wn3 = 0.001/ft2m;
Wn = ss(diag([wn1,wn2,wn3]));
Wn.OutputName = 'Wn';
Wn.InputName = 'n';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LQG design %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assemble the physical system
G = Gn*[Ga Wd];

[A,B,C,~] = ssdata(G); % Retrive model matrices
% Number of states, control inputs, process noise, measurement noise
nx = length(A); nu = 2; nw = 1; ny = 3;

Glqi = ss(A,B(:,1:2),C(1,:),0); % ss-object is used for the LQI design

% LQI design - weigthing matrices
Q = blkdiag(eye(nx),10000);
R = eye(nu);

% Compute the feedback gain K
[K,~,~] = lqi(Glqi,Q,R);

% Kalman filter design w. 'kalman'
Gkalman = ss(A,B(:,3),C,0);% ss-object for filter design
QN = eye(nw);  % Process noise covariance
RN = eye(ny);  % Measurement noise covariance
[~,L,~] = kalman(Gkalman,QN,RN); % Compute the filter

% Split the gain into the state feedback gain and the integral gain
Kx = K(:,1:end-1);  % state feedback
Ki = K(:,end);      % integral gain
Bu = B(:,1:2);      % this is the part of B that multiplies with u

% Construct the system matrices
Ac = [A-L*C-Bu*Kx -Bu*Ki;-C(1,:) 0];
Bc = [L zeros(nx,1);0 0 0 1];
Cc = -K;
Dc = 0;

% Represent the LQG controller as a dynamical system
% Its input is [ytilde;r]
% Its output is [u]
LQG = ss(Ac,Bc,Cc,Dc);
LQG.InputName={'ytilde(1)','ytilde(2)','ytilde(3)','r'};
LQG.Outputname = {'u(1)','u(2)'};
'You can ignore the above name conflict'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% A1/Ex1 - Open loop analysis
clc
G_F16 = Gn*Ga;                              % Airplane model

% Is the airplane model G_F16 is stable?


% Is the airplane model G_F16 minimal order?


% Plot the singular values of the airplane model G_F16 (see the command
% 'sigma')

% Compute the H-infinity norm of the airplane model
% (This can be done with the command 'hinfnorm'.)

% Compute the H-2 norm of the airplane model 
% (The easiest way to do this is with the command 'norm'.)
% example: norm(ss-object of interest,2)

% Compute the H-2 norm of 'Gn' 

%% A1/Ex3 - Uncertain actuator dynamics
% MATLAB-functions to use, 'tf', 'ureal', 'usample', 'step', 'sigma',
% 'blkdiag'.
clc
% Define the uncertainty parameters with 'ureal'
%t1   =  ureal('t1',2.5,'percentage',30);
%t2   =  ureal(...);
%t3   =  ureal(...);
%gT   =  ureal(...);
%ge   =  ureal(...);


% Define the uncertain dynamics for thrust
%GTu=tf(...);                % Thrust

% Define the uncertain dynamics for elevator deflection
%Geu=...                                 % Elevator deflection

% Define the uncertain actuator dynamics (block diagonal)
%Gau = blkdiag(GTu,Geu);
%Gau.InputName = 'u';
%Gau.OutputName = 'ydelta';

% Create 100 samples of Gau1 and Gau2 with 'usample'
samples = 100; % number of samples
% it takes time to generate 100 samples so feel free to change this
% the value of 'samples' when you are figuring things out.

%Geu_samples = usample(...); 
%GTu_samples = ...
%
% Plot the step responses and bode diagrams (sigma plot)
% Feel free to use the following plotting routine:
% figure(1)
% subplot(2,2,1)
% step(GTu_samples)
% hold on
% step(GTu.NominalValue,'red')
% title('Thrust: step response')
% legend('random samples','nominal')
% 
% subplot(2,2,2)
% step(Geu_samples)
% hold on
% step(Geu.NominalValue,'red')
% title('Elevator deflection: step response')
% legend('random samples','nominal')
% 
% subplot(2,2,3)
% sigma(GTu_samples)
% hold on
% sigma(GTu.NominalValue, 'red')
% title('Thrust: singular values')
% legend('random samples','nominal')
% 
% subplot(2,2,4)
% sigma(GTu_samples)
% hold on
% sigma(GTu.NominalValue, 'red')
% title('Elevator deflection: singular values')
% legend('random samples','nominal')
% hold off
%% Ex2/A2 - (1)
% Simulation of LQG controller w. nominal actuator dynamics
% Consult the documentation for 'lsim' and 'connect' to better understand
% the code
clc
% Simulation parameters:
N = 1000;                  % number of time steps in simulation
T = linspace(0,50,N);      % time vector
flag_noise = 1;            % set to zero to remove noise
flag_x0 = 1;               % set to zero to initialize system at the origin

% Inputs:
r = zeros(N,1);            % reference signal
r(1:200)=-0.5; r(201:400)=1; r(401:600)=5; r(601:800)=-2; r(801:end)=0;
noise = randn(N,4);        % disturbance and measurement noise
U = [r noise*flag_noise];  % input vector

% Define the closed loop system using 'connect':
% Provide appropriate input/output names for each ss-object
Wd.InputName = 'd';
Wd.OutputName = 'Wd';

Ga.InputName = 'u';
Ga.OutputName = 'ydelta';

Wn.InputName = 'n';
Wn.OutputName = 'Wn';

Gn.InputName = 'utilde';
Gn.OutputName = 'y';

LQG.InputName = {'ytilde(1)','ytilde(2)','ytilde(3)','r'};
LQG.OutputName = 'u';

% Summation blocks:
Sum1 = sumblk('utilde = ydelta+Wd',2);
Sum2 = sumblk('ytilde = y+Wn',3);

% Choose inputs and outputs for the resulting closed loop system
inputs = {'r','n','d'};
outputs = {'y(1)'};     % y(1) corresponds with angle of attack [deg]

% Define the closed loop system:
LQG_clp = connect(Ga,Gn,Wd,Wn,LQG,Sum1,Sum2,inputs,outputs);

% Number of states in the closed loop system
nx_LQG = length(LQG_clp.A);
x0 = randn(nx_LQG,1)*flag_x0;   % initial state

% Simulating and plotting
Y_LQG = lsim(LQG_clp,U,T,x0);
figure(2)
plot(T,Y_LQG,T,r,'r--')
legend('angle of attack','reference signal')
title('Reference tracking w. LQG on the angle of attack')
ylim([-10, 10])

%% A1/Ex3 - (2)
% Closed loop simulation w. LQG and uncertain actuator dynamics
% Repeat the above simulation with the ss-object 'Ga' replaced by a sample
% of the uncertain ss-object for the actuator dynamics 'Gau'(c.f A1/Ex2).
% The input and output names of the samples actuator dynamics should
% correspond with the input and output names of 'Ga'.

% Take a random sample of the uncertain actuator dynamics
%Ga_random = usample(...);
%Ga_random.InputName = 'u';
%Ga_random.OutputName = 'ydelta';

% Close the loop with Ga_random instead of Ga
clp_LQGu = connect(Ga_random,Gn,LQG,Wn,Wd,Sum1,Sum2,inputs,outputs);
% simulate and plot
Y_LQGu = lsim(clp_LQGu,U,T,x0);
figure(3)
plot(T,Y_LQGu,T,r,'r--')
legend('angle of attack','reference signal')
title('Reference tracking w. LQG on the angle of attack (random actuator dynamics)')
ylim([-10, 10])

%% A2/Ex1
clc
% Notes:
% Be mindful of the input/output dimensions when defining the filters.
% Feel free to use 'makeweight' when apropriate
% 
% Wn and Wd are defined in the first code section of this file.
% (Be aware of their input and output names!)
% To define the filter Wm (consult the documentation for 'ucover').
%
% Use 'connect' to contstruct the 'P' matrix.
% For convenience, order the inputs/outputs as follows
% Input: [udelta;r;n;d;u]
% Output: [ydelta;z;v] (v=[y+n;r])
%
% To compute the Hinf controller use 'hinfsyn'
% To compute the H2 controller use 'h2syn'

% Define the filters Wra, We, Wp, Wu


% Define WmT, Wme, and Wm (read the documentation for 'ucover')
% Feel free to use the samples taken in A1/Ex2


%WmT = ...             % Thrust
%Wme = ...             % Elevator deflection
%Wm = blkdiag(WmT,Wme);

% Provide appropriate input/output names


% Define the summation blocks (there are three):

% Define the appropriate inputs and outputs (the order matters!)
%inputs = ...
%outputs = ...

%P = connect(...)

%% A2/Ex2
% Compute the H-infinity controller

% plot the singular values
figure(4)

%% A2/Ex3
% Compute the H2-controller

% plot the singular values
figure(5)



%% A3/Ex1



%% A3/Ex2 Simulation
% - Construct with 'connect' the  closed loop system between the LQG, Hinf,
%   and H2 controllers and the open loop plant.
% The open loop plant consists of:
% - Ga (nominal value or a sample of the uncertain actuator dynamics)
% - Wn, Wd, and, Gn
% - The delta block and the filters Wra, We, Wp, Wu, Wm are not included!

% Take a random sample of the uncertain actuator dynamics

% Define the appropriate summation blocks (there are two)


% Choose the apropriate inputs and outputs:


% Define the closed loop systems:
% LQG closed loop system
%LQG_clp = connect(...);

% H-infinity closed loop system
%Hinf_clp = connect(...);

% H-2 closed loop system
%H2_clp = connect(...);


% Feel free to use the following plotting routine:
% Simulation parameters:
% N = 1000;                  % number of time steps in simulation
% T = linspace(0,50,N);      % time vector
% flag_noise = 1;            % set to zero to remove noise
% flag_x0 = 1;               % set to zero to initialize system at the origin
% % Inputs:
% r = zeros(N,1);            % reference signal
% r(1:200)=-0.5; r(201:400)=1; r(401:600)=5; r(601:800)=-2; r(801:end)=0;
% noise = randn(N,4);        % Disturbance and measurement noise
% U = [r noise*flag_noise];  % Input vector
% 
% % Initial state(s)
% % Number of states is the LQG closed loop system
% nx_LQG = length(LQG_clp.A);
% % Number of states in the H_(2,infinity) closed loop systems
% nx_H = length(Hinf_clp.A);
% 
% x0_LQG = randn(nx_LQG,1)*flag_x0;   % initial state (LQG)
% x0_H = randn(nx_H,1)*flag_x0;       % initial state (Hinf,H2)
% 
% % Simulate:
% Yinf =  lsim(Hinf_clp,U,T,x0_H);
% Y2   =  lsim(H2_clp,U,T,x0_H);
% YLQG =  lsim(LQG_clp,U,T,x0_LQG);
% 
% % Feel free to use the following plotting routine.
% figure(8)
% plot(T,Yinf,T,Y2,T,YLQG,T,r,'r--')
% title('Reference tracking on angle of attack [deg]')
% legend('Hinf','H2','LQG','reference')
% ylim([-3,6])