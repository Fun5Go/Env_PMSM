%% ------- SYSTEM PARAMETERS ------- %%
Vdc = 1200; % DC-Link Voltage       [V]
Fs  = 10e3; % PWM Carrier Frequency [Hz]
Ts  = 1/Fs; % Control Sampling Time [s]

%% ------- MOTOR PARAMETERS ------- %%
% Electrical
p     = 4;		    % Pairs of poles        [-]
Rs    = 29.0808e-3; % Stator resistance     [Ohm]
Ld    = 0.91e-3;    % Inductance in d-axis  [H]
Lq    = 1.17e-3;    % Inductance in q-axis  [H]
lambda_PM = 0.172312604; % Flux-linkage due to permanent magnets [Wb]
% Mechanical
J     = 0;	        % Inertia               [kg*m2]
F     = 0;          % Friction coefficient

%% ------- REFERENCE TORQUE ------- %%
% Reference Torque for non-saturated operation
Tref = 60;

%% ------- DESIGN OF PI ------- %%
Kp_Id = 3.9932;             % Proportional gain in d-frame
Ki_Id = Kp_Id*(1-0.8776)/Ts;% Integrative gain in d-frame
Kp_Iq = 3.9932;             % Proportional gain in d-frame
Ki_Iq = Kp_Iq*(1-0.8776)/Ts;% Integrative gain in d-frame

z = tf('z',Ts);
s = tf('s');

% -----  TORQUE LOOP ----- %
Gd_c = 1/(Rs + Ld*s);   % Continuous plant TF in d-frame
Gd   = c2d(Gd_c,Ts)/z;  % Discrete plant TF in d-frame (including delay)
Cd   = Kp_Id + Ki_Id*Ts/(z-1);  % Controller in d-frame

Gq_c = 1/(Rs + Lq*s);   % Continuous plant TF in q-frame
Gq   = c2d(Gq_c,Ts)/z;  % Discrete plant TF in q-frame (including delay)
Cq   = Kp_Iq + Ki_Iq*Ts/(z-1);  % Controller in q-frame

% sisotool(Gd,Cd);
% sisotool(Gq,Cq);
% % Tranfer functions of the closed loop in d and q frame
% H2d = Gd*Cd/(1+Gd*Cd);
% H2q = Gq*Cq/(1+Gq*Cq);
% % Continuous tranfer functions of the closed loop in d and q frame
% H2d_c = d2c(H2d,'tustin');
% H2q_c = d2c(H2q,'tustin');
