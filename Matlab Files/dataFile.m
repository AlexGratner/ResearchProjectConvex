%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data required for running the optimizations                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cm = 91350;                         %Machine type specific constant, torque [N/m^2.5]
Cmj = 1056;                         %Machine type specific constant, intertia [kg/m^3]
Cgr = 1.2;                          %Outer radius by reference radious [~]
Jload=1.1;                          %Inertia load [kg*m^2]
inertia_ratio = 10/9;               %(Jm+J0)/Jm
Trms = 60.9338/0.95;                %RMS Torque [Nm]
Tpeak = 221.937/0.95;               %Peak torque [Nm]
delta_hlim = 1200e6;                %Maximum allowed Hertzian pressure
SF = 1;                             %Safety factor
delta_hmax = delta_hlim/SF;         %Maximum allowed flank pressure [Pa]
lb = 0.5;                           %Form factor, lower limit [~]
ub = 5;                             %Form factor, upper limit [~]

k = 4e10*Cgr^2*Tpeak/(delta_hmax^2);%Introduced variable for gearbox [...]
k_inv = 1/k;                        %Inverse of k [k^-1]

G = 80e9;                           %Shaft shear modulus
damp = 20;                          %Shaft damping
