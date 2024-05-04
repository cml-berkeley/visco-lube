%% Create a fine grid in appropriate domains based on simulation parameters
% Non-dimensionalize appropriately
% ============= Scales ==============
% Length (x,y)  : lrad
% Thickness (z)   : h0
% Temperature (C) : T0
% Time (s)      : lrad/ux

%% Inputs
lrad = 20e-9; % NFT Laser beam (on disk) full width half maximum [m]
dx_dim = 1e-9/1; 
dy_dim = 1e-9/1;
Lend_x = 3;
Lend_y = 3;
Lend_x_disk = 3;
Lend_y_disk = 3;
h0 = 1e-9; % Initial lubricant thickness [m]
T0_disk = 25; % Ambient disk temperature [deg C]
T0_slider = 25; % Ambient slider temperature [deg C]
p0 = 101325;
Tmax_disk = 500; % Maximum temperature on disk [deg C]
Tmax_slider = 300; % Maximum temperature on slider [deg C]
ux = 10; % Disk speed [m/s]
M = 2.7; % Lubricant molecular weight [kg/mol]
rho = 1600; % Lubricant density [kg/m3]
tf = 1e-6; % Final time [s]
dt = 1e-10; % Time step [s]
disp(['Times steps to cover distance dx: ',num2str(dx_dim/(ux*dt))])
if abs(dx_dim/(ux*dt)-1) > 1e-9
    disp('Error: Wrong Time step!')
end
%disp([num2str(round(tf/dt)),' time steps needed.'])
disp(['s_interval = ',num2str(round(tf/dt/200)),' for 200 recorded time steps.'])
s_interval = input('Choose a s_interval:                    ');
delx = ux*tf;
G0_lube = 0.5e6; % Lubricant thin-film shear modulus [MPa]
fh = 4e-9/h0; % Head-disk spacing, normalized by h0
lube_name = 2; % 1 for Zdol, 2 for Ztetraol

%% Grid generation - Control Volume center points

% Create x vector
x = (-Lend_x*lrad:dx_dim:Lend_x*lrad);% uniform grid
x_disk = (-Lend_x_disk*lrad:dx_dim:Lend_x_disk*lrad);% uniform grid

% Create y vector
y = (-Lend_y*lrad:dy_dim:Lend_y*lrad);% uniform grid
y_disk = (-Lend_y_disk*lrad:dx_dim:Lend_y_disk*lrad);% uniform grid

nx=length(x);
ny=length(y);
nx_d=length(x_disk);
ny_d=length(y_disk);

%% Pre-process: NON-DIMENSIONALIZE
% scaling
ts = lrad/ux; 
x=x/lrad; 
y=y/lrad;
x_disk=x_disk/lrad; 
y_disk=y_disk/lrad;

T0_disk=273.15+T0_disk; Tmax_disk=273.15+Tmax_disk;  % convert to Kelvin
T0_slider=273.15+T0_slider; Tmax_slider=273.15+Tmax_slider;  % convert to Kelvin
dT_disk=Tmax_disk-T0_disk; % [C] maximum temperature increase at laser center, used in non-dim below
dT_slider=Tmax_slider-T0_slider;

tf=tf/ts; dt=dt/ts;

c=6e-5; % [N/m/K] negative gradient of surface tension with temperature, -d(gamma)/dT