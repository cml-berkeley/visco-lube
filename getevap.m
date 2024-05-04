%------------------------------------------------------------

function mdot=getevap(T,pd,M,rho,p_vapor,lube_name)

if lube_name == 1
%% Zdol

% Calculate thin film evaporation flux at local temperature [K] and 
% thickness [m] for specified molecular weight [kg/mol]
% All input quantities dimensional: K,Pa,Pa,M
% OUTPUT: 
% mdot ::  Mass flux [kg/m^2/s]
% Evalutate bulk vapor pressure according to Karis Ch 22 (2009)
% Thin film effect is the Gibbs-Duhem equation L Wu (2007)
R=8.314;       % universal gas constant [J/mol*K]
hP=6.626e-34;  % Plank's constant [J*s]
N=6.022e23;    % Avogadro's number [molecules/mol]
Pa = 101325;   % [Pa] ambient pressure of 1 atm
sigma_lube = 0.05*sqrt(M*1000); % [nm] vapor phase molecular diamater 
Sliq = 107; % [J/(mol*K)] liquid entropy for Zdol, indep of M
% vapor phase translational entropy for oils (Sackur-Tetrode eqn)
Svap_trans = R*(5/2+log((2*pi)^1.5/hP^3*(R*T).^2.5/N^4/Pa)...
    +3/2*(log(M))); % [J/(mol*K)] Eqn 22.16
% vapor phase rotational entropy Eqn 22.17
I = (sigma_lube*1e-9)^2*M/16/N^2; % [kg*m^2*mol]
a = 1; % number of independent rotation axes
q = 1; % degeneracy, internal rotations not included
Svap_rot = R*(1+log(1/(pi*q)*(8*pi^3*I*R*T/hP^2).^(a/2))); % [J/(mol*K)]
% vaporization entropy dSvap = Svap - Sliq
dSvap = Svap_trans + Svap_rot - Sliq;
% vaporization activation energy
dEvap = 50e3+29*M*1000;  % [J/mol]
% pure component BULK vapor pressure, Eqn 22.15
Pvap_bulk = Pa*exp(dSvap/R)*exp(-1).*exp(-dEvap./(R*T));   
% Mass flux theoretical maximum according to kinetic theory (collision
% theory) Hertz-Knudsen-Langmuir equation. alpha = 1
% Thin film effect: Gibbs-Duhem equation equating chemical potential of
% vapor and thin film [kg/m^2/s]

% mdot = Pvap_bulk.*sqrt(M./(2*pi*R*T)).*exp(M./(rho*R*T).*(-pd-pl));

Pvap_film = Pvap_bulk.*exp(M./(rho*R*T).*(-pd));
mdot = (Pvap_film-p_vapor).*sqrt(M./(2*pi*R*T));

elseif lube_name == 2
%% ZT from Jones et. al.
R = 8.314; % universal gas constant [J/mol*K]
P0 = 685415*133.32; % [J/mol]
dH = 68.112e3; 
Pvap_bulk = P0.*exp((-dH/R)./T);
Pvap_film = Pvap_bulk.*exp(M./(rho*R.*T).*(-pd));
mdot = (Pvap_film-p_vapor).*sqrt(M./(2*pi*R.*T));

else
    mdot = zeros(size(T));
end

end