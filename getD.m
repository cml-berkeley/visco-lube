%------------------------------------------------------------
%   subroutine getmu
%------------------------------------------------------------
function D = getD(T,P)
% Calculate Binary Diffusion Coefficient of lubricant vapour at local temperature and
% pressure 
% Returns D [m2/s] 
% P [Pa] is ambient pressure
% T [K] is ambient temperature
% Model: Hirschfelder approximation used in Karis "Lubricants for the Disk
% Drive Industry"
M = 2;
M_air = 28.97/1000;
omega = 1.2; % collision integral
sigma_lube = 0.05*sqrt(M*1000); % [nm] vapor phase molecular diamater 
sigma_air = 0.315; % [nm] nitrogen molecular diamater 
sigma = (sigma_lube+sigma_air)/2;
D = (1.858e-4*sqrt(1/(M*1000)+1/(M_air*1000))).*(T.^1.5)./P./(sigma^2*omega);

end