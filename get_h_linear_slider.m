function [h_lube_slider_try,m_dot_slider] = get_h_linear_slider(p_vapor,Text_x_lube_slider,Text_y_lube_slider,dText_x_dt,dText_y_dt,T_lube_slider,G0_lube,n_h_lube_slider,h_lube_slider_old,dt,ts,h0,c,dT_slider,lrad,dx_ND,dy_ND,dx_dim,dy_dim,M,rho,lube_name,nx,ny,ux)

f = zeros(nx,ny);

%% Material Properties
if lube_name == 1
    % Zdol
    lambda_lube = getlambda_ZD(T_lube_slider)/2/pi;
    n_mu_lube = 20/13*G0_lube.*lambda_lube;
    G_lube = 20/13*G0_lube*ones(size(n_h_lube_slider));
elseif lube_name == 2
    % Ztetraol
    lambda_lube = getlambda_ZT(T_lube_slider)/2/pi;
    n_mu_lube = 20/13*G0_lube.*lambda_lube;
    G_lube = 20/13*G0_lube*ones(size(n_h_lube_slider));
else
    problem = 1;
end

% Slip Length (nm)
b = 1;

%% Viscous Shear Stress

A_x = ( Text_x_lube_slider(3:end,2:end-1).*n_h_lube_slider(3:end,2:end-1).^2./n_mu_lube(3:end,2:end-1).*(1/2 + b./n_h_lube_slider(3:end,2:end-1)) ...
    -  Text_x_lube_slider(1:end-2,2:end-1).*n_h_lube_slider(1:end-2,2:end-1).^2./n_mu_lube(1:end-2,2:end-1).*(1/2 + b./n_h_lube_slider(1:end-2,2:end-1)) )/(2*dx_ND);
A_x = A_x.*(c*dT_slider*h0^2/lrad^2);

A_y = ( Text_y_lube_slider(2:end-1,3:end).*n_h_lube_slider(2:end-1,3:end).^2./n_mu_lube(2:end-1,3:end).*(1/2 + b./n_h_lube_slider(2:end-1,3:end)) ...
    -  Text_y_lube_slider(2:end-1,1:end-2).*n_h_lube_slider(2:end-1,1:end-2).^2./n_mu_lube(2:end-1,1:end-2).*(1/2 + b./n_h_lube_slider(2:end-1,1:end-2)) )/(2*dy_ND);
A_y = A_y.*(c*dT_slider*h0^2/lrad^2);

f(2:end-1,2:end-1) = f(2:end-1,2:end-1) - A_x - A_y;

%% Elastic Shear Stress

B_x = ( dText_x_dt(3:end,2:end-1)./(2.*G_lube(3:end,2:end-1)).*(n_h_lube_slider(3:end,2:end-1)).^2 ...
    - dText_x_dt(1:end-2,2:end-1)./(2.*G_lube(1:end-2,2:end-1)).*(n_h_lube_slider(1:end-2,2:end-1)).^2 )./(2*dx_ND);
B_x = B_x.*(c*dT_slider*h0^2/lrad^2/ts);

B_y = ( dText_y_dt(2:end-1,3:end)./(2.*G_lube(3:end,2:end-1)).*(n_h_lube_slider(2:end-1,3:end)).^2 ...
    - dText_y_dt(2:end-1,1:end-2)./(2.*G_lube(1:end-2,2:end-1)).*(n_h_lube_slider(2:end-1,1:end-2)).^2 )./(2*dy_ND);
B_y = B_y.*(c*dT_slider*h0^2/lrad^2/ts);

f(2:end-1,2:end-1) = f(2:end-1,2:end-1) - B_x - B_y;

%% Evap
A_ham = 5.3e-21;

% A_ham = 3e-20/6/pi;

% % temp dependent hammaker
% A_ll = 3.78e-20*(1-1.72e-3*(T_lube_slider-273.15));
% A_ss = 1.16e-19;
% A_ham = (sqrt(A_ll*A_ss)-A_ll)/6/pi;


d0 = 0.3e-9;
p_disj = A_ham./(n_h_lube_slider*h0 + d0).^3;

E = getevap(T_lube_slider,p_disj,M,rho,p_vapor,lube_name)/rho;
f = f - E;
m_dot_slider = E*rho;

%% DP Elastic

dp_dh = (3*A_ham./(n_h_lube_slider*h0 + d0).^4);
dp_dh_old = (3*A_ham./(h_lube_slider_old*h0 + d0).^4);

dh_dt = (n_h_lube_slider-h_lube_slider_old)*h0/dt/ts;

Diff = ((1/3./G_lube/dt/ts).*(n_h_lube_slider*h0).^3.*dp_dh).*h0 + (dh_dt.*dp_dh.*(n_h_lube_slider.*h0).^2./2./G_lube).*h0;
Diff_old = ((1/3./G_lube/dt/ts).*(n_h_lube_slider*h0).^3.*dp_dh_old).*h0;

%% DP Viscous
dp_dh = (3*A_ham./(n_h_lube_slider*h0 + d0).^4);
Diff = Diff + (1./n_mu_lube.*(n_h_lube_slider*h0).^3.*dp_dh.*(1/3 + b./n_h_lube_slider)).*h0;
%Diff_old = 0.*Diff;

% No Slider evolution
%Diff = 0.*Diff;

%% Set Matrix Equation
NT = nx*ny;
B = zeros(NT,1);
Ael = zeros(5, 1, NT);
ix = zeros(5, 1, NT);
for i = 1:5
    ix(i,1,:) = 1:NT;                              
end
iy = zeros(5, 1, NT);

% Mapping from grid indices to solution vector indices
umap = reshape(1:NT, nx, ny);

for ii = 1:nx
    for jj = 1:ny
            %% Define aE, aW, aN, aS, aP, b
            if ii == nx
                aE = -1/(dx_dim)^2*(Diff(ii,jj)+Diff(ii-1,jj))/2;               % Neumann BC
                aE_old = -1/(dx_dim)^2*(Diff_old(ii,jj)+Diff_old(ii-1,jj))/2;   % Neumann BC
            else
                aE = -1/(dx_dim)^2*(Diff(ii,jj)+Diff(ii+1,jj))/2;
                aE_old = -1/(dx_dim)^2*(Diff_old(ii,jj)+Diff_old(ii+1,jj))/2;
            end

            if ii == 1
                aW = -1/(dx_dim)^2*(Diff(ii,jj)+Diff(ii+1,jj))/2;               % Neumann BC
                aW_old = -1/(dx_dim)^2*(Diff_old(ii,jj)+Diff_old(ii+1,jj))/2;   % Neumann BC 
            else
                aW = -1/(dx_dim)^2*(Diff(ii,jj)+Diff(ii-1,jj))/2;
                aW_old = -1/(dx_dim)^2*(Diff_old(ii,jj)+Diff_old(ii-1,jj))/2;
            end
            
            if jj == ny
                aN = -1/(dy_dim)^2*(Diff(ii,jj)+Diff(ii,jj-1))/2;               % Neumann BC
                aN_old = -1/(dy_dim)^2*(Diff_old(ii,jj)+Diff_old(ii,jj-1))/2;   % Neumann BC
            else
                aN = -1/(dy_dim)^2*(Diff(ii,jj)+Diff(ii,jj+1))/2;
                aN_old = -1/(dy_dim)^2*(Diff_old(ii,jj)+Diff_old(ii,jj+1))/2;
            end

            if jj == 1
                aS = -1/(dy_dim)^2*(Diff(ii,jj)+Diff(ii,jj+1))/2;               % Neumann BC
                aS_old = -1/(dy_dim)^2*(Diff_old(ii,jj)+Diff_old(ii,jj+1))/2;   % Neumann BC
            else
                aS = -1/(dy_dim)^2*(Diff(ii,jj)+Diff(ii,jj-1))/2;
                aS_old = -1/(dy_dim)^2*(Diff_old(ii,jj)+Diff_old(ii,jj-1))/2;
            end   
            
            aP = h0/dt/ts -aE -aW -aN - aS;

            if ii == 1
                old_x =  -(aE_old*h_lube_slider_old(ii+1,jj) + aW_old*h_lube_slider_old(ii+1,jj) + (-aE_old -aW_old)*h_lube_slider_old(ii,jj)); % Neumann BC 
            elseif ii == nx
                old_x =  -(aE_old*h_lube_slider_old(ii-1,jj) + aW_old*h_lube_slider_old(ii-1,jj) + (-aE_old -aW_old)*h_lube_slider_old(ii,jj)); % Neumann BC 
            else
                old_x =  -(aE_old*h_lube_slider_old(ii+1,jj) + aW_old*h_lube_slider_old(ii-1,jj) + (-aE_old -aW_old)*h_lube_slider_old(ii,jj));
            end

            if jj == 1
                old_y =  -(aN_old*h_lube_slider_old(ii,jj+1) + aS_old*h_lube_slider_old(ii,jj+1) + (-aN_old -aS_old)*h_lube_slider_old(ii,jj)); % Neumann BC 
            elseif jj == ny
                old_y =  -(aN_old*h_lube_slider_old(ii,jj-1) + aS_old*h_lube_slider_old(ii,jj-1) + (-aN_old -aS_old)*h_lube_slider_old(ii,jj)); % Neumann BC 
            else
                old_y =  -(aN_old*h_lube_slider_old(ii,jj+1) + aS_old*h_lube_slider_old(ii,jj-1) + (-aN_old -aS_old)*h_lube_slider_old(ii,jj));
            end

            b = f(ii,jj) - old_x - old_y + h0/dt/ts*h_lube_slider_old(ii,jj);

            %% Populate row umap(ii,jj) of coefficient matrix A and RHS vector B(umap(ii,jj))

            if ii == 1 % Neuman BC at x = 0
                 aE = aE + aW;  % Set hW = hE
                 Ael(1, 1, umap(ii,jj)) = aE;   iy(1, 1, umap(ii,jj)) = umap(ii+1, jj);  
                 Ael(2, 1, umap(ii,jj)) = 0;    iy(2, 1, umap(ii,jj)) = 1;
            elseif ii == nx % Neuman BC at x = xend
                 aW = aW + aE;  % Set hE = hW
                 Ael(1, 1, umap(ii,jj)) = 0;    iy(1, 1, umap(ii,jj)) = 2; %Dummy insertion   
                 Ael(2, 1, umap(ii,jj)) = aW;   iy(2, 1, umap(ii,jj)) = umap(ii-1, jj);
            else
                Ael(1, 1, umap(ii,jj)) = aE;    iy(1, 1, umap(ii,jj)) = umap(ii+1, jj); 
                Ael(2, 1, umap(ii,jj)) = aW;    iy(2, 1, umap(ii,jj)) = umap(ii-1, jj);
            end

            if jj == 1 % Neuman BC at y = 0
                 aN = aN + aS;  % Set hS = hN
                 Ael(4, 1, umap(ii,jj)) = aN;   iy(4, 1, umap(ii,jj)) = umap(ii, jj+1);
                 Ael(5, 1, umap(ii,jj)) = 0;    iy(5, 1, umap(ii,jj)) = 5; %Dummy insertion
            elseif jj == ny % Neumann BC at y = yend
                aS = aS + aN;  % Set hN = hS
                Ael(4, 1, umap(ii,jj)) = 0;     iy(4, 1, umap(ii,jj)) = 4; %Dummy insertion
                Ael(5, 1, umap(ii,jj)) = aS;    iy(5, 1, umap(ii,jj)) = umap(ii, jj-1);
            else
                Ael(4, 1, umap(ii,jj)) = aN;    iy(4, 1, umap(ii,jj)) = umap(ii, jj+1);
                Ael(5, 1, umap(ii,jj)) = aS;    iy(5, 1, umap(ii,jj)) = umap(ii, jj-1);
            end

            Ael(3, 1, umap(ii,jj)) = aP;    
            B(umap(ii,jj)) = b;        
    end
end    

iy(3, 1, :) = 1:NT;
A = sparse(ix(:), iy(:), Ael(:), NT, NT);

%% Determine h
h_lube_slider_try = reshape(A\B,nx,ny);

end