function [n_p_vapor, n_rho_vapor] = get_p_vapor(T_lube_slider,T_lube_disk,T_vapor,n_h_lube_disk,n_h_lube_slider,h_a,h_a_old,rho_vapor_old,D,M,dt,ts,h0,rho,nx,ny,dx_dim,dy_dim,lube_name)

A_ham_disk = 5.3e-21;
A_ham_slider = 5.3e-21;

% A_ham = 3e-20/6/pi;

% temp dependent hammaker
% A_ll = 3.78e-20*(1-1.72e-3*(T_lube_disk-273.15));
% A_ss = 1.16e-19;
% A_ham_disk = (sqrt(A_ll*A_ss)-A_ll)/6/pi;
% 
% A_ll = 3.78e-20*(1-1.72e-3*(T_lube_slider-273.15));
% A_ss = 1.16e-19;
% A_ham_slider = (sqrt(A_ll*A_ss)-A_ll)/6/pi;

d0 = 0.3e-9;
p_disj_disk = A_ham_disk./(n_h_lube_disk*h0 + d0).^3;
p_disj_slider = A_ham_slider./(n_h_lube_slider*h0 + d0).^3;

p_film_disk = getevap_vapor(T_lube_disk,p_disj_disk,M,rho,lube_name);
p_film_slider = getevap_vapor(T_lube_slider,p_disj_slider,M,rho,lube_name);

R = 8.314;
inv = h_a.*h0/dt/ts + sqrt(M./(2*pi*R.*T_lube_slider)).*(R.*T_vapor./M) + sqrt(M./(2*pi*R.*T_lube_disk)).*(R.*T_vapor./M);
rhs = ( rho_vapor_old.*h_a_old.*h0/(dt*ts) ...
                                + sqrt(M./(2*pi*R.*T_lube_disk)).*p_film_disk ...
                                + sqrt(M./(2*pi*R.*T_lube_slider)).*p_film_slider );
%%
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
            aE = -1/(dx_dim)^2*(D(ii,jj)*h_a(ii,jj)+D(ii-1,jj)*h_a(ii-1,jj))/2*h0; % Neumann BC
        else
            aE = -1/(dx_dim)^2*(D(ii,jj)*h_a(ii,jj)+D(ii+1,jj)*h_a(ii+1,jj))/2*h0;
        end

        if ii == 1
            aW = -1/(dx_dim)^2*(D(ii,jj)*h_a(ii,jj)+D(ii+1,jj)*h_a(ii+1,jj))/2*h0; % Neumann BC    
        else
            aW = -1/(dx_dim)^2*(D(ii,jj)*h_a(ii,jj)+D(ii-1,jj)*h_a(ii-1,jj))/2*h0;    
        end

        if jj == ny
            aN = -1/(dy_dim)^2*(D(ii,jj)*h_a(ii,jj)+D(ii,jj-1)*h_a(ii,jj-1))/2*h0; % Neumann BC
        else
            aN = -1/(dy_dim)^2*(D(ii,jj)*h_a(ii,jj)+D(ii,jj+1)*h_a(ii,jj+1))/2*h0;
        end

        if jj == 1
            aS = -1/(dy_dim)^2*(D(ii,jj)*h_a(ii,jj)+D(ii,jj+1)*h_a(ii,jj+1))/2*h0; % Neumann BC
        else
            aS = -1/(dy_dim)^2*(D(ii,jj)*h_a(ii,jj)+D(ii,jj-1)*h_a(ii,jj-1))/2*h0;
        end   
        aP = inv(ii,jj) -aE - aW - aS - aN;
        b = rhs(ii,jj);

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
n_rho_vapor = reshape(A\B,nx,ny);

n_p_vapor = n_rho_vapor.*R.*T_vapor./M;

end