tic;
%% Initialize
dbstop if error
AA_init_inputs_3D
Nf = floor(tf/dt/s_interval);
h0_slider = 0.3e-9/h0;
h_lube_disk = ones(nx,ny,Nf);
h_lube_disk_big = ones(nx_d,ny_d,Nf);
h_lube_slider = h0_slider*ones(nx,ny,Nf);
pres_vapor = zeros(nx,ny,Nf);

total_mass_lost = zeros(Nf,1);
timestep_mass_lost = zeros(Nf,1);

n_h_lube_disk_moving_big = ones(nx_d,ny_d);
n_h_lube_disk_moving = ones(nx,ny);
n_h_lube_slider_stationary = h0_slider*ones(nx,ny);
n_h_lube_slider_moving = h0_slider*ones(nx,ny);
n_p_vapor_moving = zeros(nx,ny);
p_vapor_moving = zeros(nx,ny);
n_rho_vapor_moving = zeros(nx,ny);

m_evap_disk = 0;
m_evap_slider = 0;

ifile = 0;
%% Parameters
dx_ND = dx_dim/lrad;
dy_ND = dy_dim/lrad;
tol_lube = 1e-5;
tol_vapor = 1e-3;
%c = 6e-5; 
t_ramp = 2e-9;
x_init = x;
x0_init = 0;

%% Loop
for tstep=1:floor(tf/dt)
	t = dt*tstep;  % current time [s]
    if mod(tstep,s_interval)==0
        disp(['Time = ',num2str(t*ts),' s'])
    end
    
    x0=ux*t*ts; % laser center location [m]
    x0_old=ux*(t-dt)*ts;
    y0=0;
    y0_old=0;

    if mod(tstep,round(dx_dim/(ux*dt*ts))) == 0
        % Create 'new' x vector
        x_temp = zeros(size(x));
        x_temp(1:end-1) = x(2:end);
        x_temp(end) = x_temp(end-1) + dx_ND;
        x = x_temp;
        
        x_temp = zeros(size(x_disk));
        x_temp(1:end-1) = x_disk(2:end);
        x_temp(end) = x_temp(end-1) + dx_ND;
        x_disk = x_temp;
        
        % Initialize
        % Disk
        h_lube_disk_old_moving_big = ones(size(n_h_lube_disk_moving_big));
        h_lube_disk_old_moving_big(1:end-1,:) = n_h_lube_disk_moving_big(2:end,:);
        h_lube_disk_old_moving_big(end,:) = n_h_lube_disk_moving_big(end-1,:); % Neumann BC

        h_lube_disk_old_moving = ones(size(n_h_lube_disk_moving));
        h_lube_disk_old_moving(1:end-1,:) = n_h_lube_disk_moving(2:end,:);
        h_lube_disk_old_moving(end,:) = n_h_lube_disk_moving(end-1,:); % Neumann BC
               
        % Slider
        h_lube_slider_old_stationary = n_h_lube_slider_stationary;
        h_lube_slider_old_moving = h0_slider*ones(size(n_h_lube_slider_moving));
        h_lube_slider_old_moving(1:end-1,:) = n_h_lube_slider_moving(2:end,:);
        h_lube_slider_old_moving(end,:) = n_h_lube_slider_moving(end-1,:); % Neumann BC
        
        % Vapor Pressure
        p_vapor_old_moving = zeros(size(n_p_vapor_moving));
        p_vapor_old_moving(1:end-1,:) = n_p_vapor_moving(2:end,:);
        p_vapor_old_moving(end,:) = n_p_vapor_moving(end-1,:);  % Neumann BC
                
        % Vapor Density
        rho_vapor_old_moving = zeros(size(n_rho_vapor_moving));
        rho_vapor_old_moving(1:end-1,:) = n_rho_vapor_moving(2:end,:);
        rho_vapor_old_moving(end,:) = n_rho_vapor_moving(end-1,:);  % Neumann BC
        
    else
        problem = 1;
        %I choose mod(tstep,round(dx_dim/(ux*dt*ts))) = 0, so give error
        disp('Error: Wrong time step!!!!!!!!')        
    end
        
%%  Shear Stress
    % Set Temperature: Gaussian profile centered at (ux*t,0)
    T_lube_disk_moving_big = T0_disk*ones(nx_d,ny_d);
    T_lube_disk_moving = T0_disk*ones(nx,ny);
    T_lube_slider_stationary = T0_slider*ones(nx,ny);
    T_lube_slider_moving = T0_slider*ones(nx,ny);
    
    FWHM=lrad;    
    sigmax=FWHM/(2*sqrt(2*log(2))); 
    sigmay=FWHM/(2*sqrt(2*log(2))); 
    
    Text_x_lube_disk_moving_big = zeros(size(T_lube_disk_moving_big));
    Text_x_lube_disk_old_moving_big = zeros(size(T_lube_disk_moving_big));
    
    Text_y_lube_disk_moving_big = zeros(size(T_lube_disk_moving_big));
    Text_y_lube_disk_old_moving_big = zeros(size(T_lube_disk_moving_big));
    
    Text_x_lube_slider_stationary = zeros(size(T_lube_slider_stationary));
    Text_y_lube_slider_stationary = zeros(size(T_lube_slider_stationary));
    
    Text_x_lube_slider_old_stationary = zeros(size(T_lube_slider_stationary));
    Text_y_lube_slider_old_stationary = zeros(size(T_lube_slider_stationary));
    
    for i=2:nx-1
        for j = 2:ny-1
           T_lube_disk_moving(i,j)=min(t*ts/t_ramp,1)*dT_disk*exp(-((x(i)*lrad-x0)^2/(2*sigmax^2)+(y(j)*lrad-y0)^2/(2*sigmay^2)))+T_lube_disk_moving(i,1);
           T_lube_slider_moving(i,j)=min(t*ts/t_ramp,1)*dT_slider*exp(-((x(i)*lrad-x0)^2/(2*sigmax^2)+(y(j)*lrad-y0)^2/(2*sigmay^2)))+T_lube_slider_moving(i,1);
           T_lube_slider_stationary(i,j)=min(t*ts/t_ramp,1)*dT_slider*exp(-((x_init(i)*lrad-x0_init)^2/(2*sigmax^2)+(y(j)*lrad-y0)^2/(2*sigmay^2)))+T_lube_slider_stationary(i,1);
           
           Text_x_lube_slider_stationary(i,j) = min(t*ts/t_ramp,1)*c/sigmax^2*(x_init(i)*lrad-x0_init)*dT_slider*exp(-((x_init(i)*lrad-x0_init)^2/(2*sigmax^2)+(y(j)*lrad-y0)^2/(2*sigmay^2)));
           Text_y_lube_slider_stationary(i,j) = min(t*ts/t_ramp,1)*c/sigmay^2*(y(j)*lrad-y0)*dT_slider*exp(-((x_init(i)*lrad-x0_init)^2/(2*sigmax^2)+(y(j)*lrad-y0)^2/(2*sigmay^2)));

           if tstep>1
                Text_x_lube_slider_old_stationary(i,j) = min((t-dt)*ts/t_ramp,1)*c/sigmax^2*(x_init(i)*lrad-x0_init)*dT_slider*exp(-((x_init(i)*lrad-x0_init)^2/(2*sigmax^2)+(y(j)*lrad-y0)^2/(2*sigmay^2)));
                Text_y_lube_slider_old_stationary(i,j) = min((t-dt)*ts/t_ramp,1)*c/sigmay^2*(y(j)*lrad-y0)*dT_slider*exp(-((x_init(i)*lrad-x0_init)^2/(2*sigmax^2)+(y(j)*lrad-y0)^2/(2*sigmay^2)));
           end
        end
    end
    
    for i=2:nx_d-1
        for j = 2:ny_d-1
           T_lube_disk_moving_big(i,j)=min(t*ts/t_ramp,1)*dT_disk*exp(-((x_disk(i)*lrad-x0)^2/(2*sigmax^2)+(y_disk(j)*lrad-y0)^2/(2*sigmay^2)))+T_lube_disk_moving_big(i,1);

           Text_x_lube_disk_moving_big(i,j) = min(t*ts/t_ramp,1)*c/sigmax^2*(x_disk(i)*lrad-x0)*dT_disk*exp(-((x_disk(i)*lrad-x0)^2/(2*sigmax^2)+(y_disk(j)*lrad-y0)^2/(2*sigmay^2)));
           Text_y_lube_disk_moving_big(i,j) = min(t*ts/t_ramp,1)*c/sigmay^2*(y_disk(j)*lrad-y0)*dT_disk*exp(-((x_disk(i)*lrad-x0)^2/(2*sigmax^2)+(y_disk(j)*lrad-y0)^2/(2*sigmay^2)));
           
           if tstep>1
                Text_x_lube_disk_old_moving_big(i,j) = min((t-dt)*ts/t_ramp,1)*c*dT_disk/sigmax^2*(x_disk(i)*lrad-x0_old)*exp(-((x_disk(i)*lrad-x0_old)^2/(2*sigmax^2)+(y_disk(j)*lrad-y0_old)^2/(2*sigmay^2)));
                Text_y_lube_disk_old_moving_big(i,j) = min((t-dt)*ts/t_ramp,1)*c*dT_disk/sigmay^2*(y_disk(j)*lrad-y0_old)*exp(-((x_disk(i)*lrad-x0_old)^2/(2*sigmax^2)+(y_disk(j)*lrad-y0_old)^2/(2*sigmay^2)));
           end
        end
    end
    
    T_vapor_moving = 0.5*(T_lube_disk_moving + T_lube_slider_moving);
    
    Text_x_lube_disk_moving_big=Text_x_lube_disk_moving_big/(c*dT_disk/lrad);
    Text_x_lube_disk_old_moving_big=Text_x_lube_disk_old_moving_big/(c*dT_disk/lrad);
    Text_y_lube_disk_moving_big=Text_y_lube_disk_moving_big/(c*dT_disk/lrad);
    Text_y_lube_disk_old_moving_big=Text_y_lube_disk_old_moving_big/(c*dT_disk/lrad);
    
    Text_x_lube_slider_stationary=Text_x_lube_slider_stationary/(c*dT_slider/lrad);
    Text_y_lube_slider_stationary=Text_y_lube_slider_stationary/(c*dT_slider/lrad);
    Text_x_lube_slider_old_stationary=Text_x_lube_slider_old_stationary/(c*dT_slider/lrad);
    Text_y_lube_slider_old_stationary=Text_y_lube_slider_old_stationary/(c*dT_slider/lrad);
    
    dText_x_dt_disk_moving_big = (Text_x_lube_disk_moving_big - Text_x_lube_disk_old_moving_big)/dt;
    dText_y_dt_disk_moving_big = (Text_y_lube_disk_moving_big - Text_y_lube_disk_old_moving_big)/dt;
    
    dText_x_dt_slider_stationary = (Text_x_lube_slider_stationary - Text_x_lube_slider_old_stationary)/dt;
    dText_y_dt_slider_stationary = (Text_y_lube_slider_stationary - Text_y_lube_slider_old_stationary)/dt;
    
    %% Loop
    go = true;
    count = 0;
    while go
        count = count + 1;
        %% Solve for Disk Lube: Get n_h_lube_disk_moving using p_vapor_moving 
        if tstep*dt*ts < 20e-9
            go_disk = true;
            count_disk = 0;
            % Convert from "p_vapor_moving" to "p_vapor_moving_big"
            diff_x_d = (nx_d - nx)/2;
            diff_y_d = (ny_d - ny)/2;
            p_vapor_moving_big = zeros(nx_d,ny_d);
            p_vapor_moving_big(1+diff_x_d:nx_d-diff_x_d,1+diff_y_d:ny_d-diff_y_d) = p_vapor_moving;
            
            while go_disk
                count_disk = count_disk + 1;

                [h_lube_disk_try,m_dot_disk] = get_h_linear_disk(p_vapor_moving_big,Text_x_lube_disk_moving_big,Text_y_lube_disk_moving_big,dText_x_dt_disk_moving_big,dText_y_dt_disk_moving_big,T_lube_disk_moving_big,G0_lube,n_h_lube_disk_moving_big,h_lube_disk_old_moving_big,dt,ts,h0,c,dT_disk,lrad,dx_ND,dy_ND,dx_dim,dy_dim,M,rho,lube_name,nx_d,ny_d);
                
                error_lube_disk = max(max(abs(n_h_lube_disk_moving_big - h_lube_disk_try)));
    %            disp(['Iteration number: ',num2str(count_disk),' Error for h Disk: ',num2str(error_lube_disk)])

                if error_lube_disk < tol_lube %&& count_disk > 1
                    go_disk = false;
                    n_h_lube_disk_moving_big = h_lube_disk_try;
    %                disp(['Iteration number: ',num2str(count_disk),' Error for h Disk: ',num2str(error_lube_disk)])
                else
                    n_h_lube_disk_moving_big = h_lube_disk_try;
                    if error_lube_disk > 1 || count_disk > 100
                        problem = 1;
                        disp(['Iteration number: ',num2str(count_disk),' Error for h Disk: ',num2str(error_lube_disk)])
                        disp('Error Disk')
                        go_disk = false;
                    end
                end    
            end
            % Convert from "n_h_lube_disk_moving_big" to "n_h_lube_disk_moving"
            diff_x_d = (nx_d - nx)/2;
            diff_y_d = (ny_d - ny)/2;
            n_h_lube_disk_moving = n_h_lube_disk_moving_big(1+diff_x_d:nx_d-diff_x_d,1+diff_y_d:ny_d-diff_y_d);
        else
            % Do nothing, "n_h_lube_disk_moving" stays the same
            error_lube_disk = 0;
            count_disk = 0;
            go_disk = false;
        end

        %% Solve for Slider: Get n_h_lube_slider_stationary using p_vapor_stationary 
        go_slider = true;
        count_slider = 0;
        
        % Convert from p_vapor_moving (moving frame fixed to disk) to p_vapor_stationary (frame fixed to slider)
        diff = (x0 - x(ceil(nx/2))*lrad)/dx_dim;
        if abs(diff)<1e-7
            diff = 0; 
        end
        p_vapor_stationary = zeros(size(p_vapor_moving));
        p_vapor_stationary(1:end-1,:) = (1-diff)*p_vapor_moving(1:end-1,:) + (diff)*p_vapor_moving(2:end,:);        
        p_vapor_stationary(end,:) = p_vapor_moving(end,:);

%         % Diff is 0 since I choose mod(tstep,round(dx_dim/(ux*dt*ts))) = 0
        
        while go_slider
            count_slider = count_slider + 1;

            [h_lube_slider_try,m_dot_slider] = get_h_linear_slider(p_vapor_stationary,Text_x_lube_slider_stationary,Text_y_lube_slider_stationary,dText_x_dt_slider_stationary,dText_y_dt_slider_stationary,T_lube_slider_stationary,G0_lube,n_h_lube_slider_stationary,h_lube_slider_old_stationary,dt,ts,h0,c,dT_slider,lrad,dx_ND,dy_ND,dx_dim,dy_dim,M,rho,lube_name,nx,ny,ux);


            error_lube_slider = max(max(abs(n_h_lube_slider_stationary - h_lube_slider_try)));
%            disp(['Iteration number: ',num2str(count_slider),' Error for h Slider: ',num2str(error_lube_slider)])

            if error_lube_slider < tol_lube %&& count_slider > 1
                go_slider = false;
                n_h_lube_slider_stationary = h_lube_slider_try;
%                disp(['Iteration number: ',num2str(count_slider),' Error for h Slider: ',num2str(error_lube_slider)])
            else
                n_h_lube_slider_stationary = h_lube_slider_try;
                if error_lube_slider > 1 || count_slider > 100
                    problem = 1;
                    disp(['Iteration number: ',num2str(count_slider),' Error for h Slider: ',num2str(error_lube_slider)])
                    disp('Error slider')
                    go_slider = false;
                end
            end    
        end
        
        % Convert from n_h_lube_slider_stationary (stationary frame fixed to slider) to n_h_lube_slider_moving (moving frame fixed to disk)
        diff = (x0 - x(ceil(nx/2))*lrad)/dx_dim;
        if abs(diff)<1e-7
            diff = 0; 
        end
        n_h_lube_slider_moving = h0_slider*ones(size(n_h_lube_slider_stationary));
        n_h_lube_slider_moving(2:end,:) = (diff)*n_h_lube_slider_stationary(1:end-1,:) + (1-diff)*n_h_lube_slider_stationary(2:end,:);        
        n_h_lube_slider_moving(1,:) = n_h_lube_slider_stationary(1,:);
        
%         % Diff is 0 since I choose mod(tstep,round(dx_dim/(ux*dt*ts))) = 0
        
        %% Solve for Vapor: Get n_p_vapor_moving using n_h_lube_slider_moving & n_h_lube_disk_moving
        h_a_moving = fh - n_h_lube_disk_moving - n_h_lube_slider_moving; 
        h_a_old_moving = fh - h_lube_disk_old_moving - h_lube_slider_old_moving; 
        D = getD(T_vapor_moving,2.2e6*ones(size(T_vapor_moving)));
        [n_p_vapor_moving, n_rho_vapor_moving] = get_p_vapor(T_lube_slider_moving,T_lube_disk_moving,T_vapor_moving,n_h_lube_disk_moving,n_h_lube_slider_moving,h_a_moving,h_a_old_moving,rho_vapor_old_moving,D,M,dt,ts,h0,rho,nx,ny,dx_dim,dy_dim,lube_name);
        
        %% Check for Convergence
        error_vapor = max(max(abs(n_p_vapor_moving - p_vapor_moving)));
        if error_vapor < tol_vapor 
            go = false;
            p_vapor_moving = n_p_vapor_moving;
            if mod(tstep,s_interval)==0
                disp(['Iteration number: ',num2str(count_disk),' Error for h Disk: ',num2str(error_lube_disk)])
                disp(['Iteration number: ',num2str(count_slider),' Error for h Slider: ',num2str(error_lube_slider)])
                disp(['Iteration number: ',num2str(count),' Error for vapor: ',num2str(error_vapor)])
            end
        else
            p_vapor_moving = n_p_vapor_moving;
            if count > 100
                problem = 1;
                disp(['Iteration number: ',num2str(count),' Error for vapor: ',num2str(error_vapor)])
                disp('Error Vapor')
                go = false;
            end
        end 
    end
    
    %% Mass Conservation
    m_evap_disk = m_evap_disk + int(m_dot_disk)*(dx_dim*dy_dim)*(dt*ts);
    m_evap_slider = m_evap_slider + int(m_dot_slider)*(dx_dim*dy_dim)*(dt*ts);
    
    aa = h_a_moving.*n_rho_vapor_moving*h0;
    aa_old = h_a_old_moving.*rho_vapor_old_moving*h0;
    m_vapor = int(aa)*(dx_dim*dy_dim);

    %% Save
    if mod(tstep,s_interval)==0
        ifile = ifile+1;
        h_lube_disk(:,:,ifile) = n_h_lube_disk_moving;
        h_lube_disk_big(:,:,ifile) = n_h_lube_disk_moving_big;
        h_lube_slider(:,:,ifile) = n_h_lube_slider_moving;
        pres_vapor(:,:,ifile) = n_p_vapor_moving;
        
        total_mass_lost(ifile) = m_vapor - m_evap_disk - m_evap_slider;
        timestep_mass_lost(ifile) = int(aa-aa_old)*(dx_dim*dy_dim) - int(m_dot_disk)*(dx_dim*dy_dim)*(dt*ts) - int(m_dot_slider)*(dx_dim*dy_dim)*(dt*ts);
    end
end

toc;

    
