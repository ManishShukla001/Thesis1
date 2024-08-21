
clear; clc;

%-------------------------------------------------------------------------%
% Read dataaxis([-0.96 0.48 min(y) max(y) min(z) max(z)])

%for g = 20 %3000 %for loop for time steps
g=399;   
disp(g)
   % for p = 1 % for loop for planes
        
        p=1;
        disp(p)


        
    for kk=1:2
        %g=700+2*kk;
        start = tic;
        % data_files_path = 'data_on_planes/XY_plane/plane1';
        % data_files_prefix = 'vel_XYplane1_';
    
        data_files_path = strcat('extracted_data');  %,num2str(p));
       % data_files_prefix = strcat('selected_data_');%,num2str(p));
        data_files_prefix = strcat('ts')
        saving_file_path = strcat('Result');%,num2str(p));
   
    
        file_numbers =   sort(200:1:g,'descend');
        % file_numbers =   sort(1:1:50,'ascend');
        
        data_files_postfix = '.csv';
        
        %data_files_postfix = '.xlsx';
        number_of_files = length(file_numbers);
        dt=10;
        
        % T=number_of_files*dt;        % Total integration /advection time 
        T = number_of_files*dt;
        
        % Column numbers of Data in excel files.
        x_coord_colm = 1;
        y_coord_colm = 2;
        x_vel_colm = 3;
        y_vel_colm = 4;
        
        %% Part 1 - Initialize grid of particles through vector field
        
        fileindex = 1; 
        filename = strcat(data_files_path, '/', data_files_prefix, num2str(file_numbers(fileindex)), data_files_postfix);
        table = readmatrix(filename);
        x_coord = table(:,x_coord_colm);
        y_coord = table(:,y_coord_colm);
        x_vel = table(:,x_vel_colm);
        y_vel = table(:,y_vel_colm);
        
        x_min = min(x_coord);
        y_min = min(y_coord);
        x_max = max(x_coord);
        y_max = max(y_coord);
        
        %% Intial grid of particles
        % yin = yIC;
        resolution1 = 300; 
        resolution2 = 300;
        
        x_interp = linspace(x_min, x_max, resolution1);
        y_interp = linspace(y_min, y_max, resolution2);
        [XXb, YYb] = meshgrid(x_interp, y_interp);
        XX1 = reshape(XXb, resolution1*resolution2, 1);
        YY1 = reshape(YYb, resolution1*resolution2, 1);
    
        save(strcat(saving_file_path,'/XXb','_plane',num2str(p)),'XXb')
        save(strcat(saving_file_path,'/YYb','_plane',num2str(p)),'YYb')
        
        %% Compute_trajectories    % 
        x_in = XX1;
        y_in = YY1;
        
        for i = 1:1:number_of_files  % i=1:T/dt. % One file is already read/Initialized.
        
            for_tic = tic;
            filename = strcat(data_files_path, '/', data_files_prefix, num2str(file_numbers(fileindex)), data_files_postfix);
            disp(filename)
            table = readmatrix(filename);
            fileindex = fileindex + 1;
            x_coord = table(:,x_coord_colm);
            y_coord = table(:,y_coord_colm);
            x_vel = table(:,x_vel_colm);
            y_vel = table(:,y_vel_colm);
        
            x_vel_Int = scatteredInterpolant(x_coord,y_coord,x_vel);
            y_vel_Int = scatteredInterpolant(x_coord,y_coord,y_vel);
        
            % Flow Map interpolate new positions
        
            % RK-4 integrator
        
            k1=dt*x_vel_Int(x_in, y_in);
            l1=dt*y_vel_Int(x_in, y_in);
        
            k2=dt*x_vel_Int(x_in + k1/2, y_in + l1/2);
            l2=dt*y_vel_Int(x_in + k2/2, y_in + l1/2);
        
            k3=dt*x_vel_Int(x_in + k2/2, y_in + l2/2);
            l3=dt*y_vel_Int(x_in + k2/2, y_in + l2/2);
        
            k4=dt*x_vel_Int(x_in + k3, y_in + l3);
            l4=dt*y_vel_Int(x_in + k3, y_in + l3);
        
            % Flow Map interpolate new positions
        
            % for backward LCS
            x_out = x_in-1/6*(k1+2*k2+2*k3+k4);
            y_out = y_in-1/6*(l1 + 2*l2 + 2*l3 + l4);
        
            % % For forward LCS calculation
%           x_out = x_in+1/6*(k1+2*k2+2*k3+k4);
%           y_out = y_in+1/6*(l1 + 2*l2 + 2*l3 + l4);
        
            x_in = x_out;
            y_in = y_out;
             
            clear x_coord y_coord x_vel y_vel x_vel_Int y_vel_Int 
            for_loop_time = toc(for_tic);
            disp(for_loop_time)
        
        end
        
        %% Part 3 -  Compute the finite-time Lyapunov exponent (sigma)
        xT=reshape(x_out, resolution2,resolution1);
        yT=reshape(y_out, resolution2,resolution1);       
        
        dx=(x_max - x_min)/(resolution1-1);
        dy=(y_max - y_min)/(resolution2-1);
                
        % Finite difference to compute the gradient
        [dxTdx0,dxTdy0] = gradient(xT,dx,dy);
        [dyTdx0,dyTdy0] = gradient(yT,dx,dy);
        
        % Calculation of  Finite time Lyapnuov exponent (FTLE): Forward and Backward FTLE
        sigmab = zeros(resolution2,resolution1);
        for i=1: resolution2
            for j=1:resolution1
                
                D(1,1) = dxTdx0(i,j);
                D(1,2) = dxTdy0(i,j);
                D(2,1) = dyTdx0(i,j);
                D(2,2) = dyTdy0(i,j);
                
                lamda=max(eig(D'*D)); % cOMPUTE MAXIMUM EIGEN VALUES OF cAUCHY gREEN TENSOR
                sigmab(i,j)=1/abs(2*T)*log(lamda);
                
                %sigmab(i,j) = (1/T)*sqrt(max(eig(D'*D)))
            end
        end
        % % Defining a velocity field which will be zero inside the cylinders
        a = 0.46/2; %Radius of the cylinders
        %
        % Location coordinate (CX and CY) of the cylinders (C)
        CX = 0; % X-coordinate of cylinder center
        CY = 0; % Y-coordinate of cylinder center
        % Set values of X and Y inside the cylinder to zero
        
    %     [I1, J1]=find( sqrt((XXb-CX).^2+(YYb-CY).^2) < a );
    %     
    %     T1 =isempty(I1);
    %     sigmab(I1)=NaN;
        
    %     writematrix(sigmab,sprintf('data_on_planes/XY_plane/plane1/LCS_results/sigmab_XYplane1_%d%',g))
        writematrix(sigmab,strcat(saving_file_path,'/NSigma',num2str(p),sprintf('_%d%',g)));
        total_time = toc(start);
         
    end
        disp(total_time)
    beep
  %  end
%end
