%%
% Simulations for a geometrically exact beam clamped at one end 
% and with velocity feedback control applied at the other end
% We work with the Intrinsic formulation of the GEB model (Hodges 2003)
% We use P2 elements for the spatial discretization

%%
clear all
close all
clc
fs = 14; % font size
size_line = 2; % line size for plots;


%% To be chosen
linearized = false;     % true: we take the IGEB system without g(y)
                        % false: we take the real (semilinear) IGEB system
toyIniData = false;     % true: we take toy initial data y^0 (bump function)
                        % false: we take meaningful initial data
BC = 1;   % 0: we take transparent boundary conditions (BC)
          % 1: we take something close to transparent BC
          % 2: we take somthing far from transparent BC   (does not work)
plot_y = true;     % true: we plot the solution y
                   % false: we don't plot the solution y
plot_y0 = true;  % true: we plot the initial data y0
                 % false: we don't plot the initial data y0
plot_centerline = true;   % true: we plot the centerline of the beam
                          % false: we don't plot the centerline  
centerline_t_scheme = 0;     % 0: mid point rule
                             % 1: explicit Euler      (does not work properly)
                             % 2: ode45 (RK)
                             % 3: implicit Euler
                             % 4: zupan paper scheme  (does not work properly)
centerline_x_scheme = 0; % 0: mid point rule
                         % 1: explicit Euler      (probably does not work properly either)     
zeroOrderBC = false; % true: we choose the initial velocities so that the 
                    % the zero-order compatibility conditions of the IGEB
                    % model are fulfilled
                    % false: the initial velocities are set to zero
                         
%%% (meaninful initial data) beam of length 1
ell = 1;      % length of the space interval
Ne = 15*ell;      % number of elements
Tfactor = ell*5;  % the final time will be equal to :
                %%% then T = 0.2 * factor and Nt = 100 * Tfactor
           
% display options:
if ell == 3
    view_centerline = 2;
    p_plot = 15;   % plot centerline every p_plot time instances
else
    view_centerline = 1;
    p_plot = 3;   % plot centerline every p_plot time instances
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         SPACE DISCRETIZATION                            %
% nodes:      1   2   3   4   5     2e-1 2e  2e+1               Nx-1  Nx  %
%             |---o---|---o---|  ...  |---o---|---o---|---o---|---o---|   %
% elements:       1       2               e              Ne-1     Ne      %
%                                                                         %
%                        TIME DISCRETIZATION                              %
% time instances:     1       2                     Nt-1     Nt           %
%                     |-------|-------|  ...  |-------|-------|           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% curvature before deformation
kap = [0; 0; 0];  % curvature before deformation

%% The beam parameters: toy parameters
if toyIniData == true
    ell = 1;      % length of the space interval
    p_a = 1;      % area of cross section
    p_rho = 1;    % density
    p_I2 = 1;     % area moment of inertia
    p_I3 = 1;     % area moment of inertia
    p_k1 = 1;     % factor correcting the moment of area
    p_k2 = 1;     % shear correction factor
    p_k3 = 1;     % shear correctionf actor
    p_E = 1;      % young modulus
    p_G = 1;      % shear modulus
    
    % Mass, flexibility, feedback and initial strain matrices
    p_J = diag([(p_I2 + p_I3)*p_k1, p_I2, p_I3]);   % inertia matrix
    p_S1 = p_a*diag([p_E, p_k2*p_G, p_k3*p_G]);
    p_S2 = p_J*diag([p_G, p_E, p_E]);
    flexMat = inv(blkdiag(p_S1, p_S2)); 
    massMat = p_rho*blkdiag(p_a*eye(3), p_J); 
end
% ---------------------------------------- %
% parameters from the paper:
% Consistent structural linearisation in flexible-body dynamics with large rigid-body motion
% Henrik Hesse, Rafael Palacios, 2012
EA = 10^4; GAs = 10^4;
EI = 500; GJ = 500;
rhoA = 1; rhoJ = diag([20, 10, 10]);
massMat = blkdiag(rhoA*eye(3), rhoJ);            % the MASS matrix
flexMat = inv(diag([EA, GAs, GAs, GJ, EI, EI])); % the FLEXIBILITY matrix
%massMat = diag([1, 1, 1, 20, 10, 10]);
%flexMat = inv(diag([10^4, 10^4, 10^4, 500, 500, 500]));
% ---------------------------------------- %

boldE = zeros(6, 6); boldE(5, 3) = -1; boldE(6, 2) = 1;   % initial strain matrix

%%
%%% THE PDE SYSTEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We consider the IGEB model.                                            %
% The beam is clamped at x=0.                                            %
% A velocity feedback is applied at x=\ell.                              %
%                                                                        %   
% The unknown is:   y = (v^T, z^T)^T    where:                           %
% v = linear and angular velocities (expressed in body-attached basis),  %
% z = forces and moments (expressed in body-attached basis).             %
%                                                                        %
% The system reads:                                                      %
% J(x) y_t + A y_x + B(x)y = g(x,y)        in (0, \ell) x (0, T)         %
% v(0, t) = 0                              t in (0, T)                   %
% z(\ell, t) = - kappa*v(\ell, t)          t in (0, T)                   %
% y(x,0) = y^0(x)                          x in (0, \ell).               %
%                                                                        %
% M,C are the mass and flexibility matrices.                             %
% kappa is a positive definite symmetric matrix.                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% The feedback matrix
kappa = eye(6)*flexMat^(-1/2)*massMat^(1/2); % transparent BC
diag_kappa = diag(kappa);
mu1 = sqrt(max(diag_kappa(1:3)*min(diag_kappa(1:3))));
mu2 = sqrt(max(diag_kappa(4:6)*min(diag_kappa(4:6))));
new_kappa_diag = [mu1, mu1, mu1, mu2, mu2, mu2];
if BC == 1          % close to transparent BC
    kappa = diag(new_kappa_diag);
elseif BC == 2      % far from transparent BC
    kappa = kappa +  diag(rand(1, 6))*diag([0.001, 0.001, 0.001, 0.1, 0.1, 0.1])*5;
end
%kappa = mean(diag(kappa))*eye(6);  % test
%kappa = eye(6)                     % test
%kappa = zeros(6);                  % test

%% Coeffficients of the PDE system
J = blkdiag(massMat, flexMat);
B = zeros(12, 12); B(1:6, 7:12) = - boldE; B(7:12, 1:6) = transpose(boldE);
A = zeros(12, 12); A(1:6, 7:12) = -eye(6); A(7:12, 1:6) = -eye(6);
Ni = 12;      % number of PDEs

%% Some dependent variables
Nx = 2*Ne + 1;               % number of nodes
Nt = 100*Tfactor;            % number of time instances
T = 0.2*Tfactor;             % end of the time interval
he = ell/Ne;                 % length of one element
x = linspace(0,ell,Nx);      % spatial grid with the node positions
Ntot = Ni*Nx;                % number of unknowns without BC
Nf = Ntot-6;                 % degree of freedom
ht = T/(Nt-1);               % time step
t = linspace(0, T, Nt);      % time instances

if view_centerline == 1
    %orient_cent = [90, -90; -37.5, 30];
    orient_cent = [96, 10; -37.5, 30];
    locLegend = ["southeast", "northwest"];
else
    orient_cent = [0, 0; 10, 40];
    locLegend = ["northwest", "northwest"];
end

if linearized == true
    linNonlin = 'lin';
else
    linNonlin = 'nonlin';
end
if BC == 0
    type_K = 'transp';
elseif BC == 1
    type_K = 'closeTransp';
else
    type_K = 'NoTransp';
end
if zeroOrderBC == true
    zero_ord = 'zero-ord';
else
    zero_ord = 'no-zero-ord';
end
file_name_y = ['fig/SOLY_len', num2str(ell), '_', linNonlin, '_', type_K, '_', zero_ord, '.pdf'];
name_ini = ['fig/Y0_len', num2str(ell), '_', type_K, '_', zero_ord, '.pdf'];

%% Element matrices
syms XI;
Ntild =  [(1-XI)*(1-2*XI), 4*XI*(1-XI), XI*(2*XI - 1)];
Dxi_Ntild = diff(Ntild, XI);

Me_syms = int((Ntild')*Ntild, 0, 1);
Ke_syms = int((Dxi_Ntild')*Ntild, 0, 1);
P1e_syms = int(Ntild(1)*(Ntild')*Ntild, 0, 1);
P2e_syms = int(Ntild(2)*(Ntild')*Ntild, 0, 1);
P3e_syms = int(Ntild(3)*(Ntild')*Ntild, 0, 1);

Me = double(Me_syms);
Ke = double(Ke_syms);
P1e = double(P1e_syms);
P2e = double(P2e_syms);
P3e = double(P3e_syms);

%% Node numbers
NNB = reshape(1:Ntot, Ni, Nx);        % i.e. NNB(i, k) = 2*(k-1)+i
NNBc = NNB-6; NNBc = subplus(NNBc);   % with the Dirichlet BC

%% Assemble the mass and stiffness matrices
disp('Assembling the mass and stiffness matrices..')
tic

M = sparse(Ntot,Ntot);    % initialize zero mass matrix 
K = sparse(Ntot,Ntot);    % initialize zero stiffness matrix

for ii = 1:Ni
    for jj = 1:Ni
        for ee = 1:Ne    
            idxR = [NNB(ii, 2*ee-1), NNB(ii, 2*ee), NNB(ii, 2*ee+1)];
            idxC = [NNB(jj, 2*ee-1), NNB(jj, 2*ee), NNB(jj, 2*ee+1)];
            M(idxR, idxC) = M(idxR, idxC) + J(ii, jj)*he*Me;
            K(idxR, idxC) = K(idxR, idxC) + B(ii, jj)*he*Me - A(ii, jj)*Ke;
        end
    end
end

% boundary terms
for ii = 1:6
    for jj = 1:6
        K(NNB(ii, Nx), NNB(jj, Nx)) = K(NNB(ii, Nx), NNB(jj, Nx)) + kappa(ii, jj);
    end
    K(NNB(ii+6, Nx), NNB(ii, Nx)) = K(NNB(ii+6, Nx), NNB(ii, Nx)) -1;
end

toc

%% Assemble the matrices that define the nonlinearity
disp('Assembling the matrices defining the map Q(y)..')
tic

P1 = cell(Ni, Ne);
P2 = cell(Ni, Ne);
P3 = cell(Ni, Ne);
Pdagger1 = cell(Ni, Ne);
Pdagger2 = cell(Ni, Ne);
Pdagger3 = cell(Ni, Ne);
for pp = 1:Ni
    for ee = 1:Ne
        P1pe = sparse(Ntot, Ntot);   % Initialize zero matrices
        P2pe = sparse(Ntot, Ntot);
        P3pe = sparse(Ntot, Ntot);
        Pdagger1pe = sparse(Ntot, Ntot);
        Pdagger2pe = sparse(Ntot, Ntot);
        Pdagger3pe = sparse(Ntot, Ntot);
        Gp = G(pp, massMat, flexMat);
        Gdaggerp = Gdagger(pp, massMat, flexMat);
        for ii=1:Ni
            for jj = 1:Ni
                idxR = [NNB(ii, 2*ee-1), NNB(ii, 2*ee), NNB(ii, 2*ee+1)];
                idxC = [NNB(jj, 2*ee-1), NNB(jj, 2*ee), NNB(jj, 2*ee+1)];
                P1pe(idxR, idxC) = P1pe(idxR, idxC) + Gp(ii, jj)*he*P1e;
                P2pe(idxR, idxC) = P2pe(idxR, idxC) + Gp(ii, jj)*he*P2e;
                P3pe(idxR, idxC) = P3pe(idxR, idxC) + Gp(ii, jj)*he*P3e;
                Pdagger1pe(idxR, idxC) = Pdagger1pe(idxR, idxC) + Gdaggerp(ii, jj)*he*P1e;
                Pdagger2pe(idxR, idxC) = Pdagger2pe(idxR, idxC) + Gdaggerp(ii, jj)*he*P2e;
                Pdagger3pe(idxR, idxC) = Pdagger3pe(idxR, idxC) + Gdaggerp(ii, jj)*he*P3e;
            end
        end
        P1(pp, ee) = {P1pe};
        P2(pp, ee) = {P2pe};
        P3(pp, ee) = {P3pe};  
        Pdagger1(pp, ee) = {Pdagger1pe};
        Pdagger2(pp, ee) = {Pdagger2pe};
        Pdagger3(pp, ee) = {Pdagger3pe}; 
    end
end

toc

%% Apply Dirichlet boundary conditions
dof = 7:Ntot;                         % degrees of fredom
M = M(dof, dof); K = K(dof, dof); 
for pp = 1:Ni
    for ee = 1:Ne
        temp1 = cell2mat(P1(pp, ee)); P1(pp, ee) = {temp1(dof, dof)};
        temp2 = cell2mat(P2(pp, ee)); P2(pp, ee) = {temp2(dof, dof)};
        temp3 = cell2mat(P3(pp, ee)); P3(pp, ee) = {temp3(dof, dof)};
        temp4 = cell2mat(Pdagger1(pp, ee)); Pdagger1(pp, ee) = {temp4(dof, dof)};
        temp5 = cell2mat(Pdagger2(pp, ee)); Pdagger2(pp, ee) = {temp5(dof, dof)};
        temp6 = cell2mat(Pdagger3(pp, ee)); Pdagger3(pp, ee) = {temp6(dof, dof)};
    end
end

%% Initial data
disp('Building the initial data..')
tic

if toyIniData == true
    %%% meaningless initial data %%%
    x0 = 0.5*ell;   % center
    a = 0.2*ell;    % width
    c0 = 0.1;       % magnitude
    f = @(x) c0 * exp( -1/a^2*((x-x0).^2) ) ; % bump function
    Y0 = zeros(Ntot, 1);
    for ii = 1:Ni
        for kk = 1 : Nx
            Y0(NNB(ii, kk), 1) = f(x(kk));
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    %%% meningful initial data %%%
    % we first choose the initial position of the centerline p0(x)
    %             and the initial rotation matrix R0(x)
    % from this we build z0(x)
    % then we choose v0 in such a way that it fulfills the compatibility
    % conditions at x=0 and x=ell
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% building the initial data (p0, R0) for the GEB model %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    syms eta;
     
    p0_syms = 1/sqrt(2)*[eta; (1-cos(eta)); sin(eta)];
    R0_syms = [[1/sqrt(2), 0, -1/sqrt(2)]; 
               [sin(eta)/sqrt(2), cos(eta), sin(eta)/sqrt(2)];
               [cos(eta)/sqrt(2), -sin(eta), cos(eta)/sqrt(2)]];
%     dp0_syms = diff(p0_syms, eta); % test
%     dR0_syms = diff(R0_syms, eta); % test


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Below are some tests to find other initial data %%%
    %%% --------------------------------------------------------------- %%%
%     p0_syms = 1/sqrt(2)*[eta; -cos(eta); sin(eta)];
% %     dp0_syms = diff(p0_syms, eta);
%     R0_syms = [[1/sqrt(2), 0, -1/sqrt(2)]; 
%                [sin(eta)/sqrt(2), cos(eta), sin(eta)/sqrt(2)];
%                [cos(eta)/sqrt(2), -sin(eta), cos(eta)/sqrt(2)]];
% %     dR0_syms = diff(R0_syms, eta);
    %%% --------------------------------------------------------------- %%%
%     p0_syms = 1/(sqrt(2)*3)*[eta*3; (1-cos(eta*3)); sin(eta*3)];
%     dp0_syms = diff(p0_syms, eta);
%     R0_syms = [[1/sqrt(2), 0, -1/sqrt(2)]; 
%                [sin(eta*3)/sqrt(2), cos(eta*3), sin(eta*3)/sqrt(2)];
%                [cos(eta*3)/sqrt(2), -sin(eta*3), cos(eta*3)/sqrt(2)]];
%     dR0_syms = diff(R0_syms, eta);
    %%% --------------------------------------------------------------- %%%
%     p0_syms = 1/(sqrt(2)*3)*[eta*3; -cos(3*eta); sin(3*eta)];
% %     dp0_syms = diff(p0_syms, eta);
%     R0_syms = [[1/sqrt(2), 0, -1/sqrt(2)]; 
%                [sin(3*eta)/sqrt(2), cos(3*eta), sin(3*eta)/sqrt(2)];
%                [cos(3*eta)/sqrt(2), -sin(3*eta), cos(3*eta)/sqrt(2)]];
% %     dR0_syms = diff(R0_syms, eta);

    % not param by arc length
%     p0_syms = [eta; 1-cos(3*eta); sin(3*eta)];
%     R0_syms = 1/sqrt(2)*[[1, 0, -1];...
%         [sin(3*eta), sqrt(2)*cos(3*eta), sin(3*eta)];...
%         [cos(3*eta), -sqrt(2)*sin(3*eta), cos(3*eta)]];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% building the initial data y0 for the IGEB model %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Gamma0 = [0; 0; 0];    % there is no initial shear
    W0hat = [0, -1/sqrt(2), 0; 1/sqrt(2), 0, 1/sqrt(2); 0, -1/sqrt(2), 0];
    Upsilon0 = func_vec(W0hat);
    
    % not param by arc length
%     Gamma0 = (3*sqrt(2)-1)*[1; 0; 0]; % test
%     Upsilon0 = 3/sqrt(2)*[-1; 0; 1];  % test 
    
    z0 = [Gamma0; Upsilon0];  % strains
    z0 = flexMat\z0;          % corresponding stresses
    v0_ell = - kappa\z0;      % velocities at ell fulfilling 
                              % compatibility conditions at x = ell

    % to also fulfill compatibility conditions at x=0 we do as follows:
    
    %%% ------- imposes zero order comp. cond. ------- %%%
%     y0Mat = zeros(12, Nx);
%     for kk = 1:Nx
%         for ii = 1:6
%             y0Mat(ii, kk) = v0_ell(ii)*x(kk)/ell;
%         end
%         y0Mat(7:12, kk) = z0;
%     end
    %%% ---------------------------------------------- %%%
    
    %%% imposes zero-order comp. cond + null acceleration at x=0 %%%
    y0Mat = zeros(12, Nx);
    if zeroOrderBC == true
        for ii = 1:6
            x_constr = [0, 0.05, ell-0.05, ell];
            y_constr = [0, 0, v0_ell(ii), v0_ell(ii)];
            v0_interp = pchip(x_constr, y_constr, x);
            for kk = 1:Nx
                y0Mat(ii, :) = v0_interp;
            end
        end
    end
    for kk = 1:Nx
        y0Mat(7:12, kk) = z0;
    end
    %%% ---------------------------------------------- %%%
    
    %%% --------------------------------------------------------------- %%%
    
    %%% transform y0Mat to a vector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y0 = zeros(Ntot, 1);
    for ii = 1:Ni
        for kk = 1 : Nx
            Y0(NNB(ii, kk), 1) = y0Mat(ii, kk);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%% plot of the initial data y0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plot_y0
        y0_max = max(max(Y0));
        y0_min = min(min(Y0));
                
        f_y0 = figure();
        ii2subplot = [1, 5, 9, 2, 6, 10, 3, 7, 11, 4, 8, 12];
        ii2title = ["Linear velocity", "Angular velocity", "Forces", "Moments"];
        ii2label = ["$V_1$","$V_2$","$V_3$","$W_1$","$W_2$","$W_3$",...
            "$\Phi_1$","$\Phi_2$","$\Phi_3$","$\Psi_1$","$\Psi_2$","$\Psi_3$"];
        for ii = 1:12

            subplot(3, 4, ii2subplot(ii));
            s = plot(x, y0Mat(ii, :), 'LineWidth', size_line);
            ylabel(ii2label(ii),'Interpreter','latex');
            xlabel('$x$','Interpreter','latex'); 

            if ii2subplot(ii) <= 4
                title(ii2title(ii2subplot(ii)),'Interpreter','latex', 'fontsize',11);
            end
            %ylim([y0_min, y0_max]);
            grid on;
         end
        exportgraphics(f_y0,name_ini,'ContentType','vector')
    end        
    %%% --------------------------------------------------------------- %%%
    
    %%%% p0 and R0 evaluated over spatial points %%%%%%%%%%%%%%%%%%%%%%%%%%
    p0 = zeros(3, Nx);       % position of the centerline at t = 0
    R0 = zeros(3, 3, Nx);    % orientation of cross sections at t = 0
    for kk=1:Nx
        p0(:, kk) = subs(p0_syms, x(kk));
        R0(:, :, kk) = subs(R0_syms, x(kk));
    end
%     q0_true = zeros(4, Nx); % test
%     for kk=1:Nx
%         q0_true(:, kk) = rotm2quat(R0(:, :, kk));
%     end
    %%% --------------------------------------------------------------- %%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% test if we can recover p0, R0 from z0 via transfo %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hx = x(2) - x(1);    % spatial step
    z0 = y0Mat(7:12, :); % forces and moments
    s0 = flexMat*z0;     % strains
    z01 = s0(1:3, :);    % strains first 3 components
    z02 = s0(4:6, :);    % strains last 3 components

    q0_b = zeros(4, Nx);                    % quaternion
    q0_b(:, 1) = rotm2quat(R0(:, :, 1));    % 'ini' data at x=0
    R0_b = zeros(3, 3, Nx);                 % rotation
    R0_b(:, :, 1) = quat2rotm(q0_b(:, 1)'); % 'ini' data at x=0
    
%     for kk=1:Nx-1          % euler explicit: find the quaternions
%         q0_b(:, kk+1) = ( eye(4) + hx * func_U( z02(:, kk)) )*q0_b(:, kk);  
%         R0_b(:, :, kk+1) = quat2rotm(q0_b(:, kk+1)'); % corresp. rotation
%     end

    for kk=1:Nx-1          % mid point: find the quaternions
        Axt = func_U( z02(:, kk)/2 + z02(:, kk+1)/2 + kap ); 
        q0_b(:, kk+1) = ( eye(4) - hx/2 * Axt )\(( eye(4) + hx/2*Axt)*q0_b(:, kk));  
        R0_b(:, :, kk+1) = quat2rotm(q0_b(:, kk+1)'); % corresp. rotation
    end
    %sum(sum(sum(R0 - R0_b)))  % test
    
    p0_b = zeros(3, Nx);       % position of centerline
    p0_b(:, 1) = p0(:, 1);     % 'ini' data at x=0
    MM = zeros(3, Nx); 
    MM(:, 1) = R0_b(:, :, 1)*(z01(:, 1)+[1; 0; 0]); 
    for kk = 2:Nx         % integrate to find position of centerline      
        MM(:, kk) = R0_b(:, :, kk)*(z01(:, kk)+[1; 0; 0]);
        idx = 1:kk;
        p0_b(1, kk) = p0_b(1, 1) + trapz(x(idx), MM(1, idx), 2);
        p0_b(2, kk) = p0_b(2, 1) + trapz(x(idx), MM(2, idx), 2);
        p0_b(3, kk) = p0_b(3, 1) + trapz(x(idx), MM(3, idx), 2);
    end

    %%% plot of p0 and p0_b:
%     figure();
%     plot3(p0(1, :), p0(2, :), p0(3, :), 'lineWidth', 2); hold on;
%     grid on;
%     plot3(p0_b(1, :), p0_b(2, :), p0_b(3, :), ':r', 'lineWidth', 2);

    %%% plot of p0:
    f_p0 = figure();
    plot3(p0(1, :), p0(2, :), p0(3, :), 'lineWidth', 2);
    view(gca, 96,10);
    axis equal;
    grid on;
    xlabel('X', 'Interpreter', 'latex');
    ylabel('Y', 'Interpreter', 'latex');
    zlabel('Z', 'Interpreter', 'latex');
    title('Initial position of the centerline', 'Interpreter', 'latex', 'fontsize', 12);
    exportgraphics(f_p0,'fig/initial_centerline.pdf','ContentType','vector')
    %%% --------------------------------------------------------------- %%%
    
end

toc

%% Initialization of the state
Y_exp = zeros(Nf, Nt);     % zero matrix
Y_exp(:, 1) = Y0(dof, 1);  % set the initial data at time t1

%% Solving the ODE via implicit midpoint rule

if linearized
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        linearized IGEB           %
    % Approx equation: M y_t + K y = 0 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    M1 = M + 1/2 * ht * K; 
    M2 = M - 1/2 * ht * K; 
    [LL,UU,PP,QQ] = lu(M1); % LU decomposition of sparse matrix M %
                            % i.e. PMQ=LU where P and Q are unitary
                            % Thus, Mx = y equiv x = Q (U\(L\(Py)))
    disp('Solving the linearized system..')
    tic
    for kk=1:Nt-1
        Y_exp(:, kk+1) = QQ * ( UU\( LL\( PP * M2 * Y_exp(:, kk) ) )); 
    end
    toc   

else 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                  IGEB                    %
    % Approx equation: M y_t + K y + Q(y)y = 0 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Yk=Y_exp(:, 1);
    %[Q_yk, Qdagger_yk] = Q(Yk, P1, P2, P3, Pdagger1, Pdagger2, Pdagger3, Ni, Ne, Nf, NNBc);
    
    %%% M.A. code %%%
    tolzero=1e-12;  %-12
    reltolX=1e-6;   %-6
    tolF=1e-15;     %-15
    %%%%%%%%%%%%%%%%%
    
    disp('Solving the semilinear system..')
    tic
    for kk=1:Nt-1  % loop over time
        Yk = Y_exp(:, kk);
        [Q_yk, Qdagger_yk] = Q(Yk, P1, P2, P3, Pdagger1, Pdagger2, Pdagger3, Ni, Ne, Nf, NNBc);
    %     Q_yk = zeros(Nf); Qdagger_yk = zeros(Nf);   % for test 
        
        zm = Yk; % initialize zm
        while 1 
            %%% Newton method: the scheme reads %%%
            % zm1 = zm - (JacFk(zm))^{-1} Fk(zm)  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [Q_zm, Qdagger_zm] = Q(zm, P1, P2, P3, Pdagger1, Pdagger2, Pdagger3, Ni, Ne, Nf, NNBc);
    %         Q_zm = zeros(Nf);          % for test
    %         Qdagger_zm = zeros(Nf);    % for test

            Fk_zm = (M + ht/2*K)*zm - (M - ht/2*K)*Yk + ht/4*(Q_yk*Yk + (Q_yk + Qdagger_yk)*zm + Q_zm*zm);
            JacFk_zm = M + ht/2*K + ht/4*(Q_yk + Qdagger_yk) + ht/4*(Q_zm + Qdagger_zm);
            zm1 = zm - JacFk_zm\Fk_zm; 
            
            %%% M.A. code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            rel_err = (zm1 - zm)./zm;                                     %
            nan_or_inf = find( isnan(rel_err) + isinf(rel_err) ...        %
                + (abs(zm)<=tolzero) );                                   %
            rel_err(nan_or_inf) = 0;                                      %
            if (norm(rel_err,inf) <= reltolX) && (norm(Fk_zm,inf) <= tolF)%
                break                                                     %
            end                                                           %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
            
            zm = zm1;
        end
        Y_exp(:, kk+1) = zm1;
    end
    toc
end

%% Full state including the Dirichlet nodes
Y = zeros(Ntot, Nt);
Y(dof, :) = Y_exp(:, :);
ymax = max(max(Y)); ymin = min(min(Y));

%% True approximation of y ?
% disp('Building the approximation of y..')
% tic
% 
% %xNew = linspace(0, ell, 2000); % different spatial discretization
% xNew = x; % same spatial discretization
% y = zeros(Nx, Nt, Ni);
% for ii = 1:Ni
%     for nn = 1:Nt
%         for kk = 1:numel(xNew)
%             y(kk, nn, ii) = func_yi(xNew(kk), ii, nn, Y, Nx, NNB,Ne,he, x);
%         end
%     end
% end
% 
% truc = 0;
% for ii=1:12
%     truc = truc + sum(sum(Y(NNB(ii, :), :)-y(:, :, ii)));
% end
% disp(['truc = ', num2str(truc)])
% 
% toc

%% Plot the solution
if plot_y
    disp('Plots of the solution y..')

    %[t_mat, x_mat] = meshgrid(t, x);

    f = figure();
    set(gcf,'Position',[100 100 1200 600])
    ii2subplot = [1, 5, 9, 2, 6, 10, 3, 7, 11, 4, 8, 12];
    ii2title = ["Linear velocity", "Angular velocity", "Forces", "Moments"];
    ii2label = ["$V_1$","$V_2$","$V_3$","$W_1$","$W_2$","$W_3$",...
        "$\Phi_1$","$\Phi_2$","$\Phi_3$","$\Psi_1$","$\Psi_2$","$\Psi_3$"];
    for ii = 1:12

        subplot(3, 4, ii2subplot(ii));
        s = surf(x, t, Y(NNB(ii, :), :)');
        %contourf(t_mat, x_mat, Y(NNB(ii, :), :));
    % % %     s = surf(t, xNew, y(:, :, ii));        % 'true approx' of y 
        s.EdgeColor = 'none'; 
        ylabel([ii2label(ii),'\ \ $t$'],'Interpreter','latex'); 
        xlabel('$x$','Interpreter','latex'); 
        %zlabel(ii2label(ii),'Interpreter','latex');
        colorbar
        view(2)
        
        if ii2subplot(ii) <= 4
            title(ii2title(ii2subplot(ii)),'Interpreter','latex', 'fontsize', 12);
        end
    end
    orient(f,'landscape')
    exportgraphics(f,file_name_y,'ContentType','vector')
end

%% recover position of the centerline using velocities
space_solve = false;  % true: we recover the centerline via the strains
                      % false: we recover the centerline via the velocities
if space_solve == true
    type_centerline = 'XSolve';
    title_centerline = 'Centerline: by space integration using $z$';
    title_arclenght = 'Arclength: by space integration using $z$';
else
    type_centerline = 'TSolve';
    title_centerline = 'Centerline: by time integration using $v$';
    title_arclength = 'ArcLength: by time integration using $v$';
end
file_name_centerline = ['CENTERL_len', num2str(ell), '_', linNonlin, '_', type_K, '_', type_centerline, '_', zero_ord, '.pdf'];
file_name_arclength = ['fig/ARCLEN_len', num2str(ell), '_', linNonlin, '_', type_K, '_', type_centerline,'_', zero_ord, '.pdf'];

%%% WITH TIME : using the first 6 eq. of the transformation %%%
disp('Recovering the position of the centerline (using the velocities)..')
tic

H = zeros(4, Nt, Nx);     % quaternions
R = zeros(3, 3, Nt, Nx);  % rotation matrices
p = zeros(3, Nt, Nx);     % positions
for kk = 1:Nx                                 % initial data at t=0
    H(:, 1, kk) = rotm2quat(R0(:, :, kk))';   % quat angle at t=0
    R(:, :, 1, kk) = R0(:, :, kk);            % angle at t=0
    p(:, 1, kk) = p0(:, kk);                  % position at t=0
end
% we have y(x,t) containing the velocities and strains
% here we will just use the velocities
% the velocoties are the first six components of y
velocities = zeros(6, Nt, Nx); 
for kk = 1:Nx               % for all x
    for nn = 1:Nt-1         % for all t
        for ii = 1:6
            velocities(ii, nn, kk) = Y(NNB(ii, kk), nn); % extract velocities
        end
    end
end

if centerline_t_scheme == 0       %%% mid point rule
    for kk = 1:Nx                 % for all x
        for nn = 1:Nt-1           % for all t
            Axt = func_U( (velocities(4:6, nn, kk)+velocities(4:6, nn+1, kk) )/2);  % uses angular velocities  
            H(:, nn+1, kk) = (eye(4)-ht/2*Axt)\((eye(4) + ht/2*Axt)*H(:,nn, kk)); % GOOD
            %disp(norm(H(:, nn+1, kk)))  % test
            %R(:, :, nn+1, kk) = quat2rotm(H(:, nn+1, kk)');
            R(:, :, nn+1, kk) = func_quat2rotm(H(:, nn+1, kk));
        end
    end
elseif centerline_t_scheme == 1   %%% explicit Euler (does not work properly)
    for kk = 1:Nx                 % for all x
        for nn = 1:Nt-1           % for all t
            Axt = func_U( velocities(4:6, nn, kk));  % uses angular velocities  
            H(:, nn+1, kk) = (eye(4) + ht*Axt)*H(:, nn, kk);
            %R(:, :, nn+1, kk) = quat2rotm(H(:, nn+1, kk)');
            R(:, :, nn+1, kk) = func_quat2rotm(H(:, nn+1, kk));
        end
    end
elseif centerline_t_scheme == 2   %%% ode 45
    tspan = t;
    opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
for kk = 2:Nx
    ic = H(:, 1, kk); 
    [t_ode,y_ode] = ode45(@(tau,q) myode(tau,q,velocities(4:6, :, kk), tspan), tspan, ic, opts);
    H(:, :, kk) = transpose(y_ode);
    for nn = 1:Nt
        %R(:, :, nn, kk) = quat2rotm(H(:, nn, kk)');
        R(:, :, nn, kk) = func_quat2rotm(H(:, nn, kk));
    end
end
elseif centerline_t_scheme == 3   %%% implicit euler
    for kk = 1:Nx                 % for all x
        for nn = 1:Nt-1           % for all t
            Axt = func_U( velocities(4:6, nn+1, kk));  % uses angular velocities  
            %abs(eig((eye(4) - ht*Axt)\eye(4)))   % test
            H(:, nn+1, kk) = (eye(4) - ht*Axt)\H(:, nn, kk);
            %R(:, :, nn+1, kk) = quat2rotm(H(:, nn+1, kk)');
            R(:, :, nn+1, kk) = func_quat2rotm(H(:, nn+1, kk));
        end
    end
else                              %%% Z. paper (does not work properly)
    for kk = 1:Nx                 % for all x
        for nn = 1:Nt-1           % for all t
            Wm = (velocities(4:6, nn, kk)+velocities(4:6, nn+1, kk) )/2;
            AAA = (func_U_right(Wm) + eye(4)*4/ht)*4*ht^2/(16+ht^2*norm(Wm)^2);            
            qqq0 = H(1, nn, kk);
            qqq = H(2:4, nn, kk);
            bbb = [1/ht*qqq0 - 1/4*transpose(qqq)*Wm;...
                1/ht*qqq + 1/4*(qqq0*Wm + hat(Wm)*qqq)];
            H(:, nn+1, kk) = AAA*bbb/norm(AAA*bbb);
            %disp(norm(H(:, nn+1, kk)))   % test
            R(:, :, nn+1, kk) = func_quat2rotm(H(:, nn+1, kk));
        end
    end    
end

Bxt = zeros(3, Nt, Nx); 
for kk = 1:Nx
    Bxt(:, 1, kk) = R(:, :, 1, kk)*(velocities(1:3, 1, kk));
    for nn = 2:Nt
        idx = 1:nn;
        Bxt(:, nn, kk) = R(:, :, nn, kk)*(velocities(1:3, nn, kk));
        p(1, nn, kk) = p(1, 1, kk) + trapz(t(idx), Bxt(1, idx, kk), 2);
        p(2, nn, kk) = p(2, 1, kk) + trapz(t(idx), Bxt(2, idx, kk), 2);
        p(3, nn, kk) = p(3, 1, kk) + trapz(t(idx), Bxt(3, idx, kk), 2);
    end
end
p = permute(p, [1, 3, 2]); % for plotting we need to change the order of the space and time indexes

toc
    
%%% ----------- plot arclength ----------- %%%
arcLength_time = zeros(Nt, 1);
for nn = 1:Nt
    arcLength_time(nn, 1) = trapz(x, sqrt(sum((gradient(p(:, :, nn), hx)).^2))); % arclength
end
f_arc = figure();
plot(t, arcLength_time, 'lineWidth', 2);
xlabel('$t$','Interpreter','latex');
grid on
title(title_arclength, 'Interpreter', 'latex');
%print(f_arc, file_name_arclength,'-dpdf')
exportgraphics(f_arc,file_name_arclength,'ContentType','vector');
%%% -------------------------------------- %%%

%%% ----------- plot centerline ----------- %%%
pf = zeros(3, Nx); % the undeformed configuration fulfilling the clamped BC
temp = R0(:, 1, 1);
for kk = 1:Nx
    pf(:, kk) = p0(:, 1) + x(kk)*temp;
end
if plot_centerline
    for ii = 1: size(orient_cent, 1)
        f_c = figure();
        h = plot3(p0(1, :), p0(2, :), p0(3, :), 'b', 'lineWidth', 3); hold on;
        %h.Color(4) = 0.5;
        grid on;
        h = plot3(p(1, :, 1), p(2, :, 1), p(3, :, 1), ':b', 'lineWidth', 2);
        h.Color(4) = 0.5;
        nn_disp = 2:p_plot:Nt;
        for nn = nn_disp
            h = plot3(p(1, :, nn), p(2, :, nn), p(3, :, nn), ':b', 'lineWidth', 2);
            h.Color(4) = 0.5;
            set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
        plot3(p(1, :, end), p(2, :, end), p(3, :, end), 'r', 'lineWidth', 3);
        plot3(pf(1, :), pf(2, :), pf(3, :), '--g', 'lineWidth', 3);
        axis equal;
        view(gca, orient_cent(ii, 1), orient_cent(ii, 2));
        legend('$t=0$', '$t$', '$t = T$',...
            'undeformed','Interpreter','latex', 'Location', locLegend(ii), 'fontsize', 12);
        xlabel('X', 'Interpreter', 'latex');
        ylabel('Y', 'Interpreter', 'latex');
        zlabel('Z', 'Interpreter', 'latex');
        title(title_centerline, 'Interpreter', 'latex', 'fontsize', 12);
        exportgraphics(f_c,['fig/Orient', num2str(ii), '_', file_name_centerline],'ContentType','vector')
    end
end
%%% -------------------------------------- %%%

%% recover position of the centerline using internal forces and moments
disp('Recovering the position of the centerline (using the internal forces and moments / strains)..')
tic

space_solve = true;  % true: we recover the centerline via the strains

if space_solve == true
    type_centerline = 'XSolve';
    title_centerline = 'Centerline: by space integration using $z$';
    title_arclenght = 'Arclength: by space integration using $z$';
else
    type_centerline = 'TSolve';
    title_centerline = 'Centerline: by time integration using $v$';
    title_arclength = 'ArcLength: by time integration using $v$';
end
file_name_centerline = ['CENTERL_len', num2str(ell), '_', linNonlin, '_', type_K, '_', type_centerline, '_', zero_ord, '.pdf'];
file_name_arclength = ['fig/ARCLEN_len', num2str(ell), '_', linNonlin, '_', type_K, '_', type_centerline, '_', zero_ord, '.pdf'];

%%% WITH SPACE: using the last 6 eq. of the transformation  %%%
H = zeros(4, Nx, Nt);     % quaternions
R = zeros(3, 3, Nx, Nt);  % rotation matrices
p = zeros(3, Nx, Nt);     % positions
for nn = 1:Nt                                 % clamped at all times
    H(:, 1, nn) = rotm2quat(R0(:, :, 1))';    % quat angle at x=0
    R(:, :, 1, nn) = R0(:, :, 1);             % angle at x=0
    p(:, 1, nn) = p0(:, 1);                   % position at x=0
end
% we have y(x,t) containing the velocities and strains
% here we will just use the strains
% the strains are the last six components of y left-multiplied by flexMat
strains = zeros(6, Nx, Nt); 
for nn = 1:Nt            % for all t
    for kk = 1:Nx-1      % for all x
        forces = zeros(6, 1);  
        for ii = 7:12
            forces(ii-6, 1) = Y(NNB(ii, kk), nn);    % extract the forces
        end
        strains(:, kk, nn) = flexMat*forces;         % left-mult by flexMat
    end
end

if centerline_x_scheme == 0  %%% mid point rule
    for nn = 1:Nt            % for all t
        for kk = 1:Nx-1      % for all x
            Axt = func_U(strains(4:6, kk, nn)/2 + strains(4:6, kk+1, nn)/2 + kap );
            H(:, kk+1, nn) = (eye(4) - hx/2*Axt)\((eye(4) + hx/2*Axt)*H(:, kk, nn));
            %R(:, :, kk+1, nn) = quat2rotm(H(:, kk+1, nn)');
            R(:, :, kk+1, nn) = func_quat2rotm(H(:, kk+1, nn));
        end
    end
else                         %%% euler explicit
    for nn = 1:Nt            % for all t
        for kk = 1:Nx-1      % for all x
            Axt = func_U( strains(4:6, kk, nn) + kap );  
            H(:, kk+1, nn) = (eye(4) + hx*Axt)*H(:, kk, nn);
            %R(:, :, kk+1, nn) = quat2rotm(H(:, kk+1, nn)');
            R(:, :, kk+1, nn) = func_quat2rotm(H(:, kk+1, nn));
        end
    end
end

Bxt = zeros(3, Nx, Nt); 
for nn = 1:Nt
    Bxt(:, 1, nn) = R(:, :, 1, nn)*(strains(1:3, 1, nn) + [1; 0; 0]);
    for kk = 2:Nx
        idx = 1:kk;
        Bxt(:, kk, nn) = R(:, :, kk, nn)*(strains(1:3, kk, nn) + [1; 0; 0]);
        p(1, kk, nn) = p(1, 1, nn) + trapz(x(idx), Bxt(1, idx, nn), 2);
        p(2, kk, nn) = p(2, 1, nn) + trapz(x(idx), Bxt(2, idx, nn), 2);
        p(3, kk, nn) = p(3, 1, nn) + trapz(x(idx), Bxt(3, idx, nn), 2);
    end
end

toc

%%% ----------- plot arclength ----------- %%%
arcLength_time = zeros(Nt, 1);
for nn = 1:Nt
    %arcLength_time(nn, 1) = trapz(x(1:end-1),sqrt(sum((p(:, 2:end, nn)-p(:,1:end-1, nn)).^2)));
    %arcLength_time(nn, 1) = sum(sqrt(sum((p(:, 2:end, nn)-p(:,1:end-1, nn)).^2)));
    %arcLength_time(nn, 1) = trapz(x, sqrt(sum((strains(1:3, :, nn)+[1; 0; 0]).^2)));
    arcLength_time(nn, 1) = trapz(x, sqrt(sum((gradient(p(:, :, nn), hx)).^2))); % arclength
end
f_arc = figure();
plot(t, arcLength_time, 'lineWidth', 2);
xlabel('$t$','Interpreter','latex');
grid on
title(title_arclength, 'Interpreter', 'latex')
%print(f_arc, file_name_arclength,'-dpdf')
exportgraphics(f_arc,file_name_arclength,'ContentType','vector')
%%% -------------------------------------- %%%

%%% ----------- plot centerline ----------- %%%
pf = zeros(3, Nx); % the undeformed configuration fulfilling the clamped BC
temp = R0(:, 1, 1);
for kk = 1:Nx
    pf(:, kk) = p0(:, 1) + x(kk)*temp;
end
if plot_centerline
    for ii = 1: size(orient_cent, 1)
        f_c = figure();
        h = plot3(p0(1, :), p0(2, :), p0(3, :), 'b', 'lineWidth', 3); hold on;
        grid on;
        h = plot3(p(1, :, 1), p(2, :, 1), p(3, :, 1), ':b', 'lineWidth', 2);
        h.Color(4) = 0.5;
        nn_disp = 2:p_plot:Nt;
        for nn = nn_disp
            h = plot3(p(1, :, nn), p(2, :, nn), p(3, :, nn), ':b', 'lineWidth', 2);
            h.Color(4) = 0.5;
            set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
        plot3(p(1, :, end), p(2, :, end), p(3, :, end), 'r', 'lineWidth', 3);
        plot3(pf(1, :), pf(2, :), pf(3, :), '--g', 'lineWidth', 3);
        axis equal;
        xlabel('X', 'Interpreter', 'latex');
        ylabel('Y', 'Interpreter', 'latex');
        zlabel('Z', 'Interpreter', 'latex');
        title(title_centerline, 'Interpreter', 'latex', 'fontsize', 12);
        view(gca, orient_cent(ii, 1), orient_cent(ii, 2));
        legend('$t=0$', '$t$', '$t = T$',...
            'undeformed','Interpreter','latex', 'Location','northwest', 'fontsize', 12);
        exportgraphics(f_c,['fig/Orient', num2str(ii), '_', file_name_centerline],'ContentType','vector');
    end
end
%%% -------------------------------------- %%%

disp('End.')

%% test Nk definition
% % % figure();
% % % %kk=1;
% % % x_new = linspace(0, ell, 2000);
% % % for kk=1:Nx
% % %     Nk_vec = zeros(Nx, 1);
% % %     for mm = 1:numel(x_new)
% % %         Nk_vec(mm) = Nk(x_new(mm), kk,Ne,he, x);
% % %     end
% % %     plot(x_new, Nk_vec); hold on;
% % % end
