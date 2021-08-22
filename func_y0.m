function res = func_y0(dp0, R0, dR0, p1, w0, kap, Nx, flexMat, kappa)
% func_yO creates the initial data y0 for the IGEB model
% from the initial data of the GEB model.
    e1 = [1; 0; 0];
    res = zeros(12, Nx);
    
    for kk = 1:Nx        
        % velocities:
        %res(1:3, kk) = R0(:, :, kk)'*p1(:, kk);
        %res(4:6, kk) = R0(:, :, kk)'*w0(:, kk);
        
        % set strains:        
        res(7:9, kk) = R0(:, :, kk)'*dp0(:, kk) - e1;       
        res(10:12, kk) = func_vec(R0(:, :, kk)'*dR0(:, :, kk)) - kap;
        
        % convert strains to stresses:
        res(7:12, kk) = flexMat\res(7:12, kk);       
    end    
    for kk = 1:Nx
        % velocities constant and fulfilling C0 BC
        res(1:6, kk) = - kappa\res(7:12, Nx);
    end
end