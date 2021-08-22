function res = G(ii, M, C)
% G computes the matrix G(\hat{e}_{ii})
    y = zeros(12, 1); y(ii, 1) = 1;
    v = y(1:6, 1); z = y(7:12, 1);
    res = zeros(12, 12); 
    res(1:6, 1:6) = L1(v)*M;
    res(1:6, 7:12) = L2(z)*C;
    res(7:12, 7:12) = -transpose(L1(v))*C;
end