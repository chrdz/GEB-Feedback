function res = Gdagger(ii, M, C)
% Gdagger computes the matrix G_\dagger(\hat{e}_{ii})
    y = zeros(12, 1); y(ii, 1) = 1;
    v = y(1:6, 1); z = y(7:12, 1);
    res = zeros(12, 12); 
    res(1:6, 1:6) = -L2(M*v);
    res(1:6, 7:12) = -L1(C*z);
    res(7:12, 1:6) = transpose(L1(C*z));
end