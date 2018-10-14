function x_set = engn8535_omp_test(D, b_set, error)
l2_norm = @(x) sqrt(x'*x);
p = size(D, 2);
if (nargin < 2)
    error('not enough inputs');
elseif (nargin < 3)
    error = 1e-3;
end

for tmpPoint = 1:size(b_set,2)
    b = b_set(:,tmpPoint);
    r = b;
    V = [ ];

    while (l2_norm(r) >= error)
        temp = D'*r;
        [~, tmp_idx] = max(abs(temp));
        V = [V,tmp_idx(1)];
        x(V) = D(:, V) \ b;
        r = b - D*x;
    end
    x = zeros(p, 1);
    x_set(:,tmpPoint) = x;
%     err_set(:,tmpPoint) = err;
end
end