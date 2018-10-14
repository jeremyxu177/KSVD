function x_set = engn8535_omp(D, b_set, error)
p = size(D, 2);
if (nargin < 2)
    error('not enough inputs');
elseif (nargin < 3)
    error = 1e-3;
end
num = size(b_set,2);
x_set = zeros(p, num);
for tmpPoint = 1:num
    b = b_set(:,tmpPoint);
    r = b;
    V = [];
    while ((r'*r) > error)
        proxy = D'*r;
        [~, tmp_idx] = max(abs(proxy));
        V = [V,tmp_idx(1)];
        x = D(:, V) \ b;
        r = b - D(:, V) * x;
    end
    x_final = zeros(p, 1);
    x_final(V) = x;
    x_set(:,tmpPoint) = x_final;
end
end
