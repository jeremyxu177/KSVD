function Dictionary = DCT_dictionary(n, K)
V = zeros(sqrt(n), sqrt(K));
Pn = sqrt(K);
Dictionary = zeros(n, K);
for k = 0:(Pn-1)
    temp = (0:(sqrt(n)-1)) + 1;
    V(:, k + 1) = cos(temp * k * pi/ (Pn)); %compute cos term for each k
end
count = 0;
for i = 1:sqrt(K)
    for j = 1:sqrt(K)
%         tem = kron(V(:,j),V(:,j)');
        temp = V(:, j) * V(:, i)'; %kroneker delta
        Dictionary(:, count + 1) = temp(:);
        count = count + 1;
    end
end
Dictionary = normc(Dictionary);
end

