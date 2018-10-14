function displayDict(D, n, k)
sqrt_n = sqrt(n);
image_patches = zeros(sqrt(n) * sqrt(k));
count = 1;
for i = 1:sqrt(k)
for j = 1:sqrt(k)
V = D(:, count);
V = (V - min(V)) / (max(V) - min(V));
image_patches((sqrt_n * (i - 1) + 1): (sqrt_n * i), (sqrt_n * (j - 1) + 1): (sqrt_n * j)) = reshape(V, sqrt(n), sqrt(n));
count = count + 1;
end
end
imshow(image_patches)
imwrite(image_patches, 'Dict.jpg')
end