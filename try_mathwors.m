imageName = 'lena.png';
I = imread(imageName);
I = im2double(I);
n = size(I,1);
sigma = 0.06;
M = I + sigma*randn(n);

w = 12;
m = 40*w^2;
x = floor( rand(1,1,m)*(n-w) )+1;
y = floor( rand(1,1,m)*(n-w) )+1;

[dY,dX] = meshgrid(0:w-1,0:w-1);
Xp = repmat(dX,[1 1 m]) + repmat(x, [w w 1]);
Yp = repmat(dY,[1 1 m]) + repmat(y, [w w 1]);
P = M(Xp+(Yp-1)*n);

P = P - repmat( mean(mean(P)), [w w] );

P = reshape(P, [w^2 m]);
clf;
plot_dictionnary(P, [], [8 12]);

% sigma = .06;
% n0 = 256;
% name = 'barb';
% f0 = rescale( crop(load_image(name),n0) );
% clf;
% imageplot(f0);
% f = f0 + sigma*randn(n0);
% clf;
% imageplot(f);
% if not(exist('D'))
%     sparsity_4_dictionary_learning;
% end
% w = sqrt(size(D,1));