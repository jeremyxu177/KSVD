%% image pre-processing and parameters setting
imageName = 'barbara.png';
I = imread(imageName);
I = imresize((I), [256 256]);
sigma = 100; % standard deviation of Gaussian white noise
variance = (sigma / 255) ^ 2;
I_noisy = imnoise(I, 'gaussian', 0, variance);
Y = double(I_noisy); %rescale to 0-1
% I_noisy = I + sigma*randn(size(I));
% imshow(I_noisy);
PSNRIn = psnr(I_noisy, I);
% PSNRIn = 10*log10(255/sqrt(mean((I_noisy(:)-I(:)).^2)));
bb = 8; %(patch width/height) block size
n = bb^2; %patch size
k = 256; % dictionary size
J = 10; % number of training iterations %however,smaller J will not affect much
lambda = 30 / sigma; % lagrange multiplier, given in the paper
C = 1.15; % noise gain
error = (C * sigma) ^ 2 * n; % error tolerance in OMP
Dict = DCT_dictionary(n, k); % initialize the dictionary
% max_iter = 60000;
[a b] = size(Y) % total number of pixels
total_N = a*b;
row_num = a - bb + 1; % number of patches in a row,should be 249
col_num = b - bb +1; % number of patches in a column,should be 249
patch_num = row_num * col_num; %total number of patches,should be 62001
[x_ind, y_ind] = meshgrid(1:row_num, 1:row_num); % x,y location of patch
y_hat = zeros(n, patch_num); % collect all patches where D * x_set = y_hat
M_all = zeros(size(Y)); % record the repetition of the patch at each pixel

for i = 1:patch_num
    x = x_ind(i); %x location
    y = y_ind(i); %y location
    tem = Y(x:(x + bb - 1), y:(y + bb - 1)); % record each 8*8 patch
    M_all(x:(x + bb - 1), y:(y + bb - 1)) = M_all(x:(x + bb - 1), y:(y + bb - 1)) + 1;
    y_hat(:, i) = tem(:); % collect the patch
end

%% K-SVD PART
%initialization
X = Y(:);
for iter = 1:J
    % Sparse coding stage
%     alpha = engn8535_omp(Dict, y_hat, error);
    alpha = engn8535_omp_test(Dict, y_hat, error);
    % Dictionary update stage
    R = y_hat - Dict * alpha; % residual
    for l = 1:k
        index = find(alpha(l, :) ~= 0); %find non-zero sparse coefficient
        if ~isempty(index) %if has non-zero entries
            E_l = R(:,index) + Dict(:,l)*alpha(l,index);
            % Compute svd for E_l
            [U,S,V] = svds(E_l);
            Dict(:,l) = U(:,1); %use only first column
            alpha(l,index) = S(1)*V(:,1)';
            % update the residual
            R(:,index) = E_l - Dict(:,l)*alpha(l,index);
        end
    end
end
%% compute clean image
M_all = lambda + M_all(:); % add weighting to the noisy image
recons_x = Dict * alpha; % reconstruct image patches using dictionary
x_hat = lambda * Y;
% combining patches together
for i = 1:patch_num
    x = x_ind(i);
    y = y_ind(i);
    block = reshape(recons_x(:, i), bb, bb);
    x_hat(x:(x + bb - 1), y:(y + bb - 1)) = x_hat(x:(x + bb - 1), y:(y + bb - 1)) + block;
end
result = x_hat(:) ./ M_all; % average the contructed image patches with the noisy image
% result = pinv(M_all)* x_hat(:);
% display the denoised image
clean_image = reshape(uint8(result), size(I_noisy));
% PSNROut = 10*log10(255/sqrt(mean((clean_image(:)-I(:)).^2)));
PSNROut= psnr(clean_image, I);
figure;
subplot(1,3,1); imshow(I,[]); title('Original image');
subplot(1,3,2); imshow(I_noisy,[]); title(strcat(['Noisy image,sigma=100,',num2str(PSNRIn),'dB']));
subplot(1,3,3); imshow(clean_image,[]); title(strcat(['Clean Image by Adaptive dictionary, ',num2str(PSNROut),'dB']));

% imshow(clean_image)
% imwrite(I_noise, 'noise_100.jpg')
% imwrite(denoised_image, 'denoised_100.jpg')
% Compare the PSNR before and after the denoising


%% DCT
% [IoutDCT,output] = denoiseImageDCT(I_noisy, sigma, k);
% PSNROut = 20*log10(255/sqrt(mean((IoutDCT(:)-I(:)).^2)));
% figure;
% subplot(1,3,1); imshow(I,[]); title('Original clean image');
% subplot(1,3,2); imshow(I_noisy,[]); title(strcat(['Noisy image, ',num2str(PSNRIn),'dB']));
% subplot(1,3,3); imshow(IoutDCT,[]); title(strcat(['Clean Image by DCT dictionary, ',num2str(PSNROut),'dB']));