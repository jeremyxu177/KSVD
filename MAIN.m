

clear
bb=8; % block size
% RR=4; % redundancy factor
K=256; % number of atoms in the dictionary
sigma = 25; 
imageName = 'barbara.png';
I = imread(imageName);
I = im2double(I);
I = I*255;
I_noisy = I + sigma*randn(size(I));

SNR_in = 20*log10(255/sqrt(mean((I_noisy(:)-I(:)).^2)));

%% DCT
[IoutDCT,output] = denoiseImageDCT(I_noisy, sigma, K);
SNR_out = 20*log10(255/sqrt(mean((IoutDCT(:)-I(:)).^2)));
figure;
subplot(1,3,1); imshow(I,[]); title('Original clean image');
subplot(1,3,2); imshow(I_noisy,[]); title(strcat(['Noisy image, ',num2str(SNR_in),'dB']));
subplot(1,3,3); imshow(IoutDCT,[]); title(strcat(['Clean Image by DCT dictionary, ',num2str(SNR_out),'dB']));
% figure;
% I = displayDictionaryElementsAsImage(output.D, floor(sqrt(K)), floor(size(output.D,2)/floor(sqrt(K))),bb,bb,0);
% title('The DCT dictionary');

%% Global dictionary
[IoutGlobal,output] = denoiseImageGlobal(I_noisy, sigma,K);

SNR_out = 20*log10(255/sqrt(mean((IoutGlobal(:)-I(:)).^2)));
figure;
subplot(1,3,1); imshow(I,[]); title('Original clean image');
subplot(1,3,2); imshow(I_noisy,[]); title(strcat(['Noisy image, ',num2str(SNR_in),'dB']));
subplot(1,3,3); imshow(IoutGlobal,[]); title(strcat(['Clean Image by Global Trained dictionary, ',num2str(SNR_out),'dB']));
% figure;
% I = displayDictionaryElementsAsImage(output.D, floor(sqrt(K)), floor(size(output.D,2)/floor(sqrt(K))),bb,bb);
% title('The dictionary trained on patches from natural images');

%% Adaptive
[IoutAdaptive,output] = denoiseImageKSVD(I_noisy, sigma,K);

SNR_out = 20*log10(255/sqrt(mean((IoutAdaptive(:)-I(:)).^2)));
figure;
subplot(1,3,1); imshow(I,[]); title('Original clean image');
subplot(1,3,2); imshow(I_noisy,[]); title(strcat(['Noisy image, ',num2str(SNR_in),'dB']));
subplot(1,3,3); imshow(IoutAdaptive,[]); title(strcat(['Clean Image by Adaptive dictionary, ',num2str(SNR_out),'dB']));

% figure;
% I = displayDictionaryElementsAsImage(output.D, floor(sqrt(K)), floor(size(output.D,2)/floor(sqrt(K))),bb,bb);
% title('The dictionary trained on patches from the noisy image');