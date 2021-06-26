% Demo of SAED
%
% Reference:
% X. Wang, Y. Zhong, C. Cui, L. Zhang and Y. Xu, "Autonomous Endmember
% Detection via an Abundance Anomaly Guided Saliency Prior for
% Hyperspectral Imagery," in IEEE Transactions on Geoscience and Remote
% Sensing, vol. 59, no. 3, pp. 2336-2351, March 2021, doi:
% 10.1109/TGRS.2020.3001353.
% 
% X. Wang, Y. Zhong, Y. Xu, L. Zhang and Y. Xu, "Saliency-Based Endmember
% Detection for Hyperspectral Imagery," in IEEE Transactions on Geoscience
% and Remote Sensing, vol. 56, no. 7, pp. 3667-3680, July 2018, doi:
% 10.1109/TGRS.2018.2805340.
%
% Copyright (C) 2020, Xinyu Wang (wangxinyu@whu.edu.cn)
%                     Yanfei Zhong (zhongyanfei@whu.edu.cn)
%                     Wuhan university
%                     All rights reserved.%
%
% Last Modified:
% 24 June, 2021


%% Load pathes
clear
close all
clc
addpath(genpath('../code'));
addpath(genpath('../data'));
if exist('..\results','dir')==1
    rmdir('..\results','s');
end


%% Load simulated data
load 'D_Simu_A'; % Endmember data
load 'D_Simu_S'; % Abundance data

L = size(A,1); 
[row, col, M] = size(S); 
N = row * col; 
S = reshape(S, N, M)'; % colomnwise abundance matrix, where M is the ground-truth VD.
X = A * S; %  colomnwise data matrix
Y = reshape(X',row,col,L);  % 3D HSI row*col*L
%please normalize the radiances to a range of 0 - 1.0 for real data.
SNR = 30; % add noise
noise_type = 'additive'; eta = 0;
[X, n, Cn] = addNoise (Y,noise_type,SNR, eta, 1);
X = max (X,eps); %  colomnwise data matrix
Y = reshape(X',row,col,L);  % 3D HSI row*col*L


%% Compared method: ATGP(OSP)
A_Atgp = hyperAtgp(X, M);
S_Atgp = fcls(A_Atgp, X);
fprintf('\n The SAD accuracy of Atgp:\n');
SAD_Atgp = sam(A, A_Atgp); 


%% Input parameters of SAED
para.Y = Y;           % 3D HSI row*col*L
para.X = X;           % 2D HSI L*N_pixels
para.Maximumend = 20; % the maximum number of endmembers
para.r = 5;           % the size of superpixels

para.verbose = 1;     % for visualization
para.saveimage = 1;  % for saving intermediate results 
fprintf('\n The parameter setting of SAED:\n');disp(para);


%% Main program of SAED
A_Saed = SAED(para);     % main function 
S_Saed = fcls(A_Saed,X); % Abundances
fprintf('\n The SAD accuracy of SAED:\n');
SAD_Saed = sam(A,A_Saed); 

if para.saveimage
    S_Saed = reshape(S_Saed',row,col,size(A_Saed,2));
    for i = 1 : size(A_Saed,2)
       figure;axis off;axis image; set(gcf,'visible','off');
       imagesc(S_Saed(:,:,i),[0,1]);drawnow;
       saveas_center(gcf,['..\results','\','S_',num2str(i),'.tif'], row, col); 
    end
    figure; plot(A_Saed); set(gcf,'visible','off');drawnow;
    saveas(gcf,['..\results','\','A_Saed.tif'])
    save(['..\results\','A_Saed.mat'],'A_Saed');
    save(['..\results\','S_Saed.mat'],'S_Saed');
end


















