function A = SAED (para)
%==========================================================================
% Input:  para.X                -  L * N hyperspectrtal data;
%         para.Y                -  row * col * N hyperspectrtal data; 
%         para.Maximumend       -  the maximum number of endmembers
%         para.verbose          -  1- display mid-results, and 0 - otherwise  
%         para.saveimage        -  1- save mid-results, and 0 - otherwise
%         %
% 
% Output: A  -  L * () estimated endmember matrix obtained by SAED
% 
%
%
% Reference:
% X. Wang, Y. Zhong, C. Cui, L. Zhang and Y. Xu, "Autonomous Endmember
% Detection via an Abundance Anomaly Guided Saliency Prior for
% Hyperspectral Imagery," in IEEE Transactions on Geoscience and Remote
% Sensing, vol. 59, no. 3, pp. 2336-2351, March 2021, doi:
% 10.1109/TGRS.2020.3001353.
%
%X. Wang, Y. Zhong, Y. Xu, L. Zhang and Y. Xu, "Saliency-Based Endmember
%Detection for Hyperspectral Imagery," in IEEE Transactions on Geoscience
%and Remote Sensing, vol. 56, no. 7, pp. 3667-3680, July 2018, doi:
%10.1109/TGRS.2018.2805340.
%
%
% Tips:
%
% Copyright (C) 2020, Xinyu Wang (wangxinyu@whu.edu.cn)
%                     Yanfei Zhong (zhongyanfei@whu.edu.cn)
%                     Wuhan university
%                     All rights reserved.%
%
% Last Modified:
% 24 June, 2021



    X = para.X; 
    Y = para.Y;
    [row,col,L] =size(Y);
    [~,N] =size(X);

    %% 1: Generate superpixels
    [l, nk ,~] = slic_HSI(Y, ceil(N /para.r^2), 0.1);
    l = reshape(l,1,N); % labels of SLIC superpixels
    P = length(nk);     % the number of SLIC superpixels
    IR = -ones(size(l));
    fprintf('    iter      tol       C_max        D_win      C/D_ratio     nk\n');

    for k = 1 : para.Maximumend
        

       %% 2: Initialize A 
        if k == 1
            [~,disindex]= sort((sum(X.^2)),'descend');
             A(:,1) = X(:,disindex(1));
             A(:,2) = X(:,disindex(N));
        end


       %% 3: Generate AA matrix
            S = (A.'*A) \A.' * X;
            V_k = max(S-1,0)+ max(-S,0); % P1(10)
            V_s = abs(ones(1,size(X,2))-sum(S,1)); % P2(11)
            V_r = sqrt(sum((X - A*S).^2,1)); % P3(12)


       %% 4: Generating location-based and superpixel-based saliency maps
        % location-based saliency map
        V = [mat2gray(V_k); mat2gray(V_s); mat2gray(V_r)]; 
        [V, C] = maxNormlizeLocalmax(V, row, col);
        V (:, IR==1) = 0; 
        C (IR==1) = 0;
        % object-based saliency map
        Dp = zeros(1,P); D = zeros(size(C));
        for i = 1 : N
            Dp(l(i)) = Dp(l(i)) + C(i);  
        end
        Dp = Dp./nk;
        for i = 1 : N
            D(i) = Dp(l(i));
        end
        % find the winner superpxiel
        D_win = max(Dp);
        l_win = find(Dp == D_win);
        nk_win = nk(l_win);  


       %% 5: Stop conditions
       % #1  c_win/d_win
        C_nk = sort(C,'descend'); 
        C_win = sum(C_nk(1:nk_win))/nk_win;
        t1(k) = C_win/D_win;

       % #2  tol
        t2(k) = norm(V_r(l == l_win),'fro')/sqrt(nk_win*L)/max(X(:)); 
        fprintf('%6g\t    %6.4f\t    %6.4f\t    %6.4f\t     %6.4f\t    %4g\t \n', k, t2(k), C_win, D_win, t1(k), nk_win);

        if k > 4
            % #1
            if t1(k) > 2 
                fprintf('\n There is no salient object in scene. The number of endmembers =%2g.\n', k-1);
                if para.verbose 
                    figure('units','normalized','position',[0.7,0.2,0.2,0.3]); bar(t1); hold on;
                    xlabel('iter','fontsize',9); hold on;
                    ylabel('tol-1','fontsize',9); hold on;
                    title(['Stop condition 1, VD = ', num2str(k-1)],'fontsize',9);  hold on;  
                    plot([0,k+2],[2, 2],'r-','LineWidth',1); drawnow;
                    if para.saveimage
                        saveas(gcf,['..\results','\','tol_plot.tif']); 
                    end
                end
                break;
            end
            % #2
            if t2(k) < 1e-3
                fprintf('\n The HSI is well reconstructed by current endmembers. The number of endmembers =%2g.\n', k-1);
                if para.verbose 
                    figure;bar(-log(t2));
                    xlabel('iter','fontsize',9); hold on;
                    ylabel('tol-2','fontsize',9); hold on;
                    title(['Stop condition 2, VD = ', num2str(k-1)],'fontsize',9);  hold on;  
                    plot([0,k+2],[3, 3],'r-','LineWidth',1); drawnow; 
                    if para.saveimage
                        saveas(gcf,['..\results','\','tol_plot.tif']); 
                    end
                end
                break;
            end
        end


       %% 6: Purify endmember by using X_win
       [a_win,X_win]  = purend (X(:,l==l_win)); 
       IR (l==l_win) = 1; 
       A(:,k) = a_win;


       %% Show intermediate results    
        if para.verbose 
            % Plot saliency maps
            C_2D = reshape(C,row,col);
            D_2D = reshape(D,row,col);
            figure('units','normalized','position',[0.25,0.55,0.5,0.35]);
            subplot_tight(1, 2, 1,[.05 .05]); 
            imagesc(C_2D); axis image;axis off;hold on;
            title(['C No. ', num2str(k)],'fontsize',8);
            [r_i,c_i] = find(C_2D == max(C_2D(:))); plot(c_i(1),r_i(1),'*r');
            subplot_tight(1, 2, 2,[.05 .05]);
            imagesc(D_2D); axis image;axis off;hold on;
            title(['D No. ', num2str(k)],'fontsize',8);
            [r_i,c_i] = find(D_2D == max(D_2D(:))); plot(c_i(1),r_i(1),'*r'); 
            drawnow;        
            % Plot spectra
            figure('units','normalized','position',[0.25,0.1,0.5,0.35]);
            subplot_tight(1, 2, 1,[.05 .05]);plot(X_win); axis([0 L 0 1]);       
            title(['Xwin No. ', num2str(k)],'fontsize',8);     
            subplot_tight(1, 2, 2,[.05 .05]);plot(a_win);axis([0 L 0 1]);
            title(['awin No. ', num2str(k)],'fontsize',8);
            drawnow;
        end

        if  para.saveimage
            % Save saliency maps
            Path_str = ['..\results','\',num2str(k)];  mkdir(Path_str);
            K_v = size(V,1);
            V_2D = reshape(V.',row,col,K_v);
            C_2D = reshape(C,row,col);
            D_2D = reshape(D,row,col);
            for i = 1 : K_v
               figure;axis off;axis image; set(gcf,'visible','off');
               imagesc(V_2D(:,:,i),[0,max(V_2D(:))]);drawnow;
               saveas_center(gcf,[Path_str,'\','V',num2str(i),'.tif'], row, col); 
            end
            % Location-based salency map
            figure;axis off;axis image; set(gcf,'visible','off');
            imagesc(C_2D);hold on;[r_i,c_i] = find(C_2D == max(C_2D(:)));  
            plot(c_i(1),r_i(1),'o','Color','red','LineWidth',1.4); drawnow;
            saveas_center(gcf,['..\results','\','C_',num2str(k),'.tif'], row, col);      
            % Object-based salency map
            figure;axis off;axis image; set(gcf,'visible','off');
            imagesc(D_2D);hold on; [r_i,c_i] = find(D_2D == max(D_2D(:))); 
            plot(round(mean(c_i)),round(mean(r_i)),'o','Color','red','LineWidth',1.4); drawnow;
            saveas_center(gcf,['..\results','\','D_',num2str(k),'.tif'], row, col);
        end 
    end   
end

function [a, A]  = purend (X)
    %   find the winner point in winner superpixel;
    [u,~,~] = svds(X,1); 
    %   delete outliers
    SAD = acos(min(max(X'* abs(u)./(sqrt(sum(X.^2))'.*norm(u,2)),0),1));
    [dis, disindex] = sort(SAD,'ascend'); 
    A = X(:,disindex(dis <= mean(SAD))); %% remove outliers
    % update endmember a
    [u,~,~] = svds(A,1); 
    factor = sum(sqrt(sum(A.^2)))/size(A,2);
    a = abs(u)*factor;    
end

