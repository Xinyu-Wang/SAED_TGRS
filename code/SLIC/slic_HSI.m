% Preprocessing of SGSNMF:
% Simple Linear Iterative Clustering SuperPixels (SLIC) for HSIs
%
% Usage:   [l, Am, S ,C, Cj] = slic_HSI(im, k, m, seRadius, nItr)
%==========================================================================
% Input:      im - 3D hyperspectral image to be segmented.
%              P - the number of desired superpixels.
%              m - Weighting factor between spectral and spatial differences. {0.5}
%       seRadius - Regions morphologically smaller than this are merged with adjacent regions. {1}
%           nItr - the number of iterations
%
% Output:      l - Labeled image of superpixels. Labels range from 1 to k.
%             Am - Adjacency matrix of segments.  Am(i, j) indicates whether
%                  segments labeled i and j are connected/adjacent
%              S - the averange size of superpixels
%              C - cluster center of each superpixel
%             Cj - the confidence index
%==========================================================================
% 
% Copyright (C) 2018, Xinyu Wang (wangxinyu@whu.edu.cn)
%                     Yanfei Zhong (zhongyanfei@whu.edu.cn)
%                     Wuhan university
%                     All rights reserved.
%
% Reference:
%
% R. Achanta, A. Shaji, K. Smith, A. Lucchi, P. Fua and S. Susstrunk. "SLIC
% Superpixels Compared to State-of-the-Art Superpixel Methods"  PAMI. Vol
% 34 No 11. November 2012. pp 2274-2281. 
%
% X. Wang, Y. Zhong, L. Zhang, and Y. Xu, ¡°Spatial Group Sparsity
% Regularized Nonnegative Matrix Factorization for Hyperspectral Unmixing,¡±
% IEEE Transactions on Geoscience and Remote Sensing, vol. 55, no. 11, pp.
% 6287-6304, NOV 2017, 2017.
%
% X. Wang, Y. Zhong, C. Cui, L. Zhang and Y. Xu, "Autonomous Endmember
% Detection via an Abundance Anomaly Guided Saliency Prior for
% Hyperspectral Imagery," in IEEE Transactions on Geoscience and Remote
% Sensing, vol. 59, no. 3, pp. 2336-2351, March 2021, doi:
% 10.1109/TGRS.2020.3001353..
%
% Last Modified:
% 27 Feb, 2018

function [l, nk ,C] = slic_HSI(im, P, m, seRadius, nItr)

    if ~exist('m','var')   || isempty(m),     m = 0.5;     end
    if ~exist('seRadius','var')   || isempty(seRadius),     seRadius = 1;     end
    if ~exist('nItr','var')   || isempty(nItr),     nItr = 10;     end
    
    [rows, cols, L] = size(im);
    N = rows*cols;
    % DR before SLIC 
    
    if L > 30
        DR = 10;
        im_2D = reshape(im, N, L)';
        [Ud, ~, ~] = svds((im_2D*im_2D.')/N,DR);
        Xd = Ud.'* im_2D;
        im_2D = Xd ./ repmat(sum( Xd .* repmat(mean(Xd, 2),1,N)), DR, 1);
        im = reshape(im_2D',rows,cols, DR);
        L = DR;
    end
%     im_2D = Ud * (Xd ./ repmat(sum( Xd .* repmat(u,1,rows*cols)), 10, 1));
%     im = reshape(im_2D',rows,cols,L);
    
    % averange size of superpxiels
    S = sqrt(N / (P * sqrt(3)/2));
    countCols = round(cols/S - 0.5);
    S = cols/(countCols + 0.5); 
    countRows = round(rows/(sqrt(3)/2*S));
    vSpacing = rows/countRows;
    % Recompute the number of superpixels P
    P = countRows * countCols;
    % Allocate memory and initialise clusters, labels and distances.
    C = zeros(L+3,P);        % Cluster centre data  1:L is mean spectral value,
                             % L+1: L+2 is row, col of centre, 
                             % L+3 is the number of pixels in each clustering
    l = -ones(rows, cols);   % Pixel labels.
    d = inf(rows, cols);     % Pixel distances from cluster centres.
       
    % Initialise clusters on a hexagonal grid
    k = 1;
    r = vSpacing/2;
    for ri = 1:countRows
        if mod(ri,2), c = S/2; else, c = S;  end     
        for ci = 1:countCols
            cc = round(c); rr = round(r);
            C(1:L+2, k) = [squeeze(im(rr,cc,:)); cc; rr];
            c = c + S;
            k = k + 1;
        end
        r = r + vSpacing;
    end
    Ss = ceil(S);
 %  Upadate cluster centers   
    for n = 1:nItr
       for k = 1:P
           % Get subimage around cluster
           rmin = max(C(L+2,k)- Ss, 1);   rmax = min(C(L+2,k)+ Ss, rows); 
           cmin = max(C(L+1,k)- Ss, 1);   cmax = min(C(L+1,k)+ Ss, cols); 
           subim = im(rmin:rmax, cmin:cmax, :);  
           assert(numel(subim) > 0);
           % Calculate distance
           D = dist(C(:, k), subim, rmin, cmin, vSpacing, m);
           % If any pixel distance from the cluster centre is less than its
           % previous value, update its distance and label
           subd =  d(rmin:rmax, cmin:cmax);
           subl =  l(rmin:rmax, cmin:cmax);
           updateMask = D < subd;
           subd(updateMask) = D(updateMask);
           subl(updateMask) = k;
           d(rmin:rmax, cmin:cmax) = subd;
           l(rmin:rmax, cmin:cmax) = subl;           
       end
       
       % Update cluster centres
       C = zeros(L+3,P);
       for r = 1:rows
           for c = 1:cols
              C(:, l(r,c)) = C(:, l(r,c)) + [reshape(im(r,c,:),L,1); c; r; 1];
           end
       end
       for k = 1:P 
           C(1:L+2,k) = C(1:L+2,k)/C(L+3,k); 
       end
       C(L+1:L+2,:) = round(C(L+1:L+2,:));
    end
    
    % Clean up small orphaned regions.     
    if seRadius
        l = mcleanupregions(l, seRadius);
    end
    
    % Recalculate the center
    P = max(unique(l(:)));
    C = zeros(L+3,P); 
    for r = 1:rows
       for c = 1:cols
          C(:, l(r,c)) = C(:, l(r,c)) + [reshape(im(r,c,:),L,1); c; r; 1];
       end
    end
    for k = 1:P
       C(1:L+2,k) = C(1:L+2,k)/C(L+3,k); 
    end
    %Outputs
    nk = C(L+3,:); % 1*N the number of pixels in each superpxiel 
    C = C(1:L,:);  % (L+2) * N clustering center matrix

%     figure;set(gcf,'visible','off');axis image;axis off;
%     imshow(showsegresults(im,l));
%     saveas_center(gcf,'Img_seg.tif', rows, cols);

end
       
    
%%%%%%%%%%%%%%%%%%%%%%  dist  %%%%%%%%%%%%%%%%%%%%%%
%
% Usage:  D = dist(C, im, r1, c1, S, m)
% 
% Arguments:   C - Cluster center
%             im - sub-image surrounding cluster center
%         r1, c1 - row and column of top left corner of sub image within the
%                  overall image.
%              S - grid spacing
%              m - weighting factor between spectral and spatial distances.
%
% Returns:     D - Distance image giving distance of every pixel in the
%                  subimage from the cluster center
%
% Distance = sqrt( dc^2 + (ds/S)^2*m^2 )
% where:
% dc = arccos(ims*imc/sqrt(sum(ims.^2,2)*sum(imc.^2)))  % Spectral distance
% ds = sqrt(dx^2 + dy^2)                                % Spatial distance

function D = dist(C, im, r1, c1, S, m) 
    [rows, cols, L] = size(im);  
    % Squared colour difference
    ims = reshape(im,rows*cols,L);imc = C(1:L); 
    dcos = (ims * imc)./(sqrt(sum(ims.^2,2)) * norm(imc));
    dc2 = acos(min(max(reshape(dcos,rows,cols),0.0000),1.0000));
    % Squared spatial distance
    [x,y] = meshgrid(c1:(c1+cols-1),r1:(r1+rows-1));
    ds2 = (x-C(L+1)).^2 + (y-C(L+2)).^2;
    D = sqrt(dc2 + m^2 * ds2 ./ S^2 );
end


function Img_seg = showsegresults(Y,labels)
[p,q,L] = size(Y);
Img=zeros(p,q,3); lim = [round(L/4),round(L/2),round(3*L/4)];
for itr= 1:3
    Sp = Y(:,:,lim(itr));
    X_r = sort(reshape(Sp,p*q,1)); 
    a=X_r(round(0.02*p*q));
    b=X_r(round(0.98*p*q));
    for i=1:p
        for j=1:q
            if Sp(i,j) < a
                Img(i,j,itr) = a;
            else
                if Sp(i,j)< b
                Img(i,j,itr) = Sp(i,j);
                else
                Img(i,j,itr) = b;
                end
            end
        end    
    end
    Img(:,:,itr)  =(Sp -a).*(255/(b-a));
end
Img=uint8(Img);
Img_seg = drawregionboundaries(labels, Img, [255 0 0]);
end

