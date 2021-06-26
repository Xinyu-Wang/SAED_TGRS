function [V, C] = maxNormlizeLocalmax(m, row, col)
[k,~] = size(m);
H = fspecial('gaussian', [7 7], 1);
c = zeros(row, col);
v = zeros(row, col, k);

for i = 1:k
    
  vi = reshape(m(i,:),row, col);  
  M = imfilter(vi,H,'symmetric');
  % find local maximums
  refData = M(2:end-1,2:end-1);
  localMax = (refData >= M(1:end-2,2:end-1)) & ...
             (refData >= M(3:end,2:end-1)) & ...
             (refData >= M(2:end-1,1:end-2)) & ...
             (refData >= M(2:end-1,3:end)) & ...
             (refData >= 0.05);
  maxData = refData(localMax(:));
  % calculate the mean value and number of local maximums
  
  globmax = max(M(:));
  locmax_sum = sum(maxData) - globmax;
  locmax_num = numel(maxData) - 1;
  locmax_avg = locmax_sum / max(locmax_num,1);
  
  % remove visual redundancy
  if locmax_num > 1
     vi = max(vi - locmax_avg, 0);
  end
  
  if k > 5
      vi = globmax(1) * vi;
  end
  
  % the location-based salency map
  v(:,:,i) = vi;
  c = c + vi;
  
end
V = reshape(v,row*col,k)';
C = reshape(c,1,row*col);
