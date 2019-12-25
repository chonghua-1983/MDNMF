function module = select_modlemember(H)
% adapting activethreold for each disease or microbe
% the row is microbe, the column is module id
% activett(m) = mean(m) + sigma(m) or 3sigma rule
% H: the coefficent matrix obtained by NMF
[n, k] = size(H);
% max_H = max(H, [], 2);

% delete smaller element to obtain variance
% H((H < 0.5 * repmat(max_H,1,k))) = 0;
% b = (H ~= 0);
% num = sum(b,2);
% H_mean = sum(H,2) ./ num;

H_mean = sum(H,2) / k;
for i = 1:n
%     hh = H(i,find(H(i,:) ~= 0));
%     len = length(hh);
    H_std(i) = sqrt(sum((H(i,:) - H_mean(i)) .^2) / (k-1));   
end

module = cell(k,1);
H_std = H_std';
% F = 1 ./ (1 + H_std.^2);
% tt = H_mean .* F + (1-F) .* (H_mean + 3 * H_std); 
% for i = 1:k
%     c1 = find(H(:,i) ./ F >= 1);
%     module{i,1} = c1; 
% end
tt = H_mean + 1.5 * H_std;
 for i = 1:k
     c1 = find(H(:,i) ./ tt >= 1);
     module{i,1} = c1; 
 end



