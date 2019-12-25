function [W, H, residual] = NMF(V, parameters)
% V = W * H
k = parameters.k;
tol = parameters.tol;
max_iter = parameters.max_iter;

[nr, nc] = size(V);

W = abs(randn(nr,k)); 
H = abs(randn(k,nc));

iter = 1;
residual = zeros(max_iter,1);

res = Inf;
ForRes = Inf;

while ((res > tol) && (iter < max_iter))
   % step1:fix W, update H
   numerator = W' * V;
   denominator = W' * W * H + eps;
   H = H .* (numerator ./ denominator);
   clear denominator numerator
    
   % step2: Fix H, update W
   numerator = V * H';
   denominator = W * H * H'  + eps;
   W = W .* (numerator ./ denominator);
   clear denominator numerator
     
   residual(iter,1) = norm((V - W * H), 'fro') .^2;
   res = abs(residual(iter,1) - ForRes);
   ForRes = residual(iter,1);
   iter = iter + 1;
   
    val = max(max(W));
    W = W / val;
    H = H * val;
    clear val

end
final_iter = iter - 1;
residual = residual(1:final_iter,1);

