function [H, S, residual] = TriNMF(A, parameters)
% A = HSH'
k = parameters.k;
tol = parameters.tol;
max_iter = parameters.max_iter;

[nr, ~] = size(A);

H = abs(randn(nr,k));
S = abs(randn(k,k)); 
S = (S + S') / 2;
S = sum(S) - diag(S); 
S = S / max(max(S));

iter = 1;
residual = zeros(max_iter, 1);

res = Inf;
ForRes = Inf;

while ((res > tol) && (iter < max_iter))
   % step1:fix H, update S
   numerator = H' * A * H;
   denominator = H' * (H * S * H') * H + eps;
   S = S .* (numerator ./ denominator);
   clear denominator numerator
   
   % step1:fix S, update H
   numerator = A * (H * S);
   denominator = H * (S * (H' * H) * S) + eps;
   H = H .* (numerator ./ denominator);
   clear denominator numerator
    
   residual(iter,1) = norm((A - H * S * H'), 'fro') .^2;
   
   res = abs(residual(iter,1) - ForRes);
   ForRes = residual(iter,1);
   iter = iter + 1;
   
    val = max(max(S));
    sqrtval = sqrt(val);
    S = S / val;
    H = H * sqrtval;
    clear val sqrtval

end
final_iter = iter - 1;
residual = residual(1:final_iter,1);

