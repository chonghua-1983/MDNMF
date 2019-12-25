function [H1, H2, S1, S2, residual] = MDNMF(A1, A2, A12, A11, A22, parameters)
% ||A1 - H1S1H1'|| + lamda1*||A12 - H1H2'|| + lamda2 * ||A2 - H2S2H2'||
k = parameters.k;
mu = parameters.mu;
tol = parameters.tol;
max_iter = parameters.max_iter;

[nr1, nc1] = size(A1);
[nr2, nc2] = size(A2);
[nr, nc] = size(A12);
D1 = diag(sum(A11));
L1 = D1 - A11;
D2 = diag(sum(A22));
L2 = D2 - A22;

if (nr1 ~= nc1) || (nr2 ~= nc2) || (nr1 ~= nr) || (nr2 ~= nc) 
   error('The input data matrices can not be applied to the algorithm!')
end

H1 = abs(randn(nr1,k)); 
H2 = abs(randn(nr2,k));
S1 = abs(randn(k,k)); S1 = (S1 + S1') / 2;
S1 = sum(S1) - diag(S1); S1 = S1 / max(max(S1));
S2 = abs(randn(k,k)); S2 = (S2 + S2') / 2;
S2 = sum(S2) - diag(S2); S2 = S2 / max(max(S2));

lamda1 = nr1 / nr2;
lamda2 = (nr1 / nr2).^2;

iter = 1;
residual = zeros(max_iter, 4);

res = Inf;
ForRes = Inf;

while ((res > tol) && (iter < max_iter))
   % step1:fix H1 H2, update S1 S2
   numerator = H1' * A1 * H1;
   denominator = (H1' * H1) * S1 * (H1' * H1);
   S1 = S1 .* (numerator ./ max(denominator,1e-10));
   clear denominator numerator
   
   H2A2 = H2' * A2;
   numerator =  H2A2 * H2;
   HTH2 = H2' * H2; 
   HTHS = HTH2 * S2;
   denominator = HTHS * HTH2;
   S2 = S2 .* (numerator ./ max(denominator,1e-10));
   clear denominator numerator
   
   % step 2a: Fix H2, S1, S2, update H1
   numerator = 2 * A1 * (H1 * S1) + lamda1 * (A12 * H2) + mu * D1 * H1;
   denominator = 2 * H1 * (S1 * (H1' * H1) * S1) + lamda1 * H1 * (H2' * H2) + mu * A11 * H1 + eps;
   H1 = H1 .* (numerator ./ denominator);
   clear denominator numerator
   
   % step 2b: Fix H1, S1, S2, update H2
   numerator = 2 * lamda2 * A2 * (H2 * S2) + lamda1 * (A12' * H1) + mu * D2 * H2;
   denominator = 2 * lamda2 * H2 * (S2 * (H2' * H2) * S2) + lamda1 * H2 * (H1' * H1)  + eps + mu * A22 * H2;
   H2 = H2 .* (numerator ./ denominator);
   clear denominator numerator
   
   residual(iter,1) = norm((A1 - H1 * S1 * H1'), 'fro') .^2;
   residual(iter,2) = lamda1 * norm((A12 - H1 * H2'), 'fro') .^2;
   residual(iter,3) = lamda2 * norm((A2 - H2 * S2 * H2'), 'fro') .^2;
   residual(iter,4) = residual(iter,1) + residual(iter,2) + residual(iter,3);
   
   res = abs(residual(iter,4) - ForRes);
   ForRes = residual(iter,4);
   iter = iter + 1;
   
    val = max(max(S1));
    sqrtval = sqrt(val);
    S1 = S1 ./ val;
    H1 = H1 .* sqrtval;
    H2 = H2 ./ sqrtval;
    S2 = S2 .* val;
    clear val sqrtval

end
final_iter = iter - 1;
residual = residual(1:final_iter,:);

