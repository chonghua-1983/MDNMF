%%
A1 = xlsread('data/Gaussian_disease.xlsx');
A2 = xlsread('data/Gaussian_microbe.xlsx');
mic_dis = xlsread('data/microbe-disease.xlsx');
A12 = accumarray(mic_dis,1);

% introduce symptoms-based disease network
dis_symptom = xlsread('data/Symptom-based disease similarity.xlsx');
dis_symptom_sim = accumarray(dis_symptom(:,1:2), dis_symptom(:,3));
A11 = add0_into(dis_symptom_sim, 39);
A11 = (A11 + A11')/2;

A22 = A22 / sum(sum(A22));
clear mic_dis dis_symptom dis_symptom_sim
parameters.k = 10;
parameters.mu = 0.001;
parameters.tol = 1e-4;
parameters.max_iter = 2000;

%% running MDNMF on microbe-disease dataset
% randomly initialize 50 times
bestobj = 10000000;
for i = 1:50
   [H1, H2, S1, S2, residual] = MDNMF(A1, A2, A12, A11, A22, parameters);
   newobj = residual(end,4);
   if newobj < bestobj
        bestobj = newobj;
        bestH1 = H1;
        bestH2 = H2;
        bestS1 = S1;
        bestS2 = S2;
   end
   clear H1 H2 S1 S2
end
H1 = bestH1; H2 = bestH2; S1 = bestS1; S2 = bestS2;

% selecting member for each modual i=1,2,...,k
% H1_norm = zscore(H1);
% H2_norm = zscore(H2);

% setting threadthold tt for selecting member(microbes, disease) of moduals
% by Z-score for each row of H1(H2)
% tt1 = 1.4;
% tt2 = 1.5;
%  for i = 1:parameters.k
%    c1 = find(H1_norm(:,i) > tt1);
%    module1{i,1} = c1'; 
% 
%    c2 = find(H2_norm(:,i) > tt2);
%    module2{i,1} = c2'; 
%   
%    Co_module{i,1} = c1;
%    Co_module{i,2} = c2;
%  end
% clear c1 c2
% selecting module members using 3 sigma rule or mean + 1.5std
module1 = select_modlemember(H1);
module2 = select_modlemember(H2);
Co_module(:,1) = module1;
Co_module(:,2) = module2;

% model selection: k, tt1, tt2
% quantile(H1_norm, 0.25);

for i = 1:parameters.k    % compute confusion matrix of overlapping elements
    for j = 1:parameters.k
        cm1(i,j) = length(intersect(Co_module{i,1}, Co_module{j,1}));
        cm2(i,j) = length(intersect(Co_module{i,2}, Co_module{j,2}));
    end  
end   
clear i j tt

% judege the connection between moudules
H1_norms = sqrt(sum(H1.^2,1)); S1_star = (H1_norms' * H1_norms) .* S1;
S1_norm = zscore(S1_star);
H2_norms = sqrt(sum(H2.^2,1)); S2_star = (H2_norms' * H2_norms) .* S2;
S2_norm = zscore(S2_star);

microbe = readtable('data/microbe.xlsx','ReadVariableName',false);
for i = 1:parameters.k
    microbe_module{i} = microbe(Co_module{i,2},2);
    microbe_module_del{i} = table2array(microbe_module{i});  % delete ''
    
    if isempty(microbe_module_del{i})
        continue;
    else
        xlrange = ['A', num2str(i)];
        xlswrite('case study/microbe_moudles_mdnmf_k10',microbe_module_del{i}',1,xlrange);
    end
end

%% NMF
parameters.k = 10;
parameters.tol = 1e-4;
parameters.max_iter = 500;

bestobj = 1000000;
for i = 1:50
   [W, H, residual] = NMF(A12, parameters);
   newobj = residual(end,1);
   if newobj < bestobj
        bestobj = newobj;
        bestH = H;
        bestW = W;
   end
   clear H W
end
H = bestH'; W = bestW;
% selecting member for each modual i=1,2,...,k
H_norm = zscore(H);
W_norm = zscore(W);
  
tt1 = 1.4;
tt2 = 1.5;
% tt2 = [1 2];
 for i = 1:parameters.k
   c1 = find(H_norm(:,i) > tt1);
   module1{i,1} = c1'; 

   c2 = find(W_norm(:,i) > tt2);
   module2{i,1} = c2'; 
  
   Co_module{i,1} = c1;
   Co_module{i,2} = c2;
 end
clear c1 c2

microbe = readtable('data/microbe.xlsx','ReadVariableName',false);
for i = 1:parameters.k
    microbe_module{i} = microbe(Co_module{i,1},2);
    microbe_module_del{i} = table2array(microbe_module{i});  % delete ''
    
    if isempty(microbe_module_del{i})
        continue;
    else
        xlrange = ['A', num2str(i)];
        xlswrite('result/microbe_moudles_nmf',microbe_module_del{i}',1,xlrange);
    end
end

 %% NetNMF
 bestobj = 1000000;
for i = 1:50
   [H1, H2, S1, S2, residual] = NetNMF(A1, A2, A12, parameters);
   newobj = residual(end,4);
   if newobj < bestobj
        bestobj = newobj;
        bestH1 = H1;
        bestH2 = H2;
        bestS1 = S1;
        bestS2 = S2;
   end
   clear H1 H2 S1 S2
end
H1 = bestH1; H2 = bestH2; S1 = bestS1; S2 = bestS2;

% selecting member for each modual i=1,2,...,k
H1_norm = zscore(H1);
H2_norm = zscore(H2);

% setting threadthold tt for selecting member(microbes, disease) of moduals
tt1 = 1.4;
tt2 = 1.5;
% tt2 = [1 2];
 for i = 1:parameters.k
   c1 = find(H1_norm(:,i) > tt1);
   module1{i,1} = c1'; 

   c2 = find(H2_norm(:,i) > tt2);
   module2{i,1} = c2'; 
  
   Co_module{i,1} = c1;
   Co_module{i,2} = c2;
 end
clear c1 c2

% model selection: k, tt1, tt2
% quantile(H1_norm, 0.25);

for i = 1:parameters.k    % compute confusion matrix of overlapping elements
    for j = 1:parameters.k
        cm1(i,j) = length(intersect(Co_module{i,1}, Co_module{j,1}));
        cm2(i,j) = length(intersect(Co_module{i,2}, Co_module{j,2}));
    end  
end   
clear i j tt

% judege the connection between moudules
H1_norms = sqrt(sum(H1.^2,1)); S1_star = (H1_norms' * H1_norms) .* S1;
S1_norm = zscore(S1_star);
H2_norms = sqrt(sum(H2.^2,1)); S2_star = (H2_norms' * H2_norms) .* S2;
S2_norm = zscore(S2_star);

microbe = readtable('data/microbe.xlsx','ReadVariableName',false);
for i = 1:parameters.k
    microbe_module{i} = microbe(Co_module{i,2},2);
    microbe_module_del{i} = table2array(microbe_module{i});  % delete ''
    
    if isempty(microbe_module_del{i})
        continue;
    else
        xlrange = ['A', num2str(i)];
        xlswrite('result/microbe_moudles_netNMF',microbe_module_del{i}',1,xlrange);
    end
end



