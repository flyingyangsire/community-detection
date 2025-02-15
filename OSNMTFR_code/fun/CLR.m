function [y, S] = CLR(alpha, A00, c, lambda,S0)
% This function is a modification of the code provided by the following work.
% Ref:
% Feiping Nie, Xiaoqian Wang, Michael I. Jordan, Heng Huang.
% The Constrained Laplacian Rank Algorithm for Graph-Based Clustering.
% The 30th Conference on Artificial Intelligence (\textbf{AAAI}), Phoenix, USA, 2016.

% S是一致性图。A是相似性矩阵
viewnum = size(alpha,2);
if nargin<5
   S0 = zeros(size(A00{1},1));
   for v = 1:viewnum
       S0 = S0+alpha(1,v)*A00{v};%根据A初始化S
   end   
end

S0 = S0-diag(diag(S0));
num = size(S0,1); %样本数目
S10 = (S0+S0')/2;
D10 = diag(sum(S10));
L0 = D10 - S10; %初始化拉普拉斯矩阵

NITER = 30;%迭代次数
zr = 10e-11;%收敛阈值

[F0, ~, evs]=eig1(L0, num, 0);  % 求拉普拉斯矩阵的所有特征向量

F = F0(:,2:c+1);  %对F0取前2到c+1列得到F

for iter = 1:NITER   %开始迭代
    dist = L2_distance_1(F',F');  %对F求二范数
    S = zeros(num);  %negative part is all 0
    for i=1:num   %for every i  %遍历样本
        a0 = zeros(1,num);  %a0为1xN的行向量
        for v = 1:viewnum  %遍历视图
            temp = A00{v};
            a0 = a0+alpha(1,v)*temp(i,:);% 计算a_i^{v}
        end
        
        idxa0 = find(a0>0);
        ai = a0(idxa0);  %找出大于0的索引
        di = dist(i,idxa0);  %计算v_i

        ad = (ai-0.5*lambda*di)/sum(alpha);  %subtractor
        S(i,idxa0) = EProjSimplex_new(ad);
   end;
   
   if lambda ==0
       y = zeros(size(A00{1},1),1);
       A = S;
       break;
   end
       
    A = S;
    A = (A+A')/2;
    D = diag(sum(A));
    L = D-A;
    F_old = F; % store F temporaly
    [F, ~, ev]=eig1(L, c, 0);
    evs(:,iter+1) = ev;
    
    fn1 = sum(ev(1:c));
    fn2 = sum(ev(1:c+1));
    if fn1 > zr
        lambda = 2*lambda;
    elseif fn2 < zr
        lambda = lambda/2;
         F = F_old;
    else
        break;
    end;
end;
 
[clusternum, y]=graphconncomp(sparse(A)); y = y';
if clusternum ~= c
    %sprintf('Can not find the correct cluster number: %d', c)
end;


