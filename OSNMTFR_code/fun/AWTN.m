function [Y1, U, Z, C]=AWTN(X,lambda1,lambda2,lambda3,omega,gt)

epson = 1e-7;
C = size(unique(gt),1); % number of clusters
V = size(X,2); % number of views
N = size(X{1},2);% number of data points
A = X;
%normalized X
%for i=1:V
    %X{i} = X{i}./repmat(sqrt(sum(X{i}.^2,1)),size(X{i},1),1);
%end
[W_star, W1_star, L_star] = get_consensus_graph(X, V);
%Initilize Z,E,tensor G,multiplier Y,W
for i = 1:V
    Z{i} = constructW_PKN(X{i},10);
    E{i} = zeros(size(X{i},1),N);
    Y{i} = zeros(size(X{i},1),N);
    G{i} = zeros(N,N);
    W{i} = zeros(N,N);  %multiplier
    U{i} = zeros(N,N);
end

weight = 1/V*ones(1,V);  % 初始化自适应权重
sX = [N, N, V];
pho = 0.1;
mu = 0.1;
pho_max = 10e10;
mu_max = 10e10;
eta1 = 2;
eta2 = 2;
converge_Z=[]; converge_Z_G=[];
I = eye(N);

Isconverg = 0;
iter = 1;
%%  iteration
while(Isconverg == 0)
    % == update U{i} ==
    for i=1:V
        B3 = E{i} - Y{i}./mu;
        U{i} = inv(2*I+2*mu-2*A{i}-2*A{i}'+2*lambda3*L_star-mu*Z{i}-mu*Z{i}'+mu*Z{i}*Z{i}')*(mu*B3-mu*Z{i}*B3);
    end

    % == update E{i} ==
    for i=1:V
        F1 = U{i}-U{i}*Z{i}+Y{i}./mu;
        E{i} = solve_l1l2(F1,lambda2./mu);
    end


    % == update Z{i} ==
    for i = 1:V
        %B1 = U{i}-E{i}+Y{i}./mu;
        %B2 = G{i}-W{i}./pho;
        Z{i} = inv(pho*I+mu*I)*(pho*G{i}+mu*I+U{i}'*Y{i}-mu*U{i}'*E{i}-W{i});
    end

    % == update G{i} ==
    Z_tensor = cat(3, Z{:,:});
    W_tensor = cat(3, W{:,:});
    z = Z_tensor(:);
    w = W_tensor(:);
    [g, ~] = wshrinkObj(z+1/pho*w,lambda1/pho,sX,0,3,omega);
    G_tensor = reshape(g, sX);
    for i=1:V
        G{i} = G_tensor(:,:,i);
    end

    % == update W{i} ==
    for i=1:V
        W{i} = W{i}+pho*(Z{i}-G{i});
    end

    % == update Y{i} ==
    for i=1:V
        Y{i} = Y{i}+mu*(U{i}-U{i}*Z{i}-E{i});
    end

    max_Z=0;
    max_Z_G=0;
    Isconverg = 1;
    for k = 1:V

        if (norm(X{k} - X{k} * Z{k} - E{k}, inf) > epson)
            history.norm_Z = norm(X{k} - X{k} * Z{k} - E{k}, inf);
            Isconverg = 0;
            max_Z = max(max_Z, history.norm_Z);

        end

        if (norm(Z{k} - G{k}, inf) > epson)
            history.norm_Z_G = norm(Z{k} - G{k}, inf);
            Isconverg = 0;
            max_Z_G = max(max_Z_G, history.norm_Z_G);
        end

    end
    converge_Z=[converge_Z max_Z];
    converge_Z_G=[converge_Z_G max_Z_G];


    % == update pho  mu ==
    pho = min(pho_max,pho*eta1);
    mu  = min(mu_max,mu*eta2);
    iter = iter + 1;
    if (iter==50)
        Isconverg = 1;
    end

end
%affinity = (abs(A)+abs(A'))/2;
plot(converge_Z,'LineWidth', 2, 'Color', 'red')
affinity1 = zeros(N,N);
for i = 1:V
    affinity1 = abs(Z{i})+abs(Z{i}');
end

affinity1 = affinity1/V;

Y1 = affinity1;
[~, Y1] = max(Y1, [], 2);

