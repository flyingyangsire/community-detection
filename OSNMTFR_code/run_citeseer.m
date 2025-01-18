clear all;
tic;
addpath('fun/');
addpath('datasets/');
test_num = 1;
load('citeseer.mat');
%lambda1=[0.001,0.01,0.1,1,10,100,1000];
%lambda2=[0.001,0.01,0.1,1,10,100,1000];
lambda1=0.001;
lambda2=0.001;
lambda3=0.001;
omega = [5,30];  %1*V
K1 = 20;
%%转换成邻接矩阵
%假设x是一个NxM的矩阵，其中N是点的数量，D是维度
X1=W_cell{1};
X2=W_cell{2};

%计算距离矩阵
D1 = pdist2(X1,X1);
D2 = pdist2(X2,X2);

%构建邻接矩阵
A1 = zeros(size(X1,1));
A2 = zeros(size(X2,1));


for i=1:size(X1,1)
[sortedDistances, sortedIndices]= sort(D1(i,:));
A1(i,sortedIndices(2:K1+1))=1; %从2开始，因为最近的是自己
end

for i=1:size(X2,1)
[sortedDistances, sortedIndices]= sort(D2(i,:));
A2(i,sortedIndices(2:K1+1))=1; %从2开始，因为最近的是自己
end

A{1}=A1;
A{2}=A2;


%for ii=1:7
    %for jj=1:7
for test = 1:test_num
 [Y1, U, Z, C] = AWTN(A,lambda1,lambda2,lambda3,omega,labels);
  Y1 = litekmeans(Y1, C, 'Replicates', 20);
gt = sort(labels);
%[gt,Y1] = bestmap(gt,Y1);
%[ACC,NMI,PUR] = ClusteringMeasure(gt,Y1); %ACC NMI Purity
%[Fscore,Precision,R] = compute_f(gt,Y1);
%[AR,~,~,~]=RandIndex(gt,Y1);
%result = [ACC NMI PUR Fscore Precision R];
%result = [ACC NMI PUR];
%result1 = result;
result(test,:) = EvaluationMetrics(gt, Y1);
    %result(test,:) = result1;
    %ACC_repeat = result1(1);
    %ACC_mea = mean(ACC_repeat);
    %z_acc(ii,jj) = ACC_mea;
    %fprintf("  AC = %5.4f，NMI = %5.4f，purity = %5.4f\n lambda_1:%5.4f, lambda_2:%5.4f, lambda_3:%5.4f",result1(1),result1(2),result1(3),lambda1(ii),lambda2(jj),lambda3(ll));

end
%end
%end

%x_label =string(lambda1);%将数值数组转换为字符串数组
%y_label = string(lambda2);%将数值数组转换为字符串数组兴
%plot_surf(x_label,y_label,'\lambda','\beta',z_acc,"acc","acc NMTF_yaleA.eps");
      
me = mean(result,1);
st = std(result,1);
result(test_num+1,:) = me;
result(test_num+2,:) = st;
   

record_time = toc;
mea_time = record_time/test_num;
