function [] = acc_nmi_omega(z_acc,z_nmi,omega1,text)
figure('Position', [200, 200, 650, 350]);
y = [z_acc;z_nmi];
y = y';
Go = bar(y,1);


% 原始矩阵
originalMatrix = omega1;

% 获取矩阵的行数
numRows = size(originalMatrix, 1);

% 初始化一个空的字符串数组来存储结果
resultArray = strings(numRows, 1);

% 循环遍历每一行并将其转换为所需的格式
for i = 1:numRows
    rowStr = sprintf('[%s]', strjoin(string(originalMatrix(i, :)), ','));
    resultArray(i) = rowStr;
end
legend({'ACC','NMI'},'fontsize',14);
set(gca,'Xticklabel',resultArray,'fontsize',14);
xlabel('The value of the weights','fontsize',16);

text = strcat(text,"_omega_sen_AWTN_",".eps");


print(text, '-depsc');

end