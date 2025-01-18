function [] = affinity_tsne_fun(W,gt,fun_name,data_name,method_name)
%% 相似性矩阵，真实标签，选用的功能，数据集名称，模型名称

if fun_name =="affinity"
    imagesc(W);
    ax = gca;  % 获取当前坐标轴句柄
    ax.FontSize = 24;  % 设置字体大小
    text = strcat(method_name,'_',fun_name,'_',data_name,".eps");
    print(text,'-depsc');
end

if fun_name =="tsne"
    Y2 = tsne(W);
    gscatter(Y2(:,1), Y2(:,2),gt);
    legend('off');
    ax = gca;
    % 设置坐标数字的大小
    ax.FontSize = 24;
    text = strcat(method_name,'_',fun_name,'_',data_name,".eps");
    print(text,'-depsc');
end

end

