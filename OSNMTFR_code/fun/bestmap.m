function [p1, p2] = bestmap(labels1, labels2)
    % labels1 和 labels2 分别是两个聚类结果的标签向量

    n = length(labels1);
    % 创建匹配矩阵，计算标签之间的共现次数
    match_matrix = zeros(max(labels1), max(labels2));
    for i = 1:n
        match_matrix(labels1(i), labels2(i)) = match_matrix(labels1(i), labels2(i)) + 1;
    end

    % 使用匈牙利算法来找到最佳匹配
   
    [p2, ~] = munkres(-match_matrix);
    % 创建从 labels1 到 labels2 的映射 p1
    p1 = zeros(size(labels1));
    for i = 1:n
        p1(i) = p2(labels2(i));
    end
end