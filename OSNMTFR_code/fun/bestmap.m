function [p1, p2] = bestmap(labels1, labels2)
    % labels1 �� labels2 �ֱ��������������ı�ǩ����

    n = length(labels1);
    % ����ƥ����󣬼����ǩ֮��Ĺ��ִ���
    match_matrix = zeros(max(labels1), max(labels2));
    for i = 1:n
        match_matrix(labels1(i), labels2(i)) = match_matrix(labels1(i), labels2(i)) + 1;
    end

    % ʹ���������㷨���ҵ����ƥ��
   
    [p2, ~] = munkres(-match_matrix);
    % ������ labels1 �� labels2 ��ӳ�� p1
    p1 = zeros(size(labels1));
    for i = 1:n
        p1(i) = p2(labels2(i));
    end
end