function [output] = tenjudge(A)
%ʮ������ȥ��һ�����ֵ��һ����Сֵ����ƽ��ֵ
    maxx = max(A);
    minn = min(A);
    output = (sum(A) - maxx - minn) / (length(A) - 2);
end

