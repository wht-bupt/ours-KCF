function [output] = tenjudge(A)
%十个数，去掉一个最大值和一个最小值，求平均值
    maxx = max(A);
    minn = min(A);
    output = (sum(A) - maxx - minn) / (length(A) - 2);
end

