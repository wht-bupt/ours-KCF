function [output] = qlt(response)
%平均能量判据，来自论文《基于改进KCF的室内人体追踪方法》
fmax = max(max(response));
fmin = min(min(response));
[x, y] = size(response);
sumdis = 0;
for i = 1 : x
    for j = 1 : y
        sumdis = sumdis + (response(i, j) - fmin) ^ 2;
    end
end
output = x * y * (fmax - fmin) ^ 2 / sumdis;
end

