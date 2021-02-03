function [expansion] = avgline(response, chazhi)
%均值插值，输入为行向量和插成原先的几倍，输出为插值后的行向量
%   三次样条插值运行速度过慢，用简易插值法代替
%    response = [1 2 3 4]; chazhi = 4;
    [x, y] = size(response);
    if(x == 1) t = y;
    else t = x;
    end
    expansion = zeros(1, chazhi * (t - 1) + 1);
    for i = 1 : (t - 1)
        expansion(chazhi * (i - 1) + 1) = response(i);
        for j = 2 : chazhi 
        expansion(chazhi * (i - 1) + j) = ((j - 1) * response(i + 1) + (chazhi - j + 1) * response(i)) / chazhi;
        end
    expansion(1, chazhi * (t - 1) + 1) = response(t);
    end
end

