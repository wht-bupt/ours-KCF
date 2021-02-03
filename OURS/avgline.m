function [expansion] = avgline(response, chazhi)
%��ֵ��ֵ������Ϊ�������Ͳ��ԭ�ȵļ��������Ϊ��ֵ���������
%   ����������ֵ�����ٶȹ������ü��ײ�ֵ������
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

