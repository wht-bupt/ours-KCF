function [logi] = neighsum(pool,frame)
%��֤ģ�������ģ��������5֡
%   ����5֡����1�����򷵻�0
if(frame <= 5) 
%     if(sum(pool) == 0) logi = 1;
%     else logi = 0;
%     end
    logi = 0;
elseif (sum(pool((frame - 5) : (frame - 1))) ~= 0) logi = 0;
else logi = 1;
end
end

