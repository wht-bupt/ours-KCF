function [logi] = neighsum(pool,frame)
%保证模板池中两模板间隔超过5帧
%   超过5帧返回1，否则返回0
if(frame <= 5) 
%     if(sum(pool) == 0) logi = 1;
%     else logi = 0;
%     end
    logi = 0;
elseif (sum(pool((frame - 5) : (frame - 1))) ~= 0) logi = 0;
else logi = 1;
end
end

