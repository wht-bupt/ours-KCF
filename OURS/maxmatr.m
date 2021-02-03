function [C] = maxmatr(A, B)
%
if(max(max(A)) > max(max(B)))
    C = A;
else
    C = B;
end
end

