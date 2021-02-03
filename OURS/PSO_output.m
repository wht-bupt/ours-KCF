function [distavg] = PSO_output(response1)
%PSO算法（粒子群算法）判定原图是否发生遮挡
%   response1为重新拼接后的响应矩阵。本函数通过粒子群算法，判定粒子群中各粒子到种群最优解
%   的距离并对其归一化，并由本函数来决定是否发生了异常。
%   值得注意的是，粒子群算法的目的并不是要优化的很好来找到异常情况下的最优解，而是要找到一个对
%   异常情况抵抗力很差的粒子群，来区分出异常与正常情况。
%   粒子群迭代完成后，依据种群最优解是否陷入了局部极大值做一次强制判定，若陷入局部极大值则直接
%   将distavg置为-1，此时说明原图状态为峰值接近的多峰，需要增强追踪；若未陷入局部极大值，但
%   distavg仍然很大，此时原图状态为峰值很低的单峰，说明目标发生形变/光照变化等，可以作为一个
%   提供异常判断的依据。
%    clear; clc;
%    response1 = [2 6 5 2 5; 3 3 4 4 5; 0 1 8 -1 3; 4 2 5 6 6; 2 7 3 4 4];
%    response1 = rand(12,35);
    [x, y] = size(response1);
    rand('seed', sum(100*clock));
    
    %设定粒子群初值
    seednum = 10;
    pso_x(1 : seednum) = randi([floor(0.2 * x), ceil(0.8 * x)], [1, seednum]);
    pso_y(1 : seednum) = randi([floor(0.2 * y), ceil(0.8 * y)], [1, seednum]);
%      pso_x(1 : 15) = randi([1, floor(0.4*x)], [1, 15]);
%      pso_x(16 : 30) = randi([ceil(0.6*x), x], [1, 15]);
%      pso_y(1 : 7) = randi([1, floor(0.4*y)], [1, 7]);
%      pso_y(8 : 15) = randi([ceil(0.6*y), y], [1, 8]);
%      pso_y(16 : 22) = randi([1, floor(0.4*y)], [1, 7]);
%      pso_y(23 : 30) = randi([ceil(0.6*y), y], [1, 8]);
     pso_x(1) = ceil(x / 2); pso_y(1) = ceil(y / 2); 
%      pso_x(2) = ceil(x / 2) + 1; pso_y(1) = ceil(y / 2) + 1; 
%      pso_x(3) = ceil(x / 2) + 1; pso_y(8) = ceil(y / 2) - 1;
%      pso_x(4) = ceil(x / 2) - 1; pso_y(16) = ceil(y / 2) - 1;
%      pso_x(5) = ceil(x / 2) - 1; pso_y(23) = ceil(y / 2) + 1;
%     pso_x(2) = ceil(3*x / 4); pso_y(2) = ceil(3*y / 4);
%     pso_x(3) = ceil(x / 4); pso_y(3) = ceil(y / 4);
%     pso_x(4) = ceil(3*x / 4); pso_y(4) = ceil(y / 4);
%     pso_x(5) = ceil(x / 4); pso_y(5) = ceil(3*y / 4);
    
    c = 2; w = 0.8; c2 = 2; max0 = 1;
    gen = 90;
      
    for j = 0 : gen
        for i = 1 : seednum
            pso_result(i) = response1(pso_x(i), pso_y(i));
            
        end
        pso_maxp = find(pso_result == max(pso_result(:)), 1);
        pso_maxx = pso_x(pso_maxp);
        pso_maxy = pso_y(pso_maxp);
        pso_postx = pso_x;
        pso_posty = pso_y;
        for i = 1 : seednum
            if(j == 0)
                if(floor(c * (pso_maxx - pso_x(i))) <= max0)
                pso_x(i) = pso_x(i) + floor(c * (pso_maxx - pso_x(i)));
                else
                    pso_x(i) = pso_x(i) + max0 + 1;
                end
                if(pso_x(i) > x) pso_x(i) = x; end
                if(pso_x(i) < 1) pso_x(i) = 1; end
                if(floor(c * (pso_maxy - pso_y(i))) <= max0)
                pso_y(i) = pso_y(i) + floor(c * (pso_maxy - pso_y(i)));
                else
                    pso_y(i) = pso_y(i) + max0 + 1;
                end
                if(pso_y(i) > y) pso_y(i) = y; end
                if(pso_y(i) < 1) pso_y(i) = 1; end
                pso_histx(i) = pso_x(i);
                pso_histy(i) = pso_y(i);
                vx(i) = pso_x(i) - pso_postx(i);
                vy(i) = pso_y(i) - pso_posty(i);
            else
                aaax = floor(w * vx(i)) + floor(c2 * (pso_histx(i) - pso_x(i))) + floor(c * (pso_maxx - pso_x(i)));
                pso_x(i) = pso_x(i) + min([aaax, max0]);
                if(pso_x(i) > x) pso_x(i) = x; end
                if(pso_x(i) < 1) pso_x(i) = 1; end
                if(pso_x(i) > pso_histx(i)) pso_histx(i) = pso_x(i); end
                aaay = floor(w * vy(i)) + floor(c2 * (pso_histy(i) - pso_y(i))) + floor(c * (pso_maxy - pso_y(i)));
                pso_y(i) = pso_y(i) + min([aaay, max0]);
                if(pso_y(i) > y) pso_y(i) = y; end
                if(pso_y(i) < 1) pso_y(i) = 1; end
                if(pso_y(i) > pso_histy(i)) pso_histy(i) = pso_y(i); end
                vx(i) = pso_x(i) - pso_postx(i);
                vy(i) = pso_y(i) - pso_posty(i);
            end
        end
    end
    for i = 1 : seednum
        pso_result(i) = response1(pso_x(i), pso_y(i));
    end
    distavg = 0;
    pso_maxp = find(pso_result == max(pso_result(:)), 1);
    pso_maxx = pso_x(pso_maxp);
    pso_maxy = pso_y(pso_maxp);
   % biggest = response1(pso_maxx, pso_maxy);
    for i = 1 : seednum
        distavg = distavg + ((pso_maxx - pso_x(i))^2 + (pso_maxy - pso_y(i))^2);
    end
    distavg = 100 * distavg / (seednum * x * y);
    
    neigh = 6;
    
    if(max(max(response1)) > max(pso_result))
        if((pso_maxx - neigh) > 0) && ((pso_maxx + neigh) < x) && (pso_maxy - neigh > 0) && (pso_maxy + neigh < y)
            pso_neighbor(1:(2*neigh+1), 1:(2*neigh+1)) = response1((pso_maxx - neigh) : (pso_maxx + neigh), (pso_maxy - neigh) : (pso_maxy + neigh));
            if(max(max(response1)) > max(max(pso_neighbor)))
                distavg = -1;                 %强制判定是否陷入了局部极大值
            end
        else
            distavg = 1;
        end
    end
end