function [bg_stable] = bg_variety(im_series)
%计算五幅图片的相似程度
%   计算方法：为不影响快速性（每五帧计算一次），令相邻图片循环移位后相减，找到最小的结果，计算五个结果平均值
%     im1 = rand(48, 48);                           %测试用
%     im2(2 : 48, 2 : 48) = im1(1 : 47, 1 : 47);
%     im2(1, 1 : 48) = zeros(1, 48);
%     im2(1 : 48, 1) = zeros(48, 1);
%     im3 = rand(48, 48);
%     im4 = im1;
%     im4(8 : 10, 8 : 10) = ones(3, 3) + im1(8 : 10, 8 : 10);
%     im5 = im2;
%     im5(8 : 10, 8 : 10) = ones(3, 3) + im2(8 : 10, 8 : 10);
%     im_series = {im1, im2, im3, im4, im5};
    for t = 1 : 4
        img1 = im_series{t};
        img2 = im_series{t + 1};
        img12=imresize(img1,[32,32]);
        img22=imresize(img2,[32,32]);
        imgdct1=dct2(img12);    %计算二维dct
        imgdct2=dct2(img22);
        imgdct11=imgdct1(1:8,1:8);  %截取左上角8*8
        imgdct21=imgdct2(1:8,1:8);
        mean1=sum(imgdct11(:))/64;  %计算均值
        mean2=sum(imgdct21(:))/64;
        imghash1=zeros(8,8);
        imghash2=zeros(8,8);
        for i=1:8   %遍历生成hash指纹
            for j=1:8
                if(imgdct11(i,j)>=mean1)
                    imghash1(i,j)=1;end
                if(imgdct21(i,j)>=mean2)
                    imghash2(i,j)=1;end
            end
        end
        cyjz=xor(imghash1,imghash2);    %求异或
        hanming=sum(cyjz(:));   %求汉明距离
        diff(t)=(64-hanming)/64;
    end
    bg_stable = sum(diff);
%     [x, y] = size(im_series{1});
%     move_x = floor(x / 40);
%     move_y = floor(y / 40);
%     up_x = x - move_x;
%     up_y = y - move_y;
%     down_x = move_x + 1;
%     down_y = move_y + 1;
%     size_x = x - 2 * move_x;
%     size_y = y - 2 * move_y;
%     diff = 1000 * ones(1, 4);
%     for t = 1 : 4
%         for i = 1 : 2 * move_x + 1
%             for j = 1 : 2 * move_y + 1
%                 cha = sum(sum(abs(im_series{t}(down_x : up_x, down_y : up_y)  - im_series{t + 1}(i : size_x + i - 1, j : size_y + j - 1)))) / 255;
%                 cha = floor(cha) / (size_x * size_y);
%                 if(cha < diff(t)) diff(t) = cha; end
%             end
%         end
%     end
%     bg_stable = sum(diff);

end

