function [positions, time] = tracker(video_path, img_files, pos, target_sz, ...
	padding, kernel, lambda, output_sigma_factor, interp_factor, cell_size, ...
	features, show_visualization)
%TRACKER Kernelized/Dual Correlation Filter (KCF/DCF) tracking.
%   This function implements the pipeline for tracking with the KCF (by
%   choosing a non-linear kernel) and DCF (by choosing a linear kernel).
%
%   It is meant to be called by the interface function RUN_TRACKER, which
%   sets up the parameters and loads the video information.
%
%   Parameters:
%     VIDEO_PATH is the location of the image files (must end with a slash
%      '/' or '\').
%     IMG_FILES is a cell array of image file names.
%     POS and TARGET_SZ are the initial position and size of the target
%      (both in format [rows, columns]).
%     PADDING is the additional tracked region, for context, relative to 
%      the target size.
%     KERNEL is a struct describing the kernel. The field TYPE must be one
%      of 'gaussian', 'polynomial' or 'linear'. The optional fields SIGMA,
%      POLY_A and POLY_B are the parameters for the Gaussian and Polynomial
%      kernels.
%     OUTPUT_SIGMA_FACTOR is the spatial bandwidth of the regression
%      target, relative to the target size.
%     INTERP_FACTOR is the adaptation rate of the tracker.
%     CELL_SIZE is the number of pixels per cell (must be 1 if using raw
%      pixels).
%     FEATURES is a struct describing the used features (see GET_FEATURES).
%     SHOW_VISUALIZATION will show an interactive video if set to true.
%
%   Outputs:
%    POSITIONS is an Nx2 matrix of target positions over time (in the
%     format [rows, columns]).
%    TIME is the tracker execution time, without video loading/rendering.
%
%   Joao F. Henriques, 2014
%    лл���ߣ�����˵�ú�
%    �������û�ͼ�������£�
%    figure4��������Ⱥ�����figure5����ԭ������ά���棬figure6�����ض�֡���棬figure7���Ʋ�ֵ����ά����
%    figure8����qlt��figure9�鿴��ֵ��figure10�鿴lock, figure11�鿴

    draw4 = 0; draw5 = 0; draw6 = 0; draw7 = 0; draw8 = 0; draw9 = 0; draw10 = 0; draw11 = 0; draw12 = 0;
    frame6 = 42;   %��Ҫ�鿴��֡��
    poolnum = 10;   %ģ��ش�С
    allowc = 0;    %�Ƿ������ֵ
    first_judge = 0.35;
%% ���ߴ����1

	%if the target is large, lower the resolution, we don't need that much
	%detail
	resize_image = (sqrt(prod(target_sz)) >= 100);  %diagonal size >= threshold
	if resize_image
		pos = floor(pos / 2);
		target_sz = floor(target_sz / 2);
    end

	%window size, taking padding into account
	window_sz = floor(target_sz * (1 + padding));
 %   window_sz = 2 .^ nextpow2(window_sz);
    %����KCF��Χ
    lock(1) = 0;
     if(sum(lock) ~= 0)
%         padding = padding * 4;
%        window_sz = 2 .^ nextpow2(window_sz);
%        window_sz = floor(target_sz * (1 + padding));
     end
	
% 	%we could choose a size that is a power of two, for better FFT
% 	%performance. in practice it is slower, due to the larger window size.
% 	window_sz = 2 .^ nextpow2(window_sz);

	
	%create regression labels, gaussian shaped, with a bandwidth
	%proportional to target size
	output_sigma = sqrt(prod(target_sz)) * output_sigma_factor / cell_size;
	yf = fft2(gaussian_shaped_labels(output_sigma, floor(window_sz / cell_size)));

	%store pre-computed cosine window
	cos_window = hann(size(yf,1)) * hann(size(yf,2))';	
	
	
	if show_visualization  %create video interface
		update_visualization = show_video(img_files, video_path, resize_image);
	end
	
	
	%note: variables ending with 'f' are in the Fourier domain.

	time = 0;  %to calculate FPS
	positions = zeros(numel(img_files), 2);  %to calculate precision

	for frame = 1:numel(img_files)
		%load image
		im = imread([video_path img_files{frame}]);
		if size(im,3) > 1
			im = rgb2gray(im);
		end
		if resize_image
			im = imresize(im, 0.5);
		end

		tic()

		if frame > 1
			%obtain a subwindow for detection at the position from last
			%frame, and convert to Fourier domain (its size is unchanged)
			patch = get_subwindow(im, pos, window_sz);
			zf = fft2(get_features(patch, features, cell_size, cos_window));
			
%% �����Ķ�������ڵ�����ʱ���������ڵ����ߵ�����
            %����������ֵĲ��ԣ�
            %ʹ�ñ���ı�׼ģ����е�����������ģ�����ͶƱ����ģ��response�����趨����ֵ��Ϊ�Ͽ��ڵ�������
            %�������λ��Ϊ����ģ��response����λ��
            if(sum(lock) == 0 || sum(pool) == 0)
                %calculate response of the classifier at all shifts
                switch kernel.type
                case 'gaussian'
                    kzf = gaussian_correlation(zf, model_xf, kernel.sigma);
                case 'polynomial'
                    kzf = polynomial_correlation(zf, model_xf, kernel.poly_a, kernel.poly_b);
                case 'linear'
                    kzf = linear_correlation(zf, model_xf);
                end
                response = real(ifft2(model_alphaf .* kzf));  %equation for fast detection
                lock(frame) = 0;
            else
                response = 0;
                vote = 0;
                for i = 1 : length(alphaf_pool)
                    model_xf = xf_pool{i};
                    model_alphaf = alphaf_pool{i};
                    switch kernel.type
                    case 'gaussian'
                        kzf = gaussian_correlation(zf, model_xf, kernel.sigma);
                    case 'polynomial'
                        kzf = polynomial_correlation(zf, model_xf, kernel.poly_a, kernel.poly_b);
                    case 'linear'
                        kzf = linear_correlation(zf, model_xf);
                    end
                    response3 = real(ifft2(model_alphaf .* kzf));  %equation for fast detection
                    if(max(max(response3)) >= lock_threshold) vote = vote + 1;end
                    response = maxmatr(response3, response);
                end
                %sprintf("pool = %d, vote = %d, max = %.4f", length(alphaf_pool), vote, max(max(response)))
%                if(vote >= floor(length(alphaf_pool) / 2))
                if(vote >= 1)
                    lock(frame) = -1;
                else
                    lock(frame) = 0;
                end
            end
            
            
%% �˴���ʼΪ��Ӳ���
            
            %�����仯����
            if(mod(frame, 5) == 1)
                if(frame < 22)
                    im_series{floor(frame / 5) + 1} = im;
                else
                    im_series{mod(floor(frame / 5), 5) + 1} = im;
                end
            end
            if(frame < 25)
                bg_stable(frame) = 0;
            elseif(mod(frame, 25) == 1)
                bg_stable(frame) = bg_variety(im_series);
            else
                bg_stable(frame) = bg_stable(frame - 1);
            end
            
            %��������
            [xnum, ynum] = size(response);
            [xim, yim] = size(im);
            xcmp = xnum / xim;
            ycmp = ynum / yim;
            if(frame >= 22)
                velocity(frame) = 0;
                for i = 1 : 19   %����Manhattan�������
                    velocity(frame) = velocity(frame) + sum(abs(pos_post{i} - pos_post{i + 1}));
                end
                velocity(frame) = 0.001 * velocity(frame) / (19 * xcmp * ycmp);
                %sprintf("bg_stable = %.4f, velocity = %.4f", bg_stable, velocity(frame))
            end
            
            %����֡������
            if(frame >= 22)
                velo2 = 1 / (velocity(frame) + 0.001);
                bg2 = 5 * (bg_stable(frame) - 3);
                upper = floor((velo2 ^ 2) * (1 + bg2));
                upper(upper > 20) = 20; 
                upper(upper < 2) = 0;
            else
                upper = 0;
            end
            
            %ʶ����������
            if(frame >= 26)
                response_threshold = 0.35 + 0.2 * (bg_stable(frame) - 3);
                response_threshold(response_threshold > 0.5) = 0.5;
                response_threshold(response_threshold < 0.35) = 0.35;
                lock_threshold = 0.23 + 0.14 * (bg_stable(frame) - 3);
                lock_threshold(lock_threshold > 0.35) = 0.35;
                lock_threshold(lock_threshold < 0.22) = 0.22;
                pool_threshold = response_threshold;
            else
                response_threshold = 0.35;
                lock_threshold = 0.23;
                pool_threshold = 0.35;
            end
            
            %������ά����
            xx = linspace(1 , xnum, xnum);
            yy = linspace(1 , ynum, ynum);
            xxc = xx;
            yyc = yy;
            
            %ת��ƴ���Ľ�����
            response1 = [];
            a = floor(size(zf,1) / 2);
            c = floor(size(zf,2) / 2);
            response1(1 : a, 1 : c) = response((size(zf,1) - a + 1) : size(zf,1), (size(zf,2) - c + 1) : size(zf,2));
            response1((a + 1) : size(zf,1), 1 : c) = response(1 : (size(zf,1) - a), (size(zf,2) - c + 1) : size(zf,2));
            response1(1 : a, (c + 1) : size(zf,2)) = response((size(zf,1) - a + 1) : size(zf,1), 1 : (size(zf,2) - c));
            response1((a + 1) : size(zf,1), (c + 1) : size(zf,2)) = response(1 : (size(zf,1) - a), 1 : (size(zf,2) - c));

            % ԭ�����������
            if(draw5 == 1)
                zzz = response1(xx, yy);
                if(xnum > ynum)
                    for i = ynum+1 : 1 : xnum
                        yy(i) = i;
                end
                else
                    for i = xnum+1 : 1 : ynum
                        xx(i) = i;
                    end
                end
%                [xxx, yyy] = meshgrid(xx, yy);
                if(xnum > ynum)
                    zzz(:,ynum+1:xnum)=0;
                else
                    zzz(xnum+1:ynum,:)=0;
                end
                figure(5),
                mesh(xx, yy, zzz);
                zlim([-0.1, 0.6]);
            end

            % �ض�֡�鿴
            if(draw6 == 1)
                if(frame == frame6)
                    zzz = response1(xx, yy);
                if(xnum > ynum)
                    for i = ynum+1 : 1 : xnum
                        yy(i) = i;
                end
                else
                    for i = xnum+1 : 1 : ynum
                        xx(i) = i;
                    end
                end
%                [xxx, yyy] = meshgrid(xx, yy);
                if(xnum > ynum)
                    zzz(:,ynum+1:xnum)=0;
                else
                    zzz(xnum+1:ynum,:)=0;
                end
                figure(6),
                mesh(xx, yy, zzz);
                zlim([-0.1, 0.6]);
                end
            end
            
            %��ԭ������в�ֵ
            if(allowc == 1)
                chazhi = 4;
                response2 = [];
                response_temp = [];
                for i = 1 : ynum
                    response_temp(:, i) = avgline(response1(:, i), chazhi);
                end
                for i = 1 : (chazhi * xnum - (chazhi - 1))
                    response2(i, :) = avgline(response_temp(i, :), chazhi);
                end
            else
                response2 = response1;
            end
            
           
            
            %��ֵ�����������
            if(draw7 == 1)
                figure(7);
                [xnum2, ynum2] = size(response2);
                xx2 = 1 : xnum2;
                yy2 = 1 : ynum2;
                zzz = response2(xx2, yy2);
                if(xnum2 > ynum2)
                    for i = ynum2+1 : 1 : xnum2
                        yy2(i) = i;
                    end
                else
                    for i = xnum2+1 : 1 : ynum2
                        xx2(i) = i;
                    end
                end
                if(xnum2 > ynum2)
                    zzz(:,ynum2+1:xnum2)=0;
                else
                    zzz(xnum2+1:ynum2,:)=0;
                end
                mesh(xx2, yy2, zzz);
                zlim([-0.1, 0.6]);
            end

            distavg(frame) = PSO_output(response2);
           % qltt(frame) = qlt(response1);


            % ����Ⱥ���
            if(draw4 == 1)
                figure(4),
                plot(1:frame, distavg);
            end
            
            %qlt
            if(draw8 == 1)
                figure(8),
                plot(1:frame, qltt);
            end
            
            %qlt
            if(draw11 == 1)
                figure(11),
                plot(1:frame, bg_stable);
            end
            
            %qlt
            if(draw12 == 1)
                figure(12),
                plot(1:frame, velocity);
            end
            


            
%% ���ߴ����2
			%target location is at the maximum response. we must take into
			%account the fact that, if the target doesn't move, the peak
			%will appear at the top-left corner, not at the center (this is
			%discussed in the paper). the responses wrap around cyclically.
			[vert_delta, horiz_delta] = find(response == max(response(:)), 1);
%             if((vert_delta == 1) && (horiz_delta == 1)) disp(0); 
%             else disp(vert_delta); disp(horiz_delta); end
			if vert_delta > size(zf,1) / 2  %wrap around to negative half-space of vertical axis
				vert_delta = vert_delta - size(zf,1);
			end
			if horiz_delta > size(zf,2) / 2  %same for horizontal axis
				horiz_delta = horiz_delta - size(zf,2);
			end
			pos = pos + cell_size * [vert_delta - 1, horiz_delta - 1];
		end

		%obtain a subwindow for training at newly estimated target position
		patch = get_subwindow(im, pos, window_sz);
		xf = fft2(get_features(patch, features, cell_size, cos_window));

		%Kernel Ridge Regression, calculate alphas (in Fourier domain)
		switch kernel.type
		case 'gaussian'
			kf = gaussian_correlation(xf, xf, kernel.sigma);
		case 'polynomial'
			kf = polynomial_correlation(xf, xf, kernel.poly_a, kernel.poly_b);
		case 'linear'
			kf = linear_correlation(xf, xf);
		end
		alphaf = yf ./ (kf + lambda);   %equation for fast training
        
%% �����޸Ĳ�����ӣ�ģ�͸��²��Բ���

	    if frame == 1   %first frame, train with a single image
			model_alphaf = alphaf;
			model_xf = xf;
%             alphaf_post = {alphaf};
%             xf_post = {xf};
            pos_post = {pos};
            alphaf_pool = {alphaf};
            xf_pool = {xf};
            lock(1) = 0;
            distavg(1) = 0;
            pool(1) = 1;
            maxi(1) = 0;
            response = 0;
            poslock = 0;
            im_series = {im};
            bg_stable = 0;
            velocity = 0;
            response_threshold = 0.35;
            lock_threshold = 0.23;
            pool_threshold = 0.35;
		else
			%subsequent frames, interpolate model
			model_alphaf = (1 - interp_factor) * model_alphaf + interp_factor * alphaf;
			model_xf = (1 - interp_factor) * model_xf + interp_factor * xf;
        end
        
        %����response��ֵ�����׼ģ��أ����ౣ��poolnum��ģ�壬��ģ��֮�������볬��5֡
        %ģ����������������ģ���豣֤��֮ǰ����ģ����������response���ֵ������һ��ֵ����֤�����ԣ�
        %��һ֡ģ�岻�ᱻ��������ģ����������poolnum�����ӵڶ���ģ�忪ʼ����Ϊ�µ�ģ��
%         if(max(max(response)) >= response_threshold && sum(pool) <= poolnum - 1 && neighsum(pool, frame) == 1)
%             pool(frame) = 1;
%             alphaf_pool{sum(pool)} = model_alphaf;
%             xf_pool{sum(pool)} = model_xf;
%         elseif(frame ~= 1)
%             pool(frame) = 0;
%         end
        if(frame ~= 1) pool(frame) = 0; end
        if(max(max(response)) >= response_threshold && neighsum(pool, frame) == 1)
            for i = 1 : min(sum(pool), poolnum)
                switch kernel.type
                case 'gaussian'
                    kzf_judge = gaussian_correlation(zf, xf_pool{i}, kernel.sigma);
                case 'polynomial'
                    kzf_judge = polynomial_correlation(zf, xf_pool{i}, kernel.poly_a, kernel.poly_b);
                case 'linear'
                    kzf_judge = linear_correlation(zf, xf_pool{i});
                end
                response_judge = real(ifft2(alphaf_pool{i} .* kzf_judge));
                if(max(max(response_judge)) >= pool_threshold)
                    break;
                end
                if(i == sum(pool))
                    pool(frame) = 1;   
                end
            end
        end
        if(pool(frame) == 1)
            if(pool(1) == 0)
                temp = 1 + mod(sum(pool) , 10);             %1��10
            else
                if(sum(pool) <= poolnum)
                    temp = sum(pool);
                else
                    temp = 2 + mod((sum(pool) - 11) , 9);   %2��10
                end        
            end
            alphaf_pool{temp} = model_alphaf;
            xf_pool{temp} = model_xf;
        end
        
        %���ֵ����
        maxi(frame) = max(max(response));
        if(draw9 == 1)
            figure(9),
            plot(1 : frame, maxi);
        end
        
        %���϶��ж��ڵ������ģ�͵��޸ģ��˲��ֽ���֡���Ļ���
        if(frame <= 20)
            pos_post{frame} = pos;
        else
            for i = 1 : 19
                pos_post{i} = pos_post{i + 1};
            end
            pos_post{20} = pos;
        end
        
        if(sum(lock) ~= 0 && poslock == 0 && upper > 0)
            pos = pos_post{20 - upper + 1}; 
            poslock = 1;
        end
        
        %һ֡�ܷ���������
        if(distavg(frame) == -1 && sum(lock) == 0 && frame >= 22)
            lock(frame) = 1;
            poslock = 0;
            %sprintf("velocity = %.4f, bg_stable = %.2f, upper = %d", velocity(frame), bg_stable(frame), upper)
            %sprintf("response_threshold = %.4f, lock_threshold = %.2f", response_threshold, lock_threshold)
        elseif(sum(lock(1 : frame - 1)) == 0)
            lock(frame) = 0;
        end
        
        
        
        %(����)�Ե�һ֡�������Ŷ��жϣ���ǰʮ֡response�����޴󣬷�����һ֡��Ϊ��׼
        if(frame == 11)
            first = tenjudge(maxi(2 : 11));
            %sprintf("first = %.4f", first)
            if(first < first_judge) 
                pool(1) = 0;
            end
        end
        

        
        if(sum(pool) < 2)
            lock(frame) = 0;
        end
%         for i = 1 : frame
%             if(pool(frame) == 1)
%                 disp(frame);
%             end
%         end
        %����lock
        if(draw10 == 1)
            figure(10),
            plot(1:frame, lock);
        end

		%save position and timing
		positions(frame,:) = pos;
		time = time + toc();

		%visualization
		if show_visualization
			box = [pos([2,1]) - target_sz([2,1])/2, target_sz([2,1])];
			stop = update_visualization(frame, box);
			if stop, break, end  %user pressed Esc, stop early
			
			drawnow
% 			pause(0.05)  %uncomment to run slower
		end
		
	end

	if resize_image
		positions = positions * 2;
	end
end

