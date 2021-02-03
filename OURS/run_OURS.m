  function results = run_OURS(seq, res_path, bSaveImage)

    kernel_type = 'gaussian'; 
    kernel.type = kernel_type;


    padding = 1.5;  %extra area surrounding the target
    lambda = 1e-4;  %regularization
    output_sigma_factor = 0.1;  %spatial bandwidth (proportional to target)

    interp_factor = 0.02;
    kernel.sigma = 0.5; 
    kernel.poly_a = 1;
    kernel.poly_b = 9;  
    features.hog = true;
    features.gray = false;
    features.hog_orientations = 9;
    cell_size = 4;  
    show_visualization = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%–ﬁ∏ƒ’‚¿Ô%%%%%%%%%%%%%%%%%%5
%    seq_KCF = subS;
    seq_KCF = seq;
    target_sz = seq_KCF.init_rect(1,[4,3]);
    pos = seq_KCF.init_rect(1,[2,1]) + floor(target_sz/2);
    img_files = seq_KCF.s_frames;
    video_path = [];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [positions , time] = tracker(video_path, img_files, pos, target_sz, ...
            padding, kernel, lambda, output_sigma_factor, interp_factor, ...
            cell_size, features,show_visualization);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if bSaveImage
        imwrite(frame2im(getframe(gcf)),[res_path num2str(frame) '.jpg']); 
    end

    %return results to benchmark, in a workspace variable
    rects = [positions(:,2) - target_sz(2)/2, positions(:,1) - target_sz(1)/2];
    rects(:,3) = target_sz(2);
    rects(:,4) = target_sz(1);

    fps = numel(img_files) / time;
    results.type = 'rect';
    results.res = rects;%each row is a rectangle
    results.fps = fps;

  end