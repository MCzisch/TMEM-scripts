%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Extract regressors for CompCorr method             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Log: 30.09.2016 MC from IE, - scritp was used for unified IST project
% CompCor component extraction for IST-task /RSN and DEX-RSNs -

%
warning off all
clear all
addpath /usr/local/MATLAB/R2015a/toolbox/nifti/;
basis_scripts = '/media/spectro/_media_spectro/RAID/data/TMEM_RSN/autoscripts/';

cd('/media/spectro/_media_spectro/RAID/data/TMEM_RSN');
dirmaster = dir('RSN_*');


% RSN 1, 10
%for subject = 11:master_end
for subject = 32:size(dirmaster,1)
    
    cd ('/media/spectro/_media_spectro/RAID/data/TMEM_RSN/');
    subject_dir = dirmaster(subject).name;
    
    clear file code maskfile fid 
  
    % begin user area -----------------
    cd (subject_dir);
    c = 'rm npca*'; system(c);
    
    %Convert_NIFTI_to_ANALYZE;
    tmp = dir('drwavolbet0*.nii');
    for i=1:size(tmp,1)
        cmd = ['fslchfiletype ANALYZE ', tmp(i).name]; system(cmd);  
    end
    prefix = 'drwavolbet0';
    from     = 5; % NB = drwavol0004 
    to       = from + size(dir('drwavolbet0*.img'),1)-1; %NB: -1 because of zero-based counting 
    maxc     = 50; % use 4 as default
    
    tmp = dir('rwrc3vol0000.nii');
    if size(tmp,1)
        cmd = ['fslchfiletype ANALYZE rwrc3vol0000.nii']; system(cmd);
    end
    file(1,:) = ['rwrc3vol0000.img     '];
    
    tmp = dir('rwrc2vol0000_woBG.nii');
    if size(tmp,1)
        cmd = ['fslchfiletype ANALYZE rwrc2vol0000_woBG.nii']; system(cmd);
    end
    file(2,:) = ['rwrc2vol0000_woBG.img'];
    
    tmp = dir('std98_perc.nii');
    if size(tmp,1)
        cmd = ['fslchfiletype ANALYZE std98_perc.nii']; system(cmd);
    end
    file(3,:) = ['std98_perc.img       '];
    
    code(1,:) = 'fpca_csf  '; thr(1) = 0.95; % 0.80 % more lenient threshold
    code(2,:) = 'fpca_wm   '; thr(2) = 0.99; % strict threshold as large area!
    code(3,:) = 'fpca_std98'; thr(3) = 0;    % already a binary mask
    
    cmd = ['gunzip betmask.nii.gz']; system(cmd);
    cmd = ['fslchfiletype ANALYZE betmask.nii']; system(cmd);
    maskfile  = ['betmask.img'];
    
    dummyheader = [basis_scripts,'dummy_dartel.hdr'];
    % end user area -------------------
    
    
    % calculate number of images
    n_images = to-from+1;
    %cd (dir_root);
    c1 = clock;
    
    % ---------------------------------------------------------------------
    % read in masks for the 3 voxel pools and betmask, threshhold and erode
    % ---------------------------------------------------------------------
    
    %SQU_csf = 2; SQU_wm = 3; SQU_std98 = 2;
    SQU_csf = 1; SQU_wm = 1; SQU_std98 = 1;
    
    fid = fopen(deblank(file(1,:)),'r'); mask(1,:) = fread(fid,inf,'float32') > thr(1); fclose(fid);
    tmp = reshape(mask(1,:),91,109,91);
    mask_vol(:,:,:,1) = tmp;
    tmp = imerode(tmp, strel('square',SQU_csf));
    emask_vol(:,:,:,1) = tmp;
    emask(1,:) = reshape(tmp,1,902629);
    
    fid = fopen(deblank(file(2,:)),'r'); mask(2,:) = fread(fid,inf,'float32') > thr(2); fclose(fid); tmp = reshape(mask(2,:),91,109,91);
    mask_vol(:,:,:,2) = tmp; tmp = imerode(tmp, strel('square', SQU_wm)); emask_vol(:,:,:,2) = tmp; emask(2,:) = reshape(tmp,1,902629);
    
    fid = fopen(deblank(file(3,:)),'r'); mask(3,:) = fread(fid,inf,'int16') > thr(3); fclose(fid); tmp = reshape(mask(3,:),91,109,91);
    mask_vol(:,:,:,3) = tmp; tmp = imerode(tmp, strel('square', SQU_std98)); emask_vol(:,:,:,3) = tmp; emask(3,:) = reshape(tmp,1,902629);
    
    fid = fopen(maskfile,'r');  betmask = fread(fid,inf,'int16') > 0; fclose(fid);
    for i=1:3; emask(i,:) = emask (i,:) .* betmask'; end; % only mask 1 & 2 should be betmasked >> changed to all 3
    
    % ------------------------------------
    % Write out eroded and uneroded mask
    % ------------------------------------
    mask_vol=single(mask_vol); emask_vol=single(emask_vol);
    img = make_nii(mask_vol(:,:,:,1), [2 2 2]); save_nii(img, deblank(code(1,:)));  img = make_nii(mask_vol(:,:,:,2), [2 2 2]); save_nii(img, deblank(code(2,:)));  img = make_nii(mask_vol(:,:,:,3), [2 2 2]); save_nii(img, deblank(code(3,:)));
    img = make_nii(emask_vol(:,:,:,1), [2 2 2]); save_nii(img, [deblank(code(1,:)),'_eroded']); img = make_nii(emask_vol(:,:,:,2), [2 2 2]); save_nii(img, [deblank(code(2,:)),'_eroded']); img = make_nii(emask_vol(:,:,:,3), [2 2 2]); save_nii(img, [deblank(code(3,:)),'_eroded']);
    
    % ------------------------------------
    % read in images (generally necessary)
    % ------------------------------------
    
    clear alarm; alarm(1:to-from+1) = 0;
    clear images; images(to-from+1,902629) = 0;
    
    for i=from:to
        str = num2str(i);
        if i<10;str=['0',str];end % he adds one '0' if < 10...
        if i<100;str =['0',str];end % and then another '0' when < 100! aaaah. 
        str = [prefix,str,'.img'];
        fid = fopen(str,'r');
        tmp = fread(fid,inf,'float32','l');fclose(fid); % note high endian 'b' low endian 'l'
        alarm(i-from+1) = (min(tmp)<-20000);
        fprintf(['Reading in ', str,'\n']);
        images(i-from+1,:) = tmp;
    end
    
    % check for byteswap accident
    % ---------------------------
    if sum(alarm) > 0
        fprintf(['Warning. Byte swaps detected.\n']);
        stop;
    end
    
    % ----------------
    % start masterloop
    % ----------------
    
    for master = 1:size(file,1) % loop round masks
        
        clear v_mask ve_mask
        
        v_mask = find (mask(master,:));
        s = size(v_mask,2); fprintf(['Size of v_mask is ',num2str(s),' voxels.\n']);
        
        ve_mask = find (emask(master,:));
        se = size(ve_mask,2); fprintf(['Size of ve_mask is ',num2str(se),' voxels.\n']);
        
        % write mean values into an array mean_values
        clear mean_values
        for i=from:to; tmp = images(i-from+1,:);
            mean_values(master,i-from+1) = mean(tmp(ve_mask));
        end
        
        % ----------------------------------------
        % define array for pca & store mean values
        % append demean step
        % ----------------------------------------
        
        fprintf('Build up X-array.\n');
        tmp = size(ve_mask,2); clear X X_demeaned
        for i=1:tmp
            curpos  = ve_mask(i);
            curtime = images(1:n_images,curpos); % all timepoints of one voxel within mask
            X(:,i)  = curtime;
            X_demeaned(:,i) = X(:,i) - mean(curtime);
        end
        
        % --------------------
        %      perform PCA
        % --------------------
        fprintf('Perform PCA on X_demeaned.\n');
        clear COEFF SCORE LATENT TSQUARE
        [COEFF, SCORE, LATENT, TSQUARE, EXPLAINED] = princomp(X_demeaned);
        
        % ---------------------
        % plot PCA time courses
        % ---------------------
        clear mean_corr; 
          mean_corr = (mean_values(master,:)-mean(mean_values(master,:)));
          figure, plot (mean_corr*50,'k','LineWidth',2), hold on
          hold on;
          for i=1:5
              if i == 1; plot (SCORE(:,i),'r','LineWidth',1);end
              if i == 2; plot (SCORE(:,i),'g','LineWidth',1);end
              if i>2;plot (SCORE(:,i),'b'); end
         end
         str = ['MASK ', num2str(master),': '];
         title([str,'black = mean, red = 1st, green = 2nd, blue = higher']);
          
          
        %-------------------------
%       % plot explained variance 
        %-------------------------
        clear exp_in_n_comp; clear total_expl_in_n_comp;
        
        if maxc > size(EXPLAINED,1) 
            maxc = size(EXPLAINED,1) 
        end;
        
        exp_in_n_comp = cumsum(EXPLAINED(1:maxc));
        total_expl_in_n_comp(master) = exp_in_n_comp(maxc);
        %figure, plot (exp_in_n_comp,'k','LineWidth',3), hold on; 
        %title([str,': variance explained in set components']);
        
        % ---------------------------------------
        % identify origins of first three sources and write out .txt file
        % with component time-series
        % ---------------------------------------
        clear comp_matrix; 
        
        for i = 1:maxc
            
            comp_matrix(:,i) = SCORE(:,i); %make matrix with component_timeseries
            exp_var(subject,master,i)= (EXPLAINED(i)); % samples explained variance per component across masks and subjects 
            
        end
        str = [deblank(code(master,:)),'_',num2str(maxc),'_components.txt']; save(str, 'comp_matrix', '-ascii');
        str = [deblank(code(master,:)),'_',num2str(maxc),'_explained_variance.txt']; save(str, 'exp_in_n_comp', '-ascii');
        
            
        for i = 1:9
           tmp=abs(COEFF(:,i)); dummy = images(1,:); dummy(1,:) = 0; dummy(ve_mask) = tmp;
            fid = fopen([deblank(code(master,:)),'_source',num2str(i),'.img'],'w');
            fwrite(fid,dummy, 'float32');
            fclose(fid);
            copyfile (dummyheader, [deblank(code(master,:)),'_source',num2str(i),'.hdr']); 
        end
        
    
        
        % write out mean value of mask
        % ----------------------------
        tmp2 = mean_values(master,:);
        str = [deblank(code(master,:)),'_mean.txt'];
        save (str, 'tmp2', '-ascii')
        
   
        % create and write out docu_cell
        % ------------------------------
        
        docu(master).subject = {deblank(['/media/spectro/_media_spectro/RAID/data/TMEM_RSN/',subject_dir])};
        docu(master).mask = {file(master,:)};
        docu(master).n_comp = {maxc};
        docu(master).thr = {thr(master)};
        docu(master).voxel_in_mask = {s};
        docu(master).voxel_in_eroded_mask = {se};
        save ('pca_docu.mat', 'docu');
       
    end       
     
        
        
        
end


   %-------------------------------------------
%  % Plot explained variance across subjects
   %-------------------------------------------
%{
for master = 1:3
    
    figure; title(['Component ', num2str(code(master,:))]); hold on;
    for subject = [1:16, 18:51, 53:56, 58:75, 77:110, 112:115, 117:134, 136:169, 171, 173:174, 176:177]
        
        tmp = exp_var(subject, master, 1:maxc); tmp = squeeze(tmp);
        tmp3(subject,:) = cumsum(tmp)';
        plot(tmp3(subject,:)); hold on;
        
    end
    hold off;
    
end
%}   

    
    c2 = clock;
    duration = c2-c1

%cd /mnt/sda/IST_FULLRUN/scripts;
%firstlevel_auto_IST_multisession;