% ==========================================
% Analysis 2016; TMEM: RSN
% 04.03.2016 MC on KORE
%
% requires IXI templates (I have copied vbm8 from andromeda to "toolbox")
%
% ==========================================

clear all
addpath /usr/local/MATLAB/R2015a/toolbox/nifti/;
addpath /media/spectro/_media_spectro/RAID/data/TMEM_RSN/autoscripts/;
addpath /media/spectro/_media_spectro/RAID/data/TMEM_RSN/autoscripts/filters_Brice/;

cd('/media/spectro/_media_spectro/RAID/data/TMEM_RSN');

% =====================================================
%
%   start of user input section
%
% =====================================================

restart_spm     = 0;
clean_up        = 0;

slice_timing    = 0;
realignment     = 0;
motion_scrubbing= 0;
segmentation    = 0;
dartel_exist    = 0;
dartel_warp     = 0;
dartel_warp_seg = 0;
detrend_4D      = 0;
produce_std98   = 0;
temporal_filter = 0;
BET             = 1;    %check normalized_resolution case
reslice         = 0; 
reslice_wavolbet= 0; 
reslice_segments= 0; 
normalized_resolution = 2;
smoothing       = 0; fwhm_kernel = [6 6 6];
BG_correction   = 0;

% =====================================================
%
%  end of user input section
%
% =====================================================
global basis_scripts;
basis_scripts = '/media/spectro/_media_spectro/RAID/data/TMEM_RSN/autoscripts/';

basis_dir = '/media/spectro/_media_spectro/RAID/data/TMEM_RSN/';

tr= 2.56;
nslices= 40;
ta= tr - tr/nslices;
refslice= 1;
sliceorder = [1:2:40,2:2:40];


% =====================================================

master = dir('RSN_*');

master_start = 1; 
master_end = size(master,1);


% missing 1
% not run 10
for i = 2:master_end
    
    % make sure spm is up
    if restart_spm
        spm('defaults','fMRI');
        spm_jobman('initcfg');
        spm fmri;
    end
    
    fprintf(['### Working on files ', master(i).name, ' ... ###\n\n']);
    basis_fmri = [basis_dir,'/',master(i).name,'/'];
    cd (basis_fmri);
    
    if clean_up
        c = sprintf('rm *avol*  *txt *ps *bet* *c?vol* *master* std*'); system(c);
    end
    
    tmp = dir([basis_fmri,'/vol*.nii']);
    n_images = size(tmp,1);
    
    % ==============
    % define arrays
    % ==============
    
    clear image_array image_array_a image_array_wa image_array_w
    
    % define the first volume name of each stack
    % assuming 5*2048 ms for dummy scans
    
    image_array         = {[deblank(basis_fmri),'vol0005.nii']};
    image_array_a       = {[deblank(basis_fmri),'avol0005.nii']};
    image_array_wa      = {[deblank(basis_fmri),'wavol0005.nii']};
    image_array_wa_bet  = {[deblank(basis_fmri),'wavolbet0005.nii,1']};
    
    for i=6:n_images-1
        
        tmp = sprintf('vol%04d',i);
        tmp2 = sprintf('%04d',i);
        item = [deblank(basis_fmri),tmp,'.nii'];
        item_a      = [deblank(basis_fmri),'a',tmp,'.nii'];
        item_wa     = [deblank(basis_fmri),'wa',tmp,'.nii'];
        item_wa_bet = [deblank(basis_fmri),'wavolbet',tmp2,'.nii'];
        image_array        = [image_array; item];
        image_array_a      = [image_array_a; item_a];
        image_array_wa     = [image_array_wa; item_wa];
        image_array_wa_bet  = [image_array_wa_bet; item_wa_bet];
     end
    image_array_drwa_bet  = {[deblank(basis_fmri),'drwavolbet0005.nii']}; % Need separate loop due to the differences in how the bet images are named.
    image_array_r333kdrwa_bet  = {[deblank(basis_fmri),'r333kdrwavolbet0005.nii']}; % Need separate loop due to the differences in how the bet images are named.
    image_array_kdrwa_bet  = {[deblank(basis_fmri),'kdrwavolbet0005.nii']}; % Need separate loop due to the differences in how the bet images are named.
    image_array_kdr333wa_bet  = {[deblank(basis_fmri),'kdr333wavolbet0005.nii']}; % Need separate loop due to the differences in how the bet images are named.
    for i=6:n_images-1
        tmp = sprintf('bet%04d',i);
        item_drwa_bet = [deblank(basis_fmri),'drwavol',tmp,'.nii'];
        item_dr333wa_bet = [deblank(basis_fmri),'dr333wavol',tmp,'.nii'];
        item_kdrwa_bet = [deblank(basis_fmri),'kdrwavol',tmp,'.nii'];
        item_kdr333wa_bet = [deblank(basis_fmri),'dr333wavol',tmp,'.nii'];
        item_r333kdrwa_bet = [deblank(basis_fmri),'r333kdrwavol',tmp,'.nii'];
        image_array_drwa_bet = [image_array_drwa_bet; item_drwa_bet];
        image_array_kdrwa_bet = [image_array_kdrwa_bet; item_kdrwa_bet];
        image_array_r333kdrwa_bet = [image_array_r333kdrwa_bet; item_r333kdrwa_bet];
    end
    
    % ======================
    % slice time correction
    % ======================
    
    if slice_timing
        clear matlabbatch
        load ([basis_scripts, '/dummy_slicetime_correction.mat']);
        matlabbatch{1}.spm.temporal.st.scans{1} = image_array;
        matlabbatch{1}.spm.temporal.st.nslices = nslices;
        matlabbatch{1}.spm.temporal.st.tr = tr;
        matlabbatch{1}.spm.temporal.st.ta = tr - tr/nslices;
        matlabbatch{1}.spm.temporal.st.so = sliceorder;
        matlabbatch{1}.spm.temporal.st.refslice = refslice;
        inputs = cell(0, 1);
        spm_jobman('run',matlabbatch);
    end
    
    % ============
    % realignment
    % ============
    
    if realignment
        clear matlabbatch
        load ([basis_scripts, '/dummy_realignment.mat']);
        matlabbatch{1}.spm.spatial.realign.estimate.data{1}  = image_array_a; % running on all a-images
        inputs = cell(0, 1);
        spm_jobman('run',matlabbatch);
    end
    
    % ============
    % scrubbing
    % ============
    
    % see do_fsl_scrubbing_RSN.m
    
    % ============
    % segmentation
    % ============
    
    if segmentation % run only once!
        clear matlabbatch
        load ([basis_scripts, '/newsegment_one_case12.mat']);
        matlabbatch{1}.spm.spatial.preproc.channel.vols = {[basis_fmri,'/vol0000.nii,1']};
        matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 2; %just to make sure that rc*nii are written
        matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 1];
        matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 1];
        matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 1];
        spm_jobman ('run', matlabbatch);
    end
    
    % ============
    % DARTEL exist
    % ============
    
    if dartel_exist
        clear matlabbatch
        load ([basis_scripts, '/dummy_dartel_exist.mat']);
        matlabbatch{1}.spm.tools.dartel.warp1.images{1}={[basis_fmri,'/rc1vol0000.nii']};
        matlabbatch{1}.spm.tools.dartel.warp1.images{2}={[basis_fmri,'/rc2vol0000.nii']};
        spm_jobman ('run', matlabbatch);
    end
    
    % ====================
    % DARTEL create warped
    % ====================
    
    if dartel_warp
        clear matlabbatch
        load ([basis_scripts, '/dummy_create_warped.mat']);
        clear tmp
        for i=6:n_images
            tmp{i-5} = [basis_fmri,'u_rc1vol0000.nii'];
        end
        matlabbatch{1}.spm.tools.dartel.crt_warped.flowfields = tmp';
        matlabbatch{1}.spm.tools.dartel.crt_warped.images = {image_array_a};
        spm_jobman ('run', matlabbatch);
        
        % move wa*nii files to destination
        %         c = ['mv ', deblank(Q(master,:)),'/wavol*nii ', basis_fmri];
        %         system(c);
        
    end
    
    if dartel_warp_seg
        clear matlabbatch
        load ([basis_scripts, '/dummy_create_warped.mat']);
        clear tmp
        for i=1:3
            tmp{i} = [basis_fmri,'/u_rc1vol0000.nii'];
        end
        matlabbatch{1}.spm.tools.dartel.crt_warped.flowfields = tmp';
        tmp2 = {[basis_fmri,'/rc1vol0000.nii']};
        tmp2 = [tmp2  ;[basis_fmri,'/rc2vol0000.nii']];
        tmp2 = [tmp2  ;[basis_fmri,'/rc3vol0000.nii']];
        matlabbatch{1}.spm.tools.dartel.crt_warped.images = {tmp2};
        spm_jobman ('run', matlabbatch);
    end
    
    
    if reslice_wavolbet
        clear matlabbatch
        switch normalized_resolution
            case 2
                load ([basis_scripts, '/reslice_to_222.mat']);
            case 3
                load ([basis_scripts, '/reslice_to_333.mat']);
        end
        matlabbatch{1}.spm.spatial.coreg.write.source = image_array_wa_bet;
        spm_jobman ('run', matlabbatch);
        
        %cd (basis_fmri);
        % and delete the wasting wavol*nii
        %c = 'rm wavol*.nii'; system(c);
        
    end
    
    if reslice_segments
        clear matlabbatch
        switch normalized_resolution
            case 2
                load ([basis_scripts, '/reslice_to_222.mat']);
            case 3
                load ([basis_scripts, '/reslice_to_333.mat']);
        end
        clear seg_array
        
        tmp = dir('wrc*.img');
        if size(tmp,1) > 0
           for i=1:size(tmp,1)
              cmd = ['fslchfiletype NIFTI ', tmp.name(i,:)]; system(cmd); 
           end
        end
        seg_array         = {[deblank(basis_fmri),'wrc1vol0000.nii']};
        item = [deblank(basis_fmri),'wrc2vol0000.nii'];
        seg_array        = [seg_array; item];
        item = [deblank(basis_fmri),'wrc3vol0000.nii'];
        seg_array        = [seg_array; item];
        
        matlabbatch{1}.spm.spatial.coreg.write.source = seg_array;
        spm_jobman ('run', matlabbatch);
    
    end
    
    % =================
    %   detrending (L)
    % =================
    
    if detrend_4D
        cd(basis_fmri);
        func_detrend_rwa(normalized_resolution);
    end
    
    % =================================
    %      produce std98_perc_mask
    % =================================
    
    if produce_std98
        cd(basis_fmri);
        func_produce_std98(normalized_resolution);
    end
    
    % ===========================================
    % perform bet and apply mask to all wa-images
    % ===========================================
    
    if BET
        cd (basis_fmri);
        if normalized_resolution == 2
            tmp_str = 'dr';
        elseif normalized_resolution == 3
            tmp_str = 'dr333';
        elseif normalized_resolution == 1
            tmp_str = '';
        end
%         c = sprintf('rm %swavolbet*', tmp_str); system(c);
%         for i=5:n_images-1
%             c = sprintf('/usr/lib/fsl/5.0/fslchfiletype NIFTI %swavol%04d.img',tmp_str,i); system(c);
%         end
        c = sprintf('/usr/lib/fsl/5.0/bet %swavol0005.nii %swavol0005_bet.nii -m',tmp_str, tmp_str);
        system(c);
          
        
        c = ['rm betmask', tmp_str(3:end),'.* ']; system(c);
        c = ['mv ', tmp_str, 'wavol0005_bet_mask.nii.gz betmask', tmp_str(3:end),'.nii.gz']; system(c);
%         c = ['rm ', tmp_str, 'wavol0005_bet.nii.gz']; system(c);
%         c = ['rm ', tmp_str,'wavol*.img']; system(c);
%         c = ['rm ', tmp_str, 'wavol*.hdr']; system(c);
%         c = ['/usr/lib/fsl/5.0/fslmerge -t master ', tmp_str, 'wavol????.nii']; system(c);
%         c = 'rm master_bet*'; system(c);
%         c = ['/usr/lib/fsl/5.0/fslmaths master -mas betmask', tmp_str(3:end),' master_bet']; system(c);
%         c = ['/usr/lib/fsl/5.0/fslsplit master_bet ', tmp_str, 'wavolbet']; system(c);
%         c = 'rm master*'; system(c);
%         for i = n_images-5-1:-1:0  % reverse counter to avoid overwriting
%             c = sprintf('mv %swavolbet%04d.nii.gz %swavolbet%04d.nii.gz', tmp_str, i, tmp_str, i+5);
%             system(c);
%             c = sprintf('gunzip %swavolbet%04d.nii.gz', tmp_str, i+5);
%             system(c);
%         end
    end
    
    
    % ===============================
    %   temporal filtering
    % ===============================
    
    if temporal_filter
        if normalized_resolution == 2
            tmp_str = 'dr';
        elseif normalized_resolution == 3
            tmp_str = 'dr333';
        end
        fname = [basis_fmri,'/', tmp_str, 'wavolbet*.nii'];
        outputdir = basis_fmri;
        lf  = listfiles([fname ]);
        Fs = tr;
        
        im4d_f = tfilter_nii(lf, outputdir, Fs,'N', 10, ... 
        'prefix', 'k', 'fc1', 0.01, 'fc2', 0.1, 'ftype','pass', ... 
        'ic_function', @median, 'verbose', true, 'display', true); 
    
        %cave: in MATLAB R2015a, "bandpass" is not an option, but "pass"
    end
    
    % ===============================
    %   reslice wa to rwa at 2x2x2
    % ===============================
    
    if reslice
        clear matlabbatch
        switch normalized_resolution
            case 2
                load ([basis_scripts, '/reslice_to_222.mat']);
            case 3
                load ([basis_scripts, '/reslice_to_333.mat']);
        end
        matlabbatch{1}.spm.spatial.coreg.write.source =image_array_kdrwa_bet;
        spm_jobman ('run', matlabbatch);
        
        cd (basis_fmri);
        % and delete the wasting wavol*nii
        %c = 'rm wavol*.nii'; system(c);
        
    end
    
    
    % =========
    % smoothing
    % =========
    
    if smoothing
        clear matlabbatch
        load ([basis_scripts, '/dummy_smooth.mat']);
        switch normalized_resolution
            case 2
             matlabbatch{1}.spm.spatial.smooth.data = image_array_kdrwa_bet;                
            case 3
             matlabbatch{1}.spm.spatial.smooth.data = image_array_r333kdrwa_bet;                    
        end
        matlabbatch{1}.spm.spatial.smooth.fwhm = fwhm_kernel;
        if fwhm_kernel(1) == 6
            matlabbatch{1}.spm.spatial.smooth.prefix = 's';
        elseif fwhm_kernel(1) == 8
            matlabbatch{1}.spm.spatial.smooth.prefix = 's888';            
        end
        inputs = cell(0, 1);
        jobs_all = spm_jobman('run', matlabbatch, '');
    end
        
    % =====================
    %    BG correction
    % =====================
    
    % comment: montage of thalamus/putamen into GM & WM
    
    if BG_correction
        clear matlabbatch tmp
        load ([basis_scripts, '/dummy_mean.mat']);
        calc_str = 'i1+(i2>0)+(i3>0)+(i4>0)+(i5>0)';
        switch normalized_resolution
            case 2
        tmp{1} = [basis_fmri,'/','rwrc1vol0000.nii'];
        tmp{2} = [basis_scripts,'/077_Thalamus_L.img'];
        tmp{3} = [basis_scripts,'/078_Thalamus_R.img'];
        tmp{4} = [basis_scripts,'/073_Putamen_L.img'];
        tmp{5} = [basis_scripts,'/074_Putamen_R.img'];
        tmp2 = [basis_fmri,'/rwrc1vol0000_wBG.nii'];
            case 3
        tmp{1} = [basis_fmri,'/','wr333c1vol0000.nii'];
        tmp{2} = [basis_scripts,'/r333_077_Thalamus_L.img'];
        tmp{3} = [basis_scripts,'/r333_078_Thalamus_R.img'];
        tmp{4} = [basis_scripts,'/r333_073_Putamen_L.img'];
        tmp{5} = [basis_scripts,'/r333_074_Putamen_R.img'];        
        tmp2 = [basis_fmri,'/wr333c1vol0000_wBG.nii'];
        end
        
        matlabbatch{1}.spm.util.imcalc.input = tmp';
        matlabbatch{1}.spm.util.imcalc.output = tmp2;
        matlabbatch{1}.spm.util.imcalc.expression = calc_str;
        spm_jobman ('run', matlabbatch);
        
        clear matlabbatch tmp
        load ([basis_scripts, '/dummy_mean.mat']);
        calc_str = 'i1.*(~i2).*(~i3).*(~i4).*(~i5)';
        
        switch normalized_resolution
            case 2
                tmp{1} = [basis_fmri,'/','rwrc2vol0000.nii'];
                tmp{2} = [basis_scripts,'/077_Thalamus_L.img'];
                tmp{3} = [basis_scripts,'/078_Thalamus_R.img'];
                tmp{4} = [basis_scripts,'/073_Putamen_L.img'];
                tmp{5} = [basis_scripts,'/074_Putamen_R.img'];
                tmp2 = [basis_fmri,'/rwrc2vol0000_woBG.nii'];
            case 3
                tmp{1} = [basis_fmri,'/','wr333c2vol0000.nii'];
                tmp{2} = [basis_scripts,'/r333_077_Thalamus_L.img'];
                tmp{3} = [basis_scripts,'/r333_078_Thalamus_R.img'];
                tmp{4} = [basis_scripts,'/r333_073_Putamen_L.img'];
                tmp{5} = [basis_scripts,'/r333_074_Putamen_R.img'];
                tmp2 = [basis_fmri,'/wr333c2vol0000_woBG.nii'];
        end

        matlabbatch{1}.spm.util.imcalc.input = tmp';
        matlabbatch{1}.spm.util.imcalc.output = tmp2;
        matlabbatch{1}.spm.util.imcalc.expression = calc_str;
        spm_jobman ('run', matlabbatch);
               
%         clear matlabbatch
%         load ([basis_scripts, '/reslice_to_222.mat']);
%         image_array_tmp = {[basis_fmri,'/wrc1vol0000_wBG.nii,1']};
%         image_array_tmp = [image_array_tmp; [basis_fmri,'/wrc2vol0000_noBG.nii,1']];
%         image_array_tmp = [image_array_tmp; [basis_fmri,'/wrc3vol0000.nii,1']];
%         matlabbatch{1}.spm.spatial.coreg.write.source = image_array_tmp;
%         spm_jobman ('run', matlabbatch);
        
    end
    
    % make sure spm restarts properly
    if restart_spm
        spm quit;
    end
end

%firstlevel_residualization_RSN
