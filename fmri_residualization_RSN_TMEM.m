% ==========================================
% Analysis 2016; TMEM: RSN
%
%
% 29.09.2016    Residualization
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

residualization   = 1;
transform_to_floats = 1;

% =====================================================
%
%  end of user input section
%
% =====================================================
global basis_scripts;
basis_scripts = '/media/spectro/_media_spectro/RAID/data/TMEM_RSN/autoscripts/';
basis_dir = '/media/spectro/_media_spectro/RAID/data/TMEM_RSN/';

% =====================================================

master = dir('RSN_*');

master_start = 1;
master_end = size(master,1);

% 1, 10

for i = 11:master_end
    
    % make sure spm is up
    if restart_spm
        spm('defaults','fMRI');
        spm_jobman('initcfg');
        spm fmri;
    end
    
    
    fprintf(['### Working on files ', master(i).name, ' ... ###\n\n']);
    basis_fmri = [basis_dir,'/',master(i).name,'/'];
    cd (basis_fmri);
    
    if residualization
        
        cmd = ('mkdir 1stlevel_Res'); system(cmd);
        c = sprintf('rm 1stlevel_Res/SPM.mat'); system(c);
        
        cd(basis_fmri);
        tmp = dir('drwavolbet*.img');
        n_images = size(tmp,1);
        
        % ==============
        % define arrays
        % ==============
        
        clear image_array
        
        image_array         = {[deblank(basis_fmri),'drwavolbet0005.img']};
        for i=6:n_images+5-1
            tmp = sprintf('drwavolbet%04d',i);
            item = [deblank(basis_fmri),tmp,'.img'];
            image_array        = [image_array; item];
        end
        
        nuisance = load('rp_avol0005.txt');
        tmp = load('fpca_wm_50_components.txt');
        tmp2 = tmp(:,1:3);
        nuisance = [nuisance tmp2];
        tmp = load('fpca_csf_50_components.txt');
        tmp2 = tmp(:,1:3);
        nuisance = [nuisance tmp2];
        tmp = load('fpca_std98_50_components.txt');
        tmp2 = tmp(:,1:3);
        nuisance = [nuisance tmp2];
        tmp = load('fsl_motion_outliers.txt');
        nuisance = [nuisance tmp];
        
        %nuisance =nuisance';
        save('Res_nuisance.txt','nuisance','-ascii');
        
        
        clear matlabbatch
        load ([basis_scripts, '/dummy_firstlevel.mat']);
        
        outfile = [basis_fmri,'1stlevel_Res'];      % directory
        matlabbatch{1}.spm.stats.fmri_spec.dir = {outfile};
        
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans = image_array;
        
        outfile = strcat(basis_fmri,'/Res_nuisance.txt');
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {outfile};
        
        spm_jobman ('run', matlabbatch);
        
        clear matlabbatch
        load ([basis_scripts, '/dummy_run_residualization.mat']);
        matlabbatch{1}.spm.stats.fmri_est.spmmat = {[basis_fmri,'1stlevel_Res/SPM.mat']};

        spm_jobman ('run', matlabbatch);
        
    end
    
    if transform_to_floats
        
        for i = 1:185
            
            tmp = sprintf('1stlevel_Res/Res_%04d.nii',i);
            orig_data = load_nii(tmp);
            orig_data.img = orig_data.img + 1000;
            tmp = sprintf('1stlevel_Res/floatRes_%04d.nii',i);
            orig_data.hdr.dime.datatype = 16;
            orig_data.hdr.dime.bitpix = 16;
            save_nii(orig_data,tmp);
            
            
        end
        
        cmd = ['rm 1stlevel_Res/Res_*.nii']; system(cmd);
        
    end
    
end

