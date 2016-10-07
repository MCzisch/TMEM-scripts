% do FSL scrubbing TMEM RSN
% written by MC, 06.05.2015
% requires FSL5.0

clear
warning off all

cd /media/spectro/_media_spectro/RAID/data/TMEM_RSN/autoscripts;
fid = fopen('scrubbing_logfile.log','a');
tmp = datestr(now);
fprintf(fid,'\n\n\n%s',tmp);

basis_fmri = '/media/spectro/_media_spectro/RAID/data/TMEM_RSN/';
cd(basis_fmri);

% data input is run via a file list
fid = fopen('file_list.txt','r');
tline = fgetl(fid);

while ischar(tline)

    dir = tline(1:7);   % here, only the first part of the file_list entries is relevant, e.g. RSN_002
    cd (dir);
    
    fprintf(fid, 'running %s\n', dir)
    
    out_cmd = 'rm raw4D';
    system(out_cmd);
    out_cmd = 'fslmerge -t raw4D vol*nii';  % scrubbing is done on the raw images
    system(out_cmd);
    out_cmd = 'fsl_motion_outliers -i raw4D.nii -o fsl_motion_outliers.txt -s metric_val.txt -p metric_graph --dummy=5 -v';
    system(out_cmd);
    out_cmd = 'rm raw4D.nii.gz';
    system(out_cmd);
    tmp = load('fsl_motion_outliers.txt');
    number_of_scrubbed = size(tmp,2);
    if number_of_scrubbed > 20
        fprintf(fid, 'WARNING\t')
    end
    fprintf(fid, '%s: %d\n', dir, number_of_scrubbed)
    
    cd(basis_fmri);
    tline = fgetl(fid);
end
fprintf(fid,'E N D\n\n');

fclose(fid);

