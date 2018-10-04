function DatViz_extramaps

% create mean anat image

clear variables
%% parameters 
localspm = which('spm');
if isempty(localspm)
    error('SPM is not present in your matlab path');
end
spm('defaults', 'FMRI');
spm_jobman('initcfg')

%% create mean anatomical
datapath = uigetdir(pwd,'Select the face data set directory');
cd([datapath filesep 'BIDS_processed'])
sub = dir('sub-*');

data = NaN([79 95 79 size(sub,1)]);
for subject=1:size(sub,1)
    cd([sub(subject).name filesep 'anat'])
    wm       = dir('wm*.nii');
    wc       = dir('wc*.nii');
    head     = spm_read_vols(spm_vol(wm.name));
    tissues  = spm_read_vols(spm_vol(wc(1).name)) + ...
        spm_read_vols(spm_vol(wc(2).name)) + ...
        spm_read_vols(spm_vol(wc(3).name));
    mask = tissues > 0;
    data(:,:,:,subject) = head.*mask;
    cd([datapath filesep 'BIDS_processed'])
end

avg = mean(data,4);
w = spm_vol([sub(subject).name filesep 'anat' filesep wm.name]);
w.fname = [datapath filesep 'BIDS_processed' filesep 'subjects_anat_avg.nii'];
w.descrip = 'average of warped bias corrected anatomical images';
spm_write_vol(w,avg);


%% generate subject inverse warped images for vizualization in subject space

cd([datapath filesep 'BIDS_processed'])
for subject=1:size(sub,1)
    clear matlabbatch
    T1 = dir([sub(subject).name filesep 'anat' filesep 'sub-*.nii']);
    matlabbatch{1}.spm.util.bbox.image = {[sub(subject).name filesep 'anat' filesep T1(1).name]};
    matlabbatch{1}.spm.util.bbox.bbdef.fov = 'fv';
    out = spm_jobman('run',matlabbatch);
    invf = dir([sub(subject).name filesep 'anat' filesep 'iy*.nii']);
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {[sub(subject).name filesep 'anat' filesep invf.name]};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {[sub(subject).name filesep 'stats' filesep 'con_0002.nii,1']};
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = out{1}.bb;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
    spm_jobman('run',matlabbatch);
    clear matlabbatch
end


