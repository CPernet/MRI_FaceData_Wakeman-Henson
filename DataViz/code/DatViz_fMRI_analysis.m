% This is the source code for the article Data visualization for inference in brain imaging
% by Christopher R. Madan & Cyril R. Pernet
%
% This script rerun the fMRI data analysis, including lateralization, data
% extraction and plots.

%% the data used are the fmri data of the multimodal face repetition experiment from Rik Henson
% Details about the data and experiment can be found here https://www.nature.com/articles/sdata20151
% Data were download in following BIDS from ftp://ftp.mrc-cbu.cam.ac.uk/personal/rik.henson/wakemandg_hensonrn/BIDS/
% We use SPM12, the lateralization toolbox, the anatomy toolbox:
% http://www.fil.ion.ucl.ac.uk/spm/
% https://www.medizin.uni-tuebingen.de/kinder/en/research/neuroimaging/software/
% http://www.fz-juelich.de/inm/inm-1/EN/Forschung/_docs/SPMAnatomyToolbox/SPMAnatomyToolbox_node.html
% There is another set of functions distributed here, coming from the
% robust statistical toolbox - a library of matlab function for robust
% stats and data viualization https://github.com/CPernet/Robust_Statistical_Toolbox
%
% note the way I use BIDS here is old fashion better to follow 
% https://en.wikibooks.org/wiki/SPM/BIDS#BIDS_parser_and_queries
% Cyril Pernet 26 Septembre 2018

clear variables
%% parameters for the code and libraries
localspm = which('spm');
if isempty(localspm)
    error('SPM is not present in your matlab path');
end

if exist([fileparts(localspm) filesep 'toolbox' filesep 'LI' filesep 'LI.m'],'file')
    if isempty(which('li'))
        addpath([fileparts(localspm) filesep 'toolbox' filesep 'LI'])
    end
else
    error('LI.m (lateralization toolbox) is not in your matlab path');
end

if exist([fileparts(localspm) filesep 'toolbox' filesep 'Anatomy' filesep 'PMaps' filesep 'Visual_FG1.nii'],'file')
    fusiform_mask{1}=[fileparts(localspm) filesep 'toolbox' filesep 'Anatomy' filesep 'PMaps' filesep 'Visual_FG1.nii,1'];
    fusiform_mask{2}=[fileparts(localspm) filesep 'toolbox' filesep 'Anatomy' filesep 'PMaps' filesep 'Visual_FG2.nii,1'];
    fusiform_mask{3}=[fileparts(localspm) filesep 'toolbox' filesep 'Anatomy' filesep 'PMaps' filesep 'Visual_FG3.nii,1'];
    
else
    error('The Anatomy Toolbox is not present in your SPM toolbox directory');
end

datapath = uigetdir(pwd,'Select the face data set directory');

%% get the data info and unpack
% -------------------------------------------------------------------------
spm('defaults', 'FMRI');
spm_jobman('initcfg')
cd(datapath); BIDS = spm_BIDS(datapath);
BIDS.dir = datapath; 
outdir = 'BIDS_processed';
parpool(feature('numCores')-1); % use all available cores -1

% make a conditional statement, assuming that if you have done it before,
% this is correct and folder structure is untouched
if ~exist(outdir,'dir')
    mkdir(outdir);
    disp('starting to unpack data')
    parfor s=1:size(BIDS.subjects,2)
        
        in = [BIDS.dir filesep BIDS.subjects(s).name filesep 'anat' filesep BIDS.subjects(s).anat.filename];
        subjects{s}.anat = [BIDS.dir filesep outdir filesep BIDS.subjects(s).name filesep 'anat' filesep BIDS.subjects(s).anat.filename(1:end-3)];
        fprintf('subject %g: unpacking anatomical data \n',s)
        gunzip(in, [outdir filesep BIDS.subjects(s).name filesep 'anat' ]);
        
        for frun = 1:size(BIDS.subjects(s).func,2)/2
            in = [BIDS.dir filesep BIDS.subjects(s).name filesep 'func' filesep BIDS.subjects(s).func(frun).filename];
            subjects{s}.func(frun,:) = [outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep BIDS.subjects(s).func(frun).filename(1:end-3)];
            fprintf('subject %g: unpacking functional data \n',s)
            gunzip(in, [outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun)]);
        end
    end
else % as above but no point having parallel overhaed to create a list
    for s=1:size(BIDS.subjects,2)
        in = [BIDS.dir filesep BIDS.subjects(s).name filesep 'anat' filesep BIDS.subjects(s).anat.filename];
        subjects{s}.anat = [BIDS.dir filesep outdir filesep BIDS.subjects(s).name filesep 'anat' filesep BIDS.subjects(s).anat.filename(1:end-3)];

        for frun = 1:size(BIDS.subjects(s).func,2)/2
            in = [BIDS.dir filesep BIDS.subjects(s).name filesep 'func' filesep BIDS.subjects(s).func(frun).filename];
            subjects{s}.func(frun,:) = [outdir filesep BIDS.subjects(s).name filesep 'run' num2str(frun) filesep BIDS.subjects(s).func(frun).filename(1:end-3)];
        end
    end
end

%% prepare a batch job per suject then run them in parallel
batch = cell(1,size(BIDS.subjects,2));
for s=1:size(BIDS.subjects,2)
    % each subject build a job structure around matlabbatch
    N_run = size(BIDS.subjects(s).func,2)/2;
    if s==10; N = N_run-1; else N = N_run; end

    %%%% step 1 realign
    for frun = 1:N % each run
        % each run input raw data
        for v=1:208
            filesin{v} = [BIDS.dir filesep subjects{s}.func(frun,:) ',' num2str(v)];
        end
        matlabbatch{1}.spm.spatial.realign.estwrite.data{frun}   = filesin';
    end
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep     = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm    = 5;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm     = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp  = 2;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap    = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight  = '';
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which   = [0 1];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp  = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap    = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask    = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix  = 'r';
    
    %%%% step 2 segment
    matlabbatch{2}.spm.spatial.preproc.channel.vols     = {subjects{s}.anat};
    matlabbatch{2}.spm.spatial.preproc.channel.biasreg  = 0.001;
    matlabbatch{2}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{2}.spm.spatial.preproc.channel.write    = [0 1];
    matlabbatch{2}.spm.spatial.preproc.tissue(1).tpm    = {[ fileparts(localspm) filesep 'tpm' filesep 'TPM.nii,1']};
    matlabbatch{2}.spm.spatial.preproc.tissue(1).ngaus  = 1;
    matlabbatch{2}.spm.spatial.preproc.tissue(1).native = [1 0];
    matlabbatch{2}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch{2}.spm.spatial.preproc.tissue(2).tpm    = {[ fileparts(localspm) filesep 'tpm' filesep 'TPM.nii,2']};
    matlabbatch{2}.spm.spatial.preproc.tissue(2).ngaus  = 1;
    matlabbatch{2}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{2}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{2}.spm.spatial.preproc.tissue(3).tpm    = {[ fileparts(localspm) filesep 'tpm' filesep 'TPM.nii,3']};
    matlabbatch{2}.spm.spatial.preproc.tissue(3).ngaus  = 2;
    matlabbatch{2}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{2}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{2}.spm.spatial.preproc.tissue(4).tpm    = {[ fileparts(localspm) filesep 'tpm' filesep 'TPM.nii,4']};
    matlabbatch{2}.spm.spatial.preproc.tissue(4).ngaus  = 3;
    matlabbatch{2}.spm.spatial.preproc.tissue(4).native = [0 0];
    matlabbatch{2}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{2}.spm.spatial.preproc.tissue(5).tpm    = {[ fileparts(localspm) filesep 'tpm' filesep 'TPM.nii,5']};
    matlabbatch{2}.spm.spatial.preproc.tissue(5).ngaus  = 4;
    matlabbatch{2}.spm.spatial.preproc.tissue(5).native = [0 0];
    matlabbatch{2}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{2}.spm.spatial.preproc.tissue(6).tpm    = {[ fileparts(localspm) filesep 'tpm' filesep 'TPM.nii,6']};
    matlabbatch{2}.spm.spatial.preproc.tissue(6).ngaus  = 2;
    matlabbatch{2}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{2}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{2}.spm.spatial.preproc.warp.mrf         = 1;
    matlabbatch{2}.spm.spatial.preproc.warp.cleanup     = 1;
    matlabbatch{2}.spm.spatial.preproc.warp.reg         = [0 0.001 0.5 0.05 0.2];
    matlabbatch{2}.spm.spatial.preproc.warp.affreg      = 'mni';
    matlabbatch{2}.spm.spatial.preproc.warp.fwhm        = 0;
    matlabbatch{2}.spm.spatial.preproc.warp.samp        = 3;
    matlabbatch{2}.spm.spatial.preproc.warp.write       = [1 1];
    
    %%%% step 3: coregister mean to T1 apply to all EPI (ie upsample)
    
    matlabbatch{3}.spm.spatial.coreg.estimate.ref(1)            = cfg_dep('Segment: Bias Corrected (1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
    matlabbatch{3}.spm.spatial.coreg.estimate.source(1)         = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
    matlabbatch{3}.spm.spatial.coreg.estimate.other(1)          = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','cfiles'));
    matlabbatch{3}.spm.spatial.coreg.estimate.other(2)          = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 2)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{2}, '.','cfiles'));
    matlabbatch{3}.spm.spatial.coreg.estimate.other(3)          = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 3)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{3}, '.','cfiles'));
    matlabbatch{3}.spm.spatial.coreg.estimate.other(4)          = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 4)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{4}, '.','cfiles'));
    matlabbatch{3}.spm.spatial.coreg.estimate.other(5)          = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 5)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{5}, '.','cfiles'));
    matlabbatch{3}.spm.spatial.coreg.estimate.other(6)          = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 6)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{6}, '.','cfiles'));
    matlabbatch{3}.spm.spatial.coreg.estimate.other(7)          = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 7)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{7}, '.','cfiles'));
    matlabbatch{3}.spm.spatial.coreg.estimate.other(8)          = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 8)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{8}, '.','cfiles'));
    if s~=10
        matlabbatch{3}.spm.spatial.coreg.estimate.other(9)      = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 9)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{9}, '.','cfiles'));
    end
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.sep      = [4 2];
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.fwhm     = [7 7];
    
    %%%% step 4 normalize EPI
    
    matlabbatch{4}.spm.spatial.normalise.write.subj.def(1)      = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
    matlabbatch{4}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
    matlabbatch{4}.spm.spatial.normalise.write.woptions.bb      = [-78 -112 -70; 78 76 85];
    matlabbatch{4}.spm.spatial.normalise.write.woptions.vox     = [2 2 2];
    matlabbatch{4}.spm.spatial.normalise.write.woptions.interp  = 4;
    matlabbatch{4}.spm.spatial.normalise.write.woptions.prefix  = 'w';
    
    %%%% step 5 smooth
    matlabbatch{5}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
    matlabbatch{5}.spm.spatial.smooth.fwhm    = [8 8 8];
    matlabbatch{5}.spm.spatial.smooth.dtype   = 0;
    matlabbatch{5}.spm.spatial.smooth.im      = 0;
    matlabbatch{5}.spm.spatial.smooth.prefix  = 's';
    
    matlabbatch{6}.cfg_basicio.file_dir.file_ops.file_move.files(1)      = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','cfiles'));
    matlabbatch{6}.cfg_basicio.file_dir.file_ops.file_move.files(2)      = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 2)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{2}, '.','cfiles'));
    matlabbatch{6}.cfg_basicio.file_dir.file_ops.file_move.files(3)      = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 3)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{3}, '.','cfiles'));
    matlabbatch{6}.cfg_basicio.file_dir.file_ops.file_move.files(4)      = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 4)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{4}, '.','cfiles'));
    matlabbatch{6}.cfg_basicio.file_dir.file_ops.file_move.files(5)      = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 5)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{5}, '.','cfiles'));
    matlabbatch{6}.cfg_basicio.file_dir.file_ops.file_move.files(6)      = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 6)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{6}, '.','cfiles'));
    matlabbatch{6}.cfg_basicio.file_dir.file_ops.file_move.files(7)      = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 7)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{7}, '.','cfiles'));
    matlabbatch{6}.cfg_basicio.file_dir.file_ops.file_move.files(8)      = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 8)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{8}, '.','cfiles'));
    if s~=10
        matlabbatch{6}.cfg_basicio.file_dir.file_ops.file_move.files(9)  = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 9)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{9}, '.','cfiles'));
    end
    matlabbatch{6}.cfg_basicio.file_dir.file_ops.file_move.files(10)     = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
    matlabbatch{6}.cfg_basicio.file_dir.file_ops.file_move.files(11)     = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
    matlabbatch{6}.cfg_basicio.file_dir.file_ops.file_move.action.delete = false;
    
    % normalize T1 and tissue classes
    matlabbatch{7}.spm.spatial.normalise.write.subj.def(1)      = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
    matlabbatch{7}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Segment: Bias Corrected (1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
    matlabbatch{7}.spm.spatial.normalise.write.subj.resample(2) = cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
    matlabbatch{7}.spm.spatial.normalise.write.subj.resample(3) = cfg_dep('Segment: c2 Images', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
    matlabbatch{7}.spm.spatial.normalise.write.subj.resample(4) = cfg_dep('Segment: c3 Images', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{3}, '.','c', '()',{':'}));
    matlabbatch{7}.spm.spatial.normalise.write.woptions.bb      = [-78 -112 -70; 78 76 85];
    matlabbatch{7}.spm.spatial.normalise.write.woptions.vox     = [2 2 2];
    matlabbatch{7}.spm.spatial.normalise.write.woptions.interp  = 4;
    matlabbatch{7}.spm.spatial.normalise.write.woptions.prefix  = 'w';
    
    % simply save the batch into a cell array
    save([datapath filesep 'code' filesep 'preproc_batch_' BIDS.subjects(s).name],'matlabbatch')
    batch{s} = matlabbatch; clear matlabbatch
end

% run all the batches in parallel
all_jobs = cell(1,length(batch));
parfor s=1:length(batch)
    try
       %  spm_jobman('initcfg');
        all_jobs{s} = spm_jobman('run',batch{s});
        test{s} = 'done';
    catch
        test{s} = 'failed';
    end
end
cd([datapath filesep 'code']); save all_jobs all_jobs

batch = cell(1,size(BIDS.subjects,2));
for s=1:size(BIDS.subjects,2)
    % each subject build a job structure around matlabbatch
    N_run = size(BIDS.subjects(s).func,2)/2;
    
    matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.parent = {[BIDS.dir filesep fileparts(fileparts(subjects{s}.func(1,:)))]};
    matlabbatch{1}.cfg_basicio.file_dir.dir_ops.cfg_mkdir.name = 'stats';
    
    %%%% step 6 build the model
    matlabbatch{2}.spm.stats.fmri_spec.dir(1)         = cfg_dep('Make Directory: Make Directory ''stats''', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','dir'));
    matlabbatch{2}.spm.stats.fmri_spec.timing.units   = 'secs';
    matlabbatch{2}.spm.stats.fmri_spec.timing.RT      = 2;
    matlabbatch{2}.spm.stats.fmri_spec.timing.fmri_t  = 16;
    matlabbatch{2}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
        
    %%%% step 6 stats
    if s==10; N = N_run-1; else N = N_run; end

    for frun = 1:N
        for v=1:208
            [filepath,name,ext] = fileparts(subjects{s}.func(frun,:));
            filesin{v} = [BIDS.dir filesep  filepath filesep 'sw' name ext ',' num2str(v)];
        end
        matlabbatch{2}.spm.stats.fmri_spec.sess(frun).scans = filesin';
        
        N_events = length(BIDS.subjects(s).func(frun+N_run).meta.trial_type);
        cond = unique(BIDS.subjects(s).func(frun+N_run).meta.trial_type);
        all_cond{frun} = cond;
        N_cond = length(cond);
               
        onsets = NaN(N_events,N_cond);
        durations = NaN(N_events,N_cond);
        for C = 1:N_cond
            for n=1:N_events
                if strcmp(cell2mat(BIDS.subjects(s).func(frun+N_run).meta.trial_type(n,:)),cond{C})
                    onsets(n,C) = BIDS.subjects(s).func(frun+N_run).meta.onset(n);
                    durations(n,C) = BIDS.subjects(s).func(frun+N_run).meta.duration(n);
                end
            end
            matlabbatch{2}.spm.stats.fmri_spec.sess(frun).cond(C).name = cond{C};
            matlabbatch{2}.spm.stats.fmri_spec.sess(frun).cond(C).onset = onsets(~isnan(onsets(:,C)),C);
            matlabbatch{2}.spm.stats.fmri_spec.sess(frun).cond(C).duration = durations(~isnan(durations(:,C)),C);
            matlabbatch{2}.spm.stats.fmri_spec.sess(frun).cond(C).tmod = 0;
            matlabbatch{2}.spm.stats.fmri_spec.sess(frun).cond(C).pmod = struct('name', {}, 'param', {}, 'poly', {});
            matlabbatch{2}.spm.stats.fmri_spec.sess(frun).cond(C).orth = 1;
        end
        matlabbatch{2}.spm.stats.fmri_spec.sess(frun).multi        = {''};
        matlabbatch{2}.spm.stats.fmri_spec.sess(frun).regress      = struct('name', {}, 'val', {});
        matlabbatch{2}.spm.stats.fmri_spec.sess(frun).multi_reg(1) = all_jobs{s}{1}.sess(frun).rpfile;
        matlabbatch{2}.spm.stats.fmri_spec.sess(frun).hpf          = 128;
    end
    
    matlabbatch{2}.spm.stats.fmri_spec.fact             = struct('name', {}, 'levels', {});
    matlabbatch{2}.spm.stats.fmri_spec.bases.hrf.derivs = [1 1];
    matlabbatch{2}.spm.stats.fmri_spec.volt             = 1;
    matlabbatch{2}.spm.stats.fmri_spec.global           = 'None';
    matlabbatch{2}.spm.stats.fmri_spec.mthresh          = 0.8;
    matlabbatch{2}.spm.stats.fmri_spec.mask             = {''};
    matlabbatch{2}.spm.stats.fmri_spec.cvi              = 'AR(1)';
    
    % estimate
    matlabbatch{3}.spm.stats.fmri_est.spmmat(1)        = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{3}.spm.stats.fmri_est.write_residuals  = 0;
    matlabbatch{3}.spm.stats.fmri_est.method.Classical = 1;
    
    % contrasts
    matlabbatch{4}.spm.stats.con.spmmat(1)               = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{4}.spm.stats.con.consess{1}.fcon.name    = 'canonical effect of interest';
    matlabbatch{4}.spm.stats.con.consess{1}.fcon.weights = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;  0 0 0 1 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0];
    matlabbatch{4}.spm.stats.con.consess{1}.fcon.sessrep = 'repl';
    matlabbatch{4}.spm.stats.con.consess{2}.tcon.name    = 'Faces > Scrambled Faces';
    matlabbatch{4}.spm.stats.con.consess{2}.tcon.weights = [0.5 0 0 0.5 0 0 -1 0 0 0 0 0 0 0 0];
    matlabbatch{4}.spm.stats.con.consess{2}.tcon.sessrep = 'repl';
    matlabbatch{4}.spm.stats.con.consess{3}.tcon.name    = 'Famous';
    matlabbatch{4}.spm.stats.con.consess{3}.tcon.weights = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    matlabbatch{4}.spm.stats.con.consess{3}.tcon.sessrep = 'repl';
    matlabbatch{4}.spm.stats.con.consess{4}.tcon.name    = 'Unfamiliar';
    matlabbatch{4}.spm.stats.con.consess{4}.tcon.weights = [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0];
    matlabbatch{4}.spm.stats.con.consess{4}.tcon.sessrep = 'repl';
    matlabbatch{4}.spm.stats.con.consess{5}.tcon.name    = 'Scrambled';
    matlabbatch{4}.spm.stats.con.consess{5}.tcon.weights = [0 0 0 0 0 0 1 0 0 0 0 0 0 0 0];
    matlabbatch{4}.spm.stats.con.consess{5}.tcon.sessrep = 'repl';
    matlabbatch{4}.spm.stats.con.delete = 1;
    
    % simply save the batch into a cell array
    save([datapath filesep 'code' filesep 'stat_batch_' BIDS.subjects(s).name],'matlabbatch')
    batch{s} = matlabbatch; clear matlabbatch
end

% run all the batches in parallel
% subject 10 is trouble - could be issue in downloading ?? anyway I just
% removed it here
stat_jobs = cell(1,length(batch));
parfor s=1:length(batch)
    try
       %  spm_jobman('initcfg');
        stat_jobs{s} = spm_jobman('run',batch{s});
        test{s} = 'done';
    catch
        test{s} = 'failed';
    end
end
cd([datapath filesep 'code']); save stat_jobs stat_jobs

% close parallel computing toolbox
try delete(gcp('nocreate')); end 

%% group level ANOVA

cd([datapath filesep outdir]);
mkdir('Rep_ANOVA'); sindex = 1;
matlabbatch{1}.spm.stats.factorial_design.dir = {[pwd filesep 'Rep_ANOVA']};
for s=1:16
    for c=1:3
        subject_scans(c,:) = {[pwd filesep BIDS.subjects(s).name filesep 'stats' filesep 'con_000' num2str(c+2) '.nii']};
    end
    matlabbatch{1}.spm.stats.factorial_design.des.anovaw.fsubject(sindex).scans = subject_scans;
    matlabbatch{1}.spm.stats.factorial_design.des.anovaw.fsubject(sindex).conds = [1 2 3];
    sindex = sindex+1;
end
matlabbatch{1}.spm.stats.factorial_design.des.anovaw.dept        = 1;
matlabbatch{1}.spm.stats.factorial_design.des.anovaw.variance    = 1;
matlabbatch{1}.spm.stats.factorial_design.des.anovaw.gmsca       = 0;
matlabbatch{1}.spm.stats.factorial_design.des.anovaw.ancova      = 0;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none     = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im             = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em             = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit         = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm        = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.method.Classical               = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));

matlabbatch{3}.spm.stats.con.consess{1}.fcon.name    = 'Effects of Interest';
matlabbatch{3}.spm.stats.con.consess{1}.fcon.weights = [1 0 0 repmat(1/16,1,16)
                                                        0 1 0 repmat(1/16,1,16)
                                                        0 0 1 repmat(1/16,1,16)];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name    = 'Faces (Fam+Unf) > Scrambled';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0.5 0.5 -1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{3}.tcon.name    = 'Famous';
matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [1 0 0 ones(1,16).*1/16];
matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{4}.tcon.name    = 'Unfamiliar';
matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0 1 0 ones(1,16).*1/16];
matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{5}.tcon.name    = 'Scrambled';
matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [0 0 1 ones(1,16).*1/16];
matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;
spm_jobman('run',matlabbatch);

%% extract a-priori ROI parameters

% coordinates from Neurosynth
priors = [-40 -50 -18; 45 -50 -22; ...
    -38 -84 -10; 44 -78 -12; ...
    -26, 0, -30; 34 -8 -32; ...
    -20 -6 -18; 22 -6, -18; ...
    -56 -54 8; 54 -44 8; ...
    42 14 22];

% contrast images per subject
cd([datapath filesep 'BIDS_processed'])
all = dir('sub-*');
clear files cfiles
for s=1:size(all,1)
    files(s,:) = [pwd filesep all(s).name filesep 'stats' filesep 'con_0002.nii'];
end

% for each coordinate get the mean activation in a 4mm sphere
clear results
for r=1:size(priors,1)
    results(:,r) = spm_summarise(files,struct('def','sphere', 'spec',4, 'xyz',priors(r,:)'),@mean);
end

% for each coordinate get the mean constant term in a 4mm sphere (
clear constant
for s=1:size(all,1)
    index = 1;
    for run = 1:9
        if s==10 && run<9
            cfiles(index,:) = [pwd filesep all(s).name filesep 'stats' filesep 'beta_0' num2str(120+run) '.nii'];
            index = index +1;
        elseif s~=10
            cfiles(index,:) = [pwd filesep all(s).name filesep 'stats' filesep 'beta_0' num2str(135+run) '.nii'];
            index = index +1;
        end
    end
    
    for r=1:size(priors,1)
        constant(s,r) = mean(spm_summarise(cfiles,struct('def','sphere', 'spec',4, 'xyz',priors(r,:)'),@mean));
    end
    clear cfiles
end

% percentage signal change is 

% same design for all, so get scaling factor from subject 1
load([datapath filesep 'BIDS_processed' filesep BIDS.subjects(1).name filesep 'stats' filesep 'SPM.mat'])
xBF.dt     = SPM.xBF.dt;
xBF.name   = SPM.xBF.name;
xBF.length = SPM.xBF.length;
xBF.order  = SPM.xBF.order;
xBF        = spm_get_bf(xBF);
event      = xBF.bf*[1 1 1]';
SF         = max(event);
PSC        = ((results.*SF)./constant).*100;

% plot results
[mean_activations,CI_activations]=rst_data_plot(PSC,'estimator','mean');
csvwrite('PSC_apriori_roi.csv',PSC);

%% extract data and make rfx plot

observed_coord = [-22 -80 -12; ...
    -26 -84 16; ...
    26 -92 26; ...
    26 -56 -10];

Y1 = gp_event_plot([-22 -80 -12]',[datapath filesep outdir filesep 'Rep_ANOVA' filesep 'SPM.mat'],'mm'); 
con1 = Y1.individual_parameters*[0.5 0.5 -1]';
save Y1 Y1; % left fusiform

Y2 = gp_event_plot([-26 -84 16]',[datapath filesep outdir filesep 'Rep_ANOVA' filesep 'SPM.mat'],'mm'); 
con2 = Y2.individual_parameters*[0.5 0.5 -1]';
save Y2 Y2; % left MOG

Y3 = gp_event_plot([26 -92 26]',[datapath filesep outdir filesep 'Rep_ANOVA' filesep 'SPM.mat'],'mm'); 
con3 = Y3.individual_parameters*[0.5 0.5 -1]';
save Y3 Y3; % right MOG

Y4 = gp_event_plot([26 -56 -10]',[datapath filesep outdir filesep 'Rep_ANOVA' filesep 'SPM.mat'],'mm'); 
con4 = Y4.individual_parameters*[0.5 0.5 -1]';
save Y4 Y4; % right fusiform

[mean_conestimate,CI]=rst_data_plot([con1 con2 con3 con4],'estimator','mean');
title('Faces>Scrambled left/right MOG & Fusiform')

clear constant
for s=1:size(all,1)
    index = 1;
    for run = 1:9
        if s==10 && run<9
            cfiles(index,:) = [pwd filesep all(s).name filesep 'stats' filesep 'beta_0' num2str(120+run) '.nii'];
            index = index +1;
        elseif s~=10
            cfiles(index,:) = [pwd filesep all(s).name filesep 'stats' filesep 'beta_0' num2str(135+run) '.nii'];
            index = index +1;
        end
    end
    
    for r=1:size(observed_coord,1)
        constant(s,r) = mean(spm_summarise(cfiles,struct('def','sphere', 'spec',4, 'xyz',observed_coord(r,:)'),@mean));
    end
    clear cfiles
end

% re-express contrasct as PSC for consistancy
PSC = (([con1 con2 con3 con4].*SF)./constant).*100;
csvwrite('contrast_estimates.csv',[con1 con2 con3]);


%% run lateralization indices
mask = [];
for m=1:3
    V = spm_vol(fusiform_mask{m});
    if isempty(mask)
        mask = spm_read_vols(V);
    else
        mask = mask+spm_read_vols(V);
    end
end
V.fname = [datapath filesep 'code' filesep 'fusiform_mask.nii'];  
V.descrip = 'SPM Anatomy toolbox PMap Visual_FG1+2+3';
spm_write_vol(V,mask>0);
fname = spmup_resize(V.fname,[-90 -126 -72;90 90 108],[2 2 2]); % rise to match the LI bounding box

li_all_brain      = NaN(1,16);
li_fusiform       = NaN(1,16);
li_fus_curves     = cell(1,16);
li_stim_all_brain = NaN(16,3);
li_stim_fusiform  = NaN(16,3);
for con = 2:5
    fprintf('running lateralization index toolbox on con image %g\n',con)
    for s=1:size(all,1)
        clear matlabbatch        
        % 1 - fusiform
        matlabbatch{1}.spm.tools.LI_cfg.spmT        = {[pwd filesep all(s).name filesep 'stats' filesep 'spmT_000' num2str(con) '.nii']};
        matlabbatch{1}.spm.tools.LI_cfg.inmask.im11 = fname; % [datapath filesep 'code' filesep 'rfusiform_mask.nii'];
        matlabbatch{1}.spm.tools.LI_cfg.exmask.em1  = 1;
        matlabbatch{1}.spm.tools.LI_cfg.method.thr7 = 1;
        matlabbatch{1}.spm.tools.LI_cfg.pre         = 0;
        matlabbatch{1}.spm.tools.LI_cfg.op          = 4;
        matlabbatch{1}.spm.tools.LI_cfg.vc          = 0;
        matlabbatch{1}.spm.tools.LI_cfg.ni          = 1;
        matlabbatch{1}.spm.tools.LI_cfg.outfile     = ['li_bootfus_con' num2str(con) all(s).name '.txt'];
        
        % 2 - full brain
        matlabbatch{2}.spm.tools.LI_cfg.spmT        = {[pwd filesep all(s).name filesep 'stats' filesep 'spmT_000' num2str(con) '.nii']};
        matlabbatch{2}.spm.tools.LI_cfg.inmask.im10 = 1;
        matlabbatch{2}.spm.tools.LI_cfg.exmask.em1  = 1;
        matlabbatch{2}.spm.tools.LI_cfg.method.thr7 = 1;
        matlabbatch{2}.spm.tools.LI_cfg.pre         = 0;
        matlabbatch{2}.spm.tools.LI_cfg.op          = 4;
        matlabbatch{2}.spm.tools.LI_cfg.vc          = 0;
        matlabbatch{2}.spm.tools.LI_cfg.ni          = 1;
        matlabbatch{2}.spm.tools.LI_cfg.outfile     = ['li_boot_con' num2str(con) all(s).name '.txt'];
        
        % check LI curves
        if con == 2
            matlabbatch{3}.spm.tools.LI_cfg.spmT        = {[pwd filesep all(s).name filesep 'stats' filesep 'spmT_000' num2str(con) '.nii']};
            matlabbatch{3}.spm.tools.LI_cfg.inmask.im11 = fname;
            matlabbatch{3}.spm.tools.LI_cfg.exmask.em1  = 1;
            matlabbatch{3}.spm.tools.LI_cfg.method.thr5 = 1;
            matlabbatch{3}.spm.tools.LI_cfg.pre         = 0;
            matlabbatch{3}.spm.tools.LI_cfg.op          = 4;
            matlabbatch{3}.spm.tools.LI_cfg.vc          = 0;
            matlabbatch{3}.spm.tools.LI_cfg.ni          = 1;
            matlabbatch{3}.spm.tools.LI_cfg.outfile     = ['li_curvefus_con' num2str(con) all(s).name '.txt'];
        end
        
        spm_jobman('run',matlabbatch);
        close all
        
        % assemble the information from text files
        if con == 2
            out = importdata(['li_boot_con' num2str(con) all(s).name '.txt'],'\t');
            li_all_brain(s) = str2num(cell2mat(out.textdata(2,5)));
            out = importdata(['li_bootfus_con' num2str(con) all(s).name '.txt'],'\t');
            li_fusiform(s) = str2num(cell2mat(out.textdata(2,5)));
            out = importdata(['li_curvefus_con' num2str(con) all(s).name '.txt'],'\t');
            li_fus_curves{s} = cell2mat(cellfun(@str2num,out.textdata(:,5),'un',0)); % not the values just the steps
        else
            out = importdata(['li_boot_con' num2str(con) all(s).name '.txt'],'\t');
            li_stim_all_brain(s,con-2) = str2num(cell2mat(out.textdata(2,5)));
            out = importdata(['li_bootfus_con' num2str(con) all(s).name '.txt'],'\t');
            li_stim_fusiform(s,con-2) = str2num(cell2mat(out.textdata(2,5)));
        end
    end
    
    delete LI_curves.ps
    delete LI_masking.ps
    delete LI_boot.ps
    delete (['LI_r_spmT_000' num2str(con) '.nii'])
end

figure;
lat_curve_mat = NaN(max(cellfun(@length,li_fus_curves)),15);
for s=1:size(all,1)
    nsteps(s) = length(unique(li_fus_curves{s}));
    out = importdata(['li_curvefus_' all(s).name '.txt'],'\t');
    lat_curve_mat(1:nsteps(s),s) = cell2mat(cellfun(@str2num,out.textdata(1:nsteps(s),6),'un',0));
end
lat_curve_mat(max(nsteps)+1:end,:) = [];
plot(lat_curve_mat,'LineWidth',3);
title('Fusiform lateralization curves'); grid on, box on;
ylabel('Laterality (right-left)'); xlabel('T values');

figure; 
subplot(2,3,1); rst_boxplot(li_stim_fusiform); title('Fusiform LI')
xticklabels({'Famous','Unfamiliar','Scrambled'})
subplot(2,3,2); plot(li_stim_fusiform(:,1),li_stim_fusiform(:,3),'o','LineWidth',2)
axis square; axis([-1 1 -1 1]); grid on; box on; hold on; plot([-1 1],[-1 1],'k','LineWidth',2)
plot(li_stim_fusiform(:,2),li_stim_fusiform(:,3),'o','LineWidth',2); title('Faces vs Scrambled')
xlabel('faces'); ylabel('scrambled faces')
subplot(2,3,3); plot(li_stim_fusiform(:,1),li_stim_fusiform(:,2),'o','LineWidth',2)
axis square; axis([-1 1 -1 1]); grid on; box on; hold on; plot([-1 1],[-1 1],'k','LineWidth',2)
title('Faces only'); xlabel('famous'); ylabel('Unfamiliar')

subplot(2,3,4); rst_boxplot(li_stim_all_brain); title('All brain LI')
xticklabels({'Famous','Unfamiliar','Scrambled'})
subplot(2,3,5); plot(li_stim_all_brain(:,1),li_stim_all_brain(:,3),'o','LineWidth',2)
axis square; axis([-1 1 -1 1]); grid on; box on; hold on; plot([-1 1],[-1 1],'k','LineWidth',2)
plot(li_stim_all_brain(:,2),li_stim_all_brain(:,3),'o','LineWidth',2); title('Faces vs Scrambled')
xlabel('faces'); ylabel('scrambled faces')
subplot(2,3,6); plot(li_stim_all_brain(:,1),li_stim_all_brain(:,2),'o','LineWidth',2)
axis square; axis([-1 1 -1 1]); grid on; box on; hold on; plot([-1 1],[-1 1],'k','LineWidth',2)
title('Faces only'); xlabel('famous'); ylabel('Unfamiliar')

% plot con2 (faces>scrambled)
[mean_leateralization,ci_lat]=rst_data_plot([li_all_brain' li_fusiform'],'estimator','mean');
title('Lateralization full brain and fusiform')

csvwrite('LI.csv',[li_fusiform' li_stim_fusiform li_all_brain' li_stim_all_brain]);
csvwrite('lat_curves.csv',lat_curve_mat')

save all 

