
clear
%--------
%% Sort Directories
%--------
dropbox_dir = '/data/dian';
dir_base = '/data/dian/Working_projects/neuroImg';
spm_path = fullfile(dropbox_dir, 'Dropbox/scripts/external/spm12'); % change here
ArtRepair_path = fullfile(dropbox_dir,'Dropbox/scripts/external/ArtRepair');
code_path = fullfile(dropbox_dir,'Dropbox/scripts/Stanford/fMRI_pipeline');
%colormap_path = '~/scripts/matlab/ColorBrewer';

sublist   = dir(fullfile(dir_base, 'SHICEP*'));
sublist   = {sublist.name}';
sublist   = sublist([1 3 7 8 9]);

stat      = 'spmT';

%% Initiating
%addpath(colormap_path)
addpath(genpath(spm_path))
addpath(genpath(code_path))
addpath(ArtRepair_path)
%--------------------------
spm('defaults','fmri');
spm_jobman('initcfg');

%% what analyses
glm1_spec  = 1;
glm1_est   = 1;
glm1_contr = 1;
smoothing  = 1;
normalise  = 1;
glm2_spec  = 1;
glm2_est   = 1;
glm2_contr = 1;

%% ============================== for sub ==================================
for ss = 1:length(sublist)

    clear matlabbatch
    sub = sublist{ss};

    stat_folder = fullfile(dir_base, sub, 'MRI/func/model/GLM1_SelfHot_Contrast');

    if ~exist(stat_folder, 'dir'); mkdir(stat_folder); end

    hotFC  = dir(fullfile(dir_base, sub, ['MRI/func/model/GLM1_ArtRepaired_FC*hot_sphere/' stat '_0001.nii']));
    coldFC = dir(fullfile(dir_base, sub, ['MRI/func/model/GLM1_ArtRepaired_FC*cold_sphere/' stat '_0001.nii']));
  
    if length(hotFC) <= 1
        fprintf('The sample size for %s hot electrodes is %d... \n', sub, length(hotFC))
        continue;
    end
    
    hotFiles  = strrep(cellstr([char({hotFC.folder}), repmat('/', length(hotFC), 1), char({hotFC.name})]),' ','');
    coldFiles = strrep(cellstr([char({coldFC.folder}), repmat('/', length(coldFC), 1), char({coldFC.name})]),' ','');

    stat_img = [hotFiles; coldFiles];
    
    %% smoothing
    if smoothing
        [path, img, ext] = fileparts(stat_img);
        %normed_stat_img = strrep(cellstr([char(path), repmat('/m',length(stat_img),1), char(img),  char(ext)]), ' ', '');
        smooth_img(stat_img, 6);
    end

    %% normalise
    if normalise
        smoothed_stat_img = strrep(cellstr([char(path), repmat('/s',length(stat_img),1), char(img),  char(ext)]), ' ', '');
        if strcmp(sub, 'SHICEP_S22_183_CR')
            yFile = fullfile(dir_base, sub, 'MRI/anat/y_ravg_axial_mprage6_2.nii');
        else
        yFile = fullfile(dir_base, sub, 'MRI/anat/y_T1.nii');
        normalise_writeImg(yFile, smoothed_stat_img, 'm');
        end
    end

    %% glm1 specification
    if glm1_spec == 1
      %  if exist(fullfile(stat_folder, 'SPM.mat'), 'File')~=2
            hotFiles  = strrep(cellstr([char({hotFC.folder}), repmat('/', length(hotFC), 1), char({hotFC.name})]),' ','');
            coldFiles = strrep(cellstr([char({coldFC.folder}), repmat('/', length(coldFC), 1), char({coldFC.name})]),' ','');

            matlabbatch{1}.spm.stats.factorial_design.dir = {stat_folder};
            matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = hotFiles;
            matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = coldFiles;
            matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 1;
            matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
            matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
            matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
            matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
            matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
            %-------------------------------
            spm_jobman('run', matlabbatch);
    end

    %% ============== glm1 estimation ================
    if glm1_est == 1
        clear matlabbatch
        % Model estimation
        disp(['Estimate GLM1 for Subject: ', sub]);
        matlabbatch{1}.spm.stats.fmri_est.spmmat = {[stat_folder,'/SPM.mat']}; %Select the SPM.mat
        matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0; % take value 0 or 1
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
        % run batch
        disp('Begin to run...');
        spm_jobman('run', matlabbatch);
        disp('GLM1 estimation job completed.');
        disp('   ')
    end

    %% ============== glm1 contrast ==============
    if glm1_contr
        clear matlabbatch
        matlabbatch{1}.spm.stats.con.spmmat = {[stat_folder,'/SPM.mat']};
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'hot-cold';
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

        matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'hot';
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = 1;
        matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep =  'none'; % replicate over sessions

        matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'cold';
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [0 1];
        matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep =  'none'; % replicate over sessions

        matlabbatch{1}.spm.stats.con.delete = 1;
        % run batch
        disp('Begin to run...');
        spm_jobman('run', matlabbatch);
        disp('GLM1 estimation job completed.');
        disp('   ')
    end
    %% ==============================================================

    stat_imgs1 = dir(fullfile(stat_folder, [stat '_00*.nii']));
    stat_imgs1 = cellstr([char({stat_imgs1.folder}), repmat('/', length(stat_imgs1), 1) ,char({stat_imgs1.name})]);

    stat_imgs2 = dir(fullfile(stat_folder, [stat '_00*.nii']));
    stat_imgs2 = cellstr([char({stat_imgs2.folder}), repmat('/', length(stat_imgs2), 1) ,char({stat_imgs2.name})]);

    stat_img   = [stat_imgs1; stat_imgs2];

 

end

%% glm2
for c = 1:3
    t = sprintf('ms%s_%04d.nii', stat, c);
    if c==1
        glm2_folder = fullfile(dir_base, 'RESULTS', 'GLM2_selfContr');
    elseif c == 2
        glm2_folder = fullfile(dir_base, 'RESULTS', 'GLM2_selfHot');
    elseif c == 3
        glm2_folder = fullfile(dir_base, 'RESULTS', 'GLM2_selfCold');
    end
    %%
    if glm2_spec

        if ~exist(glm2_folder, 'dir'); mkdir(glm2_folder); end

        %--------------------------
        stats_maps = {};
        for ss = 1:length(sublist)
            sub = sublist{ss};
            stats_map = fullfile(dir_base, sub,'MRI/func/model/GLM1_SelfHot_Contrast', t);
            if exist(stats_map, 'File') ~= 2; continue;
            else
                stats_maps{end+1,1} = stats_map;
            end
        end
        clear matlabbatch
        matlabbatch{1}.spm.stats.factorial_design.dir = {glm2_folder};
        matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = stats_maps;
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {fullfile(spm_path, 'tpm/TPM_SPM_brain.nii')};
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

        % run batch
        disp('Begin to run...');
        spm_jobman('run', matlabbatch);
        disp('GLM2 specification job completed.');
        disp('   ')
    end

    %%  glm2 estimation
    if glm2_est
        clear matlabbatch
        % Model estimation
        disp(['Estimate GLM2 for contrast: ', glm2_folder]);
        matlabbatch{1}.spm.stats.fmri_est.spmmat = {[glm2_folder,'/SPM.mat']}; %Select the SPM.mat
        matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0; % take value 0 or 1
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
        % run batch
        disp('Begin to run...');
        spm_jobman('run', matlabbatch);
        disp('GLM2 estimation job completed.');
        disp('   ')

    end

    %% glm2 contrast
    if glm2_contr
        clear matlabbatch
        matlabbatch{1}.spm.stats.con.spmmat = {[glm2_folder,'/SPM.mat']};
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'up';
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1 ;
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.delete = 1;

        % run batch
        disp('Begin to run...');
        spm_jobman('run', matlabbatch);
        disp('GLM2 contrast job completed.');
        disp('   ')
    end
end