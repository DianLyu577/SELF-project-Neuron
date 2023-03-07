% This script will calculate the FC of the self-hot and self-code
% electrodes, and compare between them.
% for S19_137
clear
home_dir = getuserdir;
dir_base = '/media/lvdian/storage/data/Stanford_S19137/';
spm_path = fullfile(home_dir, 'Dropbox/scripts/external/spm12'); % change here
code_path = fullfile(home_dir,'Dropbox/scripts/Stanford/fMRI_pipeline');
% add paths
addpath(fullfile(home_dir, 'Dropbox/scripts/external/ArtRepair'));
addpath(genpath(spm_path))
addpath(genpath(code_path));
%colormap_path = '~/scripts/matlab/ColorBrewer';
ArtRepair_path = fullfile(home_dir,'Dropbox/scripts/external/ArtRepair');
result_folder = fullfile(home_dir,'Dropbox/Stanford_Matters/data/Stanford/SELF');
T = readtable(fullfile(home_dir, 'Dropbox/Stanford_Matters/data/cohort_data/table_EBS_PMC_fullcohort_anatomical_info.csv')); 
sublist = dir(fullfile(dir_base, 'SHICEP_*'));
sublist = {sublist.name}';
sublist = sublist(1);

%---------------------------
spm('defaults','fmri');
spm_jobman('initcfg');

%% contruct a seed structure

for ss = 1:length(sublist)

    seeds = struct();

    sub = sublist{ss};
    sub_keys = strsplit(sub, '_');
    sub_short = strjoin(sub_keys(2:3),'_');
    Ind = find(contains(T.SUBJECT, sub_short));
    schannels = T.SCHANNEL(Ind);
    [S,iCh,~] = unique(schannels);

    for s = 1:length(S)
        I = Ind(iCh(s));

        if T.Self_Hot(I) == 1
            id = 'hot';
        elseif T.Self_Hot(I) == -1
            id = 'cold';
        end
        seeds(s).name   = [S{s} '-' id '_sphere'];
        seeds(s).coord  = [T.Native_Coord_1(I) T.Native_Coord_2(I) T.Native_Coord_3(I)];
        seeds(s).radius = 4;
    end
    %=======================================
    struc = 'anat';
    func = 'func';
    T1_img = 'T1.nii';
    % --- subfolders
    StrucDir = fullfile(dir_base,  sub, 'MRI', struc);
    FuncDir = fullfile(dir_base, sub, 'MRI', func);
    GLM0_dir = fullfile(FuncDir, 'model', 'GLM0');
    % for debug
    %load(fullfile(FuncDir, 'preOP_pipline.mat'))

    %% class initiation (assigning properties)
    fpip = fMRINetMapNative;
    fpip.spm_path = spm_path;
    fpip.ArtRepair_path = ArtRepair_path;
    fpip.sub = sub;
    fpip.T1_img = T1_img;
    fpip.FuncDir = FuncDir;
    fpip.StrucDir = StrucDir;
    func_files = dir(fullfile(fpip.FuncDir, 'sub-S19_137*.nii'));
    fpip.epi = {func_files.name}';
    fpip.GLM0_dir = GLM0_dir ;
    fpip.denoised_epi = cellstr([repmat('c', length(fpip.epi),1), char(fpip.epi)]);

    %% apply function
    fpip.T1_segmentation(1) % HEAVY
    fpip.skullstrip(1)
    fpip.realign(1)
    fpip.coregister(1) % (if_check_done, if_reslice, if_4D)
%  %---- dummy GLM0 for extracting noise signal ------
    fpip.glm0(1)
    fpip.glm0est(1)% HEAVY
    fpip.glm0cont(1)
    
    fpip.motionScrub(1)
    fpip.smooth(1,3)
    fpip.noiseSignExtr(1)
%     
    fpip.glm1(1,0)
    fpip.glm1est(1)
    fpip.glm1cont(1)
%     %% seed-based FC
    n_comp = 1; 
    fpip.roiSignExtr(1, seeds, 1:length(seeds), 1)
% 
    %% >>>>>>>>> step 3: GLM using seed signal as the predictor/independent variable
    %% Functional Network Mapping
    % ---->>>>>>  prepare GLM1
    fpip.M = fpip.M(1); % always keep the first default GLM1
    % rewrite the new GLM1 model specification
    % seed is defined by the first few PC of the sessional data
    for i = 1:length(seeds)
        M_ = struct();
        M_.GLM1_dir = fullfile(fpip.FuncDir, 'model', strjoin({'GLM1', 'FC', seeds(i).name}, '_'));
        seed = seeds(i).name;

        VAR = cell(length(fpip.epi),1);

        for sess = 1:length(fpip.epi)
            if exist(fullfile(fpip.M(1).GLM1_dir, [seed '_PC_components_' num2str(sess) '.mat']), 'file')~=2
                d = load(fullfile(fpip.M(1).GLM1_dir, ['VOI_' seed '_' num2str(sess) '.mat']));
                VAR{sess} = d.Y;
            else
                % use the first 6 PCs of the ROI as the main factors
                d = load(fullfile(fpip.M(1).GLM1_dir, [seed '_PC_components_' num2str(sess) '.mat']));
                VAR{sess} = d.PC_components;
            end
            %         T = array2table(d.PC_components,...
            %             'VariableNames', cellstr([repmat([ROIs(i).name '-PC'],6,1) num2str((1:6)')]));

        end
        M_.main_vars = VAR;
        fpip.M(i+1) = M_;
    end
    %% %%%%%%%%%%%%%%%%%%%%% stats, heavy stuff %%%%%%%%%%%%%%%%%%%%%%%%
    for m = 2:length(fpip.M)
        fpip.glm1(1,0,m)
        fpip.glm1est(1,m)
    end

    PC_single = cellstr([repmat('PC',1,1) num2str((1:1)')]);
    CON = struct('name', PC_single', ...
        'weights', {1},...
        'sessrep', repmat({'replsc'}, 1, 1)); % SPM will fill up the zeros of the rest of the colums

    fpip.glm1cont(1, 2:length(fpip.M), CON)
    fpip.glm1result(1, 2:length(fpip.M), 2, 1)
end