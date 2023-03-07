% This script will calculate the FC of the self-hot and self-code
% electrodes, and compare between them.
n = 1;
% home_dir = char(expanduser('~'));
home_dir = '/data/dian';
spm_path = fullfile(home_dir, 'Dropbox/scripts/external/spm12'); % change here
code_path = fullfile(home_dir,'Dropbox/scripts/Stanford/fMRI_pipeline');

addpath(fullfile(home_dir, 'Dropbox/scripts/external/ArtRepair'));
addpath(genpath(code_path));
addpath(genpath(spm_path));
%colormap_path = '~/scripts/matlab/ColorBrewer';
ArtRepair_path = fullfile(home_dir,'Dropbox/scripts/external/ArtRepair');
result_folder = fullfile(home_dir,'Dropbox/Stanford_Matters/data/Stanford/SELF');
T = readtable(fullfile(home_dir, 'Dropbox/Stanford_Matters/data/cohort_data/SelfProj_allStims_fullCohort.csv'));

self_cohort = {'S19_137_AF', 'S21_166_TM', 'S22_176_LB','S20_152_HT','S21_165_WN', 'S22_178_AF', 'S22_183_CR'};
T = T(ismember(T.sbj_name, self_cohort) & strcmp(T.JP_label_std, 'PMC'),:);
schannels = strcat(T.sbj_name, '_', T.FS_label);

%% contruct a seed structure
    seeds = struct();

    [S,iCh,~] = unique(schannels);
    T = T(iCh,:);

    for s = 1:length(S)

        if T.Hot_Cold_Self(s) == 1
            id = 'hot';
        elseif T.Hot_Cold_Self(s) == 0
            id = 'cold';
        end
        seeds(s).name   = [S{s} '-' id '_sphere'];
        seeds(s).coord  = [T.MNI_coord_1(s) T.MNI_coord_2(s) T.MNI_coord_3(s)];
        seeds(s).radius = 4;
    end
%% ============== prepare sublist batches ======================
%---------------------------
dir_bases = {};
subList = {};
catgs = {'filtered_subs', 'problem_subs'};
for i = 1:2
    catg = catgs{i};
    List   = dir(fullfile('/data/dian/Working_projects/HCP/rest_HCP', catg));
    Name   = {List.name};
    matchC = regexp(Name, '^\d+', 'start') ;  % Or: '^[0-9]+\.txt'
    match  = ~cellfun('isempty', matchC);
    List   = List(match);
    
    dir_bases = [dir_bases; {List.folder}'];
    [subList,I] = sort([subList; {List.name}']);
    dir_bases = dir_bases(I);
end
%---------------------------
% nSub   = 12;
% Nbatch_floor = floor(length(subList)/nSub);
% if length(subList)/nSub - Nbatch_floor <= 0.25
%     Nbatch = Nbatch_floor;
% else
%     Nbatch = Nbatch_floor+1;
% end
% 
% if n<Nbatch
%     idlist = (1:nSub) + 12*(n-1) % n is the looped number of batch
% elseif n==Nbatch
%     idlist = (nSub + 12*(n-2) +1) : length(subList)
% else
%     disp(['n cannot be bigger than ' num2str(Nbatch) '. --> Return'])
%     return
% end
idlist = 1:length(subList);
%%=========================
sublist = subList(idlist);
dir_bases = dir_bases(idlist);

%---------------------------
spm('defaults','fmri');
spm_jobman('initcfg');

%% ======================================
%% Loop through a batch of subject lists

for ss = 33:length(sublist)
    tic
    dir_base = dir_bases{ss};
     sub = sublist{ss};
     if ismember(sub, {'103111', '154734','192540', '756055'})
         continue;
     end
     func = 'MNINonLinear/Results/rfMRI_REST2_LR';
    % --- subfolders
    FuncDir = fullfile(dir_base, sub, func);
    fprintf ('\n Processing for %s...' , sub)
    % for debug
    %load(fullfile(FuncDir, 'preOP_pipline.mat'))

    %% class initiation (assigning properties)
    fpip                = fMRINetMapNative;
    fpip.TR             = 0.72;
    fpip.spm_path       = spm_path;
    fpip.ArtRepair_path = ArtRepair_path;
    fpip.noise          = 'none';
    fpip.sub            = sub;
    fpip.FuncDir        = FuncDir;

    denoised_epi   = 'rfMRI_REST2_LR_hp2000_clean_smooth.nii';
    fpip.denoised_epi   = cellstr(denoised_epi);

    if ~exist(fullfile(FuncDir,denoised_epi), 'file')
        if exist(fullfile(FuncDir,[denoised_epi '.gz']), 'file')
            gunzip(fullfile(FuncDir, denoised_epi), FuncDir)
        else
            warning(['Denoised file is not found for ' sub '...'])
            continue
        end
    end

    %% apply function
     fpip.M = struct(); % use default, denoise without GSR
    fpip.M.GLM1_dir = fullfile(fpip.FuncDir, 'model', 'GLM1');
    fpip.M.main_vars = [];
    
    fpip.M.multiregressorTXT = 'allconfounds.txt';
    if ~exist(fullfile(fpip.FuncDir, 'allconfounds.txt'), 'file')
        fpip.M.multiregressorTXT = 'rfMRI_REST2_CSF_WM.txt';
    end
    
    fpip.glm1(1,0)
    fpip.glm1est(1)
    fpip.glm1cont(1)
    % %% seed-based FC
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
 
        main_vars = cell(length(fpip.denoised_epi),1);
        
        for sess = 1:length(fpip.denoised_epi)
            to_skip = 0;
            if exist(fullfile(fpip.M(1).GLM1_dir, ['VOI_' seed '_' num2str(sess) '.mat']), 'file')==2
                d = load(fullfile(fpip.M(1).GLM1_dir, ['VOI_' seed '_' num2str(sess) '.mat']));
                VAR.vec = d.Y;
                VAR.name = [seed '_PC'];
                main_vars{sess} = VAR;
            elseif exist(fullfile(fpip.M(1).GLM1_dir, [seed '_PC_components_' num2str(sess) '.mat']), 'file')==2
                % use the first 6 PCs of the ROI as the main factors
                d = load(fullfile(fpip.M(1).GLM1_dir, [seed '_PC_components_' num2str(sess) '.mat']));
                main_vars{sess} = d.PC_components;
                
            else
                disp(['Vector file for ' seed ' is not found, skipping...'])
                to_skip = 1;
                continue
            end
            %         T = array2table(d.PC_components,...
            %             'VariableNames', cellstr([repmat([ROIs(i).name '-PC'],6,1) num2str((1:6)')]));

        end
          
        if to_skip == 1
            continue;
        end
        M_.main_vars = main_vars;
        M_.multiregressorTXT = fpip.M(1).multiregressorTXT;
        fpip.M(end+1) = M_;
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
   % fpip.glm1result(1, 2:length(fpip.M), 2, 1)
    toc
end