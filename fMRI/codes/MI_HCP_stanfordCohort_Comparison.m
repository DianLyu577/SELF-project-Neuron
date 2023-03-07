clear
%% Sort Directories
%--------
dropbox_dir = '/data/dian';
%dir_base = '/data/dian/Working_projects/neuroImg';
%dropbox_dir = 'C:\Users\lvdia';
%dir_base = 'E:\Stanford\tmp';
brainspace_path = fullfile(dropbox_dir,'Dropbox/scripts/external/brainspace-0.1.10/brainspace/matlab');
gifti_lib = fullfile(dropbox_dir,'Dropbox/scripts/external/matlab_GIfTI-master');
code_path = fullfile(dropbox_dir,'Dropbox/scripts/Stanford/fMRI_pipeline');

addpath(genpath(brainspace_path))
addpath(genpath(gifti_lib))
addpath(genpath('/data/dian/Dropbox/scripts/external/surfstat'))

%% ---- FC Similarity Test ----
%=======================================================
HCP_stat_folder = '/data/dian/Working_projects/HCP/rest_HCP/RESULTS/GLM2_selfHot_conFiles';
SELF_stat_folder = '/data/dian/Working_projects/neuroImg/RESULTS/GLM2_selfHot';

surfsboth = {'/data/dian/Working_projects/ROIs/S900.L.midthickness_MSMAll.32k_fs_LR.surf.gii';...
    '/data/dian/Working_projects/ROIs/S900.R.midthickness_MSMAll.32k_fs_LR.surf.gii'};
surf_out = convert_surface(surfsboth,'format','MATLAB');

%----------------- creating null distribution ---------------------------
featureL= gifti(fullfile(HCP_stat_folder, 'con_0001inL.shape.gii'));
featureR= gifti(fullfile(HCP_stat_folder, 'con_0001inR.shape.gii'));

plot_hemispheres([featureL.cdata; featureR.cdata], surfsboth);

n_perm = 1000;
% spin permutation
%features_spin = spin_permutations({featureL.cdata, featureR.cdata},...
%    surfsboth, n_perm);
%save(fullfile(HCP_stat_folder, 'nulldist.mat'), 'features_spin')
load(fullfile(HCP_stat_folder, 'nulldist.mat'));

%% plot sphere (demonstration)
FEATURE = [featureL.cdata; featureR.cdata];
conn = reshape(surf_out{1}.faces, [], 1);
textureL = reshape(featureL.cdata(conn), size(surf_out{1}.faces));
conn = reshape(surf_out{2}.faces, [], 1);
textureR = reshape(featureR.cdata(conn), size(surf_out{2}.faces));
texture = [textureL textureR];

surfs_spheres = {'~/scripts/external/brainspace-0.1.10/brainspace/matlab/datasets/surfaces/conte69_32k_left_sphere.gii';...
    '~/scripts/external/brainspace-0.1.10/brainspace/matlab/datasets/surfaces/conte69_32k_right_sphere.gii'};
% impute NaN
shifts = 60;
T_imp = zeros(size(texture,1),(shifts*2+1).*size(texture,2));
for sft = -shifts:shifts
    if sft == 0
        T_imp(:,(sft+shifts+1):(sft+shifts)+size(texture,2)) = texture;
    else % rotate 
        tail = flip((sft+shifts+1):size(T_imp,1));
        head = 1:(sft+shifts);
        timp = [texture(tail,:); texture(head,:)];
        T_imp(:,(sft+shifts+1):((sft+shifts)+size(texture,2))) = timp; 
    end
end
Texture_imp = knnimpute(T_imp); % impute nan based on two hemispheres and neighbours
textureL_imp = Texture_imp(:,shifts+1:shifts+3);
textureR_imp = Texture_imp(:,shifts+4:shifts+6);

featureL_imp = nan(size(featureL.cdata,1),1);
featureR_imp = nan(size(featureL.cdata,1),1);
for id = 1:size(featureL.cdata,1)
    featureL_imp(id) = mean(textureL_imp(surf_out{1}.faces==id));
    featureR_imp(id) = mean(textureR_imp(surf_out{2}.faces==id));
end

feature_imp = [featureL_imp; featureR_imp];
feature_imp(~isnan(FEATURE)) = FEATURE(~isnan(FEATURE));

features_spin_imp = spin_permutations({feature_imp(1:length(featureL_imp)),...
   feature_imp(1+length(featureL_imp) : end) },...
    surfsboth, n_perm);

% plot sphere
close all
nspin = 669;
plot_hemispheres(mean([feature_imp(1:length(featureL_imp)),...
   feature_imp(1+length(featureL_imp):end)],2), surfs_spheres(1));caxis([-0.015,0.19])
plot_hemispheres(mean([features_spin_imp{1}(:,1,nspin) features_spin_imp{2}(:,1,nspin)],2), surfs_spheres(1) ); caxis([-0.02,0.12])
colormap(cmap.cm)
% plot brain
plot_hemispheres(mean([feature_imp(1:length(featureL_imp)),...
   feature_imp(1+length(featureL_imp):end)],2), surf_out(1));caxis([-0.02,0.15])
plot_hemispheres(mean([features_spin_imp{1}(:,1,nspin) features_spin_imp{2}(:,1,nspin)],2), surf_out(1) );caxis([-0.02,0.12])

%% -------------------Spearman Correlation---------------------------------
featurel= gifti(fullfile(SELF_stat_folder, 'con_0001inL.shape.gii'));
featurer= gifti(fullfile(SELF_stat_folder, 'con_0001inR.shape.gii'));
feature = [featurel.cdata; featurer.cdata];

plot_hemispheres(featurel.cdata, surf_out(1));
plot_hemispheres(featurel.cdata, surfs_spheres(1));
caxis([-0.5 0.8])

%--- stat ---
rho = corr(feature, [featureL.cdata; featureR.cdata], 'Type','Spearman', 'rows','complete');

RHO = nan(n_perm,1);
for i = 1:1000
    FEATURE = [features_spin{1}(:,:,i); features_spin{2}(:,:,i)];
    RHO(i) = corr(feature, FEATURE, 'Type','Spearman', 'rows','complete');
end

% --- t stat ---
[h,p,ci,stats] = ttest(RHO-rho);
%{
h =

  single

     1


p =

  single

     0


ci =

  1×2 single row vector

   -0.2631   -0.2474


stats = 

  struct with fields:

    tstat: -63.8317
       df: 999
       sd: 0.1265
%}

%% histogram
figure('Position',[1147        1153         390         130]);
h=histfit(RHO, 20, 'normal');
h(1).FaceColor = [1 1 1];
h(2).Color = [.2 .2 .2];
hold on
x1=xline(rho,'r--', {'Our cohort', 't = 63.83', 'p = 0.00'},'LineWidth',1);
box off
xlabel('spatial similarities (rho)')
ylabel('density')

fig_dir = '/data/dian/Dropbox/Stanford_Matters/data/SELF/Plots/Figure2';
saveas(gcf,fullfile(fig_dir, 'hist_spinPerm.fig'))
exportgraphics(gcf,fullfile(fig_dir, 'hist_spinPerm.png'),'Resolution',300)
%x1.LabelOrientation = 'horizontal';

%% -------------------mutual information---------------------------------
muInfo = mi(feature, [featureL.cdata; featureR.cdata]);

MI = nan(size(featureL));

for i = 1:1000
    FEATURE = [features_spin{1}(:,:,i); features_spin{2}(:,:,i)];
    MI(i) = mi(feature, FEATURE) ;
end
%% --- t stat ---
[h,p,ci,stats] = ttest(MI-muInfo);


%% --- FC Generalisability test ----
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
idlist = ~ismember(subList, {'103111', '154734','192540', '756055'});
%%=========================
sublist = subList(idlist);
dir_bases = dir_bases(idlist);

%% =========================
stat      = 'con';
dir_base  = '/data/dian/Working_projects/HCP/rest_HCP/filtered_subs';
glm2_dir_base = '/data/dian/Working_projects/HCP/rest_HCP';

%% ============================== prepare images ==================================
%% sample pool: individual-level FC file of the HCP cohort
%--------------------------
hcpFC = [];
for ss = 1:length(sublist)
    sub = sublist{ss};
    dir_base = dir_bases{ss};

    if exist(fullfile(dir_base, sub, 'MRI/model/GLM1_SelfHot_Contrast_conFiles', 'con_0002inL.shape.gii'), 'File') ~= 2
        continue;
    end

    fc_L = gifti(fullfile(dir_base, sub, 'MRI/model/GLM1_SelfHot_Contrast_conFiles', 'con_0002inL.shape.gii'));
    fc_R = gifti(fullfile(dir_base, sub, 'MRI/model/GLM1_SelfHot_Contrast_conFiles', 'con_0002inR.shape.gii'));
    hcpFC = [hcpFC, [fc_L.cdata; fc_R.cdata]];

end

%% sample pool2: individual level FC file fo the SELF cohort
self_cohort = {'SHICEP_S19_137_AF';...
    'SHICEP_S21_166_TM';...
    'SHICEP_S22_176_LB';...
    'SHICEP_S22_178_AF';...
    'SHICEP_S22_183_CR'};

selfFC = [];
for ss = 1:length(self_cohort)
    sub = self_cohort{ss};
    dir_base = '/data/dian/Working_projects/neuroImg';

    if exist(fullfile(dir_base, sub, 'MRI/func/model/GLM1_SelfHot_Contrast', 'mscon_0002inL.shape.gii'), 'File') ~= 2
        continue;
    end

    fc_L = gifti(fullfile(dir_base, sub, 'MRI/func/model/GLM1_SelfHot_Contrast', 'mscon_0002inL.shape.gii'));
    fc_R = gifti(fullfile(dir_base, sub, 'MRI/func/model/GLM1_SelfHot_Contrast', 'mscon_0002inR.shape.gii'));
    selfFC = [selfFC, [fc_L.cdata; fc_R.cdata]];

end

plot_hemispheres(mean(reshape(hcpFC,[],2, size(hcpFC,2)), [2 3], 'omitnan'), surf_out(1)); caxis([-0.02,0.12])
plot_hemispheres(mean(reshape(selfFC,[],2, size(selfFC,2)), [2 3], 'omitnan'), surf_out(1));caxis([-0.5,0.7])
plot_hemispheres(mean(reshape(selfFC,[],2, size(selfFC,2)), [2 3], 'omitnan'), surfs_spheres(1));caxis([-0.5,0.6])

% plot self cohort
plot_hemispheres(mean(selfFC,2, 'omitnan'), surf_out);
cmap = load('~/scripts/PlottingTools/ColourMap_matlab/red_blue_stanford.mat');
colormap(cmap.cm)

%% Monte Carlo sampling
n_perm = 1000;
ind1 = 1:size(hcpFC, 2);
pair1 = ind1(randi(length(ind1), 2, n_perm))';
pair1(pair1(:, 1) == pair1(:,2),:) = []; % remove repetitive draw
pair1 = unique(pair1,'rows'); % index of [hcpFC, hcpFC]

RHO1 = nan(size(pair1,1),1);
for i = 1:size(RHO1,1)
    RHO1(i) = corr(hcpFC(:,pair1(i,1)), hcpFC(:,pair1(i,2)), 'Type','Spearman', 'rows','complete');
end


%% one-sample t test
[h,p,ci,stats] = ttest(RHO1-rho);

%{
h =

  single

     1


p =

  single

     0


ci =

  2×1 single column vector

   -0.0333
   -0.0233


stats = 

  struct with fields:

    tstat: -11.2142
       df: 932
       sd: 0.0771
%}



%% histogram
figure('Position',[1147        1153         390         130]);
h=histfit(RHO1, 20, 'normal');
h(1).FaceColor = [0.9 0.9 0.9];
h(1).EdgeColor = [1 1 1];
h(2).Color = [.75 .75 .75];
hold on
xline(rho,'r--', {'Our cohort'},'LineWidth',0.7);
xline([mean(RHO1) + std(RHO1), mean(RHO1) - std(RHO1) ],...
    '-', {'+1 STD','-1 STD'},'LineWidth',0.5);
box off
xlabel('individual differences (rho)')
ylabel('density')

fig_dir = '/data/dian/Dropbox/Stanford_Matters/data/SELF/Plots/Figure2';
saveas(gcf,fullfile(fig_dir, 'hist_indDiff_motecarlo.fig'))
exportgraphics(gcf,fullfile(fig_dir, 'hist_indDiff_motecarlo.png'),'Resolution',300)



