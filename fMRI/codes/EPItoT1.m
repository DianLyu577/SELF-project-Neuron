function EPItoT1(ref_img, source_img, others, is4D)
% a coregistration tool
% epi coregistered to t1, without reslicing to t1 resolution

if nargin < 4
    is4D = [];
end

if ischar(others)
    others = cellstr(others);
end

[FuncDir, epi, ext]  = fileparts(others);
FuncDir = cellstr(FuncDir);
epi = cellstr([char(epi) char(ext)]);

disp('**********************************************')
disp(['Coregistration to ' ref_img ]);

%% step 1: coregister mean fMRI to T1
if isempty(is4D)
    V = spm_vol(fullfile(FuncDir{1}, epi{1}));
    if length(V)>1 % 4d file
        is4D  = 1;
    else
        is4D  = 0;
    end
end
%% ---------------------
clear files
if is4D == 1
    files = {};
    for ie = 1:length(others) % if there are multiple sessions, loop through them
        eachfile = others{ie};
        files = [files; cellstr(spm_select('expand', eachfile))]; %#ok<AGROW>
    end
else
    files = others;
end

cd(FuncDir{1})
%% -------------------- Coregistration  ------------------------------------
clear matlabbatch
matlabbatch{1}.spm.spatial.coreg.estimate.ref  = {ref_img};
matlabbatch{1}.spm.spatial.coreg.estimate.source = {source_img};
matlabbatch{1}.spm.spatial.coreg.estimate.other = files;
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [4 4];
%     % run batch
disp('Step 1: to coregister the mean functional image to T1...');
spm_jobman('run',matlabbatch)


%% realign to the coregistered file
%     %% step 2: reslice to the source_Img
%     V = spm_vol(source_img);
%     imat = spm_imatrix(V(1).mat);
%     vsizemm = abs(imat(7:9));
%     %--------------------------------------
%     resize_img(vsizemm, others); % prefix is r (by default)
cd(FuncDir{1})
%% step 2: reslice to the reference image
%     P = [{source_img}; files];
%     flag.warp = [1 1 0];
%     flag.prefix = 'c';
%     flag.mean = 0;
clear matlabbatch
matlabbatch{1}.spm.spatial.realign.estwrite.data = {[{source_img}; files]};
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [1 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'c';
%     % run batch
disp('Step 2: to realign and reslice other functional images to the coregistered mean...');
spm_jobman('run',matlabbatch)
%--------------------
disp('**********************************************')
disp('Job completed.');
disp('   ')
disp('done')
end
