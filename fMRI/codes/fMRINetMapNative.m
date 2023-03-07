classdef fMRINetMapNative < handle
    properties      % variables that can be accessed from the object; manage them in the same way as 'structure'.
        code_path
        spm_path        % the path of the spm12
        ArtRepair_path % the path of ArtRepair for motion correction
        sub % subject ID
        debug_batch % save the problematic batch outside the functions for debugging
        %====== structural stuff ======
        StrucDir       % structural scans folder
        T1_img         % the name of T1 image
        deformPar % the path to the deformation matrix (for the inverse normalisation from MNI to native space)
        %====== functional stuff ======
        TR
        smooth_fwhm
        FuncDir % functional imgaging folder
        GLM0_dir % dummy GLM folder for noise extraction
        movement_parameters
        noise
        M % first-level analysis models, a structure with GLM1_dir, and main vectors (for main effect)
        epi            % the name of the epi, functional image; can be a string vector or a cell
        exclusion % the scanning sessions to exclude because of excessive motions
        denoised_epi % denoised files to be used for analyses; can be a string vector or a cell
        % ====== ROI stuff ======
        %--> overlay_masks are in the native space
        overlay_masks %  a cell array with size of (n,1). They are multiple masks to overlay. The default procedure is to take an intersection among them;
        %--> ROIs are in the standard space
        ROIs % regions of interest
    end

    %% =========================================================
    methods

        %% %%%%%%%%%%%%%%%%%%% Preprocessing %%%%%%%%%%%%%%%%%%%%%%%
        %% Segmentation
        function T1_segmentation(obj, check)

            outputfile = fullfile(obj.StrucDir,['c5' obj.T1_img]);

            if check == 1 && nargin > 1
                if exist(outputfile, 'file') ~= 2
                    do_it = 1;
                else
                    do_it = 0;
                    disp('T1_segmentation is already done.')
                    disp(['Output is here: ' outputfile])
                end
            else
                do_it = 1;
            end

            if do_it == 1
                %% run matlabbatch
                disp('**********************************************')
                disp(['Processing structural file for ' obj.sub ])
                clear matlabbatch
                matlabbatch{1}.spm.spatial.preproc.channel.vols = {fullfile(obj.StrucDir, [obj.T1_img ',1'])} ;
                matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
                matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
                matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
                for i = 1:6
                    matlabbatch{1}.spm.spatial.preproc.tissue(i).tpm = cellstr(fullfile(obj.spm_path, 'tpm', ['TPM.nii,' num2str(i)]));
                    if i < 3
                        matlabbatch{1}.spm.spatial.preproc.tissue(i).ngaus = 1;
                    else
                        matlabbatch{1}.spm.spatial.preproc.tissue(i).ngaus = i-1;
                    end
                    matlabbatch{1}.spm.spatial.preproc.tissue(i).native = [1 0];
                    matlabbatch{1}.spm.spatial.preproc.tissue(i).warped = [0 0];
                    if i == 6
                        matlabbatch{1}.spm.spatial.preproc.tissue(i).ngaus = 2;
                        matlabbatch{1}.spm.spatial.preproc.tissue(i).native = [0 0];
                    end
                end
                matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
                matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
                matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
                matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
                matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
                matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
                matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];
                matlabbatch{1}.spm.spatial.preproc.warp.vox = NaN;
                matlabbatch{1}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
                    NaN NaN NaN];
                % run
                obj.debug_batch = matlabbatch;
                spm_jobman('run', matlabbatch);
                % display
                disp('Segment done.')
                obj.deformPar = fullfile(obj.StrucDir, ['y_' obj.T1_img]);
            end
        end % segmentation


        %% Skull stripping
        function skullstrip(obj, check)
            outputfile = fullfile(obj.StrucDir, ['brain_' obj.T1_img]);

            if check == 1 && nargin > 1
                if exist(outputfile, 'file') ~= 2
                    do_it = 1;
                else
                    do_it = 0;
                    disp('Skull stripping is already done.')
                    disp(['Output is here: ' outputfile])
                end
            else
                do_it = 1;
            end

            if do_it == 1
                clear matlabbatch

                % configure matlabbatch
                outputs = {'brain_mask'; 'brain'};
                expressions = {'((i1 +i2 +i3) > 0.3)&((i4 +i5) < 0.1)';...
                    'i6.*(((i1 +i2 +i3) > 0.3)&((i4 +i5) < 0.1))'};
                %
                input_files = {fullfile(obj.StrucDir,['c1'  obj.T1_img]);... % i1
                    fullfile(obj.StrucDir,['c2'  obj.T1_img]);...  % i2
                    fullfile(obj.StrucDir,['c3'  obj.T1_img]);...  % i3
                    fullfile(obj.StrucDir,['c4'  obj.T1_img]); ... % i4
                    fullfile(obj.StrucDir,['c5'  obj.T1_img])};... % i5
                    %----
                %
                for i = 1:length(outputs)
                    output = outputs{i};
                    expression = expressions{i};
                    if isequal(output, 'brain')
                        input_files{6,1} = fullfile(obj.StrucDir,['m' obj.T1_img, ',1']); %i6 m*.nii has better homogeneity in contrast
                    end
                    %----------------------------
                    matlabbatch{i}.spm.util.imcalc.input = input_files;
                    matlabbatch{i}.spm.util.imcalc.output = [output '_'  obj.T1_img];
                    matlabbatch{i}.spm.util.imcalc.outdir = {obj.StrucDir};
                    matlabbatch{i}.spm.util.imcalc.expression = expression; % skull stripping-> whole-brain mask (including CSF)
                    matlabbatch{i}.spm.util.imcalc.var = struct('name', {}, 'value', {});
                    matlabbatch{i}.spm.util.imcalc.options.dmtx = 0;
                    matlabbatch{i}.spm.util.imcalc.options.mask = 0;
                    matlabbatch{i}.spm.util.imcalc.options.interp = 1;
                    matlabbatch{i}.spm.util.imcalc.options.dtype = 4;
                    %---
                end

                % run batch
                disp('Begin to run...');
                spm_jobman('run',matlabbatch)
                disp('Job completed.');
                disp('   ')
            end
        end % skullstrip


        %% Warp ROIs from MNI to native space
        function MNI2native(obj, check, i_rois)
            disp('>>>>>>>>>>>>>>> Projecting MNI --> native >>>>>>>>>>>>>>>>>>>>>>>>>')
            % i_rois is the index/indices of the given ROIs
            % it may be useful for debugging; i.e. only do one out of all
            if isempty(obj.deformPar)
                % go to the individual structural folder to find the file
                yfile = fullfile(obj.StrucDir, ['y_' obj.T1_img]);
                if exist(yfile, 'file') == 2
                    obj.deformPar = yfile;
                else % apply the previous step
                    disp('No deformation parameter is found, creating now...')
                    obj.T1_segmentation(0)
                end
            end
            %-----------------------------
            if isempty(obj.ROIs)
                error('No ROI is given.')
            end

            if nargin < 3
                i_rois = 0;
            end

            %-------------------------------------
            if i_rois == 0
                roi_loopin = 1:length(obj.ROIs);
            else
                roi_loopin = i_rois;
            end
            %============loop though all rois================
            for ir = roi_loopin
                paths = obj.ROIs(ir).path;
                for n = 1:length(paths)
                    [~, roi, ext] = fileparts(paths{n}) ;
                    outputfile = fullfile(obj.StrucDir, ['w' roi ext]);
                    %--------------------------------
                    if check == 1 && nargin > 1
                        if exist(outputfile, 'file') ~= 2
                            do_it = 1;
                        else
                            do_it = 0;
                            disp(['MNI2native transformation is already done for ' roi '.'])
                            disp(['Output is here: ' obj.StrucDir])
                        end
                    else
                        do_it = 1;
                    end
                    %---------------------------------------------
                    if do_it == 1
                        cd(obj.StrucDir)
                        clear matlabbatch
                        matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = {obj.deformPar};
                        matlabbatch{1}.spm.util.defs.comp{1}.inv.space = {fullfile(obj.StrucDir, obj.T1_img)};
                        matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = obj.ROIs(ir).path;
                        matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = cellstr(obj.StrucDir);
                        matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 1;
                        matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
                        matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
                        matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'w';
                        % run batch
                        obj.debug_batch = matlabbatch;
                        spm_jobman('run',matlabbatch);
                        disp(['Standard ' roi ' is mapped to the individual space.']);

                    end
                end
            end
        end

        %% %%%%%%%%%%%   Functional Images Processing   %%%%%%%%%%%%%%%%%%%%%
        %% reorient images
        function reorient(obj, check, reorient_mat)
            % to correct for wrong origin; once manually set the origin in
            % anterior commisure, click set origin, and reorient to write
            % file, save the matrix for applying for the functional images
            % --- without reorient the functional images, the
            % coregistration might fail
            outputfile = fullfile(obj.FuncDir, ['o' obj.epi{end}]);
            if check == 1 && nargin > 1
                if exist(outputfile, 'file') ~= 2
                    do_it = 1;
                else
                    do_it = 0;
                    disp('Functional images are already reset origin.')
                end
            else
                do_it = 1;
            end

            if do_it == 1
                % select files
                files = {};
                for ie = 1:length(obj.epi) % if there are multiple sessions, loop through them
                    func_file = obj.epi{ie};
                    fs = cellstr(spm_select('expand', fullfile(obj.FuncDir, func_file)));
                    files = [files; fs];
                end% for epi
                %% run matlabbatch
                clear matlabbatch
                disp('Reorienting functional images based on the new origin....')
                matlabbatch{1}.spm.util.reorient.srcfiles = cellstr(files);
                %%
                m = load(fullfile(obj.StrucDir, reorient_mat));
                matlabbatch{1}.spm.util.reorient.transform.transM = m.M;
                matlabbatch{1}.spm.util.reorient.prefix = 'o';

                % run matlabbatch
                disp('Begin to run...');
                spm_jobman('run',matlabbatch);
                disp('Job completed.');
                disp('   ')
                %---- update fMRI image -----
                if obj.epi{1}(1) ~= 'o'
                    obj.epi = reshape(obj.epi, [],1);
                    obj.epi = cellstr([repmat('o', length(obj.epi),1) char(obj.epi)]);
                end
            end

        end

        %% Realignment
        % estimate and reslice
        % need to estimate to get the headmovement file (rp_restxx.txt ),
        % need to reslice to get the mean image for the later
        % coregistration.
        function realign(obj, check)

            if ischar(obj.epi)
                obj.epi = cellstr(obj.epi);
            elseif isempty(obj.epi) % assuming default file names
                fs = dir(fullfile(obj.FuncDir, 'rest*.nii'));
                obj.epi = {fs.name}';
            end

            %============================================
            outputfile = fullfile(obj.FuncDir, ['mean' obj.epi{1}]);

            if check == 1 && nargin > 1
                if exist(outputfile, 'file') ~= 2
                    do_it = 1;
                else
                    do_it = 0;
                    disp('Realignment is already done.')
                    disp(['Output is here: ' outputfile])
                end
            else
                do_it = 1;
            end

            if do_it == 1
                cd(obj.FuncDir)
                disp('Realign: Estimate & reslice...');
                disp(['Realignment for  ' obj.sub ]);
                % select files
                files = cell(1,length(obj.epi));
                for ie = 1:length(obj.epi) % if there are multiple sessions, loop through them
                    func_file = obj.epi{ie};
                    fprintf('Ready to process for %s, session: %d....\n', obj.sub, ie)
                    files{ie} = cellstr(spm_select('expand', fullfile(obj.FuncDir, func_file)));
                end% for epi
                %-------------------------------
                matlabbatch{1}.spm.spatial.realign.estwrite.data = files;
                % specify estimate options
                matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
                matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
                matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
                matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
                matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
                matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
                matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
                matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [0 1]; % estimate & reslice; reslice only the mean
                matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
                matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
                matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
                matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r'; % if not reslicing images here, there is no prefix
                % run batch
                disp('Begin to run...');
                spm_jobman('run',matlabbatch);
                disp('Job completed.');
                disp('   ')

            end % if do it
        end % realign

        %% Coregister
        function coregister(obj, check, reslice, is4D)

            if ischar(obj.epi)
                obj.epi = cellstr(obj.epi);
                %             elseif isempty(obj.epi) % assuming default file names
                %                 fs = dir(fullfile(obj.FuncDir, 'rest*.nii'));
                %                 obj.epi = {fs.name}';
            end

            if nargin < 4
                is4D = 1; % default
            end


            if nargin < 3
                reslice = 1;
            end


            outputfile = fullfile(obj.FuncDir, 'coregister.log'); % header info

            if check == 1 && nargin > 1
                if exist(outputfile, 'file') ~= 2
                    do_it = 1;
                else
                    do_it = 0;
                    disp('Coregistration is already done.')
                    disp(['Output is here: ' outputfile])
                end
            else
                do_it = 1;
            end

            if do_it == 1
                cd(obj.FuncDir)
                disp('**********************************************')
                disp(['Coregistration for ' obj.sub ]);

                %% step 1: coregister mean fMRI to T1
                % prepare the brain image to coregister to
                stFile = fullfile(obj.StrucDir, ['brain_' obj.T1_img ',1']);
                mean1  = fullfile(obj.FuncDir, ['mean' obj.epi{1}]);
                files = {};
                %---------------------
                for ie = 1:length(obj.epi) % if there are multiple sessions, loop through them
                    func_file = obj.epi{ie};
                    files = [files; cellstr(spm_select('expand', fullfile(obj.FuncDir, func_file)))]; %#ok<AGROW>
                end
                % to reslice for GLM modelling:
                if reslice == 1
                    EPItoT1(stFile, mean1, files, is4D)
                    % not reslice, only estimate to visually inspect images:
                else
                    %-----------------------------------------------------------------------
                    clear matlabbatch
                    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {stFile};
                    matlabbatch{1}.spm.spatial.coreg.estimate.source = {mean1};
                    %%
                    matlabbatch{1}.spm.spatial.coreg.estimate.other = files;
                    %%
                    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
                    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
                    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

                    %-------------
                    % run batch
                    disp('Begin to run...');
                    spm_jobman('run',matlabbatch);

                    % rename the coregistered files
                    for ie = 1:length(obj.epi)
                        func_file = obj.epi{ie};
                        copyfile(fullfile(obj.FuncDir, func_file), fullfile(obj.FuncDir, ['c' func_file]));
                    end
                end
                %--> register if the coregistration is done
                save(fullfile(obj.FuncDir, 'coregister.log'), 'do_it')

                disp('Job completed.');
                disp('   ')
            end

        end


        %% %%%%%%%%%%%%%%%%%%% Denoising %%%%%%%%%%%%%%%%%%%%%
        %%  Model Specification: GLM0
        function glm0(obj, check)

            if isempty(obj.GLM0_dir)
                obj.GLM0_dir = fullfile(obj.FuncDir, 'model', 'GLM0');
            end
            % Dummy Model Specification: GLM0
            outputfile = fullfile(obj.GLM0_dir,'SPM.mat');

            if check == 1 && nargin > 1
                if exist(outputfile, 'file') ~= 2
                    do_it = 1;
                else
                    do_it = 0;
                    disp('GLM0 specification is already done.')
                    disp(['Output is here: ' outputfile])
                end
            else
                do_it = 1;
            end

            if do_it == 1
                mkdir(obj.GLM0_dir)
                %% run matlabbatch
                %--------------------- set up matlabbatch -------------------
                clear matlabbatch
                % Specify the output directory where the SPM.mat file goes.
                matlabbatch{1}.spm.stats.fmri_spec.dir = {obj.GLM0_dir };
                % Scanner info.
                matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
                if isempty(obj.TR)
                    % extract the tr from the nifiti header file
                    niihdr = spm_vol(fullfile(obj.FuncDir, obj.epi{1})); % read a nifti to get the TR: each session has to have the same TR!!
                    obj.TR = niihdr(1).private.timing.tspace;
                end
                %------------------------
                matlabbatch{1}.spm.stats.fmri_spec.timing.RT = obj.TR;
                matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
                matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
                %===== loop through sessions =====
                % ------------------------------------
                if ischar(obj.epi)
                    obj.epi = cellstr(obj.epi);
                elseif isempty(obj.epi) % assuming default file names
                    fs = dir(fullfile(obj.FuncDir, 'rest*.nii'));
                    obj.epi = {fs.name}';
                end

                for ie = 1:length(obj.epi) % if there are multiple sessions, loop through them
                    %crfunc_file = spm_select('list', obj.FuncDir, '^cr.*\.nii'); % auto select aligned and coregistered images
                    func_file = obj.epi{ie};
                    % ==== check the coregistration has been done ===
                    if exist(fullfile(obj.FuncDir, ['c' func_file]), 'file') ~= 2
                        error(['Check if the ' func_file 'has been coregistered to the T1.'])
                    end
                    %----------------
                    files = cellstr(spm_select('expand', fullfile(obj.FuncDir, ['c' func_file]))); % select the coregistered images
                    % ---------
                    % load the scans
                    matlabbatch{1}.spm.stats.fmri_spec.sess(ie).scans = files;
                    matlabbatch{1}.spm.stats.fmri_spec.sess(ie).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
                    matlabbatch{1}.spm.stats.fmri_spec.sess(ie).multi = {''};
                    matlabbatch{1}.spm.stats.fmri_spec.sess(ie).regress.name = 'linear trend'; % detrend
                    matlabbatch{1}.spm.stats.fmri_spec.sess(ie).regress.val = 1:size(files,1);
                    % select headmovement file
                    if isempty(obj.movement_parameters) || length(obj.movement_parameters)~=length(obj.epi)
                        obj.movement_parameters{ie} = ['rp_' strrep(func_file, '.nii', '.txt')];
                    end
                    matlabbatch{1}.spm.stats.fmri_spec.sess(ie).multi_reg = {fullfile(obj.FuncDir,  obj.movement_parameters{ie})};% head movement
                    matlabbatch{1}.spm.stats.fmri_spec.sess(ie).hpf = 128; % high pass filter
                end
                % Specify experiment parameters
                matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
                matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
                matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
                matlabbatch{1}.spm.stats.fmri_spec.global = 'None';

                % select structural file
                %brain_mask = spm_select('list', obj.StrucDir, '^brain_mask.*\.nii');
                brain_mask = ['brain_mask_' obj.T1_img];


                if ~exist(fullfile(obj.StrucDir, brain_mask), 'file')
                    disp('No explicily defined mask is found, use SPM implicit mask...')
                    matlabbatch{1}.spm.stats.fmri_spec.mthresh = -Inf;
                    matlabbatch{1}.spm.stats.fmri_spec.mask = {fullfile(obj.spm_path, 'tpm/TPM_SPM_brain.nii')};
                else
                    matlabbatch{1}.spm.stats.fmri_spec.mthresh = -Inf;
                    matlabbatch{1}.spm.stats.fmri_spec.mask = {fullfile(obj.StrucDir, brain_mask)};
                end

                matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
                %
                disp('**********************************************')
                disp('Begin to run GLM0....')
                % obj.debug_batch = matlabbatch;
                spm_jobman('run',matlabbatch);
                disp('GLM0 specification completed.');
                disp('   ')

            end
        end % GLM0 specification

        %% GLM0 estimation
        function glm0est(obj,check)
            %------
            if isempty(obj.GLM0_dir)
                obj.GLM0_dir = fullfile(obj.FuncDir, 'model', 'GLM0');
                mkdir(obj.GLM0_dir)
            end
            outputfile = fullfile(obj.GLM0_dir, 'mask.nii');
            %--------------------------------
            if check == 1 && nargin > 1
                if exist(outputfile, 'file') ~= 2
                    do_it = 1;
                else
                    do_it = 0;
                    disp('GLM0 estimation is already done.')
                    disp(['Output is here: ' obj.GLM0_dir])
                end
            else
                do_it = 1;
            end
            %--------------------------------
            if do_it == 1
                clear matlabbatch
                % Model estimation
                disp(['Estimate GLM0 for Subject: ', obj.sub]);
                matlabbatch{1}.spm.stats.fmri_est.spmmat = {[obj.GLM0_dir,filesep, 'SPM.mat']}; %Select the SPM.mat
                matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
                matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
                % run batch
                disp('Begin to run...');
                spm_jobman('run',matlabbatch);
                disp('GLM0 estimation job completed.');
                disp('   ')

            end % if do_it
        end % function

        %% GLM0 contrast
        function glm0cont(obj, check)

            if isempty(obj.GLM0_dir)
                obj.GLM0_dir = fullfile(obj.FuncDir,'model','GLM0');
            end

            if ischar(obj.epi)
                obj.epi = cellstr(obj.epi);
            elseif isempty(obj.epi) % assuming default file names
                fs = dir(fullfile(obj.FuncDir, 'rest*.nii'));
                obj.epi = {fs.name}';
            end

            outputfile = fullfile(obj.GLM0_dir, sprintf('ess_%04d.nii', length(obj.epi)));
            %--------------------------------
            if check == 1 && nargin > 1
                if exist(outputfile, 'file') ~= 2
                    do_it = 1;
                else
                    do_it = 0;
                    disp('GLM0 contrast is already done.')
                    disp(['Output is here: ' obj.GLM0_dir])
                end
            else
                do_it = 1;
            end
            %--------------------------------
            if do_it == 1
                clear matlabbatch
                % Model estimation
                disp(['GLM0 Contrast for Subject: ', obj.sub]);
                % load SPM for info
                s = load(fullfile(obj.GLM0_dir,'SPM.mat'));
                %-------------------------------
                matlabbatch{1}.spm.stats.con.spmmat = {fullfile(obj.GLM0_dir, 'SPM.mat')};
                %-------------------------
                matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = 'all';
                matlabbatch{1}.spm.stats.con.consess{1}.fcon.weights = eye(length(s.SPM.Sess(1).C.name)); %require same model for all sessions
                matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'sess';
                matlabbatch{1}.spm.stats.con.delete = 1;
                % run batch
                disp('Begin to run...');
                spm_jobman('run',matlabbatch);
                disp('Contrasts job completed.');
                disp('   ')

            end % if do_it
        end % function

        %% ArtRepair Motion Scrub
        function motionScrub(obj, check)
            % ---------------------------------------------------------------------
            if ischar(obj.epi)
                obj.epi = cellstr(obj.epi);
            elseif isempty(obj.epi) % assuming default file names
                fs = dir(fullfile(obj.FuncDir, 'rest*.nii'));
                obj.epi = {fs.name}';
            end

            for ie = 1:length(obj.epi) % if there are multiple sessions, loop through them
                func_file = obj.epi{ie};
                fprintf('Ready to process for %s, session: %d....\n', obj.sub, ie)

                outputfile = fullfile(obj.FuncDir , ['vc' func_file]);
                if isempty(obj.denoised_epi) || length(obj.denoised_epi)~= length(obj.epi)
                    obj.denoised_epi{ie,1} = ['vc' func_file];
                end
                %--------------------------------
                if check == 1 && nargin > 1
                    if exist(outputfile, 'file') ~= 2
                        do_it = 1;
                        disp('Ready to do motion scrub')
                    else
                        do_it = 0;
                        disp('Motion scrub is already done.')
                        disp(['Output is here: ' obj.GLM0_dir])
                    end
                else
                    do_it = 1;
                    disp('Ready to do motion scrub')
                end
                %---------------------------------------------
                if do_it % check if file has been processed
                    % files = cellstr(spm_select('expand', fullfile(obj.FuncDir,func_file)));
                    files = spm_select('expand', fullfile(obj.FuncDir, ['c' func_file]));
                    rp = ['rp_' strrep(func_file, '.nii','.txt')];
                    movement_txt = fullfile(obj.FuncDir, rp);

                    % apply art_global
                    art_global(files, movement_txt, 4, 1) %4 = automask, 1 for ArtifactRepair alone (0.5 movement and add margin)
                    % rename to avoid overwriting when processing different
                    % sessions separately
                    movefile(fullfile(obj.FuncDir, 'art_deweighted.txt'), fullfile(obj.FuncDir, ['art_deweighted' num2str(ie) '.txt']));
                    movefile(fullfile(obj.FuncDir, 'art_repaired.txt'), fullfile(obj.FuncDir, ['art_repaired' num2str(ie) '.txt']));
                    listing = dir(fullfile(fileparts(obj.FuncDir), 'artglobal*.jpg'));
                    fnames = {listing.name};
                    fdates = [listing.datenum];
                    movefile(fullfile(fileparts(obj.FuncDir), fnames{fdates==max(fdates)}), fullfile(fileparts(obj.FuncDir), ['artglobal_' obj.sub '_func' num2str(ie) '.jpg']))
                    movefile(fullfile(obj.FuncDir, 'ArtifactMask.nii'), fullfile(obj.FuncDir, ['ArtifactMask' num2str(ie) '.nii']))
                    disp('Default file names are renamed to correspond to sessions.')
                    disp(['Art repaired fast-motion file for ' obj.sub '.'])
                end
            end
        end

        %% Smoothing
        function smooth(obj, check, fwhm)
            % spatial smoothing: deliberate downgrade the spatial
            % resolution to increase signal-to-noise ratio
            % fwhm is the full width at half maximum of the Gaussian
            % smoothing kernel in mm. A single value is taken, the funtion
            % will make repeat it to a 1-by-3 array for the 3 dimensions of the kernel.
            if nargin < 3
                fwhm = 4; % the default is 4mm
            elseif nargin < 2
                check = 0; fwhm = 4;
            end
            obj.smooth_fwhm = fwhm;
            %--------
            if isempty(obj.GLM0_dir)
                obj.GLM0_dir = fullfile(obj.FuncDir, 'model', 'GLM0');
            end
            % ------------------------------------
            if ischar(obj.epi)
                obj.epi = cellstr(obj.epi);
            elseif isempty(obj.epi) % assuming default file names
                fs = dir(fullfile(obj.FuncDir, 'rest*.nii'));
                obj.epi = {fs.name}';
            end

            for ie = 1:length(obj.epi) % if there are multiple sessions, loop through them
                func_file = obj.epi{ie};
                fprintf('Ready to process for %s, session: %d....\n', obj.sub, ie)

                % ---------
                if isempty(obj.denoised_epi) % default is to choose artrepaired files
                    obj.denoised_epi{ie} =  ['vc' func_file]; % select the realigned, coregistered and repaired images
                end

                outputfile = fullfile(obj.FuncDir, ['s' obj.denoised_epi{ie}]);

                if check == 1 && nargin > 1
                    if exist(outputfile, 'file') ~= 2
                        do_it = 1;
                    else
                        do_it = 0;
                        disp('smoothing is already done.')
                        disp(['Output is here: ' outputfile])
                    end
                else
                    do_it = 1;
                end
                %-----------------------------------------------
                if do_it == 1
                    disp(['Spatial smoothing for ' obj.sub ]);
                    clear matlabbatch
                    %--
                    if exist(fullfile(obj.FuncDir, obj.denoised_epi{ie}),'file') ~=2
                        error('No denoised files are found.')
                    end
                    %--
                    denoised_epis = cellstr(spm_select('expand', fullfile(obj.FuncDir, obj.denoised_epi{ie})));

                    matlabbatch{1}.spm.spatial.smooth.data =  denoised_epis;
                    matlabbatch{1}.spm.spatial.smooth.fwhm = repmat(fwhm, 1,3);
                    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
                    matlabbatch{1}.spm.spatial.smooth.im = 0;
                    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
                    % run batch
                    disp('Begin to run...');
                    spm_jobman('run',matlabbatch)
                    disp('Job completed.');
                    disp('   ')
                end
            end
        end
        %% Extract timeseries of VOIs
        function noiseSignExtr(obj, check, i_option, n_comp)
            % Variable Explained:
            % (1) 'i_option': an index of denoising options (to extract signals from
            % aCompCor or GSR or both)
            % a signle number from 0, 1, 2 is taken
            % 1 = aCompCor, 2 = GSR, 0 = both
            % (2) 'n_comp': the number of components to be saved (used for follow-up analysis)
            if nargin < 3
                i_option = 0; n_comp = 6;
            elseif nargin < 4
                n_comp = 6;
            end
            %----------
            if isempty(obj.GLM0_dir)
                obj.GLM0_dir = fullfile(obj.FuncDir, 'model', 'GLM0');
            end
            %-------------------------------------
            denoise_options = {'aCompCor', 'GSR'};
            if i_option == 0
                i_loopin = 1:length(denoise_options);
            else
                i_loopin = i_option;
            end
            %---
            if ischar(obj.epi)
                obj.epi = cellstr(obj.epi);
            elseif isempty(obj.epi) % assuming default file names
                fs = dir(fullfile(obj.FuncDir, 'rest*.nii'));
                obj.epi = {fs.name}';
            end
            %==========================================
            for ie = 1:length(obj.epi)
                disp(['Noise extraction, ready for ' obj.sub ', Session: ' obj.epi{ie}])
                %==========================================
                for i = i_loopin
                    option = denoise_options{i};
                    % ---------------------------------------------------------------------
                    outputfile = fullfile(obj.GLM0_dir, ['VOI_' option '_' num2str(ie) '.mat']);
                    PC_file =  fullfile(obj.GLM0_dir,[option  '_PC_components_' num2str(ie) '.mat']);
                    %--------------------------------
                    if check == 1 && nargin > 1
                        if exist(outputfile, 'file') ~= 2
                            do_it = 1;
                        else
                            do_it = 0;
                            disp('Signal is already extracted')
                        end
                    else
                        do_it = 1;
                    end
                    %---------------------------------------------
                    if do_it

                        clear matlabbatch
                        disp(['Applying ' option 'for subject ', obj.sub '.']);

                        matlabbatch{1}.spm.util.voi.spmmat = {[obj.GLM0_dir filesep 'SPM.mat']};
                        matlabbatch{1}.spm.util.voi.adjust = ie; % adjust for each session
                        matlabbatch{1}.spm.util.voi.session = ie;
                        matlabbatch{1}.spm.util.voi.name = option;

                        matlabbatch{1}.spm.util.voi.roi{1}.mask.image = {fullfile(obj.GLM0_dir,  'mask.nii,1')};
                        matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.2;

                        if isequal(option, 'aCompCor')
                            wm_mask = ['c2' obj.T1_img];
                            matlabbatch{1}.spm.util.voi.roi{2}.mask.image = {fullfile(obj.StrucDir, wm_mask)};
                            matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.99;
                            csf_mask = ['c3' obj.T1_img];
                            matlabbatch{1}.spm.util.voi.roi{3}.mask.image = {fullfile(obj.StrucDir, csf_mask)};
                            matlabbatch{1}.spm.util.voi.roi{3}.mask.threshold = 0.99;
                            matlabbatch{1}.spm.util.voi.expression = '(i2+i3) & i1'; % csf,wm conbined AND within the whole brain mask
                        elseif isequal(option, 'GSR')
                            gm_mask = ['c1' obj.T1_img];
                            matlabbatch{1}.spm.util.voi.roi{2}.mask.image = {fullfile(obj.StrucDir, gm_mask)};
                            matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.99;
                            matlabbatch{1}.spm.util.voi.expression = 'i2 & i1';
                        end
                        %
                        % save batch file for review
                        spm_jobman('run',matlabbatch);
                        disp([option ' extraction job completed.']);
                        disp('   ')
                    end % if do_it

                    %% check if the PCs are written

                    if exist(PC_file, 'file') ~=2  || check == 0
                        disp ('Calculating PCA now.....')
                        % --------------------------------------------------------
                        % get the first few components (default = first 6)
                        f = load(outputfile);
                        [coeff,score,latent,tsquared,explained,mu] = pca(f.xY.y); % princomp(f.xY.y,'econ')
                        PC_components = zscore(score(:, 1:n_comp),0, 1); % 0 meaning using the sample sd; take the first 6 components as suggested in the paper

                        %% plot for extracted timeseries - quality check
                        % ------> plot1
                        % plot for variance explained for the first 6
                        % components
                        perc_explained = explained/sum(explained);
                        close
                        fig1 = figure;
                        plot(perc_explained, 'k')
                        ylabel('Percentage of variance explained')
                        xlabel('Number of PC')
                        %hold on
                        line([n_comp n_comp], [0 max(perc_explained)], 'Color','red');
                        x_lim = 100;
                        xlim([0 x_lim])
                        text(x_lim/2 ,0.7*max(perc_explained), ...
                            [num2str(round(sum(explained(1:n_comp))/sum(explained),4)*100) '% variance explained'])
                        % ===== save ========
                        save(fullfile(obj.GLM0_dir, ['pca_' option '_' num2str(ie) '.mat']), 'coeff','score','latent','tsquared','explained','mu')
                        saveas(fig1, fullfile(obj.GLM0_dir, [option '_variance_explained_' num2str(ie) '.png']))
                        save(fullfile(obj.GLM0_dir,[option  '_PC_components_' num2str(ie) '.mat']), 'PC_components')
                        disp('PCA saved.')
                    end
                end % for sessions
            end % loop in options
        end % function


        %% %%%%%%%%%%%%%%%%%%%%%%% Signal Extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Signals will be extracted from the canonical networks that have been registered in the native space
        % There are several network altases to use. To make better use of these
        % multiple sources of information,
        % I will generate an intersection map among them, extract the
        % signals from the indivial atlases as well as the intersection map, and
        % then examine their correlations. This setp severs as a quality check: if
        % the maps are valid, the first (few) eigenvariate(s) should be highly
        % correlated.

        %% Image Calc
        % define common network (can be used for other calculation than intersect)

        function masking(obj, check, input_masks, calculations, out_names, outdir)
            %===> variables explained:
            % (1) check: logical variable, it checks if the file has been written, if
            % so, skip this step. It can be useful for debugging.
            % (2) input_masks: masks on which mathmatic operations will be
            % performed on, it can be a cell or a structure (if structure, there are multiple ROIs with different images for operation)
            % (3) calculations
            % it is char or cell (if cell, each element is a character array)
            % i1 is the first input image, i2 is the second etc.
            % i1&i2 = take an intersect betewen them
            % i1+/-/*//i2 = voxel-wise arithmatic calculation
            % etc...To find more expressions, read SPM12 manual
            %-----------------------------------------------------------------------
            if nargin < 3
                if isempty(obj.overlay_masks)
                    obj.overlay_masks = struct();
                end
                input_masks = obj.ROIs; % just to borrow the structure, will rewrite
                calculations = cellstr(repmat('intersection', length(obj.ROIs),1)); % default
                out_names = cellstr([repmat('wCOMM_', length(obj.ROIs),1 ), char(reshape({obj.ROIs.name},[],1))]); % naming convention: w prefix indicates the mask in the native space
                outdir = obj.StrucDir;
            elseif nargin < 4 % define default values: default is to take the common area/intersection
                calculations = {repmat('intersection', length(obj.ROIs),1)}; % default
                out_names = cellstr([repmat('wCOMM_', length(obj.ROIs),1 ), char(reshape({obj.ROIs.name},[],1))]); % naming convention: w prefix indicates the mask in the native space
                outdir = obj.StrucDir;
            elseif nargin < 5
                out_names = cellstr([repmat('wCOMM_', length(obj.ROIs),1 ), char(reshape({obj.ROIs.name},[],1))]);
                outdir = obj.StrucDir;
            elseif nargin < 6
                outdir = obj.StrucDir;
            end
            % --------
            if ischar(calculations)
                calculations = cellstr(calculations);
            end
            %-------
            if ischar(out_names)
                out_names = cellstr(out_names);
            end
            %-----------
            if isempty(obj.ROIs)
                error('No ROI/individual mask is given.')
            end

            %%
            for ir = 1:length(obj.ROIs)
                %===========================================================
                outputfile = fullfile(outdir, [out_names{ir} '.nii']);
                if check == 1 && nargin > 1
                    if exist(outputfile, 'file') ~= 2
                        do_it = 1;
                    else
                        do_it = 0;
                        disp('Image calc is already done.')
                        disp(['File is here: ' outputfile])
                    end
                else
                    do_it = 1;
                end
                %---------------------------------------------------
                if do_it == 1
                    %---------------------------------------
                    if isstruct(input_masks)
                        if isequal(input_masks, obj.ROIs)
                            % rewrite for future record
                            mask_paths = obj.ROIs(ir).path;
                            candidate_masks = cell(length(mask_paths),1);
                            for ip = 1:length(mask_paths)
                                [~, roi_name, ext] = fileparts(mask_paths{ip});
                                candidate_masks{ip,1} = fullfile(obj.StrucDir, ['w' roi_name ext]);
                                disp('Assuming that the ROIs have been warpped to the navtive space and have a prefix of ''w''...')
                            end
                            obj.overlay_masks(ir).name = obj.ROIs(ir).name;
                            obj.overlay_masks(ir).path = candidate_masks;
                        else
                            candidate_masks = obj.overlay_masks(ir).path;
                        end % is equal

                        %---------------------------------------
                    elseif iscell(input_masks)
                        candidate_masks = input_masks;
                        disp('Use external masks for ImageCalc....')
                    else
                        error('The variable overlay_masks should either be a cell or a structure.')
                    end
                    %============================================
                    n = length(candidate_masks);
                    calculation = calculations{ir};
                    if isequal(calculation, 'intersection') % pipeline default
                        calculation = strjoin(cellstr([repmat('i',n,1) num2str((1:n)')]), '&'); % need modification if n >= 10, as the string matrix won't have equal length anymore
                    end
                    %=====================================================
                    clear matlabatch
                    disp(['Processing for ROI: ' obj.overlay_masks(ir).name '.....'])
                    matlabbatch{1}.spm.util.imcalc.input =  candidate_masks;
                    matlabbatch{1}.spm.util.imcalc.output = out_names{ir};
                    matlabbatch{1}.spm.util.imcalc.outdir = {outdir};

                    matlabbatch{1}.spm.util.imcalc.expression = calculation;
                    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
                    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
                    matlabbatch{1}.spm.util.imcalc.options.mask = -1;
                    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
                    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
                    % feedback
                    disp(' ')
                    disp('Ready to do Image Calc ....')
                    disp(' ')
                    % run
                    obj.debug_batch = matlabbatch;
                    spm_jobman('run', matlabbatch);
                    % --
                    disp('Image Calc done!')
                end % if do_it
                %                 % ===== assuming overlay masks and the common mask are all ROIs
                %                 if ~ismember(outputfile, obj.overlay_masks(ir).path)
                %                     obj.overlay_masks(ir).path{end+1} = outputfile;
                %                 end
            end
        end % IMageCalc


        %% ==============Signal extract: GLM1 specification, estimation, contrast and extraction=====================

        %%  Model Specification: GLM1
        function glm1(obj, check, smooth_option, index, bad_volumes)
            %% GLM1 Model Specification
            % -- can be used for seed-based functional connectivity
            % --- Property obj.M is a structure variale that specify the GLM1_dir and main
            %     variales in the GLM. It has fields of "GLM1_dir" and "main_vars":
            %     the main_vars is a table.mat, which contains main variables for statistical testing
            %     of the main effects of interest. If the main_vars is empty,
            %     then the GLM1 will only have the default non-neuronal
            %     covariates, which is a procedure equal to denoising.
            %     index points to which model to run, can be used for
            %     debugging.
            % -- smooth_option: 1 = 'smoothed'
            %    or 0 ='unsmoothed'; the default is to use 'smoothed' images
            %    to conform to the standard SPM pipeline; however,
            %    for functional connectivity studies on individual space,
            %    using unsmoothed images might be a better option.
            % --bad_volumes (take a numeric vector/ cell/the specific option: 'repaired'):
            %   -- (1) if vector (assuming same for each session) or a cell (each vector for each session),
            %  they are indices to the bad volumes, which will be
            %   masked out from the statistical testing. The artrepair auto-detected bad volumes will be ignored.
            %   -- if cell, each cell corresponds to a session (in order),
            %   and each cell contains a vector which are the indices to
            %   the bad volumes
            %   -- (2) if use 'none', it will use all the data, regardless of
            %   the motion warning, no exclusion will be performed.
            %   -- (3) if use 'repaired' (default, i.e. when nothing is given), then use the  bad volumes
            %    which were automatically detected by the Artrepair. When
            %   -- (4) if empty; The default option is to mask out the bad volumes from the
            %     following analysis.
            % --- Note that for the default denoising 'GLM1', i.e. when the obj.M.main_vars == [];
            %     no scans will be masked out. This is important for
            %     extracting timeseries with every frame in time.
            % --- A special case is reserved for the case when the user want to
            %     use the artrepaired volumes (a interpolated scan from the
            % --- adjacent volumes, an average of before and after) for the
            %     following analysis. To do that, specify the bad_volumes as 'repaired'.
            % --- By default, if the bad volumes exceed 1/3 of all the volumes,
            %     a warning will be generated. It has been suggested in the
            %     literature to exclude this subject. An exclude marker will be
            %     saved to the FuncDir to facilitate a future search.
            % --- The first 5 volumes will always be masked out as an fMRI
            %     convention (the time for the scanner to be stabilized)
            %------------------------------------------

            % define default
            if isempty(obj.M)
                obj.M = struct(); % use default, denoise without GSR
                obj.M.GLM1_dir = fullfile(obj.FuncDir, 'model', 'GLM1');
                obj.M.main_vars = [];
                obj.M.multiregressorTXT = '';
            end

            %--------
            if nargin < 3
                smooth_option = 1;
                disp('By default, using smoothed images.')
            end
            %--------
            if nargin < 4 || isempty(index)
                if length(obj.M)>1
                    index = 1:length(obj.M);
                else
                    index = 1;
                end
            end
            %--------
            if nargin < 5 || isempty(bad_volumes)
                bad_volumes = [];
            end

            %====================
            if ~ismember(smooth_option, [0 1] ) || isempty(smooth_option)
                error('The input var ''smooth_option'' should either be 1 (smoothed) or 0 (unsmoothed)....!')
            end
            smoothing = {'unsmoothed', 'smoothed'};
            %-----
            fprintf('Input %s images are specified if use default SPM naming convention.\n', smoothing{smooth_option+1})
            %-----
            if isempty(obj.exclusion)
                obj.exclusion = struct('GLMdir', {}, 'epi', {});
            end
            %------------------------------------------
            for i = index
                M_ = obj.M(i);
                %--------
                if ~iscell(M_.main_vars) && ~isempty(M_.main_vars)
                    M_.main_vars{1} = M_.main_vars;
                end
                [~, glm1folder] = fileparts(M_.GLM1_dir);
                %------------
                outputfile = fullfile(M_.GLM1_dir, 'SPM.mat');

                if check == 1 && nargin > 1
                    if exist(outputfile, 'file') ~= 2
                        do_it = 1;
                    else
                        do_it = 0;
                        disp('GLM1 specification is already done.')
                        disp(['Output is here: ' outputfile])
                    end
                else
                    do_it = 1;
                end
                %----------------------------------------------------------------------
                if do_it == 1

                    mkdir(M_.GLM1_dir)
                    [~,GLM,~] = fileparts(M_.GLM1_dir);
                    disp(['GLM1 specification for ' GLM])
                    % prepare exclusion variables
                    obj.exclusion(end+1).GLMdir = M_.GLM1_dir;
                    obj.exclusion(end).epi      = {};

                    %% run matlabbatch
                    %--------------------- set up matlabbatch -------------------
                    clear matlabbatch
                    % Specify the output directory where the SPM.mat file goes.
                    matlabbatch{1}.spm.stats.fmri_spec.dir = {M_.GLM1_dir};
                    % Scanner info.
                    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
                    if isempty(obj.TR)
                        % extract the tr from the nifiti header file
                        niihdr = spm_vol(fullfile(obj.FuncDir, obj.epi{1}));  % assume the sessions have the same scanning TR
                        obj.TR = niihdr(1).private.timing.tspace;
                    end
                    %------------------------
                    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = obj.TR;
                    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
                    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

                    % ------------------------------------------------------------------
                    if isempty(obj.denoised_epi)  %if not given, the default is to choose artrepaired files
                        % ------------------------------------
                        if ischar(obj.epi)
                            obj.epi = cellstr(obj.epi);
                        elseif isempty(obj.epi) % assuming default file names
                            fs = dir(fullfile(obj.FuncDir, 'rest*.nii'));
                            obj.epi = {fs.name}';
                        end
                        %-----------------------------------
                        default_denoised_epi = cellstr([repmat('c',length(obj.epi),1), char(reshape(obj.epi,[],1))]); % select the realigned, coregistered images
                        obj.denoised_epi = default_denoised_epi;
                        use_default = 1;
                    else
                        if ischar(obj.denoised_epi)
                            obj.denoised_epi = cellstr(obj.denoised_epi);
                        end
                        use_default = 0;
                    end


                    %% ===== loop through sessions =====
                    ie_clean = 0; % <- a counter to actually write into the GLM (in case of excluding sessions)

                    for ie = 1:length(obj.denoised_epi) % if there are multiple sessions, loop through them
                        % ---------
                        % load the scans
                        %------------------------------------------------------------------------------------------------------------
                        if smooth_option == 1 && use_default == 1
                            %======================================================================%
                            if obj.denoised_epi{ie}(1) ~= 's'
                                obj.denoised_epi{ie} = ['s' obj.denoised_epi{ie}];
                            end
                        end
                        %======================================================================%
                        files = cellstr(spm_select('expand', fullfile(obj.FuncDir, [obj.denoised_epi{ie}]))); % % select the denoised files
                        %======================================================================%

                        if use_default == 1
                            %% sort out bad volumes
                            % The default is to mask out the ArtRepaired volumes
                            % and the first 5 volumns (scanner unstable)
                            % read artrepaired volumn index
                            if isempty(bad_volumes) || isequal(bad_volumes, 'repaired')
                                artfile = fullfile(obj.FuncDir, ['art_repaired' num2str(ie) '.txt']);
                                if exist(artfile, 'file') ~= 2
                                    warning('ArtRepair has not been run...')
                                else
                                    fileID = fopen(artfile,'r');
                                    A = fscanf(fileID, '%f');
                                end
                                %----
                                mask_volumes = unique([(1:5)' ; A]); % discard the first 5 volumes (scanner instability)

                            elseif isequal(bad_volumes, 'none')
                                mask_volumes = [];
                            else
                                if isnumeric(bad_volumes)
                                    mask_volumes = bad_volumes;
                                elseif iscell(bad_volumes)
                                    mask_volumes = bad_volumes{ie};
                                end
                            end

                            %----------
                            if (length(mask_volumes) > length(files)/3) && ~strcmp(glm1folder, 'GLM1') % if it is a denoising GLM (GLM1), always specify every session, even if the session is heavily corrupted
                                disp('==================================')
                                warning (['The session ' obj.epi{ie} ' contains too many bad volumes. It will be excluded for further analysis by default.'])
                                obj.exclusion(end).epi        = [obj.exclusion(end).epi; obj.epi{ie}];
                                fprintf('Excluding the No.%.d session: %s\n', ie, obj.epi{ie})
                                disp(['Located in the folder: ' obj.FuncDir])
                                disp('==================================')
                                continue; % will go back to the loop immediately
                            else
                                ie_clean = ie_clean + 1;
                                [~, f_check, ~] = fileparts(files{1});
                                fprintf('\nLoading the No.%.d session: %s to the GLM.\n', ie_clean, [f_check, '.nii'])
                            end
                            V_good = 1:length(files);
                            % NOTE: masking does not apply to a denoising GLM1,
                            % or using explicitly specified denoised files
                            if ~isempty(M_.main_vars) || size(M_.main_vars, 1)==length(files)
                                V_good(mask_volumes) = []; % knock out the bad volumes
                            end

                        elseif use_default==0 % if not use default
                            ie_clean = ie_clean + 1;
                            % >>>>>>>>---indices to the good volumes--->>>>>>>>
                            V_good = 1:length(files);
                        end

                        %-----------------------------------------------------------------------
                        matlabbatch{1}.spm.stats.fmri_spec.sess(ie_clean).scans = files(V_good);
                        matlabbatch{1}.spm.stats.fmri_spec.sess(ie_clean).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
                        matlabbatch{1}.spm.stats.fmri_spec.sess(ie_clean).multi = {''};
                        %====>> COVARIATES LOADED HERE:
                        %--% total length = length(main_vars) (Main effects) + length(noise)
                        %------------------ Main effects (seeds) --------------
                        if ~isempty(M_.main_vars)
                            % M_.main_vars should be a cell containing
                            % structure(s); each cell correspondes to a
                            % session, while the structure inside corresponds to the main
                            % effect(s), with two fields: name and vector
                            VAR = M_.main_vars{ie}; % a structure
                            if isstruct(VAR)
                                main_vec = VAR.vec;
                                main_name = VAR.name;
                            elseif isnumeric(VAR)
                                main_vec = VAR;
                                main_name = [strrep(obj.denoised_epi{ie}, '.nii', '_') 'PC'];
                            end

                            for i_m = 1:size(VAR,2) % can contain multiple PCs
                                matlabbatch{1}.spm.stats.fmri_spec.sess(ie_clean).regress(i_m).name = [main_name '_' num2str(i_m)];%[M_.main_vars.Properties.VariableNames{i_m} '_' num2str(ie)]
                                matlabbatch{1}.spm.stats.fmri_spec.sess(ie_clean).regress(i_m).val = main_vec(V_good, i_m);%table2array(M_.main_vars(V_good,i_m));
                            end
                        else
                            i_m = 0;
                        end
                        %--------------------Covariates to denoise ----------------------
                        % 1. detrend
                        matlabbatch{1}.spm.stats.fmri_spec.sess(ie_clean).regress(1+i_m).name = ['linear trend_sess' num2str(ie)]; % detrend
                        matlabbatch{1}.spm.stats.fmri_spec.sess(ie_clean).regress(1+i_m).val = V_good; % it is a vector
                        %--% 2:length(noise) (assuming GLM0 is already done)
                        i_n = 0;
                        if ~strcmp(obj.noise, 'none') % unless obj.noise is explicitly specified to be none
                            if isempty(obj.noise) || length(obj.noise)~=length(obj.epi) % noise vectors for each session, if not given, use default, i.e. no GSR
                                % load aCompCor
                                clear f
                                f = load(fullfile(obj.GLM0_dir, ['aCompCor_PC_components_' num2str(ie) '.mat']));
                                obj.noise{ie} = f.PC_components;
                            end
                            %-----
                            for i_n = 1:size(obj.noise{ie}, 2)
                                matlabbatch{1}.spm.stats.fmri_spec.sess(ie_clean).regress(1 + i_m + i_n).name = ['physio_' num2str(i_n) '_' num2str(ie)];
                                matlabbatch{1}.spm.stats.fmri_spec.sess(ie_clean).regress(1 + i_m + i_n).val = obj.noise{ie}(V_good,i_n);
                            end

                        end
                        %---------- confounds in a text file-----------------
                        if isempty(M_.multiregressorTXT) % when the multiple regressor text file is not explicitly provided, use headmotion file generated from the realignment step
                            % select headmovement file
                            if isempty(obj.movement_parameters) || length(obj.movement_parameters) ~= length(obj.epi)
                                obj.movement_parameters{ie} = ['rp_' strrep(obj.epi{ie}, '.nii', '.txt')];
                            end
                            rpfile = fullfile(obj.FuncDir, obj.movement_parameters{ie});
                            M_.multiregressorTXT{ie} = rpfile;
                        elseif ischar(M_.multiregressorTXT)
                            M_.multiregressorTXT = cellstr(M_.multiregressorTXT);
                        elseif strcmp(obj.movement_parameters, 'none')
                            M_.multiregressorTXT = cellstr('');

                        end

                        multiregressorTXT = fullfile(obj.FuncDir, M_.multiregressorTXT{ie});
                        % load rp data
                        multidat = readmatrix(multiregressorTXT);
                        for i_p = 1:size(multidat,2)
                            matlabbatch{1}.spm.stats.fmri_spec.sess(ie_clean).regress(1 + i_m + i_n + i_p).name = ['multi_' num2str(i_p) '_' num2str(ie)];
                            matlabbatch{1}.spm.stats.fmri_spec.sess(ie_clean).regress(1 + i_m + i_n + i_p).val = multidat(V_good, i_p);
                        end
                        %----------------------------------------------------------------------------
                        matlabbatch{1}.spm.stats.fmri_spec.sess(ie_clean).hpf = 128; % high pass filter (default)

                    end
                    % Specify experiment parameters
                    matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
                    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
                    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
                    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';

                    % select structural file
                    %brain_mask = spm_select('list', obj.StrucDir, '^brain_mask.*\.nii');
                    % if use implicit mask
                    brain_mask = ['brain_mask_' obj.T1_img];
                    if ~exist(fullfile(obj.StrucDir, brain_mask), 'file')
                        disp('No explicily defined mask is found, use SPM implicit mask...')
                        matlabbatch{1}.spm.stats.fmri_spec.mthresh = -Inf;
                        matlabbatch{1}.spm.stats.fmri_spec.mask = {fullfile(obj.spm_path, 'tpm/TPM_SPM_brain.nii')};
                    else
                        matlabbatch{1}.spm.stats.fmri_spec.mthresh = -Inf;
                        matlabbatch{1}.spm.stats.fmri_spec.mask = {fullfile(obj.StrucDir, brain_mask)};
                    end
                    matlabbatch{1}.spm.stats.fmri_sspec.cvi = 'AR(1)';
                    %
                    disp('   ')
                    disp('>>>>>>>>>>>>>>>>>>>>>>>')
                    disp('Begin to run...')
                    spm_jobman('run',matlabbatch);
                    disp('GLM1 specification completed.');
                    disp('   ')
                end % for different GLM1 models
            end
        end % GLM1 specification

        %% GLM1 estimation
        function glm1est(obj, check, index, save_residuals)
            % 'index' refers to which GLM1 model to use; it is useful when
            % there are multiple first-level analyses to run
            % if save_residual == 1, residual files from the GLM1 will be
            % written
            %------------------------------------------
            % define default
            if isempty(obj.M)
                obj.M = struct(); % use default, denoise without GSR
                obj.M.GLM1_dir = fullfile(obj.FuncDir, 'model', 'GLM1');
                obj.M.main_vars = [];
            end

            if nargin < 4
                save_residuals = zeros(length(obj.M),1);  % default is not to save residual
            elseif length(save_residuals) == 1
                save_residuals = ones(length(obj.M),1).*save_residuals;
            end

            if nargin < 3
                index = [];
            end

            if isempty(index)
                index = 1:length(obj.M) ;
            end
            %------------------------------------------
            for i = index
                M_ = obj.M(i);
                outputfile = fullfile(M_.GLM1_dir, 'mask.nii');
                save_residual =  save_residuals(i);
                %--------------------------------
                if check == 1 && nargin > 1
                    if exist(outputfile, 'file') ~= 2
                        do_it = 1;
                    else
                        do_it = 0;
                        disp('GLM1 estimation is already done.')
                        disp(['Output is here: ' M_.GLM1_dir])
                    end
                else
                    do_it = 1;
                end
                %--------------------------------
                if do_it == 1
                    clear matlabbatch
                    % Model estimation
                    disp(['Estimate GLM1 for Subject: ', obj.sub]);
                    matlabbatch{1}.spm.stats.fmri_est.spmmat = {[M_.GLM1_dir,'/SPM.mat']}; %Select the SPM.mat
                    matlabbatch{1}.spm.stats.fmri_est.write_residuals = save_residual; % take value 0 or 1
                    matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
                    % run batch
                    disp('Begin to run...');
                    spm_jobman('run', matlabbatch);
                    disp('GLM1 estimation job completed.');
                    disp('   ')

                end % if do_it
            end % for index of the GLM1
        end % function

        %% GLM1 contrast
        function glm1cont(obj,check, index, CON)
            % CON is a structure with fields "name", "weights" and "sessrep"
            % CON.name specifies a operation name,
            % CON.weights is a multiplier to the beta, which
            % assigns weights to column/factors in the indexed GLM:
            % e.g., [1, -1] means the first two colums are to be contrasted
            % Note: this function has a nested looping of 'CON' through the given 'index',
            % i.e. it will repeat the same set of CON to each of the
            % given GLM1 (sepcified in the index)
            % CON.sessrep takes a character array, indicating the contrast
            % replication options over session. If there are multiple
            % sessions with identical conditions, one might wnat to specify
            % contrasts which are identical over sessions. This can be done
            % automatically based on the contrast sepc for one session. See
            % SPM12 manual for more information. In this pipeline, the
            % following three options are used: 'replsc'
            % (1)Replicate & Scale (for the same main effect estimated over
            % sessions): 'replsc' --> a single set of stats files will be written
            % (2)Per session (for the session-specific main effect): 'sess'
            % --> multiple sets of stats files will be written, each for
            % one session
            % (3)none (not replicate over session): 'none'--> a single set
            % of stats files will be written; index strictly corresponds to
            % the GLM collumn
            %------------------------------------------
            % define default
            if isempty(obj.M)
                obj.M = struct(); % use default, denoise without GSR
                obj.M.GLM1_dir = fullfile(obj.FuncDir, 'model', 'GLM1');
                obj.M.main_vars = []; CON = [];
            end

            if nargin < 3
                index = 1:length(obj.M) ;
                CON = [];
            elseif nargin < 4
                CON = [];
            end
            %------------------------------------------
            if ischar(obj.epi)
                obj.epi = cellstr(obj.epi);
            elseif isempty(obj.epi) % assuming default file names
                fs = dir(fullfile(obj.FuncDir, 'rest*.nii'));
                obj.epi = {fs.name}';
            end

            %------------------------------------------
            for i = index
                M_ = obj.M(i);

                s = load(fullfile(M_.GLM1_dir, 'SPM.mat'));

                %----------
                if isempty(CON) && index == 1
                    CON.name = 'all';
                    CON.sessrep = 'sess';
                    CON.weights = eye(length(s.SPM.Sess(1).col));
                elseif isempty(CON) && index ~= 1
                    error('Contrasts are not defined...')
                end

                %-------
                if ~isfield(s.SPM.xCon, 'name')
                    do_it = 1;
                else % check if the contrasts have all been written
                    con_written = reshape({s.SPM.xCon.name}, [],1);
                    %-----------
                    Nsess = length(s.SPM.Sess);
                    con_towrite = {};
                    for icon = 1:length(CON)
                        name = CON(icon).name;
                        switch CON(icon).sessrep
                            case 'sess'
                                conname = cellstr([repmat([name ' - Session '], Nsess ,1) num2str((1:Nsess)')]);
                            case 'replsc'
                                conname = cellstr([name ' - All Sessions']);
                            case 'none'
                                conname = cellstr(name);
                        end
                        con_towrite = [con_towrite; conname];
                    end
                    %--------------------------------
                    if check == 1 && nargin > 1
                        if ~isequal(con_written, con_towrite)
                            do_it = 1;
                        else
                            do_it = 0;
                            disp('GLM1 contrasts are already done.')
                            disp(['Output is here: ' M_.GLM1_dir])
                        end
                    else % if do not check
                        do_it = 1;
                    end
                end

                %--------------GLM1 CONTRAST----------------
                if do_it == 1

                    clear matlabbatch
                    [~, GLM, ~] = fileparts(M_.GLM1_dir);
                    fprintf('\n GLM1 Contrast for Subject: %s, model: %s \n', obj.sub, GLM);

                    %-------------------------------
                    matlabbatch{1}.spm.stats.con.spmmat = {fullfile(M_.GLM1_dir, 'SPM.mat')};
                    %====== To run contrasts ======
                    for ic = 1:length(CON)
                        con_name = CON(ic).name;
                        con_weights = CON(ic).weights;
                        con_sessrep = CON(ic).sessrep;

                        if size(con_weights,1)==1 || size(con_weights,2)==1  % it is a T test
                            matlabbatch{1}.spm.stats.con.consess{ic}.tcon.name = con_name;
                            matlabbatch{1}.spm.stats.con.consess{ic}.tcon.weights = con_weights;
                            matlabbatch{1}.spm.stats.con.consess{ic}.tcon.sessrep =  con_sessrep; % replicate over sessions

                        else % it is an F test

                            matlabbatch{1}.spm.stats.con.consess{ic}.fcon.name = con_name;
                            matlabbatch{1}.spm.stats.con.consess{ic}.fcon.weights = con_weights;
                            matlabbatch{1}.spm.stats.con.consess{ic}.fcon.sessrep =  con_sessrep; % replicate over sessions
                        end
                    end
                    %----
                    matlabbatch{1}.spm.stats.con.delete = 1; % delete all previously run contrasts
                    % run batch
                    disp('Begin to run...');
                    obj.debug_batch = matlabbatch;
                    spm_jobman('run',matlabbatch);
                    disp('Contrasts job completed.');
                    disp('   ')

                end % if do_it
            end % for index of the GLM1
        end % function

        %% result report
        function glm1result(obj, check, index, primary_thresh, extent, report_nii, i_CON, atlas_indices)
            % NOTE: this function can be run twice, the first time to
            % get the significance testing, so that we know what cluster size to
            % threshold (coded as 'extent'). The second time it will save the significant
            % clusters.
            % variable explained: index is the index of the models
            % specified in obj.M for GLM1
            % extent: the mininum number of cluster size (no. of voxels)
            % to show (usually it takes the value of the smallest
            % significant cluster), default is 2 voxels
            % save nifti will save the nifti that exceed the threshold, a
            % logical value is taken, 0/1
            %------------------------------------------
            % define default
            if nargin < 2
                check = 0; % default is to redo the analysis without checking if done
            end
            %------------------define default-----------------------
            if nargin < 3 && isempty(obj.M)
                obj.M = struct(); % use default, denoise without GSR
                obj.M.GLM1_dir = fullfile(obj.FuncDir, 'model', 'GLM1');
                obj.M.main_vars = [];
                index = 1; extent = 2; report_nii = 0;
                i_CON = [];
            elseif nargin < 3 && ~isempty(obj.M)
                index = 1:length(obj.M) ;
                extent = 2;  report_nii = 0; i_CON = [];
            elseif nargin < 4
                primary_thresh = 0.001; extent = 2;  report_nii = 0; i_CON = [];
            elseif nargin < 5
                extent = 2; report_nii = 0; i_CON = [];
            elseif nargin < 6
                report_nii = 0; i_CON = [];
            elseif nargin < 7
                i_CON = []; atlas_indices = [1 2];
            elseif nargin < 8 % atlas index refers to the input to the external function mni2atlas.m
                atlas_indices = [1 2];
               %{'Juelich Histological Atlas'              }
               %{'Harvard-Oxford Cortical Structural Atlas'}
            end
            %------------------------------------------
            for i = index
                M_ = obj.M(i);
                s = load(fullfile(M_.GLM1_dir, 'SPM.mat'));
                %-----------check contrast------------------
                if isempty(s.SPM.xCon)
                    [~, GLM,~] = fileparts(M_.GLM1_dir);
                    error([ GLM ' has no contrasts defined!'])
                else
                    con_written = reshape({s.SPM.xCon.name}, [],1);
                    if isempty(i_CON)
                        i_CON = 1:length(con_written);
                    end
                end
                %-------------------------------------------------
                for ic = i_CON
                    cfile = dir(fullfile(M_.GLM1_dir, ['spm*_' sprintf('%04d', ic) '.nii']));
                    cfile = cfile(1).name(1:end-4);
                    contrast_name = s.SPM.xCon(ic).name;
                    sigTfile = fullfile(M_.GLM1_dir, ['SignificanceTable_' cfile '_' contrast_name '.csv']);
                    if report_nii == 1
                        outputfile = fullfile(M_.GLM1_dir, [cfile '_thr.nii']);
                    end

                    %--------------------------------
                    if check == 1 && nargin > 1
                        if exist(sigTfile, 'file') ~= 2
                            report_sigT = 1;
                        else
                            report_sigT = 0;
                            if exist(outputfile, 'file') ~= 2 && report_nii == 1
                                report_nii = 1;
                            else
                                report_nii = 0;
                            end
                            disp(['Significance table is already written for ' cfile])
                            disp(['Result is written as: ' outputfile])
                        end
                    else
                        report_sigT = 1;
                        if report_nii == 1
                            report_nii = 1;
                        else
                            report_nii = 0;
                        end
                    end
                    %========= Result report ===================
                    clear matlabbatch
                    % ==================
                    matlabbatch{1}.spm.stats.results.spmmat = {fullfile(M_.GLM1_dir, 'SPM.mat')};
                    matlabbatch{1}.spm.stats.results.conspec(1).titlestr = '';
                    matlabbatch{1}.spm.stats.results.conspec(1).contrasts = ic;
                    matlabbatch{1}.spm.stats.results.conspec(1).threshdesc = 'none'; % no family wise correction at the voxel level comparison
                    matlabbatch{1}.spm.stats.results.conspec(1).thresh = primary_thresh;
                    matlabbatch{1}.spm.stats.results.conspec(1).extent = extent; % do not show cluster that is smaller than 2 voxels
                    matlabbatch{1}.spm.stats.results.conspec(1).conjunction = 1;
                    matlabbatch{1}.spm.stats.results.conspec(1).mask.none = 1;
                    matlabbatch{1}.spm.stats.results.units = 1; % volumetric
                    matlabbatch{1}.spm.stats.results.export{1}.csv = true;
                    %--------------------------------
                    if report_sigT == 1

                        % run batch
                        disp('Begin to calculate signifcance ... ');
                        obj.debug_batch = matlabbatch;
                        spm_jobman('run',matlabbatch);
                        %-------------------------------------------
                        % a workaround to save csv with customized name:
                        cF = dir(fullfile(M_.GLM1_dir, 'spm_*.csv'));
                        dates = cellfun(@datetime, {cF.date}');
                        [~,I] = sort(dates,'descend');
                        new_csv = cF(I(1)).name;
                        csv = importdata(fullfile(M_.GLM1_dir,new_csv));
                        %                         opts = detectImportOptions((fullfile(M_.GLM1_dir,new_csv)), 'NumHeaderLines', 1);
                        %                         opts.VariableNamesLine = 2;
                        %                         sigT = readtable(fullfile(M_.GLM1_dir,new_csv), opts);
                        if isfield(csv, 'data')
                            VariableNames = {'set p', 'set c', 'Cluster p(FWE-corr)', 'Cluster p(FDR-corr)', 'Cluster size (equivk)', ...
                                'Cluster p(unc)', 'Peak (FWE-corr)', 'Peak (FDR-corr)', 'Peak F', 'Peak equivZ',...
                                'Peak p(unc)', 'x (mm)', 'y (mm)', 'z (mm)'};
                            sigT = array2table(csv.data,...
                                'VariableNames', VariableNames);

                           % try
                                atlas_labels = cell(size(sigT,1),numel(atlas_indices));
                               
                                for ir = 1:size(sigT,1)
                                atlas = mni2atlas(table2array(sigT(ir,end-2:end)), atlas_indices);
                                
                                for il = 1:numel(atlas_indices)
                                    if isempty(atlas(il).label)
                                        atlas_labels{ir,il} = 'NaN';
                                    else
                                         atlas_labels{ir,il} = strjoin(atlas(il).label);
                                    end
                                end

                                end
                                sigT_label = array2table(atlas_labels,...
                                'VariableNames', {atlas.name});
                                sigT = [sigT sigT_label];
                                
%                             catch
%                                 fprintf('\n Unable to print out automatic labeling. \n Check if MNI2Atlas toolbox (and its dependencies) is added to path.\n \n')
%                             end

                            writetable(sigT, sigTfile)
                            disp('Significance table is saved.');
                        else
                            disp('There is no significant result to be saved.')
                            report_nii = 0;
                        end
                        delete(fullfile(M_.GLM1_dir, new_csv))
                    end

                    %-----------------------
                    if report_nii == 1
                        % if save significant clusters, the function will be
                        % run twice, the first to generate statistical results,
                        % the  second time is to threshold the result and save
                        % the image.
                        % ---->>>>>>>>>> filter significance table >>>>>>>>>>>
                        % default is to use FWE-corrected significant clusters
                        if exist(sigTfile, 'file') == 2
                            sigT = readtable(sigTfile);
                            K = table2array(sigT(sigT.ClusterP_FWE_corr_ < 0.05, 5)); % select the column for the cluster size
                            if isempty(K)
                                disp('There is no significant cluster to be saved.')
                            else
                                Kmin = min(K);

                                %------ modify the previous matlabbatch -----------
                                matlabbatch{1}.spm.stats.results.conspec(1).extent = Kmin;
                                matlabbatch{1}.spm.stats.results.export{1} = struct();% clear the previous one
                                matlabbatch{1}.spm.stats.results.export{1}.tspm.basename = 'thr';
                                matlabbatch{1}.spm.stats.results.export{2}.binary.basename = 'bin';
                                % run batch
                                disp('Saving thresholded results...');
                                spm_jobman('run',matlabbatch);
                                %obj.debug_batch = matlabbatch;
                            end
                        else
                            disp('No significance table was found.')
                        end % is exit significance table
                    end % if save nii
                end % for ic
            end % for i = index (GLM1)
        end % function

        %% Extract timeseries of VOIs
        function roiSignExtr(obj, check, masks, i_masks, n_comp, index)
            % Variable Explained:
            % a denoising GLM1 is assumed to have existed. Then all the covariates are noise vectors
            % (1) index is the index to the denoised GLM1 to use to extract
            % signals
            % (2) masks where the signals are extracted; if masks is a
            % structure with a mask name and a vectorised
            % coordinates, then the sepcified sphere centered around the coordinate is used.
            % (3) 'i_masks': an index/indices of overlay_masks from which to extract signal
            %  0 meaning looping though all ROIs
            % (4) 'n_comp': the number of components to be saved (used for follow-up analysis)

            %-----------
            if isempty(obj.M)
                obj.M = struct(); % use default, i.e. denoise without GSR
                obj.M.GLM1_dir = fullfile(obj.FuncDir, 'model', 'GLM1');
                obj.M.main_vars = [];
                index = 1;
            end
            %-----------------------------
            if nargin < 3
                if isempty(obj.overlay_masks) && isempty(obj.ROIs)
                    error('No overlay_mask or ROI is provided.')
                end

                f = dir(fullfile(obj.StrucDir, 'wCOMM_*.nii'));
                masks = cellstr([repmat([char(obj.StrucDir) filesep], length(f), 1),  char(reshape({f.name}', [], 1))]);
            end

            if nargin < 4 || isempty(i_masks)
                i_masks = 1:length(masks);
            end

            if nargin < 5
                n_comp = 6;
            end

            if nargin < 6
                index = 1:length(obj.M);
            end
            %-------------------------------------

            %% ----  loop though models ------
            for i = index
                M_ = obj.M(i);
                % load SPM to read pars
                p = load(fullfile(M_.GLM1_dir, 'SPM.mat'));
                nsess = length(p.SPM.Sess);
                %============loop though all rois================
                for ir = i_masks

                    if iscell(masks)
                        [~, roi, ~] = fileparts(masks{ir}) ;
                    elseif isstruct(masks)
                        roi = masks(ir).name;
                        coord = masks(ir).coord;
                        radius = masks(ir).radius;
                    end

                    %-----------------------
                    for i_sess = 1:nsess
                        outputfile = fullfile(M_.GLM1_dir, ['VOI_' roi '_' num2str(i_sess) '.mat']);
                        %--------------------------------
                        if check == 1 && nargin > 1
                            if exist(outputfile, 'file') ~= 2
                                do_it = 1;
                            else
                                do_it = 0;
                                fprintf('Session %d: Signal extraction is already done for %s\n.' ,i_sess, roi)
                                disp(['Output is here: ' M_.GLM1_dir])
                            end
                        else
                            do_it = 1;
                        end
                        % --------- part 1: VOI extraction from SPM ------------------
                        if do_it

                            clear matlabbatch
                            disp(['Applying ' roi ' for subject ', obj.sub '.']);

                            matlabbatch{1}.spm.util.voi.spmmat = {fullfile(M_.GLM1_dir, 'SPM.mat')};
                            matlabbatch{1}.spm.util.voi.adjust = NaN; % adjust for everything
                            matlabbatch{1}.spm.util.voi.session = i_sess;
                            matlabbatch{1}.spm.util.voi.name = roi;

                            matlabbatch{1}.spm.util.voi.roi{1}.mask.image = {fullfile(M_.GLM1_dir,  'mask.nii,1')};
                            matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.2;

                            if isstruct(masks)
                                matlabbatch{1}.spm.util.voi.roi{2}.sphere.centre = coord;
                                matlabbatch{1}.spm.util.voi.roi{2}.sphere.radius = radius;
                                matlabbatch{1}.spm.util.voi.roi{2}.sphere.move.fixed = 1;

                            elseif iscell(masks)
                                matlabbatch{1}.spm.util.voi.roi{2}.mask.image = masks(ir);
                                matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.5;
                            end

                            matlabbatch{1}.spm.util.voi.expression = 'i2 & i1';

                            % run batch
                            obj.debug_batch = matlabbatch;
                            spm_jobman('run',matlabbatch);
                            fprintf('%s signal extraction job completed for session %d.\n', roi, i_sess);
                            disp('   ')
                        end  % if do_it

                        if n_comp == 1
                            continue;
                        end
                        % ---------part 2: to manually extract the the principle components---------
                        if (exist(fullfile(M_.GLM1_dir,[roi '_PC_components_' num2str(i_sess) '.mat']), 'file')~=2 && check == 1) || check == 0
                            % --------------------------------------------------------
                            % get the first few components (default = first 6)
                            V = load(fullfile(M_.GLM1_dir, ['VOI_' roi '_' num2str(i_sess) '.mat']));
                            if size(V.xY.y,2) > 6  % sometimes the warpped VOI is too small

                                [coeff,score,latent,tsquared,explained,mu] = pca(V.xY.y);
                                PC_components = zscore(score(:, 1:n_comp),0, 1); % 0 meaning using the sample sd; take the first 6 components as suggested in the paper

                                %% plot for extracted timeseries - quality check

                                % ------> plot1
                                % plot for variance explained for the first 6
                                % components
                                close
                                fig1 = figure;
                                perc_explained = explained/sum(explained);
                                plot(perc_explained, 'k')

                                %hold on
                                line([n_comp n_comp], [0 max(perc_explained)], 'Color','red');
                                x_lim = 100;
                                xlim([0 x_lim])
                                text(x_lim/2 ,0.7*max(perc_explained), ...
                                    [num2str(round(sum(explained(1:n_comp))/sum(explained),4)*100) '% variance explained'])
                                xlabel('Number of PC'); ylabel('Percentage of variance explained')
                                %===== save ======
                                save(fullfile(M_.GLM1_dir, ['pca_' roi '_' num2str(i_sess) '.mat']), 'coeff','score','latent','tsquared','explained','mu')
                                saveas(fig1, fullfile(M_.GLM1_dir, [roi '_variance_explained_' num2str(i_sess) '.png']))
                                save(fullfile(M_.GLM1_dir,[roi '_PC_components_' num2str(i_sess) '.mat']), 'PC_components')
                            else
                                warning('The VOI is too small: less than 6 voxels exist. Suggest to use only the first PC (i.e. SPM default output).')
                            end
                        end % if do PCA
                    end % for ie
                end % loop in rois
            end % loop in GLM1s
        end % function

    end

end

