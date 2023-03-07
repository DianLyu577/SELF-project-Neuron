% Import electrode LEPTOVOX coordinates
home_dir = expanduser('~');
subj = 'S21_166';
subDir = dir(sprintf('/Volumes/neurology_jparvizi$/SHICEP_%s*/**/elec_recon', subj));

Fcoord = fullfile(subDir(1).folder,[subj '.LEPTOVOX']);
chan_leptovox = dlmread(Fcoord, ' ',2,0);

Fname  = fullfile(subDir(1).folder,[subj '.electrodeNames']);
txt = readlines(Fname);
chan_nam   = cellstr(split(txt(3:length(chan_leptovox)+2), ' '));

gunzip(fullfile(subDir(1).folder, 'T1.nii.gz'), fullfile(home_dir, 'Downloads') );
t1 = fullfile(home_dir, 'Downloads', 'T1.nii');
chan_nativecoord = leptoVox2ScannerRAS(chan_leptovox, t1);

%% write table
T = table();
T.SUBJECT = cellstr(repmat(subj, size(chan_nam, 1), 1));
T.FS_label = chan_nam(:,1);
T.NativeCoord = chan_nativecoord;