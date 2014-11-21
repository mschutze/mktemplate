%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPM8 custom PET template creation script
% By Manuel Schutze - november 2013 - mschutze@gmail.com
% 
% First change b.Tmp on line 44!
% Then run script from inside folder containing aligned images:
% 1. Nifti MRI images named MRI_$.nii, where $ is the patient ID
% 2. Nifti PET images named PET_$.nii, where $ is the patient ID
% Example: MRI_PT001.nii and PET_PT001.nii
% The result is a template file called s08customPET.nii
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = mktemplate

clear all;
curdir = pwd;
 
%% Find MRI .nii images
files = dir('MRI_*.nii');
numFiles = length(files);

if(numFiles == 0)
    fprintf('No nifti files found in current folder!\n');
else
    message = ['Found ',int2str(numFiles),' subjects in current folder. Starting processing...\n'];
    fprintf(message);
end

b.wrPET = cell(1,numFiles);
 
%% For each MRI/PET pair run corregistration, segmentation and normalization
for i=1:numFiles
    %define variables for individual subjects
    curFile = files(i).name; % "MRI_PT001.nii"
    temp = regexp(curFile,'_','split');
    indiv = regexp(temp{2},'\.','split');
    ID = indiv{1}; % "PT001"
    b.MRI = strcat(curdir,'/MRI_',ID,'.nii');
    b.PET = strcat(curdir,'/PET_',ID,'.nii');
    b.MRIseg = strcat(curdir,'/MRI_',ID,'_seg_sn.mat');
    b.rPET = strcat(curdir,'/rPET_',ID,'.nii');
    cur_wrPET = strcat(curdir,'/wrPET_',ID,'.nii');
    b.wrPET{1,i} = cur_wrPET; %will be used later to create soft mean
    b.Tpm = { %change this!
           '/Users/manuel/Documents/MATLAB/spm8/tpm/grey.nii'
           '/Users/manuel/Documents/MATLAB/spm8/tpm/white.nii'
           '/Users/manuel/Documents/MATLAB/spm8/tpm/csf.nii'
           };
    
    if(0==0) %change one 0 to 1 to skip this part
        %specify matlabbatch variable with subject-specific inputs
        matlabbatch = batch_job(b);
        %run matlabbatch job
        message = ['Running functions for subject ',ID,' (',int2str(i),' of ',int2str(numFiles),')...\n'];
        fprintf(message);
        spm_jobman('initcfg');
        spm('defaults', 'PET');
        spm_jobman('run', matlabbatch);
    end
 
end

%% Create soft mean between normalized PET images and then smooth result
%Generate soft mean expression to be used by SPM8 ImCalc
global exp;
exp = '(';
for k=1:numFiles
    if k==1
        exp = strcat(exp,'i',int2str(k));
    else
        exp = strcat(exp,'+i',int2str(k));
    end
end
exp = strcat(exp,')./(eps');
for l=1:numFiles
    exp = strcat(exp,'+(i',int2str(l),'~=0)');
end
exp = strcat(exp,')');

b.Expression = exp; %something like 'i1+iN./(exp+(i1~=0)+(iN~=0))'
b.customPET = strcat(curdir,'/customPET.nii'); %path to soft mean img

%run ImCalc and then 8mm smooth in matlabbatch job
matlabbatch2 = batch_job2(b);
message = ['Calculating soft mean of ',int2str(numFiles),' subjects...\n'];
fprintf(message);
spm_jobman('initcfg');
spm('defaults', 'PET');
spm_jobman('run', matlabbatch2);

fprintf('Finished! Please rename and copy s08customPET to the templates folder.\n');

end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [matlabbatch]=batch_job(b)
%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4252 $)
%-----------------------------------------------------------------------

% PET/MRI Corregistration (Coregister: Estimate & Reslice)
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {b.MRI};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {b.PET};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 1;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

% MRI Segmentation (Segment)
matlabbatch{2}.spm.spatial.preproc.data = {b.MRI};
matlabbatch{2}.spm.spatial.preproc.output.GM = [0 0 0];
matlabbatch{2}.spm.spatial.preproc.output.WM = [0 0 0];
matlabbatch{2}.spm.spatial.preproc.output.CSF = [0 0 0];
matlabbatch{2}.spm.spatial.preproc.output.biascor = 1;
matlabbatch{2}.spm.spatial.preproc.output.cleanup = 0;
matlabbatch{2}.spm.spatial.preproc.opts.tpm = b.Tpm;
matlabbatch{2}.spm.spatial.preproc.opts.ngaus = [2
                                                 2
                                                 2
                                                 4];
matlabbatch{2}.spm.spatial.preproc.opts.regtype = 'mni';
matlabbatch{2}.spm.spatial.preproc.opts.warpreg = 1;
matlabbatch{2}.spm.spatial.preproc.opts.warpco = 25;
matlabbatch{2}.spm.spatial.preproc.opts.biasreg = 0.0001;
matlabbatch{2}.spm.spatial.preproc.opts.biasfwhm = 60;
matlabbatch{2}.spm.spatial.preproc.opts.samp = 3;
matlabbatch{2}.spm.spatial.preproc.opts.msk = {''};

% PET Normalisation (Normalise: Write)
matlabbatch{3}.spm.spatial.normalise.write.subj.matname = {b.MRIseg};
matlabbatch{3}.spm.spatial.normalise.write.subj.resample = {b.rPET};
matlabbatch{3}.spm.spatial.normalise.write.roptions.preserve = 0;
matlabbatch{3}.spm.spatial.normalise.write.roptions.bb = [-90 -126 -72
                                                          90 90 108];
matlabbatch{3}.spm.spatial.normalise.write.roptions.vox = [2 2 2];
matlabbatch{3}.spm.spatial.normalise.write.roptions.interp = 1;
matlabbatch{3}.spm.spatial.normalise.write.roptions.wrap = [0 0 0];
matlabbatch{3}.spm.spatial.normalise.write.roptions.prefix = 'w';
end

function [matlabbatch2] = batch_job2(b)
% PET soft mean (SPM > Image Calculator)
matlabbatch2{1}.spm.util.imcalc.input = b.wrPET;
matlabbatch2{1}.spm.util.imcalc.output = 'customPET.nii';
matlabbatch2{1}.spm.util.imcalc.outdir = {''};
matlabbatch2{1}.spm.util.imcalc.expression = b.Expression;
matlabbatch2{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch2{1}.spm.util.imcalc.options.mask = 0;
matlabbatch2{1}.spm.util.imcalc.options.interp = 1;
matlabbatch2{1}.spm.util.imcalc.options.dtype = 16; %saving as 8bit UINT
% Smooth result (SPM > Smooth)
matlabbatch2{2}.spm.spatial.smooth.data = {b.customPET};
matlabbatch2{2}.spm.spatial.smooth.fwhm = [8 8 8];
matlabbatch2{2}.spm.spatial.smooth.dtype = 0;
matlabbatch2{2}.spm.spatial.smooth.im = 0;
matlabbatch2{2}.spm.spatial.smooth.prefix = 's08';
end