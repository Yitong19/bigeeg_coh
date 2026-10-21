function CM_DICS_MEEG(param,CM_inv,CM_CSD)

%  CM_DICS_MEEG(param,inv,CM_CSD)
%  or CM_DICS(param,inv,CM_CSD,inv_ref)
%   Compute coherence at the source level
%
%   input : 
%    
%     param :        structure specifying the parameters to compute the overlay
%      .subjects_dir : freesurfer subjects' directory 
%                      (default : '/usr/local/freesurfer/subjects/')
%      .subject :      name given to the subject in "subjects_dir"
%      .exp :          folder in freesurfer where overlay will be written
%                      (default : new_fold)
%      .overly_mane :  name of the overlay 
%                      (default : name of the MEG fif file present in CM)
%      .dxyz :         size of the voxels (default: [1 1 1])
%      .w :            wether to normlize the overlay to MNI or back to
%                      individual subject's mri (default : 1)
%
%     CM_inv :       inverse solution computed with CM_sf_...
%     CM_CSD :       cross-spectral density as obtained using 
%                    CM_coh_MEG_ref_one_pass
%
%  Examples :
% 

%   % -> Load the inverse solution
% fwdfile = fullfile(exp_dir,'subjects',subjects{n_sub}(1:3),'meg',[subjects{n_sub}(1:3) '_from_MNI-AVSE-5-src-fwd.fif']);
% matfile = fullfile(exp_dir,'matfiles',[subjects{n_sub}(1:3) '_noise_cov.mat']);
% load(matfile,'CM_noise')
% 
% % fwd = mne_read_forward_solution(fwdfile);
% load([fwdfile(1:end-4) '-2comp.mat'])
% 
% CM_inv = MNE_inverse_MEEG(fwd,CM_noise.COV);
% CM_inv.fwdfile = fwdfile;
%


if ~isfield(param,'subjects_dir')
    param.subjects_dir = '/usr/local/freesurfer/subjects/';
end
if ~isfield(param,'fold_fs')
    param.fold_fs = 'new_fold';
end
if ~isfield(param,'exp')
    [temp1,param.exp,temp3,temp4] = fileparts(CM_CSD.infile);
end
if ~isfield(param,'w')
    param.w = 1;
end


% Load or create the correctly reshaped MRI
if ~exist(fullfile(param.subjects_dir , param.subject , 'norm'))
    CM_configure_norm(param.subject,param.subjects_dir);
end
if isempty(regexp(CM_inv.fwdfile,'from_MNI'))
    S = spm_vol(fullfile(param.subjects_dir , param.subject , 'norm' , 'orig.nii'));
else
    S = spm_vol(fullfile(param.subjects_dir , param.subject , 'norm' , 'MNI.nii'));
end


% create directory to store the overlays if missing
if ~exist(fullfile(param.subjects_dir , param.subject , param.fold_fs))
    mkdir(fullfile(param.subjects_dir , param.subject , param.fold_fs));
end


% Check whether the beamformer depends on the frequency or not.
if isfield(CM_inv,'Nf') % beamformer differs for different frequencies
    Lbf = CM_inv.Nf;
else                 % beamformer independennt of the frequency
    Lbf = 1;
    CM_inv.Nf = 1;
end
Nvox = double(CM_inv.nsource);
Nmeeg = size(CM_inv.Tinv,2);
Nref = size(CM_CSD.CSD,1)-Nmeeg;


% parameters used in the maximization process
Nori = size(CM_inv.Tinv,1)/CM_inv.nsource;
if Nori == 3
    ndir = 100;
    phi = -pi/2:pi/ndir:pi/2-pi/ndir;
    Ntheta = ceil(cos(phi)*ndir/2+0.000001);
    inds = [0 cumsum(Ntheta)];
    vec = zeros(sum(Ntheta),3);
    for n_phi = 1:length(phi)
        theta = (1:Ntheta(n_phi))/Ntheta(n_phi)*pi;
        vec(inds(n_phi)+1:inds(n_phi+1),1) = cos(theta)*cos(phi(n_phi));
        vec(inds(n_phi)+1:inds(n_phi+1),2) = sin(theta)*cos(phi(n_phi));
        vec(inds(n_phi)+1:inds(n_phi+1),3) = sin(phi(n_phi));
    end
    % plot3(vec(:,1),vec(:,2),vec(:,3),'.')
else
    Ngrid = 1000;
    teta = [1:Ngrid]'/Ngrid*pi;
    vec = [cos(teta) sin(teta)];
end


coh = zeros(length(CM_CSD.F),Nvox);
for f = 1:length(CM_CSD.F)
    if Lbf == 1      % ff = index of the third row of the inverse solution.
        ff = 1;      %      If the beamformer does not depend on the 
    else             %      frequency, ff = 1 all the time.
        ff = f;
    end
    disp(['starting coherence calculation at frequency F=' num2str(CM_CSD.F(f)) 'Hz with ' CM_inv.Tinfo])
    
    
    % Compute the coherence for all sources using DICS approach
    % computed with optimal direction
    Cyy = abs(CM_CSD.CSD(Nmeeg+1,Nmeeg+1,f));
    T_CSD = CM_inv.Tinv(:,:,ff)*CM_CSD.CSD(1:Nmeeg,1:Nmeeg,f);
    Cxy = CM_inv.Tinv(:,:,ff)*CM_CSD.CSD(1:Nmeeg,Nmeeg+1,f);
    if Nref == 2
        Czz = abs(CM_CSD.CSD(Nmeeg+2,Nmeeg+2,f));
        Cxz = CM_inv.Tinv(:,:,ff)*CM_CSD.CSD(1:Nmeeg,Nmeeg+2,f);
        Rzy = CM_CSD.CSD(Nmeeg+2,Nmeeg+1,f)/sqrt(Cyy*Czz);
    end
    for n = 1:Nvox
        int = (-Nori+1:0)+Nori*n; % indices of the source under consideration
        
        Cxx = T_CSD(int,:)*CM_inv.Tinv(int,:,ff)';
        % regularize Cxx
        [U,Sd,V] = svd(Cxx);
        Sd(1:Nori+1:end) = max(Sd(1:Nori+1:end),Sd(1:Nori+1:end)/50);
        Cxx = U*Sd*V';
        
        
        Ax = sqrt(abs(sum(vec*Cxx.*vec,2))); % sqrt of power in all directions
        Rxy = (vec*Cxy(int))./Ax/sqrt(Cyy); % Coherencecy
        if Nref == 1
            coh(f,n) = max(abs(Rxy).^2);
            % coh(n) = max(abs(eig(Cxx\Cxy*(Cyy\Cxy'))));  % cca
        elseif Nref == 2
            Rxz = (vec*Cxz(int))./Ax/sqrt(Czz); % Coherency
            cohtest = abs(Rxy-Rxz*Rzy).^2./((1-abs(Rxz).^2)*(1-abs(Rzy)^2)); % Partial coherence
            coh(f,n) = max(cohtest);
        end
    end
end
coh = mean(coh,1);


% interpolate on a 1-mm grid
param.fwdfile = CM_inv.fwdfile;
[map_mri,S] = CM_overlay(coh,param);


% Create the overlay normalized to MNI template
if isempty(regexp(CM_inv.fwdfile,'from_MNI'))
    S.fname = fullfile(param.subjects_dir , param.subject , param.fold_fs , [param.exp '_' CM_inv.Tinfo '_' num2str(CM_CSD.F(f)) '_Hz_Cmax_' num2str(max(map_mri(:))) '.nii']);
    S = spm_write_vol(S,map_mri);
    if param.w == 1
        load(fullfile(param.subjects_dir , param.subject , 'norm' , 'forward_norm.mat'))
        spm_write_sn(S.fname,prm,dnrm.write);
    end
else
    S.fname = fullfile(param.subjects_dir , param.subject , param.fold_fs , ['w' param.exp '_' CM_inv.Tinfo '_' num2str(CM_CSD.F(f)) '_Hz_Cmax_' num2str(max(map_mri(:))) '.nii']);
    S = spm_write_vol(S,map_mri);
    if param.w == 1
        load(fullfile(param.subjects_dir , param.subject , 'norm' , 'inverse_norm.mat'))
        dnrm.write.prefix = 'i';
        spm_write_sn(S.fname,iprm,dnrm.write);
    end
end


% inform that the overlay is now written in a file
disp(['overlay calculation finished at this frequency for ' CM_inv.Tinfo])
