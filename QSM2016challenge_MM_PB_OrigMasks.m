addpat%
%	qsm2016_script_recon_evaluation.m
%
%	Evaluation script for the QSM reconstruction challenge 2016
%
%	Find the results and report at http://qsm.neuroimaging.at
%
%	2016 BB, CL, FS
%
% Adapted to show results from different masking approaches


%%-------------------------------------------------------------------------
%% load ground truth Challenge data
%%-------------------------------------------------------------------------
chi_33=spm_read_vols(spm_vol('./FromChallenge/chi_33.nii'));                % chi_33 
chi_cosmos=spm_read_vols(spm_vol('./FromChallenge/chi_cosmos.nii'));   
N=size(chi_33);


		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%% LOAD AND ADD YOUR QSM TO THE LISTS BELOW %%%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %============== BMMR University Tuebingen ===============%
        %========= High Field MRI Max-Planck-Institute ==========%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                  Supplementary DATA                    %  
        %   how to generate masks including brain surface voxels %
        %                  described in                          %
        %'Phase-based masking for QSM of the human brain at 9.4T'%
        %Hagberg GE et al submitted September 2021, MagnResonMed %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TE=25;
field=3;
B0=[0 0 1]

%Needed toolboxes to run this code are:
%STISuite_V3.0
%MEDI
%SPM12
%fsl 6.0.0; 6.0.4
%input images: phs_wrap.nii, phs_unwrap.nii, msk.nii and magn_raw.nii chi_33 from
%the challenge data set

addpath(genpath('PathToSPM12GoesHere'));
addpath(genpath('PathToSTI_SuiteGoesHere'));
addpath(genpath('PathToMEDI_toolboxGoesHere'));

%load the wrapped phase
v=spm_vol(fullfile('./FromChallenge','phs_wrap.nii'));%load the wrapped phase
wr=spm_read_vols(v);
P=spm_imatrix(v.mat);
res=abs(P(7:9));
mtx=v.dim;

%for writing Results
a=v;
a=rmfield(a,'private');
a=rmfield(a,'pinfo');
[pth,nam,ext]=spm_fileparts(a.fname);
opth='./Results';


%generate magnitude based mask using BET
AMfile=spm_vol('./FromChallenge/magn_raw.nii');
[pth,nam,ext]=spm_fileparts(AMfile.fname);
unix(['bet ' AMfile.fname ' ' fullfile(opth,[nam '_brain' ext]) ' -R -f 0.1 -g 0 -m'])
unix(['fslchfiletype NIFTI ' fullfile(opth,[nam '_brain'])]);
unix(['fslchfiletype NIFTI ' fullfile(opth,[nam '_brain_mask'])]);
unix(['gunzip '  fullfile(opth,[nam '_brain_mask' ext])]);
fmask=(spm_vol(fullfile(opth,[nam '_brain_mask.nii'])));
sfmask=fmask;
sfmask.fname=fullfile('./Results/',[nam '_MMask' ext]);
spm_smooth(fmask,sfmask,4*res);
mask{1}=spm_read_vols(sfmask);
mask{1}=logical(mask{1});
strm{1}='MM';
disp('Done Magnitude-based Masking\n');

%generate phase based mask by further processing, see manuscrips for
%different steps and motivation
[pth,nam,ext]=spm_fileparts(v.fname);
se=strel('sphere',6);
L=del2(sign(wr));
test=convn(abs(L),se.Neighborhood,'same');
PB=mask{1}.*(test<500);
PB=imclose(PB,se);
mask{2}=round(imopen(PB,se));
strm{2}='PB';
mask{2}=logical(mask{2});
a.fname=fullfile(opth,[  nam '_PBMask' ext]);
spm_write_vol(a,mask{2});
disp('Done Phase-Based Masking\n');

%load original mask from the Challenge
mask{3}=logical(spm_read_vols(spm_vol('./FromChallenge/msk.nii')));%read the provided challenge mask
strm{3}='orig';
disp('Done loading original challenge mask\n');


%show the masks
M=spm_read_vols(AMfile);
figure(1)
imagesc(rot90(squeeze(M(end/2,:,:))));
axis off
hold on
col=hsv(3);
col(3,:)=[1 1 0];
for mval=1:3
    contour(rot90(squeeze(mask{mval}(end/2,:,:))),[1 1 ],'color',col(mval,:));
end
colormap gray
brighten(0.4)
f=legend('MM','PB','Orig');
f.Color=[0 0.5 1];



%% perform dipole inversion using three methods
%%closed form L2 solution provided in the challenge
%%QSM_iLSQR from STI_studio
%%MEDI from the MEDI toolbox

%%needed  for closed form solution, from Challenge
%%-------------------------------------------------------------------------

[ky,kx,kz] = meshgrid(-N(1)/2:N(1)/2-1, -N(2)/2:N(2)/2-1, -N(3)/2:N(3)/2-1);

kx = (kx / max(abs(kx(:)))) / res(1);
ky = (ky / max(abs(ky(:)))) / res(2);
kz = (kz / max(abs(kz(:)))) / res(3);
k2 = kx.^2 + ky.^2 + kz.^2;
R_tot = eye(3);     % orientation matrix for transverse acquisition
kernel = fftshift( 1/3 - (kx * R_tot(3,1) + ky * R_tot(3,2) + kz * R_tot(3,3)).^2 ./ (k2 + eps) );    

%% closed-form L2 recon
%%-------------------------------------------------------------------------

[k1, k2, k3] = ndgrid(0:N(1)-1,0:N(2)-1,0:N(3)-1);

E1 = 1 - exp(2i .* pi .* k1 / N(1));
E2 = 1 - exp(2i .* pi .* k2 / N(2));
E3 = 1 - exp(2i .* pi .* k3 / N(3));

EtE = abs(E1).^2 + abs(E2).^2 + abs(E3).^2;
DtD = abs(kernel).^2;
reg_param = 9e-2;   % gradient regularization parameter

%needed variables for MEDI:
iMag=spm_read_vols(spm_vol('./FromChallenge/magn_raw.nii'));
N_std=1;
delta_TE=0.025;
B0_dir=[0 0 1];
matrix_size=[160 160 160];
voxel_size=res;
CF=42.576*field*1E6;

%Load the tissue_phase, that does not need any unwrapping or background correction
%get QSM using the original mask only
pt=spm_vol('./FromChallenge/phs_tissue.nii');
uw2=spm_read_vols(pt);

[pthPT,namPT,extPT]=spm_fileparts(pt.fname);
X{1}  = real( ifftn(conj(kernel) .* fftn(uw2) ./ (DtD + reg_param * EtE))) .* mask{mval} ;
a.fname=fullfile(opth,[ 'qsmL2_' strm{3} '_' namPT extPT]);
spm_write_vol(a,X{1});

uw2=42.576*field*2*pi*TE/1000*uw2;%42.576 in MHz/T
X{2} = QSM_iLSQR(uw2,mask{3},'H',B0,'voxelsize',res,'padsize',round(mtx.*[0 0 res(3)*2]),'niter',50,'TE',TE,'B0',field); %STI-Suite function (Duke University - Dr. Chunlei Liu)
a.fname=fullfile(opth,[ 'qsmLSQR_' strm{3} '_' namPT extPT]);
spm_write_vol(a,X{2});
iFreq=uw2;
iFreq_raw=uw2;
RDF=uw2;
Mask=mask{3};
save RDF.mat RDF iFreq iFreq_raw iMag N_std Mask matrix_size...
    voxel_size delta_TE B0_dir CF;
X{3}=MEDI_L1('lambda',100,'merit');
a.fname=fullfile(opth,[ 'qsmMEDI_' strm{3} '_' namPT extPT]);
spm_write_vol(a,X{3});


%unwrap the phase using the Laplacian approach
[uw, Laplacian] = MRPhaseUnwrap(wr, 'voxelsize', res);%STI_suite
iFreq=uw;
iFreq_raw=wr;

for mval=[1: 3] %three masks
    [b0,mask_ero]=V_SHARP(uw-median(uw(mask{mval}==1)),mask{mval},'voxelsize',res, 'smvsize',20);
    a.fname=fullfile(opth,[ 'uwVSH_' strm{mval}  '_' nam ext]);
    disp('Done Backgroud field removal\n');
    spm_write_vol(a,b0);
    
    phs_tissue=b0/(42.576*field*2*pi*TE/1000);
    X{mval+3}  = real( ifftn(conj(kernel) .* fftn(phs_tissue) ./ (DtD + reg_param * EtE))) .* mask{mval} ;
    a.fname=fullfile(opth,[ 'qsmVshL2_' strm{mval} '_' nam ext]);
    spm_write_vol(a,X{mval+3});

    X{mval+6} = QSM_iLSQR(b0,mask_ero,'H',B0,'voxelsize',res,'padsize',round(mtx.*[0 0 res(3)*2]),'niter',50,'TE',TE,'B0',field); %STI-Suite function (Duke University - Dr. Chunlei Liu)
    a.fname=fullfile(opth,[ 'qsmVshLSQR_' strm{mval} '_' nam ext]);
    spm_write_vol(a,X{mval+6});
    
    RDF=b0;
    Mask=mask_ero;
    
    save RDF.mat RDF iFreq iFreq_raw iMag N_std Mask matrix_size...
        voxel_size delta_TE B0_dir CF;
    X{mval+9}=MEDI_L1('lambda',100,'merit');
    a.fname=fullfile(opth,[ 'qsmVshMEDI_' strm{mval} '_' nam ext]);
    spm_write_vol(a,X{mval+9});

end



%%

% create 4d volume of all results, and get a mask where all methods can be compared with the chi_33 and the chi_cosmos images

X{13}=chi_33;
X{14}=chi_cosmos;

mask_ero=ones(N);
for k=1:length(X)    
    mask_ero= mask_ero.*(abs(X{k}*1000)>eps);
end

for k=1:length(X)
    qsm_results(:,:,:,k) = X{k}.*mask_ero;
end

qsm_names = {'L2';'iLSQR';'MEDI';'L2vshp MM';'L2vshpPB';'L2vshpOrig';'iLSQRvshpMM';'iLSQRvshpPB';'iLSQRvshpOrig';'MEDIvshpMM';'MEDIvshpPB';'MEDIvshpOrig';'chi33'; 'cosmos'};
qsm_names = {'chi33';'L2';'iLSQR';'MEDI';'L2rshp MM';'L2rshpPB';'L2rshporig';'L2vshp MM';'L2vshpPB';'L2vshpOrig';'iLSQRrshpMM';'iLSQRrshpPB';'iLSQRrshpOrig';'iLSQRvshpMM';'iLSQRvshpPB';'iLSQRvshpOrig';'MEDIvshpMM';'MEDIvshpPB';'MEDIvshpOrig';'MEDIvshpMM';'MEDIvshpPB';'MEDIvshpOrig'};
for k=23:41
    qsm_names{k}=['ro' qsm_names{k-21}];
end

%%-------------------------------------------------------------------------
%% compute and print performance metrics
%%-------------------------------------------------------------------------
chi_33=chi_33.*mask_ero;
evaluation_mask=spm_read_vols(spm_vol('./FromChallenge/evaluation_mask.nii'));
%%
wm_mask = (evaluation_mask>6) &  mask_ero>0;
wm_px = sum(wm_mask(:));
gm_mask = (evaluation_mask>0) & (evaluation_mask<=6) &mask_ero>0;
gm_px = sum(gm_mask(:));

disp('   Algorithm: RSME  : HFEN  : SSIM   : WM ERROR : GM ERROR : W+GM ERROR : (all relative to chi33)')
for ccc = 1:(size(qsm_results, 4))
    gm = qsm_results(:,:,:,ccc);
    wm = qsm_results(:,:,:,ccc);    
    gm_error(ccc) = sum(abs((gm(gm_mask) - chi_33(gm_mask)))) / gm_px ;
    wm_error(ccc) = sum(abs((wm(wm_mask) - chi_33(wm_mask)))) / wm_px ;   
    rmse(ccc) = compute_rmse(qsm_results(:,:,:,ccc), chi_33);
    hfen(ccc) = compute_hfen(qsm_results(:,:,:,ccc), chi_33);
    ssim(ccc) = compute_ssim(qsm_results(:,:,:,ccc), chi_33);
    mean_error(ccc) = (gm_error(ccc) + wm_error(ccc))/2;
    
    disp([ sprintf('%15s ', char(qsm_names(ccc))), ...
        ' : ', num2str(rmse(ccc), '%-5.6f'), ...
        ' : ', num2str(hfen(ccc), '%5.6f'), ...
        ' : ', num2str(ssim(ccc), '%5.6f'), ...
        '  : ', num2str(wm_error(ccc), '%5.3f'), ...
        '    : ', num2str(gm_error(ccc), '%5.3f'), ...
        '    : ', num2str((gm_error(ccc) + wm_error(ccc))/2 , '%5.6f')]);
end

for k=1:length(X)
    Q{k} = X{k}.*mask_ero;
end

[kFe, pval]=ExtractROIdat_kFe(Q,mask_ero)
