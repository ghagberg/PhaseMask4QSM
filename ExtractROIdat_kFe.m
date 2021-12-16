function [kFe, pval]=ExtractROIdat_kFe(X,mask)
%median QSM values in grey matter cortical areas with tissue probability> 98% as defined by SPM12 tissue segmentation using enhanced tissue probability maps
%are extracted from ROIs defined by HardvardOxford cortex atlas, after removal of the cerebellum (SUIT) and subcortical areas (AALrois 71:78)
%median QSM values are extracted from ROIs defined by HardvardOxford sub-cortical atlas (aka Desikan-Killianey), after a 1mm erosion of the Thalamus, Caudate N, Putamen and Globus Pallidus) are  eroded in 3D 

% ***********************************************************
%Method To evaluate the k_Fe parameter (in ppb/microg Iron/mg wet weight
%tissue for the QSMchallenge 2016 data, based on the age-dependent iron content published
%by Hallgren and Sourander in 1958
%
% Confidential until publication
% Hagberg GE, University  and Max-Planck-Institute, Tuebingen, Germany
%
%*************************************************************

mpth=eval('pwd')

%load regions-of-interest from the harvard-oxford (HO) atlas
roi=dir('./kFE/HOctx_HOsc.nii');
DK=spm_read_vols(spm_vol(fullfile(roi.folder,roi.name)));
%1==Prefr Ctx: HO idx=3,4 (superior and middle frontal gyrus)
%2==Temporal Ctx: HOidx=9,10,12,13 (superior and middle temporal gyrus)
%3=Primary sensory Ctx: HOidx=17 (postcentral gyrus)
%4=Parietal Ctx: HOidx=18 19 (Superior Parietal Lobule, anterior division of Supramarginal Gyrus)
%5=Occipital Ctx: HOidx=22 48 (occipital pole, lateral occ. ctx, superior
%part)
%6=Primary motor cortex HOidx=7 (precentral gyrus)
%7=Caudate Nucleus HO subcortex idx=5, 16
%8=Putamen HO subcortex idx=6,17
%9=Globus pallidum HO subcortex id =7, 18

%age of challenge subject is 30y, here is the expected iron content for the
%different ROIs in microgram Iron per g wet weight tissue
finIR=[27.12 29.19 39.74 33.63 40.84 41.21 78.35 106.77 203.41]; %
qidx=length(X)
figure
col=hsv(9);
for tp=1:qidx
    subplot(6,7,tp)
    hold on
    Qmap=X{tp}*1000;
    for k=1:9
        A=(round(DK)==k).*mask;
        y=Qmap(A>0);
        plot(finIR(k)*ones(length(y),1),y,'.','color',col(k,:))
        QS(tp,k)=median(Qmap(A>0));
        QSst(tp,k)=std(Qmap(A>0));
        nvoxROI(tp,k)=sum(A(A>0));
    end
    
    clear  Qmap*
    
    
    
    [B,BINT,R,RINT,STATS] = regress(QS(tp,:)',[finIR' ones(length(finIR),1)]);
    kFe(tp)=B(1);
    pval(tp)=STATS(3);
    subplot(6,7,tp),hold on,plot(finIR,QS(tp,:),'o');
end
fprintf('ready');

%
