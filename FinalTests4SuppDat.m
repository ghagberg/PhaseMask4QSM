clear all
load v
for k=1:41
X{k}=spm_read_vols(v(k));
end

N=v(1).dim;
mask_ero=ones(N);
for k=1:length(X)    
    mask_ero= mask_ero.*(abs(X{k}*1000)>eps);
end

for k=1:length(X)
    qsm_results(:,:,:,k) = X{k}.*mask_ero;
end

qsm_names = {'chi33';'chi_cosmos';'L2';'iLSQR';'MEDI';'L2rshp MM';'L2rshpPB';'L2rshporig';'L2vshp MM';'L2vshpPB';'L2vshpOrig';'iLSQRrshpMM';'iLSQRrshpPB';'iLSQRrshpOrig';'iLSQRvshpMM';'iLSQRvshpPB';'iLSQRvshpOrig';'MEDIrshpMM';'MEDIrshpPB';'MEDIrshpOrig';'MEDIvshpMM';'MEDIvshpPB';'MEDIvshpOrig'};
for k=24:41
    qsm_names{k}=['ro' qsm_names{k-18}];
end
qsm_names{42}='chi_cosmos'
%%-------------------------------------------------------------------------
%% compute and print performance metrics
%%-------------------------------------------------------------------------
chi_33=X{1}.*mask_ero;
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

T=table(qsm_names,rmse',hfen',ssim',mean_error',kFe',pval','VariableNames',{'qsm_names' 'rmse' 'hfen' ' ssim' 'mean_error' 'kFe' 'pval'})
writetable(T,'TableSupplDat.xls')
