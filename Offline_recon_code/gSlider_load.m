% %%%%%%load gslider_nii data
%  gSlider Reconstruction with B1 and T1 corrections
%  Congyu Liao,PhD. cliao2@mgh.harvard.edu
%$ 05/21/2018
clear
addpath ('library')
addpath ('Nifti_Analyze')
addpath ('imagine')

    disp('Load .nii mag data')
    [f_name, f_path]=uigetfile('*.nii','select mag data');
    f_raw=strcat(f_path,f_name);  
    nifti_info=load_nifti(f_raw);
mag_data=nifti_info.vol;

    disp('Load .nii phase data')
    [f_name, f_path]=uigetfile('*.nii','select phase data');
    f_raw=strcat(f_path,f_name);
    nifti_info=load_nifti(f_raw);
phase_data=nifti_info.vol;
phase_data=(phase_data-2^11)./2^11.*pi;

tps=min(size(mag_data,4),size(phase_data,4));
cplx_data=mag_data(:,:,:,1:tps).*exp(1i*phase_data(:,:,:,1:tps));
% cplx_data=mag_data.*exp(1i*phase_data);
clear nifti_info clear ans clear f_path f_name