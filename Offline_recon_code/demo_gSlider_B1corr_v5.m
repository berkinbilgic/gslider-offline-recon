%  gSlider Reconstruction with B1 and T1 corrections
%  Congyu Liao,PhD. cliao2@mgh.harvard.edu
%$ 05/21/2018
%%
addpath ('library')
addpath ('imagine')
disp('Loading data ....');
B1_index=B1_index*Voltfac;
nrf_d=5;
convGrappa_mag=cplx_data;
complex_data_conv=reshape(convGrappa_mag,[size(convGrappa_mag,1),size(convGrappa_mag,2),size(convGrappa_mag,3),nrf_d,size(convGrappa_mag,4)/nrf_d]);
complex_data_conv=permute(complex_data_conv,[1 2 3 5 4]);
%%

PF = 1;   % for 960um PF is 7/8 partial fourier
ExtraCropPElines = 10; % number of extra k-line to not trust next to p.f. zero filled space, due to grappa smoothing around this region. 
POCS_flag = 0; % do POCS for p.f.
avgs = 1;
resolution =0.86;% 1.00;
%% complex_data: x,y,z,diffdir,RFenc,
complex_data = complex_data_conv;
disp('Load slice profile')
    [s_name, s_path]=uigetfile('*.mat','select slice profile');
    s_raw=strcat(s_path,s_name); 
    load (s_raw)

% clear phase_data mag_data
% use phase of filtered low Res image as estimate for background phase and remove it. (reasonable for low bvalues)
Img_lowRes = permute(RealDiffusion_lowRes(permute(complex_data,[2,1,3,4,5]),PF,ExtraCropPElines,POCS_flag),[2,1,3,4,5]);

%%
if nrf_d==10
    [n1,n2,nz,ndir0,nrf]=size(Img_lowRes);
    Img_lowRes = reshape(permute(reshape(Img_lowRes,[n1,n2,nz,2,ndir0/2,nrf]),[1,2,3,5,6,4]),[n1,n2,nz,ndir0/2,nrf*2]);
end
[nx,ny,nz_l,ndir,nrf_d]=size(Img_lowRes);
%% load slice profile 


nrf=nrf_d;

mxy=MOut;
z=zOut;
B1=B1map_de;
clear B1map B1map_intp B1map MOut zOut
% SliderShift = [0 0 0 0 0];
SliderShift = zeros(1, nrf);
SidelobesSlices = 6; 
CropEdge = 1;

[A_encoding,psf_matrix] = B1corr_SliderPSFmatrix_pShift_weightDownEdges(mxy,z,nz_l,nrf,SliderShift,SidelobesSlices,CropEdge,B1,B1_index);

%%
lambdaTikPercent=0.1;
Img_SuperRes = B1corr_gSliderReconstruction(Img_lowRes,A_encoding,lambdaTikPercent);
PhaseFac=pi/2;

Img_SuperRes = real(Img_SuperRes*exp(1i*PhaseFac));
%%
Img_SuperRes = permute(Img_SuperRes,[2,1,3,4]);
%%
disp('Done');
% meanDWI=mean(Img_SuperRes(:,:,:,2:31),4);
% save('meanDWI','meanDWI','-v7.3')
%% save Img_SuperRes as nifti, header from the gSlider magnitude
% in-plane: undo the x/y swap above on a copy, so voxel axes match the magnitude header.
% slices: slab s fills slices (s-1)*nrf+1:s*nrf (B1corr_SliderPSFmatrix_pShift_weightDownEdges),
% so spacing = slab spacing/nrf and slice 1 is centred (nrf-1)/(2*nrf) slab below slab 1.
% assumes contiguous slabs (thickness = spacing), as the reconstruction does.
% float32, values as reconstructed (real, incl. negatives); does not overwrite.
disp('Load .nii mag data (header for Img_SuperRes)')
[f_name, f_path]=uigetfile('*.nii','select mag data');
nii=load_untouch_nii(strcat(f_path,f_name),1);
[~,f_base]=fileparts(f_name);
f_out=fullfile(f_path,[f_base '_SuperRes.nii']);
if exist(f_out,'file'), error('%s already exists, not overwriting',f_out); end
img_out=single(permute(Img_SuperRes,[2,1,3,4]));
if ~isequal([size(img_out,1) size(img_out,2) nz_l],double(nii.hdr.dime.dim(2:4))), error('Img_SuperRes does not match the grid of %s',f_name); end
nsub=size(img_out,3)/nz_l; % slices per slab (= nrf)
k0=-(nsub-1)/(2*nsub); % centre of slice 1, in slab voxel units
hdr=nii.hdr;
hdr.dime.dim=[4 size(img_out,1) size(img_out,2) size(img_out,3) size(img_out,4) 1 1 1];
hdr.dime.datatype=16; hdr.dime.bitpix=32; hdr.dime.scl_slope=1; hdr.dime.scl_inter=0;
% qform: rotation unchanged, shift origin by k0 slabs along the slice axis
qb=hdr.hist.quatern_b; qc=hdr.hist.quatern_c; qd=hdr.hist.quatern_d; qa=sqrt(max(0,1-qb^2-qc^2-qd^2));
qfac=hdr.dime.pixdim(1); if qfac==0, qfac=1; end
q_shift=[2*(qb*qd+qa*qc); 2*(qc*qd-qa*qb); qa*qa+qd*qd-qc*qc-qb*qb]*qfac*hdr.dime.pixdim(4)*k0;
hdr.hist.qoffset_x=hdr.hist.qoffset_x+q_shift(1); hdr.hist.qoffset_y=hdr.hist.qoffset_y+q_shift(2); hdr.hist.qoffset_z=hdr.hist.qoffset_z+q_shift(3);
hdr.dime.pixdim(4)=hdr.dime.pixdim(4)/nsub;
% sform: same shift, slice column scaled by 1/nsub
srow=[hdr.hist.srow_x; hdr.hist.srow_y; hdr.hist.srow_z];
srow(:,4)=srow(:,4)+k0*srow(:,3); srow(:,3)=srow(:,3)/nsub;
hdr.hist.srow_x=srow(1,:); hdr.hist.srow_y=srow(2,:); hdr.hist.srow_z=srow(3,:);
nii.hdr=hdr; nii.img=img_out;
save_untouch_nii(nii,f_out);
disp(['Saved ' f_out]);
clear nii hdr img_out srow qa qb qc qd qfac q_shift nsub k0 f_base f_out f_name f_path
