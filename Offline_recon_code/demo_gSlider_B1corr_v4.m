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
meanDWI=mean(Img_SuperRes(:,:,:,2:31),4);
% save('meanDWI','meanDWI','-v7.3')
