%%%%%load B1 maps
%  gSlider Reconstruction with B1 and T1 corrections
%  Congyu Liao,PhD. cliao2@mgh.harvard.edu
%$ 05/21/2018
% clear
flag_b1mask=0;
Voltfac=1.0;
B1_index=[0.7:0.05:1.2];

nxi=max(size(cplx_data,1),size(cplx_data,2));nyi=nxi;nzi=size(cplx_data,3);
nxv=nxi;nyv=nyi;nzv=nzi;

addpath ('library')
addpath ('Nifti_Analyze')
addpath ('imagine')
addpath ('tools')
    disp('Load .nii B1 data')
    [f_name, f_path]=uigetfile('*.nii','select B1 data');
    f_raw=strcat(f_path,f_name);  
    nifti_info=load_nifti(f_raw);
B1_data=nifti_info.vol;

%%

B1_vox2ras = nifti_info.vox2ras; % keep B1 geometry before nifti_info is reused for the mask


%%

if(flag_b1mask==1)
    disp('Load .nii B1_mask data')
    [f_name, f_path]=uigetfile('*.nii.gz','select mask data');
    f_raw=strcat(f_path,f_name);  
    nifti_info=load_nifti(f_raw);
    mask=nifti_info.vol;
    
    B1_data=B1_data.*mask;
end
B1_data_even=B1_data(:,:,1:floor(size(B1_data,3)/2));
B1_data_odd=B1_data(:,:,floor(size(B1_data,3)/2)+1:end);
[nx,ny,nz]=size(B1_data);
B1map=zeros(size(B1_data));
B1map(:,:,1:2:end)=B1_data_odd;
B1map(:,:,2:2:end)=B1_data_even;

%%

if(flag_b1mask==0)
    disp('Load .nii B1_mask data')
    [f_name, f_path]=uigetfile('*.nii.gz','select mask data');
    f_raw=strcat(f_path,f_name);  
    nifti_info=load_nifti(f_raw);
    mask=nifti_info.vol;
    mask=mask(:,:,:,1);
    mask=crop(mask,[nxv,nyv,nzv]);
end

 x=linspace(-1,1,nx);
 y=linspace(-1,1,ny);
 z=linspace(-1,1,nz);
 [X,Y,Z] = meshgrid(x,y,z);
 
 xi = linspace(-1,1,nxi);
 yi = linspace(-1,1,nyi);
 zi = linspace(-1,1,nzi);
 [Xi,Yi,Zi] = meshgrid(xi,yi,zi);
 B1map_intp = interp3(X,Y,Z,B1map,Xi,Yi,Zi,'cubic');
%  B1map_intp = crop(B1map_intp,[size(cplx_data,1),size(cplx_data,2),size(cplx_data,3)]);



B1map_intp=B1map_intp./10/80;
% B1map_intp=B1map_intp./53;

% bb hack
% resample B1 onto the gSlider grid in world coordinates (sform vox2ras of B1 and of the
% gSlider mask, which shares the gSlider geometry; assumes flag_b1mask==0), instead of
% imresize3, which ignores the 220 vs 224 mm FOV and FOV-center differences.
% linear interpolation, nearest-value extrapolation beyond the B1 FOV edge
[vi,vj,vk] = ndgrid(0:size(cplx_data,1)-1, 0:size(cplx_data,2)-1, 0:size(cplx_data,3)-1);
vb = B1_vox2ras \ nifti_info.vox2ras * [vi(:)'; vj(:)'; vk(:)'; ones(1,numel(vi))];
B1_interp = griddedInterpolant(B1_data, 'linear', 'nearest');
B1_data_resize = reshape(B1_interp(vb(1,:)'+1, vb(2,:)'+1, vb(3,:)'+1), size(vi));
clear vi vj vk vb B1_interp B1_vox2ras
B1map_intp = B1_data_resize .* mask / 800;

if (flag_b1mask==1)
    B1map_intp(find(B1map_intp<B1_index(1)))=1.*Voltfac;
end

clear B1map_de
for ii=1:size(B1map_intp,1)
    for jj=1:size(B1map_intp,2)
        for kk=1:size(B1map_intp,3)
            temp_a=abs(squeeze(B1map_intp(ii,jj,kk))-B1_index);
            temp_b=find(temp_a(:)==min(temp_a(:)));
            B1map_de(ii,jj,kk)=B1_index(temp_b(1));
        end
    end
end
B1map_de=B1map_de.*Voltfac;
if (flag_b1mask==0)
    B1map_de=B1map_de.*mask;
end

B1map_de(find(B1map_de==0))=1.*Voltfac;
B1map_de=crop(B1map_de,[nxi,nyi,nzi]);



clear nifti_info 
clear ans 
clear f_path f_name
clear B1_data_even B1_data_odd  B1_data ii jj kk nx ny nz 
clear x X xi Xi y Y yi Yi z Z zi Zi
clear temp_a temp_b
clear nxi nyi nzi nxv nyv nzv