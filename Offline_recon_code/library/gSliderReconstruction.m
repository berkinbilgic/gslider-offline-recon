function Img_SuperRes = gSliderReconstruction(Img_lowRes,A_recon,lambdaTikPercent)

A_tik = InverseAmatrix(A_recon,lambdaTikPercent);
[nx,ny,nz_l,ndir,nrf]=size(Img_lowRes);
n_basis = size(A_recon,1)/size(A_recon,2);

Img_SuperRes = zeros(nx,ny,nz_l*nrf/n_basis,ndir);

for diff_dir = 1:ndir
    disp(['SuperRes: ', num2str(diff_dir), ' / ', num2str(ndir)])
    
    Img_lowResCurrentGroup = squeeze(Img_lowRes(:,:,:,diff_dir,:));
    Img_SuperResCurrent = zeros(nx,ny,nrf*nz_l/n_basis);
   
    for x = 1:nx
        for y = 1:ny
            vox = squeeze(Img_lowResCurrentGroup(x,y,:,:));
            rhs = permute(vox, [2,1]);
            rhs = rhs(:);
            res = A_tik * rhs;            
            Img_SuperResCurrent(x,y,:) = res;
        end
    end
    Img_SuperRes(:,:,:,diff_dir) = Img_SuperResCurrent;
end
end

