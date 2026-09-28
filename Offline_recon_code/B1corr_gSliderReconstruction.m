function Img_SuperRes = B1corr_gSliderReconstruction(Img_lowRes,A_recon,lambdaTikPercent)

[nx,ny,nz_l,ndir,nrf]=size(Img_lowRes);
n_basis = size(A_recon,1)/size(A_recon,2);
Img_SuperRes = zeros(nx,ny,nz_l*nrf/n_basis,ndir);
for x = 1:nx
    for y = 1:ny
        A_recon_current=squeeze(A_recon(:,:,x,y));
        A_tik = InverseAmatrix(A_recon_current,lambdaTikPercent);
%         vox=squeeze(Img_lowRes(x,y,:,:,:));
%         rhs=permute(vox, [3,1,2]);
%         rhs=reshape(rhs, size(rhs,1)*size(rhs,2),size(rhs,3));
%         res=A_tik *rhs;
%         Img_SuperRes(x,y,:,:)=res;
        for diff_dir = 1:ndir  
            vox= squeeze(Img_lowRes(x,y,:,diff_dir,:));
            rhs = permute(vox, [2,1]);
            rhs = rhs(:);
            res = A_tik * rhs;            
            Img_SuperRes(x,y,:,diff_dir) = res;
        end       
    end
    disp(['SuperRes: nx=', num2str(x), ' / ', num2str(nx)])
 
end




end

