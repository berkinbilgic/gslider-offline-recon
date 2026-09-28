function [ A_mtxb, psf_matrix] = B1corr_SliderPSFmatrix_pShift_weightDownEdges(mxy,z,nAcqSlc,SliderFac,SliderShift,SidelobesSlices,CropEdge,B1,B1_index)

subgroups_extended = SliderFac + SidelobesSlices; % (5+6=11) include within slide and sidelobes (sidelobes has to be even number)
psf_matrix = zeros(length(SliderShift),subgroups_extended,size(mxy,3));%size(5,11,b1)
group_distance = 1/SliderFac; %1mm/5=0.2mm
Zlength = subgroups_extended*group_distance;% 11*0.2mm = 2.2mm
SliderShiftLength = SliderShift*group_distance; %[0,0,0,0,0];
for iB1=1:size(mxy,3)
    z_b1=z(:,:,iB1);
    mxy_curb1=squeeze(mxy(:,:,iB1));
    for iRF=1:length(SliderShift)
        mxy_b1{iRF}=mxy_curb1(:,iRF);
    end
    for ShiftCount = 1:length(SliderShift)   
        SliderShiftLength_current = SliderShiftLength(ShiftCount); % 0 
        z_cut = z_b1(z_b1>=-Zlength/2+SliderShiftLength_current & z_b1<=Zlength/2+SliderShiftLength_current); % slice profile has been normalized to one so set width here as 2
        mxy_cut = mxy_b1{ShiftCount}(z_b1>=-Zlength/2+SliderShiftLength_current & z_b1<=Zlength/2+SliderShiftLength_current);
        PSF = [];
        for GroupCount = 1:subgroups_extended
            StartZ = -Zlength/2 + SliderShiftLength_current + group_distance *(GroupCount-1);
            EndZ = -Zlength/2 + SliderShiftLength_current + group_distance *(GroupCount);
            PSF(GroupCount) = sum(mxy_cut(z_cut>=StartZ & z_cut<=EndZ))./sum(z_cut>=StartZ & z_cut<=EndZ); %/length(mxy_cut(z_cut>=StartZ & z_cut<=EndZ));
            %sum(z_cut>=StartZ & z_cut<=EndZ)
        end

        psf_matrix(ShiftCount,:,iB1) = PSF; 
    end
end


% Flip angle is slightly lowered than simulated and the edge of the psf
% need slight modification as follow:
if (0)
psf_matrix(:,4) = psf_matrix(:,4)*1.04;
psf_matrix(:,8) = psf_matrix(:,8)*1.04;
psf_matrix(:,3) = psf_matrix(:,3)*0.2;
psf_matrix(:,9) = psf_matrix(:,9)*0.2;
end

psf_matrix(:,4,:) = psf_matrix(:,4,:)*0.90;
psf_matrix(:,8,:) = psf_matrix(:,8,:)*0.90;


psf_matrix = round(psf_matrix*1000)/1000;
nx=size(B1,1);ny=size(B1,2);
A_mtxb = zeros(nAcqSlc * SliderFac,nAcqSlc * SliderFac + (subgroups_extended-SliderFac),nx,ny);
for ix=1:nx
    for iy=1:ny     
        for count = 1:nAcqSlc   
            B1_col=B1(ix,iy,count);
            B1_count=find(B1_index==B1_col);
            beginCol = 1+ (count-1)*SliderFac;
            endCol = beginCol+ subgroups_extended -1; 
            A_mtxb(1+(count-1)*SliderFac:count*SliderFac,beginCol:endCol,ix,iy) = psf_matrix(:,:,B1_count);
        end
    end
end

if CropEdge == 1
   A_mtxb = A_mtxb(:,1+SidelobesSlices/2:end-SidelobesSlices/2,:,:);
end





