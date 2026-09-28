function [img_cropfilt, img_crop] = lowResAndFilter(img,Crop,FilterRatio)

 kdata = mrir_fDFT_freqencode(mrir_fDFT_phasencode(img));
 
 
 RO_CropLength = round(size(kdata,1)*Crop);
 RO_CropIndex = [-floor(RO_CropLength/2):floor(RO_CropLength/2)] + floor(size(kdata,1)/2)+1;
 PE_CropLength = round(size(kdata,2)*Crop);
 PE_CropIndex = [-floor(PE_CropLength/2):floor(PE_CropLength/2)] + floor(size(kdata,2)/2)+1;
 kdata_crop = kdata(RO_CropIndex,PE_CropIndex,:,:,:);
 img_crop = mrir_iDFT_freqencode(mrir_iDFT_phasencode(kdata_crop));
 
 
 if FilterRatio >0
     % assume data is same resolution in PE and RO but FOV can be different.- do same amount of Hamming in RO -so get higher SNR phase estimate by
     % smooting this way as well
     
     [Nf, Np, Nslc, Ndiff, Ns]  = size(kdata_crop);
     
     hamming_filterPE  = zeros(Nf, Np);
     PE_filterLength = Np*FilterRatio;
     h_PE = hamming(PE_filterLength);
     PE_filterIndex = [-floor(Np/2):floor(Np/2)] + floor(PE_filterLength/2)+1;
     hamming_filterPE = repmat(h_PE(PE_filterIndex).', [Nf, 1]);
     
     hamming_filterRO  = zeros(Nf, Np);
     RO_filterLength = Nf*FilterRatio;
     h_RO = hamming(RO_filterLength);
     RO_filterIndex = [-floor(Nf/2):floor(Nf/2)] + floor(RO_filterLength/2)+1;
     hamming_filterRO = repmat(h_RO(RO_filterIndex), [1, Np]);
     
     hamming_filter =  hamming_filterPE.* hamming_filterRO;
     img_cropfilt    = mrir_iDFT_freqencode(mrir_iDFT_phasencode(repmat(hamming_filter,[ 1 1 Nslc Ndiff Ns]).*kdata_crop));
 else
     img_cropfilt = img_crop;
 end