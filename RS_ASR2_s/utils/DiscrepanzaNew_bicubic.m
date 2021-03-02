function dis = DiscrepanzaNew_bicubic(U_image, g_image, filterGauss,up_factor)

U_imageBlurred=imfilter(U_image,filterGauss,'conv','circular');
D_U_imageBlurred=imresize(U_imageBlurred,1/up_factor);
dis=0.5*norm(D_U_imageBlurred-g_image,'fro')^2;

end