function Fx_lambda = funzionaleNew_bicubic(U_image,g_image,filterGauss,up_factor,lambda)

U_imageBlurred=imfilter(U_image,filterGauss,'conv','circular');
D_U_imageBlurred=imresize(U_imageBlurred,1/up_factor,'Bicubic');

discrepanza=0.5*norm(D_U_imageBlurred-g_image,'fro')^2;
[gx,gy]=gradient(U_image); 
Grad=sqrt(gx.^2+gy.^2);

Fx_lambda=discrepanza+lambda*sum(Grad(:));
end
