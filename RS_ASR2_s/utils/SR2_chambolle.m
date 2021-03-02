function [reconstructed_image,err,it_globale,time,PSNR]=SR2_chambolle(A, g_image, filterGauss, up_factor, iter, tol)
% Super Resolution - Compressed Sensing-Metodo SR2
% L'algoritmo SR2 risolve il modello di Super-Resolution:
% min (|| SHx - g || + lambda TV(x))
%
% S = down-sampling
% H = blur gaussiano
% lambda = parametro di regolarizzazione
% TV = variazione totale
%
% utilizzando l'algoritmo di Chambolle, l'accelerazione del metodo FISTA 
% e una tecnica adattiva per la scelta del parametro di regolarizzazione
% 
%DESCRIZIONE:
%Input:     
%   A : immagine reale
%   g_image :  immagine sottocampionata manualmente
%   blur : parametri di blur
%   up_factor : fattore di up_sampling
%   iter(1) : iterazioni esterne massime
%   iter(2) : iterazioni interne massime
%   
%Output:
%   reconstructed_image :  immagine ricostruita con fattore up 
%   err : valore di RMSE
%   it_globale : numero di iterazioni effettuate
%   time : tempo
%   PSNR : valore PSNR


tic

fista =1; 
image_syntetic=0;
etaM=0;
con_eta=0;

% iterazioni
outer_iter = iter(1);
inner_iter = iter(2);

% parametri chambolle
tau = 1/8;
filter_it = 1;


beta=[0.02]; 


weight =6; 


Const=0.0000008;
gamma=0.0000001; 


criterio_arresto=tol;  


scala=2.5*10^-4; 
s=1;  

      
[m,n]=size(A);
A=modcrop(A,[up_factor,up_factor]);


[rows,columns]=size(g_image);

channels=1;



Rows = (rows * up_factor); 
Columns = (columns * up_factor);



U_image = my_upsampling(g_image, rows, columns, up_factor);
U_image = imfilter(U_image, filterGauss,'conv','circular');

U_image_old = U_image;

V_old = U_image;


lambda(1)=((DiscrepanzaNew(U_image, g_image, filterGauss,up_factor))/(2*tlv(U_image,'iso')*weight));
Fx_lambda(1) = funzionaleNew(U_image,g_image,filterGauss,up_factor,lambda(1),etaM);
newF = Fx_lambda(1);
        


it_globale = 0; 
stopping=100;
outk=1;

while (stopping > criterio_arresto && outk<=outer_iter)                 
       
        controllo=100;
        it_globale=it_globale+1;
        

        %Innesco il ciclo while interno 
        ink=1;
        t_interno(1)=1;
     
        while((controllo > gamma && (ink <= inner_iter)) ) 
            
            it_globale=it_globale+1;
            
            %% Fase di Forward Splitting
            %Filtraggio con Gaussiana
            Bg=imfilter(U_image,filterGauss,'conv','circular');
            %Downsamplig
            DBg=Bg(2:up_factor:end,2:up_factor:end);
            tmp=DBg-g_image;
            %upsampling
            update=my_upsampling(tmp,rows, columns, up_factor);
            update=update(1:Rows,1:Columns);
            updateB=imfilter(update,filterGauss,'conv','circular');
    
            %Laplaciano
            lapU_image=div(grad(U_image));
            %aggiornamento
            V_image=U_image-beta*(updateB+etaM.*lapU_image);
      
            %% Fase di Backward
            xi = zeros(Rows,Columns,2); 
        
            V_image = chamb(xi,V_image,lambda(outk)*beta,filter_it,tau);
            
            %Proiezione in [0,max(max(A))]           
             V_image=min(max(V_image,min(min(A))), max(max(A)));
              
           %Fista
            t_interno(ink+1)=(1+sqrt(1+4*t_interno(ink)^2))/2;
              
            if(fista==0)
                U_image=V_image; % NO FISTA 
            else
                U_image=V_image+(t_interno(ink)-1)/t_interno(ink+1)*(V_image-V_old);
            end

            V_old=V_image;
                 
            oldF = newF; 
            newF = funzionaleNew(U_image,g_image,filterGauss,up_factor,lambda(outk),etaM);
            controllo = abs(oldF-newF)/newF;
 
        
            ink = ink+1; 
       end
    
       
       Fx_lambda(outk+1)=funzionaleNew(U_image,g_image,filterGauss,up_factor,lambda(outk),etaM);

    if outk==1
          lambda(outk+1)=lambda(1)*Fx_lambda(2)/Fx_lambda(1); 
    else 
          lambda(outk+1)=lambda(outk)*Fx_lambda(outk)/Fx_lambda(outk-1);
    end
       if (image_syntetic)
           reconstructed_image=U_image;
           stopping=abs(Fx_lambda(outk) - Fx_lambda(outk+1))/Fx_lambda(outk + 1);
           err(outk) = sqrt(((norm(A(:)-reconstructed_image(:),2))^2)/prod(size(A)));
  
           
       else
            stopping = norm(U_image-U_image_old,'fro')/norm(U_image,'fro');
            reconstructed_image=U_image;
            err(outk) = sqrt(((norm(A(:)-reconstructed_image(:),2))^2)/prod(size(A)));
       end

if con_eta==1
     if(outk>1)
            [nablaUx,nablaUy]= grad(U_image);
         
            etaM = 1./(1+(nablaUx.^2+nablaUy.^2)/(s^2))*scala;      
           
      end
 end

        U_image_old = U_image;
        outk=outk+1;
end

time=toc;

%% taglio i bordi e calcolo RMSE e PSNR
u_psnr=reconstructed_image(5:end-5,5:end-5);
YH_psnr=A(5:end-5,5:end-5);
RMSE=sqrt((norm(u_psnr(:)-YH_psnr(:),2)^2) /(prod(size(YH_psnr))));
PSNR= 20*log(max(max(YH_psnr))/RMSE);

%% Display High Resolution Image and all Data

 %figure;imagesc(abs(A-reconstructed_image)); colormap(gray); axis equal; axis off; caxis([min(min(abs(A-reconstructed_image))) max(max(abs(A-reconstructed_image)))]); title('Immagine differenza');
 %figure;imagesc(A);colormap(gray);axis equal;axis off;caxis([min(min(A)) max(max(A))]); title('immagine reale HR');
 %figure;imagesc(reconstructed_image);colormap(gray);axis equal;axis off;caxis([min(min(reconstructed_image)) max(max(reconstructed_image))]); title('Immagine ricostruita HR'); 
 %figure;imagesc(g_image);colormap(gray);axis equal;axis off;caxis([min(min(g_image)) max(max(g_image))]); title('Immagine iniziale LR');    
 %figure; loglog(err); title('RMSE al variare dei sottoproblemi');xlabel('itr');ylabel('RMSE');
 %figure;loglog(Fx_lambda);title('funzione obiettivo');

end