%% Super Resolution L2+TV method

% Thermal SISR algorithm solves:
%
% argmin (1/2)*(|| SHx - g || + lambda TV(x))
%
% S = down-sampling
% H = Gaussian blur
% lambda = reg. parameter 
% TV = Total Variation
%
% Please consider to cite: 
% "Super-resolution of thermal images using an automatic total variation based method"
% Cascarano P., Corsini F. et. al. "Remote Sensing", 2020.

% if you use this code. The algorithm is described in the aformentioned paper and it is applied to thermal imaging. 




%Input:     
%   A : GT image.
%   g_image :  blurred and downsample image.
%   blur : blur parameters
%   up : up_sampling factor 
%   
%Output:
%   reconstructed_image 

clc;clear all;close all;

addpath('data');
addpath('utils');

%% Parameters Input
up_factor = 4;              %%%% up_sampling factor
fista =1;                   %%%% 1 FISTA (suggested), 0 NO FISTA.
image_syntetic=1;           % 1 if GT is available,  0 if GT is not available
etaM=0;                     % 0 suggested
con_eta=0;                  % 0 suggested

%% Iteration parameters:
outer_iter = 15;            % 10 suggested
inner_iter = 100;           % [50,200] suggested

%% Chambolle parameters:
tau = 1/8;                  % 1/8 suggested
filter_it = 1;              % [1,2,3,4,5] one of this suggested
%% Step-lenght  
beta=[0.002];               % 0.002 suggested
%% Starting lambda parameter
weight =0.005;              % must be tuned to increase the performances

%% Tolerance FBS
Const=0.0000008;
gamma=0.0000001; 

%% Stopping tolerance
criterio_arresto=10^-3;  

%% Altri parametri
scala=2.5*10^-4;        %%only when con_eta=1;
s=1;                    %%only when con_eta=1;

%% GT image if available      
A=importdata('HR_Accursio_target.txt'); A=A - min(min(A)); A = A/max(max(A));
%A=A + 273.15;   %Kelvin degrees 
[m,n]=size(A);
A=modcrop(A,[up_factor,up_factor]);

%% Gaussian parameters
value=0.005;   % if no-blur value<<<1
dimension=15;  % dimension kernel 

% Gaussian assumption
filterGauss=fspecial('gaussian',[dimension dimension],value);

% Low resolution handcrafted 
g_image=imfilter(A,filterGauss,'conv','circular');
g_image=imresize(g_image,1/up_factor,'Bicubic');
%g_image=g_image + 273.15;
[rows,columns]=size(g_image);


% Compute High Resolution Image Size 
Rows = (rows * up_factor); 
Columns = (columns * up_factor);

% Compute Starting Image HR : H^TS^T
U_image = imresize(g_image,up_factor,'Bicubic');%my_upsampling(g_image, rows, columns, up_factor);
U_image = imfilter(U_image, filterGauss,'conv','circular');

U_image_old = U_image;
V_old = U_image;

%% Initialization
lambda(1)=((DiscrepanzaNew_bicubic(U_image, g_image, filterGauss,up_factor))/(2*tlv(U_image,'iso')*weight));
Fx_lambda(1) = funzionaleNew_bicubic(U_image,g_image,filterGauss,up_factor,lambda(1));
newF = Fx_lambda(1);
        
% Algorithm SR

it_globale = 0; 
stopping=100;
outk=1;

while (stopping > criterio_arresto && outk<=outer_iter)                 
       
    tol_in=100;
    it_globale=it_globale+1;
        

    %inner while
    ink=1;
    t_interno(1)=1;
     
    while((tol_in > gamma && (ink <= inner_iter)) ) 
            
        it_globale=it_globale+1;
            
        % Fase di Forward Splitting
        Bg=imfilter(U_image,filterGauss,'conv','circular');
        
        % Downsamplig
        DBg=imresize(Bg,1/up_factor,'Bicubic');
        tmp=DBg-g_image;
        % Upsampling
        update=my_upsampling_bicubic(tmp,rows, columns, up_factor);
        update=update(1:Rows,1:Columns);
        updateB=imfilter(update,filterGauss,'conv','circular');
    
        % Laplacian
        lapU_image=div(grad(U_image));
        V_image=U_image-beta*(updateB+etaM.*lapU_image);
      
        % Backward
        xi = zeros(Rows,Columns,2); 
        
        V_image = chamb(xi,V_image,lambda(outk)*beta,filter_it,tau);
              
        % Fista
        t_interno(ink+1)=(1+sqrt(1+4*t_interno(ink)^2))/2;
              
        if(fista==0)
            U_image=V_image; % NO FISTA 
        else
            U_image=V_image+(t_interno(ink)-1)/t_interno(ink+1)*(V_image-V_old);
        end

        V_old=V_image;
                 
        oldF = newF; 
        newF = funzionaleNew_bicubic(U_image,g_image,filterGauss,up_factor,lambda(outk));
        tol_in = abs(oldF-newF)/newF;
 
        %DEBUG
        ink = ink+1; 
        
    end
      
    fprintf('Num outer interation %i , lambda %f, Fx %f,  Stopping %f \n',outk,lambda(outk),Fx_lambda(outk),stopping);
    
        
    %Lambda update
    Fx_lambda(outk+1)=funzionaleNew_bicubic(U_image,g_image,filterGauss,up_factor,lambda(outk));

        if outk==1
            lambda(outk+1)=lambda(1)*Fx_lambda(2)/Fx_lambda(1); 
        else 
            lambda(outk+1)=lambda(outk)*Fx_lambda(outk)/Fx_lambda(outk-1);
        end
        
        if (image_syntetic)
           reconstructed_image=U_image;
           stopping=abs(Fx_lambda(outk) - Fx_lambda(outk+1))/Fx_lambda(outk);
           err(outk) = sqrt(((norm(A(:)-reconstructed_image(:),2))^2)/prod(size(A)));

        else
            stopping = norm(U_image-U_image_old,'fro')/norm(U_image,'fro');
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

figure;imagesc(g_image);colormap gray; axis equal; axis off; title('Starting');
figure;imagesc(reconstructed_image);colormap gray; axis equal; axis off; title('Rec');