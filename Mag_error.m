close all
clc
clear all


load Ir_Chirp2
load Ir_Fourier2
Ir_Chirp2= permute(Ir_Chirp2,[2 1 3 4]);
Ir_Fourier2=permute(Ir_Fourier2,[2 1 3 4]);
% Ir_Chirp2=Ir_Chirp2(:,end:-1:1,:,:);
% Ir_Fourier2=Ir_Fourier2(:,end:-1:1,:,:);

ACC=[4,2,1];     
acc=1;
%%%%%%%
% MODE='Fourier';
% I=Ir_Fourier2(:,:,4,end);
% Ir=Ir_Fourier2(:,:,4,acc);
% %%%%%%

%%%%%%%
MODE='Chirp';
I=Ir_Chirp2(:,:,4,end);
Ir=Ir_Chirp2(:,:,4,acc);
% %%%%%%%


        I = abs(I)/max(abs(I(:)));
        Ir = abs(Ir)/max(abs(Ir(:)));
        Mean_Error= mean(mean(abs(I(50:200,50:200)-Ir(50:200,50:200))))
        [y] = relative_error(I(50:200,50:200), (Ir(50:200,50:200)))
%         
        figure(5); imshow((Ir),[]);colormap gray;  
        FILENAME=strcat(MODE,' Mag Image with acc = ',int2str(ACC(acc)));
        title([MODE, ' Mag Image with acc = ',num2str(ACC(acc))]) 
        saveas(figure(5),[pwd '/SWI_Recon_T2/' FILENAME '.jpg']);
        saveas(figure(5),[pwd '/SWI_Recon_T2/' FILENAME '.fig']);   

%         figure(5); imshow(5*abs(Ir-I),[0 1]);colormap gray; 
%         FILENAME=strcat(MODE,' Mag Error with acc= ',int2str(ACC(acc)));
%         title([MODE, ' Mag Error with acc = ',num2str(ACC(acc))]) 
%         saveas(figure(5),[pwd '/SWI_Recon_T2/' FILENAME '.jpg']);
%         saveas(figure(5),[pwd '/SWI_Recon_T2/' FILENAME '.fig']); 