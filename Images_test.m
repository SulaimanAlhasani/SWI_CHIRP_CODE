clc
 close all
clear all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load Ir_Chirp_full_sampled.mat; 
load Ir_Fourier_full_sampled.mat; 
Ir_Chirp_Full=Ir_Chirp;
Ir_Fourier_Full=Ir_Fourier;

load Ir_Chirp_total.mat; 
load Ir_Fourier_total.mat; 
Ir_Chirp(:,:,:,3)=Ir_Chirp_Full;
Ir_Fourier(:,:,:,3)=Ir_Fourier_Full;
Ir_Chirp = permute(Ir_Chirp,[2 1 3 4]);
Ir_Fourier=permute(Ir_Fourier,[2 1 3 4]);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ACC=[4,2,1];
SLICE=1;
for s=1
    for acc=1
        %%%%%%%%%%%%%
        MODE='Chirp';
        I=Ir_Chirp(:,:,:,end);
        kimg_full=fft2c(I);
        kimgsos_full=fft2c(abs(I));
        [pha_full, phasemask_n,phasemask_p, mag_full, swi_n_full,swi_p_full] = phaserecon(kimg_full,kimgsos_full,2,1,0);
        im_res=Ir_Chirp(:,:,s,acc);
        %%%%%%%%%%%%%%
        MODE='Fourier';
        I=Ir_Fourier(:,:,:,end);
        kimg_full=fft2c(I);
        kimgsos_full=fft2c(abs(I));
        [pha_full, phasemask_n,phasemask_p, mag_full, swi_n_full,swi_p_full] = phaserecon(kimg_full,kimgsos_full,2,1,0);
        im_res=Ir_Fourier(:,:,s,acc);
%         %%%%%%%%%%%%%%
        % MODE='Noiselet';
        % I=Ir_Noiselet_Full;
        %kimg_full=fft2c(I);
        %kimgsos_full=fft2c(abs(I));
        %[pha_full, phasemask_n,phasemask_p, mag_full, swi_n_full,swi_p_full] = phaserecon(kimg_full,kimgsos_full,2,2,0);
        % im_res=Ir_Noiselet2(:,:,s,acc);
        %%%%%%%%%%%%%%%%%
        kimg=fft2c(im_res);
        kimgsos=fft2c(abs(im_res));
        %%
        [pha, phasemask_n,phasemask_p, mag, swi_n,swi_p] = phaserecon(kimg,kimgsos,2,1,0);
        
%         figure(5); imshow(mag,[]);colormap gray;  axis xy;
%         FILENAME=strcat(MODE,' Mag Image with acc = ',int2str(ACC(acc)));
%         title([MODE, ' Mag Image with acc = ',num2str(ACC(acc))]) 
%         saveas(figure(5),[pwd '/SWI_Recon_T2/' FILENAME '.jpg']);
%         saveas(figure(5),[pwd '/SWI_Recon_T2/' FILENAME '.fig']);      
        
%         figure(5); imshow(swi_p,[]);colormap gray; axis xy;
%         FILENAME=strcat(MODE,'SWIp with acc = ',int2str(ACC(acc)));
%         title([MODE, ' SWIp with acc = ',num2str(ACC(acc))]) 
%         saveas(figure(5),[pwd '/SWI_Recon_T2/' FILENAME '.jpg']);
%         saveas(figure(5),[pwd '/SWI_Recon_T2/' FILENAME '.fig']);  
        
%         figure(5); imshow(swi_n,[]);colormap gray; axis xy;
%         FILENAME=strcat(MODE,' SWIn with acc = ',int2str(ACC(acc)));
%         title([MODE, ' SWIn with acc = ',num2str(ACC(acc))]) 
%         saveas(figure(5),[pwd '/SWI_Recon_T2/' FILENAME '.jpg']);
%         saveas(figure(5),[pwd '/SWI_Recon_T2/' FILENAME '.fig']);        
%        
        figure(5); imagesc(pha,[-0.4 0.4]) ;colormap gray;colorbar; axis xy;
        FILENAME=strcat(MODE,' Phase Image with acc = ',int2str(ACC(acc)));
        title([MODE, ' Phase Image with acc = ',num2str(ACC(acc))]) 
        saveas(figure(5),[pwd '/SWI_Recon_T2/' FILENAME '.jpg']);
        saveas(figure(5),[pwd '/SWI_Recon_T2/' FILENAME '.fig']); 
      VAR1= var(var(abs(pha(70:170,70:170)-pha_full(70:170,70:170))))
      [y] = relative_error(pha_full(70:170,70:170), pha(70:170,70:170))
       [y1] = relative_error(abs(I(70:170,70:170)), abs(im_res(70:170,70:170)))
%         figure(5); imagesc((pha-pha_full),[-0.4 0.4]);colormap gray;colorbar;  axis xy;
%         FILENAME=strcat(MODE,' Phase Error with acc= ',int2str(ACC(acc)));
%         title([MODE, ' Phase Error with acc = ',num2str(ACC(acc))]) 
%         saveas(figure(5),[pwd '/SWI_Recon_T2/' FILENAME '.jpg']);
%         saveas(figure(5),[pwd '/SWI_Recon_T2/' FILENAME '.fig']); 
%         
%         I = abs(I)/max(abs(I(:)));
%         Ir = abs(im_res)/max(abs(im_res(:)));
%         
%         figure(5); imshow(Ir,[]);colormap gray;
%         FILENAME=strcat(MODE,' Mag Image with acc = ',int2str(ACC(acc)));
%         title([MODE, ' Mag Image with acc = ',num2str(ACC(acc))]) 
%         saveas(figure(5),[pwd '/SWI_Recon_T2/' FILENAME '.jpg']);
%         saveas(figure(5),[pwd '/SWI_Recon_T2/' FILENAME '.fig']);   
%         
%         
%        
%         
%         
%         
%         figure(5); imshow(4*abs(Ir-I),[]);colormap gray; 
%         FILENAME=strcat(MODE,' Mag Error with acc= ',int2str(ACC(acc)));
%         title([MODE, ' Mag Error with acc = ',num2str(ACC(acc))]) 
%         saveas(figure(5),[pwd '/SWI_Recon_T2/' FILENAME '.jpg']);
%         saveas(figure(5),[pwd '/SWI_Recon_T2/' FILENAME '.fig']); 
         
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I = abs(I)/max(abs(I(:)));
% Ir = abs(im_res)/max(abs(im_res(:)));
% error = relative_error(I,Ir)
% 
% figure;
% subplot(1,2,2); imshow(abs(Ir),[]);
% title('Reconstructed Image'); 
% subplot(1,2,1); imshow(abs(I),[]);toc
% title('Orignal Image'); 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
% MODE='Chirp';
% figure (5)
% imshow(abs(Ir),[]);
% FILENAME=strcat(MODE,'_90d data with acc = ',int2str(N/M));
% title([MODE, ' data with acc = ',num2str(N/M),',  TVWeight = ',num2str(TVWeight),' xfmWeight = ',num2str(xfmWeight)]) 
% saveas(figure(5),[pwd '/Chirp72_Recon_T2/' FILENAME '.jpg']);
% saveas(figure(5),[pwd '/Chirp72_Recon_T2/' FILENAME '.fig']);
% 
% figure(6); imshow(4*abs(I-Ir),[0 1]);
% FILENAME=strcat(MODE, ',ERROR = ',num2str(error),', acc = ',num2str(N/M));
% title([MODE, ',ERROR = ',num2str(error),', acc = ',int2str(N/M)]) 
% saveas(figure(6),[pwd '/Chirp72_Recon_T2/' FILENAME '.jpg']);
% saveas(figure(6),[pwd '/Chirp72_Recon_T2/' FILENAME '.fig']);

    end 
end 

