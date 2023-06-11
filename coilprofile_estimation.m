function [c] = coilprofile_estimation(ref, coil_img, filtertype, filterparameter)
%%%% ----------------------------------------------
%%%% estimation of coil sensitivity profiles
%%%% inputs:
%%%%        ref -- reference input, eg body coil image or sos of channel
%%%%        images Nfe * Npe
%%%%        coil_img -- channel coil images Nfe * Npe * ports
%%%%        filtertype -- g : gaussian
%%%%                      m : median
%%%%        filterparameter -- if g, sigma
%%%%                        -- if m, window size
%%%%                        -- if wavelet, wavelet level
%%%% outputs:
%%%%        c -- sensitivity profiles Nfe * Npe * ports
%%%% 
%%%% Zhaolin Chen @ HFI, Univ of Melbourne
%%%%
%%%% log: 
%%%%     July 2008 add waveltet
%%%%     Aug 2015 4D matrices
%%%% ----------------------------------------------

[Nfe, Npe] = size(ref);
[ports] = size(coil_img,3);
c = zeros(Nfe,Npe,ports);

switch filtertype
    case 'g'
         ck = c; %kspace ck
         for k=1:ports
             win = gausswin(Nfe,filterparameter)*gausswin(Npe,filterparameter)'; % 10 to 12 is a good balance  NRI data 0.5 -2 for mag 
             c(:,:,k) = coil_img(:,:,k)./ref;      % elimating the phase from bc or sos
             ck(:,:,k) = fftshift(ifft2(fftshift(c(:,:,k)))); % kspace
             [a,b] = find(ck(:,:,k) == max(max(ck(:,:,k))));  % find the center of kspace
             win = circshift(win, [a,b] - [floor(Nfe/2), floor(Npe/2)]); % shift the corresponding win
             c(:,:,k) = fftshift(fft2(fftshift(ck(:,:,k).*win)));
             %c(:,:,k) = c(:,:,k)./max(max(abs(c(:,:,k))));
         end
         
         %%% sensitivity estimation using median filter
    case 'm'
         for k=1:ports
             c(:,:,k) = coil_img(:,:,k)./ref;
             mag = medfilt2(abs(c(:,:,k)), [filterparameter filterparameter]); %magnitude
             pha = angle(c(:,:,k)); % phase
             c(:,:,k) = mag.*exp(sqrt(-1)*pha);
             %c(:,:,k) = c(:,:,k)./max(max(abs(c(:,:,k))));
         end
         
    case 'mg'
          ck = c;
         for k=1:ports
             win = gausswin(Nfe,filterparameter)*gausswin(Npe,filterparameter)'; % 10 to 12 is a good balance  NRI data 0.5 -2 for mag         
             c(:,:,k) = coil_img(:,:,k)./ref;      % elimating the phase from bc or sos
             mag = medfilt2(abs(c(:,:,k)), [10 10]); %magnitude
             ck(:,:,k) = fftshift(ifft2(fftshift(c(:,:,k)))); % kspace
             [a,b] = find(ck(:,:,k) == max(max(ck(:,:,k))));  % find the center of kspace
             win = circshift(win, [a,b] - [floor(Nfe/2), floor(Npe/2)]); % shift the corresponding win
             c(:,:,k) = fftshift(fft2(fftshift(ck(:,:,k).*win)));
             c(:,:,k) = mag.*exp(sqrt(-1)*angle(c(:,:,k)));
             %c(:,:,k) = c(:,:,k)./max(max(abs(c(:,:,k))));
         end
        
    case 'poly'
        for k = 1:ports
            c(:,:,k) = coil_img(:,:,k)./abs(ref);
            c(:,:,k) = img_polyfit(c(:,:,k),filterparameter,10,10);
            %c(:,:,k) = c(:,:,k)./max(max(abs(c(:,:,k))));
        end
        
    case 'wavelet'
        for k = 1:ports
            c(:,:,k) = coil_img(:,:,k)./abs(ref);
        end
        [c]=sensitivity_wavelet(c, 'bior5.5', filterparameter); 
%         for k = 1:ports
%             c(:,:,k) = c(:,:,k)./max(max(abs(c(:,:,k))));
%         end

end
