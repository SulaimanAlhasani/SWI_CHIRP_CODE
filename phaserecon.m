function [pha, phasemask_n,phasemask_p, mag, swi_n,swi_p] = phaserecon(kimg,kimgsos,a,intpl,thr)
% -------------------------------------------
% Usage: calculate phase images from complex inputs
% Inputs: 
%        kimg -- reconstructed kspace signals
%        kimgsos -- sos of square on complex images
%        a    -- gaussian filter size 0 to 20, 10 is good balance
%        intpl -- interplation factor
%        thr -- thresholding: 0 no thresholding
%                             0.05 a good choice
% Outputs
%        pha -- phase images
%        mag -- magnatitue images
%        swi -- phase weighted magnatitue images
%
% Zhaolin Chen @ Howard Florey Institute
% -------------------------------------------


[Nfe,Npe] = size(kimg);

% creat a Gaussina LPF
%win = gausswin(Nfe,10)*gausswin(Npe,8)';
win = gausswin(Nfe,a)*gausswin(Npe,a)';

[a,b] = find(kimg == max(kimg(:)));
win = circshift(win,[a,b]-[floor(Nfe/2), floor(Npe/2)]);

% creat a rectangular LPF
% win = zeros(Nfe,Npe);
% L = 32;W = 32;
% win(round(Nfe/2-L/2):round(Nfe/2+L/2-1),round(Npe/2-W/2):round(Npe/2+W/2-1)) = ones(L,W);


% creat a Kaiser LPF
%win = hann(Nfe)*hann(Npe)';
%win = kaiser(Nfe,10)*kaiser(Npe,100)';

img = fftshift(fft2(kimg,intpl*Nfe,intpl*Npe));
imgsos = fftshift(fft2(kimgsos,intpl*Nfe,intpl*Npe));
img_lpf = fftshift(fft2(kimg.*win,intpl*Nfe,intpl*Npe));


img_hpf = img ./ img_lpf;

mag = abs(img);

thd = (thr/sqrt(sqrt(intpl)))*max(abs(imgsos(:))); %%thr = 0.05


pha = angle(img_hpf);


if (thr ~= 0)
   for i = 1:intpl*Nfe
       for j = 1:intpl*Npe
           if (abs(imgsos(i,j)) <= thd)
              pha(i,j) = pha(i,j)/10;
           end
       end
   end
end

%mag = sqrt(mag);

if pha >=0
   phasemask_n = 1;
else
   phasemask_n = (pha+pi)./pi;
end

if pha <=0
   phasemask_n = 1;
else
   phasemask_p = (pi-pha)./pi;
end

swi_n = phasemask_n.^3.*mag;
swi_p = phasemask_p.^3.*mag;

