
function [im_res]=T2_of_SparseMRI_Parallel_Fourier_invivo_SWI(data,M)

%  data = double(data./max(abs(data(:))));
data = double(data);
N=256;
No = FFTMTX(N);
[r c ch] = size(data); 
for k =1:ch
  %   Is(:,:,k) = fftshift(fft2(ifftshift(data(:,:,k))));
    Is(:,:,k) = fftshift(fft(ifftshift(data(:,:,k),2),[],2),2);
    Is(:,:,k) = No'*squeeze(Is(:,:,k)); k
end
%%% Calculating Sensitivity maps %%%%%%
imgss = (sum(Is(:,:,:),3)); %not sum-of-sqrt
imgsos_mag = sqrt(sum(abs(Is(:,:,:)).^2,3));
I = imgsos_mag.*exp(sqrt(-1)*angle(imgss));
S=coilprofile_estimation(I,Is,'mg',5);
S = mri_sensemap_denoise(S,'bodycoil', ones(size(I)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[128=3, 64=8,  32=20, 16=200] 
if M==64
        RED=8;
    elseif M==128
        RED=3;
end 
    
%%%% Fourier encoding %%%%%%%%%%%%%%%%
F = generate_fourier(N);
F = FFTMTX(N);
% mask=Sampling_Fourier(N,M/N);


if M==256
    mask=ones(256);
else 
    red=RED
    mask = generate_mask(N,N,red);
    r1 = size(find(mask(:,1)),1)
    while ( r1~=M)
        mask = generate_mask(N,N,red);
        r1 = size(find(mask(:,1)),1)
    end
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
r = find(mask(:,1));
Fu = F(r,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%generate Fourier sampling operator
A = @(x) FT_coil(x, S, Fu);
A_= @(x) A_FT_coil (x, S, Fu);
FT = A_operator(@(x) A(x), @(x) A_(x));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1 Recon Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data = FT*I; 
data =[];
for k = 1:ch
    data = [data; Fu*Is(:,:,k)];
end
    
TVWeight = 0.002; 	% Weight for TV penalty
xfmWeight = 0.005;	% Weight for Transform L1 penalty
Itnlim = 50;		% Number of iterations

% scale data
im_dc = FT'*(data);
data = data/max(abs(im_dc(:)));
im_dc = im_dc/max(abs(im_dc(:)));


%generate transform operator
% wav = daubcqf(4);
% W1 = @(x) idwtcplx(x,wav);
% WT1 = @(x) dwtcplx(x,wav);
% XFM = A_operator(@(x) WT1(x), @(x) W1(x));
%  XFM=1;
XFM = Wavelet('Daubechies',4,4);	% Wavelet

% initialize Parameters for reconstruction
param = init;
param.FT = FT;
param.XFM = XFM;
param.TV = TVOP;
param.data = data;
param.TVWeight =TVWeight;     % TV penalty 
param.xfmWeight = xfmWeight;  % L1 wavelet penalty
param.Itnlim = Itnlim;


res = XFM*im_dc;

% figure; 
% subplot(2,2,1);imshow(log(abs(XFM*I)),[]); drawnow
% subplot(2,2,2);imshow(log(abs(res)),[]); drawnow
% subplot(2,2,3);imshow((abs(I)),[]); drawnow
% subplot(2,2,4);imshow((abs(im_dc)),[]); drawnow
% do iterations
% tic
for n1=1:10
    n1
	res = fnlCg(res,param);
	im_res = XFM'*res;
%     figure(n1);
% 	  subplot(2,2,1);imshow(log(abs(XFM*I)),[]); drawnow
%     subplot(2,2,2);imshow(log(abs(res)),[]); drawnow
%     subplot(2,2,3);imshow((abs(I)),[]); drawnow
%     subplot(2,2,4);imshow((abs(im_res)),[]); drawnow
%     pause
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I = abs(I)/max(abs(I(:)));
% Ir = abs(im_res)/max(abs(im_res(:)));
% error = relative_error(I,Ir)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
% MODE='Fourier';
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
% 
% % figure; imshow(abs(Ir),[]);
% toc



