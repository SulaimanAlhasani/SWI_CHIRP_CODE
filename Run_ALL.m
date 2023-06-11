
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc ; clear all; close all
tic 
addpath(strcat(pwd,'/Rice_Wavelet_Toolbox_2.4'));
addpath(strcat(pwd,'/irt'));
addpath(strcat(pwd,'/DATA'));
addpath(strcat(pwd,'/Code'));
setup
addpath(strcat(pwd,'/Wavelab850'));
WavePath
N=256;
%%
load DATA_Chirp; 
load DATA_Fourier; 

%%
acc_factor=256./[1];
for acc=1:length(acc_factor)
    s=0;
for  slice=1;%slice=[30,33,37,40,43,44,60,67]
    s=s+1;
 Ir_Chirp(:,:,s,acc)=T2_of_SparseMRI_Parallel_Chirp_invivo_SWI(squeeze(DATA_Chirp(:,:,:,slice)),acc_factor(acc));
 Ir_Fourier(:,:,s,acc)=T2_of_SparseMRI_Parallel_Fourier_invivo_SWI(squeeze(DATA_Fourier(:,:,:,slice)),acc_factor(acc));
%  Ir_Noiselet(:,:,s,acc)=T2_of_SparseMRI_Parallel_Noiselet_invivo_SWI(squeeze(DATA_Noiselet(:,:,:,slice)),acc_factor(acc));
end 
end 

%%
save Ir_Chirp_full_sampled.mat Ir_Chirp
save Ir_Fourier_full_sampled.mat Ir_Fourier


toc

