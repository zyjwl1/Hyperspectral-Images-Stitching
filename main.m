close all; clear all;
run('vlfeat-0.9.14/toolbox/vl_setup');
addpath('transforms');
imfolder = 'examples/remote';
im_n =18;
imfile = cell(im_n,1);
for ii = 1:im_n
    imfile{ii} = sprintf('%s\\g176_%02d.jpg', imfolder,ii);
end
b=zeros(1057,960,2);
im = cell(im_n,1);
for ii = 1:im_n
    im{ii} = imread(imfile{ii});
    im{ii}=cat(3,im{ii},b);
end

edge_list = zeros(im_n-1,2);
for ei = 1:im_n-1
    edge_list(ei,:) = [ei,ei+1];
end

imsize = zeros(im_n,3);

for ii = 1:im_n
    imsize(ii,:) = size(im{ii});
    if imsize(ii,1) > 720
        scale = 720/size(im{ii}, 1);
        imsize(ii,:) = size(im{ii});
    end
end
refi=0;
projection_type='equi';
save_results = true;
recomp_global_paras = false;
compute_global_results = true;
show_intermediate_results = false;
blend_output = true;
bgcolor = 1; % 0 for black, 1 for white

im_n = size(im, 1);
edge_n = size(edge_list, 1);

imsize = zeros(im_n,3);
for i = 1:im_n
    imsize(i,:) = size(im{i});
end

%% Parameters
lambda = 0.001 * imsize(1,1)*imsize(1,2); % weighting parameter to balance the fitting term and the smoothing term
intv_mesh = 10; % interval in pixels for the computing of deformation functions
K_smooth = 5; % the smooth transition width in the non-overlapping region is set to K_smooth times of the maximum bias.

%% feature detection and matching
[X,im_ch]=matching(im_n,im,edge_n,edge_list);

%% robust elastic warping 
[R,M,D,ur,vr,mosaich,mosaicw,m_v0_,m_v1_,m_u0_,m_u1_,fe,imh_,imw_,ubox,vbox,R_pair]=robust_elastic_warp(im_n,X,imfolder,imsize,edge_list);
%% mosaic
single_band_mosaic(im_n,X,imfolder,imsize,edge_n,edge_list,im_ch,im,lambda,K_smooth,intv_mesh,R,M,D,ur,vr,mosaich,mosaicw,m_v0_,m_v1_,m_u0_,m_u1_,fe,imh_,imw_,ubox,vbox,R_pair);