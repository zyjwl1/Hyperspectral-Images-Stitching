%Edit by: Yujie Zhang 
%Date:    03/5/2020
%  Implements the Robust Elastic Warps [2].
%  References:
%   [2]J. Li, Z. Wang, S. Lai, Y. Zhai, and M. Zhang, ¡°Parallax-tolerant image stitching based on robust elastic warping,¡± IEEE Trans. Multimed.,vol. 20, no. 7, pp. 1672¨C1687, 2017.
function [R,M,D,ur,vr,mosaich,mosaicw,m_v0_,m_v1_,m_u0_,m_u1_,fe,imh_,imw_,ubox,vbox,R_pair]=robust_elastic_warp(im_n,X,imfolder,imsize,edge_list)
refi=0;
projection_type='equi';
recomp_global_paras = false;
if exist([imfolder,'\global_paras.mat'],'file') && ~recomp_global_paras
    load([imfolder,'\global_paras.mat']);
else
    options = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt', 'Display','final',...
        'MaxFunEvals',1000*im_n, 'MaxIter',1e3, 'TolFun',1e-6, 'TolX',1e-6, 'Jacobian','off');
    paras_init = [1000*ones(1,im_n),zeros(1,3*(im_n-1))];
    for sigma = [1000, 100, 10]
        [paras, ~ ,~ ,exitflag] = lsqnonlin(...
            @(p)residual_all_robust(X, imsize, edge_list, p, sigma), paras_init,...
            [],[],options);
        if exitflag > 0
            paras_init = paras;
        end
    end
    save([imfolder,'\global_paras.mat'],'paras');
end

M = cell(im_n,1);
D = cell(im_n,1);
R_pair = cell(im_n, im_n); % pairwise rotation matrics
for i = 1 : im_n
    ki = paras(i);
    M{i} = [ki, 0, imsize(i,2)/2;
             0, ki, imsize(i,1)/2;
             0,  0, 1];
    D{i} = 0;
end
for i = 1 : im_n
    R_pair{i, i} = eye(3);
end
for i = 2:im_n
    theta = paras(im_n+3*(i-2)+1:im_n+3*(i-2)+3);
    theta_m = [0         -theta(3) theta(2)
               theta(3)  0         -theta(1)
               -theta(2) theta(1)  0];
    R_pair{1,i} = expm(theta_m);
    R_pair{i,1} = R_pair{1,i}';
end
for i = 2:im_n-1
    for j = i+1:im_n
        R_pair{i,j} = R_pair{1,j}*R_pair{i,1};
        R_pair{j,i} = R_pair{i,j}';
    end
end

if refi == 0
    % automatic straithtening
    xz_vecs = zeros(im_n*2,3);
    for i = 1:im_n
        xz_vecs(2*i-1:2*i,:) = R_pair{1,i}([1,3],:);
    end
    [~,~,V] = svd(xz_vecs,'econ');
    y_vec_glb = V(:,3);
    if y_vec_glb(2) < 0
        y_vec_glb = -y_vec_glb;
    end
    x_temp = zeros(3,1);
    for i = 1:im_n
        x_temp = x_temp + R_pair{1,i}(1,:)';
    end
    x_temp = x_temp ./ im_n;
    if norm(x_temp) < 0.1
        x_temp = [1;0;0];
    end
    z_vec_glb = cross(x_temp,y_vec_glb);
    z_vec_glb = z_vec_glb/norm(z_vec_glb);
    x_vec_glb = cross(y_vec_glb,z_vec_glb);
    % global rotation matrix for the reference image which is set to the first image here
    R_ref = [x_vec_glb,y_vec_glb,z_vec_glb]'; % global rotation matrix
    
    refi = 1;
else
    R_ref = eye(3);
end
R = cell(im_n, 1); % global rotation matrics
for i = 1 : im_n
    R{i} = R_pair{refi,i}*R_ref';
end


%% computing mosaic parameters
if strcmp(projection_type,'equi')
    fe = max(M{refi}(1,1),M{refi}(2,2));
else
    Mp = M{refi};
    Dp = D{refi};
end
ubox = cell(im_n,1);
vbox = cell(im_n,1);
ubox_ = cell(im_n,1);
vbox_ = cell(im_n,1);
ubox_all_ = [];
vbox_all_ = [];
for i = 1 : im_n
    ubox{i} = [1:imsize(i,2)        1:imsize(i,2)                     ones(1,imsize(i,1))  imsize(i,2)*ones(1,imsize(i,1))] ;
    vbox{i} = [ones(1,imsize(i,2))  imsize(i,1)*ones(1,imsize(i,2))  1:imsize(i,1)        1:imsize(i,1) ];
    if strcmp(projection_type,'equi')
        [ubox_{i}, vbox_{i}] =  trans_persp2equi(ubox{i}, vbox{i}, R{i}', M{i}, D{i}, fe);
    else
        [ubox_{i}, vbox_{i}] =  trans_persp2persp(ubox{i}, vbox{i}, R{i}', M{i}, D{i}, Mp, Dp);
    end
    ubox_all_ = cat(2,ubox_all_,ubox_{i});
    vbox_all_ = cat(2,vbox_all_,vbox_{i});
end
u0 = min(ubox_all_);
u1 = max(ubox_all_);
ur = u0:u1;
v0 = min(vbox_all_);
v1 = max(vbox_all_);
vr = v0:v1;
mosaicw = size(ur, 2);
mosaich = size(vr, 2);

m_u0_ = zeros(im_n,1);
m_u1_ = zeros(im_n,1);
m_v0_ = zeros(im_n,1);
m_v1_ = zeros(im_n,1);
imw_ = zeros(im_n,1);
imh_ = zeros(im_n,1);
for i = 1 : im_n
    % align the sub coordinates with the mosaic coordinates
    margin = 0.2 * min(imsize(1,1),imsize(1,2)); % additional margin of the reprojected image region considering the possilbe deformation
    u0_im_ = max(min(ubox_{i}) - margin, u0);
    u1_im_ = min(max(ubox_{i}) + margin, u1);
    v0_im_ = max(min(vbox_{i}) - margin, v0);
    v1_im_ = min(max(vbox_{i}) + margin, v1);
    m_u0_(i) = ceil(u0_im_ - u0 + 1);
    m_u1_(i) = floor(u1_im_ - u0 + 1);
    m_v0_(i) = ceil(v0_im_ - v0 + 1);
    m_v1_(i) = floor(v1_im_ - v0 + 1);
    imw_(i) = floor(m_u1_(i) - m_u0_(i) + 1);
    imh_(i) = floor(m_v1_(i) - m_v0_(i) + 1);
end