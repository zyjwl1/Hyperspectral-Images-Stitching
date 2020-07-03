function single_band_mosaic(im_n,X,imfolder,imsize,edge_n,edge_list,im_ch,im,lambda,K_smooth,intv_mesh,R,M,D,ur,vr,mosaich,mosaicw,m_v0_,m_v1_,m_u0_,m_u1_,fe,imh_,imw_,ubox,vbox,R_pair)
%% choose one band as the reference band
refi=0;
projection_type='equi';
save_results = true;
show_intermediate_results = false;
bgcolor = 1; % 0 for black, 1 for white%% mosaic

% elastic local alignment and mosaiking
Adj = zeros(im_n,im_n); % adjacent matrix describing the topological 
                        % relationship of the input images,
                        % the elements of Adj indicate the edge indexes
                        % between the two overlapping images, the 0
                        % elements indicate not overlapped
for ei = 1:edge_n
    i = edge_list(ei, 1);
    j = edge_list(ei, 2);
    Adj(i,j) = ei;
    Adj(j,i) = ei;
end

[u,v] = meshgrid(ur,vr) ;

imi_ = cell(im_n,1); % only for intermediate results
im_p = cell(im_n,1);
mask = cell(im_n,1);
mass = zeros(mosaich, mosaicw, im_ch);
mosaic = zeros(mosaich, mosaicw, im_ch);
for ki = 1:im_n
    i = mod(ki + refi - 2, im_n) + 1; % start from the reference image
    % align image i
    u_im = u(m_v0_(i):m_v1_(i),m_u0_(i):m_u1_(i));
    v_im = v(m_v0_(i):m_v1_(i),m_u0_(i):m_u1_(i));
    if strcmp(projection_type,'equi')
        [u_im_, v_im_] = trans_equi2persp(u_im, v_im, R{i}, M{i}, D{i}, fe);
    else
        [u_im_, v_im_] = trans_persp2persp(u_im, v_im, R{i}, Mp, Dp, M{i}, D{i});
    end
    
    need_deform = false;
    sub_u0_ = [];
    sub_u1_ = [];
    sub_v0_ = [];
    sub_v1_ = [];
    Pi = [];
    Pi_ = [];
    for kj = 1:ki-1 % for every image that has already been aligned
        j = mod(kj + refi - 2, im_n) + 1;
        if Adj(i,j) > 0
            need_deform = true;
            
            [ubox_ji, vbox_ji] =  trans_persp2persp(ubox{j}, vbox{j}, R_pair{j,i}, M{j}, D{j}, M{i}, D{i});
            sub_u0_ = cat(1,sub_u0_,max([1, min(ubox_ji)]) );
            sub_u1_ = cat(1,sub_u1_,min([imsize(i,2), max(ubox_ji)]) );
            sub_v0_ = cat(1,sub_v0_,max([1, min(vbox_ji)]) );
            sub_v1_ = cat(1,sub_v1_,min([imsize(i,1), max(vbox_ji)]) );
            
            ei = Adj(i,j);
            if i == edge_list(ei, 1) && j == edge_list(ei, 2)
                Xi = X{ei,1};
                Xj = X{ei,2};
            else
                Xi = X{ei,2};
                Xj = X{ei,1};
            end
            [xj_i, yj_i] = trans_persp2persp(Xj(1,:), Xj(2,:), R_pair{j,i}, M{j}, D{j}, M{i}, D{i});
            Pi = cat(2, Pi, Xi(1:2,:));
            Pi_ = cat(2, Pi_, [xj_i;yj_i]);
        end
    end
    if need_deform
        sub_u0_ = min(sub_u0_);
        sub_u1_ = max(sub_u1_);
        sub_v0_ = min(sub_v0_);
        sub_v1_ = max(sub_v1_);
        
        % merge the coincided points
        ok_Pi = false(size(Pi,2),1);
        [~, idx_Pi] = unique(round(Pi'), 'rows', 'stable');
        ok_Pi(idx_Pi) = true;
        ok_Pi_ = false(size(Pi_,2),1);
        [~, idx_Pi_] = unique(round(Pi_'), 'rows', 'stable');
        ok_Pi_(idx_Pi_) = true;
        ok_nd = ok_Pi & ok_Pi_;
        Pi_nd = Pi(:,ok_nd);
        Pi_nd_ = Pi_(:,ok_nd);

        % form the linear system
        xi = Pi_nd(1,:);
        yi = Pi_nd(2,:);
        xj_ = Pi_nd_(1,:);
        yj_ = Pi_nd_(2,:);
        gxn = xj_ - xi;
        hyn = yj_ - yi;
        
        n = size(xj_, 2);
        xx = xj_(ones(1,n),:);
        yy = yj_(ones(1,n),:);
        dist2 = (xx - xx').^2 + (yy - yy').^2;
        dist2(1:n+1:n*n) = ones(1,n);
        K = 0.5 * dist2 .* log(dist2);
        K(1:n+1:n*n) = lambda * 8*pi * ones(1,n);
        K_ = zeros(n+3,n+3);
        K_(1:n,1:n) = K;
        K_(n+1,1:n) = xj_;
        K_(n+2,1:n) = yj_;
        K_(n+3,1:n) = ones(1,n);
        K_(1:n,n+1) = xj_';
        K_(1:n,n+2) = yj_';
        K_(1:n,n+3) = ones(n,1);
        G_ = zeros(n+3,2);
        G_(1:n,1) = gxn';
        G_(1:n,2) = hyn';
        
        % solve the linear system
        W_ = K_\G_;
        wx = W_(1:n,1);
        wy = W_(1:n,2);
        a = W_(n+1:n+3,1);
        b = W_(n+1:n+3,2);
        
        % remove outliers based on the distribution of weights
        outlier = abs(wx)>3*std(wx) | abs(wy)>3*std(wy);
        
        inlier_idx = 1:size(xi, 2);
        for kiter = 1:10
%             if ~any(outlier)
            if sum(outlier) < 0.0027*n
                break;
            end
            
            ok = ~outlier;
            inlier_idx = inlier_idx(ok);
            K_ = K_([ok;true(3,1)],[ok;true(3,1)]);
            G_ = G_([ok;true(3,1)],:);
            W_ = K_\G_;
            n = size(inlier_idx,2);
            wx = W_(1:n,1);
            wy = W_(1:n,2);
            a = W_(n+1:n+3,1);
            b = W_(n+1:n+3,2);
            outlier = abs(wx)>3*std(wx) | abs(wy)>3*std(wy);
        end
        ok = false(size(xj_, 2),1);
        ok(inlier_idx) = true;
        xj_ = xj_(ok);
        yj_ = yj_(ok);
        gxn = gxn(ok);
        hyn = hyn(ok);
        
        eta_d0 = 0; % lower boundary for smooth transition area
        eta_d1 = K_smooth * max(abs([gxn, hyn])); % higher boundary for smooth transition area
%         eta_d1 = 0.3*max(imsize(1,1:2));
        sub_u0_ = sub_u0_ + min(gxn);
        sub_u1_ = sub_u1_ + max(gxn);
        sub_v0_ = sub_v0_ + min(hyn);
        sub_v1_ = sub_v1_ + max(hyn);
        
        if show_intermediate_results
%             margin = df_max;
            [u_im,v_im] = meshgrid(1:intv_mesh:imsize(i,2),1:intv_mesh:imsize(i,1));
            gx_im = zeros(ceil(imsize(i,1)/intv_mesh),ceil(imsize(i,2)/intv_mesh));
            hy_im = zeros(ceil(imsize(i,1)/intv_mesh),ceil(imsize(i,2)/intv_mesh));
            for kf = 1:n
                dist2 = (u_im - xj_(kf)).^2 + (v_im - yj_(kf)).^2;
                rbf = 0.5 * dist2 .* log(dist2);
                gx_im = gx_im + wx(kf)*rbf;
                hy_im = hy_im + wy(kf)*rbf;
            end
            gx_im = gx_im + a(1).*u_im+a(2).*v_im+a(3);
            hy_im = hy_im + b(1).*u_im+b(2).*v_im+b(3);
            gx_im = imresize(gx_im, [imsize(i,1),imsize(i,2)]);
            hy_im = imresize(hy_im, [imsize(i,1),imsize(i,2)]);
            
            [u_im,v_im] = meshgrid(1:imsize(i,2),1:imsize(i,1)) ;
            dist_horizontal_im = max(sub_u0_-u_im, u_im-sub_u1_);
            dist_vertical_im = max(sub_v0_-v_im, v_im-sub_v1_);
            dist_sub_im = max(dist_horizontal_im, dist_vertical_im);
            dist_sub_im = max(0, dist_sub_im);
            eta_im = (eta_d1 - dist_sub_im) ./ (eta_d1 - eta_d0);
            eta_im(dist_sub_im < eta_d0) = 1;
            eta_im(dist_sub_im > eta_d1) = 0;
            gx_im = gx_im .* eta_im;
            hy_im = hy_im .* eta_im;
            
            u_im = u_im - gx_im;
            v_im = v_im - hy_im;
            imi_{i} = zeros(imsize(i,:));
            for kc = 1:imsize(i,2)
                imi_{i}(:,:,kc) = interp2(im2double(im{i}(:,:,kc)),u_im,v_im);
            end
           
        end
        
        u_mesh_ = u_im_(1:intv_mesh:imh_(i),1:intv_mesh:imw_(i));
        v_mesh_ = v_im_(1:intv_mesh:imh_(i),1:intv_mesh:imw_(i));
        gx_mesh_ = zeros(ceil(imh_(i)/intv_mesh), ceil(imw_(i)/intv_mesh));
        hy_mesh_ = zeros(ceil(imh_(i)/intv_mesh), ceil(imw_(i)/intv_mesh));
        for kf = 1:n
            dist2 = (u_mesh_ - xj_(kf)).^2 + (v_mesh_ - yj_(kf)).^2;
            rbf = 0.5 * dist2 .* log(dist2);
            gx_mesh_ = gx_mesh_ + wx(kf)*rbf;
            hy_mesh_ = hy_mesh_ + wy(kf)*rbf;
        end
        gx_mesh_ = gx_mesh_ + a(1).*u_mesh_+a(2).*v_mesh_+a(3);
        hy_mesh_ = hy_mesh_ + b(1).*u_mesh_+b(2).*v_mesh_+b(3);
        gx_im_ = imresize(gx_mesh_, [imh_(i),imw_(i)]);
        hy_im_ = imresize(hy_mesh_, [imh_(i),imw_(i)]);
        
        %smooth tansition to global transform
        dist_horizontal = max(sub_u0_-u_im_, u_im_-sub_u1_);
        dist_vertical = max(sub_v0_-v_im_, v_im_-sub_v1_);
        dist_sub = max(dist_horizontal, dist_vertical);
        dist_sub = max(0, dist_sub);
        eta = (eta_d1 - dist_sub) ./ (eta_d1 - eta_d0);
        eta(dist_sub < eta_d0) = 1;
        eta(dist_sub > eta_d1) = 0;
        gx_im_ = gx_im_ .* eta;
        hy_im_ = hy_im_ .* eta;
        
        u_im_ = u_im_ - gx_im_;
        v_im_ = v_im_ - hy_im_;
        
        % update the feature locations
        for kj = ki+1:im_n % for every image that has not been aligned
            j = mod(kj + refi - 2, im_n) + 1;
            if Adj(i,j) > 0
                ei = Adj(i,j);
                if i == edge_list(ei, 1) && j == edge_list(ei, 2)
                    Xi = X{ei,1};
                else
                    Xi = X{ei,2};
                end
                newXi = Xi;
                % Iterative solver
                for kiter = 1:20
                    u_f = newXi(1,:);
                    v_f = newXi(2,:);
                    gx_f = zeros(1,size(newXi,2));
                    hy_f = zeros(1,size(newXi,2));
                    for kf = 1:n
                        dist2 = (u_f - xj_(kf)).^2 + (v_f - yj_(kf)).^2;
                        rbf = 0.5 * dist2 .* log(dist2);
                        gx_f = gx_f + wx(kf)*rbf;
                        hy_f = hy_f + wy(kf)*rbf;
                    end
                    gx_f = gx_f + a(1).*u_f+a(2).*v_f+a(3);
                    hy_f = hy_f + b(1).*u_f+b(2).*v_f+b(3);
                    dist_horizontal_f = max(sub_u0_-u_f, u_f-sub_u1_);
                    dist_vertical_f = max(sub_v0_-v_f, v_f-sub_v1_);
                    dist_sub_f = max(dist_horizontal_f, dist_vertical_f);
                    dist_sub_f = max(0, dist_sub_f);
                    eta_f = (eta_d1 - dist_sub_f) ./ (eta_d1 - eta_d0);
                    eta_f(dist_sub_f < eta_d0) = 1;
                    eta_f(dist_sub_f > eta_d1) = 0;
                    gx_f = gx_f .* eta_f;
                    hy_f = hy_f .* eta_f;
%                     disp([sum(abs(Xi(1,:) + gx_f - newXi(1,:))),sum(abs(Xi(2,:) + hy_f - newXi(2,:)))]);
                    newXi(1,:) = Xi(1,:) + gx_f;
                    newXi(2,:) = Xi(2,:) + hy_f;
                end
                if i == edge_list(ei, 1) && j == edge_list(ei, 2)
                    X{ei,1} = newXi;
                else
                    X{ei,2} = newXi;
                end
            end
        end
    end
    im_p{i} = zeros(imh_(i),imw_(i),imsize(i,3));
    for kc = 1:imsize(i,3)
        im_p{i}(:,:,kc) = interp2(im2double(im{i}(:,:,kc)),u_im_,v_im_);
    end
    mask{i} = ~isnan(im_p{i});
    mass(m_v0_(i):m_v1_(i),m_u0_(i):m_u1_(i),:)...
        = mass(m_v0_(i):m_v1_(i),m_u0_(i):m_u1_(i),:) + mask{i};
    im_p{i}(isnan(im_p{i})) = 0;
    mosaic(m_v0_(i):m_v1_(i),m_u0_(i):m_u1_(i),:)...
        = mosaic(m_v0_(i):m_v1_(i),m_u0_(i):m_u1_(i),:) + im_p{i};
end
mosaic = mosaic ./ mass;
mosaic(isnan(mosaic)) = bgcolor;

figure;
imshow(mosaic(:,:,1), 'border', 'tight') ;
drawnow;
if save_results
    imwrite(mosaic, [imfolder, '\mosaic.jpg']);
end
end
