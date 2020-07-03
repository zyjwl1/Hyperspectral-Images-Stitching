function [X,im_ch]=matching(im_n,im,edge_n,edge_list)
im_gray = cell(im_n, 1);
points = cell(im_n, 1);
features = cell(im_n, 1);
valid_points = cell(im_n, 1);

% detection
for i = 1:im_n
    im_ch = size(im{i},3);
    if im_ch > 1
        im_gray{i} = im2single(rgb2gray(im{i}));
    elseif im_ch == 1
        im_gray{i} = im2single(im{i});
    end
    if exist('vl_sift', 'file')
        [ points{i},features{i} ] = vl_sift(single(im_gray{i}),'PeakThresh', 0,'edgethresh',500);
    else
        points{i} = detectSURFFeatures(im_gray{i},...
            'MetricThreshold', 0, 'NumOctaves', 3, 'NumScaleLevels', 4);
        [features{i}, valid_points{i}] = extractFeatures(im_gray{i}, points{i});
    end
end

% matching
if edge_n == 0
    % image match verfication referring to AutoStitch
    X = cell(im_n*(im_n-1)/2, 2);
    for i = 1:im_n-1
        for j = i+1:im_n
            if exist('vl_sift', 'file')
                matches = vl_ubcmatch(features{i}, features{j});
                X_1 = [ points{i}(1:2,matches(1,:)) ; ones(1,size(matches,2)) ];
                X_2 = [ points{j}(1:2,matches(2,:)) ; ones(1,size(matches,2)) ];
            else
                indexPairs = matchFeatures(features{i}, features{j});
                matched_points_1 = valid_points{i}(indexPairs(:, 1), :);
                matched_points_2 = valid_points{j}(indexPairs(:, 2), :);
                X_1 = [matched_points_1.Location'; ones(1,size(matched_points_1, 1))];
                X_2 = [matched_points_2.Location'; ones(1,size(matched_points_2, 1))];
            end
                X{edge_n,1} = X_1(:,:);
                X{edge_n,2} = X_2(:,:);
        end
    end
    X = X(1:edge_n,:);
else
    X = cell(edge_n, 2);
    for ei = 1 : edge_n
        i = edge_list(ei, 1);
        j = edge_list(ei, 2);
        
        if exist('vl_sift', 'file')
            matches = vl_ubcmatch(features{i}, features{j});
            X_1 = [ points{i}(1:2,matches(1,:)) ; ones(1,size(matches,2)) ];
            X_2 = [ points{j}(1:2,matches(2,:)) ; ones(1,size(matches,2)) ];
        else
            indexPairs = matchFeatures(features{i}, features{j});
            matched_points_1 = valid_points{i}(indexPairs(:, 1), :);
            matched_points_2 = valid_points{j}(indexPairs(:, 2), :);
            X_1 = [matched_points_1.Location'; ones(1,size(matched_points_1, 1))];
            X_2 = [matched_points_2.Location'; ones(1,size(matched_points_2, 1))];
        end
       [X{ei,1},X{ei,2}]= mTopKRP(X_1(:,:),X_2(:,:));
    end
end
end