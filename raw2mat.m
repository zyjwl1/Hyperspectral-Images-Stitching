clc;
clear;
close all;
warning off;
addpath(genpath(cd));
X1 = Rawread_ref( '.raw', 960,1057,176 );
figure;
imshow(X1(:,:,100));
imwrite(X1(:,:,100),'gray05956.jpg');

% img_n=2;
%  imgList = dir('../*.mat');
%  in_name = cell(img_n,1);
%  for i = 1 : length(in_name)
%    in_name{1} = ['1.mat'];
%       in_name{2} = ['2.mat'];
%  end
% III = cell(img_n, 1);
% for i = 1 : img_n
%     III{i} = cell2mat(struct2cell(load(in_name{i})));
% end
% I2=zeros(1057,960,176);
% for xx=1:2
% I2=III{xx};
% for yy=13
%     I22=I2(:,:,yy);
%     if yy<=9
%     imwrite(I22, sprintf('g0%d_0%d.jpg', yy, xx));
%     else
%         imwrite(I22, sprintf('g%d_%d.jpg', yy, xx));
%     end
%             
% end
% end
