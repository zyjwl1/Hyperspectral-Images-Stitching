function [ X ] = Rawread( file, M, N, B)
fid=fopen(file,'r');
image=fread(fid, 'int16');
image=reshape(image,M,B,N);
for i=1:B
    X(:,:,i) = reshape(image(:,i,:), M,N)';
end
X = X/2^16;
end

