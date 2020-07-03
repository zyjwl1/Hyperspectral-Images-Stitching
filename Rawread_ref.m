function [ X ] = Rawread_ref( file, M, N, B)
M=960;N=1057;B=176;
fid=fopen(file,'r');
image=fread(fid, 'float');
image=reshape(image,M,B,N);
for i=1:B
    X(:,:,i) = reshape(image(:,i,:), M,N)';
end


