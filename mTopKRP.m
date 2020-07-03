%%  This is a demo for removing outliers Using mTopKRP.
%Edit by: Yujie Zhang 
%Date:    03/5/2020
%  Implements the Multiscale Top K Rank Preserving Matching [1].
%  References:
%   [1]Y. Zhang, Z. Wan, X. Jiang and X. Mei, "Automatic Stitching for Hyperspectral Images Using Robust Feature Matching and Elastic Warp," in IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, vol. 13, pp. 3145-3154, 2020, doi: 10.1109/JSTARS.2020.3001022.
%%
function [X,Y]=mTopKRP(X,Y)
    [numx1,numx2] = size(X);
    Num=numx1;
    if numx1<numx2
       Xt=X;Yt=Y;
       Num=numx2;
    else
        Xt = X';Yt = Y';
    end
 tic   
% 
idxUnique=1:Num;

% Parameters Setting
tic
    tau1 = 0.80;     numNeigh1 = 15;
    tau2 = 0.35;     numNeigh2 = 15;
    tau3 = 0.35;     numNeigh3 = 15;
    tau4 = 0.35;     numNeigh4 = 10; 

    conf=[];
    conf.iteration = 3; 
    conf.scale = 3;  
    conf.interval = 2; 
    conf.Num = Num;
    
% iteration1
    Idx = GetCorectIdx(Xt,Yt,numNeigh1,tau1, conf,idxUnique);
    
% iteration2
  if  conf.iteration>=2
      Idx = GetCorectIdx(Xt,Yt,numNeigh2,tau2, conf,Idx);
  end
  
% iteration3
  if  conf.iteration>=3
    Idx = GetCorectIdx(Xt,Yt,numNeigh3,tau3, conf,Idx);
  end

n = size(X,1);
tmp=zeros(1, n);
tmp(Idx) = 1;
tmp(Idx) = tmp(Idx)+1;
Correct = find(tmp == 2);
TruePos = Correct;   %Ture positive
tmp=zeros(1, n);
tmp(Idx) = 1;
tmp(Idx) = tmp(Idx)-1;
FalsePos = find(tmp == 1); %False positive
tmp=zeros(1, n);
tmp(Idx) = 1;
tmp(Idx) = tmp(Idx)-1;
FalseNeg = find(tmp == 1); %False negative

FP = FalsePos;
FN = FalseNeg;

NumPos = length(TruePos)+length(FalsePos)+length(FalseNeg);
    n1 = length(TruePos);
    n2 = length(FalsePos);
    n3 = length(FalseNeg);

per = randperm(length(TruePos));
TruePos = TruePos(per(1:n1));
per = randperm(length(FalsePos));
FalsePos = FalsePos(per(1:n2));
per = randperm(length(FalseNeg));
FalseNeg = FalseNeg(per(1:n3));
X=X(:,TruePos);
Y=Y(:,TruePos);