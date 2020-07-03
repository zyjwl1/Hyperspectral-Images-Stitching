function Idx = GetCorectIdx(Xt,Yt,numNeigh,tau, conf,idx)

scale = conf.scale;
interval = conf.interval;
KM = numNeigh+1; 
Num = conf.Num;
     for j=1:scale
     Km = KM-interval*(j-1);
     denom = sum(1./(1:Km));
     z = (Km - 4*floor(Km/2) + 2*(Km+1)*sum(1./(1:floor(Km/2))))/denom;
     PHI_MAX = (-2*Km + 2*z*denom);
     phi_o = PHI_MAX/(2*Km);
     conf.temp(j,:) = [Km,PHI_MAX,phi_o];
     end
     
    kdtreeX = vl_kdtreebuild(Xt(:,idx));
    kdtreeY = vl_kdtreebuild(Yt(:,idx));
    [neighborX, distX] = vl_kdtreequery(kdtreeX, Xt(:,idx), Xt, 'NumNeighbors', KM) ;
    [neighborY, distY] = vl_kdtreequery(kdtreeY, Yt(:,idx), Yt, 'NumNeighbors', KM) ;

    %% 
    D = zeros(1,Num);
    for id = 1:Num
            Di = MultiScale_Distance(neighborX(:,id),distX(:,id),neighborY(:,id), distY(:,id),conf.temp);
            D(id) = Di;
    end
    Idx = find(D <= tau);    
    
  end