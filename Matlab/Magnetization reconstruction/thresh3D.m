function thres = thresh3D(x,SORH,T)

%This function return the 3-dimensional soft thresholding of a family of vectors on R^3
%
%The "x" must be given as a vertical vector of the form [x1;x2;x3] where
%[x1(i), x2(i), x3(i)] ith vector of the family.
%The result "thresh" is given as a vector of the same form as x.
%
%SORH should be either 's' or 'h' as a character for soft of hard
%thresholding.
%
% Changelog 20231023
% 1.- Added an exeption in case the funtion wthresh is not accessible
% 2.- Commented out unnecessay if's

% if isequal(class(x),'mp')
%     rea  = reshape(x,[],3);
%     nrea = repmat(sqrt(sum(rea.^2,2)),1,3);
%     dir  = rea./nrea;
%     dir(isnan(dir) | dir==Inf | dir==-Inf) = mp(0,5);
%     norm_thres = max(0, nrea-T);
%     thres = reshape(norm_thres.*dir,[],1);
%     
% elseif isa(x,'gpuArray')
%         
%     n   = size(x,1);
%     nx  = repmat(sqrt(x(1:n/3).^2+x(n/3+1:2*n/3).^2+x(2*n/3+1:n).^2),3,1);
%     dir = x./nx;
%     nt  = wthresh(nx,SORH,T);
% 
%     thres = nt.*dir;
%     
% else
    
    rea  = reshape(x,[],3);
    nrea = repmat(sqrt(sum(rea.^2,2)),1,3);
    dir  = rea./nrea;
    dir(isnan(dir) | dir==Inf | dir==-Inf) = 0;
    
    try
      norm_thres = wthresh(nrea,SORH,T);
    catch
      norm_thres = (nrea-T).*(nrea>=T);
    end

    thres = reshape(norm_thres.*dir,[],1);
    
%end

end