function [solution,LE,DifRatios,RelResidual,Functional,steps,stop, ...
            extra_K] ...
                                    = Fistal1l2(b,AstarA,...
                                              AstarB,parh,par0,x0,model,...
                                              lambda,num_iterations,...
                                              stepsave,toler,l)
%This function approximate the solution of the problem
%argmin_x{|Ax-b|_2^2+lambda|x|_{1,2}
%using the Fista algorithm with Lipschitz constant the maximal eigenvalue 
%of AstarA.
%
%This algorithm is an adaptation of the one found in:
%Beck, Amir, and Marc Teboulle. "A fast iterative shrinkage-thresholding 
%algorithm for linear inverse problems." SIAM journal on imaging 
%sciences 2.1 (2009): 183-202
%
%x0 is the initial value for the algorithm, it must be a column vector in 
%the form [m1;m2;m3] where mi is the ith component of the magnetization
%
%The solution is a vector of the same form as x0.
%
%
% The variable "stepsave" states how often should the reconstruction is
% saved. If not given is taken as a floor(num_iterations/100)
%
% Also, this algorithm introduces some changes: 
%   - The variable "stop" takes value "SOL" when the maximum relative
%   difference between reconstructions is lower than 1e-15
%
% This version adds an extra parameter "extraparam" that tries to increase
% the speed of the algorithm by changing how much y differs from x on each
% step. It increased the difference as much as it cans at the beginning and
% starts decresing it as soon as the functional starts increasing.
%
% The variable "extra_K" denotes the number of times the algorithm went
% back because the extra parameter was too big
%
% Note that "stopping when the functional is increasing and the extra 
% parameter is zero" is valid since when the parameter is zero we go back
% to the ISTA algorithm and in such a case the article mentioned above
% mentions that the ISTA algorithm always gives a nonincreasing sequence of
% functionals.
%
% The "stop" variable is a string unlike other algorithms where it is a
% character array
%
% CHANGELOG
%   20200808
%       The when the algorithm stopped due to machine precision it was
%       returning the last reconstruction with had a bigger Functional.
%       This issue has been fixed.
%   20200813
%   1.- Lower the number of maximal iterations in K by -1 to try to avoid
%       errors.
%   2.- Commented out the original stopping criteria for SOL and FUN and
%       replaced the one for SOL for when the solution actually stopped
%       changing, and the one for FUN for when the functional had remained
%       as zero for 10 iterations.
%   20201023
%   1.- Corrected the stopping criteria when the functional stopped
%       changing. Before it stopped when one or the functional repeated and
%       now it only stops when it reaches 10 functionals that have the save
%       value
%   2.- Added a tolerance for the stopping criteria for the ratios of 
%       solution and functional difference.
%   20210302
%   1.- Changed the use of build_Ax_frac to build_Ax which supports the use
%       of vector valued parh.h.
%   20210308
%   1.- Added the display of the relative residual in the status text.
%   20210623
%   1.- Fixed an error where the indices of the variable DifRatios where
%       taken too high
%   20210705
%   1.- When eigs doesn't converge the matrix AstarA is now smoothed
%   20210806
%   1.- Made the text more legible

broke = false;


if nargin==9
    stepsave = floor(num_iterations/100);
    l        = eigs(AstarA,1)*2;
    toler    = 0;
elseif nargin == 10
    l     = eigs(AstarA,1)*2;
    toler = 0;
elseif nargin == 11
    l     = eigs(AstarA,1)*2;
end

factor = 1;
while isnan(l)
    warning('off','all')
    nA = sum(sum(AstarA.^2));
    nA = nA / numel(AstarA);
    nE = numel(find(AstarA.^2 < 1e-10*nA*factor));
    AstarA( AstarA.^2 < 1e-10*nA*factor ) = 0;
    l  = eigs(AstarA,1)*2*factor;
    factor = 2*factor;
    warning('on','all')
end
 
 
y  = x0;
xp = x0;
yp = x0;
t0 = 1;                   %auxiliary term

extra_K = 0;         %for this algorithm only
LE          = nan(1,100);
DifRatios   = nan(1,100);
Functional  = [];
stop        ="It";
extraparam  = 1;         %for this algorithm only
grow        = 2;         %for this algorithm only
OK          = 1;         %for this algorithm only
Ax0         = build_Ax2(parh,par0,x0,model);
Functional0 = norm(Ax0-b)^2+lambda*norm21(x0);
maxreldif   = [];

%n = numel(x0);
    

K = 1;

sizetext=0;

while K<floor(num_iterations/stepsave)-1
    %tic
for k=1:stepsave
    argu     = y-2/l*(AstarA*y-AstarB);   %auxiliary term
    solution = thresh3D(argu,'s',lambda/l); 
    
    t = (1+sqrt(1+4*t0^2))/2;
    d = solution-x0;
    y = solution+((t0-1)/t)*d*extraparam;
    
    x0=solution;   %setting terms for next iteration
    t0=t;
end

    ATRes = (AstarB-AstarA*solution);    %changed by experiment
    
    
argu           = solution+ATRes;    %
thresh_argu    = thresh3D(argu,'s',lambda/2);
LE(K)          = norm21(solution-thresh_argu)/norm21(solution); 
DifRatios(K)   = norm21(d)/norm21(solution); 
Asolution      = build_Ax2(parh,par0,solution,model);
RelResidual(K) = norm(Asolution-b)/norm(b); %#ok<AGROW>
Functional(K)  = (RelResidual(K)*norm(b))^2 + ...
                 lambda*norm21(solution); %#ok<AGROW>
solnon0        = x0;    % note that x0 and solution are the same now
solnon0(x0==0) = 1;
maxreldif(K)   = max(abs(d)./abs(solnon0));   %#ok<AGROW>

if K==1
    Functional_dif = (Functional0-Functional(K))/Functional(K);
else
    Functional_dif(K) = (Functional(K-1)-Functional(K))/Functional(K);
end

%Time = toc; 

    steps = K*stepsave;
    
text = sprintf(['    FISTA step %i, extra K = %i, ',...
                'Functional_dif %.2e, Extra_parameter 1+(%.2e),\n',...
                '    RelResidual %.2e'],...
                            steps, extra_K, Functional_dif(K), ...
                            extraparam-1,RelResidual(end));
eraser = repmat('\b',sizetext);
fprintf(eraser)
fprintf(text)

sizetext = size(text);

                        
if norm21(solution)==0
    fprintf('\n    Lambda is too big\n')
    stop = "lambda_big";
    broke   = true;
    break
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% added for this particular algorithm
if Functional_dif(K)<0 || isnan(Functional_dif(K))
    solution    = xp;
    if  extraparam  == 0
        fprintf('\n    Cannot continue due to Machine precision\n')
        stop = "Precision";
        broke   = true;
        break
    end
    x0          = xp;
    y           = yp;
    grow        = grow*(1-OK/4);
                                    % arbitrarily chosen way of changing 
                                    % the variable "grow". Note that  we 
                                    % don't want to decrease "grow" if it
                                    % was just decreased
    extraparam  = extraparam/(1+grow/stepsave); 
    extra_K = extra_K+1;
    OK          = OK-1;
    if OK < -4
        extraparam = 0;
        grow       = 5;
    end
else
    if isequal(xp,solution)
        fprintf('\n    The solution stopped changing\n')
        stop = "SOL";
        return
    end
    if K>10 && prod(~(Functional(end-10:end-1)-Functional(end)))
        fprintf('\n    The Functional stopped changing\n')
        stop = "FUN";
        return
    end
    xp = solution;
    yp = y;
    K  = K+1;
    if (K>2) && (9*Functional_dif(end-1) < 10*Functional_dif(end))
        % if the difference between functionals is not decreasing by at
        % least 9/10 we want to speed up by increasing extraparam
        extraparam = extraparam*(1+grow/stepsave);
        if extraparam == 0
            extraparam = 0.5;
        end
    end
    OK = 1;
    
    if K>6   %we need 6 because we are using "while" instead of "for"
    if prod(DifRatios(K-6:K-1) <= toler)
        fprintf('\n    The solution stopped changing up to tolerance\n')
        stop = "Tsol";
        return
    end
    if  prod(Functional_dif(end-5:end) < toler)
        fprintf('\n    The Functional stopped changing up to tolerance\n')
        stop = "Tfun";
        return
    end
    end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%



end

k = 0;

if ~broke

    for k = 1:stepsave-1
        argu     = y-2/l*(AstarA*y-AstarB);   %auxiliary term
        solution = thresh3D(argu,'s',lambda/l); 

        t = (1+sqrt(1+4*t0^2))/2;
        d = solution-x0;
        y = solution+((t0-1)/t)*d*extraparam;

        x0=solution;   %setting terms for next iteration
        t0=t;
    end


    %final step:
        argu     = y-2/l*(AstarA*y-AstarB);   %auxiliary term
        solution = thresh3D(argu,'s',lambda/l);

        ATRes = (AstarB-AstarA*solution);   %changed by experiment
    argu        = solution+ATRes;    %
    thresh_argu = thresh3D(argu,'s',lambda/2);
    LE(K+1)  = norm21(solution-thresh_argu)/norm21(solution);
    DifRatios(K+1)= norm21(d)/norm21(solution);

    k = k+1;
    
    if Functional_dif(K)<0 || isnan(Functional_dif(K))
        solution = xp;
        extra_K  = k;
        k        = 0;
    end

    fprintf('\n    Ended by number or iterations\n')
    
end

steps = steps+k;


end
