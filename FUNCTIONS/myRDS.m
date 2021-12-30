function [L1,R1] = myRDS(dc,dr,flag,myseed,varargin)
%Generate a Random Dot Stereogram (RDS) input. RDS is a stereo image where
%each pixel has (randomly) -1 or 1 as value. If the RDS is correlated
%then the two images (right and left) have the same seed but the right image is
%shifted of a certain value respect to the left image.
%   dc:         displacement along columns
%   dr:         displacement along rows
%   flag:       '1' -> RDS correlated,
%               '2' -> RDS correlated (background + foreground),
%               '3' -> RDS uncorrelated
%   myseed      seed of random number generator
%   varargin{1} number of columns of frames (256 by default)
%   varargin{2} number of rows of frames (256 by default)
%
%Fulvio Missoni
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CONTROLS
if nargin>4
    sx = varargin{1};
    sy = varargin{2};
else
    sx = 256;
    sy = 256;
end
    
if ~isint(sx)||~isint(sy)
    error('The size must be an integer value!')
end

%switch-case
if flag==1
    %RDS correlated 
    if ~isempty(myseed)
        rng(myseed)
    end 
%     L1 = randi([0 1],[sy*2 sx*2])*2-1;
    L1 = rand(sy*2, sx*2)*2-1;
    R1 = L1;
    R1 = circshift(R1,dr,1);
    R1 = circshift(R1,dc,2);
    L1 = L1(floor(sy/2):floor(3/2*sy),floor(sx/2):floor(3/2*sx));
    R1 = R1(floor(sy/2):floor(3/2*sy),floor(sx/2):floor(3/2*sx));
end
if flag == 2
    %Background 'static' and 'foreground' shifted
    %RDS correlated with central patch shift the central patch and the 
    %external patch are uncorrelated
    if ~isempty(myseed)
        rng(myseed)
    end 
%     L1 = randi([0 1],[sy sx])*2-1;
    L1 = rand(sy, sx)*2-1;

    R1 = L1;
    if ~isempty(myseed)
        myseed2 = myseed/2;
        rng(myseed2)    %use another seed
    end    
    cl = rand(sy/2+1,sx/2+1)*2-1;
%     cl = randi([0 1],[sy/2+1 sx/2+1])*2-1;
    cr = cl;
    
    st_indx = floor(sx/4);      st_indy = floor(sy/4);      
    end_indx = floor(3*sx/4);   end_indy = floor(3*sy/4);
    L1(st_indx:end_indx,st_indy:end_indy)=cl;
    R1(st_indx+dr:end_indx+dr,st_indy+dc:end_indy+dc)=cr;
end
if flag == 3
    %RDS uncorrelated
    if ~isempty(myseed)
        rng(myseed)
    end 
% 
%     L1 = randi([0 1],[sy*2 sx*2])*2-1;
%     R1 = randi([0 1],[sy*2 sx*2])*2-1;
    L1 = rand(sy*2,sx*2)*2-1;
    if ~isempty(myseed)
        myseed2 = myseed/2;
        rng(myseed2)    %use another seed
    end
    R1 = rand(sy*2,sx*2)*2-1;
    R1 = circshift(R1,dr,1);
    R1 = circshift(R1,dc,2);
    L1 = L1(floor(sy/2):floor(3/2*sy),floor(sx/2):floor(3/2*sx));
    R1 = R1(floor(sy/2):floor(3/2*sy),floor(sx/2):floor(3/2*sx));
end
end

function [bool,idx] = isint(x)
    % From Matlab forum
    % Check whether input is integer or not
    % Inf and NaN are not integers
    if ~isnumeric(x)
       error('Input must be a numeric, not a %s.',class(x))
    end
    bool = (mod(x,1) == 0);
%     bool = (round(x) == floor(x));  % Other approach. Fails with Inf and
%     NaN.
    idx = find(bool); 
% Manolín Sept-2019
end