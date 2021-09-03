%% imageCorrelation
% debug the cleaning of C (varargin in function)
function [xoffSet,yoffSet] = imageCorrelation(xc,yc,ACC1,ACC2,w,filterOrder,varargin)
% varargin: ,'cleanC',dxPass01,dyPass01,R);
ACC1sub = zeros(w+1,w+1,'uint8');
ACC1sub = ACC1(yc-w:yc+w,xc-w:xc+w);
C = normxcorr2(imgaussfilt(ACC1sub,filterOrder),imgaussfilt(ACC2,filterOrder));

% set C to zero above a predefined radius
% Checking varargin structure
if ( length(varargin) > 1 )
    %fprintf('cleaning C \n')
    dxPass01 = double(varargin{:,2});
    dyPass01 = double(varargin{:,3});
    R = double(varargin{:,4});
    x0 = round(xc+dxPass01 + size(ACC1sub,1)/2);
    y0 = round(yc+dyPass01 + size(ACC1sub,2)/2);
    x = 1:size(C,2);
    y = 1:size(C,1);
    [xx,yy] = meshgrid(x,y);
    %     figure
    %     imagesc(C)
    C(((xx-x0).^2+(yy-y0).^2) > R^2)=0;
    %     figure
    %     imagesc(C)
end

[ypeak,xpeak] = find(C==max(C(:)));
yoffSet = ypeak-size(ACC1sub,1) + w;
xoffSet = xpeak-size(ACC1sub,2) + w;
end
