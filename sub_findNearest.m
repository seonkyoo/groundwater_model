
function [ixC, iyC, ixFce, iyFce] = sub_findNearest(xp,yp,xFaces,yFaces)
% All inputs are supposed to be a column vector
% select the neareset cell

    [~,ixFce] = min(abs(bsxfun(@minus,xFaces,xp'))); ixFce = ixFce';
    [~,iyFce] = min(abs(bsxfun(@minus,yFaces,yp'))); iyFce = iyFce';
    
    ixC = ixFce;
    ixC(xFaces(ixFce)>xp) = ixFce(xFaces(ixFce)>xp) - 1;
    
    iyC = iyFce;
    iyC(yFaces(iyFce)>yp) = iyFce(yFaces(iyFce)>yp) - 1;

end
