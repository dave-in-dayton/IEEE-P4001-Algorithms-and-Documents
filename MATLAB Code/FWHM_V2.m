function [xpeak, xmid, fwhm, ypeak] = FWHM_V2(x, y)
%FWHM Improved center and width metrics for HS response functions.
%
%   Inputs x and y are assumed to arise from lab tests of the SPSF, LRF, or
%   SRF. The x-coordinate is either an input angle or wavelength, and y is
%   the system response in ADU's. Under the assumption of a well-planned
%   and executed data collection session, the y-values contain one or more
%   maxima relative to a noise baseline. 
%
%   Given a half-max threshold, we require at least one positive-going
%   threshold crossing to the left of the first maxima and one
%   negative-going threshold crossing to the right of the last maxima.
%   These point pairs determine half-response points taken as nearest
%   neigbor or interpolated values. 
%
%   This setup implies a minimum sequence length of 3 when there is
%   a unique maximum, and allows determination of widths as narrow as
%   one point-spacing interval.
%
%   xpeak = x-value of peak, or average of first and last peak location
%           x-values
%   ypeak = maximum y value
%   xmid  = average of first and last 50% response x coordinates
%   fwhm  = difference of first and last 50% response x coordinates
%
% D. Perry, Leidos, Dayton, OH 
% February  2023, created from FWHM.m
% September 2023, added ypeak as fourth return value

% first check inputs,
if isempty(x) || isempty(y) || length(x) ~= length(y)
    xpeak = NaN; xmid = NaN; fwhm = NaN; ypeak = NaN;
    return
end

% then locate xpeak (as location of max positive value) and calculate the
% 50% response value
[ypeak, ind1] = max(y);

[~, ind2] = max(flip(y));

ind2 = length(y)- ind2 + 1;

% detect/resolve cases where two identical maxima are present
if ind1 == ind2
    xpeak = x(ind1);
else
    xpeak = (x(ind1) + x(ind2)) / 2.0;
end

halfmax = ypeak / 2.0;

flag1 = false; flag2 = false;

flag3 = false; flag4 = false;

L = length(y);

if ind1 > 1 && ind2 < L
    
    % process lower data segment
    abovethreshind = y(2:ind1)   >   halfmax;

    atthreshind    = y(2:ind1)   ==  halfmax;

    belowthreshind = y(1:ind1-1) <   halfmax;

    % check for positive-going crossing
    pos_crossing = belowthreshind & abovethreshind;
 
    first_crossing = find(pos_crossing,1,'first');

    if ~isempty(first_crossing)
            lowxind = first_crossing;
            flag1 = true;      
    end

    % check for near-crossing terminating on threshold
    at_thresh    = belowthreshind & atthreshind;

    first_at_thresh = find(at_thresh,1,'first');

    if ~isempty(first_at_thresh)
            % note actual position of half max
            lowxind = first_at_thresh + 1;
            flag2 = true;      
    end

    
    % repeat for upper data segment
    abovethreshind = y(ind2:end-1) >  halfmax;

    atthreshind    = y(ind2+1:end) == halfmax;

    belowthreshind = y(ind2+1:end) < halfmax;

    % check for negative crossing
    neg_crossing = abovethreshind & belowthreshind; 

    last = find(neg_crossing,1,'last');

    if ~isempty(last)
            highxind = ind2 + last - 1;
            flag3 = true;      
    end

    % check for near-crossing terminating on threshold
    at_thresh    = abovethreshind & atthreshind;

    last_at_thresh = find(at_thresh,1,'first');

    if ~isempty(last_at_thresh)
            highxind = ind2 + last_at_thresh;
            flag4 = true;      
    end
    
end

% compute xmid and fwhm if the lower and upper half response points or
% crossings were found; otherwise set results to NaN
if (flag1 || flag2) && (flag3 || flag4)

    % process lower data segment
    if flag1
%         % select x,y point closest to half max
%         err1 = abs(y(lowxind)   - halfmax);
%         err2 = abs(y(lowxind+1) - halfmax);
%         
%         if err1 == err2
%             lowx = (x(lowxind) + x(lowxind+1))/2.0;
%         elseif err1 < err2
%             lowx  = x(lowxind);
%         else
%             lowx  = x(lowxind+1);
%         end
        
        % 3-9-23: new way, not as desirable due to interpolation, but much
        % better performance!
        lowx = interp1([y(lowxind+1),y(lowxind)],[x(lowxind+1),x(lowxind)],halfmax);
        
    else
        % flag2; found actual half max
        lowx  = x(lowxind);
    end

    % process upper data segment
    if flag3
%         % select x/y point closest to half max
%         err1 = abs(y(highxind)   - halfmax);
%         err2 = abs(y(highxind+1) - halfmax);
%         
%         if err1 == err2
%             highx  = (x(highxind) + x(highxind+1))/2.0;
%         elseif err1 < err2
%             highx  = x(highxind);
%         else
%             highx  = x(highxind+1);
%         end
        
        % 3-9-23: new way
        highx = interp1([y(highxind+1),y(highxind)],[x(highxind+1),x(highxind)],halfmax);  
        
    else
        % flag 4; found actual half max at
        highx = x(highxind);
    end

    fwhm  =  highx - lowx;
    xmid  = (highx + lowx)/2.0;

else
    xmid  = NaN;
    fwhm  = NaN;
end

end

