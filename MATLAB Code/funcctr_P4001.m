function [center] = funcctr_P4001(x, y, chanwidth, method)
    %FUNCCTR Computes center of function y = f(x) by specified method
    %
    %   x,y    = input data
    %
    %   method = 1  abscissa of maximum ordinate
    %   method = 2  average x-value for first/last 50% response points
    %   method = 3  1D centroid, per Mouroulis, Easton, etc.
    %   method = 4  AFRL median (the abscissa that divides y = f(x) into equal areas)
    %   method = 5  abscissa of peak of f(x) * rect(x/chanwidth)
    %   method = 6  first moment of probability distribution, with opt. noise mitig.
    %
    % Uses: FWHM_V2_P4001.m
    %
    % D. Perry, Dayton, Ohio
    % V1.0 August 2024
    verbose = false;
    
    % optional plot for method 5 - can also be used for width result
    do_plot = false;
    
    % check for valid x and y inputs
    if isempty(x) || isempty(y) || length(x) ~= length(y)
        center = NaN;
        return
    end
    
    % make sure x and y are row vectors and transpose as needed
    if ~isrow(x); x = x'; end
    
    if ~isrow(y); y = y'; end
    
    
    switch method
        
        case 1
            
            % x-value for maximum ordinate (average of first and last such
            % values when equal maxima exist)
            [center, ~, ~, ~] = FWHM_V2_P4001(x,y);
            
            return
            
            
        case 2
            
            % average of x-values at first and last 50% response points, now
            % using linear interpolation to refine earlier nearest-neighbor
            % approach
            [~, center, ~, ~] = FWHM_V2_P4001(x,y);
            
            return
            
            
        case 3
            
            % centroid of response; returns NaN when less than two points
            % are present or when the integral of the data is zero; refer
            % to Easton Eqn. 13.10 for the definition of centroid for a
            % continuous 1D function (with no restriction on the range of
            % f(x)); notice that the implementation below treats f(x) as a
            % continuous function (hence integration rather than summation)
            % while obviating the need for equally spaced x values
            if length(x) > 1
                
                denom = trapz(x,y);
                
                % set result to NaN if denominator is zero (i.e., as could
                % occur when when noise is present)
                if denom ~= 0.0
                    center = trapz(x, y .* x) / denom;
                else
                    center = NaN;
                end
                
            else
                center = NaN;
            end
            
            return
            
            
        case 4
            
            % median; (an be related to the 50% probability abcissa
            % when considering the response as a probability distribution )
            
            % first get cumulative distribution using the trapezoid rule
            % where we need not worry about unequally spaced x values;
            cumdist = cumtrapz(x,y);
            
            cumdist = cumdist / cumdist(end);
            
            %         % find nearest-neighbor x value for 50% area criteria
            %         [~, approxctrind] = min(abs(cumdist - 0.5));
            %
            %         center = x(approxctrind);
            
            % 3-9-2023: new way now involving linear interpolation rather than
            % earlier nearest-neighbor approach above
            medind = find(cumdist == 0.5,1);
            
            if isempty(medind)
                lowind = find(cumdist < 0.5,1,'last');
                center = interp1([cumdist(lowind+1),cumdist(lowind)],[x(lowind+1),x(lowind)],0.5);
            else
                center = x(medind);
            end
            
            return
            
            
        case 5
            
            % peak position of convolution with a rect() function of nominal
            % channel width
            
            % first check for equal point spacing
            deltaxvec = x(2:end) - x(1:end-1);
            
            deltax = median(deltaxvec);
            
            if (max(deltaxvec) - min(deltaxvec) ) / deltax > 0.0001
                if verbose
                    fprintf('ERROR - x-data for box center method is NOT evenly spaced!\n');
                    fprintf('Returning NaN for function center.\n\n');
                end
                center = NaN;
                return
            end
            
            % set up box/rect function
            numchansamples = round(chanwidth/deltax) + 1;
            
            boxfct = ones(1,numchansamples);
            
            % perform the convolution and find the location of its maximum;
            % returns location of first maximum if not unique
            convresult = conv(y, boxfct, 'same');
            
            [maxval, maxind] = max(convresult);
            
            center = x(maxind);
            
            % optional plot
            if do_plot
                
                boxplotdata = zeros(size(y));
                
                plotoffset = round(0.1*length(x));
                
                boxplotdata(1,plotoffset:plotoffset+length(boxfct)-1) = boxfct;
                
                convresult = convresult/maxval;
                
                figure; plot(x,y,x,boxplotdata,x,convresult); grid on;
                title('Box Function Center/Width Algorithm Components');
                legend({'Input','Box Function','Conv. of Input and Box'});
            end
            
            return
            
            
        case 6
            
            % first moment of response when it is assumed to be a
            % probability distribution (whose integral via trapezoid rule
            % s/b 1.0); identical results relative to the centroid above
            % for unity area responses, except for cases where negative
            % points are mitigated
            
            % optionally relax these conditions via the next two statements
            % to allow normalization and/or mitigation of negative values,
            % similar to what is done with funcwid.m's standard deviation
            % method
            no_mit     = false;
            
            allow_norm = true;
            
            if (min(y) < 0.0 && no_mit) || length(x) < 2
                
                center = NaN;
                
            else
                
                % perform mitigation for negative points if any are
                % present; method 3, and then 4 have been shown to provide
                % the highest flexibility scores
                mitig_mthd = 3;
                
                if min(y) < 0.0
                    
                    switch mitig_mthd
                        
                        case 1
                            % subset the data by eliminating points with negative
                            % values
                            posvalind = y >= 0.0;
                            
                            x = x(posvalind);
                            y = y(posvalind) ;
                            
                        case 2
                            % subset the data by thresholding the baseline region
                            % via the peak observed excursion below zero
                            yabovethreshind = y > abs(min(y));
                            
                            x = x(yabovethreshind);
                            y = y(yabovethreshind);
                            
                        case 3
                            % similar to case 1, but replace negative data with
                            % zeros
                            y(y < 0.0) = 0.0;
                            
                        case 4
                            % similar to case 2, but replace noisy baseline
                            % points with zeros
                            y( y <= abs(min(y)) ) = 0.0;
                            
                    end
                    
                end
                
                % check that length is still adequate after noise mitigation
                if length(x) < 2
                    
                    if verbose
                        
                        fprintf('ERROR: Noise-mitigated data length is < 2!\n');
                        
                    end
                    
                    center = NaN;
                    
                    return
                    
                else
                    
                    % compute area via trapezoid rule; should be non-negative
                    % but re-check it below
                    denom = trapz(x, y);
                    
                end
                
                % check to see that area is very close to 1.0, and re-normalize
                % for better accuracy
                if (denom > 0.9999 && denom < 1.0001) || (allow_norm && denom > 0.0)
                    
                    % normalize, or re-normalize under strict area checking
                    % when very close to unity
                    y = y/denom;
                    
                    % compute first moment
                    center = trapz(x, y .* x);
                    
                else
                    
                    if verbose
                        fprintf('ERROR: Response area must be unity for first moment computation!\n\n');
                    end
                    
                    center = NaN;
                    
                end
                
            end
            
            return
            
    end
    
    % flag error in method argument if not handled in switch statement above
    disp('ERROR: Method argument must be an integer in the range 1-6. Exiting');
    
    center = NaN;
    
    return
    
end
