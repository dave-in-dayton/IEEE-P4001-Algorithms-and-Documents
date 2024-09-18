function [width] = funcwid_P4001(x, y, chanwidth, method)
    %FUNCWID Computes width of function y = f(x) by specified method
    %
    %   x,y       = input data
	%   chanwidth = arbitrary integer
    %
    %   method = 1  FWHM
    %   method = 2  CIE 233 with peak ordinate center definition
    %   method = 6  Std. Dev. scaled to Gaussian FWHM per P4001 concept
    %   method = 8  CIE 233 with convolution peak center definition
    %   method = 10 AFRL 76% area algorithm with AFRL median
    %
    % Uses: FWHM_V2_P4001.m, and funcctr_P4001.m
    
    % D. Perry, Leidos, Dayton, OH
    % V1.0 August 2024
    
    % set verbosity - can make this a function argument if desired
    verbose = false;
    
    % check for valid x and y inputs
    if isempty(x) || isempty(y) || length(x) ~= length(y)
        width = NaN;
        return
    end
    
    % make sure x and y are row vectors and transpose as needed
    if ~isrow(x); x = x'; end
    
    if ~isrow(y); y = y'; end
    
    % compute function center or peak for width methods that require these
    % quantities
    switch method
        
        case 2
            % use peak value obtained from FWHM_V2.m
            [~,~,~,ypeak] = FWHM_V2_P4001(x,y);
                       
        case 6
            % std dev scaled to normal FWHM - prelimnary calc of scale
            % factor
            outputsf = 2.0 * sqrt(2.0 * log(2.0));
       
        case 8
            % CIE alternate - uses width method 2 algorithm, but
            % now obtains ypeak via the conv. peak center algorithm
             conv_peak_ctr = funcctr_P4001(x, y, chanwidth, 5);
            
            [~,xind] = min(abs(x-conv_peak_ctr));
             
            ypeak = y(xind);
      
            method = 2;
            
        case 10
            % AFRL 76% area algorithm with AFRL median
            center = funcctr_P4001(x, y, 1.0, 4);
            
    end
    
    % check for function valid center if it was computed above
    if exist('center', 'var')
        if isnan(center)
            width = NaN;
            return
        end
    end
    
    % now compute specified width metric
    switch method
        
        case 1
            
            % FWHM from actual response values
            %[~, ~, width] = FWHM(x,y);
            [~, ~, width, ~] = FWHM_V2_P4001(x,y);
            
            return
            
        case 2
            
            % CIE 233 method - area divided by user-supplied center ordinate
            
            % implement area algorithm via the trapezoid rule, accounting for
            % the possibility of un-equal point spacing by using the two-
            % argument form of the MATLAB trapz() function
            if ypeak ~= 0.0
                width = trapz(x,y) / ypeak;
            else
                width = NaN;
            end
            
            return
 
        case 6
            
            % standard deviation; one must use positive y values (since we
            % are computing from a distribution rather than observed data),
            % so specify the mitigation method desired, noting that method
            % 7 and then 3 have been shown to provide the highest
            % flexibility scores
            %
            % method 7 does not change the data in any way and only flags
            % detection of a negative variance result. Code for methods 5
            % and 6 has been deleted for brevity and hence function the
            % same as method 7
            mitig_mthd = 7;
            
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
                        % similar to case 2, but replace near-baseline points
                        % with zeros
                        y( y <= abs(min(y)) ) = 0.0;
                        
                    case 5
                        % find maximum symmetric domain about the mean that
                        % just avoids the first zero crossing
                        
                    case 6
                        % similar to previous method, but without symmetry
                        % requirement
                                                
                    case 7
                        % do nothing, simply check the sign of the variance
                        % below
                end
                
            end
            
            % compute normalized second central moment via one of four
            % similar definitions, and using funcctr.m method 6 for the
            % mean; note that the mean can fail if the there are less than
            % two points or the area under y is non-positive
            ctr_alg = 6;
            
            mu = funcctr_P4001(x,y,1,ctr_alg);
            
            def_code = 1;
            
            if ~isnan(mu)
                
                if def_code < 3
                    % normalize by area,
                    denom = trapz(x,y);
                else
                    % or otherwise normalize y by its sum
                    denom = sum(y);
                end
                
                if denom > 0
                    ynorm = y / denom;
                else
                    % exit if denom term is non-positive
                    if verbose
                        fprintf('ERROR: Sample area or sum must be positive when normalizing!\n\n');
                    end
                    
                    width = NaN;
                    
                    return;
                end
                
                % apply specified definition of variance
                switch def_code
                    
                    case 1
                        % probability function method
                        vari = trapz(x, (x - mu).^2 .* ynorm);
                    case 2
                        % normalized moments method
                        vari = trapz(x, x .* x .* ynorm) - mu^2;
                    case 3
                        % similar to case 1 but with summation
                        vari = sum((x - mu).^2 .* ynorm);
                    case 4
                        % similar to case 2 but with summation
                        vari = sum(x .* x .* ynorm) - mu^2;
                        
                end
                
                % check for negative variance before taking square root that
                % could occur due to numeric imprecision effects or noise
                if vari >=0
                    width = outputsf * sqrt(vari);
                else
                    
                    if verbose
                        fprintf('ERROR: Negative Variance!\n');
                    end
                    
                    width = NaN;
                    
                end
                
            else
                
                % mean calculation failed
                if verbose
                    fprintf('ERROR: Cannot compute mean for 2nd moment calc.!\n');
                end
                
                width = NaN;
                
            end
            
            return
                        
        case 10
            
            % AFRL 76% area algorithm
            
            % first create two x-y sequences that start at the median, checking
            % to see first if it is located at a current x-value within a
            % reasonable tolerance based on knowledge of the max sampling
            % rate expected  
             [minval, ctr_ind] = min(abs(x-center));
            
            if abs(minval) < 0.0001
                ctr_at_x_value = true;
                center = x(ctr_ind);
            else
                ctr_at_x_value = false;
            end

            if ctr_at_x_value
                x1 = x(1:ctr_ind);
                y1 = y(1:ctr_ind);
                
                x2 = x(ctr_ind:end);
                y2 = y(ctr_ind:end); 
            else
                % use linear interpolation to estimate y value at the median
                y_med = interp1(x,y,center);
                
                ind1 = find(x < center);
                ind2 = find(x > center);
                
                x1 = [x(ind1), center];
                y1 = [y(ind1), y_med];
                
                x2 = [center,  x(ind2)];
                y2 = [y_med,   y(ind2)];
            end
            
            % get cumulative distributions for each half of the
            % distribution above and below the median and normalize
            x1_rel = fliplr(center - x1);
            x2_rel = x2 - center;

            cum1 = cumtrapz(x1_rel,fliplr(y1));
            cum2 = cumtrapz(x2_rel,y2);
            
            totarea = cum1(end) + cum2(end);
            
            cum1 = cum1/totarea;
            cum2 = cum2/totarea;
            
            % find total distance where the 76% area sum is reached
            Lmax = min(length(cum1), length(cum2));
            
            areasum = cum1(1:Lmax) + cum2(1:Lmax);
            
            x_dist = x2_rel(1:Lmax) + x1_rel(1:Lmax);
            
            fullarea = 0.7609681085504878;
            
            % must subset the data to unique areasum values - since areasum
            % is in theory monotonically increasing, use all points that
            % are below the fullarea value plus one additional point that
            % allows interpolation
            L = find(areasum > fullarea, 1, 'first');
            
            % be sure areasum will allow interpolation, since in particular, it
            % may not quite reach the fullarea magnitude due to noise, etc.
            if fullarea < areasum(1) || fullarea > areasum(end) || L < 2 || isempty(L)
                width = NaN;
            else
                width = interp1(areasum(1:L),x_dist(1:L),fullarea);
            end
            
            return
            
    end
    
    % flag error in method argument if not handled in switch statement above
    disp('ERROR: Method argument must be 1, 2, 6, 8 or 10. Exiting');
    
    width = NaN;
    
    return
    
end