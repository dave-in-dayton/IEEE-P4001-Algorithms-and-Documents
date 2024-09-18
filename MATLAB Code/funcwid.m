function [width] = funcwid(x, y, chanwidth, method)
    %FUNCWID Computes width of function y = f(x) by specified method
    %
    %   x,y    = input data
    %
    %   method = 1  FWHM
    %   method = 2  CIE 233 with peak ordinate center definition
    %   method = 3  AFRL version of DLR's equivalent Guassian area technique
    %   method = 4  FWHM after the convolution g(x) = f(x) * rect(x/chanwidth)
    %   method = 5  Standard deviation (as second central moment of y = p(x) )
    %   method = 6  Method 5, but scaled to Gaussian FWHM per P4001 concept
    %   method = 7  DLR width similar to (3) but with 3rd-order B-Spline interp
    %   method = 8  CIE 233 with convolution peak center definition
    %   method = 9  DLR "resolution" (and median) via DLR-supplied Python module
    %   method = 10 AFRL version of DLR width, via AFRL version of DLR median
    %
    % Uses: FWHM.m or FWHM_V2.m, funcctr.m, Curve-Fitting Toolbox
    %
    % For more information, see working notes dated March-April, 2022,
    % related IEEE WG discussion notes, Easton's Fourier Methods in
    % Imaging, Frieden's Probability, Statistical Optics, and Data Testing,
    % and Palmer's 1980 and 1983 papers on LandSat radiometric bandwidth
    % normalization (key cited methods there are the equivalent of the
    % CIE233, FWHM, and standard deviation.
    %
    % D. Perry, Leidos, Dayton, OH
    %
    % V1.0 April    2022
    % V1.1 November 2022  - added method 6 (output scaling for method 5) and
    %                       edited comments
    % V1.2 December 2022  - added method 7, improved version of AFRL algorithm 3
    % V1.3 January  2023  - removed center input argument and assigned the peak
    %                       ordinate center for method 2; added
    %                       method 8, alternate CIE center
    % V1.4 February 2023  - added DLR's Python resolution calc
    % V1.5 September 2023 - updated function call to FWHM_V2
    % V1.6 February 2024  - correct coding error for CIE method that
    %                       employs conv center and associated ypeak
    
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
            % old way: CIE center is obtained from companion function
            % funcctr.m, and then the peak ordinate is found
            %peak_center = funcctr(x, y, 1.0, 1);        
            
            %[~,xind] = min(abs(x-peak_ctr));
             
            %ypeak = y(xind);
            
            % new approach, 9-25-23 - use peak value directly obtained
            % from updated form of FWHM_V2.m
            [~,~,~,ypeak] = FWHM_V2(x,y);
            
        case 3
            % AFRL version of DLR width using AFRL version of DLR median
            center = funcctr(x, y, 1.0, 4);
            
        case 5
            % default std deviation output is unscaled
            outputsf = 1.0;
            
        case 6
            % std dev scaled to normal FWHM; changes to method 5 code spec
            outputsf = 2.0 * sqrt(2.0 * log(2.0));
            
            method = 5;
            
        case 7
            % improved DLR width method, using improved DLR median, both of
            % which employ B-Spline interpolation
            center = funcctr(x, y, 1.0, 7);
            
        case 8
            % CIE alternate - uses width method 2 algorithm below, but
            % obtains ypeak via the conv. peak center
            %
            % it was noted in January 2023 trials of all other center types
            % that CIE can be sensitive to the denominator, since with
            % bi-modal functions, the function center is not at the max
            % function value, and can in fact be small, hence amplifying
            % the effect of noise and sampling errors on center
            % determination
            conv_peak_ctr = funcctr(x, y, chanwidth, 5);
            
            [~,xind] = min(abs(x-conv_peak_ctr));
             
            ypeak = y(xind);
      
            method = 2;
            
        case 10
            % AFRL version of DLR width using AFRL version of DLR median
            center = funcctr(x, y, 1.0, 4);
            
    end
    
    % check for function valid center if it was computed above
    if exist('center', 'var')
        if isnan(center)
            width = NaN;
            return
        end
    end
    
    % select optional plotting for method 4
    do_plot = false;
    
    
    switch method
        
        case 1
            
            % FWHM from actual response values
            %[~, ~, width] = FWHM(x,y);
            [~, ~, width, ~] = FWHM_V2(x,y);
            
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
            
            
        case 3
            
            % AFRL/DLR Gaussian equivalent area method - needs center abcissa to
            % solve for the delta x's that give a normalized area of 0.761
            % (i.e., the area under a gaussian curve over the range of its FWHM
            % x-axis limits); here the supplied center would typically be from
            % the companion AFRL/DLR center method (i.e., the median, or 50%-
            % probability abscissa)
            
            % get cumulative distribution, using the trapezoidal rule -
            % obviates the need for equal x-axis point spacing
            %
            % note 4-28-23 - revised the function call to correct the
            % earlier error where the x values should have been used; not
            % tested...
            cumdist = cumtrapz(x,y);
            
            cumdist = cumdist / cumdist(end);
            
            % find nearest-neighbor x index for center
            [~, xind] = min(abs(x-center));
            
            if isnan(xind)
                width = NaN;
                return
            end
            
            % now get abcissa value for the lower half-area region bounded by
            % the center (2-2023: updated to DLR's DP value)
            fullarea = 0.7609681085504878;
            
            %         % original approach - incorrect, since the area does not have to
            %         % be equal on each side of center!
            %         halfarea = fullarea/2.0;
            %         areaatctr = cumdist(xind);
            %         goalarea = areaatctr - halfarea;
            %         [~,lowxind] = min(abs(cumdist-goalarea));
            %         % and then locate upper abcissa value using full area
            %         goalarea = cumdist(lowxind) + fullarea;
            %         [~,highxind] = min(abs(cumdist-goalarea));
            
            % new 1-18-2023 - now use a simple index search method (tried fzero
            % here, but had trouble where endpoint search value never changed
            % sign, required by the method)
            myfun = @(deltind)  cumdist(xind+round(deltind)) - ...
                cumdist(xind-round(deltind)) - fullarea;
            
            maxdelind = min([xind-1,length(x)-xind]);
            
            exitflag = -1;
            
            for iii = 0:maxdelind
                if myfun(iii) > 0
                    exitflag = 1;
                    halfwidth = iii;
                    break
                end
            end
            
            % set width value if normal exit condtion is not detected
            if exitflag ~= 1
                
                width = NaN;
                
                if verbose
                    
                    fprintf('Abnormal termination during AFRL/DLR width computation!\n\n');
                    
                end
                
            else
                % assign width value
                highxind = xind + round(halfwidth);
                
                lowxind  = xind - round(halfwidth);
                
                width    = x(highxind)-x(lowxind);
            end
            
            % optional check on result (use a wide threshold to pass most
            % results for later assessments of accuracy - just trying to catch
            % gross errors here!)
            do_check = false;
            
            if do_check && ~isnan(width)
                
                estarea = cumdist(highxind) - cumdist(lowxind);
                
                if abs(estarea - fullarea) >  fullarea * 0.5
                    
                    if verbose
                        fprintf('Accuracy check for AFRL/DLR width method failed.\n');
                        fprintf('Returning NaN for function width.\n\n');
                    end
                    
                    width = NaN;
                    
                end
                
            end
            
            return
            
        case 4
            
            % width of convolution with a rect function of one spectral channel
            % width
            
            % first verify that x point spacing is OK
            deltaxvec = x(2:end) - x(1:end-1);
            
            deltax = median(deltaxvec);
            
            if (max(deltaxvec) - min(deltaxvec) ) / deltax > 0.0001
                if verbose
                    fprintf('ERROR - x-data for box center method is NOT evenly spaced.\n');
                    fprintf('Returning NaN for function width.\n\n');
                end
                width = NaN;
                return
            end
            
            % set up box/rect function
            numchansamples = (chanwidth/deltax) + 1;
            
            boxfct = ones(1,round(numchansamples));
            
            % perform convolution and find the width about its maximum
            convresult = conv(y,boxfct,'same');
            
            [~, ~, width, ~] = FWHM_V2(x,convresult);
            
            % optional plot
            if do_plot
                
                maxval = max(convresult);
                
                convresult = convresult/maxval;
                
                boxplotdata = zeros(size(y));
                
                plotoffset = round(0.1*length(x));
                
                boxplotdata(1,plotoffset:plotoffset+length(boxfct)-1) = boxfct;
                
                figure; plot(x,y,x,boxplotdata,x,convresult); grid on;
                
                title('Box Function Width Algorithm Components');
                
                legend({'Input','Box Function','Conv. of Input and Box'});
                
            end
            
            return
            
            
        case 5
            
            % standard deviation; one must use positive y values (since we are
            % computing from a distribution rather than observed data), so
            % specify the mitigation method desired, noting that choice 7 does
            % not change the data in any way and only flags detection of a
            % negative variance result
            %
            % method 7 and then 3 have been shown to provide the highest
            % flexibility scores - See 'Noise Mitigation Study - For Standard
            % Deviation Width Algorithm.xlsx', second tab, 4-13-2023; later
            % work showed that method 7 extends the SNR range when working
            % with the difficult fat tails shape
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
                        
                        % get function mean and work outward from that
                        % x-location if it exists and is positive
                        mu = funcctr(x,y,1,6);
                        
                        % check for valid center and that the ordinate is
                        % positive
                        mu_status = false;
                        
                        if ~isnan(mu)
                            
                            [~, ctr_ind] = min(abs(x-mu));
                            
                            if y(ctr_ind) > 0
                                mu_status = true;
                            end
                            
                        end
                        
                        if mu_status
                            
                            ctr_offset = length(x) - ctr_ind;
                            
                            if ctr_offset < ctr_ind - 1
                                max_search_delta = ctr_offset;
                            else
                                max_search_delta = ctr_ind - 1;
                            end
                            
                            % set default lower/upper indices then over-write
                            pos1 = ctr_ind - max_search_delta;
                            pos2 = ctr_ind + max_search_delta;
                            
                            for ii = 1:max_search_delta
                                pos1 = ctr_ind - ii;
                                pos2 = ctr_ind + ii;
                                if y(pos1) < 0.0 || y(pos2) < 0.0
                                    pos1 = pos1 + 1;
                                    pos2 = pos2 - 1;
                                    break
                                end
                            end
                            
                            % get new x, y vectors where all values are
                            % positive and contiguous
                            x = x(pos1:pos2);
                            y = y(pos1:pos2);
                            
                            % for diagnostic use only: double check that all
                            % values are positive
                            %                         if min(y) < 0
                            %                             fprintf('ERROR: Cannot eliminate all negative values from std. dev. calc.!\n');
                            %                         end
                            
                            if length(x) < 2
                                
                                if verbose
                                    fprintf('ERROR: Noise-mitigated data length is < 2!\n');
                                end
                                
                                width = NaN;
                                
                                return
                                
                            end
                            
                        else
                            
                            if verbose
                                fprintf('ERROR: Cannot compute function mean for 2nd moment!\n');
                            end
                            
                            width = NaN;
                            
                            return
                            
                        end
                        
                    case 6
                        % similar to previous method, but without symmetry
                        % requirement
                        
                        % create two search ranges below and above the mean
                        % and locate last and first zero crossing below/above
                        % it
                        mu = funcctr(x,y,1,6);
                        
                        % check for valid center and that the ordinate is
                        % positive
                        mu_status = false;
                        
                        if ~isnan(mu)
                            
                            [~, ctr_ind] = min(abs(x-mu));
                            
                            if y(ctr_ind) > 0
                                mu_status = true;
                            end
                            
                        end
                        
                        if mu_status
                            
                            [~, ctr_ind] = min(abs(x-mu));
                            
                            % check lower portion of data
                            if min(y(1:ctr_ind)) < 0.0
                                
                                negvalind = y(1:ctr_ind) < 0.0;
                                pos1 = find(negvalind,1,'last');
                                pos1 = pos1 + 1;
                                
                            else
                                pos1 = 1;
                            end
                            
                            % check upper portion of data
                            if min(y(ctr_ind:end)) < 0.0
                                
                                negvalind = y(ctr_ind:end) < 0.0;
                                pos2 = find(negvalind,1,'first');
                                pos2 = pos2 + ctr_ind - 2;
                                
                            else
                                pos2 = length(x);
                            end
                            
                            % get new x, y vectors where all values are
                            % positive and contiguous
                            x = x(pos1:pos2);
                            y = y(pos1:pos2);
                            
                            % for diagnostic use only: double check that all
                            % values are positive
                            %                         if min(y) < 0
                            %                             fprintf('ERROR: Cannot eliminate all negative values from std. dev. calc.!\n');
                            %                         end
                            
                        else
                            
                            if verbose
                                fprintf('ERROR: Cannot compute function center for 2nd moment!\n');
                            end
                            
                            width = NaN;
                            
                            return
                            
                        end
                        
                    case 7
                        % do nothing, simply check the sign of the variance
                        % below
                end
                
            end
            
            % compute normalized second central moment via one of four similar
            % definitions, and using funcctr.m method 6 for the mean; note
            % that the mean can fail if the there are less than two points or
            % the area under y is non-positive; early tests used the centroid,
            % while the current choice is the first moment
            ctr_alg = 6;
            
            mu = funcctr(x,y,1,ctr_alg);
            
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
            
        case 7
            
            % first make sure there are enough points to construct spline
            if length(x) > 2
                
                % fit data with 3rd-order B-Spline
                spline3 = spapi(3,x,y);
                
                % check to see that a structure was returned (a better check
                % based on its contents should be added?)
                if isstruct(spline3)
                    
                    % compute integral via spline coefficients
                    integ   = fnint(spline3);
                    
                    % get full area and then scaled area between FWHM's
                    full_area = fnval(integ,x(end));
                    
                    fwhm_area = 0.760968* full_area;
                    
                    % set up search limits based on
                    
                    % then iterate about center for area goal equal to that
                    % between 1-sigma points on a normal PDF; may be able to
                    % improve search time by providing "reasonable" search
                    % limits on deltax
                    myfun = @(deltax) fnval(integ,center+deltax) - ...
                        fnval(integ,center-deltax) - fwhm_area;
                    
                    [halfwidth, ~, exitflag, ~] = fzero(myfun,0);
                    
                    % set width value if abnormal exit condtion is not detected
                    if exitflag ~= 1
                        
                        width = NaN;
                        
                        if verbose
                            
                            fprintf('Abnormal termination of fzero during width computation!\n\n');
                            
                            %fprintf('fzero output message: %s\n',outputmsg.message);
                            
                        end
                        
                    else
                        % assign width value
                        width = 2* halfwidth;
                    end
                    
                else
                    
                    width = NaN;
                    
                    if verbose || true
                        
                        fprintf('Could not complete spline fit during center computation!\n\n');
                        
                    end
                    
                end
                
            else
                width = NaN;
            end
            
            return
            
        case 9
            % DLR-supplied center and width using MATLAB's Python execution
            % mode; must load Python module prior to entering this routine
            
            % first make sure there are enough points to construct the spline,
            % assuming for now that 2 points are sufficient, as seen with the
            % AFRL proxy method above
            if length(x) > 2
                
                % 3-22-23 - now calling Python interface function with revised
                % argument list
                ctr_only = false;
                
                g = py.IEEE.peakspline_afrl.fit_peakspline_via_matlab(x,y,ctr_only);
                
                % convert Python object members to MATLAB cell arrays and then
                % floating point MATLAB scalars
                g1 = g(1).cell;
                g2 = g(2).cell;
                
                DLRctr = g1{1};
                DLRwid = g2{1};
                
                width = DLRwid;
                
            else
                width = NaN;
            end
            
            return
            
        case 10
            % alternate implementation of method 3 above, now using linear
            % interpolation during the center ordinate and 76% area
            % determination processes
            
            % first create two x-y sequences that start at the median, checking
            % to see first if it is located at a current x-value within a
            % reasonable tolerance based on knowledge of the max sampling
            % rate expected
            
            % note 4-27-23: significant revision this date due to accuracy
            % sensitvity observed with minor changes in function shapes
            % such as the fat tails test function
            %
            % another (faster?) alternative to method 7 might be to use the
            % 'fit' function to approximate the y-data, integrate it with
            % 'integrate', and then linearly interpolate the integral about
            % the median to solve for the 76% area region
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
            %
            % note 4-27-23: revised to correct an error in the use of
            % cumtrapz - MUST use the (X,Y) version, NOT the (X) or unit
            % spacing basic option!! This is important when the median does
            % not fall at one of the x values
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
            
            fullarea = 0.7609681085504878; % from DLR Python code
            
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
    disp('ERROR: Method argument must be an integer in the range 1-9. Exiting');
    
    width = NaN;
    
    return
    
end