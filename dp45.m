function [t_out, y_out] = dp45( f, t_rng, y0, h, eps_abs )

    if ~isa( f, 'function_handle' )
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument f is not a function handle' ) );
    end
    if ~isscalar( eps_abs) || (eps_abs <= 0) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument eps_abs is not a positive scalar' ) );
    end
    if ~isscalar( y0 ) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument y0 is not a scalar' ) );
    end
    if ~isscalar( h ) || (h <= 0) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument h is not a positive scalar' ) );
    end
    if ~all( size( t_rng) == [2, 1] ) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument t_rng is not a 2-dimensional column vector' ) );
    end

    A = [ 0          0         0        0         0        0   0; 
          1          0         0        0         0        0   0;
         1/4        3/4        0        0         0        0   0;
        11/9      -14/3      40/9       0         0        0   0;
      4843/1458 -3170/243  8056/729  -53/162      0        0   0;
      9017/3168  -355/33  46732/5247  49/176 -5103/18656   0   0;
        35/384       0      500/1113 125/192 -2187/6784  11/84 0]';
 
    by = [5179/57600 0 7571/16695 393/640 -92097/339200 187/2100 1/40]';
    bz = [  35/384   0  500/1113  125/192  -2187/6784    11/84    0]';
 
    c = [0 1/5 3/10 4/5 8/9 1 1]';
    
    % Initialize t_out and y_out by giving an inital scalar value of y and
    % t
    n_D = 7;
    D = zeros( 1, n_D );
    
    t0 = t_rng(1);
    tf = t_rng(2);
    t_out(1) = t0;
    y_out(1) = y0;
    
    
    k = 1;
    while t_out(d) < tf
        %Use Dormand Prince to find two approximations
        %User Euler's and Heun's method to approximate y(t) at
        %t = t_out(k) + h for the current value of h
        for m = 1:n_D
            D(m) = f(t_out(d) + h*(c(m)), y_out(d) + h*(c(m))*(D)*(A(:,m)) );
        end
 
        y_tmp = y_out(d) + (h)*(D)*(by);
        z_tmp = y_out(d) + (h)*(D)*(bz);
        
        %Calculating the scaling factor to determine the value of h
        %depending on the interval. Calculated h such that the error is
        %less than the relative width of h from t-initial to t-final:
        s = ((h*(eps_abs))/((2)*(tf-t0)*abs(y_tmp - z_tmp)))^(0.25);
    
        %checks to see if the scaling factor is greater than or equal to 2
        %and will double the value of h:
        if s >= 2
        
            y_out(d+1) = z_tmp;
            t_out(+1) = t_out(k) + h;
            h = h*2;
            k = k+1;
        
        %checks to see if the the scaling factor is between 1 and 2 and
        %will leave the value of h unchanged:
        elseif s >= 1
    
            y_out(d+1) = z_tmp;
            t_out(d+1) = t_out(d) + h;
            d = d+1;
        
        %checks to see if the value of the scaling factor is less than one
        %and will divide h by 2 and try again:
        else s < 1
  
            h = h/2;
            
        end
    
        %We must make one final check before we end the loop to see if we
        %have to reduce the value of h
        if t_out(d) + h > tf
            h = tf - t_out(d);
        end
    
    end
end
