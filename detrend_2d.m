
function Z_p = detrend_2d(Z)
    %This function is written by Munther Gdeisat-The General Engineering
    %  Research Institute (GERI) at Liverpool John Moores University.
    % This program is written on 9th October 2011
    %This function is the 2D equivalent of detrend function in Matlab
    %  Z_f = DETREND(Z) removes the best plane fit trend from the
    %     data in the 2D array Z and returns the residual in the 2D array Z_f
    %Thanks for 
    %    http://www.mathworks.co.uk/support/solutions/en/data/1-1AVW5/index.html?solution=1-1AVW5
    if size(Z,2) < 2
        disp('Z must be a 2D array')
        return
    end 
    M = size(Z,2);
    N = size(Z,1);
    [X,Y] = meshgrid(1:M,1:N);
    %Make the 2D data as 1D vector
    Xcolv = X(:); % Make X a column vector
    Ycolv = Y(:); % Make Y a column vector
    Zcolv = Z(:); % Make Z a column vector
    Const = ones(size(Xcolv)); % Vector of ones for constant term
    % find the coeffcients of the best plane fit
    Coefficients = [Xcolv Ycolv Const]\Zcolv; % Find the coefficients
    XCoeff = Coefficients(1); % X coefficient
    YCoeff = Coefficients(2); % X coefficient
    CCoeff = Coefficients(3); % constant term
    % detrend the data
    Z_p = XCoeff * X + YCoeff * Y + CCoeff;

end
