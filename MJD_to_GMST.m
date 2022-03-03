function GMST = MJD_to_GMST(MJD)
    % d is the solar time in days since Jan 1, 2000 at 12h
    d = MJD - 51544.5;
    
    % Calculate GMST in degrees
    GMST = wrapTo360(280.4606 + 360.9856473*d);
end