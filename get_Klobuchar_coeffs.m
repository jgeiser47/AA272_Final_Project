function coeffs = get_Klobuchar_coeffs(ephem, tx_time, constant_coeffs)
if constant_coeffs == 1
    % Vince Thesis
    coeffs(1,:) = [0.2142e-7,0.7451e-8,-0.1192e-6,0];
    coeffs(2,:) = [0.1229e6,0,-0.2621e6,0.1966e6];
else
    % Variable Coefficients
    t = ((ephem.GPSWeek-1356)*7*24*60*60 + tx_time)/(60*60*24*365.25); %Weeks since Jan1 2006
    t = mod(t,1); %nearest beginning of year

    coeffs(1,1) = 1.02e-10*t - 1.9e-7 + 2.5e-9*cos(6.13*t/0.5+602);
    coeffs(1,2) = 1.7215e-8*cos(6.2323*t+105.4373);
    coeffs(1,3) = -5.96e-8;
    coeffs(1,4) = 1.2107e-7*cos(6.2299*t+107.0307);
    coeffs(2,1) = 490.17*t - 8.9656e5 + 0.6523e4*cos(6.3363*t-106.8294);
    coeffs(2,2) = 9.1067e4*cos(6.226*t+130.6530);
    coeffs(2,3) = 3760.1*t - 7.6904e6 - 0.6325e5*cos(6.3761*t-186.8136);
    coeffs(2,4) = 0.5457e6*cos(6.2248*t+117.5348);
end
end