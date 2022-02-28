function rho = get_carrier_phase(ephem, tx_time, ecef_sat, bias_sat, ecef_rcv, bias_rcv, N)

af1 = ephem.SVclockDrift;
L1freq = 1575.42e6;
GPSsecond = tx_time - floor(tx_time/(86400*7))*86400*7;

constant_coeffs = 1;
Klobuchar = get_Klobuchar_coeffs(ephem, tx_time, constant_coeffs);

I = GNSSionosphere(GPSsecond,r_rcv,r_sat,Klobuchar(1,:),Klobuchar(2,:));

rho = (1+af1)*(norm(ecef_sat-ecef_rcv)+bias_rcv-bias_sat) - I + L1freq*N;


end