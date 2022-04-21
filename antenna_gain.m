function out_power = antenna_gain(theta_rad,beam_angle_rad,antenna_orientation_rad,N_array)

d_ant_wrt_lambda = 0.5;

theta_rad = wrapToPi(theta_rad-antenna_orientation_rad);
beam_angle_rad = wrapToPi(beam_angle_rad-antenna_orientation_rad);

steering_vec = exp(1j*2*pi*d_ant_wrt_lambda*bsxfun(@times,...
    (0:(N_array-1)).',sin(theta_rad(:)')))/sqrt(N_array);

beamformer = exp(1j*2*pi*d_ant_wrt_lambda*(0:(N_array-1)).'.*sin(beam_angle_rad));

% element gain:

% element_gain = 8-12*(theta_rad(:)'/(65/180*pi)).^2;
% element_gain(element_gain<-30) = -30;
% 
% out_power = reshape(20*log10(abs(beamformer'*steering_vec)) + element_gain,size(theta_rad));

out_power = reshape(20*log10(abs(beamformer'*steering_vec)),size(theta_rad));