function NMISE_y_hat=fnmise(original,y_hat)


    vec_original = original(:);
	zeroes = find(vec_original == 0);
	vec_original(zeroes) = [];

	vec_y_hat = y_hat(:);
	vec_y_hat(zeroes) = [];

	% Normalized mean integrated square error of the estimate y_hat
	N = numel(vec_original);
	NMISE_y_hat = sum((vec_y_hat-vec_original).^2 ./ vec_original) / N;