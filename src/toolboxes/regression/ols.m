%% ols -- to perform ordinary least squares regression 
%% --------------------------------------------------------
function [beta, alpha, q, ci] = ols(x, y, plot_data, test_data)

	% test data ::
	if nargin == 4
		sz = [1E2, 1]; 
		x = rand(sz); 
		epsilon = normrnd(0, 1E-1, sz);
		y = x + epsilon;
	end

	% get lengths of data ::
	n_x = length(x); 
	n_y = length(y); 
	if n_x ~= n_y
		error('Input vectors are not the same length.');
	else
		n = n_x; 
	end
	dof = n - 2; 

	% calculate sample means ::
	x_mean = sum(x) / n; 
	y_mean = sum(y) / n; 

	% calculate sample (co)variances (*not* taking sample sizes to be sufficiently large, i.e. coefficient is n-1) ::
	S_xx = sum((x - x_mean) .^ 2) / (n - 1); 
	S_yy = sum((y - y_mean) .^ 2) / (n - 1); 
	S_xy = sum((x - x_mean) .* (y - y_mean)) / (n - 1);

	% calculate sample correlation coefficient ::
	rho = S_xy / sqrt(S_xx * S_yy); 

	% calculate sample variance ratio ::
	b = sqrt(S_yy / S_xx); 

	% calculate gradient ::
	beta = rho * b; 

	% calculate intercept ::
	alpha = y_mean - (beta * x_mean); 

	% get fit ::
	y_fit = beta * x + alpha;

	% calculate confidence interval ::
	% a = 1 - p;
	% t = tinv([a / 2, 1 - (a / 2)], dof);
	% se = sqrt(sum((y - y_fit).^2) ./ (n - 2)) .* sqrt((1 / n) + (x - mean(x)).^2 / sum((x - mean(x)).^2));
	% q = se * t; 
	% ci = sqrt((1 - (rho ^ 2)) / (dof * (rho ^ 2)));
	RSS = sum((y - (alpha + (beta * x))) .^ 2); 
	sSquared = RSS / dof; 
	varBeta0 = (sSquared * sum(x .^ 2)) / ((n * sum(x .^ 2)) - (sum(x) ^ 2)); 
	varBeta1 = (n * sSquared) / ((n * sum(x .^ 2)) - (sum(x) ^ 2)); 
	coVar = (- sSquared * sum(x)) / ((n * sum(x .^ 2)) - (sum(x) ^ 2));
	% q = sqrt((varBeta0 + (varBeta1 * (x .^ 2)) + (2 * x * coVar)) / n);
	q = sqrt(varBeta0 + (varBeta1 * (x .^ 2)) + (2 * x * coVar));  % error propogation
	ci = [varBeta1, varBeta0, coVar, n];

	% test ::
	if nargin > 2
		figure; 
		hold('on'); 
		scatter(x, y, 80, 'k', 'filled');
		xlabel('x');
		ylabel('y');
		title('Test Data'); 
		x_tilde = linspace(min(x), max(x), 1E3); 
		y_tilde = (beta * x_tilde) + alpha; 
		y_p = [y_tilde - q; y_tilde + q]; 
		plot(x_tilde, y_tilde, 'color', 'r', 'lineWidth', 2);
		fill([x_tilde, fliplr(x_tilde)], [y_p(1, :), fliplr(y_p(2, :))], 'r', 'faceAlpha', 0.2, 'edgeColor', 'r');
		hold('off')
		legend('Data', 'OLS', 'SD', 'location', 'southEast');
		set(gca, 'box', 'on', 'xLim', [min(x), max(x)]); 
	end

end
