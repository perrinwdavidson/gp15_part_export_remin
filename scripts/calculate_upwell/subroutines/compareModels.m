%%  compare models
%   load and get data ::
w = load([sim_output_basepath 'interpData/upwelling/w_ecco.mat'], 'gp15_w');
i = 1; 
clear('gp15_w'); 
gp15_w{i} = w.gp15_w; 
for dataProduct = {'cglo', 'foam', 'glor', 'oras', 'grep'}

	%   load ::
	i = i + 1; 
	w = load([sim_output_basepath 'interpData/upwelling/w_' dataProduct{1} 'Ave.mat'], 'gp15_wAve');
	gp15_w{i} = w.gp15_wAve;  

	%   get ppz and 100 ::
	

end

%   plot 


%%  end subroutine
