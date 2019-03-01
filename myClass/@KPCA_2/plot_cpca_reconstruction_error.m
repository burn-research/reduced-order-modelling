function plot_cpca_reconstruction_error(obj, varargin)

y1 = obj.cpca_reconstruction_error_observations;
y2 = obj.cpca_reconstruction_error_variables;

obj.plot_reconstruction_error(y1, y2, 'CPCA reconstruction error');

end