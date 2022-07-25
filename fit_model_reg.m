function [model_pars, n_dat, cov_matrix] = fit_model_reg(vel, varargin)
% This function fits the given model to the existing data

n_dat = numel(vel); % number of velocity data
model_mat = [varargin{:}];  % model matrix
n_pars = size(model_mat, 2); % number of model parameters
% handle rank deficient model matrices
model_pars = nan(n_pars, 1);
cov_matrix = nan(n_pars, n_pars);

if rank(model_mat) < n_pars
    return
end
% fit model tot data
cn = cond(model_mat);
p=0;
if cn>100
    disp(strcat('n_data = ', num2str(n_dat), 'while cond(A) = ', num2str(cn)))
    p = 2;
end
G = p*eye(n_pars);
model_pars = (model_mat'*model_mat + G'*G)\model_mat'*vel;
end



        function plot_contour(obj,var)
% Plot the mesh optionally colored with a variable
%
%   plot(obj) plots the mesh with the bed and water surface
%
%   plot(obj,var) plot the mesh and color the cells with the varibale var
%
%   see also: SigmaZetaMesh, plot3
            if ~isscalar(obj)
                for ce=1:numel(obj)
                    subplot(numel(obj),1,ce)
                    plot(obj(ce))
                end
                return
            end
            hold_stat=get(gca,'NextPlot');
            plot(obj.nb_all,obj.zb_all,'k','Linewidth',2)
            hold on
            plot(obj.nw,obj.nw*0+obj.water_level,'b','Linewidth',2)
            plot_var=nan(obj.ncells,1);
            if nargin > 1
                assert(numel(var)==obj.ncells, 'Variable to plot should have same number of elements as cells in the mesh');
                plot_var=var;
            end
            [n,z] = meshgrid(obj.n_middle(obj.col_to_cell)', obj.z_center);
            p  = griddata(obj.n_middle(obj.col_to_cell)',obj.z_center,plot_var,n,z);
%             n = obj.n_middle(obj.col_to_cell)'; % Problem: 
%             [n,z,p] = deal(obj.n_middle(obj.col_to_mat), obj.z_center(obj.cell_to_mat), plot_var(obj.cell_to_mat));
%             contourf(n,z,p,'ShowText','on')
            contourf(n,z,p, 'LineStyle', 'none', 'ShowText','on');
%             clabel(c,h, 0)

            set(gca,'NextPlot',hold_stat);
        end

        function plot_vec(obj,var)
% Plot the mesh optionally colored with a variable, with an optional
% superimposed quiver plot vector field
%
%   plot(obj) plots the mesh with the bed and water surface
%
%   plot(obj,var) plot the mesh and color the cells with the varibale var
%   plot(obj, varx, vary, varz) plot the mesh and color the cells with the
%   varibale varx, with arrows starting at mesh cells and having components
%   vary and varz
%
%   see also: SigmaZetaMesh, plot3
            if ~isscalar(obj)
                for ce=1:numel(obj)
                    subplot(numel(obj),1,ce)
                    plot(obj(ce))
                end
                return
            end
            hold_stat=get(gca,'NextPlot');
            plot(obj.nb_all,obj.zb_all,'k','Linewidth',2)
            hold on
            plot_var=nan(obj.ncells,3);
            if nargin > 1
                assert(size(var,1)==obj.ncells, 'Variable to plot should have same number of elements as cells in the mesh');
                plot_var(:,1)=var(:,1);
            end
            if size(var,2) > 1
%                 assert(numel(vary)==obj.ncells, 'Variable to plot should have same number of elements as cells in the mesh');
                plot_var(:,2)=var(:,2);
            end
            if size(var,2) > 2
%                 assert(numel(varz)==obj.ncells, 'Variable to plot should have same number of elements as cells in the mesh');
                plot_var(:,3)=var(:,3);
            end
            patch(obj.n_patch, obj.z_patch, plot_var(:,1), 'LineStyle', 'None');
            q=quiver(obj.n_middle(obj.col_to_cell)', obj.z_center, plot_var(:,2), plot_var(:,3), .2, 'k');
%             q.ShowArrowHead = 'off';
%             q.Marker = '.';
            plot(obj.nw,obj.nw*0+obj.water_level,'b','Linewidth',2)
            set(gca,'NextPlot',hold_stat);
        end



        function [us,un]=xy2sn_pars(obj, u, v)
            % Transform parameter vectors as obtained as output from the 
            % get_parameters function from xy to sn coordinates
            %
            %   [us,un]=xy2sn(obj,u,v) transform a vector with component 
            %   (u,v) both in projected coordinates to the corresponding 
            %   (s,n) locations with (us,un) components across and along 
            %   the cross-section respectively.
            %
            %   see also: XSection, sn2xy_vel
%                 validateattributes(u,{'numeric'},{});
%                 validateattributes(v,{'numeric'},{});
%                 assert(isequal(size(u),size(v)),'size of u and v should match')
                us = u * obj.direction_orthogonal(1) + v * obj.direction_orthogonal(2);
                un = u * obj.direction(1) + v * obj.direction(2);
        end


        function [u,v]=sn2xy_pars(obj, us, un)
            % Transform parameter vectors as obtained as output from the 
            % get_parameters function from xy to sn coordinates
            %
            %   [us,un]=xy2sn(obj,u,v) transform a vector with component 
            %   (u,v) both in projected coordinates to the corresponding 
            %   (s,n) locations with (us,un) components across and along 
            %   the cross-section respectively.
            %
            %   see also: XSection, sn2xy_vel
%                 validateattributes(u,{'numeric'},{});
%                 validateattributes(v,{'numeric'},{});
%                 assert(isequal(size(u),size(v)),'size of u and v should match')
                u = obj.direction(1) * un + obj.direction_orthogonal(1) * us;
                v = obj.direction(2) * un + obj.direction_orthogonal(2) * us;
        end

        function xs_mat=xyz2snz_mat(obj)
            % Transform vectors from xy to sn coordinates
            %
            %   xs_mat=get_xs_mat(obj) get matrix 3x3 such that u(s,n,z) =
            %   xs_mat*u(x,y,z)
            %
            %   see also: XSection, sn2xy_vel

                xs_mat = [obj.direction_orthogonal(1), obj.direction_orthogonal(2), 0;
                            obj.direction(1), obj.direction(2), 0;
                            0, 0, 1];
        end



