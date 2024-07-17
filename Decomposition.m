classdef Decomposition < handle & helpers.ArraySupport
    % Decomposition of functions on time-varying cross-sectional domains
    % Low-level class that decomposes function F: Omega \subset R^3 -> R
    % where Omega is a time-varying domain
    properties      
        X struct % Coordinates (t, y, sig, T, Y, Sig)

        H double % Local water depth, depending on (t,y)

        %DFf struct % Orthogonal components of the function to be decomposed
    end
    properties(Dependent)
        A
        B
        sz
    end
    methods
        function obj = Decomposition(varargin)
            obj = obj@helpers.ArraySupport(varargin{:})
        end        
        
        function A = get.A(obj)
             A = obj.wi(obj.H); % Nt by 1 vector
             A = squeeze(A(:,1,1));
        end

        function B = get.B(obj)
            B = abs((obj.X.y(end) - obj.X.y(1)));
        end
        
        function sz = get.sz(obj)
            sz = size(obj.X.T);
        end

        % Elementary functions - trapz or mean
        function tif = ti(obj, f)
            tif = trapz(obj.X.t, f, 1);
        end

        function wif = wi(obj, f)
            wif = trapz(obj.X.y, f, 2);
        end

        function dif = di(obj, f)
            dif = trapz(obj.X.sig, f, 3);
        end

        % Thickness-weighed functions %TODO look in detail.
        function f = repmat3(obj, f)
            % Reshape matrix to size sz , repeating mean values if needed
            rep = [1 1 1];
            if numel(size(f)) == 2
                rep(3) = obj.sz(3);
            end
            rep([size(f)==1]) = obj.sz([size(f)==1]);
            if ~isempty(rep)
                f = repmat(f, rep);
            end
        end


        function [twaf1, twaf2] = twa(obj, f)
            % Taking thickness-weighed average over t
            % Could be vectorized but for now I prefer readability
            % Returns squeezed AND full array

            [f, sf] = obj.handle_input_dim(f);

            if sf(2) == 1 % Laterally averaged variable
                twaf0 = obj.ti(obj.A.*f)./obj.ti(obj.A);
            else
                twaf0 = obj.ti(obj.H.*f)./obj.ti(obj.H); % Just a factor B difference. Yes.
            end
            twaf1 = twaf0(1, 1:sf(2), 1:sf(3));
            twaf2 = obj.repmat3(twaf0);
        end

        function [laf1, laf2] = la(obj,f)
            % Returns squeezed AND full array
            [f, sf] = obj.handle_input_dim(f);
            if sf(1) == 1 % Time-averaged variable
                laf0 = obj.wi(obj.ti(obj.H).*f)./obj.ti(obj.A);
            else
                laf0 = obj.wi(obj.H.*f)./obj.wi(obj.H);
            end
            laf1 = laf0(1:sf(1), 1, 1:sf(3));
            laf2 = obj.repmat3(laf0);
        end

        function [vaf1, vaf2] = va(obj,f)
            % Returns squeezed AND full array
            [f, sf] = obj.handle_input_dim(f);
            vaf0 = obj.di(f)/(obj.X.sig(end)-obj.X.sig(1));
            vaf1 = vaf0(1:sf(1), 1:sf(2), 1);
            vaf2 = obj.repmat3(vaf0);
        end

        function [f, sf] = handle_input_dim(obj,f)
            sf = size(f); % Get original size of input data
            if size(f,3)==1
                sf = [sf 1];
            end
            f = obj.repmat3(f);
        end

        function plot_components(obj, AF, q)
            figure('units','normalized','outerposition',[0 0 1 1])
            for i=1:numel(AF)
                subplot(2,4,i)
                obj.plot_component(AF{i})
%                 ncolor=100;
                colormap(gca, helpers.cmaps(q));
            end
        end

        function plot_component(obj, af) %TODO change
            sf = size(af);
            if size(af,3)==1
                sf = [sf 1];
            end
            if sf(1) == 1 % t-independent
                if sf(2) == 1 % y-independent
                    if sf(3) == 1 % f0
                        bar(af)                                      %1
                        xticklabels("f0")
                        title("F0")
                    else % only sig-dependent
                        plot(squeeze(af), obj.X.sig)    %2
                        xlabel("f(sig)")
                        ylabel("sig")
                        title("F(sig)")
                    end
                else % y-dependent
                    if sf(3) == 1 % only y-dependent
                        plot(obj.X.y, squeeze(af))     %3
                        xlabel("y")
                        ylabel("f(y)")
                        title("F(y)")
                    else % y and sig dependent
                        contourf(obj.X.y, obj.X.sig, squeeze(af)')                        %4
                        xlabel("y")
                        ylabel("sig")
                        title("F(y,sig)")
                        colorbar
                    end
                end
            else % t-dependent
                if sf(2) == 1 % y-independent
                    if sf(3) == 1 % only t-dependent
                        plot(obj.X.t, squeeze(af)')
                        xlabel("t")
                        ylabel("f(t)")%5
                        title("F(t)")
                    else % only (t,sig)-dependent
                        contourf(obj.X.t, obj.X.sig, squeeze(af)')  
                        xlabel("t")
                        ylabel("sig")                        %6
                        colorbar
                        title("F(t,sig)")
                    end
                else % y-dependent
                    if sf(3) == 1 % only (t,y)-dependent
                        contourf(obj.X.t, obj.X.y, squeeze(af)')
                        xlabel("t")
                        ylabel("y")%7
                        title("F(t,y)")
                        colorbar
                    else % t and y and sig dependent
                        histogram(af(:))
                        title("F(t,y,sig)")
                       % disp('Can only display at most 2D variables')
                    end
                end
            end
%             colormap(helpers.cmaps(q))
        end

        function [f,ff] = avg_comp(obj, f, avg_op)
            % Arbitrary composition of averaging operators
            sf0 = size(f);
            if size(f,3)==1
                sf0 = [sf0 1];
            end
            %f0 = f;
            ff = f;
            for idx = 1:numel(avg_op)
                [f, ~] = avg_op{idx}(f);
                [~, ff] = avg_op{idx}(ff);
            end
        end

        function [f, ff] = avg_0(obj, f)
            avg_op = {@obj.twa, @obj.la, @obj.va};
            [f,ff] = obj.avg_comp(f, avg_op);
        end

       % function [DF, AF] = get_decomposition(obj, f)




        function [DF, AF, DFf, AFf] = decompose_function(obj, f)

            AF{1} = f;
            AFf{1}  = f;

            [AF{2}, AFf{2}] = obj.avg_comp(f, {@obj.twa});    %twa(f)
            [AF{3}, AFf{3}] = obj.avg_comp(f, {@obj.la});     %la(f)
            [AF{4}, AFf{4}] = obj.avg_comp(f, {@obj.va});     %va(f)

            [AF{5}, AFf{5}] = obj.avg_comp(f, {@obj.twa, @obj.la});   %twa(la(f))
            [AF{6}, AFf{6}] = obj.avg_comp(f, {@obj.twa, @obj.va});   %twa(va(f))
            [AF{7}, AFf{7}] = obj.avg_comp(f, {@obj.la, @obj.va});    %la(va(f))

            [AF{8}, AFf{8}] = obj.avg_comp(f, {@obj.twa, @obj.la, @obj.va});  %f_0


            DFf{1} = AFf{8};                                DF{1} = AF{8};

            DFf{2} = AFf{7} - DFf{1};   	                DF{2} = DFf{2}(:,1,1);
            DFf{3} = AFf{6} - DFf{1};                       DF{3} = DFf{3}(1,:,1);
            DFf{4} = AFf{5} - DFf{1};                       DF{4} = DFf{4}(1,1,:);

            DFf{5} = AFf{2} - DFf{1} - DFf{3} - DFf{4};     DF{5} = DFf{5}(1,:,:);
            DFf{6} = AFf{3} - DFf{1} - DFf{2} - DFf{4};     DF{6} = DFf{6}(:,1,:);
            DFf{7} = AFf{4} - DFf{1} - DFf{2} - DFf{3};     DF{7} = DFf{7}(:,:,1);

            DFf{8} = AFf{1} - DFf{1} - DFf{2} - DFf{3} - DFf{4} - DFf{5} - DFf{6} - DFf{7};
            DF{8} = DFf{8};
        end

        function [SF, C] = decompose_product(obj, f, g)
            [DF, AF, DFf, AFf] = obj.decompose_function(f);
            [DG, AG, DGf, AGf] = obj.decompose_function(g);

            C = zeros(numel(DFf), numel(DFf));
            SF = zeros([numel(DFf),1]);
            for i = 1:numel(DFf)
                SF(i) = obj.avg_0(DFf{i}.*DGf{i}); % Since we have proven commutativity
            end
            %SF = diag(C);
        end



        function plotC(obj, C)
            figure('units','normalized','outerposition',[0 0 1 1])
            imagesc(C)
            colorbar
            pos=get(gca,'position');
            [rows,cols]=size(C);
            width= 1.05*(pos(3)-pos(1))/(cols-1);
            height = (pos(4)-pos(2))/(rows-1);
            % %create textbox annotations
            for i=1:cols
                for j=1:rows
                    annotation('textbox',[pos(1)+width*(i-1),pos(2)+height*(rows-j),width,height], ...
                        'string',num2str(C(i,j)),'LineStyle','none','HorizontalAlignment','center',...
                        'VerticalAlignment','middle');
                end
            end
            xticklabels(["$u_0$", "$\bar{u}_t^t$", "$\bar{u}_y^y$", "$\bar{u}_\sigma^\sigma$",...
                "$\widehat{u}_{y\sigma}^{y\sigma}$", "$\underline{u}_{t\sigma}^{t\sigma}$", "$[u]_{ty}^{ty}$", "$u_{t y\sigma}$"]);
            yticklabels(["$s_0$", "$\bar{s}_t^t$", "$\bar{s}_y^y$", "$\bar{s}_\sigma^\sigma$",...
                "$\widehat{s}_{y\sigma}^{y\sigma}$", "$\underline{s}_{t\sigma}^{t\sigma}$", "$[s]_{ty}^{ty}$", "$s_{t y\sigma}$"]);

        end

    end
end