function xcell = decomposeVar(x, U)

%Change this function
mm = U.mesh_mean;  
m = U.mesh; 
ri = mm.row_to_cell; nr = max(ri);
ci = mm.col_to_cell; nc = max(ci);
T = min(13, length(m)); %hour - change this! - dirty workaround

% Depth average and depth-varying.
x0t = NaN(T, nc);
xvt = NaN(T, nc, nr);



for hour = 1:T
    for c = 1:nc
        xc = x(hour, ci == c); % var at time hour, column c. Fixed time and y, varying z.
        x0t(hour, c) = mean(xc,'omitnan'); % depth-averaged var, column c. Scalar.
        xvt(hour, c, 1:length(xc)) = xc - x0t(hour,c);
    end
end
% Place in matrices
    %xDA(hour,:) = SFD{1, hour}.x0'; %t, y - place this in front of function.
    %xDV(hour,:,:) = SFD{1, hour}.uv'; %t, y


% Split tidal mean and varying
x0 = mean(x0t, 1, 'omitnan'); x1 = x0t - x0;

xcell{1,1} = mean(x0,'omitnan');        %Constant in all dimensions - correct
xcell{2,1} = x0 - xcell{1,1};               %only varying in y WA = 0 - correct
xcell{3,1} = mean(x1, 2, 'omitnan');    %only varying in t, but not tidally averaged 0.
xcell{4,1} = x1 - xcell{3,1};               %varying in t and y, but only WA = 0.
xcell{5,1} = xvt;                       %varying in t,y,z


xcell = extendu(xcell);





end