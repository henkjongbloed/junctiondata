function H = decomposeH(U)


mm = U.mesh_mean;  
m = U.mesh; 

eta = U.eta(1:13); %Dirty workaround
eta0 = mean(eta); 
eta1 = eta - eta0;
zb = mm.zb_middle;
h0 = eta0 - zb;



%Split water depth: Instantaneous depth = H{1} + H{2} + H{3}
H{1} = mean(h0);
H{2} = h0 - H{1};
H{3} = eta1;

end


