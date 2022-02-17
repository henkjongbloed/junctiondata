function msh=comp_plotQ(msh)

msh = compQ_v2(msh);

for ct=1:3
    msh{ct}.Qloc = zeros(size(squeeze(msh{ct}.pars(:,:,:,1))));
    msh{ct}.Q = nan(size(squeeze(msh{ct}.pars(1,1,:,1))));
    msh{ct}.adcp_good = 1:13;
    for hour=1:13
        %x,y,z coordinates
        vdotn = squeeze(msh{ct}.pars(:,:,msh{ct}.adcp_good(hour),1))  * msh{ct}.Nvec(1) +...
            squeeze(msh{ct}.pars(:,:,msh{ct}.adcp_good(hour),2))  * msh{ct}.Nvec(2) +...
            squeeze(msh{ct}.pars(:,:,msh{ct}.adcp_good(hour),3))  * 0;
        Qloc = msh{ct}.A(:,:,1).*vdotn;
        msh{ct}.Qloc(:,:,msh{ct}.adcp_good(hour)) = Qloc;
        msh{ct}.Q(msh{ct}.adcp_good(hour)) = nansum(Qloc,[1,2]);
        %          msh{ct}.Qloc(:,:,msh{ct}.adcp_good(hour)) = Qloc;
        Qlocrot = msh{ct}.sec.pars(:,:,msh{ct}.adcp_good(hour),1).*msh{ct}.A(:,:,1);
        msh{ct}.sec.Qloc(:,:,msh{ct}.adcp_good(hour)) = Qlocrot;
        msh{ct}.sec.Q(msh{ct}.adcp_good(hour)) = nansum(Qlocrot,[1,2]);
        Qlocrot_S = Qlocrot.*msh{ct}.CTD_int.Salinity_ex(:,:,hour);
        msh{ct}.sec.QlocS(:,:,msh{ct}.adcp_good(hour)) = Qlocrot_S;
        msh{ct}.sec.QS(msh{ct}.adcp_good(hour)) = nansum(Qlocrot_S,[1,2]);
    end
end


%sec.pars(1) is the same as dot(v(x,y,z,), n)
plotQ(msh)