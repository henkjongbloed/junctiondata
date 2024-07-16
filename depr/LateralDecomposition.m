function FT = LateralDecomposition(h,u,s)


time = 1:13;
%(time)
shearTerm = mean(u{5,1}.*s{5,1},3,'omitnan');


FT(1) = h{1}.*u{1,1}*s{1,1};
FT(2) = h{1}.*mean(u{2,1}.*s{2,1},'omitnan');
FT(3) = u{1,1}.*mean(h{2}.*s{2,1},'omitnan');
FT(4) = s{1,1}.*mean(h{2}.*u{2,1},'omitnan');
FT(5) = mean(h{2}.*u{2,1}.*s{2,1},'omitnan');

FT(6) = u{1,1}*mean(h{3}.*s{3,1},'omitnan'); %do tomorrow
FT(7) = mean(h{3}.*mean(u{2,1}.*s{4,1},2,'omitnan'),'omitnan'); %order does not matter

FT(8) = s{1,1}*mean(h{3}.*u{3,1},'omitnan'); %do tomorrow
FT(9) = mean(h{3}.*mean(s{2,1}.*u{4,1},2,'omitnan'),'omitnan');

FT(10) = h{1}.*mean(u{3,1}.*s{3,1},'omitnan'); %seems the main mechanism!
FT(11) = h{1}.*mean(mean(u{4,1}.*s{4,1},'omitnan'),'omitnan');
FT(12) = mean(u{3,1}.*mean(h{2}.*s{4,1},2,'omitnan'),'omitnan');
FT(13) = mean(s{3,1}.*mean(h{2}.*u{4,1},2,'omitnan'),'omitnan');
FT(14) = mean(mean(h{2}.*u{4,1}.*s{4,1},'omitnan'),'omitnan');
FT(15) = mean(h{3}.*u{3,1}.*s{3,1},'omitnan');
FT(16) = mean(h{3}.*mean(u{4,1}.*s{4,1},2,'omitnan'),'omitnan');

FT(17) = h{1}.*mean(mean(shearTerm,'omitnan'),'omitnan');
FT(18) = mean(h{2}.*mean(shearTerm,1,'omitnan'),'omitnan');
FT(19) = mean(h{3}.*mean(shearTerm,2,'omitnan'),'omitnan');

end