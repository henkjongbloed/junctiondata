function u = extendu(u)
u{1,2} = u{1,1} + u{2,1};
u{2,2} = u{3,1} + u{4,1};
u{3,2} = u{5,1};

u{1,3} = u{1,2} + u{2,2};
u{2,3} = u{3,2};
end