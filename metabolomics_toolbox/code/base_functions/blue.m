function c = blue(m)
if nargin <1, m=41; end
c=redblue(m*2);
c=flipud(c(m+1:end,:));