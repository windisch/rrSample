HemmeckeMatrix <- function(n)
{

I=diag(n);
O=matrix(rep(0,n*n),n,n);
v=matrix(rep(1,n),1,n);
w=matrix(rep(0,n),1,n);

rbind(t(rbind(I,I,O,O,-v,w)),t(rbind(O,O,I,I,w,-v)),t(rbind(t(w),t(w),t(w),t(w),1,1)));
}
