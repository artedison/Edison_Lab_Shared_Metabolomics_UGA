%############################################################################
%#
%#                          function DIMENSION
%#
%#      calculates the dimensionality of the data (not like ndims)
%#
%#      usage: dim=dimension(data);
%#
%#      result:
%#              dim==0  data is scalar
%#              dim==1  data is 1D vector
%#              dim==2  data is 2D
%#              dim>2   data is nD (n>2)
%#     
%#      (c) P. Blümler 1/03
%############################################################################
%----------------------------------------------------------------------------
%  version 1.1 PB 21/1/03    (please change this when code is altered)
%----------------------------------------------------------------------------



function dim=dimension(data)

data=squeeze(data);
nd=ndims(data);
ns=size(data);
nl=length(data);
if nd==2
    if nl~=prod(ns)
        dim=2;
    elseif nl==1
        dim=0;
    else
        dim=1;
    end
else
    dim=nd;
end
