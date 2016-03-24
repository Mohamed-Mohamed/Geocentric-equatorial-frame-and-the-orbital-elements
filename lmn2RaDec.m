function [ RA, Dec ] = lmn2RaDec( r )
% This function is used to get Right ascension and  Declination of a position vector
%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg
%% INPUTS:
% r      :  the position vector
%% OUTPUTS:
% RA  :  Right ascension in degree
%Dec :  Declination in degree
% ---------------------------------------------------------------------------------------------------------------------------------------------------------
mag_r=norm(r);
if mag_r==0
    mag_r=1;
end
l=r(1)/mag_r;
m=r(2)/mag_r;
n=r(3)/mag_r;
Dec=asind(n);
if m >= 0
    RA=acosd(l/cosd(Dec));
elseif m < 0
    RA=360-acosd(l/cosd(Dec));
end
end