function C = Gammacorrect( R ,lambda )
Rmax = max(R(:));
P = R/Rmax;

C = atan(lambda*P)/atan(lambda);


