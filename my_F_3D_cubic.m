function F = my_F_3D_cubic(v,RR,NN)
% define free volume

%26 vertices
%18 boxes
%24 faces

a_vec = v.^(1/3);

F = zeros(1,numel(a_vec));
for i=1:numel(a_vec)
    a = a_vec(i);
    
    if (a < 2*RR)
        F(i) = 0;
    elseif (a <= sqrt(6)*RR)
        long_double_int = pi/12*(4*(2*RR)+sqrt(2)*a).*(2*(2*RR)-sqrt(2)*a).^2;
        
        w = sqrt((2*RR)^2*12*a^4-8*a^6);
        q1 = 2*sqrt(2)*a^3;
        triple_int = w/6-3*a/sqrt(2)*(2*(2*RR)^2-a^2/3)*myatan(2*w/q1)+...
            +4*(2*RR)^3*myatan(sqrt(2)*a*w/(2*RR*q1));
        
        solid_angle = 4*asin(1/3)/(4*pi);
        dihedral_angle = acos(-1/3)/(2*pi);
        
        F(i) = 4*a^3/3-6*solid_angle*4/3*pi*(2*RR).^3+...
            +12*long_double_int*dihedral_angle-4*triple_int;
    elseif (a <= 2*sqrt(2)*RR)
        double_int = pi/12*(4*(2*RR)+a).*(2*(2*RR)-a).^2;
        long_double_int = pi/12*(4*(2*RR)+sqrt(2)*a).*(2*(2*RR)-sqrt(2)*a).^2;
        
        w = sqrt((2*RR)^2*4*a^4-2*a^6);
        q1 = 2*a^3;
        triple_int = w/6-a*(2*(2*RR)^2-a^2/6)*myatan(2*w/q1)+...
            -a*sqrt(2)/2*(2*(2*RR)^2-a^2/3)*pi/2+...
            +4*(2*RR)^3/3*(myatan(a*w/(2*RR*q1))+pi/2)+...
            +2*(2*RR)^3/3*(2*myatan(a*w/(2*RR*q1)));
        
        quad_int = 2*triple_int-long_double_int;
        
        F(i) = 8*v(i)-7*4/3*pi*(2*RR).^3+18*double_int+36*long_double_int+...
            -60*triple_int+12*quad_int;
    elseif (a <= 4*RR)
        double_int = pi/12*(4*(2*RR)+a).*(2*(2*RR)-a).^2;
        
        F(i) = 8*v(i)-7*4/3*pi*(2*RR).^3+18*double_int;
    else
%         F(i) = v(i)*NN-(NN-1)*4/3*pi*(2*RR).^3;
        F(i) = 8*v(i)-7*4/3*pi*(2*RR).^3;
    end
end
end

function y = myatan(x)
  
   y = atan(x);
   if (y < 0)
       y = y+pi;
   end
end