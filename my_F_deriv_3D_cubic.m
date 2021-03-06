% calculate derivative of free volume in SC lattice
function F = my_F_deriv_3D_cubic(v,RR)
a_vec = v.^(1/3);
da_dv_vec = 1./(3*v.^(2/3));

F = zeros(1,numel(a_vec));
for i=1:numel(a_vec)
    a = a_vec(i);
    
    if (a < 2*RR)
        F(i) = 0;
    elseif (a <= sqrt(6)*RR)
        
        da_dv = da_dv_vec(i);

        long_double_int = pi/12*(4*(2*RR)+sqrt(2)*a).*(2*(2*RR)-sqrt(2)*a).^2;
        dlong_double_int = sqrt(2)*pi/12*da_dv.*((2*(2*RR)-sqrt(2)*a).^2+...
            -(4*(2*RR)+sqrt(2)*a)*2.*(2*(2*RR)-sqrt(2)*a));

        w = sqrt((2*RR)^2*12*a^4-8*a^6);
        q1 = 2*sqrt(2)*a^3;
        
        dw = 1./(2*sqrt((2*RR)^2*12*a^4-8*a^6)).*(4*(2*RR)^2*12*a^3-48*a^5);
        dq1 = 6*sqrt(2)*a^2;
        
        triple_int = w/6-3*a/sqrt(2)*(2*(2*RR)^2-a^2/3)*myatan(2*w/q1)+...
            +4*(2*RR)^3*myatan(sqrt(2)*a*w/(2*RR*q1));
        
        dline1 = dw/6-3/sqrt(2)*(2*(2*RR)^2-a^2/3)*myatan(2*w/q1)+...
            +3*a/sqrt(2)*(2*a/3)*myatan(2*w/q1)+...
            -3*a/sqrt(2)*(2*(2*RR)^2-a^2/3)*(2*dw*q1-2*w*dq1)/q1^2/(1+(2*w/q1)^2);
        dline2 = 4*(2*RR)^3/(1+(sqrt(2)*a*w/(2*RR*q1))^2)*sqrt(2)*((w+a*dw)*(2*RR*q1)+...
            -a*w*2*RR*dq1)/(2*RR*q1)^2;
        dtriple_int = da_dv*(dline1+dline2);
        
        solid_angle = 4*asin(1/3)/(4*pi);
        dihedral_angle = acos(-1/3)/(2*pi);
        
%         F(i) = da_dv*2*(a*sqrt(2))^2+12*dlong_double_int*dihedral_angle+...
%             -4*dtriple_int;
        F(i) = da_dv*4*a^2+12*dlong_double_int*dihedral_angle+...
            -4*dtriple_int;
    elseif (a <= 2*sqrt(2)*RR)
        
        da_dv = da_dv_vec(i);

        double_int = pi/12*(4*(2*RR)+a).*(2*(2*RR)-a).^2;
        ddouble_int = pi/12*da_dv.*((2*(2*RR)-a).^2-(4*(2*RR)+a)*2.*(2*(2*RR)-a));

        long_double_int = pi/12*(4*(2*RR)+sqrt(2)*a).*(2*(2*RR)-sqrt(2)*a).^2;
        dlong_double_int = pi/12*da_dv.*(sqrt(2)*(2*(2*RR)-sqrt(2)*a).^2+...
            -sqrt(2)*(4*(2*RR)+sqrt(2)*a)*2.*(2*(2*RR)-sqrt(2)*a));
        
        w = sqrt((2*RR)^2*4*a^4-2*a^6);
        q1 = 2*a^3;
        
        dw = 1./(2*sqrt((2*RR)^2*4*a^4-2*a^6)).*(4*(2*RR)^2*4*a^3-12*a^5);
        dq1 = 6*a^2;
        
        triple_int = w/6-a*(2*(2*RR)^2-a^2/6)*myatan(2*w/q1)+...
            -a*sqrt(2)/2*(2*(2*RR)^2-a^2/3)*pi/2+...
            +4*(2*RR)^3/3*(myatan(a*w/(2*RR*q1))+pi/2)+...
            +2*(2*RR)^3/3*(2*myatan(a*w/(2*RR*q1)));
        
        dline1 = dw/6-(2*(2*RR)^2-a^2/6)*myatan(2*w/q1)+...
            +a*(2*a/6)*myatan(2*w/q1)+...
            -a*(2*(2*RR)^2-a^2/6)*(2*dw*q1-2*w*dq1)/q1^2/(1+(2*w/q1)^2);
        dline2 = -sqrt(2)/2*(2*(2*RR)^2-a^2/3)*pi/2+...
             +a*sqrt(2)/2*(2*a/3)*pi/2;
        dline3 = 4*(2*RR)^3/3/(1+(a*w/(2*RR*q1))^2)*((w+a*dw)*(2*RR*q1)+...
            -a*w*2*RR*dq1)/(2*RR*q1)^2;
        dline4 = 4*(2*RR)^3/3/(1+(a*w/(2*RR*q1))^2)*((w+a*dw)*(2*RR*q1)+...
            -a*w*2*RR*dq1)/(2*RR*q1)^2;
        dtriple_int = da_dv*(dline1+dline2+dline3+dline4);
        
        dquad_int = 2*dtriple_int-dlong_double_int;
        
        F(i) = 8+18*ddouble_int+36*dlong_double_int+...
            -60*dtriple_int+12*dquad_int;

    elseif (a <= 4*RR)
        double_int = pi/12*(4*(2*RR)+a).*(2*(2*RR)-a).^2;
        
        da_dv = da_dv_vec(i);
        
        ddouble_int = pi/12*da_dv.*((2*(2*RR)-a).^2-(4*(2*RR)+a)*2.*(2*(2*RR)-a));
        
        F(i) = 8+18*ddouble_int;
    else
        F(i) = 8;
    end
end
end

function y = myatan(x)
  
   y = atan(x);
   if (y < 0)
       y = y+pi;
   end
end