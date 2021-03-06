% calculate derivative of free volume in FCC lattice
function F = my_F_deriv_3D_fcc(v,RR)
a_vec = (sqrt(2)*v).^(1/3);
da_dv_vec = sqrt(2)./(3*a_vec.^2);

F = zeros(1,numel(a_vec));
for i=1:numel(a_vec)
    a = a_vec(i);
    
    if (a < 2*RR)
        F(i) = 0;
    elseif (a <= 2*sqrt(2)*RR)
        da_dv = da_dv_vec(i);
        
        Nt = 8; No = 6;
        Vt = a.^3/(6*sqrt(2)); Vo = sqrt(2)*a.^3/3;
        
        solid_angle_t = acos(23/27)/(4*pi);
        solid_angle_o = 4*asin(1/3)/(4*pi);

        dihedral_angle_t = acos(1/3)/(2*pi);
        dihedral_angle_o = acos(-1/3)/(2*pi);

        double_int = pi/12*(4*(2*RR)+a).*(2*(2*RR)-a).^2;
        ddouble_int = pi/12*da_dv.*((2*(2*RR)-a).^2-(4*(2*RR)+a)*2.*(2*(2*RR)-a));

        long_double_int = pi/12*(4*(2*RR)+sqrt(2)*a).*(2*(2*RR)-sqrt(2)*a).^2;
        dlong_double_int = sqrt(2)*pi/12*da_dv.*((2*(2*RR)-sqrt(2)*a).^2+...
            -(4*(2*RR)+sqrt(2)*a)*2.*(2*(2*RR)-sqrt(2)*a));
        
        w = sqrt((2*RR)^2*3*a^4-a^6);
        q1 = a^3;
         
        dw = 1./(2*sqrt((2*RR)^2*3*a^4-a^6)).*(4*(2*RR)^2*3*a^3-6*a^5);
        dq1 = 3*a^2;
        
        triple_int = w/6-3*a/2*(2*(2*RR)^2-a^2/6)*myatan(2*w/q1)+...
            +4*(2*RR)^3*myatan(a*w/(2*RR*q1));
        
        dline1 = dw/6-3/2*(2*(2*RR)^2-a^2/6)*myatan(2*w/q1)+...
            +3*a/2*(a/3)*myatan(2*w/q1)+...
            -3*a/2*(2*(2*RR)^2-a^2/6)*(2*dw*q1-2*w*dq1)/q1^2/(1+(2*w/q1)^2);
        dline2 = 4*(2*RR)^3/(1+(a*w/(2*RR*q1))^2)*((w+a*dw)*(2*RR*q1)+...
            -a*w*2*RR*dq1)/(2*RR*q1)^2;
        dtriple_int = da_dv*(dline1+dline2);
        
        w = sqrt((2*RR)^2*4*a^4-2*a^6);
        q1 = 2*a^3;
        
        dw = 1./(2*sqrt((2*RR)^2*4*a^4-2*a^6)).*(4*(2*RR)^2*4*a^3-12*a^5);
        dq1 = 6*a^2;
        
        square_triple_int = w/6-a*(2*(2*RR)^2-a^2/6)*myatan(2*w/q1)+...
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
        dsquare_triple_int = da_dv*(dline1+dline2+dline3+dline4);
        
        quad_int = 2*square_triple_int-long_double_int;
        dquad_int = 2*dsquare_triple_int-dlong_double_int;

        F(i) = Nt*(da_dv*a.^2/(2*sqrt(2))...
            +3*ddouble_int*dihedral_angle_t)...
            +No*(da_dv*a.^2/sqrt(2)...
            +2*ddouble_int*dihedral_angle_o...
            +dlong_double_int-2*dsquare_triple_int+0.5*dquad_int)...
            -4*dtriple_int;
        
    elseif (a <= 2*sqrt(3)*RR)
        da_dv = da_dv_vec(i);

        Nt = 8; No = 6;
        Vt = a.^3/(6*sqrt(2)); Vo = sqrt(2)*a.^3/3;
        
        solid_angle_t = acos(23/27)/(4*pi);
        solid_angle_o = 4*asin(1/3)/(4*pi);

        dihedral_angle_t = acos(1/3)/(2*pi);
        dihedral_angle_o = acos(-1/3)/(2*pi);

        double_int = pi/12*(4*(2*RR)+a).*(2*(2*RR)-a).^2;
        ddouble_int = pi/12*da_dv.*((2*(2*RR)-a).^2-(4*(2*RR)+a)*2.*(2*(2*RR)-a));

        w = sqrt((2*RR)^2*3*a^4-a^6);
        q1 = a^3;
         
        dw = 1./(2*sqrt((2*RR)^2*3*a^4-a^6)).*(4*(2*RR)^2*3*a^3-6*a^5);
        dq1 = 3*a^2;
        
        triple_int = w/6-3*a/2*(2*(2*RR)^2-a^2/6)*myatan(2*w/q1)+...
            +4*(2*RR)^3*myatan(a*w/(2*RR*q1));
        
        dline1 = dw/6-3/2*(2*(2*RR)^2-a^2/6)*myatan(2*w/q1)+...
            +3*a/2*(a/3)*myatan(2*w/q1)+...
            -3*a/2*(2*(2*RR)^2-a^2/6)*(2*dw*q1-2*w*dq1)/q1^2/(1+(2*w/q1)^2);
        dline2 = 4*(2*RR)^3/(1+(a*w/(2*RR*q1))^2)*((w+a*dw)*(2*RR*q1)+...
            -a*w*2*RR*dq1)/(2*RR*q1)^2;
        dtriple_int = da_dv*(dline1+dline2);
        
        F(i) = Nt*(da_dv*a.^2/(2*sqrt(2)) ...
            +3*ddouble_int*dihedral_angle_t)+No*(da_dv*sqrt(2)*a.^2 ...
            +8*ddouble_int*dihedral_angle_o)-16*dtriple_int;
        
    elseif (a <= 4*RR)
        da_dv = da_dv_vec(i);
        Nt = 8; No = 6;

        dihedral_angle_t = acos(1/3)/(2*pi);
        dihedral_angle_o = acos(-1/3)/(2*pi);
        
        double_int = pi/12*(4*(2*RR)+a).*(2*(2*RR)-a).^2;
        ddouble_int = pi/12*da_dv.*((2*(2*RR)-a).^2-(4*(2*RR)+a)*2.*(2*(2*RR)-a));

        F(i) = Nt*(da_dv*a.^2/(2*sqrt(2)) ...
            +3*ddouble_int*dihedral_angle_t) ...
            +No*(da_dv*sqrt(2)*a.^2 ...
            +8*ddouble_int*dihedral_angle_o);
    else
        da_dv = da_dv_vec(i);
        Nt = 8; No = 6;
                
        F(i) = Nt*da_dv*a.^2/(2*sqrt(2)) ...
            +No*da_dv*sqrt(2)*a.^2;
    end
end
end

function y = myatan(x)
  
   y = atan(x);
   if (y < 0)
       y = y+pi;
   end
end