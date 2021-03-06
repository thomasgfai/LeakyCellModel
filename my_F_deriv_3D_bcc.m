% compute derivative of free volume in BCC lattice
function F = my_F_deriv_3D_bcc(v,RR)
a_vec = (2*v).^(1/3);
da_dv_vec = 2^(1/3)/3*1./(v.^(2/3));

F = zeros(1,numel(a_vec));
for i=1:numel(a_vec)
    a = a_vec(i);
    
    if (a < 4/sqrt(3)*RR)
        F(i) = 0;
    elseif (a <= 8/3*RR)
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
        
        F(i) = 3*a^2*da_dv+3*ddouble_int+6*dlong_double_int+...
            -12*dtriple_int+3*dquad_int;
         
    elseif (a <= 2*sqrt(2)*RR)
        da_dv = da_dv_vec(i);

        double_int = pi/12*(4*(2*RR)+a).*(2*(2*RR)-a).^2;
        ddouble_int = pi/12*da_dv.*((2*(2*RR)-a).^2-(4*(2*RR)+a)*2.*(2*(2*RR)-a));
        
        long_double_int = pi/12*(4*(2*RR)+sqrt(2)*a).*(2*(2*RR)-sqrt(2)*a).^2;
        dlong_double_int = pi/12*da_dv*(sqrt(2)*(2*(2*RR)-sqrt(2)*a).^2 ...
            -sqrt(2)*(4*(2*RR)+sqrt(2)*a)*2*(2*(2*RR)-sqrt(2)*a));
        
        short_double_int = pi/12*(4*(2*RR)+sqrt(3)/2*a).*(2*(2*RR)-sqrt(3)/2*a).^2;
        dshort_double_int = pi/12*da_dv*(sqrt(3)/2*(2*(2*RR)-sqrt(3)/2*a).^2 ...
            -sqrt(3)/2*(4*(2*RR)+sqrt(3)/2*a)*2*(2*(2*RR)-sqrt(3)/2*a));

        w = sqrt((2*RR)^2*2*a^4-9/16*a^6);
        q1 = 0.5*a^3;
        q2 = sqrt(3)/2*a^3;
        
        dw = 1./(2*sqrt((2*RR)^2*2*a^4-9/16*a^6)).*(4*(2*RR)^2*2*a^3-27/8*a^5);
        dq1 = 1.5*a^2;
        dq2 = 3*sqrt(3)/2*a^2;
        
        triple_int = w/6-a/2*(2*(2*RR)^2-a^2/6)*myatan(2*w/q1)+...
            -a*sqrt(3)/2*(2*(2*RR)^2-3*a^2/24)*myatan(2*w/q2)+...
            +4*(2*RR)^3/3*(myatan(a*w/(2*RR*q1))+myatan(sqrt(3)/2*a*w/(2*RR*q2)))+...
            +2*(2*RR)^3/3*(2*myatan(sqrt(3)/2*a*w/(2*RR*q2)));
        
        dline1 = dw/6-1/2*(2*(2*RR)^2-a^2/6)*myatan(2*w/q1)+...
            +a/2*(a/3)*myatan(2*w/q1)+...
            -a/2*(2*(2*RR)^2-a^2/6)*(2*dw*q1-2*w*dq1)/q1^2/(1+(2*w/q1)^2);
        dline2 = -sqrt(3)/2*(2*(2*RR)^2-3*a^2/24)*myatan(2*w/q2)+...
             +a*sqrt(3)/2*(a/4)*myatan(2*w/q2)+...
             -a*sqrt(3)/2*(2*(2*RR)^2-3*a^2/24)*(2*dw*q2-2*w*dq2)/q2^2/(1+(2*w/q2)^2);
        dline3 = 4*(2*RR)^3/3*(1/(1+(a*w/(2*RR*q1))^2)*((w+a*dw)*(2*RR*q1)-a*w*2*RR*dq1)/(2*RR*q1)^2+...
            1/(1+(sqrt(3)/2*a*w/(2*RR*q2))^2)*((sqrt(3)/2*(w+a*dw))*(2*RR*q2)-sqrt(3)/2*a*w*2*RR*dq2)/(2*RR*q2)^2);
        dline4 = 4*(2*RR)^3/3/(1+(sqrt(3)/2*a*w/(2*RR*q2))^2)*((sqrt(3)/2*(w+a*dw))*(2*RR*q2)+...
            -sqrt(3)/2*a*w*2*RR*dq2)/(2*RR*q2)^2;
        dtriple_int = da_dv*(dline1+dline2+dline3+dline4);

        w = sqrt((2*RR)^2*4*a^4-2*a^6);
        q1 = 2*a^3;
        
        dw = 1./sqrt((2*RR)^2*4*a^4-2*a^6).*(4*(2*RR)^2*4*a^3-18*a^5);
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
        
        long_triple_int = long_double_int;
        dlong_triple_int = dlong_double_int;

        quad_int = square_triple_int;
        dquad_int = dsquare_triple_int;

        square_quad_int = 2*square_triple_int-long_double_int;
        dsquare_quad_int = 2*dsquare_triple_int-dlong_double_int;

        quint_int = square_quad_int;
        dquint_int = dsquare_quad_int;
        
        No = 6;
        Vo = a^3/3;
        
        solid_angle_vertex = (2*pi/3)/(4*pi);
        solid_angle_base = (pi/3)/(4*pi);
        
        dihedral_angle_base = pi/2/(2*pi);
        dihedral_angle_vertex = acos(-1/2)/(2*pi);
        
        F(i) = No*a^2*da_dv ...
             +No*(4*ddouble_int*dihedral_angle_base+4*dshort_double_int*dihedral_angle_vertex+...
             +2*dlong_double_int-2*dtriple_int ...
              -4*dsquare_triple_int-2*dlong_triple_int ...
              +4*dquad_int+dsquare_quad_int-dquint_int);
        
    elseif (a <= 8*sqrt(2)/3*RR)
        da_dv = da_dv_vec(i);

        double_int = pi/12*(4*(2*RR)+a).*(2*(2*RR)-a).^2;
        ddouble_int = pi/12*da_dv.*((2*(2*RR)-a).^2-(4*(2*RR)+a)*2.*(2*(2*RR)-a));

        short_double_int = pi/12*(4*(2*RR)+sqrt(3)/2*a).*(2*(2*RR)-sqrt(3)/2*a).^2;
        dshort_double_int = pi/12*da_dv*(sqrt(3)/2*(2*(2*RR)-sqrt(3)/2*a).^2 ...
            -sqrt(3)/2*(4*(2*RR)+sqrt(3)/2*a)*2*(2*(2*RR)-sqrt(3)/2*a));

        w = sqrt((2*RR)^2*2*a^4-9/16*a^6);
        q1 = 0.5*a^3;
        q2 = sqrt(3)/2*a^3;
        
        dw = 1./(2*sqrt((2*RR)^2*2*a^4-9/16*a^6)).*(4*(2*RR)^2*2*a^3-27/8*a^5);
        dq1 = 1.5*a^2;
        dq2 = 3*sqrt(3)/2*a^2;
        
        triple_int = w/6-a/2*(2*(2*RR)^2-a^2/6)*myatan(2*w/q1)+...
            -a*sqrt(3)/2*(2*(2*RR)^2-3*a^2/24)*myatan(2*w/q2)+...
            +4*(2*RR)^3/3*(myatan(a*w/(2*RR*q1))+myatan(sqrt(3)/2*a*w/(2*RR*q2)))+...
            +2*(2*RR)^3/3*(2*myatan(sqrt(3)/2*a*w/(2*RR*q2)));
        
        dline1 = dw/6-1/2*(2*(2*RR)^2-a^2/6)*myatan(2*w/q1)+...
            +a/2*(a/3)*myatan(2*w/q1)+...
            -a/2*(2*(2*RR)^2-a^2/6)*(2*dw*q1-2*w*dq1)/q1^2/(1+(2*w/q1)^2);
        dline2 = -sqrt(3)/2*(2*(2*RR)^2-3*a^2/24)*myatan(2*w/q2)+...
             +a*sqrt(3)/2*(a/4)*myatan(2*w/q2)+...
             -a*sqrt(3)/2*(2*(2*RR)^2-3*a^2/24)*(2*dw*q2-2*w*dq2)/q2^2/(1+(2*w/q2)^2);
        dline3 = 4*(2*RR)^3/3*(1/(1+(a*w/(2*RR*q1))^2)*((w+a*dw)*(2*RR*q1)-a*w*2*RR*dq1)/(2*RR*q1)^2+...
            1/(1+(sqrt(3)/2*a*w/(2*RR*q2))^2)*((sqrt(3)/2*(w+a*dw))*(2*RR*q2)-sqrt(3)/2*a*w*2*RR*dq2)/(2*RR*q2)^2);
        dline4 = 4*(2*RR)^3/3/(1+(sqrt(3)/2*a*w/(2*RR*q2))^2)*((sqrt(3)/2*(w+a*dw))*(2*RR*q2)+...
            -sqrt(3)/2*a*w*2*RR*dq2)/(2*RR*q2)^2;
        dtriple_int = da_dv*(dline1+dline2+dline3+dline4);
        
        No = 6;
        Vo = a^3/3;
        
        solid_angle_vertex = (2*pi/3)/(4*pi);
        solid_angle_base = (pi/3)/(4*pi);
        
        dihedral_angle_base = pi/2/(2*pi);
        dihedral_angle_vertex = acos(-1/2)/(2*pi);

        F(i) = No*a^2*da_dv+No*(4*ddouble_int*dihedral_angle_base+4*dshort_double_int*dihedral_angle_vertex+...
            -2*dtriple_int);
    elseif (a <= 4*RR)
        da_dv = da_dv_vec(i);

        double_int = pi/12*(4*(2*RR)+a).*(2*(2*RR)-a).^2;
        ddouble_int = pi/12*da_dv.*((2*(2*RR)-a).^2-(4*(2*RR)+a)*2.*(2*(2*RR)-a));

        short_double_int = pi/12*(4*(2*RR)+sqrt(3)/2*a).*(2*(2*RR)-sqrt(3)/2*a).^2;
        dshort_double_int = pi/12*da_dv*(sqrt(3)/2*(2*(2*RR)-sqrt(3)/2*a).^2 ...
            -sqrt(3)/2*(4*(2*RR)+sqrt(3)/2*a)*2*(2*(2*RR)-sqrt(3)/2*a));

        No = 6;
        Vo = a^3/3;
        
        solid_angle_vertex = (2*pi/3)/(4*pi);
        solid_angle_base = (pi/3)/(4*pi);
        
        dihedral_angle_base = pi/2/(2*pi);
        dihedral_angle_vertex = acos(-1/2)/(2*pi);

        F(i) = No*a^2*da_dv+No*(4*ddouble_int*dihedral_angle_base+4*dshort_double_int*dihedral_angle_vertex);
    elseif (sqrt(3)/2*a <= 4*RR)
        da_dv = da_dv_vec(i);

        short_double_int = pi/12*(4*(2*RR)+sqrt(3)/2*a).*(2*(2*RR)-sqrt(3)/2*a).^2;
        dshort_double_int = pi/12*da_dv*(sqrt(3)/2*(2*(2*RR)-sqrt(3)/2*a).^2 ...
            -sqrt(3)/2*(4*(2*RR)+sqrt(3)/2*a)*2*(2*(2*RR)-sqrt(3)/2*a));
    
        No = 6;
        Vo = a^3/3;
        
        solid_angle_vertex = (2*pi/3)/(4*pi);
        solid_angle_base = (pi/3)/(4*pi);
        
        dihedral_angle_vertex = acos(-1/2)/(2*pi);

        F(i) = No*a^2*da_dv+No*(4*dshort_double_int*dihedral_angle_vertex);
    else
        da_dv = da_dv_vec(i);
                
        No = 6;
        Vo = a^3/3;
        
        F(i) = No*a^2*da_dv;
    end
end
end

function y = myatan(x)
  
   y = atan(x);
   if (y < 0)
       y = y+pi;
   end
end