% compute free volume in BCC lattice
function F = my_F_3D_bcc(v,RR)

a_vec = (2*v).^(1/3);

F = zeros(1,numel(a_vec));
for i=1:numel(a_vec)
    a = a_vec(i);
    
    if (a < 4/sqrt(3)*RR)
        F(i) = 0;
    elseif (a <= 8/3*RR)
        double_int = pi/12*(4*(2*RR)+a).*(2*(2*RR)-a).^2;
        long_double_int = pi/12*(4*(2*RR)+sqrt(2)*a).*(2*(2*RR)-sqrt(2)*a).^2;
        
        w = sqrt((2*RR)^2*4*a^4-2*a^6);
        q1 = 2*a^3;
        triple_int = w/6-a*(2*(2*RR)^2-a^2/6)*myatan(2*w/q1)+...
            -a*sqrt(2)/2*(2*(2*RR)^2-a^2/3)*pi/2+...
            +4*(2*RR)^3/3*(myatan(a*w/(2*RR*q1))+pi/2)+...
            +2*(2*RR)^3/3*(2*myatan(a*w/(2*RR*q1)));
        
        quad_int = 2*triple_int-long_double_int;
        
         F(i) = a^3-4/3*pi*(2*RR).^3+3*double_int+6*long_double_int+...
             -12*triple_int+3*quad_int;
         
    elseif (a <= 2*sqrt(2)*RR)
        double_int = pi/12*(4*(2*RR)+a).*(2*(2*RR)-a).^2;
        long_double_int = pi/12*(4*(2*RR)+sqrt(2)*a).*(2*(2*RR)-sqrt(2)*a).^2;
        short_double_int = pi/12*(4*(2*RR)+sqrt(3)/2*a).*(2*(2*RR)-sqrt(3)/2*a).^2;
        
        w = sqrt((2*RR)^2*4*a^4-2*a^6);
        q1 = 2*a^3;
        square_triple_int = w/6-a*(2*(2*RR)^2-a^2/6)*myatan(2*w/q1)+...
            -a*sqrt(2)/2*(2*(2*RR)^2-a^2/3)*pi/2+...
            +4*(2*RR)^3/3*(myatan(a*w/(2*RR*q1))+pi/2)+...
            +2*(2*RR)^3/3*(2*myatan(a*w/(2*RR*q1)));
        
        w = sqrt((2*RR)^2*2*a^4-9/16*a^6);
        q1 = 0.5*a^3;
        q2 = sqrt(3)/2*a^3;
        triple_int = w/6-a/2*(2*(2*RR)^2-a^2/6)*myatan(2*w/q1)+...
            -a*sqrt(3)/2*(2*(2*RR)^2-3*a^2/24)*myatan(2*w/q2)+...
            +4*(2*RR)^3/3*(myatan(a*w/(2*RR*q1))+myatan(sqrt(3)/2*a*w/(2*RR*q2)))+...
            +2*(2*RR)^3/3*(2*myatan(sqrt(3)/2*a*w/(2*RR*q2)));
        
        long_triple_int = long_double_int;

        quad_int = square_triple_int;

        square_quad_int = 2*square_triple_int-long_double_int;
        
        quint_int = square_quad_int;

        No = 6;
        Vo = a^3/3;

        solid_angle_vertex = (2*pi/3)/(4*pi);
        solid_angle_base = (pi/3)/(4*pi);
        
        dihedral_angle_base = pi/2/(2*pi);
        dihedral_angle_vertex = acos(-1/2)/(2*pi);
        
% count all intersections
        F(i) = No*Vo-No*(4*solid_angle_base+solid_angle_vertex)*4/3*pi*(2*RR).^3 ...
        +No*(4*double_int*dihedral_angle_base+4*short_double_int*dihedral_angle_vertex+...
             +2*long_double_int-2*triple_int ...
              -4*square_triple_int-2*long_triple_int ...
              +4*quad_int+square_quad_int-quint_int);

    elseif (a <= 8*sqrt(2)/3*RR)
        double_int = pi/12*(4*(2*RR)+a).*(2*(2*RR)-a).^2;
        short_double_int = pi/12*(4*(2*RR)+sqrt(3)/2*a).*(2*(2*RR)-sqrt(3)/2*a).^2;
        
        w = sqrt((2*RR)^2*2*a^4-9/16*a^6);
        q1 = 0.5*a^3;
        q2 = sqrt(3)/2*a^3;
        triple_int = w/6-a/2*(2*(2*RR)^2-a^2/6)*myatan(2*w/q1)+...
            -a*sqrt(3)/2*(2*(2*RR)^2-3*a^2/24)*myatan(2*w/q2)+...
            +4*(2*RR)^3/3*(myatan(a*w/(2*RR*q1))+myatan(sqrt(3)/2*a*w/(2*RR*q2)))+...
            +2*(2*RR)^3/3*(2*myatan(sqrt(3)/2*a*w/(2*RR*q2)));
        
        No = 6;
        Vo = a^3/3;
        
        solid_angle_vertex = (2*pi/3)/(4*pi);
        solid_angle_base = (pi/3)/(4*pi);
        
        dihedral_angle_base = pi/2/(2*pi);
        dihedral_angle_vertex = acos(-1/2)/(2*pi);

        F(i) = No*Vo-No*(4*solid_angle_base+solid_angle_vertex)*4/3*pi*(2*RR).^3 ...
        +No*(4*double_int*dihedral_angle_base+4*short_double_int*dihedral_angle_vertex+...
            -2*triple_int);
    elseif (a <= 4*RR)
        double_int = pi/12*(4*(2*RR)+a).*(2*(2*RR)-a).^2;
        short_double_int = pi/12*(4*(2*RR)+sqrt(3)/2*a).*(2*(2*RR)-sqrt(3)/2*a).^2;
        
        No = 6;
        Vo = a^3/3;
        
        solid_angle_vertex = (2*pi/3)/(4*pi);
        solid_angle_base = (pi/3)/(4*pi);
        
        dihedral_angle_base = pi/2/(2*pi);
        dihedral_angle_vertex = acos(-1/2)/(2*pi);

        F(i) = No*Vo-No*(4*solid_angle_base+solid_angle_vertex)*4/3*pi*(2*RR).^3 ...
        +No*(4*double_int*dihedral_angle_base+4*short_double_int*dihedral_angle_vertex);
    elseif (sqrt(3)/2*a <= 4*RR)
        short_double_int = pi/12*(4*(2*RR)+sqrt(3)/2*a).*(2*(2*RR)-sqrt(3)/2*a).^2;
        
        No = 6;
        Vo = a^3/3;
        
        solid_angle_vertex = (2*pi/3)/(4*pi);
        solid_angle_base = (pi/3)/(4*pi);
        
        dihedral_angle_vertex = acos(-1/2)/(2*pi);

        F(i) = No*Vo-No*(4*solid_angle_base+solid_angle_vertex)*4/3*pi*(2*RR).^3 ...
        +No*(4*short_double_int*dihedral_angle_vertex);
    else
        No = 6;
        Vo = a^3/3;
        
        solid_angle_vertex = (2*pi/3)/(4*pi);
        solid_angle_base = (pi/3)/(4*pi);
        
        F(i) = No*Vo-No*(4*solid_angle_base+solid_angle_vertex)*4/3*pi*(2*RR).^3;
    end
end
end

function y = myatan(x)
  
   y = atan(x);
   if (y < 0)
       y = y+pi;
   end
end