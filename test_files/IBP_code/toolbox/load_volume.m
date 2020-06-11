function h = load_volume(name, n, options)

% load_volume - load a 3-D data
%
%   h = load_volume(name, n, options)
%
%   Copyright (c) 2015 Gabriel Peyre

options.null = 0;

% helpers
t = linspace(0,1,n);
[X,Y,Z] = ndgrid(t,t,t);
sphere = @(c,r)( (X-c(1)).^2 + (Y-c(2)).^2 + (Z-c(2)).^2 )<= r^2;
cube = @(c,r)  max( max(abs(X-c(1)),abs(Y-c(2))), abs(Z-c(2)) )<= r;

switch name
    case 'sphere'
        r = .35;
        c = [1 1 1]/2;
        h = sphere(c,r);        
    case 'spheres'
        r = .25;
        c1 = [1 1 1]*.3;
        c2 = [1 1 1]*.7;
        h = sphere(c1,r) + sphere(c2,r);
        
    case 'boxes'
        r = .15;
        c1 = [.25 1/2 1/2];
        c2 = [.75 1/2 1/2];
        h = cube(c1,r) + cube(c2,r);

    case 'cone'
        r = getoptions(options, 'radius', .8);
        c = getoptions(options, 'center', [0 0 0]);
        x = linspace(-1,1,n);
        [X,Y,Z] = ndgrid(x,x,x);

        ax = getoptions(options, 'axis', [1 1 1]);
        A = rand(3); A(:,1) = ax;
        [A,R] = qr(A); 
        A = A * sign( sum(A(:,1).*ax(:)) );            
        pos = [X(:),Y(:),Z(:)]';
        pos = A'*pos;
        X = reshape(pos(1,:),n,n,n);
        Y = reshape(pos(2,:),n,n,n);
        Z = reshape(pos(3,:),n,n,n);
        D = sqrt(Y.^2+Z.^2);
        eta = getoptions(options, 'eta', .95);
        h = (X>=0) .* ( D <= r*(1-X/eta));
        
    case 'cone_rotated'
        r = getoptions(options, 'radius', .8);
        c = getoptions(options, 'center', [0 0 0]);
        x = linspace(-1,1,n);
        [X,Y,Z] = ndgrid(x,x,x);

        ax = getoptions(options, 'axis', [1 1 1]);
        A = rand(3); A(:,1) = ax;
        [A,R] = qr(A); 
        A = A * sign( sum(A(:,1).*ax(:)) );            
        pos = [X(:),Y(:),Z(:)]';
        pos = A'*pos;
        X = reshape(pos(1,:),n,n,n);
        Y = reshape(pos(2,:),n,n,n);
        Z = reshape(pos(3,:),n,n,n);
        D = sqrt(Y.^2+Z.^2);
        eta = getoptions(options, 'eta', .95);
        h = (X>=0) .* ( D <= r*(1-X/eta));
        
    case 'spiky'
        options.radius = .3;
        cones_layout = getoptions(options, 'cones_layout', 'ico');
        vertex = compute_base_mesh(cones_layout);
        vertex = vertex';
        h = zeros(n,n,n);
        for i=1:size(vertex,2)
            % progressbar(i, size(vertex,2));
            options.axis = vertex(:,i);
            h = h | load_volume('cone', n, options);
        end
    otherwise
        % error('Unknown 3D shape');
        h = [];
end
h = double(h);

end



function [vertex,face] = compute_base_mesh(type, j, options)

% compute_base_mesh - generate a simple triangulation.
%
%   [vertex,face] = compute_base_mesh(type,j, options);
%
%   'type' can be one 'triangle', 'square', 'square1', 'L', 'L1', 'tetra',
%   'oct', or 'ico' or 'rand'.
%
%	j is the (optional) number of subdivision levels
%
%   Copyright (c) 2004 Gabriel Peyre

if nargin<2
    j = 0;
end

switch  lower(type)
    case 'rand'
        nverts_base = getoptions(options, 'nverts_base', 30);
        vertex = rand(3,nverts_base);
        vertex(3,:) = 0;
        face = delaunay(vertex(1,:),vertex(2,:));
        
    case 'triangle'
        vertex = [0,0,0; 1,0,0; 0.5,1/sqrt(2),0]';
        vertex = vertex - repmat(mean(vertex,2), [1 3]);
        face = [1,2,3]';
        
    case 'square'
        nbr = getoptions(options, 'n', 2^j);
        x = 1:nbr;
        [Y,X] = meshgrid(x,x);
        I = X+(Y-1)*nbr;
        a = I(1:end-1, 1:end-1);
        b = I(2:end, 1:end-1);
        c = I(1:end-1, 2:end);
        face = cat(1,a(:)',b(:)',c(:)');
        a = I(1:end-1, 2:end);
        b = I(2:end, 1:end-1);
        c = I(2:end, 2:end);
        face= [face, cat(1,a(:)',b(:)',c(:)')];
        x = linspace(0,1,nbr);
        [Y,X] = meshgrid(x,x);
        vertex = cat(1, X(:)', Y(:)');
        return;
    case 'square1'
        vertex = [0,0,0; 1,0,0; 1,1,0; 0,1,0];
        face = [1,2,3; 3,4,1];
    case 'l'
        vertex = [  0,0,0; 1,0,0; 2,0,0;
                    0,1,0; 1,1,0; 2,1,0;
                    0,2,0; 1,2,0;
                    0,3,0; 1,3,0];
        face = [1,2,5; 1,5,4;
                2,3,6; 2,6,5;
                4,5,8; 4,8,7;
                7,8,10; 7,10,9];
    case 'l1'
        vertex = [  0,0,0; 1,0,0; 2,0,0;
                    0,1,0; 1,1,0; 2,1,0;
                    0,2,0; 1,2,0;
                    0.5,0.5,0; 1.5,0.5,0; 0.5,1.5,0];
        face = [1,2,9; 2,5,9; 5,4,9; 4,1,9; 
                2,3,10; 3,6,10; 6,5,10; 5,2,10;
                4,5,11; 5,8,11; 8,7,11; 7,4,11];  
    case 'tetra',
        sqrt_3 = 0.5773502692;
        vertex = [  sqrt_3,  sqrt_3,  sqrt_3 ;
                -sqrt_3, -sqrt_3,  sqrt_3 ;
                -sqrt_3,  sqrt_3, -sqrt_3 ;
                 sqrt_3, -sqrt_3, -sqrt_3 ]; 
        face = [ 1, 2, 3;
                1, 4, 2;
                3, 2, 4;
                4, 1, 3 ]; 
    case 'oct',
        vertex = [  1,  0,  0 ;
              -1,  0,  0 ;
               0,  1,  0 ;
               0, -1,  0 ;
               0,  0,  1 ;
               0,  0, -1 ];
        face = [ 1 5 3 ;
              3 5 2 ;
              2 5 4 ;
              4 5 1 ;
              1 3 6 ;
              3 2 6 ;
              2 4 6 ;
              4 1 6 ];    
    case 'ico',
        tau = 0.8506508084;
        one = 0.5257311121;
        vertex = [  tau,  one,    0;
                -tau,  one,    0
                -tau, -one,    0;
                tau, -one,    0;
                one,   0 ,  tau;
                one,   0 , -tau;
                -one,   0 , -tau;
                -one,   0 ,  tau;
                0 ,  tau,  one;
                0 , -tau,  one;
                0 , -tau, -one;
                0 ,  tau, -one ];
        face = [  5,  9,  8 ;
               5,  8, 10 ;
               6,  7, 12 ;
               6, 11,  7 ;
               1,  5,  4 ;
               1,  4,  6 ;
               3,  8,  2 ;
               3,  2,  7 ;
               9,  1, 12 ;
               9, 12,  2 ;
              10, 11,  4 ;
              10,  3, 11 ;
               9,  5,  1 ;
              12,  1,  6 ;
               5, 10,  4 ;
               6,  4, 11 ;
               8,  9,  2 ;
               7,  2, 12 ;
               8,  3, 10 ;
               7, 11,  3 ];
    case 'cube',
        vertex = [  -1 1 1;
                    1 1 1;
                    1 -1 1;
                    -1 -1 1;
                    -1 1 -1;
                    1 1 -1;
                    1 -1 -1;
                    -1 -1 -1 ];
        face = [1 2 3; 3 4 1;
                2 6 7; 7 3 2;
                5 8 7; 7 6 5;
                5 1 4; 4 8 5;
                5 6 2; 2 1 5;
                3 4 8; 8 7 3 ];
            
            
        vertex = [ -1 -1 -1;
                    1 -1 -1;
                    1  1 -1;
                   -1  1 -1;
                    1 -1  1;
                   -1 -1  1;
                   -1  1  1;
                    1  1  1;];
        face = [1 2 3; 3 4 1;
                3 2 5; 5 8 3;
                8 5 7; 7 5 6;
                7 6 1; 4 7 1;
                7 4 3; 8 7 3;
                1 6 5; 5 2 1];
            
    case 'sphere',
        [vertex,face] = compute_base_mesh('ico');
        % subdivision
        [vertex,face] = subdivide_sphere(vertex,face,j);
		[vertex,face] = check_face_vertex(vertex,face);
        return;
    otherwise
        error('Unknown type of mesh.');
end

end

% subdivision
% [vertex,face] = perform_mesh_subdivision(vertex,face,j);
% [vertex,face] = check_face_vertex(vertex,face);
