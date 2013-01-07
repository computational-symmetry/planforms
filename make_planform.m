function pf = make_planform( varargin )
% MAKE_PLANFORM creates a planform, a sum of 4 or 6 sine wave gratings.
%   pf  :    A data structure for the planform with fields
%            
%            img_pix:           image size in pixels
%            X:                 array of X values [img_pix, img_pix]
%            Y:                 array of Y values [img_pix, img_pix]
%            cyc_per_img:       sine wave cycles per image
%            cyc_pix:           cycles per pixel
%            gaussian_space_constant: width (pix) of Gaussian
%            gaussian_mask:     image mask [img_pix, img_pix]
%            gray_scale:        Grayscale value
%            ampl:              sine wave amplitude
%            base_angle_offset: orientation relative to horizontal
%            angle_rad:         lattice angle
%            phase_rad:         for 4-component planforms, sq vs. supsq
%            class_type:        number of components
%            img:               data structure with img maps
%                   C1...C6:    Grating components
%                   P12, P34, P56:  Sums of pairs
%                   P1234, P123456: Sums
%
%   Takes a planform structure or no input.
%
%   Use cases:
%
%   >> pf1 = make_planform();
%   Assigning default value of img_pix = 600
%   Assigning default value of cyc_per_img = 12
%   Assigning default value of gaussian_space_constant = 1800
%   Assigning default value of gray_scale = 255
%   Assigning default value of ampl = 1
%   Assigning default value of base_angle_offset = 0
%   Assigning default value of angle_rad = 3.217506e-01
%   Assigning default value of phase_rad = 0
%   Assigning default value of class_type = 4
%   >> imagesc( pf1.img.P1234 ); colormap( gray );
%
%   >> pf2 = pf1;
%   >> pf2.phase_rad = pi; % make super-square
%   >> pf2 = make_planform( pf2 );
%   >> imagesc( pf2.img.P1234 );

%--------------------------------------------------------------------------
% History
% 
% 130107 rog adapted from symmetry_sim.m, cleaned up code, documentation
%


%--------------------------------------------------------------------------
% Known bugs, problems, feature wish list
%
% 130107 need to add citation to planforms paper.

%--------------------------------------------------------------------------


%----   Check parameter values

if nargin > 1
    error('Too many input arguments.');
end

if nargin == 0
    pf = check_planform_params();
else
    pf = varargin{1};
    pf = check_planform_params( pf );    
end

%----   Generate based on class type
switch pf.class_type
    case 4       
        pair_angle = pi/2; % C1 and C2, C3 and C4 are orthogonal 
        
        % C1 and C2 are rotated with respect to the origin by angle_rad
        % C3 and C4 are rotated -angle_rad
        comp_axis  = [ 0 pair_angle pair_angle 0 ];
        comp_offset = [ pf.angle_rad pf.angle_rad -pf.angle_rad -pf.angle_rad ];
       
        % Generate C1...C4
        % In theory, could generate C1 and then rotate the whole image
        C1 = grating( pf.X, pf.Y, pf.ampl, 0, pf.base_angle_offset + comp_axis(1) + comp_offset(1), pf.cyc_pix );
        C2 = grating( pf.X, pf.Y, pf.ampl, 0, pf.base_angle_offset + comp_axis(2) + comp_offset(2), pf.cyc_pix );
        
        C3 = grating( pf.X, pf.Y, pf.ampl, pf.phase_rad, pf.base_angle_offset + comp_axis(3) + comp_offset(3), pf.cyc_pix );
        C4 = grating( pf.X, pf.Y, pf.ampl, pf.phase_rad, pf.base_angle_offset + comp_axis(4) + comp_offset(4), pf.cyc_pix );
        
        % Create grayscale images as sums of components, gray-scaled
        
        P1234=((C1+C2+C3+C4).*pf.gaussian_mask)/4;
        P1234=(P1234*pf.gray_scale/2) + ones( size(P1234) )*pf.gray_scale/2;
        img.P1234 = P1234;

        P12 = (C1 + C2).*pf.gaussian_mask/2;
        P12=(P12*pf.gray_scale/2) + ones( size(P12) )*pf.gray_scale/2;
        img.P12 = P12;
        
        P34 = (C3 + C4).*pf.gaussian_mask/2;
        P34=(P34*pf.gray_scale/2) + ones( size(P34) )*pf.gray_scale/2;
        img.P34 = P34;
                
    case 6
        % Origins of 3, 6 component patterns are offset by 2*pi/3
        pair_angle = 2*pi/3;
        comp_axis  = [ 0 pair_angle pair_angle 0 2*pair_angle 2*pair_angle ];
        comp_offset = [ pf.angle_rad pf.angle_rad -pf.angle_rad -pf.angle_rad pf.angle_rad -pf.angle_rad];
                
        % phase_rad argument doesn't make sense with n=3, 6?
        C1 = grating( pf.X, pf.Y, pf.ampl, 0, pf.base_angle_offset + comp_axis(1) + comp_offset(1), pf.cyc_pix );
        C2 = grating( pf.X, pf.Y, pf.ampl, 0, pf.base_angle_offset + comp_axis(2) + comp_offset(2), pf.cyc_pix );
        
        C3 = grating( pf.X, pf.Y, pf.ampl, 0, pf.base_angle_offset + comp_axis(3) + comp_offset(3), pf.cyc_pix );
        C4 = grating( pf.X, pf.Y, pf.ampl, 0, pf.base_angle_offset + comp_axis(4) + comp_offset(4), pf.cyc_pix );
        
        C5 = grating( pf.X, pf.Y, pf.ampl, 0, pf.base_angle_offset + comp_axis(5) + comp_offset(5), pf.cyc_pix ); 
        C6 = grating( pf.X, pf.Y, pf.ampl, 0, pf.base_angle_offset + comp_axis(6) + comp_offset(6), pf.cyc_pix );      
    
         % Create grayscale images as sums of components, gray-scaled
        
        P123456 = ((C1+C2+C3+C4+C5+C6).*pf.gaussian_mask)/6;
        P123456=(P123456*pf.gray_scale/2) + ones( size(P123456) )*pf.gray_scale/2;
        img.P123456 = P123456;
        
        P12 = (C1 + C2).*pf.gaussian_mask/2;
        %scale=max(max(P12));
        P12=(P12*pf.gray_scale/2) + ones( size(P12) )*pf.gray_scale/2;
        img.P12 = P12;
        
        P34 = (C3 + C4).*pf.gaussian_mask/2;
        %scale=max(max(P34));
        P34=(P34*pf.gray_scale/2) + ones( size(P34) )*pf.gray_scale/2;
        img.P34 = P34;
        
        P56 = (C5 + C6).*pf.gaussian_mask/2;
        %scale=max(max(P56));
        P56=(P56*pf.gray_scale/2) + ones( size(P56) )*pf.gray_scale/2;
        img.P56 = P56;
        
        C5=C3.*pf.gaussian_mask;
        C5=(C5*pf.gray_scale/2) + ones( size(C5) )*pf.gray_scale/2;
        img.C5 = C5;
        
        C6=C4.*pf.gaussian_mask;
        C6=(C6*pf.gray_scale/2) + ones( size(C6) )*pf.gray_scale/2;
        img.C6 = C6;
    
    otherwise
        error('Planform class_type not allowed.');
end

%----   C1...C4 in common across both class 4, 6 so can generate once.
C1=C1.*pf.gaussian_mask;
C1=(C1*pf.gray_scale/2) + ones( size(C1) )*pf.gray_scale/2;
img.C1 = C1;

C2=C2.*pf.gaussian_mask;
C2=(C2*pf.gray_scale/2) + ones( size(C2) )*pf.gray_scale/2;
img.C2 = C2;

C3=C3.*pf.gaussian_mask;
C3=(C3*pf.gray_scale/2) + ones( size(C3) )*pf.gray_scale/2;
img.C3 = C3;

C4=C4.*pf.gaussian_mask;
C4=(C4*pf.gray_scale/2) + ones( size(C4) )*pf.gray_scale/2;
img.C4 = C4;

%----   Add img to planform struct
pf.img = img;

return;
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function pf = check_planform_params( varargin )

%----   Create data structure with default values for comparison
d = planform_assign_defaults(); 

if nargin == 0
    pf = [];
end

if nargin == 1
    arg = varargin{1};
    if ~isstruct( arg ) % input not a structure
        error('Input not a data structure.\n');
    end
    pf = arg;
end

%----   Pixels in image (img_pix x img_pix)
if ~isfield( pf, 'img_pix' )
    pf.img_pix = d.img_pix;
    fprintf('Assigning default value of img_pix = %d\n', pf.img_pix);
end

%----   Make X, Y grid
[ pf.X, pf.Y ] = meshgrid( 0:pf.img_pix-1 );

%----   Cycles per image
if ~isfield( pf, 'cyc_per_img' )
    pf.cyc_per_img = d.cyc_per_img;
    fprintf('Assigning default value of cyc_per_img = %d\n', pf.cyc_per_img);
end

%----   Calculate cyc/pix
pf.cyc_pix = pf.cyc_per_img/pf.img_pix;

%----   Gaussian space constant
if ~isfield( pf, 'gaussian_space_constant' )
    pf.gaussian_space_constant = d.gaussian_space_constant;
    fprintf('Assigning default value of gaussian_space_constant = %d\n', pf.gaussian_space_constant);
end

%----   Calculate Gaussian mask
pf.gaussian_mask = make_gaussian_mask( pf );

%----   Gray scale maximum value 
if ~isfield( pf, 'gray_scale' )
    pf.gray_scale = d.gray_scale;
    fprintf('Assigning default value of gray_scale = %d\n', pf.gray_scale);
end

%----   Amplitude
if ~isfield( pf, 'ampl' )
    pf.ampl = d.ampl;
    fprintf('Assigning default value of ampl = %d\n', pf.ampl);
end

%----   Base angle offset from horizontal
if ~isfield( pf, 'base_angle_offset' )
    pf.base_angle_offset = d.base_angle_offset; % Should read from some defaults file
    fprintf('Assigning default value of base_angle_offset = %d\n', pf.base_angle_offset);
end

%----   Angle 
if ~isfield( pf, 'angle_rad' )
    pf.angle_rad = d.angle_rad; % Should read from some defaults file
    fprintf('Assigning default value of angle_rad = %d\n', pf.angle_rad);
end

%----   Phase of 2nd checkerboard pair (4 component planforms only)
if ~isfield( pf, 'phase_rad' )
    pf.phase_rad = d.phase_rad;
    fprintf('Assigning default value of phase_rad = %d\n', pf.phase_rad);
end

%----   Class type (n components)
if ~isfield( pf, 'class_type' )
    pf.class_type = d.class_type;
    fprintf('Assigning default value of class_type = %d\n', pf.class_type);
end


return
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
function defaults = planform_assign_defaults()

defaults.img_pix                 = 600;
defaults.cyc_per_img             = 12;
defaults.gaussian_space_constant = 3 * defaults.img_pix; % No edge blurring
defaults.ampl                    = 1;
defaults.base_angle_offset       = 0;
defaults.angle_rad               = atan2(1,3); % 3,1 lattice
defaults.gray_scale              = 255;
defaults.phase_rad               = 0; % Square; pi for super square with 3:1 lattice
defaults.class_type              = 4; % 4 components
defaults.pair_angle_offset       = 0;

return
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
function g = make_gaussian_mask( pf )

[x, y] = meshgrid( linspace(-pf.img_pix/2, pf.img_pix/2, pf.img_pix), linspace(-pf.img_pix/2, pf.img_pix/2, pf.img_pix) );

g = exp(-((x .^ 2) + (y .^ 2)) / (pf.gaussian_space_constant ^ 2));
return
%--------------------------------------------------------------------------




%--------------------------------------------------------------------------
function cc = grating( X, Y, ampl, phase_rad, angle_rad, cyc_pix )
% Generates gray scale grating

f  = cyc_pix*2*pi;
aa = cos( angle_rad )*f;
bb = sin( angle_rad )*f;

cc = ampl*cos( aa*X + bb*Y +  phase_rad );
return
%--------------------------------------------------------------------------

