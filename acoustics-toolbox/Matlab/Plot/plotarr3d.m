function plotarr3d( ARRFile, irtheta, irr, ird, isd, Arr, Pos )

% plot the arrivals calculated by BELLHOP3D
%
% Usage: plotarr3d( ARRFile, irtheta, irr, ird, isd, Arr, Pos )
%
%   ARRFile - name of the Bellhop3D Arrivals File
%   irtheta - index of receiver bearing
%   irr     - index of receiver range
%   ird     - index of receiver depth
%   isd     - index of source   depth
%   Arr     - Arrival structure array (output of read_arrivals_xxx, OPTIONAL)
%   Pos     - Src/Rcv Position structure (output of read_arrivals_xxx, OPTIONAL)
%
% mbp, Apr 2009
% jcp, Jul 2018

% read the Bellhop Arrivals File (unless the Arr, Pos structs were given)

if nargin < 7
  [ Arr, Pos ] = read_arrivals_asc( ARRFile );
  %[ Arr, Pos ] = read_arrivals_bin( ARRFile );
end

%%
% amp-time stem plot

figure
Narr = Arr( irr, ird, irtheta, isd ).Narr;
Amps = abs( Arr( irr, ird, irtheta, isd ).A( 1:Narr  ) );
Delays = real( Arr( irr, ird, irtheta, isd ).delay( 1:Narr ) );

stem( Delays, Amps )

% [ foo, ii ] = sort( Delays );
% doo = Amps( ii );
% plot( foo, doo )

xlabel( 'Time (s)' )
ylabel( 'Amplitude' )
title( ['Sd = ',      num2str( Pos.s.z( isd ) ), ...
  ' (m)   Rtheta = ', num2str( Pos.r.theta( irtheta ) ), ...
  ' (deg)   Rd = ',   num2str( Pos.r.z( ird ) ), ...
  ' (m)   Rr = ',     num2str( Pos.r.r( irr ) ), ' (m)' ] )

%%
% depth-time stem plot

figure
for ird1 = 1 : size( Arr, 2 )
   Narr = Arr( irr, ird1, irtheta, isd ).Narr;
   stem3( real( Arr( irr, ird1, irtheta, isd ).delay( 1:Narr ) ), ...
          Pos.r.z( ird1 ) * ones( Narr , 1 ), ...
          abs( Arr( irr, ird1, irtheta, isd ).A( 1:Narr ) ) )
   hold on;
end

xlabel( 'Time (s)' )
ylabel( 'Depth (m)' )
zlabel( 'Amplitude' )
title( ['Sd = ',      num2str( Pos.s.z( isd ) ), ...
  ' (m)   Rtheta = ', num2str( Pos.r.theta( irtheta ) ), ...
  ' (deg)   Rr = ',   num2str( Pos.r.r( irr ) ), ' (m)' ] )

%%
% range-time stem plot

figure
for irr1 = 1 : size( Arr, 1 )
   Narr = Arr( irr1, ird, irtheta, isd ).Narr;
   stem3( real( Arr( irr1, ird, irtheta, isd ).delay( 1:Narr ) ), ...
          Pos.r.r( irr1 ) * ones( Narr , 1 ), ...
          abs( Arr( irr1, ird, irtheta, isd ).A( 1:Narr ) ) )
   hold on
end

xlabel( 'Time (s)' )
ylabel( 'Range (m)' )
zlabel( 'Amplitude' )
title( ['Sd = ',      num2str( Pos.s.z( isd ) ), ...
  ' (m)   Rtheta = ', num2str( Pos.r.theta( irtheta ) ), ...
  ' (deg)   Rd = ',   num2str( Pos.r.z( ird ) ), ' (m)' ] )
