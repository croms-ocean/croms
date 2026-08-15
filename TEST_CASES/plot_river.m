%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Make 1 animation from the results of the RIVER test case
% 
%  Further Information:  
%  http://www.croco-ocean.org
%  
%  This file is part of CROCOTOOLS
%
%  CROCOTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  CROCOTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Copyright (c) 2002-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%
step=2;
vname='salt';
makemovie=0;
makepdf=0;
%
% Read data
%
nc=netcdf('river_his.nc','r');
tis=nc{'scrum_time'}(:);
h=nc{'h'}(:);
x=(nc{'x_rho'}(:))/1000;
y=(nc{'y_rho'}(:))/1000;
N=length(nc('s_rho'));
mask=nc{'mask_rho'}(:);

mask(mask==0)=NaN;
[M,L]=size(x);
[I,J]=meshgrid([0:L-1],[0:M-1]);
hmax=max(max(h));
map=colormap(jet);
NCOL=length(map(:,1));

if makemovie==1
  tstart=1;
else 
  tstart=length(tis);  
end

for tndx=1:length(tis)
  s=squeeze(nc{vname}(tndx,N,:,:));
  u=squeeze(nc{'u'}(tndx,N,:,:));
  v=squeeze(nc{'v'}(tndx,N,:,:));
  pcolor(x,y,mask.*s)
  caxis([18 36]) 
  shading flat
  colorbar
  hold on
  [C1,h1]=contour(x,y,h,[0:25:300],'k');
  u=u2rho_2d(u);
  v=v2rho_2d(v);
  quiver(x(1:step:end,1:step:end),y(1:step:end,1:step:end),...
         u(1:step:end,1:step:end),v(1:step:end,1:step:end),1,'k');
  axis image
  axis([0 40 0 80]) 

  hold off
  title(['RIVER: ',vname,' - day = ',num2str(tis(tndx)/(24*3600))])

  if makemovie==1
    MOV(tndx) = getframe;
  end  
end

close(nc)

if makemovie==1
  movie(MOV,1);
end
if makepdf
 print -dpdf river_salt.pdf
 eval('!pdfcrop river_salt.pdf river_salt.pdf')
end

