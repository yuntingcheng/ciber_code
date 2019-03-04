function mypaths=get_paths(flight)
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  function mypaths=get_paths.m
%%  May 7, 2009
%%  Mike Zemcov
%%  This function sets up the user defined paths in the Matlab branch of
%%   the ciber code tree.
%%  Inputs: flight - flight number
%%  Outputs: mypaths - a structure defining the different paths we need.
%%
%%  IMPORTANT!!!
%%  You must copy this file to a new file get_paths.m and edit that file
%%  to match your system.  Please DO NOT add your local copy of 
%%  get_paths.m to the repository.
%%  IMPORTANT!!!  
%%  
%%  The following default to MZ's path structure on ciber0, which is 
%%  a reasonable one if you want to copy it.  
%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %% lamp
  if (flight == 10000)
    mypaths.in_data   = '/data/LAMP';
    mypaths.fitsdata  = 'data/lamp/fits/';
    mypaths.framedata = 'data/lamp/framedata/';
    mypaths.cat       = 'data/lamp/cat/';
    mypaths.slopedata = 'data/lamp/slopedata/';
    mypaths.scratch   = 'scratch/';
    mypaths.lookup    = 'Matlab/lookup/';
    mypaths.plots      = 'plots/';
    mypaths.supporting = 'data/lamp/supporting/';
    mypaths.fielddata  = 'data/lamp/fielddata/';
    mypaths.flatfield  = 'data/lamp/flatfield/';
    mypaths.labdata    = '/data/CIBER/blackbrant/';
  end
  %% ciber 1st flight
  if (flight == 36226)
    mypaths.in_data   = '/data/CIBER/';
    mypaths.fitsdata  = 'data/36226/fits/';
    mypaths.framedata = 'data/36226/framedata/';
    mypaths.cat       = 'data/36226/cat/';
    mypaths.slopedata = 'data/36226/slopedata/';
    mypaths.scratch   = 'scratch/';
    mypaths.lookup    = 'Matlab/lookup/';
    mypaths.plots      = 'plots/';
    mypaths.supporting = 'data/36226/supporting/';
    mypaths.fielddata  = 'data/36226/fielddata/';
    mypaths.flatfield  = 'data/36226/flatfield/';
    mypaths.labdata    = '/data/CIBER/blackbrant/';
  end
  % ciber 2nd flight
  if (flight == 36265)
    mypaths.in_data    = '/data/CIBER/';
    mypaths.fitsdata   = 'data/36265/fits/';
    %mypaths.framedata  = 'data/36265/framedata/' ;
    mypaths.framedata  = '/Volumes/HD1TB/CIBER/data/36265/framedata/';%yt
    mypaths.cat        = 'data/36265/cats/';
    %mypaths.slopedata  =  'data/36265/slopedata/';
    mypaths.slopedata =  '/Users/ytcheng/ciber/mz/36265/slopedata/';%yt
    %mypaths.scratch    = 'scratch/';
    mypaths.scratch    = '/Volumes/HD1TB/CIBER/scratch/';%yt
    mypaths.lookup     = 'Matlab/lookup/';
    %mypaths.plots      = 'plots/';
    mypaths.plots      = '/Users/ytcheng/ciber/mz/plots/';%yt
    %mypaths.supporting = 'data/36265/supporting/';
    mypaths.supporting = '/Users/ytcheng/ciber/mz/36265/supporting/';%yt
    mypaths.fielddata  = 'data/36265/fielddata/';
    mypaths.flatfield  = 'data/36265/flatfield/';
    %mypaths.labdata    = '/data/CIBER/blackbrant/';
    mypaths.labdata    = '/Volumes/HD1TB/CIBER/data/blackbrant/';%yt
    %mypaths.release    = 'data/dr/20120625/36265/';
    mypaths.release    = '/Users/ytcheng/ciber/mz/36265/dr/';%yt
    %mypaths.noisemodel = 'data/36265/noisemodel/';
    mypaths.noisemodel = '/Volumes/HD1TB/CIBER/data/36265/noisemodel/';%yt
  end
  % ciber 3rd flight
  if (flight == 36277)
    mypaths.in_data    = '/data/CIBER/';
    mypaths.fitsdata   = 'data/36277/fits/';
    %mypaths.framedata  = 'data/36277/framedata/' ;
    mypaths.framedata  = '/Volumes/HD1TB/CIBER/data/36277/framedata/';%yt
    mypaths.cat        = 'data/36277/cats/';
    %mypaths.slopedata =  'data/36277/slopedata/';
    mypaths.slopedata =  '/Users/ytcheng/ciber/mz/36277/slopedata/';%yt
    %mypaths.scratch    = 'scratch/';
    mypaths.scratch    = '/Volumes/HD1TB/CIBER/scratch/';%yt
    mypaths.lookup     = 'Matlab/lookup/';
    %mypaths.plots      = 'plots/';
    mypaths.plots      = '/Users/ytcheng/ciber/mz/plots/';%yt
    %mypaths.supporting = 'data/36277/supporting/';
    mypaths.supporting = '/Users/ytcheng/ciber/mz/36277/supporting/';%yt
    mypaths.fielddata  = 'data/36277/fielddata/';
    mypaths.flatfield  = 'data/36277/flatfield/';
    %mypaths.labdata    = '/data/CIBER/blackbrant/';
    mypaths.labdata    = '/Volumes/HD1TB/CIBER/data/blackbrant/';%yt
    %mypaths.release    = 'data/dr/20120625/36265/';
    mypaths.release    = '/Users/ytcheng/ciber/mz/36277/dr/';%yt
    %mypaths.noisemodel = 'data/36277/noisemodel/';
    mypaths.noisemodel = '/Volumes/HD1TB/CIBER/data/36277/noisemodel/';%yt
  end
  % ciber 4th flight
  if (flight == 40030)
    mypaths.in_data    = '/data/CIBER/';
    mypaths.fitsdata   = 'data/40030/fits/';
    %mypaths.framedata  = 'data/40030/framedata/' ;
    %mypaths.framedata  = '/Volumes/HD1TB/CIBER/data/40030/framedata/';%yt
    mypaths.framedata  = '/Users/ytcheng/ciber/data/framedata/';
    mypaths.cat        = 'data/40030/cats/';
    %mypaths.slopedata =  'data/40030/slopedata/';
    mypaths.slopedata =  '/Users/ytcheng/ciber/mz/40030/slopedata/';%yt
    %mypaths.scratch    = 'scratch/';
    mypaths.scratch    = '/Volumes/HD1TB/CIBER/scratch/';%yt
    mypaths.lookup     = 'Matlab/lookup/';
    %mypaths.plots      = 'plots/';
    mypaths.plots      = '/Users/ytcheng/ciber/mz/plots/';%yt
    %mypaths.supporting = 'data/40030/supporting/';
    mypaths.supporting = '/Users/ytcheng/ciber/mz/40030/supporting/';%yt
    mypaths.fielddata  = 'data/40030/fielddata/';
    mypaths.flatfield  = 'data/40030/flatfield/';
    %mypaths.labdata    = '/data/CIBER/blackbrant/';
    mypaths.labdata    = '/Volumes/HD1TB/CIBER/data/blackbrant/';%yt
    %mypaths.release    = 'data/dr/';
    mypaths.release    = '/Users/ytcheng/ciber/mz/40030/dr/';%yt
    %mypaths.noisemodel = 'data/40030/noisemodel/';
    mypaths.noisemodel = '/Volumes/HD1TB/CIBER/data/40030/noisemodel/';%yt
    mypaths.Dark_Raw = '/Volumes/HD1TB/CIBER/data/Dark_Raw/40030/';%yt
    mypaths.filtmap = '/Volumes/HD1TB/CIBER/data/40030/filtmap/';
    mypaths.alldat = '/Users/ytcheng/ciber/doc/20170325_alldat/';
    mypaths.ciberdir = '/Users/ytcheng/ciber/';
  end
  
return
