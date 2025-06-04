% Code to load one given OCM dataset (data format version V2).
%
% Bruno Madore, Harvard Medical School, Advanced Lab for MRI and Acoustics (ALMA),
% Brigham and Women's Hospital, Radiology, 2023.

function [S_cplx S dphi dT tstamp tags par] = load_OCMdata(fname);

eps = 1e-6;                            % Arbitrarily small number

% Define file names
fname_raw = sprintf('%s%s', fname, '.ocm');
fname_proc = sprintf('%s%s', fname, '.ocm.proc');
fname_params = sprintf('%s%s', fname, '.params');

% The acquisition parameters in fname_params were saved in Scilab
% format and cannot be readily read from Matlab, for this reason
% they are hardcoded below instead. The name of each variable
% corresponds to a field in the structure params_saved, for example,
% bytesize should have been set as bytesize = params_saved.bytesize.
%load(fname_params);
par.Nocm = 1;			 % Number of sensors
par.c = 1540;                    % Should most likely be 1540 m/s
%%%%%%%%%%%%%%%%%%%%%
par.attn = 0.5;                  % Should most likely be 0.5 db/MHz/cm
par.attn = 0;                  % Should most likely be 0.5 db/MHz/cm
%%%%%%%%%%%%%%%%%%%%%
par.dt = 1e-7;                   % 1/(sampling frequency)
par.F0_all = 1;                  % Nominal frequency of each sensors, in MHz
par.bytesize = 8;                % Should most likely be 8
par.npts = 2000;                 % Number of data points per readout
par.Ndepths = 3;                 % Number of reconstructed depths, 'advanced' screen
par.max_Tmatch = 5;		 % 'Advanced' screen

% Open the raw and processed data files.
fpraw = fopen(fname_raw, 'rb');   % Open the raw data file for reading
fpproc = fopen(fname_proc, 'rb'); % Open the processed data file for reading
fseek(fpproc,0,'bof');
vproc = fread(fpproc,1,'double'); % Read the file version
if (sum(vproc==[2.0 2.1]) ~= 1)   % V 2.0 and 2.1 are previous and current versions
   error('Wrong version number %0.1f for OCM software/hardware version', vproc);
end
   
% Figure out the length of a readout in bytes, rl_raw and rl_proc, and the number
% of readouts, NT.
rl_raw = (2+par.npts)*par.bytesize;% Record length per readout, record+timestamp+readout
if (vproc == 2.0)                  % ctrig was just zeros in V 2.0
   rl_proc_short = (3+4*(par.Ndepths+1)+(par.max_Tmatch+1)+2*par.npts)*par.bytesize;% record+timestamp+freq+v+z+zc+Tmatch+ctrig+cplx
   rl_proc_long = rl_proc_short;
   Nctrig1 = par.Ndepths+1;
   Nctrig2 = par.Nocm;
elseif (vproc == 2.1)              % The cardiac application was implemented in V 2.1
   rl_proc_short = (3+3*(par.Ndepths+1)+(par.max_Tmatch+1)+2*par.npts)*par.bytesize; % Without ctrig
   rl_proc_long = (3+3*(par.Ndepths+1)+3+(par.max_Tmatch+1)+2*par.npts)*par.bytesize;% With ctrig
   Nctrig1 = 3;
   Nctrig2 = 1;
end

fseek(fpraw,0,'eof');
flength = ftell(fpraw);
NT_1 = floor(flength / (par.Nocm*rl_raw));
fseek(fpproc,0,'eof');
flength = ftell(fpproc);
NT_2 = floor(flength / ((par.Nocm-1)*rl_proc_short + rl_proc_long));
% Make sure the file lengths make sense
if (NT_1 == NT_2)
   NT = NT_1;
else
   error('The raw and processed files do not seem to have compatible lengths');
end

% Define the arrays
readout = zeros(par.npts,par.Nocm,NT);
readout_real = zeros(par.npts,par.Nocm,NT);
readout_imag = zeros(par.npts,par.Nocm,NT);
v = zeros(NT,par.Ndepths+1,par.Nocm); % The +1 is because of the 'whole range' option
z = zeros(NT,par.Ndepths+1,par.Nocm);
zc = zeros(NT,par.Ndepths+1,par.Nocm);
Tmatch = zeros(NT,par.max_Tmatch+1,par.Nocm); % The +1 is to store Tnow
ctrig = zeros(NT,Nctrig1,Nctrig2);   
tstamp = zeros(par.Nocm,NT);              % Time stamp for individual readouts (oscilloscope time)
tstamp_proc = zeros(par.Nocm,NT);         % Time stamp for individual readouts (oscilloscope time)
recnum_raw = zeros(par.Nocm);		  % Record number from raw data file

% Set filters/weights associated with the transmit gain correction (TGC)
tgc = zeros(par.npts,par.Nocm);           % Weights for TGC, initialize to zero
t = ((0:par.npts-1)*par.dt)';             % t axis, in s
for (iocm = 1:par.Nocm)
   tgc_dB = par.attn*par.F0_all(iocm)*(par.c*100)*t; % TGC, in dB; c*100 is in cm/s
   tgc(:,iocm) = 10.^(tgc_dB./10);        % TGC, as a multiplicative factor
end

% Read the data.
for (Tnow = 1:NT)
   % Read the raw OCM entry: a timestamp and a readout per OCM channel.
   fseek(fpraw,(Tnow-1)*par.Nocm*rl_raw,'bof');
   for (iocm = 1:par.Nocm)
      recnum_raw(iocm) = fread(fpraw,1,'double');% Record number
      tmp = fread(fpraw,1,'double');             % Time stamp
      tstamp(iocm,Tnow) = tmp;
      tmp = fread(fpraw,par.npts,'double');
      readout(1:length(tmp),iocm,Tnow) = tmp;
      readout(:,iocm,Tnow) = readout(:,iocm,Tnow).*tgc(:,iocm); % Apply TGC
   end

   % Read the processed OCM data: timestamp, frequency, v, z and cplx readout.
   % The first 'bytesize' entry was the file version number.
   fseek(fpproc,par.bytesize+(Tnow-1)*((par.Nocm-1)*rl_proc_short + rl_proc_long),'bof');
   for (iocm = 1:par.Nocm)
      recnum_proc(iocm) = fread(fpproc,1,'double');% Record number
      tmp = fread(fpproc,1,'double');    % Timestamp, read as 'double'
      tstamp_proc(iocm,Tnow) = tmp;      % Should be the same as tstamp
      tmp = fread(fpproc,1,'double');
      f(iocm) = tmp;                     % Best current guess at OCM freq
      tmp = fread(fpproc,(par.Ndepths+1),'double');
      v(Tnow,:,iocm) = tmp;              % Tissue velocity
      tmp = fread(fpproc,(par.Ndepths+1),'double');
      z(Tnow,:,iocm) = tmp;              % Tissue displacement
      tmp = fread(fpproc,(par.Ndepths+1),'double');
      zc(Tnow,:,iocm) = tmp;             % (Corrected) tissue displacement
      tmp = fread(fpproc,(par.max_Tmatch+1),'double'); % These are integers, not doubles
      Tmatch(Tnow,:,iocm) = tmp;         % T index for best past matches
      if (vproc==2.0)
         tmp = fread(fpproc,Nctrig1,'double');
         ctrig(Tnow,:,par.Nocm) = tmp;          % No cardiac functionality in V 2.0
      end
      tmp = fread(fpproc,2*par.npts,'double');
      readout_real(:,iocm,Tnow) = tmp(1:2:2*par.npts);
      readout_imag(:,iocm,Tnow) = tmp(2:2:2*par.npts);
      readout_real(:,iocm,Tnow) = readout_real(:,iocm,Tnow).*tgc(:,iocm); % Apply TGC
      readout_imag(:,iocm,Tnow) = readout_imag(:,iocm,Tnow).*tgc(:,iocm); % Apply TGC
      if (vproc==2.1)
         if (iocm==par.Nocm)
            tmp = fread(fpproc,Nctrig1,'double');
            ctrig(Tnow,:,1) = tmp;          % Cardiac-related information
         end
      end
   end

   % Check that everything seems in sync between the two files (raw and proc), and
   % that no NaNs corrupt the processed parameters (v, z and zc).
   if (sum(abs(tstamp_proc(:,Tnow) - tstamp(:,Tnow))) > eps)
      fprintf('WARNING: Time mismatch between time stamps from raw and proc files\n')
   end
   if (sum(abs(recnum_raw(:) - recnum_proc(:))) > eps)
      fprintf('WARNING: Record number mismatch between raw and proc files\n')
   end
   % Make sure there are no NaNs in the processed parameters that were just read
   tmp = v(Tnow,:,:) + z(Tnow,:,:) + zc(Tnow,:,:);
   if (sum(sum(isnan(tmp))) > 0)
      error('WARNING: NaNs were found in processed parameters');
   end
end
% Close the files from which OCM readouts and processed data were read
fclose(fpraw);
fclose(fpproc);

% The raw data may be a little more 'raw' than necessary, process the US signals
% to get the magnitude signal S, the phase difference dphi, and dT.
S_cplx = readout_real + 1i.*readout_imag; 
S = abs(S_cplx);
dphi = zeros(par.npts,par.Nocm,NT);
dphi(:,:,2:NT) = angle(S_cplx(:,:,2:NT).*conj(S_cplx(:,:,1:NT-1))); % dphi along T
dT = zeros(1,1,NT);
dT(2:NT) = tstamp(2:NT) - tstamp(1:NT-1);
par.NT = NT;

% Time tags, for syncronization with PET data, were made by double-firing the
% system every 3.0 +/- 0.2 s. Find the timing of these tags.
dT_typ = median(dT(:)); 
dT(1) = dT_typ;		% To avoid divisions by zero
% Look for firings that occurred about twice to quickly, using <60% as
% the criterion here.
tag_locs = find(dT(:)<0.6*dT_typ);
tags = zeros(floor(length(tag_locs)/2),1);
% Make sure the tag_locs make sense; double firings should always occur
% in, well, double, so tag_locs should be even and consist of a series
% of two consecutive events.
if (mod(length(tag_locs),2)==0 & min(tag_locs(1:2:end)==(tag_locs(2:2:end)-1)))
   tags(:,1) = tag_locs(1:2:end);	% Index location for the sync tags
   tags(:,2) = tstamp(tags(:,1));	% Time stamp for the sync tags
else
   error('Ran into problems figuring out the PET synchronization tags\n');
end
