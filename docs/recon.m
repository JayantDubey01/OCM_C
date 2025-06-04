% Code to load OCM data and generate respiratory bins, associated
% with different motion states.
%
% Bruno Madore, Harvard Medical School, Advanced Lab for MRI and Acoustics (ALMA),
% Brigham and Women's Hospital, Radiology, 2025.
clear
% Define variable that determine the behavior of the algorithm
d_min = 2;      % Min depth to look at, in cm
d_max = 8;      % Max depth to look at for phase info, in cm
Nb_raw = 5000;	% Initial number of bins, a relatively large number
Nb_targ1 = 100;	% Temp number of motion states, or bins, after step 1
Nb_targ2 = 30;  % Final number of motion states, or bins, after step 2
% Set the seed algorithm for the random number generator, so that
% the 'randperm' function used here reliable always gives the
% same random selection of numbers.
rng(1,'twister');
% Variables related to saving and displaying the results
fname_templ = 'dataset';	% Template for output file names
Nrow = 4;	% Number of rows to be used in tiled display
Nx_fig = 1500;	% For display purposes, size of figures along x
Ny_fig = 1000;	% For display purposes, size of figures along y

[sets stemdir] = inputdatasets();
isets = 1:8;	% List of datasets to reconstruct
for (iset = isets)
   tic		% Let's keep track of how long it takes
   % Pick, and then load, the dataset to be processed
   S0 = sets(iset);
   fprintf('Loading data from %s, %s%s\n', S0.tag, S0.dir, S0.fname);
   fname = sprintf('%s%s', S0.dir, S0.fname);
   [S_cplx Smag dphi_raw dT tstamp tags par] = load_OCMdata(fname);
   par.period = 1/(par.F0_all*1e6);
   NT = par.NT;
   Nt = par.npts;
   Nocm = par.Nocm;
   % Compute the phase increment along t
   phi = zeros(Nt,Nocm,NT);
   phi(2:Nt,:,:) = angle(S_cplx(2:Nt,:,:).*conj(S_cplx(1:Nt-1,:,:)));

   % Generate a de-modulated version of the complex signals
   fprintf('Generating de-modulated complex signals\n');
   median_phi = median(phi(:));
   par.F0_act = -par.F0_all*median_phi*par.period/par.dt/(2*pi);
   % Demodulate
   phi = phi - median_phi;
   % Select a subset of the t axis. Early t points represent signals from the
   % hardware and capsule themselves and are not relevant to patient motion, while
   % signals may get very weak and noisy for later/deeper data.
   lambda = 1e3*par.c/(par.F0_all*1e6); % Wavelength, in mm 
   t_min = round(2*(d_min/100)/par.c/par.dt);   % t point corresponding to d_min
   t_max = round(2*(d_max/100)/par.c/par.dt);   % t point corresponding to d_max
   trange = (t_min:t_max)'; % column vector with subset of t indeces being considered
   Nt_ = length(trange);
   % Make a demodulated complex entity that binning (below) will be based on.
   S = abs(S_cplx(trange,:,:)) .* exp(1i*phi(trange,:,:));

   % Pick Nb_raw random time points as an initial set of motion states
   % and associate all other time points to these motion states based
   % on similarity.
   fprintf('Associating all time points to an initial, large set of motion states\n');
   fprintf('    Processing time point (out of %6d) #      ', NT);
   states_raw = sort(randperm(NT,Nb_raw)); % Initial set of motion states
   states_list = zeros(round(Nt_/2),1,Nb_raw);
   states_list(1,1,:) = states_raw;
   states_v = S(:,:,states_raw);	% Representative vectors for each motion state
   norm = sqrt(sum(abs(states_v).^2,1)); % 'Vector length' for each motion state
   states_u = states_v./repmat(norm,[Nt_ 1 1]); % Unit vector version
   N_list = ones(Nb_raw,1);         % Number of time points per state
   for (iT = 1:NT)
      fprintf('\b\b\b\b\b\b%6d', iT);
      % Check whether this time point has already been assigned a state
      if (sum(states_raw == iT) == 0)
         % This time point has not yet been assigned, test where it belongs.   
         test_v = S(:,:,iT);	% Test vector, not normalized yet
         test_u = test_v./sqrt(sum(abs(test_v).^2,1)); % Unit vector version
         test = abs(sum(repmat(test_u,[1 1 Nb_raw]).*conj(states_u),1));
%test = sum(abs(repmat(test_u,[1 1 Nb_raw]).*conj(states_u)),1);
         % Find which motion state this time point was most similar to
         [val loc] = max(test);
         % Add this time point to the list of time point(s) associated with
         % this motion state.
         N_list(loc) = N_list(loc) + 1; % One more time point for this state
         states_list(N_list(loc),1,loc) = iT;
      end
   end
   fprintf('\n');
   states_list_safe = states_list;	% Make a safe copy
   N_list_safe = N_list;		% Make a safe copy
   states_u_safe = states_u;		% Make a safe copy

   % Consolidate motion states by repeatedly picking the smallest
   % one and merging it with the one most similar to it, until we
   % reach the intermediary target for the number of bins, Nb_targ1.
   % desired number.
   Nb = Nb_raw;			% Initialize
   fprintf('Consolidation step 1: Going from %d to %d motion states\n', ...
   Nb_raw, Nb_targ1);
   fprintf('    Number of motion states left:     ');
   while (Nb > Nb_targ1)
      % Identify the motion state with the least members, to
      % be assimilated by another.
      [val loc1] = min(N_list);
      % Get all members of this motion state
      n = N_list(loc1);
      to_move = states_list(1:n,1,loc1);
      % Get the unit vector descriptive of this particular motion state
      test_u = states_u(:,:,loc1);
      % Remove this motion state from existence, so that it could
      % not be compared with itself in 'test' below.
      N_list(loc1) = [];
      states_list(:,:,loc1) = [];
      states_u(:,:,loc1) = [];
      Nb = Nb - 1;
      % Find where the member(s) of this removed motion state
      % should best go.
      test = sum(abs(repmat(test_u,[1 1 Nb]).*conj(states_u)),1);
      [val loc2] = max(test);
      states_list(N_list(loc2)+1:N_list(loc2)+n,:,loc2) = to_move;
      N_list(loc2) = N_list(loc2) + n;
      fprintf('\b\b\b\b\b%5d', Nb);
   end
   fprintf('\n');
   % A sanity checks.
   if (Nb ~= Nb_targ1)
      error('The number of motion states is %d but should have been %d\n', ...
      Nb, Nb_targ1);
   end

   % The current ordering of the bins is random; instead, now re-
   % order the bins according to the number of time points associated
   % with them, in decreasing order. This ordering is needed for
   % the next step, so that when two similar bins will get merged,
   % the 'states_u' vector of the larger bin will be the one that
   % will be kept to describe/characterize the new, merged bin.
   [N_list bin_order] = sort(N_list,'descend');
   states_list = states_list(:,:,bin_order);
   states_u = states_u(:,:,bin_order);

   % Continue consolidating states, but now instead of focusing on
   % eliminating states with the least time points, the focus here is
   % now instead on combining states that are most similar, regardless
   % of size. Keep doing so until we reach the final target, Nb_targ2.
   fprintf('Consolidation step 2: Going from %d to %d motion states\n', ...
   Nb, Nb_targ2);
   % Calculate the correlation matrix between the different remaining
   % motion states.
   fprintf('    Calculating the correlation matrix\n');
   corr_matrix = zeros(Nb_targ1,Nb_targ1); % Correlation matrix
   for i1 = 1:Nb_targ1
      for i2 = i1+1:Nb_targ1
         test_u1 = states_u(:,:,i1);
         test_u2 = states_u(:,:,i2);
         corr12 = abs(sum(test_u1.*conj(test_u2),1));
         corr_matrix(i1,i2) = corr12;
      end
   end
   fprintf('    Number of motion states left:     ');
   while (Nb > Nb_targ2)
      % Find the maximum location in the correlation matrix.
      [val loc] = max(corr_matrix(:));
      row = mod(loc, Nb);
      col = floor(loc/Nb) + 1;
      % A sanity checks.
      if (val ~= corr_matrix(row,col))
         error('The max value %0.2f does not match the corr value %0.2f\n', ...
         val, corr_matrix(row,col));
      end
      % The row and column of the correlation matrix directly
      % represent the bins/motion states that should be merged.
      i1 = row;
      i2 = col;
      % Get all members of this motion state i2, to be eliminated.
      n = N_list(i2);
      to_move = states_list(1:n,1,i2);
      % Remove the motion state i2 from existence; i1 is the one
      % being preserved, because it has a higher number of associated
      % time points.
      N_list(i2) = [];
      states_list(:,:,i2) = [];
      states_u(:,:,i2) = [];
      Nb = Nb - 1;
      % The motion state i2 also needs to be removed from the
      % correlation matrix.
      corr_matrix(:,i2) = [];
      corr_matrix(i2,:) = [];
      % Give all the time points that were associated with bin i2
      % to bin i1 instead.
      states_list(N_list(i1)+1:N_list(i1)+n,:,i1) = to_move;
      N_list(i1) = N_list(i1) + n;
      fprintf('\b\b\b\b\b%5d', Nb);
   end
   fprintf('\n');

   % Make a couple sanity checks.
   if (Nb ~= Nb_targ2)
      error('The number of motion states is %d but should have been %d\n', ...
      Nb, Nb_targ2);
   end
   if (sum(N_list) ~= NT)
      error('The number of assigned time points is %d but should have been %d\n', ...
      sum(N_list), NT);
   end

   % Once again, order them according to the number of time points
   % associated with them, it could have changed since last time due
   % to the consolidation step above.
   [N_list bin_order] = sort(N_list,'descend');
   states_list = states_list(:,:,bin_order);
   states_u = states_u(:,:,bin_order);

   % Save to file the main results for this particular dataset.
   fname = sprintf('%s%s%03d', S0.dir, fname_templ, iset);
   fprintf('Saving results to file %s\n', fname);
   save(fname, 'S', 'states_list', 'states_u', 'N_list', 'tstamp', 'tags', 'par');

   % Generate a figure with all Nb motion states shown in a tiled
   % layout.
   figure(iset);
   [A] = get(gcf, 'Position');
   set(gcf, 'Position', [A(1) A(2) Nx_fig Ny_fig]);
   set(gca,'position',[0 0 1 1]);
   % Find a reasonable scaling
   colormap jet;
   map = colormap;
   sc = 1.0*size(map,1)/max(max(abs(S)));
   % Assuming we use Nrow rows, find how many columns are needed.
   Ncol = ceil(Nb/Nrow);
   % Initiate a tiled display
   t = tiledlayout(Nrow, Ncol);
   % In the first tile, put an example of the free-breathing, unbinned
   % data.
%  nexttile(t);
%  image(sc*squeeze(abs(S(:,1,1:5000))))
   % Go through all motion states and display each one of them
   for (ib = 1:Nb)
      nexttile(t);
      imagesc(squeeze(abs(S(:,1,states_list(1:N_list(ib),1,ib)))))
   end

   % Get ready for next pass
   clear N_list S S_cplx Smag dphi_raw norm phi states_list states_list_safe;
   clear states_raw states_u states_u_safe states_v;
   pause(0.1);	% Just to make sure the figure does get displayed
   t_end = toc;
   fprintf('It took %02dmin:%02ds to reconstruct dataset #%d\n\n', ...
   floor(t_end/60), round(mod(t_end,60)), iset);

end % Loop on datasets, i.e., patients 

