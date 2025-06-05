#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <mtwist.h>
#include <string.h>

#define PI 3.14159

int d_min = 2;      // Min depth to look at, in cm
int d_max = 8;      // Max depth to look at for phase info, in cm
int Nb_raw = 5000;	// Initial number of bins, a relatively large number
int Nb_targ1 = 100;	// Temp number of motion states, or bins, after step 1
int Nb_targ2 = 30; // Final number of motion states, or bins, after step 2
// Random function generator

char fname_templ[] = "dataset";
int Nrow = 4;	// Number of rows to be used in tiled display
int Nx_fig = 1500;	// For display purposes, size of figures along x
int Ny_fig = 1000;	// For display purposes, size of figures along y


int recon() {
    // Generate random number with Mersenne Twister Generator
    unsigned long rand_int = genrand_int32();
    
    // An array of structs and filepath strings returned by the 'inputdatasets' script
    struct S {
        char dir[10];
        char fname[10];
        char tag[10];
    };
    struct S sets[8];
    char stemdir[] = "/datadrive/OCM/PET/";

    // A struct returned by the 'load_OCMdata' for S0 (maybe it inputs loaddatasets struct and generates a new array of structs)
    struct pars {
        int NT;
        int Nocm;
        int c;
        int attn;
        int dt;
        int F0_all;
        int F0_act;
        int bytsize;
        int npts;
        int Ndepths;
        int max_Tmatch;
        int period;
    };

    struct OCMdata {
        int S_cplx[10];
        int Smag;
        int dphi_raw;
        int dT;
        int tstamp;
        int tags;
        struct pars par;
        
    };
    struct OCMdata O[8];

    int isets = sizeof(sets) / sizeof(sets[0]); // Total bytes of sets / 1 struct (assuming all structs are of the same size)
    
    for (int i = 1; i<=isets; i++) {
        struct S S0 = sets[i];
        printf('Loading data from %s, %s%s\n', S0.tag, S0.dir, S0.fname);
        char fname[] = sprintf('%s%s', S0.dir, S0.fname);

        struct OCMdata OCM0[i];
        OCM0->par.period = 1/(OCM0->par.F0_all*1000000);
        int NT = OCM0->par.NT;
        int Nt = OCM0->par.npts;
        int Nocm = OCM0->par.Nocm;
        // Compute the phase increment along t
        int phi[Nt][Nocm][NT];
        memset(phi, 0, sizeof(phi)); // Initialize to 0
        
        // This line slices the Array for all timestamps and multiplies the angle with its conjugate
        // Need to include our own impl in 'utility.c' (hint: use pointers for slicing)
        //phi(2:Nt,:,:) = angle(S_cplx(2:Nt,:,:).*conj(S_cplx(1:Nt-1,:,:)));

        // Generate a de-modulated version of the complex signals
        printf('Generating de-modulated complex signals\n');
        // Need to include our own median impl in 'utility.c'
        int median_phi = median(phi);

        OCM0->par.F0_act = ((-1 * OCM0->par.F0_all*median_phi*OCM0->par.period) / OCM0->par.dt) / (2*PI);
        
        // Demodulate
        int phi = phi - median_phi;
        // Select a subset of the t axis. Early t points represent signals from the
        // hardware and capsule themselves and are not relevant to patient motion, while
        // signals may get very weak and noisy for later/deeper data.
        double lambda = 1e3*OCM0->par.c/(OCM0->par.F0_all*1e6); // Wavelength, in mm 
        double t_min = round(2*(d_min/100)/OCM0->par.c/OCM0->par.dt);   // t point corresponding to d_min
        double t_max = round(2*(d_max/100)/OCM0->par.c/OCM0->par.dt);   // t point corresponding to d_max
    
        // Need to include our own [start:stop] range function impl in 'utility.c'
        //trange = (t_min:t_max)
        int trange[10];
        int Nt_ = sizeof(trange) / sizeof(trange[0]);
        // Make a demodulated complex entity that binning (below) will be based on.
        // Slice S_cplx with trange, and perform element-wise multiplication with e^-i(phi[trange,:,:])
        //S = abs(S_cplx(trange,:,:)) .* exp(1i*phi(trange,:,:));
        int S[10];

        // Pick Nb_raw random time points as an initial set of motion states
        // and associate all other time points to these motion states based
        //  on similarity.
        fprintf('Associating all time points to an initial, large set of motion states\n');
        fprintf('    Processing time point (out of %6d) #      ', NT);
        
        // randperm(n,k) returns a row vector containing k unique integers selected randomly from 1 to n. Impl. in 'utility.c'
        double states_raw[] = sort(randperm(NT,Nb_raw)); // Initial set of motion states
        int Nt_round = round(Nt_/2);
        double states_list[Nt_round][1][Nb_raw];
        memset(states_list, 0, sizeof(states_list));

        // Replace original values (use pointer arithmetic)
        //states_list(1,1,:) = states_raw;
        //states_v = S(:,:,states_raw);	% Representative vectors for each motion state
        double states_list[10];
        double states_v[10];
        
        // Write a norm function in 'utility.c'
        //norm = sqrt(sum(abs(states_v).^2,1)); % 'Vector length' for each motion state
        double norm[10];

        //states_u = states_v./repmat(norm,[Nt_ 1 1]); % Unit vector version
        double states_u[10];

        double N_list[Nb_raw];
        memset(N_list, 1, sizeof(N_list));

        for (int iT = 1;iT<=NT,iT++)
        {        
            fprintf('\b\b\b\b\b\b%6d', iT);
            // Check whether this time point has already been assigned a state
            if (sum(states_raw == iT) == 0)
            {
                // This time point has not yet been assigned, test where it belongs.   
                double test_v[10];
                // test_v = S(:,:,iT);	% Test vector, not normalized yet
                double test_u[10];
                //test_u = test_v./sqrt(sum(abs(test_v).^2,1)); % Unit vector version
                double test[10];
                //test = abs(sum(repmat(test_u,[1 1 Nb_raw]).*conj(states_u),1));
            }
            //test = sum(abs(repmat(test_u,[1 1 Nb_raw]).*conj(states_u)),1);
            // Find which motion state this time point was most similar to
            
            // Find max value and their index, store in a 2D array 
            //[val loc] = max(test);
            typedef struct {
                double value;
                int index;
            } MaxResult;

            MaxResult max;

            // Add this time point to the list of time point(s) associated with
            // this motion state.
            N_list[max.index] = N_list(max.index) + 1; // One more time point for this state
            states_list(N_list(loc),1,loc) = iT;
        }
        printf('\n');
    }



}



