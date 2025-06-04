#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <mtwist.h>

int d_min = 2;      // Min depth to look at, in cm
int d_max = 8;      // Max depth to look at for phase info, in cm
int Nb_raw = 5000;	// Initial number of bins, a relatively large number
int Nb_targ1 = 100;	// Temp number of motion states, or bins, after step 1
int Nb_targ2 = 30; // Final number of motion states, or bins, after step 2

// Random function generator
uint32_t rand = mt_seed();

char fname_templ[] = "dataset";
int Nrow = 4;	// Number of rows to be used in tiled display
int Nx_fig = 1500;	// For display purposes, size of figures along x
int Ny_fig = 1000;	// For display purposes, size of figures along y


int recon() {
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
        int bytsize;
        int npts;
        int Ndepths;
        int max_Tmatch;
        int period;
    };

    struct OCMdata {
        int S_cplx;
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
        // Write and test this line in the utility.c 
        
        //phi(2:Nt,:,:) = angle(S_cplx(2:Nt,:,:).*conj(S_cplx(1:Nt-1,:,:)));

        // Generate a de-modulated version of the complex signals
        printf('Generating de-modulated complex signals\n');
    }



}



