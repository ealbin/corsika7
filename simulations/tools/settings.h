
#ifndef tools_settings_h
#define tools_settings_h

struct Colors {
    Colors() :
        photon   (kRed),
        electron (kMagenta),
        muon     (kGreen),
        proton   (kBlue),
        neutron  (kCyan),
        ocharged (kOrange),
        oneutral (kGray),
        nuclei   (kBlack) {}

    EColor photon;
    EColor electron;
    EColor muon;
    EColor proton;
    EColor neutron;
    EColor ocharged;
    EColor oneutral;
    EColor nuclei;
};

struct DensityOpts {
    DensityOpts() :
        name   ("density_"),
        title  (""),
        xtitle ("Distance from Shower Core [km]"),
        ytitle ("Number Density [1 / cm**2]"),
        stats  (false),
        xbins  (1000),
        xmin   (0.),  /* [km] */
        xmax   (10.), /* [km] */ 
        width  (800),
        height (800) {}
    
    TString  name;
    TString  title;
    TString  xtitle;
    TString  ytitle;
    Bool_t   stats;
    Int_t    xbins;
    Double_t xmin;
    Double_t xmax;
    Int_t    width;
    Int_t    height;
};


struct SliceOpts {
    SliceOpts() :
        name   ("slice_"),
        title  (""),
        xtitle ("West-to-East [km]"),
        ytitle ("South-to-North [km]"),
        stats  (false),
        xbins  (100),
        ybins  (100),
        xmin   (-10.), /* [km] */
        ymin   (-10.), /* [km] */
        xmax   ( 10.), /* [km] */
        ymax   ( 10.)  /* [km] */ {}

    TString  name;
    TString  title;
    TString  xtitle;
    TString  ytitle;
    Bool_t   stats;
    Int_t    xbins;
    Int_t    ybins;
    Double_t xmin;
    Double_t ymin;
    Double_t xmax;
    Double_t ymax;
};

struct SpectrumOpts {
    SpectrumOpts() :
        name   ("spectrum_"),
        title  (""),
        xtitle ("Energy [eV]"),
        ytitle ("dN/dE"),
        stats  (false),
        xbins  (1000), 
        xmin   (3.),  /* 10^(xmin) [eV] */
        xmax   (21.)  /* 10^(xmax) [eV] */ {}
    
    TString  name;
    TString  title;
    TString  xtitle;
    TString  ytitle;
    Bool_t   stats;
    Int_t    xbins;
    Double_t xmin;
    Double_t xmax;
};

struct ContentOpts {
    ContentOpts() :
        name   ("content"),
        title  (""),
        ytitle ("Counts"),
        stats  (false),
        xbins  (37) {}

    TString name;
    TString title;
    TString ytitle;
    Bool_t  stats;
    Int_t   xbins;

    const char* GetLabel(Int_t vBin) {
        if (vBin == 1)
            return "#gamma";
        else if (vBin == 2)
            return "e^{+}";
        else if (vBin == 3)
            return "e^{-}";
        else if (vBin == 4)
            return "#mu^{+}";
        else if (vBin == 5)
            return "#mu^{-}";
        else if (vBin == 6)
            return "p";
        else if (vBin == 7)
            return "#bar{p}";
        else if (vBin == 8)
            return "n";
        else if (vBin == 9)
            return "#bar{n}";
        else if (vBin == 10)
            return "o^{0}";
        else if (vBin == 11)
            return "o^{#pm}";
        else if (vBin == 12)
            return "H";
        else if (vBin == 13)
            return "He";
        else if (vBin == 14)
            return "Li";
        else if (vBin == 15)
            return "Be";
        else if (vBin == 16)
            return "B";
        else if (vBin == 17)
            return "C";
        else if (vBin == 18)
            return "N";
        else if (vBin == 19)
            return "O";
        else if (vBin == 20)
            return "F";
        else if (vBin == 21)
            return "Ne";
        else if (vBin == 22)
            return "Na";
        else if (vBin == 23)
            return "Mg";
        else if (vBin == 24)
            return "Al";
        else if (vBin == 25)
            return "Si";
        else if (vBin == 26)
            return "P";
        else if (vBin == 27)
            return "S";
        else if (vBin == 28)
            return "Cl";
        else if (vBin == 29)
            return "Ar";
        else if (vBin == 30)
            return "K";
        else if (vBin == 31)
            return "Ca";
        else if (vBin == 32)
            return "Sc";
        else if (vBin == 33)
            return "Ti";
        else if (vBin == 34)
            return "V";
        else if (vBin == 35)
            return "Cr";
        else if (vBin == 36)
            return "Mn";
        else if (vBin == 37)
            return "Fe";
        else
            return "?";
    }

    Int_t GetBin(Int_t vID) {
        if (vID == 1)      // photon
            return 1;
        else if (vID == 2) // e+
            return 2;
        else if (vID == 3) // e-
            return 3;
        else if (vID == 5) // mu+
            return 4;
        else if (vID == 6) // mu-
            return 5;
        else if (vID == 14) // p
            return 6;
        else if (vID == 15) // p_bar
            return 7;
        else if (vID == 13) // n
            return 8;
        else if (vID == 25) // n_bar
            return 9;
        else if (
            vID == 7   || // pi0
            vID == 10  || // KL0
            vID == 16  || // KS0
            vID == 17  || // eta
            vID == 18  || // Lambda0
            vID == 20  || // Sigma0
            vID == 22  || // Xi0
            vID == 28  || // Sigma0_bar
            vID == 30  || // Xi0_bar
            vID == 48  || // eta'
            vID == 49  || // Phi
            vID == 50  || // omega0
            vID == 51  || // rho0
            vID == 56  || // Delta0
            vID == 60  || // Delta0_bar
            vID == 62  || // K*0
            vID == 65  || // K*0_bar
            vID == 66  || // nu_e
            vID == 67  || // nu_e_bar
            vID == 68  || // nu_mu
            vID == 69  || // nu_mu_bar
            vID == 116 || // D0
            vID == 119 || // D0_bar
            vID == 123 || // D*0
            vID == 126 || // D*0_bar
            vID == 130 || // J/psi
            vID == 133 || // nu_tau
            vID == 134 || // nu_tau_bar
            vID == 139 || // Xic0
            vID == 142 || // Sigmac0
            vID == 144 || // Sigmac0'
            vID == 145 || // Omegac0
            vID == 151 || // Xic0_bar
            vID == 154 || // Sigmac0_bar
            vID == 156 || // Xic0'_bar
            vID == 157 || // Omegac0_bar
            vID == 163 || // Sigmac*0
            vID == 173 || // Sigmac*0_bar
            vID == 176 || // B0
            vID == 179 || // B0_bar
            vID == 180 || // Bs0
            vID == 181 || // Bs0_bar
            vID == 184 || // Lambdab0
            vID == 187 || // Xib0
            vID == 190 || // Lambdab0_bar
            vID == 193 )  // Xib0_bar
            return 10;    // (other neutral)
        else if (vID >= 200)
            Int_t intPart;
            Int_t Z = (Int_t)(100. * modf(vID / 100., &intPart));
            if (Z == 1)       // H
                return 12;
            else if (Z == 2)  // He
                return 13;
            else if (Z == 3)  // Li
                return 14;
            else if (Z == 4)  // Be
                return 15;
            else if (Z == 5)  // B
                return 16;
            else if (Z == 6)  // C
                return 17;
            else if (Z == 7)  // N
                return 18;
            else if (Z == 8)  // O
                return 19;
            else if (Z == 9)  // F
                return 20;
            else if (Z == 10) // Ne
                return 21;
            else if (Z == 11) // Na
                return 22;
            else if (Z == 12) // Mg
                return 23;
            else if (Z == 13) // Al
                return 24;
            else if (Z == 14) // Si
                return 25;
            else if (Z == 15) // P
                return 26;
            else if (Z == 16) // S
                return 27;
            else if (Z == 17) // Cl
                return 28;
            else if (Z == 18) // Ar
                return 29;
            else if (Z == 19) // K
                return 30;
            else if (Z == 20) // Ca
                return 31;
            else if (Z == 21) // Sc
                return 32;
            else if (Z == 22) // Ti
                return 33;
            else if (Z == 23) // V
                return 34;
            else if (Z == 24) // Cr
                return 35;
            else if (Z == 25) // Mn
                return 36;
            else if (Z == 26) // Fe
                return 37;
            else 
                return 0;
        else
            return 11; // (other charged) 
    }

};

struct ImpactOpts {
    ImpactOpts() :
        name   ("impact"),
        title  (""),
        ytitle ("Altitude [km]"),
        stats  (false),
        xbins  (3),
        ybins  (1000),
        xmin   (0),
        ymin   (-3),
        xmax   (3),
        ymax   (2.3) {}

    TString  name;
    TString  title;
    TString  ytitle;
    Bool_t   stats;
    Int_t    xbins;
    Int_t    ybins;
    Double_t xmin;
    Double_t ymin;
    Double_t xmax;
    Double_t ymax;

    const char* GetLabel(Int_t vBin) {
        if (vBin == 1)
            return "Nitrogen";
        else if (vBin == 2)
            return "Oxygen";
        else if (vBin == 3)
            return "Argon";
        else 
            return "?";
    }
};

#endif

