#ifndef SimTracker_TrackHistory_CMSProcessTypes_h
#define SimTracker_TrackHistory_CMSProcessTypes_h

//! Struct holding legacy CMS convention for process types
struct CMS
{
    enum Process
    {
        Undefined = 0,
        Unknown,
        Primary,
	Ionisation,
	Bremsstrahlung,
	EPairProd,
        Annihilation,
	/*	MultipleScattering, */
        Photon,
        Compton,
        Conversions,
        SynchrotronRadiation,
        Hadronic,
	Neutron,
        Decay
        /* EIoni, */
        /* HIoni, */
        /* MuIoni, */
        /* MuPairProd, */
        /* EBrem, */
        /* MuBrem, */
        /* MuNucl */
    };
};

#endif
