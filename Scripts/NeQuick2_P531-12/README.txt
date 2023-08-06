=========================================================================
      NeQuick2 P.531-12 electron density model 
        for International Telecommunication Union (ITU-R) Study Group 3
	associated to ITU-R Recommendation P.531.       
      Release date 22 October 2011
      (Internal developers reference 2.1.0)
      This version of NeQuick is also applicable to version 14 of Recommendation ITU-R P.531

=========================================================================
> > DISCLAIMER < <
=========================================================================
     This software (NeQuick2 P.531) is meant for scientific use only.
     Please acknowledge the Aeronomy and Radiopropagation Laboratory
        of the Abdus Salam International Centre for Theoretical Physics, Trieste, Italy.
     
     This software is provided to The User "as is". ITU-R and the authors assume no responsibility whatsoever for its use by other parties, and makes no guarantees, expressed or implied, about the quality, reliability, or any other characteristic of this software. 
	Under no circumstances shall ITU, the author(s) or contributor(s) be liable for damages resulting directly or indirectly from the use, misuse or inability to use this software.
	
    The user may not modify this software in content and then present it or results derived from it as ITU-R or author(s) material. 

=========================================================================
> > Description < <
=========================================================================

 Detailed documentation is provided by the ITU Radiocommunication Bureau in the ITU-R Study Group 3 ITU-R Report ÒElectron density models and data for transionospheric radio propagationÓ

 The NeQuick was developed at the Abdus Salam International Centre for Theoretical Physics, Trieste, Italy and at the University of Graz, Austria.

 The NeQuick formulation is based on the Di Giovanni - Radicella (DGR) model which was modified 
to the requirements of the COST 238 Action PRIME to give vertical electron content from ground to 
1000 km consistent with the COST 238 regional electron content model.

 The NeQuick is a quick-run ionospheric electron density model particularly designed for 
transionospheric propagation applications. To describe the electron density of the ionosphere up to the peak of  the F2 layer, the NeQuick uses a profile formulation which includes  five semi-Epstein 
layers with modelled thickness parameters. Three profile anchor points are used: the E layer peak, the F1 peak and the F2 peak, that are modelled in terms of the ionosonde parameters foE, foF1, foF2 and M(3000)F2. 

 These values can be modelled (e.g. ITU- R coefficients for foF2, M3000) or experimentally derived. 
A semi-Epstein layer represents the model topside with a height- dependent thickness parameter empirically determined. The NeQuick2 package includes routines to evaluate the electron density along any ray-path and the corresponding Total Electron Content (TEC) by numerical integration.

 The sub-models contained in NeQuick2 use monthly average values of solar activity in two forms: 
twelve month running mean sunspot number R12 and 10.7 cm wavelength solar radio flux F107. 
The latter is considered to be the primary input parameter. The following relation between R12 and 
F107 is used:
          F107= 63.7+8.9D-4*R12**2+0.728*R12  (ITU-R recommendation)
   or   R12 = sqrt(167273.0+(flx-63.7)*1123.6)-408.99

 NeQuick2 calculates modip by third order interpolation in geographic latitude and longitude using grid 
values contained in the file modip.asc.

=========================================================================
> > Usage < <
=========================================================================

This version (NeQuick2 P.531) of the package contains the subroutines and two driver 
 programs written in FORTRAN 77, and ITUR (CCIR) coefficients, modip and R12 index files needed 
 for the model.
 Driver program:
	slQu_2.for            	usage: slQu_2
					provides a text-based interface step-by-step procedure to calculate
					electron density and total electron content along any 				            
					arbitrarily chosen straight line from a number of input parameters.
					The output file slQu.dat is created, which contains the electron 
					density in units of [m-3] along the profile (if required) and TEC 
					in units 1015 [m-2].
	NeQVal.for 			usage: NeQVal F10.7 InputFile (> OutputFile) 
					calculates Slant Total Electron Content for a given Solar Flux at 
					10.7 cm (F10.7) all lines in InputFile, where 
					each line defines the epoch and the coordinates of begin and end 
					points, the year. The results is provided in the Standard Output 
					with a line including the same line of the input file plus an addtional 
					column with the Slant Total Electron Content along the ray.

 Data files needed:
 
1)	12 CCIR coefficients files :   ccir11.asc (for January),...,ccir22.asc(for December)
2)	modip grid file:  modip.asc
 
Additional historical data file (usage is optional): 
1)	R12.dat (twelve-month smoothed mean sunspot number values from 1931 to 2009)


=========================================================================
> > Compilation < <
=========================================================================
A Makefile is provided in the package.
Information on compilation is available in COMPILING.txt


=========================================================================
> > Testing and Validation < <
=========================================================================
Information on testing and validation input/output data are available in TESTING_VALIDATION.txt

 
=========================================================================
> > Valid range for input solar indices < <
=========================================================================
 
The model provided in this package NeQuick2 P.531 and recommended by ITU-R 
Recommendation P.531 does not allow Solar Flux at 10.7 cm (F10.7) input below 63 F.U. (R12=0) 
and saturates the F10.7 at 193 F.U (R12=150) if solar flux input exceeds 193 F.U (R12=150).
In the case the NeQuick functions and subroutines are called in another programs, if a F10.7 
(or R12) value exceeds the upper limits of 193 F.U. (R12=150), the subroutine will saturate F10.7 at 
193 F.U. (R12 at 150). If a value below 63 F.U (R12=0) is input, the subroutine will stop the 
program.

It is recognised that several versions of NeQuick and NeQuick2 model exists in the public domain 
and they have been adapted for different applications and published in scientific literature. 
Application ranges from monthly-mean climatology to daily and quasi-real-time assimilation. Some of 
those applications require to extend the range of the F10.7 (or R12) limits.

For those alternative uses, the limits on F10.7 (or R12) input might be removed 
by commenting the lines 197 to 208 in the source file NeQuick_2.for (prepNeq subroutine) as follows:

c     *** flux saturation above 193 FU and blocked below 63 FU to avoid
c     unrealistic or undefined electron density values! ***
c      if (flx1 .gt. 193.0D0) then
c      flx1=193.0D0
c      write(*,'(2A/A/A)')'***WARNING! Solar flux limit F=193 (R12=150)',
c     & 		' exceeded.',
c     & 		' NeQuick saturates F to 193 units',
c     & 		' (ITU-R P.1239 reccomendation).'
c      endif
c      if (flx1 .lt. 63.0D0) then
c      write(*,'(2A/A)')'***WARNING! Solar flux below 63 FU (R12 <0)',
c     &                   'program stopped!'
c      stop
c      endif

This change would be applicable at the user's discretion in which case, if values outside the range 
[63,193] for F10.7 (or [0,150] for R12) are used, unrealistic or undefined electron densities 
may be obtained.
WARNING: USAGE OF THE MODEL OUTSIDE THE EXPECTED RANGE HAS NOT BEEN VALIDATED BY 
ITU-R STUDY GROUP 3 AND THEREFORE MAY BE USED ONLY AT THE USER'S RISK.


=========================================================================
> > Authors/Acknowledgements < <
=========================================================================
Authors of the modified Di Giovanni - Radicella (DGR) model code:
	Man-Lian Zhang and S.M. Radicella, Abdus Salam ICTP Trieste.
Author of modifications of DGR code and adaptations of CCIR coefficients code of the original 
NeQuick package: 
	Reinhart Leitinger.

Contributor authors for this version:
	Bruno Nava, Pierdavide Coisson, Johanna Wagner and Yenca Migoya Orue'

	the Abdus Salam International Centre for Theoretical Physics
	Strada Costiera 11, 34014 Trieste (TS), Italy
	E-mail: bnava@ictp.it, yenca@ictp.it          
	       
Validation/testing examples, validation driver/script, version control, documentation revisions for 
ITU-R provided by:

	Roberto Prieto-Cerdeira, Raul Orus-Perez, 
	European Space Agency 
	Noordwijk, The Netherlands
	E-mail: Roberto.Prieto.Cerdeira@esa.int, Raul.Orus.Perez@esa.int
