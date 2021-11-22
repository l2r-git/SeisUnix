/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUMUTECSV: $Revision: 1.2 $ ; $Date: 2021/10/13 03:57:26 $		*/
 
#include "su.h"
#include "segy.h" 
#include "gridread.h"
#include "gridxy.h"

/*********************** self documentation ******************************/
char *sdoc[] = {
"									     ",
" SUMUTECSV - MUTE above (or below) bilinearly interpolated polygonal curves ",
"									     ",
"  sumutecsv <stdin >stdout [required parameters] [optional parameters]      ",
"									     ",
" Required Parameters:							     ",
" cdp=		   list of CDPs for which tims & offs are specified.         ",
"                  Example cdp=1,17,333. (Not needed if just 1 in list).     ",
" offs=            list of offset values. This parameter must be repeated    ",
"                  for each number in cdp= list. Values in this list must    ",
"                  be in increasing order.                                   ",
" tims=      	   list of mute time values (ms.) Note: MILLISECONDS.        ",
"                  This parameter must be repeated for each number in the    ",
"                  cdp= list. There must be the same number of tims and offs ",
"                  values for a cdp (but can be a different number than for  ",
"                  other cdps).                                              ",
"									     ",
" Optional Parameters:							     ",
"									     ",
" rfile=           If set, read a K-file containing 3D grid definition.      ", 
"                  Assume 2D if no K-file and no grid definition is          ", 
"                  supplied via command line parameters. The required        ", 
"                  command line parameters are: grid_xa,grid_ya,             ", 
"                  grid_xb,grid_yb, grid_xc,grid_yc, grid_wb, grid_wc.       ", 
"                  (See program SUBINCSV for 3D grid details).               ", 
"                  If this is a 3D then the input CDPs for the mute          ", 
"                  locations and the seismic trace CDP key should be         ", 
"                  3D grid cell numbers as produced by SUBINCSV.             ", 
"                  A 3D also forces restrictions on the locations of         ", 
"                  the input mute locations. Their CDP locations must        ", 
"                  form aligned rectangles (see Notes).                      ", 
"                  Note: no need to input a grid if cdp= list has            ", 
"                        only 1 number (or is not specified).                ", 
"									     ",
" offkey=offset    Key header word specifying trace offset                   ",
" abs=1            use the absolute value of offkey.                         ",
"               =0 do not use absolute value of offkey.                      ",
" ntaper=0         number of samples to taper (sine-squared) from            ",
"                  computed full-mute sample                                 ",
" mode=0           mute ABOVE the polygonal curves                           ",
"               =1 to zero BELOW the polygonal curves                        ",
"               =2 to mute below AND above a straight line. In this case     ",
"                  offs,tims describe the total time length of the muted     ",
"                  zone as a function of offs. the slope of the line is      ",
"                  given by vel=                                             ",
"               =3 to mute below AND above a constant velocity hyperbola     ",
"                  as in mode=2 offs,tims describe the total time length of  ",
"                  mute zone as a function of offs, the velocity is given by ",
"                  the value of vel=                                         ",
" vel=330          constant velocity for linear or hyperbolic mute           ",
" tzero=0          time shift (ms.) of linear or hyperbolic mute at the      ",
"                  offkey value of 0. Note: MILLISECONDS.                    ",
"									     ",
" extrapi=0        do not extrapolate at ends in igi direction.              ",
"                  (Mute times beyond ends are held constant).               ",
"               =1 extrapolate at both igi ends                              ",
"               =2 extrapolate only at lower igi end                         ",
"               =3 extrapolate only at higher igi end                        ",
"									     ",
" extrapc=0        do not extrapolate at ends in igc direction.              ",
"                  (Mute times beyond ends are held constant).               ",
"               =1 extrapolate at both igc ends                              ",
"               =2 extrapolate only at lower igc end                         ",
"               =3 extrapolate only at higher igc end                        ",
"									     ",
" extrapt=0        do not extrapolate at ends in offset direction.           ",
"                  (Mute times beyond ends are held constant).               ",
"               =1 extrapolate at both offset ends                           ",
"               =2 extrapolate only at lower offset end                      ",
"               =3 extrapolate only at higher offset end                     ",
"									     ",
" check=0          Do not print grid checking and function locations.        ",
"               =1 If grid definition input, run some grid functions on      ",
"                  the 4 grid corner points and print results. Also print    ",
"                  input mute function location information by using the     ",
"                  grid definition and the input cdp= values.                ",
"                  The output values are:                                    ",
"                     G,cdp,igi,igc,xgrid,ygrid,xworld,yworld.               ",
"                  (The format is intended to make it easy to review and     ",
"                  also to make it easy to load to spreadsheets).            ",
"                  This information can be written to a file by              ",
"                  putting 2>yourfile on the command line.                   ",
"									     ",
" Notes:								     ",
"									     ",
" For muting with one function only, specify the arrays                      ",
"          offs=o1,o2,... tims=t1,t2,...                                     ",
" where t1 is the time at offset o1, t2 is the time at offset o2, ...        ",
" The offsets specified in the offs array must be monotonically increasing.  ",
" Linear interpolation and constant extrapolation of the specified times     ",
" is used to compute the times at offsets not specified.                     ",
"                                                                            ",
" For muting with a function of offset and CDP, specify the array            ",
"          cdp=cdp1,cdp2,...                                                 ",
" and, for each CDP specified, specify another offs and tims array as above. ",
"									     ",
" For 3D, user needs to input mute locations (cdp numbers) which form aligned",
" rectangles. That is, how-ever-many grid inlines the user chooses to put    ",
" mute locations on, there must always be the same number of mute locations  ",
" on each inline and those functions must be located at same grid crosslines.",
" For instance, if user inputs locations for inline 7 at crosslines 15,25,40 ",
" then the user must input locations at crosslines 15,25,40 for any other    ",
" inlines that the user wants to supply locations for. (If user is lucky     ",
" enough that the grid has 100 inlines, then the input locations could be at ",
" CDPs 715,725,740 and 1115,1125,1140 and 2015,2025,2040 and 2515,2525,2540. ",
" Note that the CDP numbers specified on the cdp= parameter and also in the  ",
" input seismic trace cdp key are translated to inline and crossline numbers ",
" using the input 3D grid definition - so those cdp numbers need to          ",
" correspond to the input 3D grid definition.                                ",
"									     ",
" For trace cdps that are not explicitly in the input cdp= list, bilinear    ",
" interpolation is done if the trace cdp location is surrounded by 4 mute    ",
" functions specified in the cdp= list. If the trace cdp is not surrounded   ",
" by 4 input mute functions, the result depends on the extrapi and extrapc   ",
" options. The result can be any combination of linear interpolation and     ",
" linear extrapolation and constant extrapolation. If input mute functions   ",
" are only located along 1 inline or 1 crossline the result is linear        ", 
" interpolation in that direction (and linear or constant extrapolation      ",
" the outer ending functions).                                               ",
"									     ",
" The interpolation related to cdp is actually done after the interpolation  ",
" related to offset. That is, first the trace offset is used to compute times",
" from the offs,tims arrays for 4, 2, or 1 mute functions and then weighting ",
" based on the cdp locations of those mute functions is applied. Note also   ",
" that restricting the mute to the earliest and latest trace time is done    ",
" just before application to each trace. A consequence of this is that both  ",
" negative offs= and negative tims= values are allowed and honored even if   ",
" the abs= option is 1 (the default).                                        ",
"									     ",
NULL};

/* Amalgamated Credits:
 *	SEP: Shuki Ronen, Chuck Sword
 *	CWP: Shuki Ronen, Jack K. Cohen, Dave Hale, Bjoern Rommel, 
 *           Carlos E. Theodoro, Sang-yong Suh, John Stockwell                                 
 *      DELPHI: Alexander Koek.
 *      USBG: Florian Bleibinhaus. 
 *
 *      Merged/Modified: Oct 2021: Andre Latour   
 *	 1. This program started as a copy of sunmocsv.c which itself
 * 	    started from sunmo.c (yes, NMO). The reason to do this is    
 *	    that sumute.c just has input parameters of tmute[],xmute[].     
 *	    For 3d muting, this program needs  cdp[], tims[][], offs[][].
 *          And, as it happens, sunmocsv.c has cdp[], tnmo[][], vnmo[][]. 
 *          So the easiest thing was to start from sunmocsv.c and 
 *	    rename input parameters. Parts of sumute.c were then copied
 *	    into this program (the code for ntaper= and mode=0,1,2,3).       
 *	    The mode=4 option/code was not copied to here since it does  
 *	    not seem to make sense combined with bilinear interpolation
 *	    of CDP mute function locations (but I could be wrong).          
 *	 2. Changed to expect milliseconds for all parameter inputs.               
 *	 3. Put in error checks to stop users from accidentally         
 *	    trying to use sumute parameter names.                       
 */
/**************** end self doc *******************************************/

segy tr;

struct  VelInfo { /* Structure for mute function information */
     int *kinf;
     int nto;
     float *tims;
     float *offs;
};

struct VelInfo *RecInfo; /* Storage for all mute function location value pointers */

/* Find the 4 igi,igc locations (near) a cdp and compute their (spatial) weights. */
static void binterpfind(int kigi, int *mgi, int mgi_tot, int mgiextr,
                        int kigc, int *mgc, int mgc_tot, int mgcextr,
                        int *mgix, int *mgcx, float *wi, float *wc) ;

/* Use the 4 igi,igc near cdp, compute mute time from offs,tims arrays and (spatial) weights.*/
static void binterpvalue(float offset, int mgtextr,
                         struct VelInfo lwi, struct VelInfo hii, int mgi_tot, float wi, 
                         struct VelInfo lwc, struct VelInfo hic, int mgc_tot, float wc,
                         float *timeout) ;

/* Input offset value, and compute mute time from the offs,tims arrays at one location. */
static void linterpmute(float offset,float *offs,float *tims,int nto,int mgtextr,float *time) ;

int compSort2 (const void * q1, const void * q2) ; /* comparison function for qsort  */

int bhighi(int *all, int last, int iguy) ;             /* binary search */

int bhighf(float *all, int last, float iguy) ;         /* binary search */

#define SQ(x) ((x))*((x))

int
main(int argc, char **argv)
{

        char *key=NULL;         /* header key word from segy.h          */
        char *type=NULL;        /* ... its type                         */
        int index;              /* ... its index                        */
        Value val;              /* ... its value                        */
        float fval;             /* ... its value cast to float          */

	int ncdp;		/* number of cdps specified */
  	int *cdp=NULL;	        /* array[ncdp] of cdps */
	int icdp;		/* index into cdp array */
	int jcdp;		/* index into cdp array */
	int oldcdp;     	/* cdp of previous trace */

	int k;    		/* just a general integer               */
	int ntims;		/* number of tims specified */
	int noffs;		/* number of offs specified */

        cwp_String Rname=NULL;  /* text file name for values            */
        FILE *fpR=NULL;         /* file pointer for Rname input file    */
        double gvals[999];      /* to contain the grid definition       */
	
        float linvel;           /* linear velocity                      */
        float tm0;              /* time shift of mute=2 or 3 for key=0  */
        float *taper=NULL;      /* ...          taper values            */
        int ntaper;             /* ...          taper values            */

        int mgiextr = 0;        /* for igi extrapolation option         */
        int mgcextr = 0;        /* for igc extrapolation option         */
        int mgtextr = 0;        /* for offset extrapolation option      */

        int mode;               /* kind of mute (top, bottom, linear)   */
        int iabsoff;            /* Take absolute value of offkey        */
        cwp_Bool seismic;       /* cwp_true if seismic, cwp_false not seismic */

	/* hook up getpar */
	initargs(argc, argv);
	requestdoc(1);

/* NOTE: In this program code I make an effort to conform to the       */
/*       original indentation style (even tho I find it awkward).      */

        int nerr = 0;
	if (countparname("linvel") > 0) {
          nerr++;
          warn ("error: linvel is not a parameter. You probably meant: vel");
        }
	if (countparname("tm0")    > 0) {
          nerr++;
          warn ("error: tm0 is not a parameter. You probably meant: tzero (in MILLISECONDS)");
        }
	if (countparname("key")    > 0) {
          nerr++;
          warn ("error: key is not a parameter. You probably meant: offkey");
        }
	if (countparname("xmute")  > 0) {
          nerr++;
          warn ("error: xmute is not a parameter. You probably meant: offs");
        }
	if (countparname("tmute")  > 0) {
          nerr++;
          warn ("error: tmute is not a parameter. You probably meant: tims (in MILLISECONDS)");
        }
	if (countparname("nmute")  > 0) {
          nerr++;
          warn ("error: nmute is not a parameter.");
        }
	if (countparname("xfile")  > 0) {
          nerr++;
          warn ("error: xfile is not a parameter.");
        }
	if (countparname("tfile")  > 0) {
          nerr++;
          warn ("error: tfile is not a parameter.");
        }
	if (countparname("hmute")  > 0) {
          nerr++;
          warn ("error: hmute is not a parameter.");
        }
	if (countparname("twindow")  > 0) {
          nerr++;
          warn ("error: twindow is not a parameter.");
        }
	if (countparname("twfile")  > 0) {
          nerr++;
          warn ("error: twfile is not a parameter.");
        }
        if(nerr>0) err ("error: parameter name for SUMUTE was specified (see above).");

        if (!getparint("mode", &mode))          mode = 0;
        if(mode<0 || mode>3) err ("error: mode must be 0,1,2, or 3");

        if (!getparint("ntaper", &ntaper))      ntaper = 0;
        if(ntaper<0) err ("error: ntaper cannot be negative.");

        if (!getparint("abs", &iabsoff)) iabsoff = 1;

        if (!getparint("extrapi", &mgiextr)) mgiextr = 0;
        if(mgiextr<0 || mgiextr>3) err ("error: extrapi option not in range ");

        if (!getparint("extrapc", &mgcextr)) mgcextr = 0;
        if(mgcextr<0 || mgcextr>3) err ("error: extrapc option not in range ");

        if (!getparint("extrapt", &mgtextr)) mgtextr = 0;
        if(mgcextr<0 || mgcextr>3) err ("error: extrapt option not in range ");

/* Set up taper weights if tapering requested */

        if (ntaper) {
                taper = ealloc1float(ntaper);
                for (k = 0; k < ntaper; ++k) {
                        float s = sin((k+1)*PI/(2*ntaper));
                        taper[k] = s*s;
                }
        }       

        if (!getparfloat("vel", &linvel))    linvel = 330;
        if (linvel==0) err ("error: vel cannot be 0");
        if (!getparfloat("tzero", &tm0))          tm0 = 0;
        tm0 /= 1000.; /* so as not to change code copied from sumute, convert to seconds */

        if (!getparstring("offkey", &key)) key = "offset";
        type = hdtype(key);
        index = getindex(key);

/* ------------------------------------------------------------------- */
/* Process and set the grid definition values?                         */

        getparstring("rfile", &Rname);
      
        int maygrid;;
        gridcommand(&maygrid);
      
        int is3d = 1;
        if(maygrid==1  && Rname != NULL) err("error: input k-file not allowed when full grid on command line.");
        if(maygrid==-1 && Rname == NULL) err("error: input k-file required when partial grid on command line.");
        if(maygrid==0  && Rname == NULL) is3d = 0;

        int icheck;
        if (!getparint("check", &icheck)) icheck = 0;

        if(is3d==1) {

                if(maygrid!=1) { /* open if not full grid on command line (else pass fpR still NULL) */
                        fpR = fopen(Rname, "r");
                        if(fpR==NULL) err("error: input K-file did not open correctly.");
                }

                int errwarn = 1; /* print if error or unusual thing inside gridread */
                gridread(fpR,gvals,&errwarn); 
                if(errwarn>0) err("error reading grid (from K-file or command line parameters)");

                gridset(gvals,&errwarn); 

                if(errwarn==1) err ("gridset error: grid_wb cell width must be positive.");
                else if(errwarn==2) err ("gridset error: grid_wc cell width must be positive.");
                else if(errwarn==3) err ("gridset error: corner B is within grid_wb cell width of corner A.");
                else if(errwarn>0) err ("gridset error: returned with some unrecognized error code.");
                else if(errwarn==-1) warn ("gridset warning: corner C is near A and is reset to A.");
 
                gridcheck(gvals,icheck,&errwarn); 
                if(errwarn>0) err ("gridcheck error: returned with some unrecognized error code.");
        }

	/* get information from the first header */
	if (!gettr(&tr)) err("error: can't get first trace");

        oldcdp = tr.cdp - 1; /* make sure first trace goes thru the complete code for new cdp */

        seismic = ISSEISMIC(tr.trid);

        if (seismic) {
                if (!tr.dt) err("dt header field must be set");
        } else if (tr.trid==TRID_DEPTH) {   /* depth section */
                if (!tr.d1) err("d1 header field must be set");
        } else {
                err ("tr.trid = %d, unsupported trace id",tr.trid);
        }

	/* get tims and offs functions */
	ncdp = countparval("cdp");
	if (ncdp>0) {
		if (countparname("tims")!=ncdp)
			err("error: a tims array must be specified for each listed cdp");
		if (countparname("offs")!=ncdp)
			err("error: a offs array must be specified for each listed cdp");
	} else {
		ncdp = 1;
		if (countparname("tims")!=1)
			err("error: must be 1 tims array when 0 cdps listed");
		if (countparname("offs")!=1)
			err("error: must be 1 offs array when 0 cdps listed");
	}
  	cdp = ealloc1int(ncdp); 
        if (!getparint("cdp",cdp)) cdp[0] = tr.cdp;

/* The values from each record are going to be stored.              */
/* For quick "finding" we will sort by igi,igc (or icdp for 2d).    */
/* But we are only going to sort the pointers to the record values, */
/* not the record values themselves. The record values will stay    */
/* where they were stored during read-in.                           */
/* Allocate the memory in big chunks, then use pointer arithmatic   */
/* to divide that memory amoung the individual record pointers.     */
/* (We could, of course, simply allocate record-by-record but that  */
/*  means the memory could be spangled around - usually faster to   */
/*  keep values near each other).                                   */

        RecInfo = ealloc1(ncdp,sizeof(struct VelInfo)); 

        int *tinf = ealloc1int(ncdp*3);
	for (icdp=0; icdp<ncdp; ++icdp) RecInfo[icdp].kinf = tinf + icdp*3;

        int ierr = 0;
        if(is3d==1 && icheck>0) 
                warn("Input Mute function location information follows: G,cdp,igi,igc,  xgrid,ygrid,  xworld,yworld");

	for (icdp=0; icdp<ncdp; ++icdp) {

		ntims = countnparval(icdp+1,"tims");
		noffs = countnparval(icdp+1,"offs");

		if (noffs!=ntims || noffs==0)
			err("error: tims and offs arrays are different lengths for cdp= %d",cdp[icdp]);

                RecInfo[icdp].nto = ntims;

                RecInfo[icdp].tims = ealloc1float(ntims);
                RecInfo[icdp].offs = ealloc1float(ntims);

  		if (!getnparfloat(icdp+1,"offs",RecInfo[icdp].offs)) RecInfo[icdp].offs[0] = 0.0;
  		if (!getnparfloat(icdp+1,"tims",RecInfo[icdp].tims)) RecInfo[icdp].tims[0] = 0.0;

  		/* check that offs are increasing */
  		for(k=1; k<ntims; k++) {
                        if(RecInfo[icdp].offs[k] <= RecInfo[icdp].offs[k-1]) 
			        err("error: offs array values not in increasing order for cdp= %d",cdp[icdp]);
                }

  		/* so as not to change the code copied from sumute, convert to seconds */
  		for(k=0; k<ntims; k++) RecInfo[icdp].tims[k] /= 1000.;

                RecInfo[icdp].kinf[0] = cdp[icdp];
                if(is3d==1) {
                        gridcdpic(gvals,RecInfo[icdp].kinf[0],RecInfo[icdp].kinf+1,RecInfo[icdp].kinf+2);

                        if(RecInfo[icdp].kinf[1] < -2147483644) ierr = 1;

                        if(icheck>0) {
                                double xg;
                                double yg;
                                double xw;
                                double yw;
                                gridicgridxy(gvals,RecInfo[icdp].kinf[1],RecInfo[icdp].kinf[2],&xg,&yg);
                                gridicrawxy(gvals,RecInfo[icdp].kinf[1],RecInfo[icdp].kinf[2],&xw,&yw);
                                warn("G,%12d,%12d,%12d,  %.20g,%.20g,  %.20g,%.20g",
                                RecInfo[icdp].kinf[0],RecInfo[icdp].kinf[1],RecInfo[icdp].kinf[2],xg,yg,xw,yw);
                        }
                }
                else { /* if no 3d grid input, set inline to cdp and crossline to 1 */ 
                        RecInfo[icdp].kinf[1] = RecInfo[icdp].kinf[0];
                        RecInfo[icdp].kinf[2] = 1;
                }

	}

        if(ierr>0) err("error: At least one Input Mute function cdp is not in grid");

        checkpars(); /* I do not know what this does? Call it after all parameters are read? */

/* Sort by the 2 igi,igc values (inline and crossline grid index numbers of each cdp).       */

        qsort(RecInfo,ncdp,sizeof(struct VelInfo),compSort2);

        for (jcdp=1; jcdp<ncdp; ++jcdp) {
                if(RecInfo[jcdp-1].kinf[1] == RecInfo[jcdp].kinf[1] &&
                   RecInfo[jcdp-1].kinf[2] == RecInfo[jcdp].kinf[2]) {
                        err("error: Two mute functions input for cdp,igi,igc = %d %d %d",
                            RecInfo[jcdp].kinf[0],RecInfo[jcdp].kinf[1],RecInfo[jcdp].kinf[2]);
                }
        }

/* For bilinear interpolation, user must input mute function locations which form aligned rectangles.   */
/* That is, howevermany inlines the user chooses to put mute functions on, there must always be the     */
/* same number of functions on each inline and those functions must be located at the same crosslines.  */
/* For instance, if user inputs mute functions for inline 7 at crosslines 15,25,40 then the user must   */
/* input the functions at crosslines 15,25,40 for any other inlines that the user wants to supply mute  */
/* functions for. The following code enforces that restriction on the user input.                       */

        int mgi_tot = -1;
        for (jcdp=0; jcdp<ncdp; ++jcdp) {
                if(RecInfo[jcdp].kinf[2] != RecInfo[0].kinf[2]) {
                        mgi_tot = jcdp; /* since jcdp starts at 0 */
                        break;
                }
        }

        if(mgi_tot==-1) mgi_tot = ncdp; /* just incase all kinf[2] are the same value */

        int igc_set = RecInfo[0].kinf[2]; /* igc value of first set, so cannot match next set */
        int iset_tot = mgi_tot;
        for (jcdp=mgi_tot; jcdp<ncdp; ++jcdp) {
                if(RecInfo[jcdp-mgi_tot].kinf[1] != RecInfo[jcdp].kinf[1]) {
                        err("error: Input Mute functions are irregularly spaced (at cdp= %d)",RecInfo[jcdp].kinf[0]);
                }
                if(igc_set == RecInfo[jcdp].kinf[2]) iset_tot++;
                else {
                        if(iset_tot != mgi_tot) {
                                err("error: Not same number of input Mute functions as first set (at cpd= %d)",
                                RecInfo[jcdp].kinf[0]);
                        }
                        igc_set = RecInfo[jcdp].kinf[2];
                        iset_tot = 1;
                }
        }

        int mgc_tot = ncdp/mgi_tot;

/* OK, a brief review so as not to get confused here. We sorted on 2 values (inline and crossline */
/* numbers of each cdp). Then we checked/enforced the restriction that there always be the same   */
/* specified crossline locations ON every specified inline. That is, we made sure that we always  */
/* things like 5inlines by 3crosslines or 17inlines by 12crosslines or whatever. So, we have      */
/* effectively forced the users to specify a 2-dimensional array. Now we are going to take        */
/* advantage of that fact within binterpfind by binary-searching each dimension separately in     */
/* order to find the 4 surrounding locations that we need for bilinear interpolation.             */
/* So, allocate and copy the inline and crossline numbers from the sorted RecInfo.                */

        int *mgi = ealloc1int(mgi_tot);
        int *mgc = ealloc1int(mgc_tot);

        for (k=0; k<mgi_tot; ++k) mgi[k] = RecInfo[k].kinf[1];
        for (k=0; k<mgc_tot; ++k) mgc[k] = RecInfo[k*mgi_tot].kinf[2];

        int kigi = 0;
        int kigc = 1; /* set to 1 in case this is a 2D */ 

        int mgix = 0;
        int mgcx = 0;
        float wi = 0.;
        float wc = 0.;

/* ndxi,ndxc are the element numbers that will be computed for the stored mute functions. */
/* For the degenerate cases of just 1 inline or 1 crossline, set mgi_totdeg = 0 so that   */
/* ndxc does not get set to access elements that do not exist. For the degenerate cases   */
/* it actually does not matter which 2 extra functions are being accessed since their     */
/* weight values wi or wc will be 0.                                                      */

        int ndxi = 0; 
        int ndxc = 0;
        int mgi_totdeg = mgi_tot;
        if(mgi_tot==1 || mgc_tot==1) mgi_totdeg = 0;

	/* loop over traces */
	do {

                int nt     = (int) tr.ns;
                float tmin = tr.delrt/1000.0;
                float dt = ((double) tr.dt)/1000000.0;
                float t;
                int nmute;
                int itaper;
                int topmute;
                int botmute;
                int ntair=0;
                register int i;

                if (!seismic) {
                        tmin = 0.0;
                        dt = tr.d1;
                }

                /* get value of key and convert to float */
                gethval(&tr, index, &val);
                fval = vtof(type,val);
                if (iabsoff==1) fval = fabsf(fval);

                if(ncdp<2) { /* if just one mute function, we MUST call linterpmute directly. */ 
                        linterpmute(fval,RecInfo[0].offs,RecInfo[0].tims,RecInfo[0].nto,mgtextr,&t);
                }
                else if(tr.cdp==oldcdp) { /* just compute time for new offset */
                        binterpvalue(fval,mgtextr,RecInfo[ndxi-1],RecInfo[ndxi],mgi_tot,wi,
                                                  RecInfo[ndxc-1],RecInfo[ndxc],mgc_tot,wc,&t);
                }
                else {
                        oldcdp = tr.cdp;
                        /* compute the grid indexes igi,igc */
                        if(is3d==1) { 
                                gridcdpic(gvals,tr.cdp,&kigi,&kigc);
                                if(kigi<-2147483644) err("Error: input cdp= %d not in grid.",tr.cdp); 
                        }
                        else kigi = tr.cdp; /* if no grid, just use igi=cdp and igc=1 */
             
                        /* find input cdp (higher) locations mgix,mgcx and weights wi,wc */
                        binterpfind(kigi,mgi,mgi_tot,mgiextr,kigc,mgc,mgc_tot,mgcextr,
                                    &mgix,&mgcx,&wi,&wc);

                        /* mgix and mgcx are the locations computed for each dimension seperately.    */
                        /* From them, compute the element numbers of the stored functions in RecInfo. */
                        /* Note that for the degenerate cases of mgi_tot=1 or mgc_tot=1 the           */
                        /* mgi_totdeg=0, which results in ndxc=ndxi, which in turn means the second   */
                        /* two functions passed to binterpvalue are the same as the first two         */
                        /* (which works because either weight wi or wc will be 0.0).                  */
                        /* BUT this would not work for both mgi_tot=1 and mgc_tot=1 (just 1 function) */
                        /* which is why we MUST skip this code when there is just 1 function.         */

                        ndxi = mgix + mgi_tot * (mgcx-1); 
                        ndxc = ndxi + mgi_totdeg;

                        /* mgix and mgcx are always returned as the highest of the 2 (near) locations.*/
                        /* (So, if mgi has 10 locations, mgix is only returned from 1 to 9, not 0).   */
                        /* That means the 4 locations are ndxi-1,ndxi and ndxc-1,ndxc.                */

                        /* use the 4 locations and their offs,tims and get mute time at fval offset)  */
                        binterpvalue(fval,mgtextr,RecInfo[ndxi-1],RecInfo[ndxi],mgi_tot,wi,
                                                  RecInfo[ndxc-1],RecInfo[ndxc],mgc_tot,wc,&t);

                }

               /* do the mute */
                if (mode==0) {  /* mute above */
                        nmute = MIN(NINT((t - tmin)/dt),nt);
                        if (nmute>0) memset( (void *) tr.data, 0, nmute*FSIZE);
                        for (i = 0; i < ntaper; ++i)
                                if (i+nmute>0) tr.data[i+nmute] *= taper[i];
                        if (seismic) {
                                tr.muts = NINT(t*1000);
                        } else  {
                                tr.muts = NINT(t);
                        } 
                } else if (mode==1){    /* mute below */
                        nmute = MAX(0,NINT((tmin + nt*dt - t)/dt));
                        memset( (void *) (tr.data+nt-nmute), 0, nmute*FSIZE);
                        for (i = 0; i < ntaper; ++i)
                                if (nt>nmute+i && nmute+i>0)
                                        tr.data[nt-nmute-1-i] *= taper[i];
                        if (seismic) {
                                tr.mute = NINT(t*1000);
                        } else  {
                                tr.mute = NINT(t);
                        }
                } else if (mode==2){    /* air wave mute */
                        nmute = NINT((tmin+t)/dt);
                        ntair=NINT(tm0/dt+fval/linvel/dt);
                        topmute=MIN(MAX(0,ntair-nmute/2),nt);
                        botmute=MIN(nt,ntair+nmute/2);
                        memset( (void *) (tr.data+topmute), 0,
                                (botmute-topmute)*FSIZE);
                        for (i = 0; i < ntaper; ++i){
                                itaper=ntair-nmute/2-i;
                                if (itaper > 0) tr.data[itaper] *=taper[i];
                        }
                        for (i = 0; i < ntaper; ++i){
                                itaper=ntair+nmute/2+i;
                                if (itaper<nt) tr.data[itaper] *=taper[i];
                        }
                } else if (mode==3) {   /* hyperbolic mute */
                        nmute = NINT((tmin + t)/dt);
                        ntair=NINT(sqrt( SQ((float)(tm0/dt))+SQ((float)(fval/linvel/dt)) ));
                        topmute=MIN(MAX(0,ntair-nmute/2),nt);
                        botmute=MIN(nt,ntair+nmute/2);
                        memset( (void *) (tr.data+topmute), 0,
                                (botmute-topmute)*FSIZE);
                        for (i = 0; i < ntaper; ++i){
                                itaper=ntair-nmute/2-i;
                                if (itaper > 0) tr.data[itaper] *=taper[i];
                        }
                        for (i = 0; i < ntaper; ++i){
                                itaper=ntair+nmute/2+i;
                                if (itaper<nt) tr.data[itaper] *=taper[i];
                        }
                } /* this is where mode==4 goes in sumute */ 

	        puttr(&tr);
        } while (gettr(&tr));

	return(CWP_Exit());
}

/* Find the 4 storage locations associated with the cdp's grid index locations (kigi,kigc).  */
/* The 4 locations are: mgixo,mgcxo and mgixo-1,mgcxo and mgixo,mgcxo-1 and mgixo-1,mgcxo-1. */
/* These 4 locations do not always surround the trace cdp location. For instance, for cdps   */
/* below the minimum igi value of the input functions, mgixo is returned 1 anyway (not 0).   */
/* But the wi weight is returned as 1.0 so only the lowest igi mute function contributes to  */
/* resulting mute time from binterpvalue. Similarly, for cdps above the maximum igi value    */
/* of input mute functions, mgixo is returned as maximum BUT wi weight is returned as 0 so   */
/* only the highest mute function contributes to the resulting mute time from binterpvalue.  */
/*                                                                                           */
/* Input arguments:                                                                          */
/*                                                                                           */
/* kigi      igi number of cdp (the 3D grid inline location of the cdp)                      */
/*                                                                                           */
/* mgi       array of igi numbers of the mute functions (3D grid inline locations)           */
/*                                                                                           */
/* mgi_tot   number of values in mgi array                                                   */
/*                                                                                           */
/* mgiextr=0 no extrapolation at igi ends                                                    */
/*        =1 extrapolate both lower and higher ends                                          */
/*        =2 extrapolate only at lower end                                                   */
/*        =3 extrapolate only at higher                                                      */
/*                                                                                           */
/* kigc      igc number of cdp (the 3D grid crossline location of the cdp)                   */
/*                                                                                           */
/* mgc       array of igc numbers of the mute functions (3D grid crossline locations)        */
/*                                                                                           */
/* mgc_tot   number of values in mgc array                                                   */
/*                                                                                           */
/* mgcextr=0 no extrapolation at igc ends                                                    */
/*        =1 extrapolate both lower and higher ends                                          */
/*        =2 extrapolate only at lower end                                                   */
/*        =3 extrapolate only at higher                                                      */
/*                                                                                           */
/*                                                                                           */
/* Output arguments:                                                                         */
/*                                                                                           */
/* mgixo     mgi element number where mgi[mgixo] is usually greater than kigi (there is      */
/*           some trickyness here, read the note below)                                      */
/*                                                                                           */
/* mgcxo     mgc element number where mgc[mgcxo] is usually greater than kigc (there is      */
/*           some trickyness here, read the note below)                                      */
/*                                                                                           */
/* wi        weight in the igi direction.                                                    */
/*           This weight should be applied to the TWO mute functions associated with mgixo-1 */
/*           and (1.-wi) should be applied to the TWO mute functions associated with mgixo.  */
/*                                                                                           */
/* wc        weight in the igc direction                                                     */
/*           This weight should be applied to the TWO mute functions associated with mgcxo-1 */
/*           and (1.-wc) should be applied to the TWO mute functions associated with mgcxo.  */
/*                                                                                           */

static void binterpfind(int kigi, int *mgi, int mgi_tot, int mgiextr,
                        int kigc, int *mgc, int mgc_tot, int mgcextr,
                        int *mgixo, int *mgcxo, float *wi, float *wc) {

  if(mgi_tot==1 && mgc_tot==1) { 
    *wi    = 1.;
    *wc    = 0.;
    *mgixo = 1; 
    *mgcxo = 1;
    return;
  }

/* Note the trickyness here. We never return an mgix=0 because other  */
/* code is going to use mgix-1 for the lower "surrounding" location.  */
/* So, for kigi less than lowest, we set mgix=1. But we reset kigi so */
/* subsequent computation puts all weight on the values at mgi[0]     */
/* (unless we want to extrapolate the low end).                       */
/*                                                                    */
/* But when kigi is greater than highest, we set mgix to highest and  */
/* reset kigi so that subsequent computation puts all weight on the   */
/* values at the highest (unless we want to extrapolate the high end).*/

  int mgix;
  if(kigi<=mgi[0]) {
    mgix = 1;
    if(mgiextr==0 || mgiextr==3) kigi = mgi[0];
  }
  else if(kigi>=mgi[mgi_tot-1]) { 
    mgix = mgi_tot - 1;     
    if(mgiextr==0 || mgiextr==2) kigi = mgi[mgi_tot-1];
  }
  else {
    mgix = bhighi(mgi, mgi_tot, kigi);
  }

  if(mgc_tot==1) {
    *wi    = ((float)(mgi[mgix]-kigi)) / ((float)(mgi[mgix]-mgi[mgix-1]));
    *wc    = 0.; 
    *mgixo = mgix;
    *mgcxo = 1;
    return;
  }

/* Same trickyness next as explained above.  */

  int mgcx;
  if(kigc<=mgc[0]) {
    mgcx = 1;
    if(mgcextr==0 || mgcextr==3) kigc = mgc[0];
  }
  else if(kigc>=mgc[mgc_tot-1]) {
    mgcx = mgc_tot - 1;    
    if(mgcextr==0 || mgcextr==2) kigc = mgc[mgc_tot-1];
  }
  else {
    mgcx = bhighi(mgc, mgc_tot, kigc);
  }

  if(mgi_tot==1) { 
    *wi    = 0.;
    *wc    = ((float)(mgc[mgcx]-kigc)) / ((float)(mgc[mgcx]-mgc[mgcx-1]));
    *mgixo = 1;
    *mgcxo = mgcx;
    return;
  }

  *wi    = ((float)(mgi[mgix]-kigi)) / ((float)(mgi[mgix]-mgi[mgix-1]));
  *wc    = ((float)(mgc[mgcx]-kigc)) / ((float)(mgc[mgcx]-mgc[mgcx-1]));
  *mgixo = mgix;
  *mgcxo = mgcx;

  return;
}

/* Use the storage locations and weights found by binterpfind and compute the mute  */
/* time related to the input offset (by linear interpolation within those functions */
/* and then applying the input wi,wc weigths).                                      */
/* Note: Within SUNMOCSV.C binterpfind and binterpvalue are done within one routine,*/
/*       but here they are seperated into 2 routines because you cannot pre-compute */
/*       and store mute times for all possible offsets. The mute time for an offset */
/*       must be computed on-the-fly, but I still want to be able to skip over the  */
/*       binterpfind code when a trace belongs to same cdp as the previous trace.   */

static void binterpvalue(float offset, int mgtextr,
                         struct VelInfo lwi, struct VelInfo hii, int mgi_tot, float wi, 
                         struct VelInfo lwc, struct VelInfo hic, int mgc_tot, float wc,
                         float *timeout) {

  if(mgi_tot==1 && mgc_tot==1) { /* never get here in sumutecsv due to previous code*/ 
    linterpmute(offset,lwi.offs,lwi.tims,lwi.nto,mgtextr,timeout);
    return;
  }

  *timeout = 0.;
  float time = 0.;

  if(mgc_tot==1) {
    if(wi != 0.) { /* because of extrapolation options, check exactly 0 */
      linterpmute(offset,lwi.offs,lwi.tims,lwi.nto,mgtextr,&time);
      *timeout += wi*time;  
    }
    if(wi != 1.) { /* because of extrapolation options, check exactly 1 */
      linterpmute(offset,hii.offs,hii.tims,hii.nto,mgtextr,&time);
      *timeout += (1.0-wi)*time;  
    }
    return;
  }

  if(mgi_tot==1) { 
    if(wc != 0.) {
      linterpmute(offset,lwc.offs,lwc.tims,lwc.nto,mgtextr,&time);
      *timeout += wc*time;  
    }
    if(wc != 1.) {
      linterpmute(offset,hic.offs,hic.tims,hic.nto,mgtextr,&time);
      *timeout += (1.0-wc)*time;  
    }
    return;
  }

/* The 4 point weighting equation looks like this:           */  
/*  *timeout =  wc      * (wi*timea + (1.0-wi)*timeb)        */  
/*           + (1.0-wc) * (wi*timec + (1.0-wi)*timed);       */
/*                                                           */
/* But reduce some brackets and it looks like this:          */  
/*  *timeout =  wc*wi*timea + wc*(1.0-wi)*timeb              */  
/*           + (1.0-wc)*wi*timec + (1.0-wc)*(1.0-wi)*timed;  */
/*                                                           */
/* So we can isolate the weight factors needed for each of   */  
/* the 4 locations, as follows:                              */  

  float aw = wc*wi;
  float bw = wc*(1.0-wi);
  float cw = (1.0-wc)*wi;
  float dw = (1.0-wc)*(1.0-wi);

/* Which means we do not have to call linterpmute when we    */  
/* know the corresponding weight is zero. This may seem like */  
/* it will only save a small amount of CPU time but remember */  
/* that most situations have many cdps outside of the area   */  
/* that is completely surrounded by input mute locations.    */  
/* When outside the surrounded area, the binterpfind routine */  
/* has produced wi=0 or 1 and/or wc=0 or 1.                  */  

  if(aw != 0.) {
    linterpmute(offset,lwi.offs,lwi.tims,lwi.nto,mgtextr,&time);
    *timeout += aw*time;
  }
  if(bw != 0.) {
    linterpmute(offset,hii.offs,hii.tims,hii.nto,mgtextr,&time);
    *timeout += bw*time;
  }

  if(cw != 0.) {
    linterpmute(offset,lwc.offs,lwc.tims,lwc.nto,mgtextr,&time);
    *timeout += cw*time;
  }
  if(dw != 0.) {
    linterpmute(offset,hic.offs,hic.tims,hic.nto,mgtextr,&time);
    *timeout += dw*time;
  }

  return;
}

/* -----------------------------------------------------------         */
/* Specify compare function for qsort.                                 */

int compSort2 (const void * q1, const void * q2) {

  struct VelInfo* p1 = (struct VelInfo*) q1;
  struct VelInfo* p2 = (struct VelInfo*) q2;

/* Note I decided so sort to igc,igi order (kinf[2], then kinf[1])     */

  if(p1->kinf[2] < p2->kinf[2]) return (-1);
  if(p1->kinf[2] > p2->kinf[2]) return (1); 
  if(p1->kinf[1] < p2->kinf[1]) return (-1);
  if(p1->kinf[1] > p2->kinf[1]) return (1); 

  return (0); 

}

/* Standard binary search. But which side includes equal value */
/* is an important detail for other code in this program.      */

int bhighi(int *all, int last, int iguy) {
  int mid;
  int low = 0;
  int high = last;
  while (low < high) {
    mid = low + (high - low) / 2;
    if (iguy >= all[mid]) low = mid +1;
    else high = mid;
  }
  return low;
}

/* Linearly interpolate time at desired offset from arrays of offsets and times.      */
/*                                                                                    */
/* Input arguments:                                                                   */
/*  offset = the offset of the desire output time                                     */
/*  offs   = array of offsets (in increasing order)                                   */
/*  tims   = a time for each input offset                                             */
/*  nto    = number of offs and tims values (must be >0, can be =1)                   */
/*  mgtextr=0 no extrapolation at offset ends                                         */
/*         =1 extrapolate both lower and higher ends                                  */
/*         =2 extrapolate only at lower end                                           */
/*         =3 extrapolate only at higher                                              */
/*                                                                                    */
/* Output argument:                                                                   */
/*  time   = linearly interpolated time value for input offset value                  */

static void linterpmute(float offset,float *offs,float *tims,int nto,int mgtextr,float *time) {

  if(nto<2) {
    *time = tims[0];
    return;
  }

  int n = 1;
  if(offset <= offs[0]) {
    if(mgtextr==0 || mgtextr==3) {
      *time = tims[0];
      return;
    }
/*  needs to be n=1 here, but already initialed 1 (to avoid compiler warnings) */
  }
  else if(offset >= offs[nto-1]) {
    if(mgtextr==0 || mgtextr==2) {
      *time = tims[nto-1];
      return;
    }
    n = nto - 1;
  }
  else {
    n = bhighf(offs, nto, offset);
  }

  float wo = (offs[n]-offset) / (offs[n]-offs[n-1]);
  *time = wo*tims[n-1] + (1.0-wo)*tims[n];

  return;

}

int bhighf(float *all, int last, float iguy) {
  int mid;
  int low = 0;
  int high = last;
  while (low < high) {
    mid = low + (high - low) / 2;
    if (iguy >= all[mid]) low = mid +1;
    else high = mid;
  }
  return low;
}