/*
     AEDI                                                   Oli E. 02/12/2011
     
     ---------------------------------------------------------------------------
     Edge Refinement for A.aeolicus ferredoxin, written for the 9WL data from
     DESY. This program combines PREP_DATA and FPP_LSQ and includes reflection
     handling with STL routines.
     
     - Reads in |Fph| and Dano for all wavelengths from an MTZ file. Reflection
       are pre-sorted by CCP4 programs. Limits are determined and an internal
       list created.
       
     - Reads in coordinates from a *.pdb file. Fcalcs are determined and added
       to the list.
       
     - Global or local scaling between |Fph| and |Fcalc|
     
     - Determines f" and f' for each iron and all sulfurs.
     
     - Calculates quality indicators
     
     - Outputs result vectors and stats
     
     In contrast to older versions, this is supposed to more strictly adhere to
     C/C++ style, wrapping all functionality in functions and classes.
     ---------------------------------------------------------------------------

     06.06.2006  First integration of the old programs into a single routine.
                 SF calculation is the time-limiting step. Major redesign of 
		 code to use STL containers for reflections and coords.
		 
     14.06.2006  v0.1.03 Include data stats, |d"|/sigD seems to correlate with
                         edge quality and correlation. 
     
     15.06.2006  v0.1.04 Include calculation of variance, sd and R-factor during
                         global scaling.
			 
     08.09.2006  v0.1.05 Include calculation of overall e.s.d.s. Add Routine
                         EdgeDisp for calculation of f'.
			 
     09.09.2006  v0.1.06 Consider a cadmium ion in the structure; requires major
                         rearrangement for another anomalous species. 

     15.09.2006  v0.1.07 Dynamic handling of anomalous species through new
                         container structure for data objects. Corrected calcu-
			 lation of e.s.d.s.

     02.12.2011  v0.1.08 Fixed an error in ... ehem ... matrix & vector 
                         multiplication. 
*/


// -----------------------------------------------------------------------------

#include "xray.h"

// -----------------------------------------------------------------------------

const string ANO = "Fe";
const int METALS = 2;
const int WAVELS = 15;

const int waves[15] = { 7080, 7100, 7110, 7112, 7114, 7116, 7118, 7120, 7122, 7124, 7126, 7130, 8000, 8010, 8020 };


// -- Data structures ----------------------------------------------------------

typedef struct
{
   string element;       	// Columns XX from the PDB file
   xCoord r;		      	// Atom coordinates
   float  q;			// Occupancy
   float  B;			// isotropic Temperature factor
} Atom;
   
typedef struct
{
   double F, sigF;              // Amplitude and error for Fobs
   double D, sigD;              // Anomalous difference and error
} FDobs;

typedef struct 
{
    double A, B;                // Real and imaginary parts of an Fcalc
} ReIm;

typedef struct
{
   int h,k,l;
   vector<FDobs> wl;            // Container wl to have [WAVELS] elements
} HKLData;

typedef struct
{
   int h, k, l;
   vector<ReIm> type;           // Container type to have a column for each atom type and
} HKLFcalc;                     // each individual iron
       

// -- Function Prototypes ------------------------------------------------------

int ReadPDB(char*, vector<Atom>&);    
int ReadMTZ(char*, vector<HKLData>&);

int ReflStats(vector<HKLData>&);
int DefSymOps(vector<SymOp>&);
int SFCalc(vector<Atom>&, vector<HKLData>&, vector<HKLFcalc>&, vector<XFormFactor>&, XCrystal&); 
int ScaleGlobal(vector<HKLData>&, vector<HKLFcalc>&);
int EdgeAno(vector<HKLData>&, vector<HKLFcalc>&, vector<XFormFactor>&, XCrystal&);
int EdgeDisp(vector<HKLData>&, vector<HKLFcalc>&, vector<XFormFactor>&, XCrystal&);
/*
int ScaleLocal(vector<HKLData>&, vector<HKLFcalc>&);
int RFactor
*/


// -- M A I N () ---------------------------------------------------------------

int main(int argc, char* argv[])
{
   // -- Say 'Hi'

   cout.setf(ios::fixed, ios::floatfield);  
   cout.setf(ios::right, ios::adjustfield);  

   cout << endl << endl
        << "   ---«« A E D I    v0.1.07    06/06/06 »»---" << endl 
	<< endl;
   
   // -- Input Files: These are hardcoded for convenience, may shift to arguments.
   
   char *pdbfile = "C1_13000_refmac4.pdb";
   char *mtzfile = "CAD_combined2.mtz";

   
   // -- Global data containers:
   
   vector<Atom> structure;               // The structure
   vector<HKLData> data;                 // The data
   vector<HKLFcalc> fcalc;               // Fcalc Re and Im parts
   vector<XFormFactor> affs;             // Form factors (and thus list of atom types)
  
    
   // -- Read the Files. All details left to the routines. Only the filenames and the
   //    names of the containers (by reference) are passed.
   
  
   ReadPDB(pdbfile, structure);
   ReadMTZ(mtzfile, data);
   
   cout << endl << "=> Read Data Files:" << endl 
        << "   PDB: [" << pdbfile << "] with " << structure.size() << " atoms." << endl
        << "   MTZ: [" << mtzfile << "] with " << data.size() << " reflections." << endl
	<< endl;


   // -- Crystal Information is taken from the PDB File, not the MTZ!
   
   XCrystal xtal(pdbfile);                  


   // -- Some stats on the reflection file

   cout << endl << "=> Statistics on Edge Data: " << WAVELS << " Wavelengths." << endl; 
   
   ReflStats(data);

 
   // -- SF Calculation
   
   cout << endl << "=> Calculating Structure Factors .....";
   cout.flush();
   if ( SFCalc(structure, data, fcalc, affs, xtal) )
   {
      cout.flush();
      cerr << " ERROR during SF calculation!" << endl;
      exit(1);
   }
   cout << " done." << endl << endl; 
   
   // -- Global Scaling
   
   cout << endl << "=> Global Scaling of Fobs and Fcalc:" << endl;
   ScaleGlobal(data, fcalc);

   // -- Edge Calculation: f"
   EdgeAno(data, fcalc, affs, xtal);
 
   // -- Edge Calculation: f'
   EdgeDisp(data, fcalc, affs, xtal);
 
   // -- .. and that was that.
   
   cout << "   ---«« All done. Exiting. »»---" << endl << endl;
 
   return(0);
}

// -- End of main() ------------------------------------------------------------


int ReadPDB(char* filename, vector<Atom>& structure)
{
   // Reads all ATOM and HETATM lines from a PDB File and copies them into 
   // an atom vector (which is assumed empty).
   
   ifstream pdbfile(filename);
   if ( !pdbfile.is_open() )
   {
       cerr << "FILE ERROR: " << filename << " not found." << endl << endl;
       exit(1);
   }
   
   structure.clear();
   string line;
   Atom atom;
      
   while ( !pdbfile.eof() )
   {
      getline(pdbfile, line);
       
      if (  (line.compare(0, 6, "ATOM  ") == 0 ) 
         || (line.compare(0, 6, "HETATM") == 0 ) )
      {
         atom.element = line.substr(12, 2);                     // Get & Format Element
         if ( atom.element[0] == ' ' ) atom.element.replace(0, 1, "");
         else if ( atom.element[1] == ' ' ) atom.element.replace(1, 1, "");
         if ( atom.element.size() == 2 && isupper(atom.element[1]) ) 
	    atom.element[1] = tolower(atom.element[1]);
         istringstream co(line);
         co.ignore(30);
	 co >> atom.r.x >> atom.r.y >> atom.r.z >> atom.q >> atom.B;
	 
         structure.push_back(atom);
      }
   }
   
   pdbfile.close();
   return 0;
}


int ReadMTZ(char* filename, vector<HKLData>& data)
{
   // Reads data from an MTZ file that MUST contain nothing but four columns (F, sigF, D, sigD)
   // for WAVELS wavelengths!

   struct mtzHead { char mid[4]; int loc; };	        // The first bits of the MTZ file
   
   ifstream mtzfile(filename, ios_base::binary);
   if ( !mtzfile.is_open() ) 
   {
       cerr << "FILE ERROR: " << filename << " not found." << endl;
       exit(1);
   }
   
   // Check whether file is MTZ

   mtzHead header;
   mtzfile.read((char*) &header, sizeof(header));
   header.mid[3] = '\0';
   if ( strcmp(header.mid, "MTZ") )
   {
       cerr << "FILE ERROR: " << filename << " does not seem to be an MTZ file." << endl;
       exit(1);
   }
   
   // Read number of reflections (header is at end of file, so we better know before)
   
   int nref;
   mtzfile.seekg( ( header.loc - 1 ) * sizeof(int) + 160 + 15);  // Where NREF should be ...
   mtzfile >> nref;
   
   
   // Now read data
   
   struct InFile{ float h, k, l, F[WAVELS * 4]; };
   InFile line;
   FDobs wave;
   mtzfile.seekg(80);         // Data start here
   
   for ( int i = 0; i < nref; i++ )
   {
       HKLData refl;
       mtzfile.read((char*) &line, sizeof(line));

       // Replace missing number flags with zeros. Is this legal?!? 
       for (int j = 0; j < ( WAVELS * 4); j++)
           if ( isnan(line.F[j]) ) line.F[j] = 0;

       // READ data and copy to the final table
       refl.h = (int)line.h; refl.k = (int)line.k; refl.l = (int)line.l;
       for ( int j = 0; j < WAVELS; j++ )
       {
           wave.F    = line.F[ (4 * j) ];
           wave.sigF = line.F[ (4 * j) + 1 ];
           wave.D    = line.F[ (4 * j) + 2 ];
           wave.sigD = line.F[ (4 * j) + 3 ];
           refl.wl.push_back(wave);
       }
       data.push_back(refl);
   }

   mtzfile.close();
   return 0;
}


int ReflStats(vector<HKLData>& data)
{
   // Calculate various useful (?) stats on the Fobs data. Hope these do reflect the
   // inconsistencies observed in the DESY data
   //
   // Ouput is:
   // Wavelength | NRefl | Av.F | Av.sigF | rmssigF | NAno | AnoPos | Av.D | Av.sigD | rmsD |
   
   // No.Refl is counted individually, skipping the 0's!
   
   vector<HKLData>::iterator  it_d;			

   vector<int> refls, anom;
   vector<double> avgf, avgsigf, fsigf, anopos, avgd, avgsigd, dovsigd, rmsD;
   
   for ( int i = 0; i < WAVELS; i++ )
   {
      int r = 0, a = 0, anop = 0;
      double ftot = 0, sigftot = 0,  
             dtot = 0, sigdtot = 0, absdtot = 0;
	     
      for ( it_d = data.begin(); it_d != data.end(); it_d++ )
      {
         if ( it_d->wl[i].F != 0 )
	 { 
	    r++;
            ftot    += it_d->wl[i].F;
            sigftot += it_d->wl[i].sigF;
	    if ( it_d->wl[i].D != 0)
	    { 
	       a++;
               dtot    += it_d->wl[i].D;
               absdtot += fabs( it_d->wl[i].D );
               sigdtot += it_d->wl[i].sigD;
	       if ( it_d->wl[i].D > 0 ) anop++;
	    }
	 }
      }   
      refls.push_back(r);
      avgf.push_back( ftot / r );
      avgsigf.push_back( sigftot / r );
      fsigf.push_back( ftot / sigftot );
      anom.push_back(a);
      anopos.push_back( anop / (double)a );
      avgd.push_back( dtot / a );
      avgsigd.push_back( sigdtot / a);
      dovsigd.push_back( absdtot / sigdtot);
   }
   
   cout << endl;
   cout << "   Energy ---- NRefl --- <F> - <sigF> - F/sigF --- NAno -- Ano>0 - <Ano> - <sigD> - |d\"|/sigD" << endl
        << "   ------------------------------------------------------------------------------------------" << endl;

   cout << setprecision(2);
   
   for ( int i = 0; i < WAVELS; i++ )
   {
      cout << "   " << waves[i] << " eV:"
           << setw(9) << refls[i]   
	   << setw(9) << avgf[i] 
	   << setw(8) << avgsigf[i] 
	   << setw(9) << fsigf[i]
           << setw(10) << anom[i]   
           << setw(7) << anopos[i]   
	   << setw(9) << avgd[i] 
	   << setw(9) << avgsigd[i] 
	   << setw(10) << dovsigd[i]
	   << endl;
   }	
   cout << "   ------------------------------------------------------------------------------------------" << endl
        << endl;
      
   return 0;
}


int DefSymOps(vector<SymOp>& symops)
{
   // Reading Symops for P21 from a file. I want all this hardwired. All 230! NOW!!
   
   ifstream symfile("p22121.symop");
   string line;
   getline(symfile, line);              // Ignore first line: SPG Name and Number.
   
   for (int i = 0; i < 4; i++)
   {
      getline(symfile, line);
      istringstream op(line);
      for (int j = 0; j < 3; j++)
         for (int k = 0; k < 3; k++)
	    op >> symops[i].r[j][k];
      op >> symops[i].t[0] >> symops[i].t[1] >> symops[i].t[2];
   }   
   symfile.close();      
   return 0;
}


int SFCalc(vector<Atom>& structure, vector<HKLData>& data, 
           vector<HKLFcalc>& fcalc, vector<XFormFactor>& affs, XCrystal& xtal)
{

   // 13.09.2006
   // Recognizing that the correlation coefficient of the edge routine strongly depends on a 
   // Complete set of anomalous scatterers, f_calcs should be determined for all atom types.
   // Only one of those will have an edge that justifies individual refinement, for all others a
   // total contribution will be calculated.
   
   // This requires a more dynamic fcalc structure, as the number of atom types is not necessarily
   // known beforehand and is only derived here.
      
   int yop;
   vector<Atom>::iterator        it_s;			
   vector<HKLData>::iterator     it_d;			
   vector<HKLFcalc>::iterator    it_c;			
   vector<XFormFactor>::iterator it_a;

   // Scan structure for atom types, read all required form factors into a 
   // container. Requires formatted element names as done in ReadPDB.
   // Make sure that iron is last in the list at this point.

   int iron = 0;
   for ( it_s = structure.begin(); it_s != structure.end(); it_s++)
   {
      if ( it_s->element == "Fe" ) iron = 1;
      else
      {
         yop = 0;
         for (it_a = affs.begin(); it_a != affs.end(); it_a++)
            if ( it_s->element == it_a->ID().substr(0, it_s->element.size()) ) yop = 1;
         if ( yop == 0 )
	     affs.push_back(XFormFactor(it_s->element));      // element at a given index.
      }
   }
   if ( iron ) affs.push_back(XFormFactor("Fe"));              // If iron was found, add as last.
   
   // The structure contains found.size() types of atoms, including METALS irons, which need to be
   // further separated, such that the type-container in fcalc will have 
   // found.size() -1 + METALS elements. 
   
   // Initialize symmetry operators

   vector<SymOp> symops(4);              // Symmetries, hardwired to P21 for now
   DefSymOps(symops);                    // Read file with ops.
   
   
   // Copy the h,k,l columns to the fcalc container and flush all fields.
   
   ReIm ttype;
   ttype.A = ttype.B = 0;
   
   for ( it_d = data.begin(); it_d != data.end(); it_d++)
   {
      HKLFcalc temp;
      temp.h = it_d->h; temp.k = it_d->k; temp.l = it_d->l;
      for (int i = 0; i < (affs.size() - 1 + METALS); i++)
	 temp.type.push_back(ttype);
      fcalc.push_back(temp);
   }
   
   // Structure factor calulation. Separate for each atom type. For now, iron is the edge element,
   // such that irons are separated as individual atoms.
   
   int col, fe_counter = 0;
   XFormFactor f_at;
   xCoord site;
   vector<xCoord> sites;                         // To hold the symmetry equvalents
   vector<xCoord>::iterator it_sit; 
   
   for ( it_s = structure.begin(); it_s != structure.end(); it_s++)
   {
      // Find correct form factor and column in fcalc container. Iron atoms will be last.
      
      for ( it_a = affs.begin(), col = 0; it_a != affs.end(); it_a++, col++)
      {
	 if ( it_s->element == it_a->ID().substr(0, it_s->element.size()) )
         {
	    f_at = *it_a;
	    break;
	 }
      }
      if ( it_s->element == "Fe" )
         col = affs.size() - 1 + fe_counter++;  // Adjust column for Fe
      
      // Fractionalize and create symmetry equivalents (unity operator is also processed!)
      
      sites.clear();
      site.x = it_s->r.x; site.y = it_s->r.y; site.z = it_s->r.z;
      site   = xtal.Deorthogonalize(site);
      xCoord temp;
      
      vector<SymOp>::iterator it_sym;
      for ( it_sym = symops.begin(); it_sym != symops.end(); it_sym++)
      {
         temp.x = it_sym->t[0] + ( it_sym->r[0][0] * site.x + it_sym->r[0][1] * site.y + it_sym->r[0][2] * site.z );    
         temp.y = it_sym->t[1] + ( it_sym->r[1][0] * site.x + it_sym->r[1][1] * site.y + it_sym->r[1][2] * site.z );    
         temp.z = it_sym->t[2] + ( it_sym->r[2][0] * site.x + it_sym->r[2][1] * site.y + it_sym->r[2][2] * site.z );
 	 sites.push_back(temp);    
      }

      // -- Go through all reflections, calculate Fcalc and sort it all.
      
      double dmin, aff;
      
      for ( it_c = fcalc.begin(); it_c != fcalc.end(); it_c++)
      {
         dmin = xtal.GetResolution(it_c->h, it_c->k, it_c->l);
         aff = f_at.f_d(dmin, it_s->B) * it_s->q;     // Is this really the place for q?
	 for ( it_sit = sites.begin(); it_sit != sites.end(); it_sit++)
	 {
	    it_c->type[col].A += aff * cos(2*PI * ( it_c->h * it_sit->x
	                                          + it_c->k * it_sit->y
	                                          + it_c->l * it_sit->z ));
	    it_c->type[col].B += aff * sin(2*PI * ( it_c->h * it_sit->x
	                                          + it_c->k * it_sit->y
	                                          + it_c->l * it_sit->z ));
	 }
      }
   }
   return 0;
}


int ScaleGlobal(vector<HKLData>& data, vector<HKLFcalc>& fcalc)
{
   // GLOBAL SCALING:

   // Original global scaling is rather coarse and probably one of the main reasons
   // for the "vertical" divergence of individual iron curves.
   // Local scaling may remedy this. For every individual reflection a scaling factor
   // is calculated from a cube surrounding this reflection in reciprocal space.
   // Variation of the cube size will be important. 
   //
   // Issue: How can we deal with missing reflections?? For now -> zero. Likely distorts 
   //        high resolution range! Fill-in may be better.
   
   // -- Calculate scaling factors between Fc_total and the Fobs
   
   double scale[WAVELS];
   vector<HKLData>::iterator  it_d;			
   vector<HKLFcalc>::iterator it_c;			
   vector<ReIm>::iterator     it_r;			

   vector<double> fc_total;
   vector<double> ti;
   double t1, t2;

   for (int i = 0; i < WAVELS; i++)
   {
      t1 = t2 = 0;
      double fctot = 0;        
      for ( it_d = data.begin(), it_c = fcalc.begin(); it_d != data.end(); it_d++, it_c++)
      {
         double atot = 0, btot = 0;           
	 for ( it_r = it_c->type.begin(); it_r != it_c->type.end(); it_r++)
	 {
	    atot += it_r->A;
	    btot += it_r->B;
	 }
         fctot = sqrt( pow( atot, 2 ) + pow( btot, 2 ) );
	 t1 += it_d->wl[i].F * fctot;  			 
	 t2 += pow( it_d->wl[i].F, 2 );
	 fc_total.push_back(fctot);          // Save Fc_total for statistics
	 ti.push_back(t2);                
      }
      scale[i] = t1 / t2;
   }

   // -- Scale Fobs and Dano with these factors, calculate variance and R-factor

   vector<double>::iterator it_fc;
   vector<double> var, rfac;
   
   for (int i = 0; i < WAVELS; i++)
   {
      double d = 0, r = 0, fosum = 0;
      for ( it_d = data.begin(), it_fc = fc_total.begin(); it_d != data.end(); it_d++, it_fc++)
      {
         it_d->wl[i].F    *= scale[i]; 
         it_d->wl[i].sigF *= scale[i]; 
	 it_d->wl[i].D    *= scale[i];   
         it_d->wl[i].sigD *= scale[i];
	 d += pow( it_d->wl[i].F - *it_fc, 2 );
	 r += fabs( it_d->wl[i].F - *it_fc );
	 fosum += it_d->wl[i].F; 
      }
      var.push_back( ( d / data.size() ) / ti[i]);
      rfac.push_back( r / fosum );
   }

   cout << setprecision(4);
   cout << endl
        << "   Energy ----- Scale --- Var ----- SD ------ Rfac -" << endl
	<< "   -------------------------------------------------" << endl;
   for (int i = 0; i < WAVELS; i++)
      cout << "   " << waves[i] << " eV: " 
           << setw(10) << scale[i] 
	   << setw(10) << var[i]
	   << setw(10) << sqrt(var[i])
	   << setw(10) << rfac[i] << endl;   
   cout << "   -------------------------------------------------" << endl
        << endl;
   cout.flush();
   
   return 0;
}


int EdgeAno(vector<HKLData>& data, vector<HKLFcalc>& fcalc, 
            vector<XFormFactor>& affs, XCrystal& xtal)
{
   // Prepares the design matrix for solving the linear equation system and
   // calls the GaussJordan solver.
   
   int n = (affs.size() - 1 + METALS);       // Number of rows in (At A), anomalous species 
   int m = WAVELS;                              // Number of columns in (At F) 
   int counter = 0;

   Matrix<double> ata(0.0,n,n);
   Matrix<double> atf(0.0,n,m);
   vector<double> arw(n);

   Matrix<double> a_full(0.0,data.size(),n);    // The full design matix 
   
   vector<HKLData>::iterator     it_d;
   vector<HKLFcalc>::iterator    it_c;
   vector<ReIm>::iterator        it_r;
   vector<XFormFactor>::iterator it_a;

   int i = 0;                                   
   
   cout << endl << endl
        << "   ---------------------------------------------------------" << endl
	<< "=> ---       REFINEMENT OF ANOMALOUS CONTRIBUTIONS       ---" << endl
        << "   ---------------------------------------------------------" << endl
        << endl;	
	
   cout << "=> The structure contains " << affs.size() << " atom types." << endl
        << "=> " << METALS << " atoms of type " << ANO << " will be refined individually." << endl
	<< endl;  
   
   for ( it_d = data.begin(), it_c = fcalc.begin(); it_d != data.end(); it_d++, it_c++)
   {
      int check = 1;    // For data cut-offs ...
      if (check)
      {  
      	  // -- Calculate the combined phase angle for this (hkl) (in radians)
    	  
          double at = 0, bt=0;
	  for ( it_r = it_c->type.begin(); it_r != it_c->type.end(); it_r++ )
	  {
	     at += it_r->A;
	     bt += it_r->B;
	  }
	  double phiprot = atan2 (bt, at);
 
          double d = xtal.GetResolution(it_c->h, it_c->k, it_c->l);
	   
   	  // -- Now calculate the ith row of the design matrix A

          int l = 0;
	  for ( it_r = (it_c->type.begin()), it_a = (affs.begin()); it_r != it_c->type.end(); it_r++ )
          {
	     arw[l++] = -2 / it_a->f_d( d ) * ( it_r->B * cos(phiprot)
	                                      - it_r->A * sin(phiprot));
 	     
	     if ( it_a->ID().substr(0,2) != ANO ) it_a++;
	  }
	  
	  // -- create the full design matrix 

	  for (int ii = 0; ii < n; ii++)
  	     a_full[i][ii] = arw[ii];

          // -- Add these values to (At A)
       
          for (int j = 0; j < n; j++)
              for (int k = 0; k < n; k++)  ata[j][k] += arw[j] * arw[k];

          // -- And to (At F) for all wavelengths
       
   	  for (int j = 0; j < n; j++)
     	  {
     	     for (int k = 0; k < m; k++)
	        atf[j][k] += arw[j] * it_d->wl[k].D;
     	  }
       counter ++;
       i++;
       }
   }
   cout << endl << "=> " << counter << " reflections out of " << data.size() << " will be used in least squares." 
        << endl << endl;
   
   // -- Now that the Matrices Are set up, solve the linear equation system
   //    (At A) * x = (At Fs) by Gauss-Jordan Elimination

   GaussJordan(ata, atf);
   
   // -- This leaves ata as (At A)^(-1) and the atf's as the solution 


   // -- Calculate correlation coefficient between observed and expected
   //    Differences.
   
   // For calculating the expected differences we need the complete
   // design matrix, a_full.
   
   vector<double> dcalc, dobs, sum;
   
   cout << endl << "=> Correlations for Least Squares of Anomalous Differences:" << endl << endl;
   
   cout << "   Energy ---- CC ---" << endl
        << "   ------------------" << endl;
   
   for (int i = 0; i < m; i++)
   {
      int j = 0;
      double s = 0;
      dcalc.clear(); dobs.clear();
      for ( it_d = data.begin(); it_d != data.end(); it_d++)
      {
         dobs.push_back( it_d->wl[i].D );
	 double clc = 0;
	 for(int k = 0; k < n; k++)
	    clc += a_full[j][k] * atf[k][i]; 
	 dcalc.push_back(clc); 
	 s += pow(it_d->wl[i].D - clc, 2);    
         j++;
      }
      sum.push_back(sqrt(s));
      cout << "   " << waves[i] << " eV:" << setw(10) << PearsonsR(dobs, dcalc) << endl ;
   }
   cout << "   ------------------" << endl << endl;


   // -- Calculate the covariance matrix 

   cout << endl << "=> Covariance Matrix:" << endl << endl << "   - ";
   
   // -- Write table header with atom names

   for ( it_a = affs.begin(); it_a != (affs.end() - 1); it_a++)
      cout << it_a->ID() << "  ----- ";
   for (int m = 0; m < METALS; m++)
      cout << ANO << (m+1) << " ---- ";   
   cout << endl <<  "   ";
   for (int m = 0; m < n; m++)
      cout << "--------";   
   cout << "-------" << endl;

   for (int i = 0; i < n; i++)
   {
      for (int j = 0; j < n; j++)
	 cout << setw(9) << ata[i][j] / sqrt(ata[i][i] * ata[j][j]);
      cout << endl;
   }
   cout <<  "   ";

   for (int m = 0; m < n; m++)
      cout << "--------";   
   cout << "-------" << endl;

   // -- output result by wavelength
   
   cout << setprecision(4);

   cout << endl << endl << "=> Refined f\" values for all sulfurs and individual metals:" << endl << endl;
   cout << "   Energy ---- Solution vectors ------------->" << endl
        << "   -[eV]-------  ";
	    
    // -- Write table header with atom names

   for ( it_a = affs.begin(); it_a != (affs.end() - 1); it_a++)
      cout << it_a->ID() << "  -----  ";
   for (int m = 0; m < METALS; m++)
      cout << ANO << (m+1) << " ----- ";   
   cout << endl <<  "   ";
   for (int m = 0; m < n; m++)
      cout << "-----------";   
   cout << endl;

   for (int i = 0; i < m; i++)
   {
       cout << "    " << waves[i] << "   ";
       for (int j = 0; j < n; j++) cout << setw(10) << atf[j][i] ;  
       cout << endl;
   }
   cout <<  "   ";
   for (int m = 0; m < n; m++)
      cout << "-----------";   
   cout << endl;

   // Output e.s.d.s: 
   // Sqrt of diagonal elements of (ATA-1) * sqrt( sum(dobs-dcalc)^2) / (Nref-nparam) )

   cout << "   e.s.d.  " << endl;
   for (int i = 0; i < m; i++)
   {
      cout << "           ";
      for (int j = 0; j < n; j++)
          cout << setw(10) << sqrt( sum[i] / (data.size() - n) ) * sqrt(ata[j][j]);
 
      cout << endl;
   }
   cout <<  "   ";
   for (int m = 0; m < n; m++)
      cout << "-----------";   
   cout << endl;

   // --------------------------------------------------------------------

   return 0;
}


int EdgeDisp(vector<HKLData>& data, vector<HKLFcalc>& fcalc,
             vector<XFormFactor>& affs, XCrystal& xtal)
{
   // Prepares the design matrix for solving the linear equation system and
   // calls the GaussJordan solver for dispersive differences
   
   // Dispersive differences are calculated against the total F_calcs. This is
   // probably the most critical parameter in this and needs to be optimized.
   // All stats are repeated code from EdgeAno(), so a separate routine to analyze
   // the results will eventually be desirable. Maybe same for setting up ATA.
   
   
   int n = (affs.size() - 1 + METALS);   // Number of rows in (At A), anomalous species 
   int m = WAVELS;                       // Number of columns in (At F) 
   int counter = 0;

   Matrix<double> ata(0.0,n,n);
   Matrix<double> atf(0.0,n,m);
   vector<double> arw(n);
   vector<double> fc_total;

   Matrix<double> a_full(0.0,data.size(),n);    // The full design matix 
   
   vector<HKLData>::iterator     it_d;
   vector<HKLFcalc>::iterator    it_c;
   vector<ReIm>::iterator        it_r;
   vector<XFormFactor>::iterator it_a;

   int i = 0;                                                 // Cutoff check flag

   cout << endl << endl
        << "   ---------------------------------------------------------" << endl
	<< "=> ---       REFINEMENT OF DISPERSIVE CONTRIBUTIONS      ---" << endl
        << "   ---------------------------------------------------------" << endl
        << endl;	
	
   
   for ( it_d = data.begin(), it_c = fcalc.begin(); it_d != data.end(); it_d++, it_c++)
   {
      int check = 1;    // For data cut-offs ...
      if (check)
      {  
	  // -- Calculate the combined phase angle for this (hkl) (in radians) and
	  //    Fc_total for the dispersive difference
    	  
    	  
          double at = 0, bt=0;
	  for ( it_r = it_c->type.begin(); it_r != it_c->type.end(); it_r++ )
	  {
	     at += it_r->A;
	     bt += it_r->B;
	  }
	  double phiprot = atan2 (bt, at);
          double fctot = sqrt( pow( at, 2 ) + pow( bt, 2 ) );
	  fc_total.push_back(fctot);          // Save Fc_total for statistics
 
          // Calculate f', the real part of anomalous scattering
           
          double dmin = xtal.GetResolution(it_c->h, it_c->k, it_c->l);
	   
   	  // -- Now calculate the ith row of the design matrix A

          int l = 0;
	  for ( it_r = it_c->type.begin(), it_a = affs.begin(); it_r != it_c->type.end(); it_r++)
          {
	     arw[l++] =  1 / it_a->f_d( dmin ) * ( it_r->A * cos(phiprot)
	                                         + it_r->B * sin(phiprot));
 	     
	     if ( it_a->ID().substr(0,2) != ANO ) it_a++;
	  }

	  // -- create the full design matrix 

	  for (int ii = 0; ii < n; ii++)
  	     a_full[i][ii] = arw[ii];

          // -- Add these values to (At A) for each metal and sulfur
       
          for (int j = 0; j < n; j++)
              for (int k = 0; k < n; k++)  ata[j][k] += arw[j] * arw[k];

          // -- And to (At F) for all wavelengths
       
   	  for (int j = 0; j < n; j++)
     	  {
     	     for (int k = 0; k < m; k++)
	        atf[j][k] += arw[j] * (it_d->wl[k].F - fctot);
     	  }
       counter ++;
       i++;
       }
   }
   cout << endl << "=> " << counter << " reflections out of " << data.size() << " will be used in least squares." 
        << endl << endl;
   
   // -- Now that the Matrices Are set up, solve the linear equation system
   //    (At A) * x = (At Fs) by Gauss-Jordan Elimination

   GaussJordan(ata, atf);
   
   // -- This leaves ata as (At A)^(-1) and the atf's as the solution 


   // -- Calculate correlation coefficient between observed and expected
   //    Differences.
   
   // For calculating the expected differences we need the complete
   // design matrix, a_full.
   
   vector<double> dcalc, dobs, sum;
   vector<double>::iterator it_fc;
   
   cout << endl << "=> Correlations for Least Squares of Dispersive Differences:" << endl << endl;
   
   cout << "   Energy ---- CC ---" << endl
        << "   ------------------" << endl;
   
   for(int i = 0; i < m; i++)
   {
      int j = 0;
      double s = 0;
      dcalc.clear(); dobs.clear();
      for ( it_d = data.begin(), it_fc = fc_total.begin(); it_d != data.end(); it_d++, it_fc++)
      {
         dobs.push_back( it_d->wl[i].F - *it_fc );
	 double clc = 0;
	 for(int k = 0; k < n; k++)
	    clc += a_full[j][k] * atf[k][i]; 
	 dcalc.push_back(clc);      
 	 s += pow(it_d->wl[i].F - clc, 2);    
        j++;
      }
      sum.push_back(sqrt(s));
      cout << "   " << waves[i] << " eV:" << setw(10) << PearsonsR(dobs, dcalc) << endl ;
   }
   cout << "   ------------------" << endl << endl;

   // -- output result by wavelength
   
   cout << setprecision(4);

   cout << endl << endl << "=> Refined f' values for all sulfurs and individual metals:" << endl << endl;
   cout << "   Energy ---- Solution vectors ------------->" << endl
        << "   -[eV]-------  ";
	    
    // -- Write table header with atom names

   for ( it_a = affs.begin(); it_a != (affs.end() - 1); it_a++)
      cout << it_a->ID() << " ----  ";
   for (int m = 0; m < METALS; m++)
      cout << ANO << (m+1) << " ----- ";   
   cout << endl <<  "   ";
   for (int m = 0; m < n; m++)
      cout << "-----------";   
   cout << endl;

   for (int i = 0; i < m; i++)
   {
       cout << "    " << waves[i] << "   ";
       for (int j = 0; j < n; j++) cout << setw(10) << atf[j][i] ;  
       cout << endl;
   }
   cout <<  "   ";
   for (int m = 0; m < n; m++)
      cout << "-----------";   
   cout << endl;


   // Output e.s.d.s: Sqrt of diagonal elements of (ATA-1)

   cout << "   e.s.d.  " << endl;
   for (int i = 0; i < m; i++)
   {
      cout << "           ";
      for (int j = 0; j < n; j++)
          cout << setw(10) << sqrt( sum[i] / (data.size() - n) ) * sqrt(ata[j][j]);
 
      cout << endl;
   }
   cout <<  "   ";
   for (int m = 0; m < n; m++)
      cout << "-----------";   
   cout << endl << endl;
 
   // --------------------------------------------------------------------

   return 0;
}

