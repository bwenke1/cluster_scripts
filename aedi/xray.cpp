/*


    XRAY.CPP                                        Oli E. 16.10.2002
    
    Methods for Xray classes
    
    16.10.2002 - Added class XCrystal member functions:
                 - XCrystal()
		 - GetResolution()
		 - GetHKL()
		 - GetCellVolume()
		 - GetMatthews()
    17.10.2002 - Added funct FileFormFactor(), FormFactor()
    18.10.2002 - Removed FormFactor functions and xglobals.h
                 Added class XFormFactor member functions:
		 - XFormFactor()
		 - f()
		 - f_d()
    20.10.2002 - Added template class Vector member functions
               - Added template class Matrix member functions
               - Added template class Matrix3D member functions
    22.10.2002 - Added function GaussJordan()
               - Added function PearsonsR()
    29.04.2006 - Added XCrystal::XCrystal(char *filename)
    16.09.2006 - Changed XFormFactor constructor to use AFFs from xglobals.h      
      
*/

#include "xray.h"
#include "xglobals.h"

using namespace std;

// -----------------------------------------------------------------------------
// -- Global functions ---------------------------------------------------------
// -----------------------------------------------------------------------------

void GaussJordan(Matrix<double> &A, Matrix<double> &b)
{
   // -- Solve the linear equation System Ax = b by Gauss-Jordan elimination
   //    with full pivoting.
   //    A contains the design matrix, the columns of b are data vectors.
   //    This routine destroys both input arrays (!). After execution, A will
   //    contain its inverse (A^(-1)), and b will contain the solution vectors.
   //
   //    -- from Numerical recipes in C++, p.42, modified to use STL vectors
    
   int n = A.nrows();
   int m = b.ncols();
   int icol, irow;
   double big, dum, pivinv;
   
   vector<int> indxc(n), indxr(n), ipiv(n);      // Working arrays for pivoting
   ipiv.assign(n, 0);                            // flush main pivoting array
   
   for (int i = 0; i < n; i++)                   // Main loop
   {
       big = 0;
       for (int j = 0; j < n; j++)               // Search for pivot element
       {
           if (ipiv[j] != 1)
	       for (int k = 0; k < n; k++)
	       {
	           if (ipiv[k] == 0)
		   {
		       if (fabs(A[j][k]) >= big)
		       {
		           big = fabs(A[j][k]);
			   irow = j;             // Position of the pivot
			   icol = k;
		       }
		   }
	       }
       }
       ++(ipiv[icol]);
       
       // -- reindex rows to put the pivot element on the diagonal
       
       if (irow != icol)              // i.e. the pivot is not on the diagonal
       {
           for (int l = 0; l < n; l++) SWAP(A[irow][l], A[icol][l]);
           for (int l = 0; l < m; l++) SWAP(b[irow][l], b[icol][l]);
       }
       
       indxr[i] = irow;
       indxc[i] = icol;
       if (A[icol][icol] == 0)
       {
           cerr << "ERROR: Gauss_Jordan: Singular matrix!\n";
	   exit(1);
       }
       
       pivinv = 1 / A[icol][icol];          // (pivot row) / (pivot element)
       A[icol][icol] = 1;
       for (int l = 0; l < n; l++) A[icol][l] *= pivinv;
       for (int l = 0; l < m; l++) b[icol][l] *= pivinv;
       
       for (int ll = 0; ll < n; ll++)         // reduce rows ...
           if (ll != icol)                    // ... except for pivot row
	   {
	       dum = A[ll][icol];
	       A[ll][icol] = 0;
               for (int l = 0; l < n; l++) A[ll][l] -= A[icol][l] * dum;
               for (int l = 0; l < m; l++) b[ll][l] -= b[icol][l] * dum;
	   } 
   }

   // -- unscramble the solution by reindexing the columns
   
   for (int l = n-1; l >= 0; l--)
   {
      if (indxr[l] != indxc[l])
         for (int k = 0; k < n; k++) SWAP (A[k][indxr[l]], A[k][indxc[l]]);
   }
}

// -----------------------------------------------------------------------------

double PearsonsR(vector<double> &x, vector<double> &y)
{
   // -- Calculate the linear correlation coefficient (Pearson's r) between
   //    arrays x and y. Prob is the significance level for disproving the 
   //    null hypothesis of zero correlation: small values indicate significant
   //    Correlation.
   //
   //    This is simplified from Numerical Recipes in C++. p.643. Modified to
   //    return r and skip the calculation of prob, Student's t and
   //    Fischer's z. Regards to Marliese.
   
   const double TINY = 1.0e-30;
   double yt, xt, r;
   double syy = 0.0, sxy = 0.0, sxx = 0.0, ay = 0.0, ax = 0.0;
   
   int n = x.size();
   for(int j = 0; j < n; j++)                        // Calculate the mean
   {
      ax += x[j];
      ay += y[j];
   }
   ax /= n;
   ay /= n;
   for(int j = 0; j < n; j++)                        // Calculate r
   {
      xt = x[j] - ax;
      yt = y[j] - ay;
      sxx += xt * xt;
      syy += yt * yt;
      sxy += xt * yt;
   }
   r = sxy / (sqrt(sxx * syy) + TINY);
   return(r);
}



// >===============--------------
// -----------------------------------------------------------------------------



// -----------------------------------------------------------------------------
// -- CLASS -- XCrystal --------------------------------------------------------
// -----------------------------------------------------------------------------

// -- Private:

int    XCrystal::GetAsymUnits()
{
    // This needs to be worked on, ehem.
    return(1);    
}


// -- Public:

XCrystal::XCrystal()
{
    a = b = c = 0;
    alpha = beta = gamma = 90;
    spg_num = 1;
    
    asym_units = GetAsymUnits();
}


XCrystal::XCrystal(double oa, double ob, double oc, 
                   double oalpha, double obeta, double ogamma,
		   int ospg)
{
    a = oa;
    b = ob;
    c = oc;
    alpha = oalpha;
    beta  = obeta;
    gamma = ogamma;
    spg_num = ospg;

    asym_units = GetAsymUnits();
}


XCrystal::XCrystal(char *filename)
{
    // Alternative Constructor from a PDB or MTZ File. So far Spacegroup is not used.
    // Other file types are not supported.
   
   a = 0;
   
   ifstream file(filename, ios_base::binary);
   if ( !file.is_open() )
   {
       cerr << "FILE ERROR: " << filename << " not found." << endl << endl;
       exit(1);
   }
   
   char mid[4];
   
   file.read((char*) &mid, sizeof(mid));
   mid[3] = '\0';
   if ( strcmp(mid, "MTZ") == 0 )                   // Is it an MTZ file?         
   {                                      
       int loc;
       file >> loc; 
       file.seekg( ( loc - 1 ) * sizeof(int) + 244 );
       file >> a >> b >> c >> alpha >> beta >> gamma;
   } 
   else                                             // Then it better be a PDB file ...
   {
       file.seekg(0);
       char *lineID = new char[7]; lineID[6] = '\0';
       char *dummy  = new char[132];
   
       while ( file >> lineID )
       {
           if (    strcmp(lineID, "CRYST1") == 0
	        || strcmp(lineID, "CRYST" ) == 0 )
           {
              file >> a >> b >> c >> alpha >> beta >> gamma;
	      break;
           }
       file.getline(dummy, 132);
       }
    }

    if ( a == 0 )
    {
        cerr << "FILE ERROR: " << filename << "is neither MTZ nor PDB with CRYST1 card." 
             << endl << endl;
             exit(1);
    }
    file.close();
    
    // A dummy. Space groups are not yet in.
    spg_num = 1;
    
    asym_units = GetAsymUnits();
}


double XCrystal::GetResolution(int oH, int oK, int oL)
{
    // If it's fine for triclinic ...

    double ssq, deta, shh, skk, sll, shk, shl, skl, dmin;
    double ca = cos(alpha * PI / 180);
    double cb = cos(beta  * PI / 180);
    double cg = cos(gamma * PI / 180);

    deta = 4 * (1 - pow(ca,2) - pow(cb,2) - pow(cg,2) + 2 * ca * cb * cg);
    shh  = (1 - ca * ca) / (deta * a * a);
    skk  = (1 - cb * cb) / (deta * b * b);
    sll  = (1 - cg * cg) / (deta * c * c);
    shk  = 2 * (ca * cb - cg) / (deta * a * b);
    shl  = 2 * (ca * cg - cb) / (deta * a * c);
    skl  = 2 * (cb * cg - ca) / (deta * b * c);

    ssq = shh * oH * oH + skk * oK * oK + sll * oL * oL
        + shk * oH * oK + shl * oH * oL + skl * oK * oL;

    dmin = sqrt(1 / (4 *ssq));
    
    return(dmin);   
}


double XCrystal::GetResolution(xHKL oHKL)
{

    double ssq, deta, shh, skk, sll, shk, shl, skl, dmin;
    double ca = cos(alpha * PI / 180);
    double cb = cos(beta  * PI / 180);
    double cg = cos(gamma * PI / 180);

    deta = 4 * (1 - pow(ca,2) - pow(cb,2) - pow(cg,2) + 2 * ca * cb * cg);
    shh  = (1 - ca * ca) / (deta * a * a);
    skk  = (1 - cb * cb) / (deta * b * b);
    sll  = (1 - cg * cg) / (deta * c * c);
    shk  = 2 * (ca * cb - cg) / (deta * a * b);
    shl  = 2 * (ca * cg - cb) / (deta * a * c);
    skl  = 2 * (cb * cg - ca) / (deta * b * c);

    ssq = shh * oHKL.H * oHKL.H + skk * oHKL.K * oHKL.K  
        + sll * oHKL.L * oHKL.L + shk * oHKL.H * oHKL.K 
	+ shl * oHKL.H * oHKL.L + skl * oHKL.K * oHKL.L;

    dmin = sqrt(1 / (4 *ssq));
    
    return(dmin);   
}


double XCrystal::GetCellVolume()
{
    double v;
    double ca = cos(alpha * PI / 180);
    double cb = cos(beta  * PI / 180);
    double cg = cos(gamma * PI / 180);
    
    v  = a * b * c * sqrt(  1 - pow(ca,2) - pow(cb,2) - pow(cg,2) 
                          + 2 * ca * cb * cg );
    return(v);
}


double XCrystal::GetMatthews(double omolweight, int on)
{
    double vm, asym, v;
    
    v  = GetCellVolume();

    if ((omolweight != 0) && (on != 0))
        return(v / (omolweight * on * asym_units));
    else
        return(0);
}


xCoord XCrystal::Orthogonalize(double ox, double oy, double oz)
{
    // Transform fractional into cartesian coordinates using the PDB
    // orthogonalization convention:
    
    xCoord result;
    double v = GetCellVolume();
    double ca = cos(alpha * PI/180);
    double cb = cos(beta  * PI/180);
    double cg = cos(gamma * PI/180);
    double sg = sin(gamma * PI/180);
    
    result.x = a * ox + b * cg * oy + c * cb * oz;
    result.y = b * sg * oy + c * ((ca - cb * cg) / sg) * oz;
    result.z = (v / (a * b * sg)) * oz;
    return (result);
}


xCoord XCrystal::Orthogonalize(xCoord oa)
{
    xCoord result;
    double v = GetCellVolume();
    double ca = cos(alpha * PI/180);
    double cb = cos(beta  * PI/180);
    double cg = cos(gamma * PI/180);
    double sg = sin(gamma * PI/180);
    
    result.x = a * oa.x + b * cg * oa.y + c * cb * oa.z;
    result.y = b * sg * oa.y + c * ((ca - cb * cg) / sg) * oa.z;
    result.z = (v / (a * b * sg)) * oa.z;
    return (result);
}


xCoord XCrystal::Deorthogonalize(double ox, double oy, double oz)
{
    // Transform cartesian into fractional coordinates using the PDB
    // deorthogonalization convention:
    
    xCoord result;
    double v = GetCellVolume();
    double ca = cos(alpha * PI/180);
    double cb = cos(beta  * PI/180);
    double cg = cos(gamma * PI/180);
    double sg = sin(gamma * PI/180);

    result.x = (1 / a) * ox - cg/(a * sg) * oy 
             + (b * c * cg * (ca - cb * cg) / sg - b * c * cb * sg) / v * oz;
    result.y = (1 / (b * sg)) * oy - (a * c * (ca - cb * cg) / (v * sg)) * oz;
    result.z = (a * b * sg) / v * oz;
    return (result);
}


xCoord XCrystal::Deorthogonalize(xCoord oa)
{
    xCoord result;
    double v = GetCellVolume();
    double ca = cos(alpha * PI/180);
    double cb = cos(beta  * PI/180);
    double cg = cos(gamma * PI/180);
    double sg = sin(gamma * PI/180);

    result.x = (1 / a) * oa.x - cg/(a * sg) * oa.y 
             + (b * c * cg * (ca - cb * cg) / sg - b * c * cb * sg) / v * oa.z;
    result.y = (1 / (b * sg)) * oa.y - (a * c * (ca - cb * cg) / (v * sg)) * oa.z;
    result.z = (a * b * sg) / v * oa.z;
    return (result);
}



// >===============--------------
// -----------------------------------------------------------------------------



// -----------------------------------------------------------------------------
// -- CLASS -- XFormFactor -----------------------------------------------------
// -----------------------------------------------------------------------------

// -- Private:

void XFormFactor::GetCCP4FF(string olabel)
{
    // -- Find the AFF corresponding to "olabel" and initialize the object
    //    AFFs are read from $CLIBD/atomsf.lib. This has been replaced by
    //    xglobals.h and could be scrapped ...
    
    string fullpath = getenv("CLIBD");
    fullpath += "/atomsf.lib" ;
    
    char *file = new char[fullpath.length()];
    fullpath.copy(file, fullpath.length());
    file[fullpath.length()] = '\0';          // string.copy() does not append \0 !
    
    ifstream lib(file);
    if (!lib)
    {
        cerr << "\nFILE ERROR: CCP4 form factor library ($CLIBD/atomsf.lib)"
	     << " not found." << endl;
	cerr << " Searched for: [" << file << "]" << endl;
	exit(1);
    }
    
    // -- Copy the input string into a 4-character working string, filling
    //    the string with spaces or clipping additional characters.
    
    int l = olabel.length();
    if ((l < 1) || (l > 4)) 
    {
        cerr << "\nERROR: FormFactor: Illegal Atom identifier." << endl;
	exit(1);
    }

    if ( l < 4 )
        for (int i = 0; i < (4 - l); i++)
	    olabel.append(" ");

    char *dummy = new char[100];
    char *fillmnt = new char[5];
    char *work    = new char[5];
    
    olabel.copy(work, 4);                  // fstream does not take string, mpf
    work[4] = '\0';                        // Add the damn Null by hand!

    // -- Skip header (31 lines until Hydrogen is encountered)
    while (lib.peek() != 'H')
        lib.getline(dummy, 90);
    
    while (lib.peek() != EOF)
    {
	lib.read(fillmnt, 4);
        fillmnt[4] = '\0';                        // Add the damn Null by hand!
	if (strcmp(fillmnt, work) == 0)
	{
	    lib >> iwt >> ielec >> c;
	    lib >> a[0] >> a[1] >> a[2] >> a[3];
	    lib >> b[0] >> b[1] >> b[2] >> b[3];
            lib >> cufp >> cufpp >> mofp >> mofpp;
	    break;
	}
	else
	    for (int i = 0; i < 5; i++) lib.getline(dummy, 90);
    }
    
    // -- This has worked, so put the atom label in the id field.
    id = olabel;
    
    lib.close();
    delete file; delete dummy; delete fillmnt; delete work;
}

// -- Public:

XFormFactor::XFormFactor()
{
    // -- Standard constructor, not to be used so far
    iwt = 0;                         // Flag for empty object
}


XFormFactor::XFormFactor(string olabel)
{
   for (int i = 0; i < TOTAL_AFFS; i++)
   {
      if ( affdata[i].id == olabel )
      {
         id    = affdata[i].id;
	 iwt   = affdata[i].i[0];
	 ielec + affdata[i].i[1];
	 for ( int j = 0; j < 4; j++ )
	 {
	    a[j] = affdata[i].a[j];
	    b[j] = affdata[i].a[j+4];
	 }
         c     = affdata[i].a[8];
	 cufp  = affdata[i].a[9];
	 cufpp = affdata[i].a[10];
	 mofp  = affdata[i].a[11];
	 mofpp = affdata[i].a[12];
         break;
      }
   }
}


double XFormFactor::f(double os, double ob)
{
    // -- Calculate f(sin(theta)/lambda)
    
    double res = 0;
    if (iwt)
    {
        for(int i = 0; i < 4; i++)
	    res +=  a[i] * exp(-b[i] * os * os);
	res += c;
	res *= exp(-ob * os * os);
	return (res);
    }
}


double XFormFactor::f_d(double od, double ob)
{
    // -- Calculate (sin(theta)/lambda) from d and return f(d)
    double res = 0;
    if (od == 0) od = 0.0001;
    double os = 1 / (2 * od);
    if (iwt)
    {
        if (od == 0) od = 0.0001;
        double os = 1 / (2 * od);
        for(int i = 0; i < 4; i++)
	    res +=  a[i] * exp(-b[i] * os * os);
	res += c;
	res *= exp(-ob * os * os);
	return (res);
    }
    return (res);
}

// >===============--------------
// -----------------------------------------------------------------------------




// -----------------------------------------------------------------------------
// -- CLASS -- XAtom -----------------------------------------------------------
// -----------------------------------------------------------------------------

// -- Private:

void XAtom::GetType()
{
}


// -- Public:

XAtom::XAtom()
{
}

XAtom::XAtom(int n, string nm, float xx, float yy, float zz, float qq, float BB)
{
    flag   = ON;                       
    num    = n;                        
    name   = nm;                    
    x      = xx;
    y      = yy;
    z      = zz;                    
    q      = qq;                          
    b      = BB;                          
}


void XAtom::SetUij(Matrix<float> &u)
{
}


void XAtom::SetType(string t)
{
   type = t;
}

  
xCoord XAtom::coord()
{
   xCoord res;
   res.x = x; res.y = y; res.z = z;
   return res;
}

  
inline bool   XAtom::Peek()
{
   return (flag);
}
	
inline void   XAtom::On()
{
   flag = ON;
}
	
inline void   XAtom::Off()
{
   flag = OFF;
}

	
XAtom::~XAtom()
{
}


// >===============--------------
// -----------------------------------------------------------------------------
