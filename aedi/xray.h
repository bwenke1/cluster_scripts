/*


    XRAY.H                                          Oli E. 16.10.2002
    
    Header file for Xray classes
    
    For the time being this should be the onlt file to include in all programs
    using the Xray classes. All functions needed go to xray.cpp and the global
    constants (form factors and symmetry operators) are in xglobals.h
    
    16.10.2002 - Added class XCrystal
    17.10.2002 - Added funct FormFactor
    18.10.2002 - Removd funct FormFactor, added class xFormFactor
    20.10.2002 - Added classes Vector, Matrix and Matrix3D
               - Added XCrystal::Orthogonalize(), XCrystal::Deorthogonalize()
    21.10.2002 - Added Biso to XFormFactor::f() and XFormFactor::f_d()	       
    22.10.2002 - Added GaussJordan()
               - Added PearsonsR()	
    29.04.2006 - Added XCrystal::XCrystal(char *filename) 
    06.06.2006 - Included SymOp structure
    12.06.2006 - Removed class Vector (STL vectors to be used)  
    16.09.2006 - Changed XFormFactor constructor to use xglobals.h   

*/

#include <algorithm>
#include <string>
#include <complex>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>

#include <cmath>     
#include <cstdio>
#include <cstring>


using namespace std;


// -- Global constants ---------------------------------------------------------


const double PI   = 3.1415926535;
const double CUKa = 1.5418;

const bool   ON   = 1;
const bool   OFF  = 0;


// -- Useful bits --------------------------------------------------------------

template<class T>
inline void SWAP(T &a, T &b)                        // From Numerical Recipes
    { T dum = a; a = b; b = dum; }


// -- Useful structures --------------------------------------------------------

// Basics

struct xHKL     {  int    H, K, L; };
struct xCoord   {  double x, y, z; };

// Waves

struct xFphi    {  double F, phi; };                // Amplitude and Phase
struct xFaFb    {  double Fa, Fb; };                // Vector components

// Symmetry Operators

struct SymOp    {  int r[3][3]; float t[3]; };      // Evtly worth a class
    
// >===============--------------




// -- CLASS -- Matrix ----------------------------------------------------------
//
//  20.10.2002 This is an implementation (copy) of the basic NRMat class from                   
//             Numerical Recipes in C++, offering a templated vector class with
//             operator overloading and basic functionality. See the book for
//             a desctiption of the member functions
//             This implementation uses indices starting from 0!
//
// -----------------------------------------------------------------------------

template<class T>
class Matrix
{
  private:
    int nn;                                   // Number of rows
    int mm;                                   // Number of columns
    T **v;                                    // ponter to the data
    
  public:
    Matrix();
    Matrix(int n, int m);                     // Zero-based array
    Matrix(const T &a, int n, int m);         // Initialize to constant
    Matrix(const T *a, int n, int m);         // Initialize to array
    Matrix(const Matrix &rhs);                // Copy constructor
    Matrix & operator=(const Matrix &rhs);    // Assignment
    Matrix & operator=(const T &a);           // Assignment to every element
    inline T* operator[](const int i);        // Subsctipting: Pointer to row i
    inline const T* operator[](const int i) const;
    inline int nrows() const;                 // Returns number of rows
    inline int ncols() const;                 // Returns number of columns
    ~Matrix();                                // Destructor
};

// -- Function templates

template<class T>
Matrix<T>::Matrix() : nn(0), mm(0), v(0) {}


template<class T>
Matrix<T>::Matrix(int n, int m) : nn(n), mm(m), v(new T*[n])
{
    v[0] = new T[n*m];
    for(int i = 1; i < n; i++)
       v[i] = v[i-1] + m;
}


template<class T>
Matrix<T>::Matrix(const T &a, int n, int m) : nn(n), mm(m), v(new T*[n])
{
    int i, j;
    v[0] = new T[n*m];
    for(int i = 1; i < n; i++)
       v[i] = v[i-1] + m;
    for(int i = 0; i < n; i++)
       for(int j = 0; j < m; j++)
         v[i][j] = a;
}


template<class T>
Matrix<T>::Matrix(const T *a, int n, int m) : nn(n), mm(m), v(new T*[n])
{
    int i, j;
    v[0] = new T[n*m];
    for(int i = 1; i < n; i++)
       v[i] = v[i-1] + m;
    for(int i = 0; i < n; i++)
       for(int j = 0; j < m; j++)
          v[i][j] = *a++;
}


template<class T>
Matrix<T>::Matrix(const Matrix &rhs) : nn(rhs.nn), mm(rhs.mm), v(new T*[nn])
{
    int i, j;
    v[0] = new T[nn*mm];
    for(int i = 1; i < nn; i++)
       v[i] = v[i-1] + mm;
    for(int i = 0; i < nn; i++)
       for(int j = 0; j < mm; j++)
          v[i][j] = rhs[i][j];
}


template<class T>
Matrix<T> & Matrix<T>::operator=(const Matrix<T> &rhs)
{
    if (this != &rhs)
    {
       int i, j;
       if (nn != rhs.nn || mm != rhs.mm)
       {
	  if (v != 0)
	  {
	     delete [] (v[0]);
	     delete [] (v);
	  }
	  nn   = rhs.nn;
	  mm   = rhs.mm;
	  v    = new T*[nn];
	  v[0] = new T[nn*mm];   
       }
       for(int i = 1; i < nn; i++)
          v[i] = v[i-1] + mm;
       for(int i = 0; i < nn; i++)
          for(int j = 0; j < mm; j++)
             v[i][j] = rhs[i][j];
    }
    return *this;
}


template<class T>
Matrix<T> & Matrix<T>::operator=(const T &a)
{
    for(int i = 0; i < nn; i++)
       for(int j = 0; j < mm; j++)
	  v[i][j] = a;
    return *this;	  
}


template<class T>
inline T* Matrix<T>::operator[](const int i)
{
    return v[i];
}


template<class T>
inline const T* Matrix<T>::operator[](const int i) const
{
    return v[i];
}


template<class T>
inline int Matrix<T>::nrows() const
{
    return nn;
}


template<class T>
inline int Matrix<T>::ncols() const
{
    return mm;
}


template<class T>
Matrix<T>::~Matrix()
{
   if (v != 0)
   {
     delete [] (v[0]);
     delete [] (v);
   }
}


// >===============--------------            




// -- CLASS -- Matrix3D --------------------------------------------------------
//
//  20.10.2002 This is an implementation (copy) of the basic NRMat3D class from                   
//             Numerical Recipes in C++, offering a templated vector class with
//             operator overloading and basic functionality. See the book for
//             a desctiption of the member functions
//             This implementation uses indices starting from 0!
//
// -----------------------------------------------------------------------------

template<class T>
class Matrix3D
{
  private:
    int nn;                                   // Dimension 1
    int mm;                                   // Dimension 2
    int kk;                                   // Dimension 3
    T ***v;                                   // pointer to the data
    
  public:
    Matrix3D();
    Matrix3D(int n, int m, int k);            // Zero-based array
    inline T** operator[](const int i);       // Subsctipting: Pointer to row i
    inline const T* const * operator[](const int i) const;
    inline int dim1() const;                  // Returns dimension 1
    inline int dim2() const;                  // Returns dimension 2
    inline int dim3() const;                  // Returns dimension 3
    ~Matrix3D();                              // Destructor
};

// -- Function templates

template<class T>
Matrix3D<T>::Matrix3D() : nn(0), mm(0), kk(0), v(0) {}


template<class T>
Matrix3D<T>::Matrix3D(int n, int m, int k) : nn(n), mm(m), kk(k), v(new T**[n])
{
    int i, j;
    v[0]    = new T*[n*m];
    v[0][0] = new T[n*m*k];
    for(int j = 1; j < m; j++)
       v[0][j] = v[0][j-1] + k;
    for(int i = 1; i < n; i++)
    {
       v[i]    = v[i-1] + m;
       v[i][0] = v[i-1][0] + m*k;
       for(int j = 1; j < m; j++)
          v[i][j] = v[i][j-1] + k;
    }	  
}


template<class T>
inline T** Matrix3D<T>::operator[](const int i)
{
    return v[i];
}


template<class T>
inline const T* const * Matrix3D<T>::operator[](const int i) const
{
    return v[i];
}


template<class T>
inline int Matrix3D<T>::dim1() const
{
    return nn;
}


template<class T>
inline int Matrix3D<T>::dim2() const
{
    return mm;
}


template<class T>
inline int Matrix3D<T>::dim3() const
{
    return kk;
}


template<class T>
Matrix3D<T>::~Matrix3D()
{
    if (v != 0)
    {
        delete [] (v[0][0]);
        delete [] (v[0]);
        delete [] (v);
    }
}

// >===============--------------            


// -- CLASS -- XCrystal --------------------------------------------------------
//
//  16.10.2002 Designed to hold information about unit cell dimensions and
//             spacegroup. Allows calculations of cell volume, Matthews coef-
//             ficient, conversions HKL -> resolution.
//
//             Implemented Member functions:
//
//           - GetAsymUnits: (incomplete)
//
//           - XCrystal: Constructor.
//           - GetResolution: calculate d from (HKL) according to Bragg's Law
//           - GetCellVolume: calculate the parallelepiped volume
//           - GetMatthews: calculate V(m) from mol.weight and number of mols.
//
//  29.04.2006 Implemented Member function>
//
//           - XCrystal(char *filename): Get cell from PDB or MTZ file
//
// -----------------------------------------------------------------------------

class XCrystal
{
    double a, b, c;                 // unit cell axes in Angstroms
    double alpha,                   // unit cell angles in degrees
           beta,
	   gamma;
    
    int    spg_num;                 // spacegroup number
    int    asym_units;              // number of asymmetric units

  private:
    int    GetAsymUnits();
    
  public:
           XCrystal();
    	   XCrystal(double oa, double ob, double oc, 
                    double oalpha, double obeta, double ogamma,
		    int ospg);
           XCrystal(char *filename);
    
    double GetResolution(int oH, int oK, int oL);
    double GetResolution(xHKL oHKL);
    double GetCellVolume();
    double GetMatthews  (double omolweight, int on);
    xCoord Orthogonalize(double ox, double oy, double oz);
    xCoord Orthogonalize(xCoord oa);
    xCoord Deorthogonalize(double ox, double oy, double oz);
    xCoord Deorthogonalize(xCoord oa);

};

// >===============--------------

	   

// -- CLASS -- XFormFactor -----------------------------------------------------
//
//  18.10.2002 Class to hold a single atom form factor which, on initialization,
//             is read from $CLIBD/atomsf.lib. Allows for the calculation of 
//             f from (sin(theta)/lambda) or from resolution.
//
//             Implemented Member functions:
//
//           - GetFormFactor(): reads aff data from $CLIBD/atomsf.lib
//
//           - XFormFactor(): Constructor, reads atom_id
//           - f(): returns f(sin(theta)/lambda)
//           - f_d(): returns f(d)
//
//  21.10.2002 included Biso into f() and f_d()
//
//  09.06.2006 fixed name transfer to id field. Names in id match atomsf.lib!
//
//             Implemented Member function:
//
//           - ID(): returns the element name (id field)
//
// -----------------------------------------------------------------------------

class XFormFactor
{
    string id;
    int iwt, ielec;
    double a[4], b[4], c;
    double cufp, cufpp, 
           mofp, mofpp;
	   
  private:
    void   GetCCP4FF(string olabel);  // Read aff from $CLIBD/atomsf.lib (outdated)

  public:
           XFormFactor();
           XFormFactor(string olabel);
    string ID()    const { return id; }
    int    IWT()   const { return iwt; }
    int    IELEC() const { return ielec; }
    double f(double os, double ob = 0);
    double f_d(double od, double ob = 0);
}; 


// >===============--------------




// -- CLASS -- XAtom -----------------------------------------------------------
//
//  23.10.2002 Without being overdemocritical, the atom is the smallest building
//             block for a structure. Basically the ATOM/HETATM and ANISOU lines
//             of a PDB file.
//
//  12.06.2006 A simple way of coding a structure would be as a container of
//             XAtom objects, more sophisticated approaches with more hierarchy.
//             Changed some members and return values. 
//
// -----------------------------------------------------------------------------

class XAtom
{
  private:
    bool   flag;                         // On/Off switch
    int    num;                          // Atom number
    string name;                         // Atom name from PDB File
    string type;                         // Atom type for form factor
    float  x, y, z;                      // Coordinates
    float  q;                            // Occupancy
    float  b;                            // Isotropic temperature factor
    Matrix<float> U[3][2];               // Anisotropic B-factor
    
    void GetType();                      // Determine atom type
    
  public:
    XAtom();
    XAtom(int n, string nm,              // Constructor from PDB file
         float xx, float yy, float zz, float qq, float BB);
    void SetUij(Matrix<float> &u);       // Set AnisoUs
    void SetType(string t);              // Define atom type (charge etc.)
    xCoord coord();                      // Get coordinates
    float X() const { return x; }        // Get x
    float Y() const { return y; }        // Get x
    float Z() const { return z; }        // Get x
    float Q() const { return q; }        // Get occupancy
    float B() const { return b; }        // Get temperature factor
    inline bool   Peek();                // Return switch status
    inline void   On();                  // Switch on
    inline void   Off();                 // Switch off
    ~XAtom();                            // The atomic destructor!!
    
};

// >===============--------------




// -- Function protoypes -------------------------------------------------------

void   GaussJordan(Matrix<double> &A, Matrix<double> &b);
double PearsonsR(vector<double> &x, vector<double> &y);

// >===============--------------


