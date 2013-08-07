//---------------------------------------------------------------------------
#include "RooFit.h"

#include "Riostream.h"
#include <TMath.h>
#include <math.h>

#include "RooRazor2DTail_SYS.h"
#include "RooRealVar.h"

ClassImp(RooRazor1DTail_SYS)
//---------------------------------------------------------------------------
RooRazor1DTail_SYS::RooRazor1DTail_SYS(const char *name, const char *title,
				       RooAbsReal &_x, 
				       RooAbsReal &_x0, 
				       RooAbsReal &_b, RooAbsReal &_n) : RooAbsPdf(name, title), 
  X("X", "X Observable", this, _x),
  X0("X0", "X Offset", this, _x0),
  B("B", "Shape parameter", this, _b),
  N("N", "Shape parameter", this, _n)
{
}
//---------------------------------------------------------------------------
RooRazor1DTail_SYS::RooRazor1DTail_SYS(const RooRazor1DTail_SYS& other, const char* name) :
   RooAbsPdf(other, name), 
   X("X", this, other.X), 
   X0("X0", this, other.X0),
   B("B", this, other.B),
   N("N", this, other.N)
{
}
//---------------------------------------------------------------------------
Double_t RooRazor1DTail_SYS::evaluate() const
{
  double myexp = B*N*pow(fabs(X-X0),1./N);
  double mycoeff = B*pow(fabs(X-X0),1./N) - 1.;

  if(myexp < -700.) {
    //std::cout << "MYEXP = "<< myexp << ", < -700 -- BAD" << std::endl;
    return  1.7e-308;}
  if(mycoeff <= 0.) {
    //std::cout << "MYCOEFF = " << mycoeff <<", IS NEGATIVE -- BAD" << std::endl;
    return  1.7e-308;}
  if(X0 >= X.min() || Y0 >= Y.min() || B <= 0. || N <= 0.) {
    //std::cout << "PARAMETERS OUT OF PHYSICAL, INTEGRABLE RANGES -- BAD" << std::endl;
    return  1.7e-308;}

  return mycoeff*exp(-myexp);
}

// //---------------------------------------------------------------------------
Int_t RooRazor1DTail_SYS::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const{
  // integral over X
  if (matchArgs(allVars, analVars, X)) return 1;
  // integrating nothing
  return 0;
}

// //---------------------------------------------------------------------------
Double_t RooRazor1DTail_SYS::analyticalIntegral(Int_t code, const char* rangeName) const{

   const Double_t xmin = X.min(rangeName); const Double_t xmax = X.max(rangeName);
   const Double_t ymin = Y.min(rangeName); const Double_t ymax = Y.max(rangeName);

   if(B <= 0. || N <= 0. || X0 >= X.min()) return 1.;

   double integral = 0.;
   if(code ==0) {
     integral = ( (xmin-X0)*exp(B*N*pow(xmax-X0,1/N)) - (xmax-X0)*exp(B*N*pow(xmin-X0,1/N)) )*exp(-B*N*(pow(xmin-X0,1/N)+pow(xmax-X0,1/N)));
     
   }
   else {
     cout << "WARNING IN RooRazor1DTail_SYS: integration code is not correct" << endl;
     cout << "                           what are you integrating on?" << endl;
     return 1.;
   }

   return integral;
}
// //---------------------------------------------------------------------------

