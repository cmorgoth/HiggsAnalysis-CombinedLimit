#ifndef HZZ2L2QROOPDFS
#define HZZ2L2QROOPDFS

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"



class RooHMDiphoton : public RooAbsPdf {
 public:
  RooHMDiphoton();
  RooHMDiphoton(const char *name, const char *title,
		RooAbsReal& _x,
		RooAbsReal& _a,
		RooAbsReal& _b
		);
  RooHMDiphoton(const RooHMDiphoton& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooHMDiphoton(*this,newname); }
  inline virtual ~RooHMDiphoton() { };
  
 protected:
  
  RooRealProxy x ;
  RooRealProxy a;
  RooRealProxy b;
  
  Double_t evaluate() const ;
  
 private:
  
  ClassDef(RooHMDiphoton,1)

};

//RooDoubleCBInterpolate
class RooDoubleCBInterpolate : public RooAbsPdf {
 public:
  RooDoubleCBInterpolate();
  RooDoubleCBInterpolate(const char *name, const char *title,
			 RooAbsReal& _x,
			  RooAbsReal& _mass
			 );
  RooDoubleCBInterpolate(const RooDoubleCBInterpolate& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooDoubleCBInterpolate(*this,newname); }
  inline virtual ~RooDoubleCBInterpolate() { }
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;
  Double_t getMean( Double_t m ) const;
  Double_t getSigma( Double_t m ) const;
  Double_t getN1( Double_t m) const;
  Double_t getAlpha1( Double_t m ) const;
  Double_t getN2( Double_t m ) const;
  Double_t getAlpha2( Double_t m ) const;
  
  
 protected:
  
  RooRealProxy x ;
  RooRealProxy mass;
  
  Double_t evaluate() const ;
  
 private:
  
  ClassDef(RooDoubleCBInterpolate,1)
};

class RooCB : public RooAbsPdf {
 public:
  RooCB();
  RooCB(const char *name, const char *title,
        RooAbsReal& _x,
        RooAbsReal& _mean,
        RooAbsReal& _width,
        RooAbsReal& _alpha,
        RooAbsReal& _n,
        RooAbsReal& _theta
	);
  RooCB(const RooCB& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooCB(*this,newname); }
  inline virtual ~RooCB() { }

 protected:

  RooRealProxy x ;
  RooRealProxy mean;
  RooRealProxy width;
  RooRealProxy alpha;
  RooRealProxy n;
  RooRealProxy theta;

  Double_t evaluate() const ;

 private:

  ClassDef(RooCB,1)
    };

 
class RooDoubleCB : public RooAbsPdf {
public:
  RooDoubleCB();
  RooDoubleCB(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _mean,
	      RooAbsReal& _width,
	      RooAbsReal& _alpha1,
	      RooAbsReal& _n1,
	      RooAbsReal& _alpha2,
	      RooAbsReal& _n2
	   );
  RooDoubleCB(const RooDoubleCB& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooDoubleCB(*this,newname); }
  inline virtual ~RooDoubleCB() { }
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

protected:

  RooRealProxy x ;
  RooRealProxy mean;
  RooRealProxy width;
  RooRealProxy alpha1;
  RooRealProxy n1;
  RooRealProxy alpha2;
  RooRealProxy n2;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooDoubleCB,1)
};
 
class RooFermi : public RooAbsPdf {
public:
  RooFermi();
  RooFermi(const char *name, const char *title,
 	    RooAbsReal& _x,
            RooAbsReal& _cutOff,
	   RooAbsReal& _beta
	   );
  RooFermi(const RooFermi& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooFermi(*this,newname); }
  inline virtual ~RooFermi() { }

protected:

  RooRealProxy x ;
  RooRealProxy cutOff ;
  RooRealProxy beta ;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooFermi,1) 
};
  
class RooRelBW : public RooAbsPdf {
public:
  RooRelBW();
  RooRelBW(const char *name, const char *title,
	   RooAbsReal& _x,
	   RooAbsReal& _mean,
	   RooAbsReal& _width,
	   RooAbsReal& _n
	   );
  RooRelBW(const RooRelBW& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooRelBW(*this,newname); }
  inline virtual ~RooRelBW() { }

protected:

  RooRealProxy x ;
  RooRealProxy mean ;
  RooRealProxy width ;
  RooRealProxy n ;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooRelBW,1)
};
 

class Triangle : public RooAbsPdf {
public:
  Triangle();
  Triangle(const char *name, const char *title,                
	   RooAbsReal& _m,
	   RooAbsReal& _start,
	   RooAbsReal& _turn,
	   RooAbsReal& _stop
	   );	
  
  Triangle(const Triangle& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { 
    return new Triangle(*this,newname); }

  inline virtual ~Triangle() { }
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

protected:

  RooRealProxy m;
  RooRealProxy start;
  RooRealProxy turn;
  RooRealProxy stop;
  
  Double_t evaluate() const;

private:
  
  ClassDef(Triangle,1)
};


class RooLevelledExp : public RooAbsPdf {
 public:
  RooLevelledExp();
  RooLevelledExp(const char *name, const char *title,
		 RooAbsReal& _x,
		 RooAbsReal& _sigma,
		 RooAbsReal& _alpha,
		 RooAbsReal& _m,
		 //RooAbsReal& _k,
		 RooAbsReal& _theta
		);

  RooLevelledExp(const RooLevelledExp& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooLevelledExp(*this,newname); }
  inline virtual ~RooLevelledExp() { }

 protected:

  RooRealProxy x ;
  RooRealProxy sigma;
  RooRealProxy alpha;
  RooRealProxy m;
  // RooRealProxy k;
  RooRealProxy theta;
  

  Double_t evaluate() const ;

 private:

  ClassDef(RooLevelledExp,1)
    };

//-------------------------------------
//RooIntepolateDSCB_W0p014_Spin0_EBEB
//--------------------------------------
class RooIntepolateDSCB_W0p014_Spin0_EBEB : public RooAbsPdf {
public:
  RooIntepolateDSCB_W0p014_Spin0_EBEB();
  RooIntepolateDSCB_W0p014_Spin0_EBEB(const char *name, const char *title,
			 RooAbsReal& _x,
			 RooAbsReal& _mass
			 );
  RooIntepolateDSCB_W0p014_Spin0_EBEB(const RooIntepolateDSCB_W0p014_Spin0_EBEB& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooIntepolateDSCB_W0p014_Spin0_EBEB(*this,newname); }
  inline virtual ~RooIntepolateDSCB_W0p014_Spin0_EBEB() { }
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;
  Double_t getMean( Double_t m ) const;
  Double_t getSigma( Double_t m ) const;
  Double_t getN1( Double_t m) const;
  Double_t getAlpha1( Double_t m ) const;
  Double_t getN2( Double_t m ) const;
  Double_t getAlpha2( Double_t m ) const;
  
  
protected:

  RooRealProxy x ;
  RooRealProxy mass;
  
  Double_t evaluate() const ;
  
private:

  ClassDef(RooIntepolateDSCB_W0p014_Spin0_EBEB,1)
};

//-----------------------------------
//RooIntepolateDSCB_W0p014_Spin0_EBEE
//----------------------------------
class RooIntepolateDSCB_W0p014_Spin0_EBEE : public RooAbsPdf {
public:
  RooIntepolateDSCB_W0p014_Spin0_EBEE();
  RooIntepolateDSCB_W0p014_Spin0_EBEE(const char *name, const char *title,
			 RooAbsReal& _x,
			 RooAbsReal& _mass
			 );
  RooIntepolateDSCB_W0p014_Spin0_EBEE(const RooIntepolateDSCB_W0p014_Spin0_EBEE& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooIntepolateDSCB_W0p014_Spin0_EBEE(*this,newname); }
  inline virtual ~RooIntepolateDSCB_W0p014_Spin0_EBEE() { }
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;
  Double_t getMean( Double_t m ) const;
  Double_t getSigma( Double_t m ) const;
  Double_t getN1( Double_t m) const;
  Double_t getAlpha1( Double_t m ) const;
  Double_t getN2( Double_t m ) const;
  Double_t getAlpha2( Double_t m ) const;
  
  
protected:

  RooRealProxy x ;
  RooRealProxy mass;
  
  Double_t evaluate() const ;
  
private:

  ClassDef(RooIntepolateDSCB_W0p014_Spin0_EBEE,1)
};
//-------------------------------------

#endif
