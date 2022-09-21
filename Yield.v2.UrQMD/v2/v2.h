#include <Math/IFunction.h>
#include <Math/SpecFuncMathMore.h>
#include <TMath.h>


using namespace ROOT::Math;

class v2 :  public ROOT::Math::IBaseFunctionMultiDim {
 public: 
  
  double DoEval(const double * x) const {
    
    double xx = x[0]; // variable 1
    double yy = x[1]; // variable 2

    double VolT = 450;
    double AnchoDisco = 0.25*TMath::Power(0.197,-1);
    double B = 0.0249;
    float mq = 0.0001;
    int g = 2;
    double alphas = TMath::Power(TMath::Pi(),-1);
    double alphaem = TMath::Power(137,-1);
    double lambdas=2;
    double eq =0.6; //{2/3, 1/3, 1/3}
    double vol= AnchoDisco*VolT*alphaem*alphas*alphas*eq*eq*TMath::Power(2,-1)*TMath::Power(2*TMath::Pi(),-6);
    float beta=0.25;
    double eta=3;
    double betaf=TMath::Power(1-beta,2)*TMath::Power(1-TMath::Power(beta,2),-1);
    double omegaq;

    double n1func = eta*TMath::Power(TMath::Exp(TMath::Sqrt(TMath::Power(xx,2)*betaf)*
						TMath::Power(lambdas,-1))-1,-1);
    double n2func = eta*TMath::Power(TMath::Exp(TMath::Sqrt(TMath::Power(0.1,2)+
							    TMath::Power(yy-xx,2))*betaf*
					       TMath::Power(lambdas,-1))-1,1);
    double argb = betaf*((TMath::Power(xx,2))-(xx*yy)+(TMath::Power(yy,2)))*
      TMath::Power(2*B*eq,-1);
   double argB1 = (1+betaf*argb)*TMath::Power(betaf*argb,-1);
   double pol =  betaf*(2*TMath::Power(xx,2)-(xx*yy)+TMath::Power(yy,2));
   
   double val = n2func;
     /*double val = TMath::Exp(-argb)*n1func*n2func*pol*
     (cyl_bessel_i(0.0,argb)-(argB1*cyl_bessel_i(1.0,argb)));
     */
return val; 
  }
  unsigned int NDim() const
  {
    return 2;
  }
  
  ROOT::Math::IBaseFunctionMultiDim* Clone() const
    {
      return new v2();
    }
  
 protected:
  
 private:
  
};
