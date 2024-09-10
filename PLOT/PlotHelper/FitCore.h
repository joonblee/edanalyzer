#ifndef FitCore_H
#define FitCore_H

Double_t Pol2Func(Double_t *x, Double_t *par) {
 return par[0]*(x[0]-par[1])*(x[0]-par[1])+par[2];
 //return par[0]+par[1]*x[0]+par[2]*x[0]*x[0];
}

Double_t Pol3Func(Double_t *x, Double_t *par) {
 return par[0]+par[1]*x[0]+par[2]*x[0]*x[0]+par[3]*x[0]*x[0]*x[0];
}

Double_t Pol4Func(Double_t *x, Double_t *par) {
 return par[0]+par[1]*x[0]+par[2]*x[0]*x[0]+par[3]*x[0]*x[0]*x[0]+par[4]*x[0]*x[0]*x[0]*x[0];
}

Double_t Pol6Func(Double_t *x, Double_t *par) {
 return par[0]+par[1]*x[0]+par[2]*x[0]*x[0]+par[3]*pow(x[0],3)+par[4]*pow(x[0],4)+par[5]*pow(x[0],5)+par[6]*pow(x[0],6);
}


Double_t PowFunc(Double_t *x, Double_t *par) {
 return par[0]*pow(x[0]-par[1], par[3])+par[2];
 //return par[0]+par[1]*x[0]+par[2]*x[0]*x[0];
}

Double_t xn_exp(Double_t *x, Double_t *par) {
  return par[0]*pow(x[0],par[1])*TMath::Exp(par[2]*x[0]);
}

Double_t exp(Double_t *x, Double_t *par) {
  return par[0]*TMath::Exp(par[1]*x[0]);
}

Double_t powx(Double_t *x, Double_t *par) {
  return par[0]*pow(par[1],x[0]);
}

Double_t xn(Double_t *x, Double_t *par) {
  return par[0]*pow(x[0]-par[2],par[1]);
}

Double_t LogNormal(Double_t *x, Double_t *par) {
  return par[0]*pow(x[0],-par[1])*TMath::Exp(-(TMath::Log(x[0])-par[2])*(TMath::Log(x[0])-par[2])/2/par[3]/par[3]);
}

Double_t xnExp2(Double_t *x, Double_t *par) {
  return par[0]*pow(x[0],par[1])*TMath::Exp(par[2]*(x[0]-par[3])*(x[0]-par[3]));
}

Double_t Gaussian(Double_t *x, Double_t *par) {
  return par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/2/par[2]/par[2]);
}
Double_t Cauchy(Double_t *x, Double_t *par) { // Breit-Wigner
  return par[0]/(1+(x[0]-par[1])*(x[0]-par[1])/par[2]/par[2]);
}
Double_t RelBreitWigner(Double_t *x, Double_t *par) {
  return par[0]/((x[0]*x[0]-par[1]*par[1])*(x[0]*x[0]-par[1]*par[1])+par[1]*par[1]*par[2]*par[2]);
}
Double_t CrystalBall(Double_t *x, Double_t *par) {
  if((x[0]-par[1])/par[2]<par[4])
    return par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/2/par[2]/par[2]);
  else
    return par[0]*pow(par[3]/par[4],par[3])*TMath::Exp(-par[4]*par[4]/2)*pow(par[3]/par[4]-par[4]+(x[0]-par[1])/par[2],-par[3]);
}
Double_t LogNormalWithExp(Double_t *x, Double_t *par) {
  if(x[0]<par[4])
    return par[0]/x[0]*TMath::Exp(-(TMath::Log(x[0])-par[1])*(TMath::Log(x[0])-par[1])/2/par[2]/par[2]);
  else
    return par[0]/par[4]*TMath::Exp(-(TMath::Log(par[4])-par[1])*(TMath::Log(par[4])-par[1])/2/par[2]/par[2])*pow(1-(1-x[0]/par[4])*(1+(TMath::Log(par[4])-par[1])/par[2]/par[2])/par[3],-par[3]);
}
Double_t GammaWithExp(Double_t *x, Double_t *par) { // k:par[1], theta:par[2]
  if(x[0]<par[4])
    return par[0]*pow(x[0],par[1]-1)*TMath::Exp(-x[0]/par[2]);
  else
    return par[0]*pow(par[4],par[1]-1)*TMath::Exp(-par[4]/par[2])*pow(1+(1/par[2]-(par[1]-1)/par[4])*(x[0]-par[4])/par[3],-par[3]);
}

Double_t test(Double_t *x, Double_t *par) {
  /*
  Double_t xx=x[0]/13000;
  return par[0]*pow((1-xx),par[1])/pow(xx,par[2]+par[3]*TMath::Log(xx));
  */

  //return par[0]*pow(x[0],par[1])*TMath::Exp(par[2]*x[0])+par[3];

  //return par[0]*pow(x[0],par[1])*(TMath::Exp(par[2]*x[0])-par[3])*(TMath::Exp(par[2]*x[0])-par[3]);

  //Double_t expx=TMath::Exp(x[0]);
  if(x[0]<par[4])
    return par[0]*pow(x[0],par[1])*TMath::Exp(par[2]*x[0]);
  else
    return par[0]*pow(par[4],par[1]-par[3])*TMath::Exp(par[2]*par[4])*pow(x[0],par[3]);


  //return par[0]*pow(x[0],par[1])*TMath::Exp(par[2]*x[0]-par[3]);
  //return par[0]*pow(x[0],par[1])*pow(par[3],par[2]*x[0]);
}

#endif
