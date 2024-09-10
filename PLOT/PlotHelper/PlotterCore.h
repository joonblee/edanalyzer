#ifndef PlotterCore_H
#define PlotterCore_H

#include "../PlotHelper/tdrstyle.C"
#include "../PlotHelper/CMS_lumi.C"

// Replace text
TString Replace(TString& title, const TString& from, const TString& to) {
  unsigned pos = 0;
  while( (pos = title.Index(from, pos)) != kNPOS ) {
    title.Replace( pos, from.Length(), to);
  }
  return title;
}

// Function to draw a TH2D histogram with title and axis labels
void drawTH2D(TH2D* hist, TString title, const TString& x_label, const TString& y_label, bool logx = false, bool logy = false) {
  TCanvas* canvas = new TCanvas("canvas","",1000,1000);
  canvas->cd();

  if(logx) canvas->SetLogx();
  if(logy) canvas->SetLogy();

  hist->Draw("COLZ");
  hist->SetMarkerSize(0.6);
  hist->Draw("TEXT SAME");
  hist->SetStats(0);
  hist->SetTitle(title);

  if(logx) {
    hist->GetXaxis()->SetMoreLogLabels();
    hist->GetXaxis()->SetNoExponent();
  }
  if(logy) {
    hist->GetYaxis()->SetMoreLogLabels();
    hist->GetYaxis()->SetNoExponent();
  }

  hist->GetXaxis()->SetLabelSize(0.035);
  hist->GetXaxis()->SetTitle(x_label);
  hist->GetXaxis()->SetTitleSize(0.035);
  hist->GetXaxis()->SetTitleOffset(1.0);

  hist->GetYaxis()->SetLabelSize(0.035);
  hist->GetYaxis()->SetTitle(y_label);
  hist->GetYaxis()->SetTitleSize(0.035);
  hist->GetYaxis()->SetTitleOffset(1.3);

  gStyle->SetPaintTextFormat(".1f");
  if( title.Index("epsilon",0)!=kNPOS || title.Index("eff",0)!=kNPOS || title.Index("SF",0)!=kNPOS ) {
    gStyle->SetPaintTextFormat(".3f");
    hist->SetMaximum(1.0); // Set maximum for efficiency plot
  }

  title = Replace( title, " ", "_" );
  title = Replace( title, "#epsilon", "eff" );
  title = Replace( title, "(", "_" );
  title = Replace( title, ")", "" );

  canvas->SaveAs("../outputs/"+title+".png");

  delete canvas;
}


/*
void SetDataStyle(TH1D *h, int bin, int MarkerStyle, double MarkerSize , int MarkerColor) {
  h->Rebin(bin);
  h->SetMarkerStyle(MarkerStyle);
  h->SetMarkerSize(MarkerSize);
  h->SetLineColor(MarkerColor);
  h->SetStats(0);
}
*/

TH1D* SetDataStyle(TH1D *h, int MarkerStyle, double MarkerSize , int MarkerColor) {
  h->SetMarkerStyle(MarkerStyle);
  h->SetMarkerSize(MarkerSize);
  h->SetLineColor(MarkerColor);
  h->SetStats(0);
  return h;
}

TH1D* SetDataStyle(TH1D *h, int bin, int MarkerStyle, double MarkerSize , int MarkerColor) {
  h->Rebin(bin);
  h->SetMarkerStyle(MarkerStyle);
  h->SetMarkerSize(MarkerSize);
  h->SetLineColor(MarkerColor);
  h->SetStats(0);
  return h;
}

TH1D* SetDataStyle(TH1D *h, vector<double> xbins, int MarkerStyle, double MarkerSize , int MarkerColor) {
  const int Nxbin = xbins.size();
  Double_t xbin[Nxbin];
  for(unsigned i=0; i<Nxbin; i++) {
    xbin[i]=xbins[i];
  }
  TH1D* out = (TH1D*)h->Rebin(Nxbin-1,h->GetName(),xbin);
  for(unsigned i=0; i<Nxbin; i++) {
    out->SetBinContent(i,out->GetBinContent(i)/((TAxis*)out->GetXaxis())->GetBinWidth(i));
    out->SetBinError(i,out->GetBinError(i)/((TAxis*)out->GetXaxis())->GetBinWidth(i));
  }
  out->SetMarkerStyle(MarkerStyle);
  out->SetMarkerSize(MarkerSize);
  out->SetLineColor(MarkerColor);
  out->SetStats(0);
  return out;
}
/*
void SetMCStyle(TH1D *h, int bin, int LineColor, int FillColor, double LineWidth, int LineStyle=1) {
  h->Rebin(bin);
  h->SetLineColor(LineColor);
  h->SetFillColor(FillColor);
  h->SetLineWidth(LineWidth);
  h->SetLineStyle(LineStyle);
}
*/
TH1D* SetMCStyle(TH1D *h, int LineColor, int FillColor, double LineWidth, int LineStyle=1) {
  h->SetLineColor(LineColor);
  if(FillColor!=0) h->SetFillColor(FillColor);
  h->SetLineWidth(LineWidth);
  h->SetLineStyle(LineStyle);
  return h;
}

TH1D* SetMCStyle(TH1D *h, int bin, int LineColor, int FillColor, double LineWidth, int LineStyle=1) {
  h->Rebin(bin);
  h->SetLineColor(LineColor);
  if(FillColor!=0) h->SetFillColor(FillColor);
  h->SetLineWidth(LineWidth);
  h->SetLineStyle(LineStyle);
  return h;
}
TH1D* SetMCStyle(TH1D* h, vector<double> xbins, int LineColor, int FillColor, double LineWidth, int LineStyle=1) {
  const int Nxbin = xbins.size();
  Double_t xbin[Nxbin];
  for(unsigned i=0; i<Nxbin; i++) {
    xbin[i]=xbins[i];
  }
  TH1D* out = (TH1D*)h->Rebin(Nxbin-1,h->GetName(),xbin);
  for(unsigned i=0; i<Nxbin; i++) {
    out->SetBinContent(i,out->GetBinContent(i)/((TAxis*)out->GetXaxis())->GetBinWidth(i));
    out->SetBinError(i,out->GetBinError(i)/((TAxis*)out->GetXaxis())->GetBinWidth(i));
  }
  out->SetLineColor(LineColor);
  if(FillColor!=0) out->SetFillColor(FillColor);
  out->SetLineWidth(LineWidth);
  out->SetLineStyle(LineStyle);
  return out;
}

/*
void SetPlotStyle(TH1D *h, int bin=1,
                  int MarkerColor=kBlack, int MarkerStyle=20, double MarkerSize=.4,
                  int LineColor=kBlack,   int LineStyle=1,    double LineWidth=1.,
                  int FillColor=kBlack,   int FillStyle=1001) {
  h->Rebin(bin);
  if(MarkerStyle!=0) {
    h->SetMarkerColor(MarkerColor);
    h->SetMarkerStyle(MarkerStyle);
    h->SetMarkerSize(MarkerSize);
  }
  if(LineStyle!=0) {
    h->SetLineColor(LineColor);
    h->SetLineStyle(LineStyle);
    h->SetLineWidth(LineWidth);
  }
  if(FillStyle!=0) {
    h->SetFillColor(FillColor);
    h->SetFillStyle(FillStyle);
  }
}
*/
TH1D* SetPlotStyle(TH1D *h, int bin=1,
                  /* Marker Style */
                  /* +:2, *:3, o:4, x:5,
                     filled o:20, rect:21, tri:22, inverted tri:23,
                     empty o:24, rect:25, tri:26, inverted tri:32 */
                  /* Line Style */
                  /* Dashed: 2, 7, 9 */
                  /* Fill Style */
                  /* hollow:0, solid:1001, pattern:3000+@ */
                  int MarkerColor=kBlack, int MarkerStyle=20, double MarkerSize=.4,
                  int LineColor=kBlack,   int LineStyle=1,    double LineWidth=1.,
                  int FillColor=kBlack,   int FillStyle=1001) {
  h->Rebin(bin);
  if(MarkerStyle!=0) {
    h->SetMarkerColor(MarkerColor);
    h->SetMarkerStyle(MarkerStyle);
    h->SetMarkerSize(MarkerSize);
  }
  if(LineStyle!=0) {
    h->SetLineColor(LineColor);
    h->SetLineStyle(LineStyle);
    h->SetLineWidth(LineWidth);
  }
  if(FillStyle!=0) {
    h->SetFillColor(FillColor);
    h->SetFillStyle(FillStyle);
  }
  return h;
}
TH1D* SetPlotStyle(TH1D *h, vector<double> xbins,
                  /* Marker Style */
                  /* +:2, *:3, o:4, x:5,
                     filled o:20, rect:21, tri:22, inverted tri:23,
                     empty o:24, rect:25, tri:26, inverted tri:32 */
                  /* Line Style */
                  /* Dashed: 2, 7, 9 */
                  /* Fill Style */
                  /* hollow:0, solid:1001, pattern:3000+@ */
                  int MarkerColor=kBlack, int MarkerStyle=20, double MarkerSize=.4,
                  int LineColor=kBlack,   int LineStyle=1,    double LineWidth=1.,
                  int FillColor=kBlack,   int FillStyle=1001) {
  const int Nxbin = xbins.size();
  Double_t xbin[Nxbin];
  for(unsigned i=0; i<Nxbin; i++) {
    xbin[i]=xbins[i];
  }
  TH1D* out = (TH1D*)h->Rebin(Nxbin-1,h->GetName(),xbin);
  for(unsigned i=0; i<Nxbin; i++) {
    out->SetBinContent(i,out->GetBinContent(i)/((TAxis*)out->GetXaxis())->GetBinWidth(i));
    out->SetBinError(i,out->GetBinError(i)/((TAxis*)out->GetXaxis())->GetBinWidth(i));
  }
  if(MarkerStyle!=0) {
    out->SetMarkerColor(MarkerColor);
    out->SetMarkerStyle(MarkerStyle);
    out->SetMarkerSize(MarkerSize);
  }
  if(LineStyle!=0) {
    out->SetLineColor(LineColor);
    out->SetLineStyle(LineStyle);
    out->SetLineWidth(LineWidth);
  }
  if(FillStyle!=0) {
    out->SetFillColor(FillColor);
    out->SetFillStyle(FillStyle);
  }
  return out;
}


void SetXRange(double& XRangeMin, double& XRangeMax, TH1D *h, TString Variable = "") {
  vector<TString> Variables = {"Eta","Phi"};
  for(unsigned iVar = 0; iVar < Variables.size(); iVar++) {
    if( Variable.Index( Variables[iVar] ) != kNPOS ) return;
  }

  XRangeMin = h->GetXaxis()->GetXmin();
  XRangeMax = XRangeMin + h->FindLastBinAbove(0)*((TAxis*)h->GetXaxis())->GetBinWidth(0);

  if( XRangeMax == h->GetXaxis()->GetXmax() ) return;

  if( XRangeMax < 1 ) return;

  TString s = TString::Itoa( XRangeMax , 10 );

  int XRangeDigit = s.Sizeof() - 1;

  TString s0 = s[0]; int i0 = s0.Atoi();
  TString s1 = s[1]; int i1 = s1.Atoi();

  if( XRangeDigit < 3 || i0 > 3 ) 
    XRangeMax = ( i0 + 1 ) * pow(10, XRangeDigit-1);
  else
    XRangeMax = i0 * pow(10, XRangeDigit-1) + (i1+1) * pow(10, XRangeDigit-2);
}

void SetXRange(double& XRangeMin, double& XRangeMax, TH1D *h1, TH1D *h2, TString Variable = "") {
  vector<TString> Variables = {"Eta","Phi"};
  for(unsigned iVar = 0; iVar < Variables.size(); iVar++) {
    if( Variable.Index( Variables[iVar] ) != kNPOS ) return;
  }

  XRangeMin = 0;
  XRangeMax = max( h1->FindLastBinAbove(0)*((TAxis*)h1->GetXaxis())->GetBinWidth(0),
                   h2->FindLastBinAbove(0)*((TAxis*)h2->GetXaxis())->GetBinWidth(0) );

  if( XRangeMax == h1->GetXaxis()->GetXmax() ) return;

  if( XRangeMax < 1 ) return;

  TString s = TString::Itoa( XRangeMax , 10 );

  int XRangeDigit = s.Sizeof() - 1;

  TString s0 = s[0]; int i0 = s0.Atoi();
  TString s1 = s[1]; int i1 = s1.Atoi();

  if( XRangeDigit < 3 || i0 > 3 ) 
    XRangeMax = ( i0 + 1 ) * pow(10, XRangeDigit-1);
  else
    XRangeMax = i0 * pow(10, XRangeDigit-1) + (i1+1) * pow(10, XRangeDigit-2);
}
void SetXRange(double& XRangeMin, double& XRangeMax, TH1D *h1, TH1D *h2, TH1D *h3, TString Variable = "") {
  vector<TString> Variables = {"Eta","Phi"};
  for(unsigned iVar = 0; iVar < Variables.size(); iVar++) {
    if( Variable.Index( Variables[iVar] ) != kNPOS ) return;
  }

  XRangeMin = 0;
  XRangeMax = max( max(
                h1->FindLastBinAbove(0)*((TAxis*)h1->GetXaxis())->GetBinWidth(0),
                h2->FindLastBinAbove(0)*((TAxis*)h2->GetXaxis())->GetBinWidth(0) ), 
                h3->FindLastBinAbove(0)*((TAxis*)h3->GetXaxis())->GetBinWidth(0) );

  if( XRangeMax == h1->GetXaxis()->GetXmax() ) return;

  if( XRangeMax < 1 ) return;

  TString s = TString::Itoa( XRangeMax , 10 );

  int XRangeDigit = s.Sizeof() - 1;

  TString s0 = s[0]; int i0 = s0.Atoi();
  TString s1 = s[1]; int i1 = s1.Atoi();

  if( XRangeDigit < 3 || i0 > 3 ) 
    XRangeMax = ( i0 + 1 ) * pow(10, XRangeDigit-1);
  else
    XRangeMax = i0 * pow(10, XRangeDigit-1) + (i1+1) * pow(10, XRangeDigit-2);
}
void SetXRange(double& XRangeMin, double& XRangeMax, TH1D *h1, TH1D *h2, TH1D *h3, TH1D *h4, TString Variable = "") {
  vector<TString> Variables = {"Eta","Phi"};
  for(unsigned iVar = 0; iVar < Variables.size(); iVar++) {
    if( Variable.Index( Variables[iVar] ) != kNPOS ) return;
  }

  XRangeMin = 0;
  XRangeMax = max( max( max(
                h1->FindLastBinAbove(0)*((TAxis*)h1->GetXaxis())->GetBinWidth(0),
                h2->FindLastBinAbove(0)*((TAxis*)h2->GetXaxis())->GetBinWidth(0) ), 
                h3->FindLastBinAbove(0)*((TAxis*)h3->GetXaxis())->GetBinWidth(0) ),
                h4->FindLastBinAbove(0)*((TAxis*)h4->GetXaxis())->GetBinWidth(0) );

  if( XRangeMax == h1->GetXaxis()->GetXmax() ) return;

  if( XRangeMax < 1 ) return;

  TString s = TString::Itoa( XRangeMax , 10 );

  int XRangeDigit = s.Sizeof() - 1;

  TString s0 = s[0]; int i0 = s0.Atoi();
  TString s1 = s[1]; int i1 = s1.Atoi();

  if( XRangeDigit < 3 || i0 > 3 ) 
    XRangeMax = ( i0 + 1 ) * pow(10, XRangeDigit-1);
  else
    XRangeMax = i0 * pow(10, XRangeDigit-1) + (i1+1) * pow(10, XRangeDigit-2);
}
void SetXRange(double& XRangeMin, double& XRangeMax, vector<TH1D*> hists, TString Variable = "") {
  vector<TString> Variables = {"Eta","Phi"};
  for(unsigned iVar = 0; iVar < Variables.size(); iVar++) {
    if( Variable.Index( Variables[iVar] ) != kNPOS ) return;
  }

  XRangeMin = 0;
  for(unsigned i=0; i<hists.size(); i++) {
    XRangeMax = max(XRangeMax, hists[i]->FindLastBinAbove(0)*((TAxis*)hists[i]->GetXaxis())->GetBinWidth(0));
  }

  if( XRangeMax == hists[0]->GetXaxis()->GetXmax() ) return;

  if( XRangeMax < 1 ) return;

  TString s = TString::Itoa( XRangeMax , 10 );

  int XRangeDigit = s.Sizeof() - 1;

  TString s0 = s[0]; int i0 = s0.Atoi();
  TString s1 = s[1]; int i1 = s1.Atoi();

  if( XRangeDigit < 3 || i0 > 3 ) 
    XRangeMax = ( i0 + 1 ) * pow(10, XRangeDigit-1);
  else
    XRangeMax = i0 * pow(10, XRangeDigit-1) + (i1+1) * pow(10, XRangeDigit-2);
}


void SetRebinSize(double XRangeMin, double XRangeMax, double OriginalBinSize, double& NewBinSize, int& RebinSize) {
  vector<double> XRangeIndex = {0.1,0.2,0.5,1.,2.,5.,10.,20.,50.,100.,200.,500.,1000.,2000.,5000.,10000.,20000.,50000.,100000.,200000.,500000.};
  double XRange = XRangeMax - XRangeMin;
  for(int i=0; i<XRangeIndex.size(); i++) {
    if( XRange <= XRangeIndex[i] ) {
      NewBinSize = max(XRangeIndex[i]*0.01, OriginalBinSize);
      RebinSize = NewBinSize/OriginalBinSize;
      if( RebinSize/2. == (int)(RebinSize/2.) ) {
        RebinSize = RebinSize/2.;
        NewBinSize = NewBinSize/2.;
      }
      return;
    }
  }
  return;
}

double GetMaxY(TH1D *h1, TH1D *h2) {
  return max( h1->GetBinContent(h1->GetMaximumBin()),h2->GetBinContent(h2->GetMaximumBin()) );
}
double GetMaxY(TH1D *h1, TH1D *h2, TH1D *h3) {
  return max( GetMaxY(h1,h2),h3->GetBinContent(h3->GetMaximumBin()) );
}
double GetMaxY(TH1D *h1, TH1D *h2, TH1D *h3, TH1D *h4) {
  return max( GetMaxY(h1,h2),GetMaxY(h3,h4) );
}
double GetMaxY(TH1D *h1, TH1D *h2, TH1D *h3, TH1D *h4, TH1D *h5) {
  return max( GetMaxY(h1,h2),GetMaxY(h3,h4,h5) );
}
double GetMaxY(vector<TH1D*> hists) {
  double out=0;
  for(unsigned i=0; i<hists.size(); i++) {
    out = max(out, hists[i]->GetBinContent(hists[i]->GetMaximumBin()));
  }
  return out;
}

void SetYRange(TH1D* h1, TH1D* h2, const bool logy=false) {
  double y_min=0; double y_max=GetMaxY(h1,h2);
  double multiplier = 1.05;
  if(logy) {y_min=.5; multiplier = 3.;}
  y_max = y_max*multiplier;
  h1->GetYaxis()->SetRangeUser(y_min,y_max);
  h2->GetYaxis()->SetRangeUser(y_min,y_max);
}

void SetYRange(TH1D* h1, TH1D* h2, double y_min=0., double multiplier=1.05, const bool logy=false) {
  double y_max=GetMaxY(h1,h2);
  y_max = y_max*multiplier;
  h1->GetYaxis()->SetRangeUser(y_min,y_max);
  h2->GetYaxis()->SetRangeUser(y_min,y_max);
}

 


TString SetLumi(int Year) {
  TString out = "???";
  if( Year == 2016 ) out = "35.9";
  else if( Year == 2017 ) out = "41.5";
  return out;
}

int CountDigit(double d) {
  int out=0;
  int n = (int)d;
  while(n>0) {
    out++;
    n/=10;
  }
  return out;
}


TString XaxisIndex(TString Object, TString Variable, TString sign) {
  TString xTitle="";
  if( Variable.Index("Mass")!=kNPOS ) {
    xTitle += "M";
    if( sign.Index("OS")!=kNPOS )
      xTitle += "(#mu^{+},#mu^{-})";
    else if(sign.Index("SS")!=kNPOS )
      xTitle += "(#mu^{#pm},#mu^{#pm})";
    else
      xTitle += "(#mu,#mu)";
    xTitle += " [GeV]";
  }
  else if( Variable.Index("Eta")!=kNPOS ) {
    xTitle += "#eta";
    if( Object.Index("0")!=kNPOS )
      xTitle += "(#mu_{Leading})";
    else if( Object.Index("1")!=kNPOS )
      xTitle += "(#mu_{Subleading})";
  }
  else if( Variable.Index("Pt")!=kNPOS ) {
    xTitle += "p_{T}";
    if( Object.Index("Lepton")!=kNPOS ) {
      if( Object.Index("0")!=kNPOS )
        xTitle += "(#mu_{Leading})";
      else if( Object.Index("1")!=kNPOS )
        xTitle += "(#mu_{Subleading})";
    }
    else if( Object.Index("Jet")!=kNPOS ) {
      if( Variable.Index("PtminusSumMuPt")!=kNPOS ) {
        xTitle += "(j-#mu_{0}-#mu_{1})";
      }
      else if( Variable.Index("Pt")!=kNPOS ) {
        xTitle += "(j)";
      }
    }
    else if( Object.Index("MET")!=kNPOS ) {
      xTitle = "MET "+xTitle;
    }
    xTitle += " [GeV]";
  }
  else if( Variable.Index("DeltaR")!=kNPOS ) {
    xTitle += "#DeltaR";
    if( sign.Index("OS")!=kNPOS )
      xTitle += "(#mu^{+},#mu^{-})";
    else if(sign.Index("SS")!=kNPOS )
      xTitle += "(#mu^{#pm},#mu^{#pm})";
    else
      xTitle += "(#mu,#mu)";
  }
  else if( Variable.Index("Phi")!=kNPOS ) {
    xTitle += "#phi(";
    if( Object.Index("0")!=kNPOS )
      xTitle += "#mu_{Leading})";
    else
      xTitle += "#mu_{Subleading})";
  }
  else if( Variable.Index("Iso")!=kNPOS ) {
    xTitle += "Isolation";
  }
  return xTitle;
}

void sort_MC_Map( vector< tuple<int,double,TString,TString,TString,TString,Color_t,int,double,Color_t,int,double,Color_t,int> >& map, vector<TH1D*>& hists ){
  vector< tuple<int,double,TString,TString,TString,TString,Color_t,int,double,Color_t,int,double,Color_t,int> > map_;
  vector<TH1D*> hists_;
  int Ntype = 0;
  for(unsigned i=0; i<map.size(); i++) {
    if( Ntype < get<0>(map[i]) ) Ntype = get<0>(map[i]);
  }
  for(unsigned i=0; i<Ntype+1; i++) {
    for(unsigned j=0; j<map.size(); j++) {
      if(get<0>(map[j])==i) {
        map_.push_back( map[j] );
        hists_.push_back( hists[j] );
      }
    }
  }
  map=map_; hists=hists_;
}

double GetSumw( TH1D* hist, double lowx, double highx ) {
  double sumw=0; double err=0;
  for(unsigned i=0; i<hist->GetNbinsX(); i++) {
    double binloc=((TAxis*)hist->GetXaxis())->GetBinCenter(i);
    if( lowx < binloc && binloc < highx ) {
      sumw += hist->GetBinContent(i);
      err = sqrt( err*err + hist->GetBinError(i)*hist->GetBinError(i) );
      //cout<<"bin loc = "<<hist->GetBinCenter(i)<<", content: "<<hist->GetBinContent(i)<<", err: "<<hist->GetBinError(i)<<endl;
    }
  }
  cout<<"bin = ["<<lowx<<", "<<highx<<"]"<<", content: "<<sumw<<" entry: "<<sqrt(sumw)/err*sumw<<", err: "<<err<<", r = "<<err/sqrt(sumw)<<endl;
  return sumw;
}

TH1D* make_ratio_plot(TH1D* h1, TH1D* h2) {
  TH1D* result = (TH1D*)h1->Clone("");
  result->Divide(h2);
  result->GetYaxis()->SetRangeUser(0.,2.);
  return result;
}

void draw_comparison_plot(map<TString, TH1D*> hist, const vector<TString> categories, map<TString, Color_t> color, TString title, TString& x_label, TString& y_label, bool logx = false, bool logy = false) {

  TH1D* hist_mc = nullptr;
  THStack* stack_mc = new THStack("","");
 
  hist["data"]=SetDataStyle(hist["data"], 20, .4, color["data"]);
  hist["data"]->GetXaxis()->SetRangeUser(1.,100.);
  for(const auto& category : categories) {
    if(hist[category]==nullptr) {
      cerr<<"Error: The histogram ("<<category<<") is null."<<endl;
      continue;
    }
    if(category=="data") continue;
    hist[category]=SetMCStyle(hist[category], color[category], color[category], 1.);
    hist[category]->GetXaxis()->SetRangeUser(1.,100.);
    if(hist_mc==nullptr) hist_mc = (TH1D*)hist[category]->Clone("");
    else hist_mc->Add(hist[category]);
    stack_mc->Add(hist[category]);
  }

  TH1D* hist_ratio = make_ratio_plot(hist["data"], hist_mc);

  SetYRange(hist["data"], hist_mc, logy);

  TCanvas *canvas = new TCanvas("canvas","",1000,1000);
  canvas->cd();

  auto pad1 = new TPad("","",0,.26,1,1,0);
  pad1->SetLeftMargin(.1);
  if(logx) pad1->SetLogx();
  if(logy) pad1->SetLogy();
  pad1->Draw();

  auto pad2 = new TPad("","",0,0,1,.33,0);
  pad2->SetLeftMargin(.1);
  pad2->SetBottomMargin(.3);
  if(logx) pad2->SetLogx();
  pad2->Draw();

  setTDRStyle();

  writeExtraText = true;
  //extraText = "Preliminary";
  lumi_13TeV = "35.9 fb^{-1}";
  lumi_sqrtS = "13 TeV";

  int iPeriod = 4;
  // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 4=13TeV, 0=free form (uses lumi_sqrtS)
  int iPos = 0;
  // second parameter in example_plot is iPos, which drives the position of the CMS logo in the plot
  // iPos=11 : top-left, left-aligned
  // iPos=33 : top-right, right-aligned
  // iPos=22 : center, centered
  // mode generally :
  //   iPos = 10*(alignement 1/2/3) + position (1/2/3 = left/center/right)

  //CMS_lumi( pad1 , iPeriod , iPos );

  // --- Upper Pad --- //
  pad1->cd();

  // Draw preliminary plots //
  hist["data"]->Draw("PE");

  if(logx) x_label = x_label + " (log)";
  if(logy) y_label = y_label + " (log)";

  hist["data"] -> GetYaxis()->SetTitle(y_label);
  hist["data"] -> GetYaxis()->SetTitleSize(0.045);
  hist["data"] -> GetYaxis()->SetLabelSize(0.040);
  hist["data"] -> GetYaxis()->SetTitleOffset(1.0);

  hist["data"] -> GetXaxis()->SetTitle(x_label);
  hist["data"] -> GetXaxis()->SetTitleSize(0);
  hist["data"] -> GetXaxis()->SetTitleOffset(0.7);
  hist["data"] -> GetXaxis()->SetLabelSize(0.025);

  TLegend* legend;
  legend = new TLegend(0.69,0.82-0.07*(categories.size()-1),0.90,0.90);
  legend->AddEntry(hist["data"], "data", "lep");
  for (int i = categories.size() - 1; i >= 0; --i) {
    const auto& category = categories[i];
		if(category=="data") continue;
    legend->AddEntry(hist[category], category, "f");
  }
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->Draw();

  stack_mc->Draw("HIST SAME");
  hist["data"]->SetStats(0);
  hist["data"]->Draw("PE SAME");

  pad1->SetGridx(); pad1->SetGridy();
  pad1->RedrawAxis();
  CMS_lumi( pad1 , iPeriod , iPos );

  // --- Lower Pad --- //
  pad2->cd();

  //hist_ratio->Divide(hist_mc);

  hist_ratio->GetXaxis()->SetTitle(x_label);
  hist_ratio->GetXaxis()->SetTitleOffset(1.1);
  hist_ratio->GetXaxis()->SetTitleSize(0.1);
  hist_ratio->GetXaxis()->SetLabelSize(0.09);
  hist_ratio->GetXaxis()->SetLabelOffset(0.01);

  hist_ratio->GetYaxis()->SetTitle("Data / MC");
  hist_ratio->GetYaxis()->SetTitleSize(0.095);
  hist_ratio->GetYaxis()->SetTitleOffset(0.40);
  hist_ratio->GetYaxis()->SetLabelSize(0.075);

  hist_ratio->SetStats(0);
  hist_ratio->Draw("P");

  pad2->SetGridx(); pad2->SetGridy();

  TF1 *hline = new TF1("","1",hist["data"]->GetXaxis()->GetXmin(),hist["data"]->GetXaxis()->GetXmax());
  hline->SetLineColor(kRed);
  hline->Draw("SAME");

  // --- Save Plot --- //
  canvas->SaveAs("../outputs/"+title+".png");
}

#endif
