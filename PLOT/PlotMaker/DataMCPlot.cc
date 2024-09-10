#include "../PlotHelper/PlotterCore.h"

void DataMCPlot(const TString target = "mass_den") {
  // --- Initial Settings --- //
  TString x_name = "M(#mu,#mu)";
  TString y_name = "Events";
  const Int_t nxbin = 17;
  const Double_t xbin[nxbin+1] = {0.,1.,2.,2.9,3.3,4.,5.,6.,8.,10.,15.,20.,25.,30.,35.,40.,50.,100.};
  double x_min = 1.0; double x_max = 100.0;
  bool logx = true, logy = true;

  double lumi=36.3*1000; // 35.9 fb-1

  bool validation_mode = 0;
  //TString input_path = "../inputs/";
  TString input_path = "../inputs/240827_PfJet450/";

  //            sample  color   title   xsec   Nevent
  std::vector< std::tuple<TString, TString, double, Color_t> > samples = { 
    {"data", "JetHT",           1.                              ,  1}, // data
    {"ST",   "ST_antitop",      lumi * 19.28       / 8681541.   , 46}, // 38.06 pb
    {"ST",   "ST_top",          lumi * 19.28       / 8681495.   , 46}, // 38.09
    {"tt",   "TTLL",            lumi * 88.29       / 79140880.  , 46}, // 76.7
    {"tt",   "TTLJ",            lumi * 365.34      / 152669800. , 46}, // 320.1
    //{"QCD",  "QCD_Pt-15to20", lumi * 3819570     / 4141251.   ,  9},
    //{"QCD",  "QCD_Pt-20to30", lumi * 2960198.4   / 31878740.  ,  9},
    {"QCD",  "QCD_Pt-30to50",   lumi * 1652471.46  / 29954815.  ,  9}, // 1662000
    {"QCD",  "QCD_Pt-50to80",   lumi * 437504.1    / 19662175.  ,  9}, // 452200
    {"QCD",  "QCD_Pt-80to120",  lumi * 106033.6648 / 23705386.  ,  9}, // 106500
    {"QCD",  "QCD_Pt-120to170", lumi * 25190.51514 / 7897731.   ,  9}, // 25700
    {"QCD",  "QCD_Pt-170to300", lumi * 8654.49315  / 17350231.  ,  9}, // 8683
    {"QCD",  "QCD_Pt-300to470", lumi * 797.35269   / 49005976.  ,  9}, // 797.3
    {"QCD",  "QCD_Pt-470to600", lumi * 79.02553776 / 19489276.  ,  9}, // 79.25
    {"QCD",  "QCD_Pt-600to800", lumi * 25.09505908 / 9981311.   ,  9}, // 25.25
    {"QCD",  "QCD_Pt-800to1000",lumi * 4.707368272 / 19940747.  ,  9}, // 4.723
    {"QCD",  "QCD_Pt-1000toInf",lumi * 1.62131692  / 13628219.  ,  9} // 1.613
  };  

  vector<TString> categories;
  map<TString, TH1D*> hist;
  map<TString, Color_t> color;
  map<TString, double> sumw;

  for(const auto& sample : samples) {
    TString category = get<0>(sample);
    if(categories.empty() || category != categories.back()) {
      categories.push_back(category);
      hist[category] = new TH1D("hist_"+category, "", nxbin, xbin);
      color[category] = get<3>(sample);
      sumw[category] = 0.;
    }
  }

  for(const auto& sample : samples) {
    TString category = get<0>(sample);
    TString sample_name = get<1>(sample);
    double wgt = get<2>(sample);

    TFile* file = TFile::Open(input_path+"/"+sample_name+".root");
    TH1D* tmp = dynamic_cast<TH1D*>(file->Get("demo/"+target));
    if(!tmp) {
      cerr<<"Error: Could not get histogram '"<<target<<"'"<<endl;
      continue;
    }

    if(validation_mode) {
      cout<<endl<<" ++ "<<category<<" (before wgt)"<<endl;
      for(int i=0; i<nxbin; i++) {
        cout<<"("<<xbin[i]<<", "<<tmp->GetBinContent(i+1)<<")  ";
      }
      cout<<endl;
    }
 
    tmp->Sumw2();
    tmp->Scale(wgt);

    if(validation_mode) {  
      cout<<"             (after wgt)"<<endl;
      for(int i=0; i<nxbin; i++) {
        cout<<"("<<xbin[i]<<", "<<tmp->GetBinContent(i+1)<<")  ";
      }
      cout<<endl;
    }

    hist[category]->Add(tmp);
    delete tmp; delete file;
  }

  int bin_low = hist["data"]->FindBin(1.);
  int bin_high = hist["data"]->FindBin(2.9)-1;

  for(const auto& category : categories) {
    sumw[category] = hist[category]->Integral(bin_low, bin_high);
  }
  double QCD_rescale = (sumw["data"]-sumw["ST"]-sumw["tt"])/sumw["QCD"];

  //if(validation_mode) 
  cout<<endl<<"New QCD Rescale Factor within 1 <= x < 2.9: "<<QCD_rescale<<endl;

  hist["QCD"]->Scale(QCD_rescale);

  draw_comparison_plot(hist, categories, color, "test", x_name, y_name, logx, logy);
  return;
}
