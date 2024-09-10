#include "../PlotHelper/PlotterCore.h"
#include "../PlotHelper/tdrstyle.C"
#include "../PlotHelper/CMS_lumi.C"

void hlteff2d() {
  // --- Plot settings ---
  const TString x_name = "M(#mu,#mu)";
  const TString y_name = "p_T(#mu_{lead})";
  const Int_t nxbin = 17; const Int_t nybin = 10;
  const Double_t xbin[nxbin+1] = {0.,1.,2.,2.9,3.3,4.,5.,6.,8.,10.,15.,20.,25.,30.,35.,40.,50.,100.};
  const Double_t ybin[nybin+1] = {31., 32., 33., 34., 35., 40., 50.,60., 80., 100., 200.};
  const double x_min = 1.0, x_max = 100.0;
  const double y_min = 31.0, y_max = 200.0;
  const bool logx = true, logy = true;

  // --- File paths and object names ---
  const TString den_name = "demo/eff_den";
  const TString num_name = "demo/eff_num";

  double lumi=35.9;

  //            sample  color   title   xsec   Nevent
  std::vector< std::tuple<TString, Color_t, TString, double> > samples = {
    {"data", 1, "JetHT",           1.                              }, // data
    {"ST",  46, "ST_antitop",      lumi * 19.28       / 8681541.   }, // 38.06
    {"ST",  46, "ST_top",          lumi * 19.28       / 8681495.   }, // 38.09
    {"tt",  46, "TTLL",            lumi * 88.29       / 79140880.  }, // 76.7
    {"tt",  46, "TTLJ",            lumi * 365.34      / 152669800. }, // 320.1
    //{"QCD",  9, "QCD_Pt-15to20", lumi * 3819570     / 4141251.   },
    //{"QCD",  9, "QCD_Pt-20to30", lumi * 2960198.4   / 31878740.  },
    {"QCD",  9, "QCD_Pt-30to50",   lumi * 1652471.46  / 29954815.  }, // 1662000
    {"QCD",  9, "QCD_Pt-50to80",   lumi * 437504.1    / 19662175.  }, // 452200
    {"QCD",  9, "QCD_Pt-80to120",  lumi * 106033.6648 / 23705386.  }, // 106500
    {"QCD",  9, "QCD_Pt-120to170", lumi * 25190.51514 / 7897731.   }, // 25700
    {"QCD",  9, "QCD_Pt-170to300", lumi * 8654.49315  / 17350231.  }, // 8683
    {"QCD",  9, "QCD_Pt-300to470", lumi * 797.35269   / 49005976.  }, // 797.3
    {"QCD",  9, "QCD_Pt-470to600", lumi * 79.02553776 / 19489276.  }, // 79.25
    {"QCD",  9, "QCD_Pt-600to800", lumi * 25.09505908 / 9981311.   }, // 25.25
    {"QCD",  9, "QCD_Pt-800to1000",lumi * 4.707368272 / 19940747.  }, // 4.723
    {"QCD",  9, "QCD_Pt-1000toInf",lumi * 1.62131692  / 13628219.  } // 1.613
  };

  // --- Loop for MC sets ---
  vector<TString> categories;
  map<TString, TH2D*> hist_den;
  map<TString, TH2D*> hist_num;
  hist_den["inclusive"] = new TH2D("hist_den_inclusive", "", nxbin, xbin, nybin, ybin);
  hist_num["inclusive"] = new TH2D("hist_num_inclusive", "", nxbin, xbin, nybin, ybin);
  for(const auto& sample : samples) {
    TString category = get<0>(sample);
    if(categories.empty() || category != categories.back()) {
      categories.push_back(category);
      hist_den[category] = new TH2D("hist_den_" + category, "", nxbin, xbin, nybin, ybin);
      hist_num[category] = new TH2D("hist_num_" + category, "", nxbin, xbin, nybin, ybin);
    }
  }

  for(const auto& sample : samples) {
    TString category = get<0>(sample);
    TString sample_name = get<2>(sample);
    double wgt = get<3>(sample);

    TFile* file = TFile::Open("../inputs/"+sample_name+".root");

    TH2D* tmp_den = dynamic_cast<TH2D*>(file->Get(den_name));
    if (!tmp_den) std::cerr << "Error: Could not get histogram '" << den_name << "'" << std::endl;
    TH2D* tmp_num = dynamic_cast<TH2D*>(file->Get(num_name));
    if (!tmp_num) std::cerr << "Error: Could not get histogram '" << num_name << "'" << std::endl;
    tmp_den->Sumw2();
    tmp_num->Sumw2();

    tmp_den->Scale(wgt);
    tmp_num->Scale(wgt);

    hist_den[category]->Add(tmp_den);
    hist_num[category]->Add(tmp_num);
    if(category!="data") {
      hist_den["inclusive"]->Add(tmp_den);
      hist_num["inclusive"]->Add(tmp_num);
    }
    delete tmp_den; delete tmp_num; delete file; 
  }

  // --- Loop for each category ---
  for(const auto& category : categories) {
    drawTH2D(hist_den[category], "Before HLT_Mu30_TkMu11 ("+category+")", x_name, y_name, logx, logy);
    drawTH2D(hist_num[category], "After HLT_Mu30_TkMu11 ("+category+")", x_name, y_name, logx, logy);

    TH2D* eff = new TH2D(*hist_num[category]);
    eff->Divide(hist_den[category]);
    drawTH2D(eff, "#epsilon ("+category+")", x_name, y_name, logx, logy);

    delete eff;
  }

  drawTH2D(hist_den["inclusive"], "Before HLT_Mu30_TkMu11 (MC)", x_name, y_name, logx, logy);
  drawTH2D(hist_num["inclusive"], "After HLT_Mu30_TkMu11 (MC)", x_name, y_name, logx, logy);
}
