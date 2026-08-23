//hijing_wrapper.cpp
#include "wrapper.h"
#include "WrapperQCorrEvent.h"
//---------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
#include <string>
#include <math.h>
#include <random>
#include <chrono>
#include <csignal>
#include <algorithm>
#include <set>
#include <unordered_map>
#include <functional>
//---------------------------------------------------------------
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TMath.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TUnixSystem.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include "TClassTable.h"
#include "TSystem.h"
#include "TDirectory.h"
#include <TClonesArray.h>
#include <TRandom3.h>
//--------------------------------------------------------------- 
//---------------------------------------------------------------
#ifdef DEBUG
#define LogDebug(message) std::cout << "[DEBUG] " << message << std::endl
#else
#define LogDebug(message)
#endif
//---------------------------------------------------------------
//---------------------------------------------------------------
bool FMpTW = false;

using namespace std;
const double RefPTcut = 2.0;
const double etaMin = -1.0;
const double etaMax = 1.0;
const double ptMin = 0.2;
const double ptMax = 5.0;

const int kNCentBins = 10;
const int kNEtaBins = 20;
const int kNPhiBins = 36;

const int kNPtaBins = 20;
const double pTbin[20] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 5.0};
//---------------------------------------------------------------
//---------------------------------------------------------------
//---------------------------------------------------------------
//---------------------------------------------------------------
const double FMtoMM = 1.0e-12;
const double Pi     = 3.141592653589793;
//---------------------------------------------------------------
bool GetMpT = true;  // true; //false;
//---------------------------------------------------------------
extern "C" {
  void hijset_(float* efrm, char* frame, char* proj, char* targ, int* iap, int* izp, int* iat, int* izt);
}
#define HIJSET hijset_

extern "C" {
  void hijing_(char* FRAME, float* BMIN, float* BMAX);
}
#define HIJING hijing_

extern "C" {
  extern struct{ 
    int natt;
    int eatt;
    int jatt;
    int nt;
    int np;
    int n0;
    int n01;
    int n10;
    int n11;
  }himain1_;
}
#define himain1 himain1_

extern "C" {
  extern struct{ 
    int   katt[4][130000];
    float patt[4][130000];
    float vatt[4][130000];
  }himain2_;
}
#define himain2 himain2_

extern "C" {
  extern struct{
    float vatt[4][130000];
  }himain3_;
}
#define himain3 himain3_

extern "C" {
  extern struct{ 
    float  hipr1[100];
    int    ihpr2[50];
    float  hint1[100];
    int    ihnt2[50];
  }hiparnt_;
}
#define hiparnt hiparnt_

extern "C" {
  extern struct{
    int   nfp[15][300];
    float pp[15][300];
    int   nft[15][300];
    float pt[15][300];   
  }histrng_;
}
#define histrng histrng_

extern "C" {
  extern struct{
    float yp[300][3];
    float yt[300][3];
  }hijcrdn_;
}
#define hijcrdn hijcrdn_

extern "C" {
  extern struct{
    int nseed;
  }ranseed_;
}
#define ranseed ranseed_
//---------------------------------------------------------------
//---------------------------------------------------------------
int getNparCentID(float nX) {
  float Npar[11] = {0.0,5.5,12.5,23.5,41.5,66.5,101.5,146.5,206.5,283.5,1000.0};

  for (int i = 0; i < 10; i++) {
    if (nX >= Npar[i] && nX < Npar[i+1]) {
      return i;
    }
  }
  return -999; // Return last bin if value is higher than all boundaries
}
//---------------------------------------------------------------
// ---------- HistoAutoFile ----------
HistoAutoFile::HistoAutoFile(const TString& fname, const char* opt)
  : mFile(TFile::Open(fname, opt))
{
  if (!mFile || mFile->IsZombie()) {
    throw std::runtime_error(Form("Cannot open output file '%s'", fname.Data()));
  }
  mFile->cd();
  // ensure histograms attach to the current directory (this file)
  TH1::AddDirectory(kTRUE);
}

HistoAutoFile::~HistoAutoFile() {
  if (mFile && mFile->IsOpen()) {
    // writes all objects attached to this file's directories
    mFile->Write();
    mFile->Close();
  }
}

TDirectory* HistoAutoFile::Mkdir(const char* name) {
  if (!mFile) return nullptr;
  TDirectory* d = mFile->mkdir(name);
  if (d) d->cd();
  return d;
}

// ---------- Wrapper ----------
void Wrapper::OpenOutput(const TString& fname) {
  mOut = std::make_unique<HistoAutoFile>(fname);
}

void Wrapper::CloseOutput() {
  // Destroying HistoAutoFile triggers Write() and Close() automatically
  mOut.reset();
}

TDirectory* Wrapper::MakeDir(const char* name) {
  return mOut ? mOut->Mkdir(name) : nullptr;
}
//---------------------------------------------------------------
//---------------------------------------------------------------
struct Params {
  int   MODE{};
  int   DECAY{};
  int   FRAMEID{};
  int   N_EVENT{};
  int   IAP{};
  int   IZP{};
  int   IAT{};
  int   IZT{};
  float EFRM{};
  float BMIN{};
  float BMAX{};
};
//---------------------------------------------------------------
static inline std::string upper(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::toupper(c); });
  return s;
}
//---------------------------------------------------------------
bool load_params_from_stream(std::istream& in, Params& p, std::string& err) {
  std::set<std::string> seen;
  std::unordered_map<std::string, std::function<void(const std::string&)>> set = {
    {"MODE",    [&](const std::string& v){ p.MODE    = std::stoi(v); seen.insert("MODE"); }},
    {"DECAY",   [&](const std::string& v){ p.DECAY   = std::stoi(v); seen.insert("DECAY"); }},
    {"FRAMEID", [&](const std::string& v){ p.FRAMEID = std::stoi(v); seen.insert("FRAMEID"); }},
    {"N_EVENT", [&](const std::string& v){ p.N_EVENT = std::stoi(v); seen.insert("N_EVENT"); }},
    {"IAP",     [&](const std::string& v){ p.IAP     = std::stoi(v); seen.insert("IAP"); }},
    {"IZP",     [&](const std::string& v){ p.IZP     = std::stoi(v); seen.insert("IZP"); }},
    {"IAT",     [&](const std::string& v){ p.IAT     = std::stoi(v); seen.insert("IAT"); }},
    {"IZT",     [&](const std::string& v){ p.IZT     = std::stoi(v); seen.insert("IZT"); }},
    {"EFRM",    [&](const std::string& v){ p.EFRM    = std::stof(v); seen.insert("EFRM"); }},
    {"BMIN",    [&](const std::string& v){ p.BMIN    = std::stof(v); seen.insert("BMIN"); }},
    {"BMAX",    [&](const std::string& v){ p.BMAX    = std::stof(v); seen.insert("BMAX"); }}
  };

  std::string line;
  size_t lineno = 0;
  while (std::getline(in, line)) {
    ++lineno;
    // strip comments
    if (auto pos = line.find('#'); pos != std::string::npos) line.erase(pos);
    // trim
    auto notspace = [](int ch){ return !std::isspace(ch); };
    line.erase(line.begin(), std::find_if(line.begin(), line.end(), notspace));
    line.erase(std::find_if(line.rbegin(), line.rend(), notspace).base(), line.end());
    if (line.empty()) continue;

    std::string key, val;
    if (auto eq = line.find('='); eq != std::string::npos) {
      key = line.substr(0, eq);
      val = line.substr(eq + 1);
    } else {
      std::istringstream iss(line);
      iss >> key >> val;
      if (val.empty()) { err = "Missing value at line " + std::to_string(lineno); return false; }
    }
    key = upper(key);
    // trim val again (in case of "KEY = value")
    val.erase(val.begin(), std::find_if(val.begin(), val.end(), notspace));
    val.erase(std::find_if(val.rbegin(), val.rend(), notspace).base(), val.end());

    auto it = set.find(key);
    if (it == set.end()) {
      err = "Unknown key '" + key + "' at line " + std::to_string(lineno);
      return false;
    }
    try {
      it->second(val);
    } catch (const std::exception& e) {
      err = "Bad value for '" + key + "' at line " + std::to_string(lineno) + ": " + e.what();
      return false;
    }
  }

  // verify all required keys present
  for (auto& kv : set) {
    if (!seen.count(kv.first)) {
      err = "Missing required key: " + kv.first;
      return false;
    }
  }
  return true;
}

bool load_params_from_file(const std::string& path, Params& p, std::string& err) {
  if (path == "-") {
    return load_params_from_stream(std::cin, p, err);
  }
  std::ifstream fin(path);
  if (!fin) { err = "Failed to open config file: " + path; return false; }
  return load_params_from_stream(fin, p, err);
}
//---------------------------------------------------------------
void drop_fraction_from_tclones(TClonesArray* arr, double dropFrac, UInt_t seed=12345)
{
  if (!arr) return;
  TRandom3 rng(seed);

  // Iterate backwards so indices stay valid while removing
  for (int i = arr->GetEntriesFast() - 1; i >= 0; --i) {
    if (rng.Uniform() < dropFrac) arr->RemoveAt(i);
  }
  arr->Compress(); // physically remove holes
}
//---------------------------------------------------------------
//---------------------------------------------------------------
//---------------------------------------------------------------
int getPtaBin(double pt) {
  int bin = -1;
  for (int i = 0; i < kNPtaBins-1; ++i) {
    if (pt >= pTbin[i] && pt < pTbin[i + 1]) bin = i;
  }
  return bin;
}
//---------------------------------------------------------------
int getEtaBin(double eta) {
  if (eta < etaMin || eta > etaMax) return -1;

  const double width = (etaMax - etaMin) / static_cast<double>(kNEtaBins);
  int bin = static_cast<int>((eta - etaMin) / width);

  if (bin < 0) return -1;
  if (bin >= kNEtaBins) bin = kNEtaBins - 1;
  return bin;
}
//---------------------------------------------------------------
//---------------------------------------------------------------
void Wrapper::BookHistos() {  
  //-------------------------------------------------------------
  TH1::AddDirectory(kTRUE);
  //-------------------------------------------------------------
  if (!mOut || !mOut->File() || !mOut->File()->IsOpen()) return;
  mOut->cd(); 
  //-------------------------------------------------------------
  char  hname[2000];
  char  histname[4000];
  char  histname1[4000];
  char  histname2[4000];
  char  histname3[4000];
  char  histname4[4000];
  char  histname5[4000];
  char  histname6[4000];
  char  histname7[4000];
  char  histname8[4000];
  //-------------------------------------------------------------
  MultA       = new TH1F ("MultA","MultA", 5000, 0.0, 5000.0) ;
  MultA->Sumw2();

  ImbaA       = new TH1F ("ImbaA","ImbaA", 5000 , 0.0, 20.0   ) ;
  ImbaA->Sumw2();

  NparA       = new TH1F ("NparA","NparA", 600 , 0.0, 600.0 ) ;
  NparA->Sumw2();
  //-------------------------------------------------------------
  for (int i = 0; i < 10; ++i) {
    sprintf(histname,"Eta_%d",i);
    Etah[i] = new TH1D(histname,histname,44,-1,1);
    Etah[i]->Sumw2();

    sprintf(histname,"Phi_%d",i);
    Phih[i] = new TH1D(histname,histname,50,-3.14159,3.14159);
    Phih[i]->Sumw2();

    sprintf(histname,"PT_%d",i);
    PTh[i] = new TH1D(histname,histname,1000,0,150);
    PTh[i]->Sumw2();
  }
  //-------------------------------------------------------------
  NChvsCent = new TProfile("NChvsCent", "NChvsCent", kNCentBins, -0.5, double(kNCentBins) - 0.5);
  NChvsCent->Sumw2();

  NpartvsCent = new TProfile("NpartvsCent", "NpartvsCent", kNCentBins, -0.5, double(kNCentBins) - 0.5);
  NpartvsCent->Sumw2();
  //-------------------------------------------------------------
  TP_d00_c11t0           =new TProfile("TP_d00_c11t0"           ,"TP_d00_c11t0"           ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d00_c11t0           ->Sumw2();
  TP_d00_c22t0           =new TProfile("TP_d00_c22t0"           ,"TP_d00_c22t0"           ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d00_c22t0           ->Sumw2();
  TP_d00_c33t0           =new TProfile("TP_d00_c33t0"           ,"TP_d00_c33t0"           ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d00_c33t0           ->Sumw2();
  TP_d00_c44t0           =new TProfile("TP_d00_c44t0"           ,"TP_d00_c44t0"           ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d00_c44t0           ->Sumw2();
  TP_d00_c55t0           =new TProfile("TP_d00_c55t0"           ,"TP_d00_c55t0"           ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d00_c55t0           ->Sumw2();
  TP_d00_c66t0           =new TProfile("TP_d00_c66t0"           ,"TP_d00_c66t0"           ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d00_c66t0           ->Sumw2();
  TP_d000_c211t0         =new TProfile("TP_d000_c211t0"         ,"TP_d000_c211t0"         ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000_c211t0         ->Sumw2();
  TP_d000_c312t0         =new TProfile("TP_d000_c312t0"         ,"TP_d000_c312t0"         ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000_c312t0         ->Sumw2();
  TP_d000_c413t0         =new TProfile("TP_d000_c413t0"         ,"TP_d000_c413t0"         ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000_c413t0         ->Sumw2();
  TP_d000_c514t0         =new TProfile("TP_d000_c514t0"         ,"TP_d000_c514t0"         ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000_c514t0         ->Sumw2();
  TP_d000_c422t0         =new TProfile("TP_d000_c422t0"         ,"TP_d000_c422t0"         ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000_c422t0         ->Sumw2();
  TP_d000_c523t0         =new TProfile("TP_d000_c523t0"         ,"TP_d000_c523t0"         ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000_c523t0         ->Sumw2();
  TP_d000_c624t0         =new TProfile("TP_d000_c624t0"         ,"TP_d000_c624t0"         ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000_c624t0         ->Sumw2();
  TP_d000_c633t0         =new TProfile("TP_d000_c633t0"         ,"TP_d000_c633t0"         ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000_c633t0         ->Sumw2();
  //---------------------------------------------------------------------------------------------------------------------------
  TP_d0000_c1111t0       =new TProfile("TP_d0000_c1111t0"       ,"TP_d0000_c1111t0"       ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c1111t0       ->Sumw2();
  TP_d0000_c2222t0       =new TProfile("TP_d0000_c2222t0"       ,"TP_d0000_c2222t0"       ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c2222t0       ->Sumw2();
  TP_d0000_c3333t0       =new TProfile("TP_d0000_c3333t0"       ,"TP_d0000_c3333t0"       ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c3333t0       ->Sumw2();
  TP_d0000_c1212t0       =new TProfile("TP_d0000_c1212t0"       ,"TP_d0000_c1212t0"       ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c1212t0       ->Sumw2();
  TP_d0000_c1313t0       =new TProfile("TP_d0000_c1313t0"       ,"TP_d0000_c1313t0"       ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c1313t0       ->Sumw2();
  TP_d0000_c1414t0       =new TProfile("TP_d0000_c1414t0"       ,"TP_d0000_c1414t0"       ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c1414t0       ->Sumw2();
  TP_d0000_c1515t0       =new TProfile("TP_d0000_c1515t0"       ,"TP_d0000_c1515t0"       ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c1515t0       ->Sumw2();
  TP_d0000_c2323t0       =new TProfile("TP_d0000_c2323t0"       ,"TP_d0000_c2323t0"       ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c2323t0       ->Sumw2();
  TP_d0000_c2424t0       =new TProfile("TP_d0000_c2424t0"       ,"TP_d0000_c2424t0"       ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c2424t0       ->Sumw2();
  TP_d0000_c2525t0       =new TProfile("TP_d0000_c2525t0"       ,"TP_d0000_c2525t0"       ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c2525t0       ->Sumw2();
  TP_d0000_c3434t0       =new TProfile("TP_d0000_c3434t0"       ,"TP_d0000_c3434t0"       ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c3434t0       ->Sumw2();
  TP_d0000_c3535t0       =new TProfile("TP_d0000_c3535t0"       ,"TP_d0000_c3535t0"       ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c3535t0       ->Sumw2();
  TP_d0000_c3111t0       =new TProfile("TP_d0000_c3111t0"       ,"TP_d0000_c3111t0"       ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c3111t0       ->Sumw2();
  TP_d0000_c6222t0       =new TProfile("TP_d0000_c6222t0"       ,"TP_d0000_c6222t0"       ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c6222t0       ->Sumw2();
  TP_d0000_c3324t0       =new TProfile("TP_d0000_c3324t0"       ,"TP_d0000_c3324t0"       ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c3324t0       ->Sumw2();
  TP_d0000_c4435t0       =new TProfile("TP_d0000_c4435t0"       ,"TP_d0000_c4435t0"       ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c4435t0       ->Sumw2();
  TP_d0000_c2534t0       =new TProfile("TP_d0000_c2534t0"       ,"TP_d0000_c2534t0"       ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c2534t0       ->Sumw2();
  //---------------------------------------------------------------------------------------------------------------------------
  TP_d00000_c21111t0     =new TProfile("TP_d00000_c21111t0"     ,"TP_d00000_c21111t0"     ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d00000_c21111t0     ->Sumw2();
  TP_d00000_c42222t0     =new TProfile("TP_d00000_c42222t0"     ,"TP_d00000_c42222t0"     ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d00000_c42222t0     ->Sumw2();
  TP_d00000_c43322t0     =new TProfile("TP_d00000_c43322t0"     ,"TP_d00000_c43322t0"     ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d00000_c43322t0     ->Sumw2();
  TP_d00000_c52232t0     =new TProfile("TP_d00000_c52232t0"     ,"TP_d00000_c52232t0"     ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d00000_c52232t0     ->Sumw2();
  TP_d00000_c53332t0     =new TProfile("TP_d00000_c53332t0"     ,"TP_d00000_c53332t0"     ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d00000_c53332t0     ->Sumw2();
  TP_d00000_c22435t0     =new TProfile("TP_d00000_c22435t0"     ,"TP_d00000_c22435t0"     ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d00000_c22435t0     ->Sumw2();
  TP_d00000_c23344t0     =new TProfile("TP_d00000_c23344t0"     ,"TP_d00000_c23344t0"     ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d00000_c23344t0     ->Sumw2();
  TP_d00000_c33345t0     =new TProfile("TP_d00000_c33345t0"     ,"TP_d00000_c33345t0"     ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d00000_c33345t0     ->Sumw2();
  TP_d00000_c24455t0     =new TProfile("TP_d00000_c24455t0"     ,"TP_d00000_c24455t0"     ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d00000_c24455t0     ->Sumw2();
  TP_d00000_c33455t0     =new TProfile("TP_d00000_c33455t0"     ,"TP_d00000_c33455t0"     ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d00000_c33455t0     ->Sumw2();
  TP_d00000_c22233t0     =new TProfile("TP_d00000_c22233t0"     ,"TP_d00000_c22233t0"     ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d00000_c22233t0     ->Sumw2();
  TP_d00000_c62233t0     =new TProfile("TP_d00000_c62233t0"     ,"TP_d00000_c62233t0"     ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d00000_c62233t0     ->Sumw2();
  //---------------------------------------------------------------------------------------------------------------------------
  TP_d000000_c111111t0   =new TProfile("TP_d000000_c111111t0"   ,"TP_d000000_c111111t0"   ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000000_c111111t0   ->Sumw2();
  TP_d000000_c222222t0   =new TProfile("TP_d000000_c222222t0"   ,"TP_d000000_c222222t0"   ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000000_c222222t0   ->Sumw2();
  TP_d000000_c333333t0   =new TProfile("TP_d000000_c333333t0"   ,"TP_d000000_c333333t0"   ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000000_c333333t0   ->Sumw2();
  //---------------------------------------------------------------------------------------------------------------------------
  TP_d000000_c112112t0   =new TProfile("TP_d000000_c112112t0"   ,"TP_d000000_c112112t0"   ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000000_c112112t0   ->Sumw2();
  TP_d000000_c113113t0   =new TProfile("TP_d000000_c113113t0"   ,"TP_d000000_c113113t0"   ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000000_c113113t0   ->Sumw2();
  TP_d000000_c221221t0   =new TProfile("TP_d000000_c221221t0"   ,"TP_d000000_c221221t0"   ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000000_c221221t0   ->Sumw2();
  TP_d000000_c331331t0   =new TProfile("TP_d000000_c331331t0"   ,"TP_d000000_c331331t0"   ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000000_c331331t0   ->Sumw2();
  TP_d000000_c123123t0   =new TProfile("TP_d000000_c123123t0"   ,"TP_d000000_c123123t0"   ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000000_c123123t0   ->Sumw2();
  //---------------------------------------------------------------------------------------------------------------------------
  TP_d000000_c223223t0   =new TProfile("TP_d000000_c223223t0"   ,"TP_d000000_c223223t0"   ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000000_c223223t0   ->Sumw2();
  TP_d000000_c224224t0   =new TProfile("TP_d000000_c224224t0"   ,"TP_d000000_c224224t0"   ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000000_c224224t0   ->Sumw2();
  TP_d000000_c225225t0   =new TProfile("TP_d000000_c225225t0"   ,"TP_d000000_c225225t0"   ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000000_c225225t0   ->Sumw2();
  TP_d000000_c332332t0   =new TProfile("TP_d000000_c332332t0"   ,"TP_d000000_c332332t0"   ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000000_c332332t0   ->Sumw2();
  TP_d000000_c334334t0   =new TProfile("TP_d000000_c334334t0"   ,"TP_d000000_c334334t0"   ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000000_c334334t0   ->Sumw2();
  TP_d000000_c442442t0   =new TProfile("TP_d000000_c442442t0"   ,"TP_d000000_c442442t0"   ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000000_c442442t0   ->Sumw2();
  TP_d000000_c234234t0   =new TProfile("TP_d000000_c234234t0"   ,"TP_d000000_c234234t0"   ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000000_c234234t0   ->Sumw2();
  TP_d000000_c235235t0   =new TProfile("TP_d000000_c235235t0"   ,"TP_d000000_c235235t0"   ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000000_c235235t0   ->Sumw2();
  TP_d000000_c345345t0   =new TProfile("TP_d000000_c345345t0"   ,"TP_d000000_c345345t0"   ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000000_c345345t0   ->Sumw2();
  //---------------------------------------------------------------------------------------------------------------------------
  TP_d00_c11sbc          =new TProfile("TP_d00_c11sbc"          ,"TP_d00_c11sbc"          ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d00_c11sbc          ->Sumw2();
  TP_d00_c22sbc          =new TProfile("TP_d00_c22sbc"          ,"TP_d00_c22sbc"          ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d00_c22sbc          ->Sumw2();
  TP_d00_c33sbc          =new TProfile("TP_d00_c33sbc"          ,"TP_d00_c33sbc"          ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d00_c33sbc          ->Sumw2();
  TP_d00_c44sbc          =new TProfile("TP_d00_c44sbc"          ,"TP_d00_c44sbc"          ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d00_c44sbc          ->Sumw2();
  TP_d00_c55sbc          =new TProfile("TP_d00_c55sbc"          ,"TP_d00_c55sbc"          ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d00_c55sbc          ->Sumw2();
  TP_d00_c66sbc          =new TProfile("TP_d00_c66sbc"          ,"TP_d00_c66sbc"          ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d00_c66sbc          ->Sumw2();
  //---------------------------------------------------------------------------------------------------------------------------
  TP_d000_c211sbc        =new TProfile("TP_d000_c211sbc"        ,"TP_d000_c211sbc"        ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000_c211sbc        ->Sumw2();
  TP_d000_c312sbc        =new TProfile("TP_d000_c312sbc"        ,"TP_d000_c312sbc"        ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000_c312sbc        ->Sumw2();
  TP_d000_c413sbc        =new TProfile("TP_d000_c413sbc"        ,"TP_d000_c413sbc"        ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000_c413sbc        ->Sumw2();
  TP_d000_c514sbc        =new TProfile("TP_d000_c514sbc"        ,"TP_d000_c514sbc"        ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000_c514sbc        ->Sumw2();
  TP_d000_c422sbc        =new TProfile("TP_d000_c422sbc"        ,"TP_d000_c422sbc"        ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000_c422sbc        ->Sumw2();
  TP_d000_c523sbc        =new TProfile("TP_d000_c523sbc"        ,"TP_d000_c523sbc"        ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000_c523sbc        ->Sumw2();
  TP_d000_c624sbc        =new TProfile("TP_d000_c624sbc"        ,"TP_d000_c624sbc"        ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000_c624sbc        ->Sumw2();
  TP_d000_c633sbc        =new TProfile("TP_d000_c633sbc"        ,"TP_d000_c633sbc"        ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d000_c633sbc        ->Sumw2();
  //---------------------------------------------------------------------------------------------------------------------------
  TP_d0000_c1111sbc      =new TProfile("TP_d0000_c1111sbc"      ,"TP_d0000_c1111sbc"      ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c1111sbc      ->Sumw2();
  TP_d0000_c2222sbc      =new TProfile("TP_d0000_c2222sbc"      ,"TP_d0000_c2222sbc"      ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c2222sbc      ->Sumw2();
  TP_d0000_c3333sbc      =new TProfile("TP_d0000_c3333sbc"      ,"TP_d0000_c3333sbc"      ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c3333sbc      ->Sumw2();
  TP_d0000_c1212sbc      =new TProfile("TP_d0000_c1212sbc"      ,"TP_d0000_c1212sbc"      ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c1212sbc      ->Sumw2();
  TP_d0000_c1313sbc      =new TProfile("TP_d0000_c1313sbc"      ,"TP_d0000_c1313sbc"      ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c1313sbc      ->Sumw2();
  TP_d0000_c1414sbc      =new TProfile("TP_d0000_c1414sbc"      ,"TP_d0000_c1414sbc"      ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c1414sbc      ->Sumw2();
  TP_d0000_c1515sbc      =new TProfile("TP_d0000_c1515sbc"      ,"TP_d0000_c1515sbc"      ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c1515sbc      ->Sumw2();
  TP_d0000_c2323sbc      =new TProfile("TP_d0000_c2323sbc"      ,"TP_d0000_c2323sbc"      ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c2323sbc      ->Sumw2();
  TP_d0000_c2424sbc      =new TProfile("TP_d0000_c2424sbc"      ,"TP_d0000_c2424sbc"      ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c2424sbc      ->Sumw2();
  TP_d0000_c2525sbc      =new TProfile("TP_d0000_c2525sbc"      ,"TP_d0000_c2525sbc"      ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c2525sbc      ->Sumw2();
  TP_d0000_c3434sbc      =new TProfile("TP_d0000_c3434sbc"      ,"TP_d0000_c3434sbc"      ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c3434sbc      ->Sumw2();
  TP_d0000_c3535sbc      =new TProfile("TP_d0000_c3535sbc"      ,"TP_d0000_c3535sbc"      ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c3535sbc      ->Sumw2();
  TP_d0000_c3111sbc      =new TProfile("TP_d0000_c3111sbc"      ,"TP_d0000_c3111sbc"      ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c3111sbc      ->Sumw2();
  TP_d0000_c6222sbc      =new TProfile("TP_d0000_c6222sbc"      ,"TP_d0000_c6222sbc"      ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c6222sbc      ->Sumw2();
  TP_d0000_c3324sbc      =new TProfile("TP_d0000_c3324sbc"      ,"TP_d0000_c3324sbc"      ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c3324sbc      ->Sumw2();
  TP_d0000_c4435sbc      =new TProfile("TP_d0000_c4435sbc"      ,"TP_d0000_c4435sbc"      ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c4435sbc      ->Sumw2();
  TP_d0000_c2534sbc      =new TProfile("TP_d0000_c2534sbc"      ,"TP_d0000_c2534sbc"      ,kNCentBins,-0.5,double(kNCentBins)-0.5);
  TP_d0000_c2534sbc      ->Sumw2();
  //----------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------
  for(int i=0; i<kNCentBins; i++){
    sprintf(histname1,"TP_d00_c00t0_00_%d", i )            ;
    TP_d00_c00t0_00[i]            =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d00_c00t0_00[i]            ->Sumw2();
      
    sprintf(histname1,"TP_d10_c11t0_v1_%d", i )            ;
    TP_d10_c11t0_v1[i]            =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d10_c11t0_v1[i]            ->Sumw2();
      
    sprintf(histname1,"TP_d10_c22t0_v2_%d", i )            ;
    TP_d10_c22t0_v2[i]            =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d10_c22t0_v2[i]            ->Sumw2();
      
    sprintf(histname1,"TP_d10_c33t0_v3_%d", i )            ;
    TP_d10_c33t0_v3[i]            =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d10_c33t0_v3[i]            ->Sumw2();
      
    sprintf(histname1,"TP_d10_c44t0_v4_%d", i )            ;
    TP_d10_c44t0_v4[i]            =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d10_c44t0_v4[i]            ->Sumw2();
      
    sprintf(histname1,"TP_d10_c55t0_v5_%d", i )            ;
    TP_d10_c55t0_v5[i]            =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d10_c55t0_v5[i]            ->Sumw2();
      
    sprintf(histname1,"TP_d10_c66t0_v6_%d", i )            ;
    TP_d10_c66t0_v6[i]            =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d10_c66t0_v6[i]            ->Sumw2();
      
    sprintf(histname1,"TP_d100_c211t0_v2_%d", i )          ;
    TP_d100_c211t0_v2[i]          =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d100_c211t0_v2[i]          ->Sumw2();
      
    sprintf(histname1,"TP_d100_c312t0_v3_%d", i )          ;
    TP_d100_c312t0_v3[i]          =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d100_c312t0_v3[i]          ->Sumw2();

    sprintf(histname1,"TP_d100_c413t0_v4_%d", i )          ;
    TP_d100_c413t0_v4[i]          =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d100_c413t0_v4[i]          ->Sumw2();
      
    sprintf(histname1,"TP_d100_c514t0_v5_%d", i )          ;
    TP_d100_c514t0_v5[i]          =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d100_c514t0_v5[i]          ->Sumw2();
      
    sprintf(histname1,"TP_d100_c422t0_v4_%d", i )          ;
    TP_d100_c422t0_v4[i]          =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d100_c422t0_v4[i]          ->Sumw2();
      
    sprintf(histname1,"TP_d100_c523t0_v5_%d", i )          ;
    TP_d100_c523t0_v5[i]          =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d100_c523t0_v5[i]          ->Sumw2();
      
    sprintf(histname1,"TP_d100_c624t0_v6_%d", i )          ;
    TP_d100_c624t0_v6[i]          =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d100_c624t0_v6[i]          ->Sumw2();
      
    sprintf(histname1,"TP_d100_c633t0_v6_%d", i )          ;
    TP_d100_c633t0_v6[i]          =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d100_c633t0_v6[i]          ->Sumw2();
      
    sprintf(histname1,"TP_d1000_c1111t0_v1_%d", i )        ;
    TP_d1000_c1111t0_v1[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c1111t0_v1[i]        ->Sumw2();
      
    sprintf(histname1,"TP_d1000_c2222t0_v2_%d", i )        ;
    TP_d1000_c2222t0_v2[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c2222t0_v2[i]        ->Sumw2();
      
    sprintf(histname1,"TP_d1000_c3333t0_v3_%d", i )        ;
    TP_d1000_c3333t0_v3[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c3333t0_v3[i]        ->Sumw2();
      
    sprintf(histname1,"TP_d1000_c1212t0_v1_%d", i )        ;
    TP_d1000_c1212t0_v1[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c1212t0_v1[i]        ->Sumw2();
      
    sprintf(histname1,"TP_d1000_c1313t0_v1_%d", i )        ;
    TP_d1000_c1313t0_v1[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c1313t0_v1[i]        ->Sumw2();
      
    sprintf(histname1,"TP_d1000_c1414t0_v1_%d", i )        ;
    TP_d1000_c1414t0_v1[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c1414t0_v1[i]        ->Sumw2();
      
    sprintf(histname1,"TP_d1000_c1515t0_v1_%d", i )        ;
    TP_d1000_c1515t0_v1[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c1515t0_v1[i]        ->Sumw2();
      
    sprintf(histname1,"TP_d1000_c2323t0_v2_%d", i )        ;
    TP_d1000_c2323t0_v2[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c2323t0_v2[i]        ->Sumw2();
      
    sprintf(histname1,"TP_d1000_c2424t0_v2_%d", i )        ;
    TP_d1000_c2424t0_v2[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c2424t0_v2[i]        ->Sumw2();
      
    sprintf(histname1,"TP_d1000_c2525t0_v2_%d", i )        ;
    TP_d1000_c2525t0_v2[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c2525t0_v2[i]        ->Sumw2();
      
    sprintf(histname1,"TP_d1000_c3434t0_v3_%d", i )        ;
    TP_d1000_c3434t0_v3[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c3434t0_v3[i]        ->Sumw2();
      
    sprintf(histname1,"TP_d1000_c3535t0_v3_%d", i )        ;
    TP_d1000_c3535t0_v3[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c3535t0_v3[i]        ->Sumw2();
      
    sprintf(histname1,"TP_d1000_c3111t0_v3_%d", i )        ;
    TP_d1000_c3111t0_v3[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c3111t0_v3[i]        ->Sumw2();
      
    sprintf(histname1,"TP_d1000_c6222t0_v6_%d", i )        ;
    TP_d1000_c6222t0_v6[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c6222t0_v6[i]        ->Sumw2();
      
    sprintf(histname1,"TP_d1000_c3324t0_v3_%d", i )        ;
    TP_d1000_c3324t0_v3[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c3324t0_v3[i]        ->Sumw2();
      
    sprintf(histname1,"TP_d1000_c4435t0_v4_%d", i )        ;
    TP_d1000_c4435t0_v4[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c4435t0_v4[i]        ->Sumw2();
      
    sprintf(histname1,"TP_d1000_c2534t0_v2_%d", i )        ;
    TP_d1000_c2534t0_v2[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c2534t0_v2[i]        ->Sumw2();
      
    sprintf(histname1,"TP_d0200_c1212t0_v2_%d", i )        ;
    TP_d0200_c1212t0_v2[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d0200_c1212t0_v2[i]        ->Sumw2();
      
    sprintf(histname1,"TP_d0200_c1313t0_v3_%d", i )        ;
    TP_d0200_c1313t0_v3[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d0200_c1313t0_v3[i]        ->Sumw2();
      
    sprintf(histname1,"TP_d0200_c1414t0_v4_%d", i )        ;
    TP_d0200_c1414t0_v4[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d0200_c1414t0_v4[i]        ->Sumw2();
      
    sprintf(histname1,"TP_d0200_c1515t0_v5_%d", i )        ;
    TP_d0200_c1515t0_v5[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d0200_c1515t0_v5[i]        ->Sumw2();
      
    sprintf(histname1,"TP_d0200_c2323t0_v3_%d", i )        ;
    TP_d0200_c2323t0_v3[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d0200_c2323t0_v3[i]        ->Sumw2();
      
    sprintf(histname1,"TP_d0200_c2424t0_v4_%d", i )        ;
    TP_d0200_c2424t0_v4[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d0200_c2424t0_v4[i]        ->Sumw2();
      
    sprintf(histname1,"TP_d0200_c2525t0_v5_%d", i )        ;
    TP_d0200_c2525t0_v5[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d0200_c2525t0_v5[i]        ->Sumw2();
      
    sprintf(histname1,"TP_d0200_c3434t0_v4_%d", i )        ;
    TP_d0200_c3434t0_v4[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d0200_c3434t0_v4[i]        ->Sumw2();

    sprintf(histname1,"TP_d0200_c3535t0_v5_%d", i )        ;
    TP_d0200_c3535t0_v5[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d0200_c3535t0_v5[i]        ->Sumw2();
      
    sprintf(histname1,"TP_d0030_c3324t0_v2_%d", i )        ;
    TP_d0030_c3324t0_v2[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d0030_c3324t0_v2[i]        ->Sumw2();
      
    sprintf(histname1,"TP_d0030_c4435t0_v3_%d", i )        ;
    TP_d0030_c4435t0_v3[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d0030_c4435t0_v3[i]        ->Sumw2();
      
    sprintf(histname1,"TP_d0004_c4435t0_v5_%d", i )        ;
    TP_d0004_c4435t0_v5[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d0004_c4435t0_v5[i]        ->Sumw2();
      
    sprintf(histname1,"TP_d0200_c2534t0_v5_%d", i )        ;
    TP_d0200_c2534t0_v5[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d0200_c2534t0_v5[i]        ->Sumw2();
      
    sprintf(histname1,"TP_d0030_c2534t0_v3_%d", i )        ;
    TP_d0030_c2534t0_v3[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d0030_c2534t0_v3[i]        ->Sumw2();
      
    sprintf(histname1,"TP_d0004_c2534t0_v4_%d", i )        ;
    TP_d0004_c2534t0_v4[i]        =new TProfile(histname1,histname1,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d0004_c2534t0_v4[i]        ->Sumw2();


    sprintf(histname2,"TP_d10000_c21111t0_v2_%d", i )      ;
    TP_d10000_c21111t0_v2[i]      =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d10000_c21111t0_v2[i]      ->Sumw2();
      
    sprintf(histname2,"TP_d10000_c42222t0_v4_%d", i )      ;
    TP_d10000_c42222t0_v4[i]      =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d10000_c42222t0_v4[i]      ->Sumw2();
      
    sprintf(histname2,"TP_d10000_c43322t0_v4_%d", i )      ;
    TP_d10000_c43322t0_v4[i]      =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d10000_c43322t0_v4[i]      ->Sumw2();
      
    sprintf(histname2,"TP_d10000_c52232t0_v5_%d", i )      ;
    TP_d10000_c52232t0_v5[i]      =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d10000_c52232t0_v5[i]      ->Sumw2();
      
    sprintf(histname2,"TP_d10000_c53332t0_v5_%d", i )      ;
    TP_d10000_c53332t0_v5[i]      =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d10000_c53332t0_v5[i]      ->Sumw2();
      
    sprintf(histname2,"TP_d10000_c22435t0_v2_%d", i )      ;
    TP_d10000_c22435t0_v2[i]      =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d10000_c22435t0_v2[i]      ->Sumw2();
      
    sprintf(histname2,"TP_d10000_c23344t0_v2_%d", i )      ;
    TP_d10000_c23344t0_v2[i]      =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d10000_c23344t0_v2[i]      ->Sumw2();
      
    sprintf(histname2,"TP_d10000_c33345t0_v3_%d", i )      ;
    TP_d10000_c33345t0_v3[i]      =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d10000_c33345t0_v3[i]      ->Sumw2();
      
    sprintf(histname2,"TP_d10000_c24455t0_v2_%d", i )      ;
    TP_d10000_c24455t0_v2[i]      =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d10000_c24455t0_v2[i]      ->Sumw2();
      
    sprintf(histname2,"TP_d10000_c33455t0_v3_%d", i )      ;
    TP_d10000_c33455t0_v3[i]      =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d10000_c33455t0_v3[i]      ->Sumw2();
      
    sprintf(histname2,"TP_d10000_c22233t0_v2_%d", i )      ;
    TP_d10000_c22233t0_v2[i]      =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d10000_c22233t0_v2[i]      ->Sumw2();
      
    sprintf(histname2,"TP_d10000_c62233t0_v6_%d", i )      ;
    TP_d10000_c62233t0_v6[i]      =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d10000_c62233t0_v6[i]      ->Sumw2();
      
    sprintf(histname2,"TP_d100000_c111111t0_v1_%d", i )    ;
    TP_d100000_c111111t0_v1[i]    =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d100000_c111111t0_v1[i]    ->Sumw2();
      
    sprintf(histname2,"TP_d100000_c222222t0_v2_%d", i )    ;
    TP_d100000_c222222t0_v2[i]    =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d100000_c222222t0_v2[i]    ->Sumw2();
      
    sprintf(histname2,"TP_d100000_c333333t0_v3_%d", i )    ;
    TP_d100000_c333333t0_v3[i]    =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d100000_c333333t0_v3[i]    ->Sumw2();
      
    sprintf(histname2,"TP_d100000_c223223t0_v2_%d", i )    ;
    TP_d100000_c223223t0_v2[i]    =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d100000_c223223t0_v2[i]    ->Sumw2();
      
    sprintf(histname2,"TP_d100000_c224224t0_v2_%d", i )    ;
    TP_d100000_c224224t0_v2[i]    =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d100000_c224224t0_v2[i]    ->Sumw2();
      
    sprintf(histname2,"TP_d100000_c225225t0_v2_%d", i )    ;
    TP_d100000_c225225t0_v2[i]    =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d100000_c225225t0_v2[i]    ->Sumw2();

    sprintf(histname2,"TP_d100000_c332332t0_v3_%d", i )    ;
    TP_d100000_c332332t0_v3[i]    =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d100000_c332332t0_v3[i]    ->Sumw2();
      
    sprintf(histname2,"TP_d100000_c334334t0_v3_%d", i )    ;
    TP_d100000_c334334t0_v3[i]    =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d100000_c334334t0_v3[i]    ->Sumw2();
      
    sprintf(histname2,"TP_d100000_c442442t0_v4_%d", i )    ;
    TP_d100000_c442442t0_v4[i]    =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d100000_c442442t0_v4[i]    ->Sumw2();
      
    sprintf(histname2,"TP_d100000_c234234t0_v2_%d", i )    ;
    TP_d100000_c234234t0_v2[i]    =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d100000_c234234t0_v2[i]    ->Sumw2();
      
    sprintf(histname2,"TP_d100000_c235235t0_v2_%d", i )    ;
    TP_d100000_c235235t0_v2[i]    =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d100000_c235235t0_v2[i]    ->Sumw2();
      
    sprintf(histname2,"TP_d100000_c345345t0_v3_%d", i )    ;
    TP_d100000_c345345t0_v3[i]    =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d100000_c345345t0_v3[i]    ->Sumw2();
      
    sprintf(histname2,"TP_d020000_c223223t0_v3_%d", i )    ;
    TP_d020000_c223223t0_v3[i]    =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d020000_c223223t0_v3[i]    ->Sumw2();
      
    sprintf(histname2,"TP_d020000_c224224t0_v4_%d", i )    ;
    TP_d020000_c224224t0_v4[i]    =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d020000_c224224t0_v4[i]    ->Sumw2();
      
    sprintf(histname2,"TP_d020000_c225225t0_v5_%d", i )    ;
    TP_d020000_c225225t0_v5[i]    =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d020000_c225225t0_v5[i]    ->Sumw2();
      
    sprintf(histname2,"TP_d020000_c332332t0_v2_%d", i )    ;
    TP_d020000_c332332t0_v2[i]    =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d020000_c332332t0_v2[i]    ->Sumw2();
      
    sprintf(histname2,"TP_d020000_c334334t0_v4_%d", i )    ;
    TP_d020000_c334334t0_v4[i]    =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d020000_c334334t0_v4[i]    ->Sumw2();
      
    sprintf(histname2,"TP_d020000_c442442t0_v2_%d", i )    ;
    TP_d020000_c442442t0_v2[i]    =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d020000_c442442t0_v2[i]    ->Sumw2();
      
    sprintf(histname2,"TP_d020000_c234234t0_v3_%d", i )    ;
    TP_d020000_c234234t0_v3[i]    =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d020000_c234234t0_v3[i]    ->Sumw2();
      
    sprintf(histname2,"TP_d003000_c234234t0_v4_%d", i )    ;
    TP_d003000_c234234t0_v4[i]    =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d003000_c234234t0_v4[i]    ->Sumw2();
      
    sprintf(histname2,"TP_d020000_c235235t0_v3_%d", i )    ;
    TP_d020000_c235235t0_v3[i]    =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d020000_c235235t0_v3[i]    ->Sumw2();
      
    sprintf(histname2,"TP_d003000_c235235t0_v5_%d", i )    ;
    TP_d003000_c235235t0_v5[i]    =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d003000_c235235t0_v5[i]    ->Sumw2();
      
    sprintf(histname2,"TP_d020000_c345345t0_v4_%d", i )    ;
    TP_d020000_c345345t0_v4[i]    =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d020000_c345345t0_v4[i]    ->Sumw2();
      
    sprintf(histname2,"TP_d003000_c345345t0_v5_%d", i )    ;
    TP_d003000_c345345t0_v5[i]    =new TProfile(histname2,histname2,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d003000_c345345t0_v5[i]    ->Sumw2();
    //------------------------------------------------------------------------------------------  
    //------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------
    sprintf(histname5,"TP_d10_c11sbc_v1_%d", i )           ;
    TP_d10_c11sbc_v1[i]           =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d10_c11sbc_v1[i]           ->Sumw2();

    sprintf(histname5,"TP_d10_c22sbc_v2_%d", i )           ;
    TP_d10_c22sbc_v2[i]           =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d10_c22sbc_v2[i]           ->Sumw2();

    sprintf(histname5,"TP_d10_c33sbc_v3_%d", i )           ;
    TP_d10_c33sbc_v3[i]           =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d10_c33sbc_v3[i]           ->Sumw2();

    sprintf(histname5,"TP_d10_c44sbc_v4_%d", i )           ;
    TP_d10_c44sbc_v4[i]           =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d10_c44sbc_v4[i]           ->Sumw2();

    sprintf(histname5,"TP_d10_c55sbc_v5_%d", i )           ;
    TP_d10_c55sbc_v5[i]           =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d10_c55sbc_v5[i]           ->Sumw2();

    sprintf(histname5,"TP_d10_c66sbc_v6_%d", i )           ;
    TP_d10_c66sbc_v6[i]           =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d10_c66sbc_v6[i]           ->Sumw2();

    sprintf(histname5,"TP_d100_c211sbc_v2_%d", i )         ;
    TP_d100_c211sbc_v2[i]         =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d100_c211sbc_v2[i]         ->Sumw2();

    sprintf(histname5,"TP_d100_c312sbc_v3_%d", i )         ;
    TP_d100_c312sbc_v3[i]         =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d100_c312sbc_v3[i]         ->Sumw2();

    sprintf(histname5,"TP_d100_c413sbc_v4_%d", i )         ;
    TP_d100_c413sbc_v4[i]         =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d100_c413sbc_v4[i]         ->Sumw2();

    sprintf(histname5,"TP_d100_c514sbc_v5_%d", i )         ;
    TP_d100_c514sbc_v5[i]         =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d100_c514sbc_v5[i]         ->Sumw2();

    sprintf(histname5,"TP_d100_c422sbc_v4_%d", i )         ;
    TP_d100_c422sbc_v4[i]         =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d100_c422sbc_v4[i]         ->Sumw2();

    sprintf(histname5,"TP_d100_c523sbc_v5_%d", i )         ;
    TP_d100_c523sbc_v5[i]         =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d100_c523sbc_v5[i]         ->Sumw2();

    sprintf(histname5,"TP_d100_c624sbc_v6_%d", i )         ;
    TP_d100_c624sbc_v6[i]         =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d100_c624sbc_v6[i]         ->Sumw2();

    sprintf(histname5,"TP_d100_c633sbc_v6_%d", i )         ;
    TP_d100_c633sbc_v6[i]         =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d100_c633sbc_v6[i]         ->Sumw2();

    sprintf(histname5,"TP_d1000_c1111sbc_v1_%d", i )       ;
    TP_d1000_c1111sbc_v1[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c1111sbc_v1[i]       ->Sumw2();

    sprintf(histname5,"TP_d1000_c2222sbc_v2_%d", i )       ;
    TP_d1000_c2222sbc_v2[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c2222sbc_v2[i]       ->Sumw2();

    sprintf(histname5,"TP_d1000_c3333sbc_v3_%d", i )       ;
    TP_d1000_c3333sbc_v3[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c3333sbc_v3[i]       ->Sumw2();

    sprintf(histname5,"TP_d1000_c1212sbc_v1_%d", i )       ;
    TP_d1000_c1212sbc_v1[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c1212sbc_v1[i]       ->Sumw2();

    sprintf(histname5,"TP_d1000_c1313sbc_v1_%d", i )       ;
    TP_d1000_c1313sbc_v1[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c1313sbc_v1[i]       ->Sumw2();

    sprintf(histname5,"TP_d1000_c1414sbc_v1_%d", i )       ;
    TP_d1000_c1414sbc_v1[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c1414sbc_v1[i]       ->Sumw2();

    sprintf(histname5,"TP_d1000_c1515sbc_v1_%d", i )       ;
    TP_d1000_c1515sbc_v1[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c1515sbc_v1[i]       ->Sumw2();

    sprintf(histname5,"TP_d1000_c2323sbc_v2_%d", i )       ;
    TP_d1000_c2323sbc_v2[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c2323sbc_v2[i]       ->Sumw2();

    sprintf(histname5,"TP_d1000_c2424sbc_v2_%d", i )       ;
    TP_d1000_c2424sbc_v2[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c2424sbc_v2[i]       ->Sumw2();

    sprintf(histname5,"TP_d1000_c2525sbc_v2_%d", i )       ;
    TP_d1000_c2525sbc_v2[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c2525sbc_v2[i]       ->Sumw2();

    sprintf(histname5,"TP_d1000_c3434sbc_v3_%d", i )       ;
    TP_d1000_c3434sbc_v3[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c3434sbc_v3[i]       ->Sumw2();

    sprintf(histname5,"TP_d1000_c3535sbc_v3_%d", i )       ;
    TP_d1000_c3535sbc_v3[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c3535sbc_v3[i]       ->Sumw2();

    sprintf(histname5,"TP_d1000_c3111sbc_v3_%d", i )       ;
    TP_d1000_c3111sbc_v3[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c3111sbc_v3[i]       ->Sumw2();

    sprintf(histname5,"TP_d1000_c6222sbc_v6_%d", i )       ;
    TP_d1000_c6222sbc_v6[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c6222sbc_v6[i]       ->Sumw2();

    sprintf(histname5,"TP_d1000_c3324sbc_v3_%d", i )       ;
    TP_d1000_c3324sbc_v3[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c3324sbc_v3[i]       ->Sumw2();

    sprintf(histname5,"TP_d1000_c4435sbc_v4_%d", i )       ;
    TP_d1000_c4435sbc_v4[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c4435sbc_v4[i]       ->Sumw2();

    sprintf(histname5,"TP_d1000_c2534sbc_v2_%d", i )       ;
    TP_d1000_c2534sbc_v2[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d1000_c2534sbc_v2[i]       ->Sumw2();

    sprintf(histname5,"TP_d0200_c1212sbc_v2_%d", i )       ;
    TP_d0200_c1212sbc_v2[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d0200_c1212sbc_v2[i]       ->Sumw2();

    sprintf(histname5,"TP_d0200_c1313sbc_v3_%d", i )       ;
    TP_d0200_c1313sbc_v3[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d0200_c1313sbc_v3[i]       ->Sumw2();

    sprintf(histname5,"TP_d0200_c1414sbc_v4_%d", i )       ;
    TP_d0200_c1414sbc_v4[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d0200_c1414sbc_v4[i]       ->Sumw2();

    sprintf(histname5,"TP_d0200_c1515sbc_v5_%d", i )       ;
    TP_d0200_c1515sbc_v5[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d0200_c1515sbc_v5[i]       ->Sumw2();

    sprintf(histname5,"TP_d0200_c2323sbc_v3_%d", i )       ;
    TP_d0200_c2323sbc_v3[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d0200_c2323sbc_v3[i]       ->Sumw2();

    sprintf(histname5,"TP_d0200_c2424sbc_v4_%d", i )       ;
    TP_d0200_c2424sbc_v4[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d0200_c2424sbc_v4[i]       ->Sumw2();

    sprintf(histname5,"TP_d0200_c2525sbc_v5_%d", i )       ;
    TP_d0200_c2525sbc_v5[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d0200_c2525sbc_v5[i]       ->Sumw2();

    sprintf(histname5,"TP_d0200_c3434sbc_v4_%d", i )       ;
    TP_d0200_c3434sbc_v4[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d0200_c3434sbc_v4[i]       ->Sumw2();

    sprintf(histname5,"TP_d0200_c3535sbc_v5_%d", i )       ;
    TP_d0200_c3535sbc_v5[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d0200_c3535sbc_v5[i]       ->Sumw2();

    sprintf(histname5,"TP_d0030_c3324sbc_v2_%d", i )       ;
    TP_d0030_c3324sbc_v2[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d0030_c3324sbc_v2[i]       ->Sumw2();

    sprintf(histname5,"TP_d0030_c4435sbc_v3_%d", i )       ;
    TP_d0030_c4435sbc_v3[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d0030_c4435sbc_v3[i]       ->Sumw2();

    sprintf(histname5,"TP_d0004_c4435sbc_v5_%d", i )       ;
    TP_d0004_c4435sbc_v5[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d0004_c4435sbc_v5[i]       ->Sumw2();

    sprintf(histname5,"TP_d0200_c2534sbc_v5_%d", i )       ;
    TP_d0200_c2534sbc_v5[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d0200_c2534sbc_v5[i]       ->Sumw2();

    sprintf(histname5,"TP_d0030_c2534sbc_v3_%d", i )       ;
    TP_d0030_c2534sbc_v3[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d0030_c2534sbc_v3[i]       ->Sumw2();

    sprintf(histname5,"TP_d0004_c2534sbc_v4_%d", i )       ;
    TP_d0004_c2534sbc_v4[i]       =new TProfile(histname5,histname5,kNPtaBins,-0.5,double(kNPtaBins)-0.5);
    TP_d0004_c2534sbc_v4[i]       ->Sumw2();
  }
  //----------------------------------------------------------------------------------------
  //-------------------------------------------------------------
}
//---------------------------------------------------------------
//---------------------------------------------------------------
//---------------------------------------------------------------
int main(int argc, char* argv[]){
  //-------------------------------------------------------------
  Params P;
  std::string err;

  if (argc >= 12) {
    // Legacy positional mode: argv[1]..argv[11]
    P.MODE    = std::stoi(argv[1]);
    P.DECAY   = std::stoi(argv[2]);
    P.FRAMEID = std::stoi(argv[3]);
    P.N_EVENT = std::stoi(argv[4]);
    P.IAP     = std::stoi(argv[5]);
    P.IZP     = std::stoi(argv[6]);
    P.IAT     = std::stoi(argv[7]);
    P.IZT     = std::stoi(argv[8]);
    P.EFRM    = std::stof(argv[9]);
    P.BMIN    = std::stof(argv[10]);
    P.BMAX    = std::stof(argv[11]);
  } else {
    // Config mode: argv[1] is a path, or "-" for stdin. Also support "-c file".
    std::string path;
    if (argc == 3 && (std::string(argv[1]) == "-c" || std::string(argv[1]) == "--config")) {
      path = argv[2];
    } else if (argc == 2) {
      path = argv[1];
    } else {
      // No args → read from stdin so `./run < params.txt` works
      path = "-";
    }

    if (!load_params_from_file(path, P, err)) {
      std::cerr << "Config error: " << err << "\n\n"
                << "Usage (positional): " << argv[0]
                << " MODE DECAY FRAMEID N_EVENT IAP IZP IAT IZT EFRM BMIN BMAX\n"
                << "or config file:     " << argv[0] << " -c params.txt\n"
                << "or from stdin:      " << argv[0] << " < params.txt\n";
      return 1;
    }
  }

  // --- use P here ---
  // e.g. Wrapper wrapper; wrapper.OpenOutput("output_histograms.root"); wrapper.BookHistos();
  // ... run with P.MODE, P.DECAY, etc. ...

  std::cout << "Loaded params: MODE=" << P.MODE << " DECAY=" << P.DECAY
            << " FRAMEID=" << P.FRAMEID << " N_EVENT=" << P.N_EVENT
            << " IAP=" << P.IAP << " IZP=" << P.IZP
            << " IAT=" << P.IAT << " IZT=" << P.IZT
            << " EFRM=" << P.EFRM << " BMIN=" << P.BMIN << " BMAX=" << P.BMAX
            << "\n";
  //-------------------------------------------------------------
  //-------------------------------------------------------------
  Wrapper wrapper;
  wrapper.OpenOutput("output_histograms.root");  // file->cd(); AddDirectory(true)
  wrapper.BookHistos();                          // creates & attaches all histos
  //-------------------------------------------------------------
  //-------------------------------------------------------------
  char FRAME[2][9];
  std::strcpy(FRAME[0], "CMS     ");
  std::strcpy(FRAME[1], "LAB     ");

  char PROJ[2][9];
  std::strcpy(PROJ[0], "P       ");
  std::strcpy(PROJ[1], "A       ");

  char TARG[2][9];
  std::strcpy(TARG[0],"P       ");
  std::strcpy(TARG[1],"A       ");
  
  int PROJID;
  if(P.IAP == 1){
    PROJID = 0;
  }else{
    PROJID = 1;
  }

  int TARGID;
  if(P.IAT == 1){
    TARGID = 0;
  }else{
    TARGID = 1;
  }

  float RIonP    = 1.2*pow(P.IAP,0.3333333333);
  float RIonT    = 1.2*pow(P.IAT,0.3333333333);
  float RPT      =  (RIonP+RIonT)/2.0;
  //----------------------------------------------------------------------------------------
  // Get the current time as a duration since the epoch
  auto currentTime = std::chrono::system_clock::now().time_since_epoch();
  auto currentTimeMs = std::chrono::duration_cast<std::chrono::milliseconds>(currentTime).count();
  // Seed the random number generator with the current time
  std::mt19937_64 generator(currentTimeMs);
  // Generate a random odd number
  std::uniform_int_distribution<int> distribution(1, 20000);
  int randomNumber = distribution(generator);
  if (randomNumber % 2 == 0) {
    randomNumber++; // Make sure it's odd
  }  
  std::cout << "Random odd number: " << randomNumber << std::endl;
  //----------------------------------------------------------------------------------------
  ranseed.nseed    = randomNumber;
  //hiparnt.ihpr2[9] =0;
  //hiparnt.ihpr2[11] =0; //For pi0 K0S D+- Lambda Cascad+-
  hiparnt.ihpr2[20]=P.DECAY;
  //----------------------------------------------------------------------------------------
  hijset_(&P.EFRM, FRAME[P.FRAMEID], PROJ[PROJID], TARG[TARGID], &P.IAP, &P.IZP, &P.IAT, &P.IZT);
  //--------------------------------------------------------------------------------------
  int Nev_tot = 0;
  for (int IE = 0; IE <= P.N_EVENT; ++IE) {
    //------------------------------------------------------------------------------------
    hijing_(FRAME[P.FRAMEID], &P.BMIN, &P.BMAX);
    //------------------------------------------------------------------------------------
    if(himain1.natt <= 1) continue;
    //------------------------------------------------------------------------------------
    if(hiparnt.hint1[18] > 1.3*(RIonP + RIonT) ) continue;
    //------------------------------------------------------------------------------------
    float Imb = hiparnt.hint1[18];
    //------------------------------------------------------------------------------------
    int   npart0     = 0;
    float vxpart[600]={0.};
    float vypart[600]={0.};
    //------------------------------------------------------------------------------------
    for(int ip=0; ip<P.IAP; ip++){
      if(histrng.nfp[4][ip] > 1){
	vxpart[npart0] = hijcrdn.yp[ip][0];
	vypart[npart0] = hijcrdn.yp[ip][1];
	npart0++;
      }
    }
    for(int it=0; it<P.IAT; it++){
      if(histrng.nft[4][it] > 1){
	vxpart[npart0] = hijcrdn.yt[it][0];
	vypart[npart0] = hijcrdn.yt[it][1];
	npart0++;
      }      
    }
    //------------------------------------------------------------------------------------
    int Npart  = npart0;
    //------------------------------------------------------------------------------------
    int CentralityID =  getNparCentID(Npart);
    //------------------------------------------------------------------------------------
    //-------------------Loop over TPC tracks---------------------------------------------
    //------------------------------------------------------------------------------------
    int Ntrak = himain1.natt;
    
    double Deta = 0.6;
    wrapper_qcorr::WrapperQCorrEvent qcorrEvent(Deta, kNPtaBins - 1);
    qcorrEvent.reserve(static_cast<std::size_t>(Ntrak));
    double Npt[kNPtaBins] = {0};
    double Wpt[kNPtaBins] = {0};
    int Nch = 0;
    //-----------------------------------------------------------------------------------
    for(int part=0; part<himain1.natt; part++){
      //---------------------------------------------------------------------------------
      int    ID = himain2.katt[0][part];
      double px = himain2.patt[0][part];
      double py = himain2.patt[1][part];
      double pz = himain2.patt[2][part];
      double en = himain2.patt[3][part];
      //----------------------------------------------------------------------------------
      // Hijing gives V0's status=4, they need to have status=1 to be decayed in geant
      // also change status=11 to status=2
      int status = himain2.katt[3][part];
      if(status <=10 && status>0 ) status = 1;
      if(status <=20 && status>10) status = 2;
      //----------------------------------------------------------------------------------
      if(status != 1) continue;
      //----------------------------------------------------------------------------------
      int pid = -999;
      if(fabs(ID) == 211)       pid = 0;
      else if(fabs(ID) == 321 ) pid = 1;
      else if(fabs(ID) == 2212) pid = 2;
      else continue;
      if(pid < 0) continue;
      //----------------------------------------------------------------------------------
      int ch = ID/fabs(ID);
      //----------------------------------------------------------------------------------
      TLorentzVector p4;
      p4.SetPxPyPzE(px, py, pz, en);          // (px,py,pz,E)
      //----------------------------------------------------------------------------------
      // Directly from TLorentzVector:
      double pT  = p4.Pt();      // transverse momentum
      double Eta = p4.Eta();     // pseudorapidity
      double Phi = p4.Phi();     // azimuth in (-pi, pi]
      //Phi = TVector2::Phi_0_2pi(Phi);  // now φ ∈ [0,2π)
      //double r = dis(gen);
      //----------------------------------------------------------------------------------
      if(fabs(Eta) < 0.5 && pT > 0.2 && pT < 3.0) Nch ++;
      //----------------------------------------------------------------------------------
      if (pT < ptMin || pT > ptMax) continue;
      if (std::fabs(Eta) > 1.0) continue;
      //----------------------------------------------------------------------------------
      float weight = 1.0;
      //----------------------------------------------------------------------------------
      wrapper.Etah[CentralityID]->Fill(Eta,weight);
      wrapper.Phih[CentralityID]->Fill(Phi,weight);
      wrapper.PTh[CentralityID]->Fill(pT,weight);
      //----------------------------------------------------------------------------------
      int etaBin = getEtaBin(Eta);
      if (etaBin < 0 || etaBin >= kNEtaBins) continue;

      int ptBin = getPtaBin(pT);
      if (ptBin < 0 || ptBin >= kNPtaBins) continue;
      //----------------------------------------------------------------------------------
      const bool isRefPt = (pT <= RefPTcut);
      qcorrEvent.add_track(Phi, Eta, pT, ptBin, weight, 1.0, isRefPt, true, true);
      Npt[ptBin] += pT;
      Wpt[ptBin] += 1.0;
      //----------------------------------------------------------------------------------
    }
    //------------------------------------------------------------------------------------
    if( Nch < 2) continue;
    //------------------------------------------------------------------------------------
    wrapper.ImbaA->Fill(Imb);
    wrapper.NparA->Fill(Npart);
    wrapper.MultA->Fill(Nch);
    //------------------------------------------------------------------------------------
    wrapper.NChvsCent->Fill(CentralityID, Nch);
    wrapper.NpartvsCent->Fill(CentralityID, Npart);
    //------------------------------------------------------------------------------------
    double Fixval = 0.00000000000001;
    //------------------------------------------------------------------------------------
    using QCorrSub = wrapper_qcorr::WrapperQCorrEvent::Subevent;
    const auto QS_A = QCorrSub::A; // eta > +gap/2
    const auto QS_B = QCorrSub::B; // eta < -gap/2
    //------------------------------------------------------------------------------------
    const auto r_d00_w00t0 = qcorrEvent.integrated_same({0,  0}, {false, false});
    const auto r_d00_c11t0 = qcorrEvent.integrated_same({1, -1}, {FMpTW, FMpTW});
    const auto r_d00_c22t0 = qcorrEvent.integrated_same({2, -2}, {false, false});
    const auto r_d00_c33t0 = qcorrEvent.integrated_same({3, -3}, {false, false});
    const auto r_d00_c44t0 = qcorrEvent.integrated_same({4, -4}, {false, false});
    const auto r_d00_c55t0 = qcorrEvent.integrated_same({5, -5}, {false, false});
    const auto r_d00_c66t0 = qcorrEvent.integrated_same({6, -6}, {false, false});
    const double d00_w00t0 = r_d00_w00t0.denominator;
    const double d00_c11t0 = r_d00_c11t0.numerator.real();
    const double d00_c22t0 = r_d00_c22t0.numerator.real();
    const double d00_c33t0 = r_d00_c33t0.numerator.real();
    const double d00_c44t0 = r_d00_c44t0.numerator.real();
    const double d00_c55t0 = r_d00_c55t0.numerator.real();
    const double d00_c66t0 = r_d00_c66t0.numerator.real();
    if(fabs(d00_w00t0) > Fixval){
      wrapper.TP_d00_c11t0           ->Fill(CentralityID, (d00_c11t0           )/(d00_w00t0           ), (d00_w00t0           ));
      wrapper.TP_d00_c22t0           ->Fill(CentralityID, (d00_c22t0           )/(d00_w00t0           ), (d00_w00t0           ));
      wrapper.TP_d00_c33t0           ->Fill(CentralityID, (d00_c33t0           )/(d00_w00t0           ), (d00_w00t0           ));
      wrapper.TP_d00_c44t0           ->Fill(CentralityID, (d00_c44t0           )/(d00_w00t0           ), (d00_w00t0           ));
      wrapper.TP_d00_c55t0           ->Fill(CentralityID, (d00_c55t0           )/(d00_w00t0           ), (d00_w00t0           ));
      wrapper.TP_d00_c66t0           ->Fill(CentralityID, (d00_c66t0           )/(d00_w00t0           ), (d00_w00t0           ));
    }
    //------------------------------------------------------------------------------------
    const auto r_d000_w000t0 = qcorrEvent.integrated_same({0,  0,  0}, {false, false, false});
    const auto r_d000_c211t0 = qcorrEvent.integrated_same({2, -1, -1}, {false, FMpTW, FMpTW});
    const auto r_d000_c312t0 = qcorrEvent.integrated_same({3, -1, -2}, {false, FMpTW, false});
    const auto r_d000_c413t0 = qcorrEvent.integrated_same({4, -1, -3}, {false, FMpTW, false});
    const auto r_d000_c514t0 = qcorrEvent.integrated_same({5, -1, -4}, {false, FMpTW, false});
    const auto r_d000_c422t0 = qcorrEvent.integrated_same({4, -2, -2}, {false, false, false});
    const auto r_d000_c523t0 = qcorrEvent.integrated_same({5, -2, -3}, {false, false, false});
    const auto r_d000_c624t0 = qcorrEvent.integrated_same({6, -2, -4}, {false, false, false});
    const auto r_d000_c633t0 = qcorrEvent.integrated_same({6, -3, -3}, {false, false, false});
    const double d000_w000t0 = r_d000_w000t0.denominator;
    const double d000_c211t0 = r_d000_c211t0.numerator.real();
    const double d000_c312t0 = r_d000_c312t0.numerator.real();
    const double d000_c413t0 = r_d000_c413t0.numerator.real();
    const double d000_c514t0 = r_d000_c514t0.numerator.real();
    const double d000_c422t0 = r_d000_c422t0.numerator.real();
    const double d000_c523t0 = r_d000_c523t0.numerator.real();
    const double d000_c624t0 = r_d000_c624t0.numerator.real();
    const double d000_c633t0 = r_d000_c633t0.numerator.real();
    if(fabs(d000_w000t0) > Fixval){
      wrapper.TP_d000_c211t0         ->Fill(CentralityID, (d000_c211t0         )/(d000_w000t0         ), (d000_w000t0         ));
      wrapper.TP_d000_c312t0         ->Fill(CentralityID, (d000_c312t0         )/(d000_w000t0         ), (d000_w000t0         ));
      wrapper.TP_d000_c413t0         ->Fill(CentralityID, (d000_c413t0         )/(d000_w000t0         ), (d000_w000t0         ));
      wrapper.TP_d000_c514t0         ->Fill(CentralityID, (d000_c514t0         )/(d000_w000t0         ), (d000_w000t0         ));
      wrapper.TP_d000_c422t0         ->Fill(CentralityID, (d000_c422t0         )/(d000_w000t0         ), (d000_w000t0         ));
      wrapper.TP_d000_c523t0         ->Fill(CentralityID, (d000_c523t0         )/(d000_w000t0         ), (d000_w000t0         ));
      wrapper.TP_d000_c624t0         ->Fill(CentralityID, (d000_c624t0         )/(d000_w000t0         ), (d000_w000t0         ));
      wrapper.TP_d000_c633t0         ->Fill(CentralityID, (d000_c633t0         )/(d000_w000t0         ), (d000_w000t0         ));
    }
    //------------------------------------------------------------------------------------
    const auto r_d0000_w0000t0 = qcorrEvent.integrated_same({0, 0,  0,  0}, {false, false, false, false});
    const auto r_d0000_c1111t0 = qcorrEvent.integrated_same({1, 1, -1, -1}, {FMpTW, FMpTW, FMpTW, FMpTW});
    const auto r_d0000_c2222t0 = qcorrEvent.integrated_same({2, 2, -2, -2}, {false, false, false, false});
    const auto r_d0000_c3333t0 = qcorrEvent.integrated_same({3, 3, -3, -3}, {false, false, false, false});
    const auto r_d0000_c1212t0 = qcorrEvent.integrated_same({1, 2, -1, -2}, {FMpTW, false, FMpTW, false});
    const auto r_d0000_c1313t0 = qcorrEvent.integrated_same({1, 3, -1, -3}, {FMpTW, false, FMpTW, false});
    const auto r_d0000_c1414t0 = qcorrEvent.integrated_same({1, 4, -1, -4}, {FMpTW, false, FMpTW, false});
    const auto r_d0000_c1515t0 = qcorrEvent.integrated_same({1, 5, -1, -5}, {FMpTW, false, FMpTW, false});
    const auto r_d0000_c2323t0 = qcorrEvent.integrated_same({2, 3, -2, -3}, {false, false, false, false});
    const auto r_d0000_c2424t0 = qcorrEvent.integrated_same({2, 4, -2, -4}, {false, false, false, false});
    const auto r_d0000_c2525t0 = qcorrEvent.integrated_same({2, 5, -2, -5}, {false, false, false, false});
    const auto r_d0000_c3434t0 = qcorrEvent.integrated_same({3, 4, -3, -4}, {false, false, false, false});
    const auto r_d0000_c3535t0 = qcorrEvent.integrated_same({3, 5, -3, -5}, {false, false, false, false});
    const auto r_d0000_c3111t0 = qcorrEvent.integrated_same({3,-1, -1, -1}, {false, FMpTW, FMpTW, FMpTW});
    const auto r_d0000_c6222t0 = qcorrEvent.integrated_same({6,-2, -2, -2}, {false, false, false, false});
    const auto r_d0000_c3324t0 = qcorrEvent.integrated_same({3, 3, -2, -4}, {false, false, false, false});
    const auto r_d0000_c4435t0 = qcorrEvent.integrated_same({4, 4, -3, -5}, {false, false, false, false});
    const auto r_d0000_c2534t0 = qcorrEvent.integrated_same({2, 5, -3, -4}, {false, false, false, false});
    const double d0000_w0000t0 = r_d0000_w0000t0.denominator;
    const double d0000_c1111t0 = r_d0000_c1111t0.numerator.real();
    const double d0000_c2222t0 = r_d0000_c2222t0.numerator.real();
    const double d0000_c3333t0 = r_d0000_c3333t0.numerator.real();
    const double d0000_c1212t0 = r_d0000_c1212t0.numerator.real();
    const double d0000_c1313t0 = r_d0000_c1313t0.numerator.real();
    const double d0000_c1414t0 = r_d0000_c1414t0.numerator.real();
    const double d0000_c1515t0 = r_d0000_c1515t0.numerator.real();
    const double d0000_c2323t0 = r_d0000_c2323t0.numerator.real();
    const double d0000_c2424t0 = r_d0000_c2424t0.numerator.real();
    const double d0000_c2525t0 = r_d0000_c2525t0.numerator.real();
    const double d0000_c3434t0 = r_d0000_c3434t0.numerator.real();
    const double d0000_c3535t0 = r_d0000_c3535t0.numerator.real();
    const double d0000_c3111t0 = r_d0000_c3111t0.numerator.real();
    const double d0000_c6222t0 = r_d0000_c6222t0.numerator.real();
    const double d0000_c3324t0 = r_d0000_c3324t0.numerator.real();
    const double d0000_c4435t0 = r_d0000_c4435t0.numerator.real();
    const double d0000_c2534t0 = r_d0000_c2534t0.numerator.real();
    if(fabs(d0000_w0000t0) > Fixval){
      wrapper.TP_d0000_c1111t0       ->Fill(CentralityID, (d0000_c1111t0       )/(d0000_w0000t0       ), (d0000_w0000t0       ));
      wrapper.TP_d0000_c2222t0       ->Fill(CentralityID, (d0000_c2222t0       )/(d0000_w0000t0       ), (d0000_w0000t0       ));
      wrapper.TP_d0000_c3333t0       ->Fill(CentralityID, (d0000_c3333t0       )/(d0000_w0000t0       ), (d0000_w0000t0       ));
      wrapper.TP_d0000_c1212t0       ->Fill(CentralityID, (d0000_c1212t0       )/(d0000_w0000t0       ), (d0000_w0000t0       ));
      wrapper.TP_d0000_c1313t0       ->Fill(CentralityID, (d0000_c1313t0       )/(d0000_w0000t0       ), (d0000_w0000t0       ));
      wrapper.TP_d0000_c1414t0       ->Fill(CentralityID, (d0000_c1414t0       )/(d0000_w0000t0       ), (d0000_w0000t0       ));
      wrapper.TP_d0000_c1515t0       ->Fill(CentralityID, (d0000_c1515t0       )/(d0000_w0000t0       ), (d0000_w0000t0       ));
      wrapper.TP_d0000_c2323t0       ->Fill(CentralityID, (d0000_c2323t0       )/(d0000_w0000t0       ), (d0000_w0000t0       ));
      wrapper.TP_d0000_c2424t0       ->Fill(CentralityID, (d0000_c2424t0       )/(d0000_w0000t0       ), (d0000_w0000t0       ));
      wrapper.TP_d0000_c2525t0       ->Fill(CentralityID, (d0000_c2525t0       )/(d0000_w0000t0       ), (d0000_w0000t0       ));
      wrapper.TP_d0000_c3434t0       ->Fill(CentralityID, (d0000_c3434t0       )/(d0000_w0000t0       ), (d0000_w0000t0       ));
      wrapper.TP_d0000_c3535t0       ->Fill(CentralityID, (d0000_c3535t0       )/(d0000_w0000t0       ), (d0000_w0000t0       ));
      wrapper.TP_d0000_c3111t0       ->Fill(CentralityID, (d0000_c3111t0       )/(d0000_w0000t0       ), (d0000_w0000t0       ));
      wrapper.TP_d0000_c6222t0       ->Fill(CentralityID, (d0000_c6222t0       )/(d0000_w0000t0       ), (d0000_w0000t0       ));
      wrapper.TP_d0000_c3324t0       ->Fill(CentralityID, (d0000_c3324t0       )/(d0000_w0000t0       ), (d0000_w0000t0       ));
      wrapper.TP_d0000_c4435t0       ->Fill(CentralityID, (d0000_c4435t0       )/(d0000_w0000t0       ), (d0000_w0000t0       ));
      wrapper.TP_d0000_c2534t0       ->Fill(CentralityID, (d0000_c2534t0       )/(d0000_w0000t0       ), (d0000_w0000t0       ));
    }
    //------------------------------------------------------------------------------------
    const auto r_d00000_w00000t0 = qcorrEvent.integrated_same({0, 0,  0,  0,  0}, {false, false, false, false, false});
    const auto r_d00000_c21111t0 = qcorrEvent.integrated_same({2, 1, -1, -1, -1}, {false, FMpTW, FMpTW, FMpTW, FMpTW});
    const auto r_d00000_c42222t0 = qcorrEvent.integrated_same({4, 2, -2, -2, -2}, {false, false, false, false, false});
    const auto r_d00000_c43322t0 = qcorrEvent.integrated_same({4, 3, -3, -2, -2}, {false, false, false, false, false});
    const auto r_d00000_c52232t0 = qcorrEvent.integrated_same({5, 2, -2, -3, -2}, {false, false, false, false, false});
    const auto r_d00000_c53332t0 = qcorrEvent.integrated_same({5, 3, -3, -3, -2}, {false, false, false, false, false});
    const auto r_d00000_c22435t0 = qcorrEvent.integrated_same({2, 2,  4, -3, -5}, {false, false, false, false, false});
    const auto r_d00000_c23344t0 = qcorrEvent.integrated_same({2, 3,  3, -4, -4}, {false, false, false, false, false});
    const auto r_d00000_c33345t0 = qcorrEvent.integrated_same({3, 3,  3, -4, -5}, {false, false, false, false, false});
    const auto r_d00000_c24455t0 = qcorrEvent.integrated_same({2, 4,  4, -5, -5}, {false, false, false, false, false});
    const auto r_d00000_c33455t0 = qcorrEvent.integrated_same({3, 3,  4, -5, -5}, {false, false, false, false, false});
    const auto r_d00000_c22233t0 = qcorrEvent.integrated_same({2, 2,  2, -3, -3}, {false, false, false, false, false});
    const auto r_d00000_c62233t0 = qcorrEvent.integrated_same({6, 2, -2, -3, -3}, {false, false, false, false, false});
    const double d00000_w00000t0 = r_d00000_w00000t0.denominator;
    const double d00000_c21111t0 = r_d00000_c21111t0.numerator.real();
    const double d00000_c42222t0 = r_d00000_c42222t0.numerator.real();
    const double d00000_c43322t0 = r_d00000_c43322t0.numerator.real();
    const double d00000_c52232t0 = r_d00000_c52232t0.numerator.real();
    const double d00000_c53332t0 = r_d00000_c53332t0.numerator.real();
    const double d00000_c22435t0 = r_d00000_c22435t0.numerator.real();
    const double d00000_c23344t0 = r_d00000_c23344t0.numerator.real();
    const double d00000_c33345t0 = r_d00000_c33345t0.numerator.real();
    const double d00000_c24455t0 = r_d00000_c24455t0.numerator.real();
    const double d00000_c33455t0 = r_d00000_c33455t0.numerator.real();
    const double d00000_c22233t0 = r_d00000_c22233t0.numerator.real();
    const double d00000_c62233t0 = r_d00000_c62233t0.numerator.real();
    if(fabs(d00000_w00000t0) > Fixval){
      wrapper.TP_d00000_c21111t0     ->Fill(CentralityID, (d00000_c21111t0     )/(d00000_w00000t0     ), (d00000_w00000t0     ));
      wrapper.TP_d00000_c42222t0     ->Fill(CentralityID, (d00000_c42222t0     )/(d00000_w00000t0     ), (d00000_w00000t0     ));
      wrapper.TP_d00000_c43322t0     ->Fill(CentralityID, (d00000_c43322t0     )/(d00000_w00000t0     ), (d00000_w00000t0     ));
      wrapper.TP_d00000_c52232t0     ->Fill(CentralityID, (d00000_c52232t0     )/(d00000_w00000t0     ), (d00000_w00000t0     ));
      wrapper.TP_d00000_c53332t0     ->Fill(CentralityID, (d00000_c53332t0     )/(d00000_w00000t0     ), (d00000_w00000t0     ));
      wrapper.TP_d00000_c22435t0     ->Fill(CentralityID, (d00000_c22435t0     )/(d00000_w00000t0     ), (d00000_w00000t0     ));
      wrapper.TP_d00000_c23344t0     ->Fill(CentralityID, (d00000_c23344t0     )/(d00000_w00000t0     ), (d00000_w00000t0     ));
      wrapper.TP_d00000_c33345t0     ->Fill(CentralityID, (d00000_c33345t0     )/(d00000_w00000t0     ), (d00000_w00000t0     ));
      wrapper.TP_d00000_c24455t0     ->Fill(CentralityID, (d00000_c24455t0     )/(d00000_w00000t0     ), (d00000_w00000t0     ));
      wrapper.TP_d00000_c33455t0     ->Fill(CentralityID, (d00000_c33455t0     )/(d00000_w00000t0     ), (d00000_w00000t0     ));
      wrapper.TP_d00000_c22233t0     ->Fill(CentralityID, (d00000_c22233t0     )/(d00000_w00000t0     ), (d00000_w00000t0     ));
      wrapper.TP_d00000_c62233t0     ->Fill(CentralityID, (d00000_c62233t0     )/(d00000_w00000t0     ), (d00000_w00000t0     ));
    }
    //------------------------------------------------------------------------------------
    const auto r_d000000_w000000t0 = qcorrEvent.integrated_same({0,  0, 0,  0,  0,  0}, {false, false, false, false, false, false});
    const auto r_d000000_c111111t0 = qcorrEvent.integrated_same({1,  1, 1, -1, -1, -1}, {FMpTW, FMpTW, FMpTW, FMpTW, FMpTW, FMpTW});
    const auto r_d000000_c222222t0 = qcorrEvent.integrated_same({2,  2, 2, -2, -2, -2}, {false, false, false, false, false, false});
    const auto r_d000000_c333333t0 = qcorrEvent.integrated_same({3,  3, 3, -3, -3, -3}, {false, false, false, false, false, false});


    const auto r_d000000_c112112t0 = qcorrEvent.integrated_same({1,  1, 2, -1, -1, -2}, {FMpTW, FMpTW, false, FMpTW, FMpTW, false});
    const auto r_d000000_c113113t0 = qcorrEvent.integrated_same({1,  1, 3, -1, -1, -3}, {FMpTW, FMpTW, false, FMpTW, FMpTW, false});
    const auto r_d000000_c221221t0 = qcorrEvent.integrated_same({2,  2, 1, -2, -2, -1}, {false, false, FMpTW, false, false, FMpTW});
    const auto r_d000000_c331331t0 = qcorrEvent.integrated_same({3,  3, 1, -3, -3, -1}, {false, false, FMpTW, false, false, FMpTW});
    const auto r_d000000_c123123t0 = qcorrEvent.integrated_same({1,  2, 3, -1, -2, -3}, {FMpTW, false, false, FMpTW, false, false});
    

    const auto r_d000000_c223223t0 = qcorrEvent.integrated_same({2,  2, 3, -2, -2, -3}, {false, false, false, false, false, false});
    const auto r_d000000_c224224t0 = qcorrEvent.integrated_same({2,  2, 4, -2, -2, -4}, {false, false, false, false, false, false});
    const auto r_d000000_c225225t0 = qcorrEvent.integrated_same({2,  2, 5, -2, -2, -5}, {false, false, false, false, false, false});
    const auto r_d000000_c332332t0 = qcorrEvent.integrated_same({3,  3, 2, -3, -3, -2}, {false, false, false, false, false, false});
    const auto r_d000000_c334334t0 = qcorrEvent.integrated_same({3,  3, 4, -3, -3, -4}, {false, false, false, false, false, false});
    const auto r_d000000_c442442t0 = qcorrEvent.integrated_same({4,  4, 2, -4, -4, -2}, {false, false, false, false, false, false});
    const auto r_d000000_c234234t0 = qcorrEvent.integrated_same({2, -2, 3, -3,  4, -4}, {false, false, false, false, false, false});
    const auto r_d000000_c235235t0 = qcorrEvent.integrated_same({2, -2, 3, -3,  5, -5}, {false, false, false, false, false, false});
    const auto r_d000000_c345345t0 = qcorrEvent.integrated_same({3, -3, 4, -4,  5, -5}, {false, false, false, false, false, false});

    
    const double d000000_w000000t0 = r_d000000_w000000t0.denominator;
    const double d000000_c111111t0 = r_d000000_c111111t0.numerator.real();
    const double d000000_c222222t0 = r_d000000_c222222t0.numerator.real();
    const double d000000_c333333t0 = r_d000000_c333333t0.numerator.real();

    const double d000000_c112112t0 = r_d000000_c112112t0.numerator.real();
    const double d000000_c113113t0 = r_d000000_c113113t0.numerator.real();
    const double d000000_c221221t0 = r_d000000_c221221t0.numerator.real();
    const double d000000_c331331t0 = r_d000000_c331331t0.numerator.real();
    const double d000000_c123123t0 = r_d000000_c123123t0.numerator.real();

    const double d000000_c223223t0 = r_d000000_c223223t0.numerator.real();
    const double d000000_c224224t0 = r_d000000_c224224t0.numerator.real();
    const double d000000_c225225t0 = r_d000000_c225225t0.numerator.real();
    const double d000000_c332332t0 = r_d000000_c332332t0.numerator.real();
    const double d000000_c334334t0 = r_d000000_c334334t0.numerator.real();
    const double d000000_c442442t0 = r_d000000_c442442t0.numerator.real();
    const double d000000_c234234t0 = r_d000000_c234234t0.numerator.real();
    const double d000000_c235235t0 = r_d000000_c235235t0.numerator.real();
    const double d000000_c345345t0 = r_d000000_c345345t0.numerator.real();

    if(fabs(d000000_w000000t0) > Fixval){
      wrapper.TP_d000000_c111111t0   ->Fill(CentralityID, (d000000_c111111t0   )/(d000000_w000000t0   ), (d000000_w000000t0   ));
      wrapper.TP_d000000_c222222t0   ->Fill(CentralityID, (d000000_c222222t0   )/(d000000_w000000t0   ), (d000000_w000000t0   ));
      wrapper.TP_d000000_c333333t0   ->Fill(CentralityID, (d000000_c333333t0   )/(d000000_w000000t0   ), (d000000_w000000t0   ));
      wrapper.TP_d000000_c112112t0   ->Fill(CentralityID, (d000000_c112112t0   )/(d000000_w000000t0   ), (d000000_w000000t0   ));
      wrapper.TP_d000000_c113113t0   ->Fill(CentralityID, (d000000_c113113t0   )/(d000000_w000000t0   ), (d000000_w000000t0   ));
      wrapper.TP_d000000_c221221t0   ->Fill(CentralityID, (d000000_c221221t0   )/(d000000_w000000t0   ), (d000000_w000000t0   ));
      wrapper.TP_d000000_c331331t0   ->Fill(CentralityID, (d000000_c331331t0   )/(d000000_w000000t0   ), (d000000_w000000t0   ));
      wrapper.TP_d000000_c223223t0   ->Fill(CentralityID, (d000000_c223223t0   )/(d000000_w000000t0   ), (d000000_w000000t0   ));
      wrapper.TP_d000000_c224224t0   ->Fill(CentralityID, (d000000_c224224t0   )/(d000000_w000000t0   ), (d000000_w000000t0   ));
      wrapper.TP_d000000_c225225t0   ->Fill(CentralityID, (d000000_c225225t0   )/(d000000_w000000t0   ), (d000000_w000000t0   ));
      wrapper.TP_d000000_c332332t0   ->Fill(CentralityID, (d000000_c332332t0   )/(d000000_w000000t0   ), (d000000_w000000t0   ));
      wrapper.TP_d000000_c334334t0   ->Fill(CentralityID, (d000000_c334334t0   )/(d000000_w000000t0   ), (d000000_w000000t0   ));
      wrapper.TP_d000000_c442442t0   ->Fill(CentralityID, (d000000_c442442t0   )/(d000000_w000000t0   ), (d000000_w000000t0   ));
      wrapper.TP_d000000_c234234t0   ->Fill(CentralityID, (d000000_c234234t0   )/(d000000_w000000t0   ), (d000000_w000000t0   ));
      wrapper.TP_d000000_c235235t0   ->Fill(CentralityID, (d000000_c235235t0   )/(d000000_w000000t0   ), (d000000_w000000t0   ));
      wrapper.TP_d000000_c345345t0   ->Fill(CentralityID, (d000000_c345345t0   )/(d000000_w000000t0   ), (d000000_w000000t0   ));
    }
    //------------------------------------------------------------------------------------
    //-------------------------------------Two-Sub----------------------------------------
    //------------------------------------------------------------------------------------
    const auto r_d00_w00sbc = qcorrEvent.integrated_two_subevent({0,  0}, {QS_A, QS_B}, {false, false});
    const auto r_d00_c11sbc = qcorrEvent.integrated_two_subevent({1, -1}, {QS_A, QS_B}, {FMpTW, FMpTW});
    const auto r_d00_c22sbc = qcorrEvent.integrated_two_subevent({2, -2}, {QS_A, QS_B}, {false, false});
    const auto r_d00_c33sbc = qcorrEvent.integrated_two_subevent({3, -3}, {QS_A, QS_B}, {false, false});
    const auto r_d00_c44sbc = qcorrEvent.integrated_two_subevent({4, -4}, {QS_A, QS_B}, {false, false});
    const auto r_d00_c55sbc = qcorrEvent.integrated_two_subevent({5, -5}, {QS_A, QS_B}, {false, false});
    const auto r_d00_c66sbc = qcorrEvent.integrated_two_subevent({6, -6}, {QS_A, QS_B}, {false, false});
    const double d00_w00sbc = r_d00_w00sbc.denominator;
    const double d00_c11sbc = r_d00_c11sbc.numerator.real();
    const double d00_c22sbc = r_d00_c22sbc.numerator.real();
    const double d00_c33sbc = r_d00_c33sbc.numerator.real();
    const double d00_c44sbc = r_d00_c44sbc.numerator.real();
    const double d00_c55sbc = r_d00_c55sbc.numerator.real();
    const double d00_c66sbc = r_d00_c66sbc.numerator.real();
    if(fabs(d00_w00sbc) > Fixval){
      wrapper.TP_d00_c11sbc          ->Fill(CentralityID, (d00_c11sbc          )/(d00_w00sbc          ), (d00_w00sbc          ));
      wrapper.TP_d00_c22sbc          ->Fill(CentralityID, (d00_c22sbc          )/(d00_w00sbc          ), (d00_w00sbc          ));
      wrapper.TP_d00_c33sbc          ->Fill(CentralityID, (d00_c33sbc          )/(d00_w00sbc          ), (d00_w00sbc          ));
      wrapper.TP_d00_c44sbc          ->Fill(CentralityID, (d00_c44sbc          )/(d00_w00sbc          ), (d00_w00sbc          ));
      wrapper.TP_d00_c55sbc          ->Fill(CentralityID, (d00_c55sbc          )/(d00_w00sbc          ), (d00_w00sbc          ));
      wrapper.TP_d00_c66sbc          ->Fill(CentralityID, (d00_c66sbc          )/(d00_w00sbc          ), (d00_w00sbc          ));
    }
    //------------------------------------------------------------------------------------
    const auto r_d000_w000sbc = qcorrEvent.integrated_two_subevent({0,  0,  0}, {QS_A, QS_A, QS_B}, {false, false, false});
    const auto r_d000_c211sbc = qcorrEvent.integrated_two_subevent({2, -1, -1}, {QS_A, QS_A, QS_B}, {false, FMpTW, FMpTW});
    const auto r_d000_c312sbc = qcorrEvent.integrated_two_subevent({3, -1, -2}, {QS_A, QS_A, QS_B}, {false, FMpTW, false});
    const auto r_d000_c413sbc = qcorrEvent.integrated_two_subevent({4, -1, -3}, {QS_A, QS_A, QS_B}, {false, FMpTW, false});
    const auto r_d000_c514sbc = qcorrEvent.integrated_two_subevent({5, -1, -4}, {QS_A, QS_A, QS_B}, {false, FMpTW, false});
    const auto r_d000_c422sbc = qcorrEvent.integrated_two_subevent({4, -2, -2}, {QS_A, QS_A, QS_B}, {false, false, false});
    const auto r_d000_c523sbc = qcorrEvent.integrated_two_subevent({5, -2, -3}, {QS_A, QS_A, QS_B}, {false, false, false});
    const auto r_d000_c624sbc = qcorrEvent.integrated_two_subevent({6, -2, -4}, {QS_A, QS_A, QS_B}, {false, false, false});
    const auto r_d000_c633sbc = qcorrEvent.integrated_two_subevent({6, -3, -3}, {QS_A, QS_A, QS_B}, {false, false, false});
    const double d000_w000sbc = r_d000_w000sbc.denominator;
    const double d000_c211sbc = r_d000_c211sbc.numerator.real();
    const double d000_c312sbc = r_d000_c312sbc.numerator.real();
    const double d000_c413sbc = r_d000_c413sbc.numerator.real();
    const double d000_c514sbc = r_d000_c514sbc.numerator.real();
    const double d000_c422sbc = r_d000_c422sbc.numerator.real();
    const double d000_c523sbc = r_d000_c523sbc.numerator.real();
    const double d000_c624sbc = r_d000_c624sbc.numerator.real();
    const double d000_c633sbc = r_d000_c633sbc.numerator.real();
    if(fabs(d000_w000sbc) > Fixval){
      wrapper.TP_d000_c211sbc        ->Fill(CentralityID, (d000_c211sbc        )/(d000_w000sbc        ), (d000_w000sbc        ));
      wrapper.TP_d000_c312sbc        ->Fill(CentralityID, (d000_c312sbc        )/(d000_w000sbc        ), (d000_w000sbc        ));
      wrapper.TP_d000_c413sbc        ->Fill(CentralityID, (d000_c413sbc        )/(d000_w000sbc        ), (d000_w000sbc        ));
      wrapper.TP_d000_c514sbc        ->Fill(CentralityID, (d000_c514sbc        )/(d000_w000sbc        ), (d000_w000sbc        ));
      wrapper.TP_d000_c422sbc        ->Fill(CentralityID, (d000_c422sbc        )/(d000_w000sbc        ), (d000_w000sbc        ));
      wrapper.TP_d000_c523sbc        ->Fill(CentralityID, (d000_c523sbc        )/(d000_w000sbc        ), (d000_w000sbc        ));
      wrapper.TP_d000_c624sbc        ->Fill(CentralityID, (d000_c624sbc        )/(d000_w000sbc        ), (d000_w000sbc        ));
      wrapper.TP_d000_c633sbc        ->Fill(CentralityID, (d000_c633sbc        )/(d000_w000sbc        ), (d000_w000sbc        ));
    }
    //------------------------------------------------------------------------------------
    const auto r_d0000_w0000sbc = qcorrEvent.integrated_two_subevent({0,  0,  0,  0}, {QS_A, QS_A, QS_B, QS_B}, {false, false, false, false});
    const auto r_d0000_c1111sbc = qcorrEvent.integrated_two_subevent({1,  1, -1, -1}, {QS_A, QS_A, QS_B, QS_B}, {FMpTW, FMpTW, FMpTW, FMpTW});
    const auto r_d0000_c2222sbc = qcorrEvent.integrated_two_subevent({2,  2, -2, -2}, {QS_A, QS_A, QS_B, QS_B}, {false, false, false, false});
    const auto r_d0000_c3333sbc = qcorrEvent.integrated_two_subevent({3,  3, -3, -3}, {QS_A, QS_A, QS_B, QS_B}, {false, false, false, false});
    const auto r_d0000_c1212sbc = qcorrEvent.integrated_two_subevent({1,  2, -1, -2}, {QS_A, QS_A, QS_B, QS_B}, {FMpTW, false, FMpTW, false});
    const auto r_d0000_c1313sbc = qcorrEvent.integrated_two_subevent({1,  3, -1, -3}, {QS_A, QS_A, QS_B, QS_B}, {FMpTW, false, FMpTW, false});
    const auto r_d0000_c1414sbc = qcorrEvent.integrated_two_subevent({1,  4, -1, -4}, {QS_A, QS_A, QS_B, QS_B}, {FMpTW, false, FMpTW, false});
    const auto r_d0000_c1515sbc = qcorrEvent.integrated_two_subevent({1,  5, -1, -5}, {QS_A, QS_A, QS_B, QS_B}, {FMpTW, false, FMpTW, false});
    const auto r_d0000_c2323sbc = qcorrEvent.integrated_two_subevent({2,  3, -2, -3}, {QS_A, QS_A, QS_B, QS_B}, {false, false, false, false});
    const auto r_d0000_c2424sbc = qcorrEvent.integrated_two_subevent({2,  4, -2, -4}, {QS_A, QS_A, QS_B, QS_B}, {false, false, false, false});
    const auto r_d0000_c2525sbc = qcorrEvent.integrated_two_subevent({2,  5, -2, -5}, {QS_A, QS_A, QS_B, QS_B}, {false, false, false, false});
    const auto r_d0000_c3434sbc = qcorrEvent.integrated_two_subevent({3,  4, -3, -4}, {QS_A, QS_A, QS_B, QS_B}, {false, false, false, false});
    const auto r_d0000_c3535sbc = qcorrEvent.integrated_two_subevent({3,  5, -3, -5}, {QS_A, QS_A, QS_B, QS_B}, {false, false, false, false});
    const auto r_d0000_c3111sbc = qcorrEvent.integrated_two_subevent({3, -1, -1, -1}, {QS_A, QS_A, QS_B, QS_B}, {false, FMpTW, FMpTW, FMpTW});
    const auto r_d0000_c6222sbc = qcorrEvent.integrated_two_subevent({6, -2, -2, -2}, {QS_A, QS_A, QS_B, QS_B}, {false, false, false, false});
    const auto r_d0000_c3324sbc = qcorrEvent.integrated_two_subevent({3,  3, -2, -4}, {QS_A, QS_A, QS_B, QS_B}, {false, false, false, false});
    const auto r_d0000_c4435sbc = qcorrEvent.integrated_two_subevent({4,  4, -3, -5}, {QS_A, QS_A, QS_B, QS_B}, {false, false, false, false});
    const auto r_d0000_c2534sbc = qcorrEvent.integrated_two_subevent({2,  5, -3, -4}, {QS_A, QS_A, QS_B, QS_B}, {false, false, false, false});
    const double d0000_w0000sbc = r_d0000_w0000sbc.denominator;
    const double d0000_c1111sbc = r_d0000_c1111sbc.numerator.real();
    const double d0000_c2222sbc = r_d0000_c2222sbc.numerator.real();
    const double d0000_c3333sbc = r_d0000_c3333sbc.numerator.real();
    const double d0000_c1212sbc = r_d0000_c1212sbc.numerator.real();
    const double d0000_c1313sbc = r_d0000_c1313sbc.numerator.real();
    const double d0000_c1414sbc = r_d0000_c1414sbc.numerator.real();
    const double d0000_c1515sbc = r_d0000_c1515sbc.numerator.real();
    const double d0000_c2323sbc = r_d0000_c2323sbc.numerator.real();
    const double d0000_c2424sbc = r_d0000_c2424sbc.numerator.real();
    const double d0000_c2525sbc = r_d0000_c2525sbc.numerator.real();
    const double d0000_c3434sbc = r_d0000_c3434sbc.numerator.real();
    const double d0000_c3535sbc = r_d0000_c3535sbc.numerator.real();
    const double d0000_c3111sbc = r_d0000_c3111sbc.numerator.real();
    const double d0000_c6222sbc = r_d0000_c6222sbc.numerator.real();
    const double d0000_c3324sbc = r_d0000_c3324sbc.numerator.real();
    const double d0000_c4435sbc = r_d0000_c4435sbc.numerator.real();
    const double d0000_c2534sbc = r_d0000_c2534sbc.numerator.real();
    if(fabs(d0000_w0000sbc) > Fixval){
      wrapper.TP_d0000_c1111sbc      ->Fill(CentralityID, (d0000_c1111sbc      )/(d0000_w0000sbc      ), (d0000_w0000sbc      ));
      wrapper.TP_d0000_c2222sbc      ->Fill(CentralityID, (d0000_c2222sbc      )/(d0000_w0000sbc      ), (d0000_w0000sbc      ));
      wrapper.TP_d0000_c3333sbc      ->Fill(CentralityID, (d0000_c3333sbc      )/(d0000_w0000sbc      ), (d0000_w0000sbc      ));
      wrapper.TP_d0000_c1212sbc      ->Fill(CentralityID, (d0000_c1212sbc      )/(d0000_w0000sbc      ), (d0000_w0000sbc      ));
      wrapper.TP_d0000_c1313sbc      ->Fill(CentralityID, (d0000_c1313sbc      )/(d0000_w0000sbc      ), (d0000_w0000sbc      ));
      wrapper.TP_d0000_c1414sbc      ->Fill(CentralityID, (d0000_c1414sbc      )/(d0000_w0000sbc      ), (d0000_w0000sbc      ));
      wrapper.TP_d0000_c1515sbc      ->Fill(CentralityID, (d0000_c1515sbc      )/(d0000_w0000sbc      ), (d0000_w0000sbc      ));
      wrapper.TP_d0000_c2323sbc      ->Fill(CentralityID, (d0000_c2323sbc      )/(d0000_w0000sbc      ), (d0000_w0000sbc      ));
      wrapper.TP_d0000_c2424sbc      ->Fill(CentralityID, (d0000_c2424sbc      )/(d0000_w0000sbc      ), (d0000_w0000sbc      ));
      wrapper.TP_d0000_c2525sbc      ->Fill(CentralityID, (d0000_c2525sbc      )/(d0000_w0000sbc      ), (d0000_w0000sbc      ));
      wrapper.TP_d0000_c3434sbc      ->Fill(CentralityID, (d0000_c3434sbc      )/(d0000_w0000sbc      ), (d0000_w0000sbc      ));
      wrapper.TP_d0000_c3535sbc      ->Fill(CentralityID, (d0000_c3535sbc      )/(d0000_w0000sbc      ), (d0000_w0000sbc      ));
      wrapper.TP_d0000_c3111sbc      ->Fill(CentralityID, (d0000_c3111sbc      )/(d0000_w0000sbc      ), (d0000_w0000sbc      ));
      wrapper.TP_d0000_c6222sbc      ->Fill(CentralityID, (d0000_c6222sbc      )/(d0000_w0000sbc      ), (d0000_w0000sbc      ));
      wrapper.TP_d0000_c3324sbc      ->Fill(CentralityID, (d0000_c3324sbc      )/(d0000_w0000sbc      ), (d0000_w0000sbc      ));
      wrapper.TP_d0000_c4435sbc      ->Fill(CentralityID, (d0000_c4435sbc      )/(d0000_w0000sbc      ), (d0000_w0000sbc      ));
      wrapper.TP_d0000_c2534sbc      ->Fill(CentralityID, (d0000_c2534sbc      )/(d0000_w0000sbc      ), (d0000_w0000sbc      ));
    }
    //------------------------------------------------------------------------------------
    //-----------------------------------Diffrential-PC-----------------------------------
    //------------------------------------------------------------------------------------
    for(int i=0; i< kNPtaBins-1; i++){
      
      if (Wpt[i] <= 0.0) {
	continue;
      }
      
      if(Wpt[i] > 0.0){
	wrapper.TP_d00_c00t0_00[CentralityID]->Fill(i, Npt[i]/Wpt[i]);
      }
      if(Wpt[i] < 0.0) continue;
      //----------------------------------------------------------------------------------
      const auto r_dd0_w00t0_w0 = qcorrEvent.differential_one_same({0,  0}, {0}, i, {false, false});
      const auto r_d10_c11t0_v1 = qcorrEvent.differential_one_same({1, -1}, {0}, i, {false, FMpTW});
      const auto r_d10_c22t0_v2 = qcorrEvent.differential_one_same({2, -2}, {0}, i, {false, false});
      const auto r_d10_c33t0_v3 = qcorrEvent.differential_one_same({3, -3}, {0}, i, {false, false});
      const auto r_d10_c44t0_v4 = qcorrEvent.differential_one_same({4, -4}, {0}, i, {false, false});
      const auto r_d10_c55t0_v5 = qcorrEvent.differential_one_same({5, -5}, {0}, i, {false, false});
      const auto r_d10_c66t0_v6 = qcorrEvent.differential_one_same({6, -6}, {0}, i, {false, false});
      const double dd0_w00t0_w0 = r_dd0_w00t0_w0.denominator;
      if(dd0_w00t0_w0 < Fixval) continue;
      const double d10_c11t0_v1 = r_d10_c11t0_v1.numerator.real();
      const double d10_c22t0_v2 = r_d10_c22t0_v2.numerator.real();
      const double d10_c33t0_v3 = r_d10_c33t0_v3.numerator.real();
      const double d10_c44t0_v4 = r_d10_c44t0_v4.numerator.real();
      const double d10_c55t0_v5 = r_d10_c55t0_v5.numerator.real();
      const double d10_c66t0_v6 = r_d10_c66t0_v6.numerator.real();
      if(fabs(dd0_w00t0_w0) > Fixval){
	wrapper.TP_d10_c11t0_v1[CentralityID]            ->Fill(i, (d10_c11t0_v1            )/(dd0_w00t0_w0            ), (dd0_w00t0_w0            ));
	wrapper.TP_d10_c22t0_v2[CentralityID]            ->Fill(i, (d10_c22t0_v2            )/(dd0_w00t0_w0            ), (dd0_w00t0_w0            ));
	wrapper.TP_d10_c33t0_v3[CentralityID]            ->Fill(i, (d10_c33t0_v3            )/(dd0_w00t0_w0            ), (dd0_w00t0_w0            ));
	wrapper.TP_d10_c44t0_v4[CentralityID]            ->Fill(i, (d10_c44t0_v4            )/(dd0_w00t0_w0            ), (dd0_w00t0_w0            ));
	wrapper.TP_d10_c55t0_v5[CentralityID]            ->Fill(i, (d10_c55t0_v5            )/(dd0_w00t0_w0            ), (dd0_w00t0_w0            ));
	wrapper.TP_d10_c66t0_v6[CentralityID]            ->Fill(i, (d10_c66t0_v6            )/(dd0_w00t0_w0            ), (dd0_w00t0_w0            ));
      }
      //----------------------------------------------------------------------------------
      const auto r_dd00_w000t0_w0 = qcorrEvent.differential_one_same({0,  0,  0}, {0}, i, {false, false, false});
      const auto r_d100_c211t0_v2 = qcorrEvent.differential_one_same({2, -1, -1}, {0}, i, {false, FMpTW, FMpTW});
      const auto r_d100_c312t0_v3 = qcorrEvent.differential_one_same({3, -1, -2}, {0}, i, {false, FMpTW, false});
      const auto r_d100_c413t0_v4 = qcorrEvent.differential_one_same({4, -1, -3}, {0}, i, {false, FMpTW, false});
      const auto r_d100_c514t0_v5 = qcorrEvent.differential_one_same({5, -1, -4}, {0}, i, {false, FMpTW, false});
      const auto r_d100_c422t0_v4 = qcorrEvent.differential_one_same({4, -2, -2}, {0}, i, {false, false, false});
      const auto r_d100_c523t0_v5 = qcorrEvent.differential_one_same({5, -2, -3}, {0}, i, {false, false, false});
      const auto r_d100_c624t0_v6 = qcorrEvent.differential_one_same({6, -2, -4}, {0}, i, {false, false, false});
      const auto r_d100_c633t0_v6 = qcorrEvent.differential_one_same({6, -3, -3}, {0}, i, {false, false, false});
      const double dd00_w000t0_w0 = r_dd00_w000t0_w0.denominator;
      const double d100_c211t0_v2 = r_d100_c211t0_v2.numerator.real();
      const double d100_c312t0_v3 = r_d100_c312t0_v3.numerator.real();
      const double d100_c413t0_v4 = r_d100_c413t0_v4.numerator.real();
      const double d100_c514t0_v5 = r_d100_c514t0_v5.numerator.real();
      const double d100_c422t0_v4 = r_d100_c422t0_v4.numerator.real();
      const double d100_c523t0_v5 = r_d100_c523t0_v5.numerator.real();
      const double d100_c624t0_v6 = r_d100_c624t0_v6.numerator.real();
      const double d100_c633t0_v6 = r_d100_c633t0_v6.numerator.real();
      if(fabs(dd00_w000t0_w0) > Fixval){
	wrapper.TP_d100_c211t0_v2[CentralityID]          ->Fill(i, (d100_c211t0_v2          )/(dd00_w000t0_w0          ), (dd00_w000t0_w0          ));
	wrapper.TP_d100_c312t0_v3[CentralityID]          ->Fill(i, (d100_c312t0_v3          )/(dd00_w000t0_w0          ), (dd00_w000t0_w0          ));
	wrapper.TP_d100_c413t0_v4[CentralityID]          ->Fill(i, (d100_c413t0_v4          )/(dd00_w000t0_w0          ), (dd00_w000t0_w0          ));
	wrapper.TP_d100_c514t0_v5[CentralityID]          ->Fill(i, (d100_c514t0_v5          )/(dd00_w000t0_w0          ), (dd00_w000t0_w0          ));
	wrapper.TP_d100_c422t0_v4[CentralityID]          ->Fill(i, (d100_c422t0_v4          )/(dd00_w000t0_w0          ), (dd00_w000t0_w0          ));
	wrapper.TP_d100_c523t0_v5[CentralityID]          ->Fill(i, (d100_c523t0_v5          )/(dd00_w000t0_w0          ), (dd00_w000t0_w0          ));
	wrapper.TP_d100_c624t0_v6[CentralityID]          ->Fill(i, (d100_c624t0_v6          )/(dd00_w000t0_w0          ), (dd00_w000t0_w0          ));
	wrapper.TP_d100_c633t0_v6[CentralityID]          ->Fill(i, (d100_c633t0_v6          )/(dd00_w000t0_w0          ), (dd00_w000t0_w0          ));
      }
      //----------------------------------------------------------------------------------
      const auto r_dd000_w0000t0_w0 = qcorrEvent.differential_one_same({ 0,  0,  0,  0}, {0}, i, {false, false, false, false});
      const auto r_d1000_c1111t0_v1 = qcorrEvent.differential_one_same({ 1,  1, -1, -1}, {0}, i, {false, FMpTW, FMpTW, FMpTW});
      const auto r_d1000_c2222t0_v2 = qcorrEvent.differential_one_same({ 2,  2, -2, -2}, {0}, i, {false, false, false, false});
      const auto r_d1000_c3333t0_v3 = qcorrEvent.differential_one_same({ 3,  3, -3, -3}, {0}, i, {false, false, false, false});
      const auto r_d1000_c1212t0_v1 = qcorrEvent.differential_one_same({ 1,  2, -1, -2}, {0}, i, {false, false, FMpTW, false});
      const auto r_d1000_c1313t0_v1 = qcorrEvent.differential_one_same({ 1,  3, -1, -3}, {0}, i, {false, false, FMpTW, false});
      const auto r_d1000_c1414t0_v1 = qcorrEvent.differential_one_same({ 1,  4, -1, -4}, {0}, i, {false, false, FMpTW, false});
      const auto r_d1000_c1515t0_v1 = qcorrEvent.differential_one_same({ 1,  5, -1, -5}, {0}, i, {false, false, FMpTW, false});
      const auto r_d1000_c2323t0_v2 = qcorrEvent.differential_one_same({ 2,  3, -2, -3}, {0}, i, {false, false, false, false});
      const auto r_d1000_c2424t0_v2 = qcorrEvent.differential_one_same({ 2,  4, -2, -4}, {0}, i, {false, false, false, false});
      const auto r_d1000_c2525t0_v2 = qcorrEvent.differential_one_same({ 2,  5, -2, -5}, {0}, i, {false, false, false, false});
      const auto r_d1000_c3434t0_v3 = qcorrEvent.differential_one_same({ 3,  4, -3, -4}, {0}, i, {false, false, false, false});
      const auto r_d1000_c3535t0_v3 = qcorrEvent.differential_one_same({ 3,  5, -3, -5}, {0}, i, {false, false, false, false});
      const auto r_d1000_c3111t0_v3 = qcorrEvent.differential_one_same({ 3, -1, -1, -1}, {0}, i, {false, FMpTW, FMpTW, FMpTW});
      const auto r_d1000_c6222t0_v6 = qcorrEvent.differential_one_same({ 6, -2, -2, -2}, {0}, i, {false, false, false, false});
      const auto r_d1000_c3324t0_v3 = qcorrEvent.differential_one_same({ 3,  3, -2, -4}, {0}, i, {false, false, false, false});
      const auto r_d1000_c4435t0_v4 = qcorrEvent.differential_one_same({ 4,  4, -3, -5}, {0}, i, {false, false, false, false});
      const auto r_d1000_c2534t0_v2 = qcorrEvent.differential_one_same({ 2,  5, -3, -4}, {0}, i, {false, false, false, false});
      const auto r_d0200_c1212t0_v2 = qcorrEvent.differential_one_same({ 2,  1, -1, -2}, {0}, i, {false, FMpTW, FMpTW, false});
      const auto r_d0200_c1313t0_v3 = qcorrEvent.differential_one_same({ 3,  1, -1, -3}, {0}, i, {false, FMpTW, FMpTW, false});
      const auto r_d0200_c1414t0_v4 = qcorrEvent.differential_one_same({ 4,  1, -1, -4}, {0}, i, {false, FMpTW, FMpTW, false});
      const auto r_d0200_c1515t0_v5 = qcorrEvent.differential_one_same({ 5,  1, -1, -5}, {0}, i, {false, FMpTW, FMpTW, false});
      const auto r_d0200_c2323t0_v3 = qcorrEvent.differential_one_same({ 3,  2, -2, -3}, {0}, i, {false, false, false, false});
      const auto r_d0200_c2424t0_v4 = qcorrEvent.differential_one_same({ 4,  2, -2, -4}, {0}, i, {false, false, false, false});
      const auto r_d0200_c2525t0_v5 = qcorrEvent.differential_one_same({ 5,  2, -2, -5}, {0}, i, {false, false, false, false});
      const auto r_d0200_c3434t0_v4 = qcorrEvent.differential_one_same({ 4,  3, -3, -4}, {0}, i, {false, false, false, false});
      const auto r_d0200_c3535t0_v5 = qcorrEvent.differential_one_same({ 5,  3, -3, -5}, {0}, i, {false, false, false, false});
      const auto r_d0030_c3324t0_v2 = qcorrEvent.differential_one_same({-2,  3,  3, -4}, {0}, i, {false, false, false, false});
      const auto r_d0030_c4435t0_v3 = qcorrEvent.differential_one_same({-3,  4,  4, -5}, {0}, i, {false, false, false, false});
      const auto r_d0004_c4435t0_v5 = qcorrEvent.differential_one_same({-5, -3,  4,  4}, {0}, i, {false, false, false, false});
      const auto r_d0200_c2534t0_v5 = qcorrEvent.differential_one_same({ 5,  2, -3, -4}, {0}, i, {false, false, false, false});
      const auto r_d0030_c2534t0_v3 = qcorrEvent.differential_one_same({-3,  2,  5, -4}, {0}, i, {false, false, false, false});
      const auto r_d0004_c2534t0_v4 = qcorrEvent.differential_one_same({-4,  5, -3,  2}, {0}, i, {false, false, false, false});
      
      
      const double dd000_w0000t0_w0 = r_dd000_w0000t0_w0.denominator;
      const double d1000_c1111t0_v1 = r_d1000_c1111t0_v1.numerator.real();
      const double d1000_c2222t0_v2 = r_d1000_c2222t0_v2.numerator.real();
      const double d1000_c3333t0_v3 = r_d1000_c3333t0_v3.numerator.real();
      const double d1000_c1212t0_v1 = r_d1000_c1212t0_v1.numerator.real();
      const double d1000_c1313t0_v1 = r_d1000_c1313t0_v1.numerator.real();
      const double d1000_c1414t0_v1 = r_d1000_c1414t0_v1.numerator.real();
      const double d1000_c1515t0_v1 = r_d1000_c1515t0_v1.numerator.real();
      const double d1000_c2323t0_v2 = r_d1000_c2323t0_v2.numerator.real();
      const double d1000_c2424t0_v2 = r_d1000_c2424t0_v2.numerator.real();
      const double d1000_c2525t0_v2 = r_d1000_c2525t0_v2.numerator.real();
      const double d1000_c3434t0_v3 = r_d1000_c3434t0_v3.numerator.real();
      const double d1000_c3535t0_v3 = r_d1000_c3535t0_v3.numerator.real();
      const double d1000_c3111t0_v3 = r_d1000_c3111t0_v3.numerator.real();
      const double d1000_c6222t0_v6 = r_d1000_c6222t0_v6.numerator.real();
      const double d1000_c3324t0_v3 = r_d1000_c3324t0_v3.numerator.real();
      const double d1000_c4435t0_v4 = r_d1000_c4435t0_v4.numerator.real();
      const double d1000_c2534t0_v2 = r_d1000_c2534t0_v2.numerator.real();
      const double d0200_c1212t0_v2 = r_d0200_c1212t0_v2.numerator.real();
      const double d0200_c1313t0_v3 = r_d0200_c1313t0_v3.numerator.real();
      const double d0200_c1414t0_v4 = r_d0200_c1414t0_v4.numerator.real();
      const double d0200_c1515t0_v5 = r_d0200_c1515t0_v5.numerator.real();
      const double d0200_c2323t0_v3 = r_d0200_c2323t0_v3.numerator.real();
      const double d0200_c2424t0_v4 = r_d0200_c2424t0_v4.numerator.real();
      const double d0200_c2525t0_v5 = r_d0200_c2525t0_v5.numerator.real();
      const double d0200_c3434t0_v4 = r_d0200_c3434t0_v4.numerator.real();
      const double d0200_c3535t0_v5 = r_d0200_c3535t0_v5.numerator.real();
      const double d0030_c3324t0_v2 = r_d0030_c3324t0_v2.numerator.real();
      const double d0030_c4435t0_v3 = r_d0030_c4435t0_v3.numerator.real();
      const double d0004_c4435t0_v5 = r_d0004_c4435t0_v5.numerator.real();
      const double d0200_c2534t0_v5 = r_d0200_c2534t0_v5.numerator.real();
      const double d0030_c2534t0_v3 = r_d0030_c2534t0_v3.numerator.real();
      const double d0004_c2534t0_v4 = r_d0004_c2534t0_v4.numerator.real();
      if(fabs(dd000_w0000t0_w0) > Fixval){
	wrapper.TP_d1000_c1111t0_v1[CentralityID]        ->Fill(i, (d1000_c1111t0_v1        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
	wrapper.TP_d1000_c2222t0_v2[CentralityID]        ->Fill(i, (d1000_c2222t0_v2        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
	wrapper.TP_d1000_c3333t0_v3[CentralityID]        ->Fill(i, (d1000_c3333t0_v3        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
	wrapper.TP_d1000_c1212t0_v1[CentralityID]        ->Fill(i, (d1000_c1212t0_v1        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
	wrapper.TP_d1000_c1313t0_v1[CentralityID]        ->Fill(i, (d1000_c1313t0_v1        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
	wrapper.TP_d1000_c1414t0_v1[CentralityID]        ->Fill(i, (d1000_c1414t0_v1        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
	wrapper.TP_d1000_c1515t0_v1[CentralityID]        ->Fill(i, (d1000_c1515t0_v1        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
	wrapper.TP_d1000_c2323t0_v2[CentralityID]        ->Fill(i, (d1000_c2323t0_v2        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
	wrapper.TP_d1000_c2424t0_v2[CentralityID]        ->Fill(i, (d1000_c2424t0_v2        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
	wrapper.TP_d1000_c2525t0_v2[CentralityID]        ->Fill(i, (d1000_c2525t0_v2        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
	wrapper.TP_d1000_c3434t0_v3[CentralityID]        ->Fill(i, (d1000_c3434t0_v3        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
	wrapper.TP_d1000_c3535t0_v3[CentralityID]        ->Fill(i, (d1000_c3535t0_v3        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
	wrapper.TP_d1000_c3111t0_v3[CentralityID]        ->Fill(i, (d1000_c3111t0_v3        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
	wrapper.TP_d1000_c6222t0_v6[CentralityID]        ->Fill(i, (d1000_c6222t0_v6        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
	wrapper.TP_d1000_c3324t0_v3[CentralityID]        ->Fill(i, (d1000_c3324t0_v3        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
	wrapper.TP_d1000_c4435t0_v4[CentralityID]        ->Fill(i, (d1000_c4435t0_v4        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
	wrapper.TP_d1000_c2534t0_v2[CentralityID]        ->Fill(i, (d1000_c2534t0_v2        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
	wrapper.TP_d0200_c1212t0_v2[CentralityID]        ->Fill(i, (d0200_c1212t0_v2        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
	wrapper.TP_d0200_c1313t0_v3[CentralityID]        ->Fill(i, (d0200_c1313t0_v3        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
	wrapper.TP_d0200_c1414t0_v4[CentralityID]        ->Fill(i, (d0200_c1414t0_v4        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
	wrapper.TP_d0200_c1515t0_v5[CentralityID]        ->Fill(i, (d0200_c1515t0_v5        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
	wrapper.TP_d0200_c2323t0_v3[CentralityID]        ->Fill(i, (d0200_c2323t0_v3        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
	wrapper.TP_d0200_c2424t0_v4[CentralityID]        ->Fill(i, (d0200_c2424t0_v4        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
	wrapper.TP_d0200_c2525t0_v5[CentralityID]        ->Fill(i, (d0200_c2525t0_v5        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
	wrapper.TP_d0200_c3434t0_v4[CentralityID]        ->Fill(i, (d0200_c3434t0_v4        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
	wrapper.TP_d0200_c3535t0_v5[CentralityID]        ->Fill(i, (d0200_c3535t0_v5        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
	wrapper.TP_d0030_c3324t0_v2[CentralityID]        ->Fill(i, (d0030_c3324t0_v2        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
	wrapper.TP_d0030_c4435t0_v3[CentralityID]        ->Fill(i, (d0030_c4435t0_v3        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
	wrapper.TP_d0004_c4435t0_v5[CentralityID]        ->Fill(i, (d0004_c4435t0_v5        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
	wrapper.TP_d0200_c2534t0_v5[CentralityID]        ->Fill(i, (d0200_c2534t0_v5        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
	wrapper.TP_d0030_c2534t0_v3[CentralityID]        ->Fill(i, (d0030_c2534t0_v3        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
	wrapper.TP_d0004_c2534t0_v4[CentralityID]        ->Fill(i, (d0004_c2534t0_v4        )/(dd000_w0000t0_w0        ), (dd000_w0000t0_w0        ));
      }
      //----------------------------------------------------------------------------------
      const auto r_dd0000_w00000t0_w0 = qcorrEvent.differential_one_same({0, 0,  0,  0,  0}, {0}, i, {false, false, false, false, false});
      const auto r_d10000_c21111t0_v2 = qcorrEvent.differential_one_same({2, 1, -1, -1, -1}, {0}, i, {false, FMpTW, FMpTW, FMpTW, FMpTW});
      const auto r_d10000_c42222t0_v4 = qcorrEvent.differential_one_same({4, 2, -2, -2, -2}, {0}, i, {false, false, false, false, false});
      const auto r_d10000_c43322t0_v4 = qcorrEvent.differential_one_same({4, 3, -3, -2, -2}, {0}, i, {false, false, false, false, false});
      const auto r_d10000_c52232t0_v5 = qcorrEvent.differential_one_same({5, 2, -2, -3, -2}, {0}, i, {false, false, false, false, false});
      const auto r_d10000_c53332t0_v5 = qcorrEvent.differential_one_same({5, 3, -3, -3, -2}, {0}, i, {false, false, false, false, false});
      const auto r_d10000_c22435t0_v2 = qcorrEvent.differential_one_same({2, 2,  4, -3, -5}, {0}, i, {false, false, false, false, false});
      const auto r_d10000_c23344t0_v2 = qcorrEvent.differential_one_same({2, 3,  3, -4, -4}, {0}, i, {false, false, false, false, false});
      const auto r_d10000_c33345t0_v3 = qcorrEvent.differential_one_same({3, 3,  3, -4, -5}, {0}, i, {false, false, false, false, false});
      const auto r_d10000_c24455t0_v2 = qcorrEvent.differential_one_same({2, 4,  4, -5, -5}, {0}, i, {false, false, false, false, false});
      const auto r_d10000_c33455t0_v3 = qcorrEvent.differential_one_same({3, 3,  4, -5, -5}, {0}, i, {false, false, false, false, false});
      const auto r_d10000_c22233t0_v2 = qcorrEvent.differential_one_same({2, 2,  2, -3, -3}, {0}, i, {false, false, false, false, false});
      const auto r_d10000_c62233t0_v6 = qcorrEvent.differential_one_same({6, 2, -2, -3, -3}, {0}, i, {false, false, false, false, false});
      const double dd0000_w00000t0_w0 = r_dd0000_w00000t0_w0.denominator;
      const double d10000_c21111t0_v2 = r_d10000_c21111t0_v2.numerator.real();
      const double d10000_c42222t0_v4 = r_d10000_c42222t0_v4.numerator.real();
      const double d10000_c43322t0_v4 = r_d10000_c43322t0_v4.numerator.real();
      const double d10000_c52232t0_v5 = r_d10000_c52232t0_v5.numerator.real();
      const double d10000_c53332t0_v5 = r_d10000_c53332t0_v5.numerator.real();
      const double d10000_c22435t0_v2 = r_d10000_c22435t0_v2.numerator.real();
      const double d10000_c23344t0_v2 = r_d10000_c23344t0_v2.numerator.real();
      const double d10000_c33345t0_v3 = r_d10000_c33345t0_v3.numerator.real();
      const double d10000_c24455t0_v2 = r_d10000_c24455t0_v2.numerator.real();
      const double d10000_c33455t0_v3 = r_d10000_c33455t0_v3.numerator.real();
      const double d10000_c22233t0_v2 = r_d10000_c22233t0_v2.numerator.real();
      const double d10000_c62233t0_v6 = r_d10000_c62233t0_v6.numerator.real();
      if(fabs(dd0000_w00000t0_w0) > Fixval){
	wrapper.TP_d10000_c21111t0_v2[CentralityID]      ->Fill(i, (d10000_c21111t0_v2      )/(dd0000_w00000t0_w0      ), (dd0000_w00000t0_w0      ));
	wrapper.TP_d10000_c42222t0_v4[CentralityID]      ->Fill(i, (d10000_c42222t0_v4      )/(dd0000_w00000t0_w0      ), (dd0000_w00000t0_w0      ));
	wrapper.TP_d10000_c43322t0_v4[CentralityID]      ->Fill(i, (d10000_c43322t0_v4      )/(dd0000_w00000t0_w0      ), (dd0000_w00000t0_w0      ));
	wrapper.TP_d10000_c52232t0_v5[CentralityID]      ->Fill(i, (d10000_c52232t0_v5      )/(dd0000_w00000t0_w0      ), (dd0000_w00000t0_w0      ));
	wrapper.TP_d10000_c53332t0_v5[CentralityID]      ->Fill(i, (d10000_c53332t0_v5      )/(dd0000_w00000t0_w0      ), (dd0000_w00000t0_w0      ));
	wrapper.TP_d10000_c22435t0_v2[CentralityID]      ->Fill(i, (d10000_c22435t0_v2      )/(dd0000_w00000t0_w0      ), (dd0000_w00000t0_w0      ));
	wrapper.TP_d10000_c23344t0_v2[CentralityID]      ->Fill(i, (d10000_c23344t0_v2      )/(dd0000_w00000t0_w0      ), (dd0000_w00000t0_w0      ));
	wrapper.TP_d10000_c33345t0_v3[CentralityID]      ->Fill(i, (d10000_c33345t0_v3      )/(dd0000_w00000t0_w0      ), (dd0000_w00000t0_w0      ));
	wrapper.TP_d10000_c24455t0_v2[CentralityID]      ->Fill(i, (d10000_c24455t0_v2      )/(dd0000_w00000t0_w0      ), (dd0000_w00000t0_w0      ));
	wrapper.TP_d10000_c33455t0_v3[CentralityID]      ->Fill(i, (d10000_c33455t0_v3      )/(dd0000_w00000t0_w0      ), (dd0000_w00000t0_w0      ));
	wrapper.TP_d10000_c22233t0_v2[CentralityID]      ->Fill(i, (d10000_c22233t0_v2      )/(dd0000_w00000t0_w0      ), (dd0000_w00000t0_w0      ));
	wrapper.TP_d10000_c62233t0_v6[CentralityID]      ->Fill(i, (d10000_c62233t0_v6      )/(dd0000_w00000t0_w0      ), (dd0000_w00000t0_w0      ));
      }
      //----------------------------------------------------------------------------------
      const auto r_dd00000_w000000t0_w0 = qcorrEvent.differential_one_same({0,  0, 0,  0,  0,  0}, {0}, i, {false, false, false, false, false, false});
      const auto r_d100000_c111111t0_v1 = qcorrEvent.differential_one_same({1,  1, 1, -1, -1, -1}, {0}, i, {false, FMpTW, FMpTW, FMpTW, FMpTW, FMpTW});
      const auto r_d100000_c222222t0_v2 = qcorrEvent.differential_one_same({2,  2, 2, -2, -2, -2}, {0}, i, {false, false, false, false, false, false});
      const auto r_d100000_c333333t0_v3 = qcorrEvent.differential_one_same({3,  3, 3, -3, -3, -3}, {0}, i, {false, false, false, false, false, false});
      const auto r_d100000_c223223t0_v2 = qcorrEvent.differential_one_same({2,  2, 3, -2, -2, -3}, {0}, i, {false, false, false, false, false, false});
      const auto r_d100000_c224224t0_v2 = qcorrEvent.differential_one_same({2,  2, 4, -2, -2, -4}, {0}, i, {false, false, false, false, false, false});
      const auto r_d100000_c225225t0_v2 = qcorrEvent.differential_one_same({2,  2, 5, -2, -2, -5}, {0}, i, {false, false, false, false, false, false});
      const auto r_d100000_c332332t0_v3 = qcorrEvent.differential_one_same({3,  3, 2, -3, -3, -2}, {0}, i, {false, false, false, false, false, false});
      const auto r_d100000_c334334t0_v3 = qcorrEvent.differential_one_same({3,  3, 4, -3, -3, -4}, {0}, i, {false, false, false, false, false, false});
      const auto r_d100000_c442442t0_v4 = qcorrEvent.differential_one_same({4,  4, 2, -4, -4, -2}, {0}, i, {false, false, false, false, false, false});
      const auto r_d100000_c234234t0_v2 = qcorrEvent.differential_one_same({2, -2, 3, -3,  4, -4}, {0}, i, {false, false, false, false, false, false});
      const auto r_d100000_c235235t0_v2 = qcorrEvent.differential_one_same({2, -2, 3, -3,  5, -5}, {0}, i, {false, false, false, false, false, false});
      const auto r_d100000_c345345t0_v3 = qcorrEvent.differential_one_same({3, -3, 4, -4,  5, -5}, {0}, i, {false, false, false, false, false, false});
      const auto r_d020000_c223223t0_v3 = qcorrEvent.differential_one_same({3,  2, 2, -2, -2, -3}, {0}, i, {false, false, false, false, false, false});
      const auto r_d020000_c224224t0_v4 = qcorrEvent.differential_one_same({4,  2, 2, -2, -2, -4}, {0}, i, {false, false, false, false, false, false});
      const auto r_d020000_c225225t0_v5 = qcorrEvent.differential_one_same({5,  2, 2, -2, -2, -5}, {0}, i, {false, false, false, false, false, false});
      const auto r_d020000_c332332t0_v2 = qcorrEvent.differential_one_same({2,  3, 3, -3, -3, -2}, {0}, i, {false, false, false, false, false, false});
      const auto r_d020000_c334334t0_v4 = qcorrEvent.differential_one_same({4,  3, 3, -3, -3, -4}, {0}, i, {false, false, false, false, false, false});
      const auto r_d020000_c442442t0_v2 = qcorrEvent.differential_one_same({2,  4, 4, -4, -4, -2}, {0}, i, {false, false, false, false, false, false});
      const auto r_d020000_c234234t0_v3 = qcorrEvent.differential_one_same({3,  2, 4, -2, -3, -4}, {0}, i, {false, false, false, false, false, false});
      const auto r_d003000_c234234t0_v4 = qcorrEvent.differential_one_same({4,  2, 3, -2, -3, -4}, {0}, i, {false, false, false, false, false, false});
      const auto r_d020000_c235235t0_v3 = qcorrEvent.differential_one_same({3,  2, 5, -2, -3, -5}, {0}, i, {false, false, false, false, false, false});
      const auto r_d003000_c235235t0_v5 = qcorrEvent.differential_one_same({5,  2, 3, -2, -3, -5}, {0}, i, {false, false, false, false, false, false});
      const auto r_d020000_c345345t0_v4 = qcorrEvent.differential_one_same({4,  3, 5, -3, -4, -5}, {0}, i, {false, false, false, false, false, false});
      const auto r_d003000_c345345t0_v5 = qcorrEvent.differential_one_same({5,  4, 3, -3, -4, -5}, {0}, i, {false, false, false, false, false, false});
      const double dd00000_w000000t0_w0 = r_dd00000_w000000t0_w0.denominator;
      const double d100000_c111111t0_v1 = r_d100000_c111111t0_v1.numerator.real();
      const double d100000_c222222t0_v2 = r_d100000_c222222t0_v2.numerator.real();
      const double d100000_c333333t0_v3 = r_d100000_c333333t0_v3.numerator.real();
      const double d100000_c223223t0_v2 = r_d100000_c223223t0_v2.numerator.real();
      const double d100000_c224224t0_v2 = r_d100000_c224224t0_v2.numerator.real();
      const double d100000_c225225t0_v2 = r_d100000_c225225t0_v2.numerator.real();
      const double d100000_c332332t0_v3 = r_d100000_c332332t0_v3.numerator.real();
      const double d100000_c334334t0_v3 = r_d100000_c334334t0_v3.numerator.real();
      const double d100000_c442442t0_v4 = r_d100000_c442442t0_v4.numerator.real();
      const double d100000_c234234t0_v2 = r_d100000_c234234t0_v2.numerator.real();
      const double d100000_c235235t0_v2 = r_d100000_c235235t0_v2.numerator.real();
      const double d100000_c345345t0_v3 = r_d100000_c345345t0_v3.numerator.real();
      const double d020000_c223223t0_v3 = r_d020000_c223223t0_v3.numerator.real();
      const double d020000_c224224t0_v4 = r_d020000_c224224t0_v4.numerator.real();
      const double d020000_c225225t0_v5 = r_d020000_c225225t0_v5.numerator.real();
      const double d020000_c332332t0_v2 = r_d020000_c332332t0_v2.numerator.real();
      const double d020000_c334334t0_v4 = r_d020000_c334334t0_v4.numerator.real();
      const double d020000_c442442t0_v2 = r_d020000_c442442t0_v2.numerator.real();
      const double d020000_c234234t0_v3 = r_d020000_c234234t0_v3.numerator.real();
      const double d003000_c234234t0_v4 = r_d003000_c234234t0_v4.numerator.real();
      const double d020000_c235235t0_v3 = r_d020000_c235235t0_v3.numerator.real();
      const double d003000_c235235t0_v5 = r_d003000_c235235t0_v5.numerator.real();
      const double d020000_c345345t0_v4 = r_d020000_c345345t0_v4.numerator.real();
      const double d003000_c345345t0_v5 = r_d003000_c345345t0_v5.numerator.real();
      if(fabs(dd00000_w000000t0_w0) > Fixval){
	wrapper.TP_d100000_c111111t0_v1[CentralityID]    ->Fill(i, (d100000_c111111t0_v1    )/(dd00000_w000000t0_w0    ), (dd00000_w000000t0_w0    ));
	wrapper.TP_d100000_c222222t0_v2[CentralityID]    ->Fill(i, (d100000_c222222t0_v2    )/(dd00000_w000000t0_w0    ), (dd00000_w000000t0_w0    ));
	wrapper.TP_d100000_c333333t0_v3[CentralityID]    ->Fill(i, (d100000_c333333t0_v3    )/(dd00000_w000000t0_w0    ), (dd00000_w000000t0_w0    ));
	wrapper.TP_d100000_c223223t0_v2[CentralityID]    ->Fill(i, (d100000_c223223t0_v2    )/(dd00000_w000000t0_w0    ), (dd00000_w000000t0_w0    ));
	wrapper.TP_d100000_c224224t0_v2[CentralityID]    ->Fill(i, (d100000_c224224t0_v2    )/(dd00000_w000000t0_w0    ), (dd00000_w000000t0_w0    ));
	wrapper.TP_d100000_c225225t0_v2[CentralityID]    ->Fill(i, (d100000_c225225t0_v2    )/(dd00000_w000000t0_w0    ), (dd00000_w000000t0_w0    ));
	wrapper.TP_d100000_c332332t0_v3[CentralityID]    ->Fill(i, (d100000_c332332t0_v3    )/(dd00000_w000000t0_w0    ), (dd00000_w000000t0_w0    ));
	wrapper.TP_d100000_c334334t0_v3[CentralityID]    ->Fill(i, (d100000_c334334t0_v3    )/(dd00000_w000000t0_w0    ), (dd00000_w000000t0_w0    ));
	wrapper.TP_d100000_c442442t0_v4[CentralityID]    ->Fill(i, (d100000_c442442t0_v4    )/(dd00000_w000000t0_w0    ), (dd00000_w000000t0_w0    ));
	wrapper.TP_d100000_c234234t0_v2[CentralityID]    ->Fill(i, (d100000_c234234t0_v2    )/(dd00000_w000000t0_w0    ), (dd00000_w000000t0_w0    ));
	wrapper.TP_d100000_c235235t0_v2[CentralityID]    ->Fill(i, (d100000_c235235t0_v2    )/(dd00000_w000000t0_w0    ), (dd00000_w000000t0_w0    ));
	wrapper.TP_d100000_c345345t0_v3[CentralityID]    ->Fill(i, (d100000_c345345t0_v3    )/(dd00000_w000000t0_w0    ), (dd00000_w000000t0_w0    ));
	wrapper.TP_d020000_c223223t0_v3[CentralityID]    ->Fill(i, (d020000_c223223t0_v3    )/(dd00000_w000000t0_w0    ), (dd00000_w000000t0_w0    ));
	wrapper.TP_d020000_c224224t0_v4[CentralityID]    ->Fill(i, (d020000_c224224t0_v4    )/(dd00000_w000000t0_w0    ), (dd00000_w000000t0_w0    ));
	wrapper.TP_d020000_c225225t0_v5[CentralityID]    ->Fill(i, (d020000_c225225t0_v5    )/(dd00000_w000000t0_w0    ), (dd00000_w000000t0_w0    ));
	wrapper.TP_d020000_c332332t0_v2[CentralityID]    ->Fill(i, (d020000_c332332t0_v2    )/(dd00000_w000000t0_w0    ), (dd00000_w000000t0_w0    ));
	wrapper.TP_d020000_c334334t0_v4[CentralityID]    ->Fill(i, (d020000_c334334t0_v4    )/(dd00000_w000000t0_w0    ), (dd00000_w000000t0_w0    ));
	wrapper.TP_d020000_c442442t0_v2[CentralityID]    ->Fill(i, (d020000_c442442t0_v2    )/(dd00000_w000000t0_w0    ), (dd00000_w000000t0_w0    ));
	wrapper.TP_d020000_c234234t0_v3[CentralityID]    ->Fill(i, (d020000_c234234t0_v3    )/(dd00000_w000000t0_w0    ), (dd00000_w000000t0_w0    ));
	wrapper.TP_d003000_c234234t0_v4[CentralityID]    ->Fill(i, (d003000_c234234t0_v4    )/(dd00000_w000000t0_w0    ), (dd00000_w000000t0_w0    ));
	wrapper.TP_d020000_c235235t0_v3[CentralityID]    ->Fill(i, (d020000_c235235t0_v3    )/(dd00000_w000000t0_w0    ), (dd00000_w000000t0_w0    ));
	wrapper.TP_d003000_c235235t0_v5[CentralityID]    ->Fill(i, (d003000_c235235t0_v5    )/(dd00000_w000000t0_w0    ), (dd00000_w000000t0_w0    ));
	wrapper.TP_d020000_c345345t0_v4[CentralityID]    ->Fill(i, (d020000_c345345t0_v4    )/(dd00000_w000000t0_w0    ), (dd00000_w000000t0_w0    ));
	wrapper.TP_d003000_c345345t0_v5[CentralityID]    ->Fill(i, (d003000_c345345t0_v5    )/(dd00000_w000000t0_w0    ), (dd00000_w000000t0_w0    ));
      }
      //----------------------------------------------------------------------------------
      //------------------------------------Two-Subevents---------------------------------
      //----------------------------------------------------------------------------------
      const auto r_dd0_w00sbc_w0 = qcorrEvent.differential_one_two_subevent({0,  0}, {1}, i, {QS_A, QS_B}, {false, false});
      const auto r_d10_c11sbc_v1 = qcorrEvent.differential_one_two_subevent({1, -1}, {1}, i, {QS_A, QS_B}, {FMpTW, false});
      const auto r_d10_c22sbc_v2 = qcorrEvent.differential_one_two_subevent({2, -2}, {1}, i, {QS_A, QS_B}, {false, false});
      const auto r_d10_c33sbc_v3 = qcorrEvent.differential_one_two_subevent({3, -3}, {1}, i, {QS_A, QS_B}, {false, false});
      const auto r_d10_c44sbc_v4 = qcorrEvent.differential_one_two_subevent({4, -4}, {1}, i, {QS_A, QS_B}, {false, false});
      const auto r_d10_c55sbc_v5 = qcorrEvent.differential_one_two_subevent({5, -5}, {1}, i, {QS_A, QS_B}, {false, false});
      const auto r_d10_c66sbc_v6 = qcorrEvent.differential_one_two_subevent({6, -6}, {1}, i, {QS_A, QS_B}, {false, false});
      const double dd0_w00sbc_w0 = r_dd0_w00sbc_w0.denominator;
      const double d10_c11sbc_v1 = r_d10_c11sbc_v1.numerator.real();
      const double d10_c22sbc_v2 = r_d10_c22sbc_v2.numerator.real();
      const double d10_c33sbc_v3 = r_d10_c33sbc_v3.numerator.real();
      const double d10_c44sbc_v4 = r_d10_c44sbc_v4.numerator.real();
      const double d10_c55sbc_v5 = r_d10_c55sbc_v5.numerator.real();
      const double d10_c66sbc_v6 = r_d10_c66sbc_v6.numerator.real();
      if(fabs(dd0_w00sbc_w0) > Fixval){
	wrapper.TP_d10_c11sbc_v1[CentralityID]           ->Fill(i, (d10_c11sbc_v1           )/(dd0_w00sbc_w0           ), (dd0_w00sbc_w0           ));
	wrapper.TP_d10_c22sbc_v2[CentralityID]           ->Fill(i, (d10_c22sbc_v2           )/(dd0_w00sbc_w0           ), (dd0_w00sbc_w0           ));
	wrapper.TP_d10_c33sbc_v3[CentralityID]           ->Fill(i, (d10_c33sbc_v3           )/(dd0_w00sbc_w0           ), (dd0_w00sbc_w0           ));
	wrapper.TP_d10_c44sbc_v4[CentralityID]           ->Fill(i, (d10_c44sbc_v4           )/(dd0_w00sbc_w0           ), (dd0_w00sbc_w0           ));
	wrapper.TP_d10_c55sbc_v5[CentralityID]           ->Fill(i, (d10_c55sbc_v5           )/(dd0_w00sbc_w0           ), (dd0_w00sbc_w0           ));
	wrapper.TP_d10_c66sbc_v6[CentralityID]           ->Fill(i, (d10_c66sbc_v6           )/(dd0_w00sbc_w0           ), (dd0_w00sbc_w0           ));
      }
      //----------------------------------------------------------------------------------
      const auto r_dd00_w000sbc_w0 = qcorrEvent.differential_one_two_subevent({ 0,  0, 0}, {2}, i, {QS_A, QS_A, QS_B}, {false, false, false});
      const auto r_d100_c211sbc_v2 = qcorrEvent.differential_one_two_subevent({-1, -1, 2}, {2}, i, {QS_A, QS_A, QS_B}, {FMpTW, FMpTW, false});
      const auto r_d100_c312sbc_v3 = qcorrEvent.differential_one_two_subevent({-1, -2, 3}, {2}, i, {QS_A, QS_A, QS_B}, {FMpTW, false, false});
      const auto r_d100_c413sbc_v4 = qcorrEvent.differential_one_two_subevent({-1, -3, 4}, {2}, i, {QS_A, QS_A, QS_B}, {FMpTW, false, false});
      const auto r_d100_c514sbc_v5 = qcorrEvent.differential_one_two_subevent({-1, -4, 5}, {2}, i, {QS_A, QS_A, QS_B}, {FMpTW, false, false});
      const auto r_d100_c422sbc_v4 = qcorrEvent.differential_one_two_subevent({-2, -2, 4}, {2}, i, {QS_A, QS_A, QS_B}, {false, false, false});
      const auto r_d100_c523sbc_v5 = qcorrEvent.differential_one_two_subevent({-2, -3, 5}, {2}, i, {QS_A, QS_A, QS_B}, {false, false, false});
      const auto r_d100_c624sbc_v6 = qcorrEvent.differential_one_two_subevent({-2, -4, 6}, {2}, i, {QS_A, QS_A, QS_B}, {false, false, false});
      const auto r_d100_c633sbc_v6 = qcorrEvent.differential_one_two_subevent({-3, -3, 6}, {2}, i, {QS_A, QS_A, QS_B}, {false, false, false});
      const double dd00_w000sbc_w0 = r_dd00_w000sbc_w0.denominator;
      const double d100_c211sbc_v2 = r_d100_c211sbc_v2.numerator.real();
      const double d100_c312sbc_v3 = r_d100_c312sbc_v3.numerator.real();
      const double d100_c413sbc_v4 = r_d100_c413sbc_v4.numerator.real();
      const double d100_c514sbc_v5 = r_d100_c514sbc_v5.numerator.real();
      const double d100_c422sbc_v4 = r_d100_c422sbc_v4.numerator.real();
      const double d100_c523sbc_v5 = r_d100_c523sbc_v5.numerator.real();
      const double d100_c624sbc_v6 = r_d100_c624sbc_v6.numerator.real();
      const double d100_c633sbc_v6 = r_d100_c633sbc_v6.numerator.real();
      if(fabs(dd00_w000sbc_w0) > Fixval){
	wrapper.TP_d100_c211sbc_v2[CentralityID]         ->Fill(i, (d100_c211sbc_v2         )/(dd00_w000sbc_w0         ), (dd00_w000sbc_w0         ));
	wrapper.TP_d100_c312sbc_v3[CentralityID]         ->Fill(i, (d100_c312sbc_v3         )/(dd00_w000sbc_w0         ), (dd00_w000sbc_w0         ));
	wrapper.TP_d100_c413sbc_v4[CentralityID]         ->Fill(i, (d100_c413sbc_v4         )/(dd00_w000sbc_w0         ), (dd00_w000sbc_w0         ));
	wrapper.TP_d100_c514sbc_v5[CentralityID]         ->Fill(i, (d100_c514sbc_v5         )/(dd00_w000sbc_w0         ), (dd00_w000sbc_w0         ));
	wrapper.TP_d100_c422sbc_v4[CentralityID]         ->Fill(i, (d100_c422sbc_v4         )/(dd00_w000sbc_w0         ), (dd00_w000sbc_w0         ));
	wrapper.TP_d100_c523sbc_v5[CentralityID]         ->Fill(i, (d100_c523sbc_v5         )/(dd00_w000sbc_w0         ), (dd00_w000sbc_w0         ));
	wrapper.TP_d100_c624sbc_v6[CentralityID]         ->Fill(i, (d100_c624sbc_v6         )/(dd00_w000sbc_w0         ), (dd00_w000sbc_w0         ));
	wrapper.TP_d100_c633sbc_v6[CentralityID]         ->Fill(i, (d100_c633sbc_v6         )/(dd00_w000sbc_w0         ), (dd00_w000sbc_w0         ));
      }
      //----------------------------------------------------------------------------------
      const auto r_dd000_w0000sbc_w0 = qcorrEvent.differential_one_two_subevent({ 0,  0,  0,  0}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, false, false, false});
      const auto r_d1000_c1111sbc_v1 = qcorrEvent.differential_one_two_subevent({ 1,  1, -1, -1}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, FMpTW, FMpTW, FMpTW});
      const auto r_d1000_c2222sbc_v2 = qcorrEvent.differential_one_two_subevent({ 2,  2, -2, -2}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, false, false, false});
      const auto r_d1000_c3333sbc_v3 = qcorrEvent.differential_one_two_subevent({ 3,  3, -3, -3}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, false, false, false});
      const auto r_d1000_c1212sbc_v1 = qcorrEvent.differential_one_two_subevent({ 1,  2, -1, -2}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, false, FMpTW, false});
      const auto r_d1000_c1313sbc_v1 = qcorrEvent.differential_one_two_subevent({ 1,  3, -1, -3}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, false, FMpTW, false});
      const auto r_d1000_c1414sbc_v1 = qcorrEvent.differential_one_two_subevent({ 1,  4, -1, -4}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, false, FMpTW, false});
      const auto r_d1000_c1515sbc_v1 = qcorrEvent.differential_one_two_subevent({ 1,  5, -1, -5}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, false, FMpTW, false});
      const auto r_d1000_c2323sbc_v2 = qcorrEvent.differential_one_two_subevent({ 2,  3, -2, -3}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, false, false, false});
      const auto r_d1000_c2424sbc_v2 = qcorrEvent.differential_one_two_subevent({ 2,  4, -2, -4}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, false, false, false});
      const auto r_d1000_c2525sbc_v2 = qcorrEvent.differential_one_two_subevent({ 2,  5, -2, -5}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, false, false, false});
      const auto r_d1000_c3434sbc_v3 = qcorrEvent.differential_one_two_subevent({ 3,  4, -3, -4}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, false, false, false});
      const auto r_d1000_c3535sbc_v3 = qcorrEvent.differential_one_two_subevent({ 3,  5, -3, -5}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, false, false, false});
      const auto r_d1000_c3111sbc_v3 = qcorrEvent.differential_one_two_subevent({ 3, -1, -1, -1}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, FMpTW, FMpTW, FMpTW});
      const auto r_d1000_c6222sbc_v6 = qcorrEvent.differential_one_two_subevent({ 6, -2, -2, -2}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, false, false, false});
      const auto r_d1000_c3324sbc_v3 = qcorrEvent.differential_one_two_subevent({ 3,  3, -2, -4}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, false, false, false});
      const auto r_d1000_c4435sbc_v4 = qcorrEvent.differential_one_two_subevent({ 4,  4, -3, -5}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, false, false, false});
      const auto r_d1000_c2534sbc_v2 = qcorrEvent.differential_one_two_subevent({ 2,  5, -3, -4}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, false, false, false});
      const auto r_d0200_c1212sbc_v2 = qcorrEvent.differential_one_two_subevent({ 2,  1, -1, -2}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, FMpTW, FMpTW, false});
      const auto r_d0200_c1313sbc_v3 = qcorrEvent.differential_one_two_subevent({ 3,  1, -1, -3}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, FMpTW, FMpTW, false});
      const auto r_d0200_c1414sbc_v4 = qcorrEvent.differential_one_two_subevent({ 4,  1, -1, -4}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, FMpTW, FMpTW, false});
      const auto r_d0200_c1515sbc_v5 = qcorrEvent.differential_one_two_subevent({ 5,  1, -1, -5}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, FMpTW, FMpTW, false});
      const auto r_d0200_c2323sbc_v3 = qcorrEvent.differential_one_two_subevent({ 3,  2, -2, -3}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, false, false, false});
      const auto r_d0200_c2424sbc_v4 = qcorrEvent.differential_one_two_subevent({ 4,  2, -2, -4}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, false, false, false});
      const auto r_d0200_c2525sbc_v5 = qcorrEvent.differential_one_two_subevent({ 5,  2, -2, -5}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, false, false, false});
      const auto r_d0200_c3434sbc_v4 = qcorrEvent.differential_one_two_subevent({ 4,  3, -3, -4}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, false, false, false});
      const auto r_d0200_c3535sbc_v5 = qcorrEvent.differential_one_two_subevent({ 5,  3, -3, -5}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, false, false, false});
      const auto r_d0030_c3324sbc_v2 = qcorrEvent.differential_one_two_subevent({-2,  3,  3, -4}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, false, false, false});
      const auto r_d0030_c4435sbc_v3 = qcorrEvent.differential_one_two_subevent({-3,  4,  4, -5}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, false, false, false});
      const auto r_d0004_c4435sbc_v5 = qcorrEvent.differential_one_two_subevent({-5, -3,  4,  4}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, false, false, false});
      const auto r_d0200_c2534sbc_v5 = qcorrEvent.differential_one_two_subevent({ 5,  2, -3, -4}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, false, false, false});
      const auto r_d0030_c2534sbc_v3 = qcorrEvent.differential_one_two_subevent({-3,  2,  5, -4}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, false, false, false});
      const auto r_d0004_c2534sbc_v4 = qcorrEvent.differential_one_two_subevent({-4,  5, -3,  2}, {0}, i, {QS_B, QS_B, QS_A, QS_A}, {false, false, false, false});
      const double dd000_w0000sbc_w0 = r_dd000_w0000sbc_w0.denominator;
      const double d1000_c1111sbc_v1 = r_d1000_c1111sbc_v1.numerator.real();
      const double d1000_c2222sbc_v2 = r_d1000_c2222sbc_v2.numerator.real();
      const double d1000_c3333sbc_v3 = r_d1000_c3333sbc_v3.numerator.real();
      const double d1000_c1212sbc_v1 = r_d1000_c1212sbc_v1.numerator.real();
      const double d1000_c1313sbc_v1 = r_d1000_c1313sbc_v1.numerator.real();
      const double d1000_c1414sbc_v1 = r_d1000_c1414sbc_v1.numerator.real();
      const double d1000_c1515sbc_v1 = r_d1000_c1515sbc_v1.numerator.real();
      const double d1000_c2323sbc_v2 = r_d1000_c2323sbc_v2.numerator.real();
      const double d1000_c2424sbc_v2 = r_d1000_c2424sbc_v2.numerator.real();
      const double d1000_c2525sbc_v2 = r_d1000_c2525sbc_v2.numerator.real();
      const double d1000_c3434sbc_v3 = r_d1000_c3434sbc_v3.numerator.real();
      const double d1000_c3535sbc_v3 = r_d1000_c3535sbc_v3.numerator.real();
      const double d1000_c3111sbc_v3 = r_d1000_c3111sbc_v3.numerator.real();
      const double d1000_c6222sbc_v6 = r_d1000_c6222sbc_v6.numerator.real();
      const double d1000_c3324sbc_v3 = r_d1000_c3324sbc_v3.numerator.real();
      const double d1000_c4435sbc_v4 = r_d1000_c4435sbc_v4.numerator.real();
      const double d1000_c2534sbc_v2 = r_d1000_c2534sbc_v2.numerator.real();
      const double d0200_c1212sbc_v2 = r_d0200_c1212sbc_v2.numerator.real();
      const double d0200_c1313sbc_v3 = r_d0200_c1313sbc_v3.numerator.real();
      const double d0200_c1414sbc_v4 = r_d0200_c1414sbc_v4.numerator.real();
      const double d0200_c1515sbc_v5 = r_d0200_c1515sbc_v5.numerator.real();
      const double d0200_c2323sbc_v3 = r_d0200_c2323sbc_v3.numerator.real();
      const double d0200_c2424sbc_v4 = r_d0200_c2424sbc_v4.numerator.real();
      const double d0200_c2525sbc_v5 = r_d0200_c2525sbc_v5.numerator.real();
      const double d0200_c3434sbc_v4 = r_d0200_c3434sbc_v4.numerator.real();
      const double d0200_c3535sbc_v5 = r_d0200_c3535sbc_v5.numerator.real();
      const double d0030_c3324sbc_v2 = r_d0030_c3324sbc_v2.numerator.real();
      const double d0030_c4435sbc_v3 = r_d0030_c4435sbc_v3.numerator.real();
      const double d0004_c4435sbc_v5 = r_d0004_c4435sbc_v5.numerator.real();
      const double d0200_c2534sbc_v5 = r_d0200_c2534sbc_v5.numerator.real();
      const double d0030_c2534sbc_v3 = r_d0030_c2534sbc_v3.numerator.real();
      const double d0004_c2534sbc_v4 = r_d0004_c2534sbc_v4.numerator.real();
      if(fabs(dd000_w0000sbc_w0) > Fixval){
	wrapper.TP_d1000_c1111sbc_v1[CentralityID]       ->Fill(i, (d1000_c1111sbc_v1       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
	wrapper.TP_d1000_c2222sbc_v2[CentralityID]       ->Fill(i, (d1000_c2222sbc_v2       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
	wrapper.TP_d1000_c3333sbc_v3[CentralityID]       ->Fill(i, (d1000_c3333sbc_v3       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
	wrapper.TP_d1000_c1212sbc_v1[CentralityID]       ->Fill(i, (d1000_c1212sbc_v1       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
	wrapper.TP_d1000_c1313sbc_v1[CentralityID]       ->Fill(i, (d1000_c1313sbc_v1       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
	wrapper.TP_d1000_c1414sbc_v1[CentralityID]       ->Fill(i, (d1000_c1414sbc_v1       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
	wrapper.TP_d1000_c1515sbc_v1[CentralityID]       ->Fill(i, (d1000_c1515sbc_v1       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
	wrapper.TP_d1000_c2323sbc_v2[CentralityID]       ->Fill(i, (d1000_c2323sbc_v2       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
	wrapper.TP_d1000_c2424sbc_v2[CentralityID]       ->Fill(i, (d1000_c2424sbc_v2       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
	wrapper.TP_d1000_c2525sbc_v2[CentralityID]       ->Fill(i, (d1000_c2525sbc_v2       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
	wrapper.TP_d1000_c3434sbc_v3[CentralityID]       ->Fill(i, (d1000_c3434sbc_v3       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
	wrapper.TP_d1000_c3535sbc_v3[CentralityID]       ->Fill(i, (d1000_c3535sbc_v3       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
	wrapper.TP_d1000_c3111sbc_v3[CentralityID]       ->Fill(i, (d1000_c3111sbc_v3       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
	wrapper.TP_d1000_c6222sbc_v6[CentralityID]       ->Fill(i, (d1000_c6222sbc_v6       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
	wrapper.TP_d1000_c3324sbc_v3[CentralityID]       ->Fill(i, (d1000_c3324sbc_v3       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
	wrapper.TP_d1000_c4435sbc_v4[CentralityID]       ->Fill(i, (d1000_c4435sbc_v4       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
	wrapper.TP_d1000_c2534sbc_v2[CentralityID]       ->Fill(i, (d1000_c2534sbc_v2       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
	wrapper.TP_d0200_c1212sbc_v2[CentralityID]       ->Fill(i, (d0200_c1212sbc_v2       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
	wrapper.TP_d0200_c1313sbc_v3[CentralityID]       ->Fill(i, (d0200_c1313sbc_v3       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
	wrapper.TP_d0200_c1414sbc_v4[CentralityID]       ->Fill(i, (d0200_c1414sbc_v4       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
	wrapper.TP_d0200_c1515sbc_v5[CentralityID]       ->Fill(i, (d0200_c1515sbc_v5       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
	wrapper.TP_d0200_c2323sbc_v3[CentralityID]       ->Fill(i, (d0200_c2323sbc_v3       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
	wrapper.TP_d0200_c2424sbc_v4[CentralityID]       ->Fill(i, (d0200_c2424sbc_v4       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
	wrapper.TP_d0200_c2525sbc_v5[CentralityID]       ->Fill(i, (d0200_c2525sbc_v5       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
	wrapper.TP_d0200_c3434sbc_v4[CentralityID]       ->Fill(i, (d0200_c3434sbc_v4       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
	wrapper.TP_d0200_c3535sbc_v5[CentralityID]       ->Fill(i, (d0200_c3535sbc_v5       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
	wrapper.TP_d0030_c3324sbc_v2[CentralityID]       ->Fill(i, (d0030_c3324sbc_v2       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
	wrapper.TP_d0030_c4435sbc_v3[CentralityID]       ->Fill(i, (d0030_c4435sbc_v3       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
	wrapper.TP_d0004_c4435sbc_v5[CentralityID]       ->Fill(i, (d0004_c4435sbc_v5       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
	wrapper.TP_d0200_c2534sbc_v5[CentralityID]       ->Fill(i, (d0200_c2534sbc_v5       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
	wrapper.TP_d0030_c2534sbc_v3[CentralityID]       ->Fill(i, (d0030_c2534sbc_v3       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
	wrapper.TP_d0004_c2534sbc_v4[CentralityID]       ->Fill(i, (d0004_c2534sbc_v4       )/(dd000_w0000sbc_w0       ), (dd000_w0000sbc_w0       ));
      }	
      //----------------------------------------------------------------------------------
      //----------------------------------------------------------------------------------
      //----------------------------------------------------------------------------------
      //----------------------------------------------------------------------------------
    }
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    Nev_tot ++ ;
  }
  //-------------------------------------------------------------
  //-------------------------------------------------------------  
  // Reuse the SAME file the wrapper opened:
  TFile* outputFile = wrapper.OutputFile();
  if (outputFile && !outputFile->IsZombie()) {
    //-----------------------------------------------------------
    // Back to file root (or wherever you want the histos)
    outputFile->cd(); 
    //-----------------------------------------------------------    
  }
  //-------------------------------------------------------------
  
  //-------------------------------------------------------------
  std::cout<<"The end of the code"<<"\n";
  //-------------------------------------------------------------
  //-------------------------------------------------------------  
  return 0;
}
