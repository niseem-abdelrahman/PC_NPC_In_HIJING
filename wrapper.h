#pragma once

#include "TString.h"

#include <vector>
#include <cmath>
#include <string>
#include <cstring>
#include <array>


#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TFile.h"
#include <iostream>
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TMath.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TList.h"
#include "TString.h"
#include "TRandom.h" 
#include "TRandom3.h"
#include "TComplex.h"
#include "TVector2.h"
#include <fstream>
#include "TNtuple.h"
#include "TTree.h"
#include <memory>
#include <TROOT.h>
#include <TChain.h>
#include <stdexcept>


// run QA
class TH1F;
class TH2F;
class TH1D;
class TH2D;
class TH3D;
class TFile;
class TH3F;
class TProfile       ;
class TProfile2D     ;
class TProfile3D     ;
class TRandom        ;
class TRandom3       ;
class StRefMultCorr  ;
class TComplex       ;
class TVector2       ;




// RAII helper: open ROOT file, make it current, auto-Write() & Close() on destruction
class HistoAutoFile {
 public:
  explicit HistoAutoFile(const TString& fname, const char* opt = "RECREATE");
  ~HistoAutoFile();

  TFile*      File() const { return mFile.get(); }
  void        cd()         { if (mFile) mFile->cd(); }
  TDirectory* Mkdir(const char* name);

  // non-copyable
  HistoAutoFile(const HistoAutoFile&)            = delete;
  HistoAutoFile& operator=(const HistoAutoFile&) = delete;

 private:
  std::unique_ptr<TFile> mFile;
};

// Your wrapper that owns the HistoAutoFile
class Wrapper {

 public:
  //--------------------------------------------------
  Wrapper() = default;
  //--------------------------------------------------
  void        OpenOutput(const TString& fname);
  void        CloseOutput();
  void        BookHistos();   // creates them
  //--------------------------------------------------
  TDirectory* MakeDir(const char* name);
  //--------------------------------------------------
  // optional helpers
  TFile*      OutputFile() const { return mOut ? mOut->File() : nullptr; }
  void        cd()               { if (mOut) mOut->cd(); }
  //--------------------------------------------------
  // non-copyable
  Wrapper(const Wrapper&)            = delete;
  Wrapper& operator=(const Wrapper&) = delete;
  //--------------------------------------------------
  TH1F* MultA;
  TH1F* MultB;
  TH1F* ImbaA;
  TH1F* NparA;
  //-------------------------------------------------------------------
  TH1D* Etah[10];
  TH1D* Phih[10];
  TH1D* PTh[10];
  //-------------------------------------------------------------------
  TProfile* NChvsCent  ;
  TProfile* NpartvsCent;
  //-------------------------------------------------------------------
  TProfile*  TP_d00_c11t0           ;
  TProfile*  TP_d00_c22t0           ;
  TProfile*  TP_d00_c33t0           ;
  TProfile*  TP_d00_c44t0           ;
  TProfile*  TP_d00_c55t0           ;
  TProfile*  TP_d00_c66t0           ;
  //-------------------------------------------------------------------
  TProfile*  TP_d000_c211t0         ;
  TProfile*  TP_d000_c312t0         ;
  TProfile*  TP_d000_c413t0         ;
  TProfile*  TP_d000_c514t0         ;
  TProfile*  TP_d000_c422t0         ;
  TProfile*  TP_d000_c523t0         ;
  TProfile*  TP_d000_c624t0         ;
  TProfile*  TP_d000_c633t0         ;
  //-------------------------------------------------------------------
  TProfile*  TP_d0000_c1111t0       ;
  TProfile*  TP_d0000_c2222t0       ;
  TProfile*  TP_d0000_c3333t0       ;
  TProfile*  TP_d0000_c1212t0       ;
  TProfile*  TP_d0000_c1313t0       ;
  TProfile*  TP_d0000_c1414t0       ;
  TProfile*  TP_d0000_c1515t0       ;
  TProfile*  TP_d0000_c2323t0       ;
  TProfile*  TP_d0000_c2424t0       ;
  TProfile*  TP_d0000_c2525t0       ;
  TProfile*  TP_d0000_c3434t0       ;
  TProfile*  TP_d0000_c3535t0       ;
  TProfile*  TP_d0000_c3111t0       ;
  TProfile*  TP_d0000_c6222t0       ;
  TProfile*  TP_d0000_c3324t0       ;
  TProfile*  TP_d0000_c4435t0       ;
  TProfile*  TP_d0000_c2534t0       ;
  //-------------------------------------------------------------------
  TProfile*  TP_d00000_c21111t0     ;
  TProfile*  TP_d00000_c42222t0     ;
  TProfile*  TP_d00000_c43322t0     ;
  TProfile*  TP_d00000_c52232t0     ;
  TProfile*  TP_d00000_c53332t0     ;
  TProfile*  TP_d00000_c22435t0     ;
  TProfile*  TP_d00000_c23344t0     ;
  TProfile*  TP_d00000_c33345t0     ;
  TProfile*  TP_d00000_c24455t0     ;
  TProfile*  TP_d00000_c33455t0     ;
  TProfile*  TP_d00000_c22233t0     ;
  TProfile*  TP_d00000_c62233t0     ;
  //-------------------------------------------------------------------
  TProfile*  TP_d000000_c111111t0   ;
  TProfile*  TP_d000000_c222222t0   ;
  TProfile*  TP_d000000_c333333t0   ;

  TProfile* TP_d000000_c112112t0;
  TProfile* TP_d000000_c113113t0;
  TProfile* TP_d000000_c221221t0;
  TProfile* TP_d000000_c331331t0;
  TProfile* TP_d000000_c123123t0;


  TProfile*  TP_d000000_c223223t0   ;
  TProfile*  TP_d000000_c224224t0   ;
  TProfile*  TP_d000000_c225225t0   ;
  TProfile*  TP_d000000_c332332t0   ;
  TProfile*  TP_d000000_c334334t0   ;
  TProfile*  TP_d000000_c442442t0   ;
  TProfile*  TP_d000000_c234234t0   ;
  TProfile*  TP_d000000_c235235t0   ;
  TProfile*  TP_d000000_c345345t0   ;
  TProfile*  TP_d000000_c111122t0_N ;
  TProfile*  TP_d000000_c111122t0_D ;
  //-------------------------------------------------------------------
  TProfile*  TP_d00_c11sbc          ;
  TProfile*  TP_d00_c22sbc          ;
  TProfile*  TP_d00_c33sbc          ;
  TProfile*  TP_d00_c44sbc          ;
  TProfile*  TP_d00_c55sbc          ;
  TProfile*  TP_d00_c66sbc          ;
  //-------------------------------------------------------------------
  TProfile*  TP_d000_c211sbc        ;
  TProfile*  TP_d000_c312sbc        ;
  TProfile*  TP_d000_c413sbc        ;
  TProfile*  TP_d000_c514sbc        ;
  TProfile*  TP_d000_c422sbc        ;
  TProfile*  TP_d000_c523sbc        ;
  TProfile*  TP_d000_c624sbc        ;
  TProfile*  TP_d000_c633sbc        ;
  //-------------------------------------------------------------------
  TProfile*  TP_d0000_c1111sbc      ;
  TProfile*  TP_d0000_c2222sbc      ;
  TProfile*  TP_d0000_c3333sbc      ;
  TProfile*  TP_d0000_c1212sbc      ;
  TProfile*  TP_d0000_c1313sbc      ;
  TProfile*  TP_d0000_c1414sbc      ;
  TProfile*  TP_d0000_c1515sbc      ;
  TProfile*  TP_d0000_c2323sbc      ;
  TProfile*  TP_d0000_c2424sbc      ;
  TProfile*  TP_d0000_c2525sbc      ;
  TProfile*  TP_d0000_c3434sbc      ;
  TProfile*  TP_d0000_c3535sbc      ;
  TProfile*  TP_d0000_c3111sbc      ;
  TProfile*  TP_d0000_c6222sbc      ;
  TProfile*  TP_d0000_c3324sbc      ;
  TProfile*  TP_d0000_c4435sbc      ;
  TProfile*  TP_d0000_c2534sbc      ;
  //-------------------------------------------------------------------
  //-------------------------------------------------------------------
  TProfile*  TP_d00_c00t0_00[10]            ;
  TProfile*  TP_d10_c11t0_v1[10]            ;
  TProfile*  TP_d10_c22t0_v2[10]            ;
  TProfile*  TP_d10_c33t0_v3[10]            ;
  TProfile*  TP_d10_c44t0_v4[10]            ;
  TProfile*  TP_d10_c55t0_v5[10]            ;
  TProfile*  TP_d10_c66t0_v6[10]            ;
  //-------------------------------------------------------------------
  TProfile*  TP_d100_c211t0_v2[10]          ;
  TProfile*  TP_d100_c312t0_v3[10]          ;
  TProfile*  TP_d100_c413t0_v4[10]          ;
  TProfile*  TP_d100_c514t0_v5[10]          ;
  TProfile*  TP_d100_c422t0_v4[10]          ;
  TProfile*  TP_d100_c523t0_v5[10]          ;
  TProfile*  TP_d100_c624t0_v6[10]          ;
  TProfile*  TP_d100_c633t0_v6[10]          ;
  //-------------------------------------------------------------------
  TProfile*  TP_d1000_c1111t0_v1[10]        ;
  TProfile*  TP_d1000_c2222t0_v2[10]        ;
  TProfile*  TP_d1000_c3333t0_v3[10]        ;
  TProfile*  TP_d1000_c1212t0_v1[10]        ;
  TProfile*  TP_d1000_c1313t0_v1[10]        ;
  TProfile*  TP_d1000_c1414t0_v1[10]        ;
  TProfile*  TP_d1000_c1515t0_v1[10]        ;
  TProfile*  TP_d1000_c2323t0_v2[10]        ;
  TProfile*  TP_d1000_c2424t0_v2[10]        ;
  TProfile*  TP_d1000_c2525t0_v2[10]        ;
  TProfile*  TP_d1000_c3434t0_v3[10]        ;
  TProfile*  TP_d1000_c3535t0_v3[10]        ;
  TProfile*  TP_d1000_c3111t0_v3[10]        ;
  TProfile*  TP_d1000_c6222t0_v6[10]        ;
  TProfile*  TP_d1000_c3324t0_v3[10]        ;
  TProfile*  TP_d1000_c4435t0_v4[10]        ;
  TProfile*  TP_d1000_c2534t0_v2[10]        ;
  TProfile*  TP_d0200_c1212t0_v2[10]        ;
  TProfile*  TP_d0200_c1313t0_v3[10]        ;
  TProfile*  TP_d0200_c1414t0_v4[10]        ;
  TProfile*  TP_d0200_c1515t0_v5[10]        ;
  TProfile*  TP_d0200_c2323t0_v3[10]        ;
  TProfile*  TP_d0200_c2424t0_v4[10]        ;
  TProfile*  TP_d0200_c2525t0_v5[10]        ;
  TProfile*  TP_d0200_c3434t0_v4[10]        ;
  TProfile*  TP_d0200_c3535t0_v5[10]        ;
  TProfile*  TP_d0030_c3324t0_v2[10]        ;
  TProfile*  TP_d0030_c4435t0_v3[10]        ;
  TProfile*  TP_d0004_c4435t0_v5[10]        ;
  TProfile*  TP_d0200_c2534t0_v5[10]        ;
  TProfile*  TP_d0030_c2534t0_v3[10]        ;
  TProfile*  TP_d0004_c2534t0_v4[10]        ;
  //-------------------------------------------------------------------
  TProfile*  TP_d10000_c21111t0_v2[10]      ;
  TProfile*  TP_d10000_c42222t0_v4[10]      ;
  TProfile*  TP_d10000_c43322t0_v4[10]      ;
  TProfile*  TP_d10000_c52232t0_v5[10]      ;
  TProfile*  TP_d10000_c53332t0_v5[10]      ;
  TProfile*  TP_d10000_c22435t0_v2[10]      ;
  TProfile*  TP_d10000_c23344t0_v2[10]      ;
  TProfile*  TP_d10000_c33345t0_v3[10]      ;
  TProfile*  TP_d10000_c24455t0_v2[10]      ;
  TProfile*  TP_d10000_c33455t0_v3[10]      ;
  TProfile*  TP_d10000_c22233t0_v2[10]      ;
  TProfile*  TP_d10000_c62233t0_v6[10]      ;
  //-------------------------------------------------------------------
  TProfile*  TP_d100000_c111111t0_v1[10]    ;
  TProfile*  TP_d100000_c222222t0_v2[10]    ;
  TProfile*  TP_d100000_c333333t0_v3[10]    ;
  TProfile*  TP_d100000_c223223t0_v2[10]    ;
  TProfile*  TP_d100000_c224224t0_v2[10]    ;
  TProfile*  TP_d100000_c225225t0_v2[10]    ;
  TProfile*  TP_d100000_c332332t0_v3[10]    ;
  TProfile*  TP_d100000_c334334t0_v3[10]    ;
  TProfile*  TP_d100000_c442442t0_v4[10]    ;
  TProfile*  TP_d100000_c234234t0_v2[10]    ;
  TProfile*  TP_d100000_c235235t0_v2[10]    ;
  TProfile*  TP_d100000_c345345t0_v3[10]    ;
  TProfile*  TP_d020000_c223223t0_v3[10]    ;
  TProfile*  TP_d020000_c224224t0_v4[10]    ;
  TProfile*  TP_d020000_c225225t0_v5[10]    ;
  TProfile*  TP_d020000_c332332t0_v2[10]    ;
  TProfile*  TP_d020000_c334334t0_v4[10]    ;
  TProfile*  TP_d020000_c442442t0_v2[10]    ;
  TProfile*  TP_d020000_c234234t0_v3[10]    ;
  TProfile*  TP_d003000_c234234t0_v4[10]    ;
  TProfile*  TP_d020000_c235235t0_v3[10]    ;
  TProfile*  TP_d003000_c235235t0_v5[10]    ;
  TProfile*  TP_d020000_c345345t0_v4[10]    ;
  TProfile*  TP_d003000_c345345t0_v5[10]    ;
  //-------------------------------------------------------------------
  //-------------------------------------------------------------------
  TProfile*  TP_d10_c11sbc_v1[10]           ;
  TProfile*  TP_d10_c22sbc_v2[10]           ;
  TProfile*  TP_d10_c33sbc_v3[10]           ;
  TProfile*  TP_d10_c44sbc_v4[10]           ;
  TProfile*  TP_d10_c55sbc_v5[10]           ;
  TProfile*  TP_d10_c66sbc_v6[10]           ;
  //-------------------------------------------------------------------
  TProfile*  TP_d100_c211sbc_v2[10]         ;
  TProfile*  TP_d100_c312sbc_v3[10]         ;
  TProfile*  TP_d100_c413sbc_v4[10]         ;
  TProfile*  TP_d100_c514sbc_v5[10]         ;
  TProfile*  TP_d100_c422sbc_v4[10]         ;
  TProfile*  TP_d100_c523sbc_v5[10]         ;
  TProfile*  TP_d100_c624sbc_v6[10]         ;
  TProfile*  TP_d100_c633sbc_v6[10]         ;
  //-------------------------------------------------------------------
  TProfile*  TP_d1000_c1111sbc_v1[10]       ;
  TProfile*  TP_d1000_c2222sbc_v2[10]       ;
  TProfile*  TP_d1000_c3333sbc_v3[10]       ;
  TProfile*  TP_d1000_c1212sbc_v1[10]       ;
  TProfile*  TP_d1000_c1313sbc_v1[10]       ;
  TProfile*  TP_d1000_c1414sbc_v1[10]       ;
  TProfile*  TP_d1000_c1515sbc_v1[10]       ;
  TProfile*  TP_d1000_c2323sbc_v2[10]       ;
  TProfile*  TP_d1000_c2424sbc_v2[10]       ;
  TProfile*  TP_d1000_c2525sbc_v2[10]       ;
  TProfile*  TP_d1000_c3434sbc_v3[10]       ;
  TProfile*  TP_d1000_c3535sbc_v3[10]       ;
  TProfile*  TP_d1000_c3111sbc_v3[10]       ;
  TProfile*  TP_d1000_c6222sbc_v6[10]       ;
  TProfile*  TP_d1000_c3324sbc_v3[10]       ;
  TProfile*  TP_d1000_c4435sbc_v4[10]       ;
  TProfile*  TP_d1000_c2534sbc_v2[10]       ;
  TProfile*  TP_d0200_c1212sbc_v2[10]       ;
  TProfile*  TP_d0200_c1313sbc_v3[10]       ;
  TProfile*  TP_d0200_c1414sbc_v4[10]       ;
  TProfile*  TP_d0200_c1515sbc_v5[10]       ;
  TProfile*  TP_d0200_c2323sbc_v3[10]       ;
  TProfile*  TP_d0200_c2424sbc_v4[10]       ;
  TProfile*  TP_d0200_c2525sbc_v5[10]       ;
  TProfile*  TP_d0200_c3434sbc_v4[10]       ;
  TProfile*  TP_d0200_c3535sbc_v5[10]       ;
  TProfile*  TP_d0030_c3324sbc_v2[10]       ;
  TProfile*  TP_d0030_c4435sbc_v3[10]       ;
  TProfile*  TP_d0004_c4435sbc_v5[10]       ;
  TProfile*  TP_d0200_c2534sbc_v5[10]       ;
  TProfile*  TP_d0030_c2534sbc_v3[10]       ;
  TProfile*  TP_d0004_c2534sbc_v4[10]       ;
  //-------------------------------------------------------------------
 private:
  //-------------------------------------------------------------------
  std::unique_ptr<HistoAutoFile> mOut;
  //-------------------------------------------------------------------
};
