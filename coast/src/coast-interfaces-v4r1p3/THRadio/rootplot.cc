/***************************************************************************
 *   Copyright (C) 2006 by Tim Huege, Sven Lafebre                         *
 *   tim.huege@ik.fzk.de, s.lafebre@astro.ru.nl                            *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/*
  - Wed Sep 23 22:42:49 CEST 2009, Ralf Ulrich
    code maintenance and cleanup
 */

/* TODO
 * ****
 * - Error bars for cuts
 * - Average Energy and offset plots
 * - Y axis ranges
 * - Boundaries for cuts
 * - Relative e+/e- normalization when cuts are applied
 */

#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TFrame.h>
#include <TGraphErrors.h>
#include <TH3D.h>
#include <TH2D.h>
#include <TPostScript.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>

#include <cmath>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>

#include <popt.h>

using namespace std;


/// global variable
bool debug = false;


/* 
 * definition of plot commands
 */
typedef enum {
  TP_GRAPH,
  TP_ERRORAREA,
  TP_ERRORLIMIT,
  TP_AXISRANGE,
} TPlotType;

const double SPEED_OF_LIGHT = 2.99792458; // m/s
const int PADSX = 2;
const int PADSY = 3;


/// ************************************************************************
/// ************************************************************************
///
/// \struct Parameters
/// \brief Run parameters
///
/// This <tt>struct</tt> describes the run parameters of a rootplot run.

struct Parameters {
  string file;
  string particle;
  double layer;
  double age;
  double depth;
  double rdepth;

  int    linecolor;
  int    linestyle;
  int    linewidth;

  Parameters();
  bool        operator==(const Parameters&) const;
  bool        operator!=(const Parameters&) const;
  bool        operator<(const Parameters&) const;
  bool        operator>(const Parameters&) const;
  bool        OnlyDiffersByParticle(const Parameters&);
  void   Print();
};

Parameters::Parameters() 
: file(""),
  particle(""),
  layer(-1),
  age(-1),
  depth(-1),
  rdepth(-1e6),
  linecolor(1),
  linestyle(1),
  linewidth(1) {
}

bool Parameters::operator>(const Parameters& that) const {
  if (file      > that.file) return true;
  if (particle  > that.particle) return true;
  if (layer     > that.layer) return true;
  if (age       > that.age) return true;
  if (depth     > that.depth) return true;
  if (rdepth    > that.rdepth) return true;
  if (linecolor > that.linecolor) return true;
  if (linestyle > that.linestyle) return true;
  if (linewidth > that.linewidth) return true;

  return false;
}

bool Parameters::operator<(const Parameters& that) const {
  if (file      < that.file) return true;
  if (particle  < that.particle) return true;
  if (layer     < that.layer) return true;
  if (age       < that.age) return true;
  if (depth     < that.depth) return true;
  if (rdepth    < that.rdepth) return true;
  if (linecolor < that.linecolor) return true;
  if (linestyle < that.linestyle) return true;
  if (linewidth < that.linewidth) return true;

  return false;
}

bool Parameters::operator==(const Parameters& that) const {
  return (file      == that.file &&
	  particle  == that.particle &&
	  layer     == that.layer &&
	  age       == that.age &&
	  depth     == that.depth &&
	  rdepth    == that.rdepth);
}

bool Parameters::operator!=(const Parameters& that) const {
  return !(*this == that);
}

bool Parameters::OnlyDiffersByParticle(const Parameters& that) {
  return (file      == that.file &&
	  layer     == that.layer &&
	  age       == that.age &&
	  depth     == that.depth &&
	  rdepth    == that.rdepth);
}

void Parameters::Print() {
  cout << "-- Parameters --------" << endl
       << "File     :" << file << endl
       << "Particle :" << particle << endl;
  if (layer > 0) {
    cout
      << "Layer    :" << layer << endl;
  } else if (age > 0) {
    cout
      << "Age      :" << age << endl;
  } else if (depth > 0) {
    cout
      << "Depth    :" << depth << endl;
  } else if (rdepth > -1e5) {
    cout
      << "Rel.Depth:" << rdepth << endl;
  }
}

/// ************************************************************************
/// ************************************************************************
///
/// \struct Offset
/// \brief Structure to hold the offset data
///

struct Offset {
  // Array of offset values
  TArrayD val;
  // Array of square of offset errors
  TArrayD err;
  // Array of number of entries for value
  TArrayI has;
  string name;
  string title;

  Offset();
  Offset(unsigned int n);
  Offset& Add(const Offset&);

  unsigned int GetSize() const {return has.GetSize();}
  unsigned int GetN()    const {return has.GetSize();}
};

Offset::Offset(unsigned int n) {
  val.Set(n);
  err.Set(n);
  has.Set(n);
  for (unsigned int i=0; i<n; i++) {
    has[i] = 0;
  }
}

Offset::Offset() {
  Offset(0);
}

Offset& Offset::Add(const Offset& that) {
  if (that.GetSize() != GetSize()) {
    cerr << "Warning: adding offset histograms of different size" << endl;
    return *this;
  }
  
  for (unsigned int i=0; i<GetSize(); ++i) {
    if (has[i] && that.has[i]) {
      val[i] += that.val[i];
      err[i] += that.err[i];
      has[i] += that.has[i];
    } else if (that.has[i]) {
      val[i] = that.val[i];
      err[i] = that.err[i];
      has[i] = that.has[i];
    }
  }
  
  return *this;
}

/// ************************************************************************
/// ************************************************************************
///
/// \struct Dataset
/// \brief Data structure to hold data from the files
///

struct Dataset {
  
  TH3D   time;
  TH3D   time_err;
  TH3D   angle;
  TH3D   angle_err;

  Offset offsetX;
  Offset offsetY;

  TH2D   energy;

  double norm;

  bool   isallocated;
  bool   hasoffset;
  bool   hasenergy;
  bool   haserror;
  int    numberofstacks;

  Dataset();
  //Dataset(const Dataset&);
  void Add(const Dataset& ds, bool stack, double weight);
};

Dataset::Dataset() {
  isallocated    = false;
  hasoffset      = false;
  hasenergy      = false;
  haserror       = false;
  numberofstacks = 1;
}

/*
 *
 */
void Dataset::Add(const Dataset& ds, bool stack = false, double weight = 1.) {
  
  if (! isallocated) {
    
    if (debug) 
      cout << "Add             : New average histogram" << endl;
    
    if (weight != 1.)
      cerr << "Warning: histograms can only be initialized with weight 1." << endl;
    
    isallocated = true;
    haserror    = true;

    time      = ds.time;
    angle     = ds.angle;
    norm      = ds.norm;
    if (ds.hasoffset) {
      offsetX   = ds.offsetX;
      offsetY   = ds.offsetY;
      hasoffset = true;
    }

    time_err  = ds.time;
    angle_err = ds.angle;
    time_err. Multiply(&ds.time);
    angle_err.Multiply(&ds.angle);
  } else {
    time. Add(&ds.time, weight);
    angle.Add(&ds.angle, weight);

    if (hasoffset && ds.hasoffset) {
      offsetX.Add(ds.offsetX);
      offsetY.Add(ds.offsetY);
    }

    if (!stack) {
      norm = 1./(1./norm+1./(weight*ds.norm));

      TH3D* tmpT = (TH3D*)ds.time.Clone("tmpThst");
      tmpT->Multiply(&ds.time);
      time_err.Add(tmpT, weight);

      TH3D* tmpA = (TH3D*)ds.angle.Clone("tmpAhst");
      tmpA->Multiply(&ds.angle);
      angle_err.Add(tmpA, weight);

      numberofstacks ++;

      tmpA->Delete();
      tmpT->Delete();
    }
  }

  if (ds.hasenergy) {
    if (! hasenergy) {
      energy = ds.energy;
      hasenergy = true;
    } else {
      energy.Add(&ds.energy);
    }
  }

  if (debug) 
    cout << "Add             : Adding histogram" << endl;

}

//Dataset::Dataset(const Dataset& that) {
//time           = that.time;
//angle          = that.angle;
//isallocated    = that.isallocated;
//hasoffset      = that.hasoffset;
//numberofstacks = that.numberofstacks;
//}

/// ************************************************************************
/// ************************************************************************
///
/// \struct Options
/// \brief Store the main programs parameters
///

struct Options {

  string outfile;
  int    particles;
  int    ecuts;
  int    lcuts;
  int    tcuts;
  int    acuts;
  int    zcuts;
  int    numfiles;
  TPostScript* psstream;
  
  bool   firstplot;
  
  bool   average;
  bool   errors;
  bool   logy;
  bool   normalize;
  bool   newpage;
  int    postscript;
  bool   textout;
};

/// ************************************************************************
/// ************************************************************************
///
/// \struct TAxisRange
/// \brief Store the min/max of an axis range
///

struct TAxisRange {
  double min;
  double max;
  bool isinitialized;

  TAxisRange() { isinitialized = false; }
};

/// -----------------------------------------------------------------------
/// -----------------------------------------------------------------------
///
///  Initialize global variables
///
/// -----------------------------------------------------------------------

Options options;
TCanvas c1("c1", "Rootplot", 400, 540);
int numberofplots = 0;
TAxisRange axisrange[PADSX*PADSY+2];
vector<Parameters> plots;


/* OpenFile
 * ********
 */
TFile* OpenFile(const string& filename) {
  
  if (debug) 
    cout << "OpenFile        : Opening file " << filename << endl;

  TFile* f = new TFile(filename.c_str());
  if (f->IsZombie()) {
    cout << "Error opening file: \"" << filename << "\"" << endl;
    exit(-1);
  }
  return f;
}

/* AgeToLayer
 * **********
 */
double AgeToLayer(TTree* t, Float_t s) {
  int    i  = 0;
  Float_t s1 = 0, s2 = -1;
  t->SetBranchStatus("*", 0); //disable all branches
  t->SetBranchStatus("age", 1);
  t->SetBranchAddress("age", &s2);

  while (s2 < s && i < t->GetEntries()) {
    s1 = s2;
    t->GetEntry(i++);
  }
  i -= 2;
  if (i < 2)
    return (double)i;

  if (debug) 
    cout << "AgeToLayer      : Layer " << i + (s1 - s)/(s1 - s2) << endl;

  // Perform linear interpolation to find good value for observation level
  return i + (s1 - s)/(s1 - s2);
}

/* DepthToLayer
 * **********
 */
double DepthToLayer(TFile* f, Float_t X) {
  TTree* t = (TTree*)f->Get("shower");
  Float_t Xmax;
  
  t->SetBranchStatus("*", 0); //disable all branches
  t->SetBranchStatus("xmax", 1);
  t->SetBranchAddress("xmax", &Xmax);
  t->GetEntry(0);
  
  return AgeToLayer((TTree*)f->Get("data_layer"), 3.*X/(X+2.*Xmax));
}

/*
 *
 */
double RdepthToLayer(TFile* f, Float_t X) {
  TTree* t = (TTree*)f->Get("shower");
  Float_t Xmax;

  t->SetBranchStatus("*",0); //disable all branches
  t->SetBranchStatus("xmax",1);
  t->SetBranchAddress("xmax", &Xmax);
  t->GetEntry(0);

  X *= 36.7;
  X += Xmax;

  return AgeToLayer((TTree*)f->Get("data_layer"), 3.*X/(X+2.*Xmax));
}

/*
 *
 */
Offset InterpolateOffsetData(TGraphErrors* l, TGraphErrors* h,
			     double dl,
			     Dataset& ds,
			     string name, string title) {
  
  Offset lo(ds.time.GetNbinsZ());
  Offset hi(ds.time.GetNbinsZ());

  hi.name  = name;
  hi.title = title;

  const int nBins = ds.angle.GetNbinsZ();
  const double xMin  = ds.angle.GetZaxis()->GetXmin();
  const double xMax  = ds.angle.GetZaxis()->GetXmax();

  for (int i=0; i<h->GetN(); ++i) {
    const int bin = int((nBins-1) * (h->GetX()[i] - xMin) / (xMax - xMin));
    hi.val[bin] = h->GetY ()[i];
    hi.err[bin] = h->GetEY()[i] * h->GetEY()[i];
    hi.has[bin] = true;
  }
  for (int i=0; i<l->GetN(); ++i) {
    const int bin = int((nBins-1) * (l->GetX()[i] - xMin) / (xMax - xMin));
    lo.val[bin] = l->GetY ()[i];
    lo.err[bin] = l->GetEY()[i] * l->GetEY()[i];
    lo.has[bin] = true;
  }
  for (unsigned int i=0; i<lo.GetSize(); ++i) {
    if (lo.has[i] && hi.has[i]) {
      hi.val[i] = lo.val[i] * (1-dl)+hi.val[i]*dl;
      hi.err[i] = lo.err[i] * (1-dl)+hi.err[i]*dl;
    } else if (lo.has[i]) {
      if (dl < 0.5) {
	hi.val[i] = lo.val[i];
	hi.err[i] = lo.err[i];
	hi.has[i] = true;
      }
    } else if (hi.has[i]) {
      if (dl < 0.5) {
	hi.has[i] = false;
      }
    }
  }
  return hi;
}

/*
 *
 */
TH2D InterpolateEnergyData(TH2D* l, TH2D* h,
			   double dl,
			   string name, string title) {
  TH2D tgt(name.c_str(), title.c_str(),
	   l->GetNbinsX(), l->GetXaxis()->GetXmin(), l->GetXaxis()->GetXmax(),
	   l->GetNbinsY(), l->GetYaxis()->GetXmin(), l->GetYaxis()->GetXmax());

  for (int i = 0; i < l->GetNbinsX() + 2; i ++) {
    for (int j = 0; j < l->GetNbinsY() + 2; j ++) {
      tgt.SetBinContent(i,j,
			(l->GetBinContent(i,j)*(1-dl) + h->GetBinContent(i,j)*dl));// /
      //pow(10,l->GetBinCenter(l->GetBin(i,j))));
    }	}
  return tgt;
}

/* InterpolateLayerData
 *
 */
TH3D InterpolateLayerData(TH3D* l, TH3D* h,
                          double dl,
                          double rml, double old_rml,
                          double rmh, double old_rmh,
                          string name, string title) {
  double dil = log10(rml/old_rml)/(l->GetXaxis()->GetXmax()-l->GetXaxis()->GetXmin())*l->GetNbinsX();
  double djl = log10(rml/old_rml)/(l->GetYaxis()->GetXmax()-l->GetYaxis()->GetXmin())*l->GetNbinsY();
  double dih = log10(rmh/old_rmh)/(h->GetXaxis()->GetXmax()-h->GetXaxis()->GetXmin())*h->GetNbinsX();
  double djh = log10(rmh/old_rmh)/(h->GetYaxis()->GetXmax()-h->GetYaxis()->GetXmin())*h->GetNbinsY();

  // Make new histogram
  TH3D tgt(name.c_str(), title.c_str(),
	   l->GetNbinsX(), l->GetXaxis()->GetXmin()-log10(old_rml*SPEED_OF_LIGHT),
	   l->GetXaxis()->GetXmax()-log10(old_rml*SPEED_OF_LIGHT),
	   l->GetNbinsY(), l->GetYaxis()->GetXmin(),
	   l->GetYaxis()->GetXmax(),
	   l->GetNbinsZ(), l->GetZaxis()->GetXmin(),
	   l->GetZaxis()->GetXmax());

  for (int i = 0; i < l->GetNbinsX() + 2; i ++) {
    for (int j = 0; j < l->GetNbinsY() + 2; j ++) {
      for (int k = 0; k < l->GetNbinsZ() + 2; k ++) {
	tgt.SetBinContent(i,j,k,
			  l->GetBinContent(i+  (int)floor(dil),j+  (int)floor(djl),k)*(1+floor(dil)-dil)*(1+floor(djl)-djl)*(1-dl)+
			  l->GetBinContent(i+  (int)floor(dil),j+1+(int)floor(djl),k)*(1+floor(dil)-dil)*  (djl-floor(djl))*(1-dl)+
			  l->GetBinContent(i+1+(int)floor(dil),j+  (int)floor(djl),k)*  (dil-floor(dil))*(1+floor(djl)-djl)*(1-dl)+
			  l->GetBinContent(i+1+(int)floor(dil),j+1+(int)floor(djl),k)*  (dil-floor(dil))*  (djl-floor(djl))*(1-dl)+
			  h->GetBinContent(i+  (int)floor(dih),j+  (int)floor(djh),k)*(1+floor(dih)-dih)*(1+floor(djh)-djh)*   dl +
			  h->GetBinContent(i+  (int)floor(dih),j+1+(int)floor(djh),k)*(1+floor(dih)-dih)*  (djh-floor(djh))*   dl +
			  h->GetBinContent(i+1+(int)floor(dih),j+  (int)floor(djh),k)*  (dih-floor(dih))*(1+floor(djh)-djh)*   dl +
			  h->GetBinContent(i+1+(int)floor(dih),j+1+(int)floor(djh),k)*  (dih-floor(dih))*  (djh-floor(djh))*   dl);
	//cout << tgt.GetBinContent(i,j,k) << " = "
	//<< lo->GetBinContent(i,j,k) << " x " << 1.-dlayer << " + "
	//<< hi->GetBinContent(i,j,k) << " x " << dlayer << endl;
      }	}	}
  return tgt;
}

/*
  Dataset Square(Dataset& ds) {	
  }*/

/* GetIntegralFromDataFile
 *
 */
double GetIntegralFromDataFile(TFile* f, double layer) {
  TTree *te = 0, *tp = 0;
  tp = (TTree*)f->Get("data_positron");
  te = (TTree*)f->Get("data_electron");

  TH3D* gete = 0;
  TH3D* getp = 0;
  te->SetBranchStatus ("*",0);
  te->SetBranchStatus ("hTelectron",1);
  tp->SetBranchStatus ("*",0);
  tp->SetBranchStatus ("hTpositron",1);
  te->SetBranchAddress("hTelectron", &gete);
  tp->SetBranchAddress("hTpositron", &getp);
  te->GetEntry((int)floor(layer));
  tp->GetEntry((int)floor(layer));
  double el = gete->Integral();
  double pl = getp->Integral();

  te->GetEntry((int)floor(layer)+1);
  tp->GetEntry((int)floor(layer)+1);
  double eh = gete->Integral();
  double ph = getp->Integral();

  return (el+pl) + (eh+ph-el-pl)*(layer-floor(layer));
}

/* GetData
 *
 */
Dataset GetData(Parameters &par) {
  
  TFile *f = OpenFile(par.file);
  TTree *h = 0;
  h = (TTree*)f->Get("data_layer");

  double layer;
  if (par.layer > 0) {
    layer = par.layer;
  } else if (par.age > 0) {
    layer = AgeToLayer(h, par.age);
  } else if (par.depth > 0) {
    layer = DepthToLayer(f, par.depth);
  } else if (par.rdepth > -1e5) {
    layer = RdepthToLayer(f, par.rdepth);
  } else {
    cerr << "No valid layer specified" << endl;
    exit(-1);
  }

  TTree *d;
  string p;

  if (par.particle[0] == 'p') {
    d = (TTree*)f->Get("data_positron");
    p = "positron";
  } else {
    d = (TTree*)f->Get("data_electron");
    p = "electron";
  }
  
  // Get the shower age if no age to plot was specified
  Float_t rml, old_rml, rmh, old_rmh, age = par.age;
  if (age < 0) {
    Float_t age1, age2;
    h->SetBranchStatus("*",0); //disable all branches
    h->SetBranchStatus("age",1);
    h->SetBranchAddress("age", &age2);
    h->GetEntry((int)floor(layer));
    age1 = age2;
    h->GetEntry((int)floor(layer)+1);
    age = age1 + (layer - floor(layer))*(age2 - age1);
  }

  // Get the correct value for the Molière adius
  Float_t tmp1, tmp2;
  h->SetBranchStatus("*",0);
  h->SetBranchStatus("density",1);
  h->SetBranchStatus("rm",1);
  h->SetBranchAddress("density", &tmp1);
  h->SetBranchAddress("rm", &tmp2);
  h->GetEntry((int)floor(layer));
  rml     = 0.096/tmp1;
  old_rml = tmp2;
  h->GetEntry((int)floor(layer)+1);
  rmh     = 0.096/tmp1;
  old_rmh = tmp2;

  stringstream tmp;
  char nameA[512], nameT[512], nameOx[512], nameOy[512], nameE[512], title[512];
  tmp << par.file << ",A," << layer << "," << p << endl;
  tmp >> nameA;
  tmp << par.file << ",T," << layer << "," << p << endl;
  tmp >> nameT;
  tmp << par.file << ",Ox," << layer << "," << p << endl;
  tmp >> nameOx;
  tmp << par.file << ",Oy," << layer << "," << p << endl;
  tmp >> nameOy;

  sprintf(title, "File %s, s = %lg (%s)", par.file.c_str(), age, p.c_str());

  if (! debug) 
    cout << title << endl;

  TH3D* getT = 0;
  TH3D* getA = 0;
  TH3D* lT; TH3D* hT;
  TH3D* lA; TH3D* hA;
  d->SetBranchStatus ("*",0);
  //d->SetBranchStatus (string("hT"+p).c_str(),1);
  d->SetBranchStatus (string("hT"+p).c_str(),1);
  d->SetBranchStatus (string("hA"+p).c_str(),1);
  //d->SetBranchAddress(string("hT"+p).c_str(), &getT);
  d->SetBranchAddress(string("hT"+p).c_str(), &getT);
  d->SetBranchAddress(string("hA"+p).c_str(), &getA);
  d->GetEntry((int)floor(layer));
  lT = (TH3D*)getT->Clone("hT_tmp1");
  lA = (TH3D*)getA->Clone("hA_tmp1");
  d->GetEntry((int)floor(layer)+1);
  hT = (TH3D*)getT->Clone("hT_tmp2");
  hA = (TH3D*)getA->Clone("hA_tmp2");
  
  Dataset ds;
  ds.norm = 1./GetIntegralFromDataFile(f,layer);
  
  ds.time  = InterpolateLayerData(lT, hT, layer-floor(layer), rml, old_rml, rmh, old_rmh, nameT, title);
  ds.angle = InterpolateLayerData(lA, hA, layer-floor(layer),  1.,      1.,  1.,      1., nameA, title);
  ds.isallocated = true;
  ds.numberofstacks = 1;

  lT->Delete();
  lA->Delete();
  hT->Delete();
  hA->Delete();

  if (options.normalize) {
    ds.time .Scale(ds.norm);
    ds.angle.Scale(ds.norm);
  }

  ds.time.GetXaxis()->SetTitle("#font[132]{log(#font[12]{ct}/#font[12]{r_{M}})}");
  ds.time.GetYaxis()->SetTitle("#font[132]{log(#font[12]{r}/#font[12]{r_{M}})}");
  ds.time.GetZaxis()->SetTitle("#font[132]{log(#font[12]{E}/GeV)}");
  ds.angle.GetXaxis()->SetTitle("#font[12]{#theta}");
  ds.angle.GetYaxis()->SetTitle("#font[12]{#phi}");
  ds.angle.GetZaxis()->SetTitle("#font[132]{log(#font[12]{E}/GeV)}");

  TH2D* lE; TH2D* hE;
  if (d->GetBranch(string("hE"+p).c_str())==NULL) {
    ds.hasenergy = false;
  } else {
    TH2D* getE = 0;
    d->SetBranchStatus ("*",0);
    d->SetBranchStatus (string("hE"+p).c_str(),1);
    d->SetBranchAddress(string("hE"+p).c_str(), &getE);
    d->GetEntry((int)floor(layer));
    lE = (TH2D*)getE->Clone("hE_tmp1");
    d->GetEntry((int)floor(layer)+1);
    hE = (TH2D*)getE->Clone("hE_tmp2");

    ds.energy = InterpolateEnergyData(lE, hE, layer-floor(layer), nameE, title);
    ds.hasenergy = true;
    ds.energy.GetXaxis()->SetTitle("#font[132]{log(#font[12]{E}_{old}/GeV)}");
    ds.energy.GetYaxis()->SetTitle("#Delta#font[12]{E}");

    lE->Delete();
    hE->Delete();
  }

  TGraphErrors *lOX = 0, *hOX = 0;
  TGraphErrors *lOY = 0, *hOY = 0;
  if (d->GetBranch(string("offsetX_"+p).c_str())==NULL ||
      d->GetBranch(string("offsetY_"+p).c_str())==NULL) {
    cerr << "Warning: no offset histogram in file" << endl;
    ds.hasoffset = false;
  } else {
    TGraphErrors *getOX = 0, *getOY = 0;
    d->SetBranchStatus ("*", 0);
    d->SetBranchStatus (string("offsetX_"+p).c_str(), 1);
    d->SetBranchStatus (string("offsetY_"+p).c_str(), 1);
    d->SetBranchAddress(string("offsetX_"+p).c_str(), &getOX);
    d->SetBranchAddress(string("offsetY_"+p).c_str(), &getOY);
    d->GetEntry((int)floor(layer));
    lOX = (TGraphErrors*)getOX->Clone("hOX_tmp1");
    lOY = (TGraphErrors*)getOY->Clone("hOY_tmp1");
    d->GetEntry((int)floor(layer)+1);
    hOX = (TGraphErrors*)getOX->Clone("hOX_tmp2");
    hOY = (TGraphErrors*)getOY->Clone("hOY_tmp2");

    ds.offsetX = InterpolateOffsetData(lOX, hOX, layer-floor(layer), ds, nameOx, title);
    ds.offsetY = InterpolateOffsetData(lOY, hOY, layer-floor(layer), ds, nameOy, title);
    ds.hasoffset = true;
    //ds.offsetX.GetXaxis()->SetTitle("#font[132]{log(#font[12]{E}_{old}/GeV)}");
    //ds.offsetX.GetYaxis()->SetTitle("#Delta#font[12]{x}");
    //ds.offsetY.GetXaxis()->SetTitle("#font[132]{log(#font[12]{E}_{old}/GeV)}");
    //ds.offsetY.GetYaxis()->SetTitle("#Delta#font[12]{y}");

    lOX->Delete();
    lOY->Delete();
    hOX->Delete();
    hOY->Delete();
  }


  f->Close();

  return ds;
}


/* CalculateErrors
 * ***************
 * Calculate the standard error of the mean, given by
 * \sigma_{\bar x}
 *    = \frac{\sigma_x}{\sqrt{N}}
 *    = \sqrt{\frac1{N-1}\left(\frac1N\sum_i x_i^2-\bar x^2\right)}
 * Multiply by three to get $3\sigma$ confidence level
 */
void CalculateErrors(Dataset& ds) {
  for (int i = 0; i < ds.time.GetNbinsX() + 2; i ++) {
    for (int j = 0; j < ds.time.GetNbinsY() + 2; j ++) {
      for (int k = 0; k < ds.time.GetNbinsZ() + 2; k ++) {
	ds.time_err.SetBinContent(i,j,k,
				  3*sqrt((   ds.time_err.GetBinContent(i,j,k)
					     - ds.time.GetBinContent(i,j,k)*ds.time.GetBinContent(i,j,k)/ds.numberofstacks)
					 /double((ds.numberofstacks-1)*ds.numberofstacks)));
      }	}	}
  for (int i = 0; i < ds.angle.GetNbinsX() + 2; i ++) {
    for (int j = 0; j < ds.angle.GetNbinsY() + 2; j ++) {
      for (int k = 0; k < ds.angle.GetNbinsZ() + 2; k ++) {
	ds.angle_err.SetBinContent(i,j,k,
				   3*sqrt((   ds.angle_err.GetBinContent(i,j,k)
					      - ds.angle.GetBinContent(i,j,k)*ds.angle.GetBinContent(i,j,k)/ds.numberofstacks)
				          /double((ds.numberofstacks-1)*ds.numberofstacks)));
      }	}	}
}

/* ParseCommaSeparatedList
 * ***********************
 * Ralf Ulrich, Wed Sep 23 22:58:15 CEST 2009
 *
 * parses a comma separated string of numbers
 * and returns a vector of these numbers
 */
template<typename T>
vector<T> ParseCommaSeparatedList(const string& str) {
  vector<T> list;
  size_t pos = 0;
  while (str != "" && pos <= str.size() && pos != string::npos) {
    const size_t comma = str.find(',', pos);
    string subStr;
    if (comma == string::npos) {
      subStr = str.substr(pos, str.length()-pos);
      pos = string::npos;
    } else {
      subStr = str.substr(pos, comma-pos);
      pos = comma + 1;
    }

    istringstream convert(subStr);
    T entry;
    convert >> entry;
    list.push_back(entry);
  }
  return list;
}
  


/* ProcessOptions
 * **************
 * Reads options using popt and processes them to fill lists of settings
 * and list of parameter cases to plot.
 */
int ProcessOptions(int argc, const char** argv) {

  char const * layerbuf     = "\0";
  char const * agebuf       = "\0";
  char const * depthbuf     = "\0";
  char const * rdepthbuf    = "\0";
  char const * particlebuf  = "e,p";
  bool  do_average   = false;
  bool  do_errors    = true;
  bool  do_noerrors    = true;
  bool  do_logy      = false;
  bool  do_newpage   = false;
  bool  do_normalize = false;
  int   do_postscript= false;
  bool  do_text      = false;
  char const * outf_name    = "rootplot.ps";
  int   t_cuts       = 1;
  int   l_cuts       = 1;
  int   a_cuts       = 1;
  int   z_cuts       = 1;
  int   e_cuts       = 1;

  // Define allowed options 
  struct poptOption optionsTable[] = {
    { "average",  'a', POPT_ARG_NONE,  &do_average,  0, "Produce one averaged histogram for all files", "\0" },
    { "debug",    'd', POPT_ARG_NONE,  &debug,       0, "Print debugging information", "\0" },
    { "noerrors", 'e', POPT_ARG_NONE,  &do_noerrors, 0, "Do not draw error margins for averaged plots", "\0" },
    { "layer",    'l', POPT_ARG_STRING,&layerbuf,    0, "Comma separated list of layers to plot", "comma-separated-list" },
    { "normalize",'n', POPT_ARG_NONE,  &do_normalize,0, "Normalize the integral of histograms to one", "\0" },
    { "newpage",  'N', POPT_ARG_NONE,  &do_newpage,  0, "Put every plot (except different particle species) on a separate page", "\0" },
    { "output",   'o', POPT_ARG_STRING,&outf_name,   0, "Filename for final output (default: 'rootplot.ps')", "file" },
    { "particle", 'p', POPT_ARG_STRING,&particlebuf, 0,
      "Comma separated list of particle types to plot; possible types are electron, positron, sum and difference. (default: e,p)", "comma-separated-list" },
    { "age",      's', POPT_ARG_STRING,&agebuf,      0, "Comma separated list of shower ages to plot", "comma-separated-list" },
    { "text",     't', POPT_ARG_NONE,  &do_text,     0, "Plain text output", "\0" },
    { "depth",    'X', POPT_ARG_STRING,&depthbuf,    0, "Comma separated list of atmospheric depths to plot", "comma-separated-list" },
    { "rdepth",   'Y', POPT_ARG_STRING,&rdepthbuf,   0, "Comma separated list of relative atmospheric depths to plot", "comma-separated-list" },
    { "logy",     'y', POPT_ARG_NONE,  &do_logy,     0, "Toggle logarithmic scale on y axes", "\0" },
    { "tcuts",   '\0', POPT_ARG_INT,   &t_cuts,      0, "Cut relevant histograms in time slices", "n" },
    { "lcuts",   '\0', POPT_ARG_INT,   &l_cuts,      0, "Cut relevant histograms in lateral slices", "n" },
    { "acuts",   '\0', POPT_ARG_INT,   &a_cuts,      0, "Cut relevant histograms in momentum azimuth angle slices", "n" },
    { "zcuts",   '\0', POPT_ARG_INT,   &z_cuts,      0, "Cut relevant histograms in momentum zenith angle slices",  "n" },
    { "ecuts",   '\0', POPT_ARG_INT,   &e_cuts,      0, "Cut relevant histograms in energy slices", "n" },
    POPT_AUTOHELP
    POPT_TABLEEND
  };
  // Use popt to parse options
  poptContext arguments = poptGetContext(NULL, argc, argv, optionsTable, 0);
  poptSetOtherOptionHelp(arguments, "<histogram file(s)>");
  
  // If no command line options were given, print usage notifications
  if (argc < 2) {
    poptPrintUsage(arguments, stderr, 0);
    return 0;
  }
  
  char c;
  while ((c = poptGetNextOpt(arguments)) >= 0) { }

  do_errors = !do_noerrors;
  
  if (t_cuts < 1) t_cuts = 1;
  if (l_cuts < 1) l_cuts = 1;
  if (a_cuts < 1) a_cuts = 1;
  if (z_cuts < 1) z_cuts = 1;
  if (e_cuts < 1) e_cuts = 1;

  // Parse layer string to make list of layers to run over
  vector<int> layers = ParseCommaSeparatedList<int>(layerbuf);
  vector<double> ages = ParseCommaSeparatedList<double>(agebuf);
  vector<double> depths = ParseCommaSeparatedList<double>(depthbuf);
  vector<double> rdepths = ParseCommaSeparatedList<double>(rdepthbuf);
  vector<string> particles = ParseCommaSeparatedList<string>(particlebuf);
  // str = "e,p"; defulat
  
  // If no histogram layer or shower age is given, do not proceed
  if (layers.size() == 0 && ages.size() == 0 &&
      depths.size() == 0 && rdepths.size() == 0) {
    cerr << "Please specify a layer to plot (or use 'rootplot -?')" << endl;
    return 0;
  }

  options.outfile = outf_name;
  options.textout = do_text;
  if (options.outfile.rfind(".ps") == options.outfile.length() - 3) {
    do_postscript = 1;
  } else if (options.outfile.rfind(".pdf") == options.outfile.length() - 4) {
    do_postscript = 2;
    options.outfile.replace(options.outfile.length() - 4, 4, ".ps");
  }
  // Fill settings object with information obtained from command line
  options.average     = do_average;
  options.errors      = do_errors;
  options.logy        = do_logy;
  options.newpage     = do_newpage;
  options.normalize   = do_normalize;
  options.postscript  = do_postscript;
  options.tcuts       = t_cuts;
  options.lcuts       = l_cuts;
  options.acuts       = a_cuts;
  options.zcuts       = z_cuts;
  options.ecuts       = e_cuts;
  options.numfiles    = 0;

  
  if (debug) {
    cout << "------------------------------------------------------------\n";
    cout << " Options: \n"
	 << " do_average:   " << do_average << "\n"
	 << " do_errors:    " << do_errors << "\n"
	 << " do_logy:      " << do_logy << "\n"
	 << " do_newpage:   " << do_newpage << "\n"
	 << " do_normalize: " << do_normalize << "\n"
	 << " do_postscript:" << do_postscript << "\n"
	 << " do_test:      " << do_text << "\n"
	 << " t_cuts:       " << t_cuts << "\n"
	 << " l_cuts:       " << l_cuts << "\n"
	 << " a_cuts:       " << a_cuts << "\n"
	 << " z_cuts:       " << z_cuts << "\n"
	 << " e_cuts:       " << e_cuts << "\n"
	 << " layers:       ";
    for (vector<int>::const_iterator i = layers.begin(); i != layers.end(); ++i) cout << *i << " "; cout << "\n";
    cout << " ages:         ";
    for (vector<double>::const_iterator i = ages.begin(); i != ages.end(); ++i) cout << *i << " "; cout << "\n";
    cout << " depths:       ";
    for (vector<double>::const_iterator i = depths.begin(); i != depths.end(); ++i) cout << *i << " "; cout << "\n";
    cout << " rdepths:      ";
    for (vector<double>::const_iterator i = rdepths.begin(); i != rdepths.end(); ++i) cout << *i << " "; cout << "\n";
    cout << " particles:    ";
    for (vector<string>::const_iterator i = particles.begin(); i != particles.end(); ++i) cout << *i << " "; cout << "\n";
    cout << "output:        " << outf_name << "\n";
    cout << "-----------------------------------------------------\n";
  }
  
  
  // Create list of parameter cases to perform
  char* a = NULL;
  int lStyle = 1;
  // Now loop over files to make an array of cases
  while ((a = (char*)poptGetArg(arguments)) != NULL) {
    
    Parameters par;
    par.file = string(a);
    
    if (debug)
      cout << " found input file: " << par.file << endl;

    //if (!options.newpage) {
    par.linestyle = lStyle;
    //}
    options.numfiles ++;

    // Define line widths for different layers, ages, depths, ...
    const int ws = layers.size() + ages.size() + depths.size() + rdepths.size();
    int lineWidth[ws];
    if (options.newpage) {
      for (int w=0; w<ws; ++w) {
	lineWidth[w] = 4;
      }
    } else {
      if (ws == 1) {
	lineWidth[0] = 4;
      } else if (ws < 5) {
	for (int w=0; w<ws; ++w) {
	  lineWidth[w] = (ws-w) * 2;
	}
      } else if (ages.size() == 5) {
	for (int w=0; w<ws; ++w) {
	  lineWidth[w] = (ws-w) * 2 - 1;
	}
      } else {
	for (int w=0; w<ws; ++w) {
	  lineWidth[w] = ws-w;
	}
      }
    }

    
    int lwCounter = 0;
    for (vector<double>::const_iterator s = ages.begin();
	 s != ages.end(); ++s) {
      par.layer     = -1;
      par.age       = *s;
      par.linewidth = lineWidth[lwCounter++];
      for (vector<string>::const_iterator p = particles.begin();
	   p != particles.end(); ++p) {
	par.particle = *p;
	// Define line colors for different particles
	if        (par.particle[0] == 'e') {
	  par.linecolor = 11;
	} else if (par.particle[0] == 'p') {
	  par.linecolor = 12;
	} else if (par.particle[0] == 's') {
	  par.linecolor = 13;
	} else {
	  par.linecolor = 14;
	}
	plots.push_back(par);
      }
    }
    
    for (vector<double>::const_iterator X = depths.begin();
	   X != depths.end(); ++X) {
      par.layer     = -1;
      par.depth     = *X;
      par.linewidth = lineWidth[lwCounter++];
      for (vector<string>::const_iterator p = particles.begin();
	   p != particles.end(); ++p) {
	par.particle = *p;
	// Define line colors for different particles
	if        (par.particle[0] == 'e') {
	  par.linecolor = 11;
	} else if (par.particle[0] == 'p') {
	  par.linecolor = 12;
	} else if (par.particle[0] == 's') {
	  par.linecolor = 13;
	} else {
	  par.linecolor = 14;
	}
	plots.push_back(par);
      }
    }

    for (vector<double>::const_iterator Y = rdepths.begin();
	 Y != rdepths.end(); ++Y) {
      par.layer     = -1;
      par.rdepth    = *Y;
      par.linewidth = lineWidth[lwCounter++];
      for (vector<string>::const_iterator p = particles.begin();
	   p != particles.end(); ++p) {
	par.particle = *p;
	// Define line colors for different particles
	if        (par.particle[0] == 'e') {
	  par.linecolor = 11;
	} else if (par.particle[0] == 'p') {
	  par.linecolor = 12;
	} else if (par.particle[0] == 's') {
	  par.linecolor = 13;
	} else {
	  par.linecolor = 14;
	}
	plots.push_back(par);
      }
    }

    for (vector<int>::const_iterator l = layers.begin();
	 l != layers.end(); ++l) {
      par.layer     = *l;
      par.age       = -1;
      par.linewidth = lineWidth[lwCounter++];
      for (vector<string>::const_iterator p = particles.begin();
	   p != particles.end(); ++p) {
	par.particle = *p;
	// Define line colors for different particles
	if        (par.particle[0] == 'e') {
	  par.linecolor = 11;
	} else if (par.particle[0] == 'p') {
	  par.linecolor = 12;
	} else if (par.particle[0] == 's') {
	  par.linecolor = 13;
	} else {
	  par.linecolor = 14;
	}
	plots.push_back(par);
      }
    }
    lStyle ++;
  } // loop input files
  
  if (options.numfiles == 1) {
    options.average = false;
  }
  
  options.firstplot = true;
  
  // We don't need the command line arguments any longer now
  poptFreeContext(arguments);

  if (debug)
    cout << "\n" << endl;
  
  cout << "Rootplot v0.724" << endl
       << "Copyright (C) 2007 by Sven Lafebre (s.lafebre@astro.ru.nl)" << endl
       << endl
       << (options.average ? "Averaging and p" : "P") << "lotting " << options.numfiles
       << (options.normalize ? " normalized" : " unnormalized") << " data file"
       << (options.numfiles > 1 ? "s" : "") << endl;
  string with = "with ";
  if (options.tcuts > 1) {
    cout << with << options.tcuts << " time cuts ";
    with = "and ";
  }
  if (options.lcuts > 1) {
    cout << with << options.lcuts << " lateral cuts ";
    with = "and ";
  }
  if (options.acuts > 1) {
    cout << with << options.acuts << " momentum azimuth angle cuts ";
    with = "and ";
  }
  if (options.zcuts > 1) {
    cout << with << options.zcuts << " momentum zenith angle cuts ";
    with = "and ";
  }
  if (options.ecuts > 1) {
    cout << with << options.ecuts << " energy cuts ";
  }
  cout << "over " << (layers.size()+ages.size()+depths.size()+rdepths.size()) << " layer"
       << (layers.size()+ages.size() == 1 ? " " : "s ")
       << "for " << options.particles << " particle type"
       << (options.particles > 1 ? "s" : "") << endl << endl;

  int num;
  if ( (num = (options.newpage ? 1 : ((options.average ? 1 : options.numfiles) * (layers.size()+ages.size())))
	* options.tcuts * options.lcuts * options.acuts * options.zcuts * options.ecuts
	* particles.size()) > 6) {
    cout << "Advice:" << endl
	 << "With your settings, graphs may contain up to " << num << " plots." << endl
	 << "Consider increasing readability by reducing the number of layers, plots," << endl
	 << (options.average ? ""  : "files, ") << "or particles"
	 << (options.newpage ? "." : ", or by using the '-N' (newpage) option.")
	 << endl << endl;
  }
  return 1;
}


/* InitPlot
 * ********
 *
 * Apply settings to make the plots look nice
 * (which is not the default seventies ROOT look)
 *  and creates the Canvas layout
 */
void InitPlot() {
  
  c1.SetFillColor(kWhite);
  gStyle->SetOptStat(0);
  gStyle->SetMarkerSize(0.5);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(2);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderSize(2);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetPadBorderSize(2);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPaperSize(29.7,21);
  gStyle->SetLineWidth(3);
  gStyle->SetLabelFont(132);
  gStyle->SetTextFont(132);
  gStyle->SetTitleFont(132);
  gStyle->SetLabelFont(132);
  gStyle->SetErrorX(0);

  //float r, g, b;
  TColor* color = (TColor*)gROOT->GetListOfColors()->At(11);
  color->SetRGB(0.90, 0.30, 0.00);
  color = (TColor*)gROOT->GetListOfColors()->At(21);
  color->SetRGB(0.95, 0.65, 0.50);
  color = (TColor*)gROOT->GetListOfColors()->At(31);
  color->SetRGB(0.99, 0.93, 0.85);
  color = (TColor*)gROOT->GetListOfColors()->At(12);
  color->SetRGB(0.30, 0.90, 0.00);
  color = (TColor*)gROOT->GetListOfColors()->At(22);
  color->SetRGB(0.65, 0.95, 0.50);
  color = (TColor*)gROOT->GetListOfColors()->At(32);
  color->SetRGB(0.93, 0.99, 0.85);
  color = (TColor*)gROOT->GetListOfColors()->At(13);
  color->SetRGB(0.00, 0.30, 0.90);
  color = (TColor*)gROOT->GetListOfColors()->At(23);
  color->SetRGB(0.50, 0.65, 0.95);
  color = (TColor*)gROOT->GetListOfColors()->At(33);
  color->SetRGB(0.85, 0.93, 0.99);
  color = (TColor*)gROOT->GetListOfColors()->At(14);
  color->SetRGB(0.30, 0.00, 0.90);
  color = (TColor*)gROOT->GetListOfColors()->At(24);
  color->SetRGB(0.65, 0.50, 0.95);
  color = (TColor*)gROOT->GetListOfColors()->At(34);
  color->SetRGB(0.93, 0.85, 0.99);

  c1.Divide(PADSX, PADSY);
  c1.cd(PADSX*PADSY)->Divide(1, 2, 0.01, 0);
  c1.cd(PADSX*PADSY)->cd(1);
  gPad->SetTopMargin(0.1);
  c1.cd(PADSX*PADSY)->cd(2);
  gPad->SetBottomMargin(0.2);
  if (options.postscript) {
    options.psstream = new TPostScript(options.outfile.c_str(), 111);
    options.psstream->Range(20., 25.9);
    options.psstream->SetLineScale(0.5);
    options.psstream->NewPage();
  } else {
    options.psstream = NULL;
  }
}


/* DrawGraph
 * *********
 * Draws a single graph
 */
void DrawGraph(Dataset& ds, TH1* h, string title, double sum, Parameters par, int style, TPlotType error, bool dotime) {

  //if (debug) cout << "DrawGraph       : Adding plot" << endl;
  TH3D* h3  = 0;
  TH3D* h3e = 0;

  if (dotime) {
    h3  = &ds.time;
    h3e = &ds.time_err;
  } else {
    h3  = &ds.angle;
    h3e = &ds.angle_err;
  }

  int mini, minj, mink, maxi, maxj, maxk;
  int minx, miny, minz, maxx, maxy, maxz;
  int  i,  j,  k;
  int *x, *y, *z;

  minx = h3->GetXaxis()->GetFirst();
  maxx = h3->GetXaxis()->GetLast();
  miny = h3->GetYaxis()->GetFirst();
  maxy = h3->GetYaxis()->GetLast();
  minz = h3->GetZaxis()->GetFirst();
  maxz = h3->GetZaxis()->GetLast();

  string s = h->GetName();
  int pos = s.rfind("_") + 1;
  if (s[pos] == 'x') {
    mini = minx; maxi = maxx;
    minj = miny; maxj = maxy;
    mink = minz; maxk = maxz;
    x = &i; y = &j; z = &k;
  } else if (s[pos] == 'y') {
    mini = miny; maxi = maxy;
    minj = minz; maxj = maxz;
    mink = minx; maxk = maxx;
    x = &k; y = &i; z = &j;
  } else if (s[pos] == 'z') {
    mini = minz; maxi = maxz;
    minj = minx; maxj = maxx;
    mink = miny; maxk = maxy;
    x = &j; y = &k; z = &i;
  } else {
    cerr << "Cannot determine errors. Skipping this plot." << endl;
    return;
  }

  double norm = 1./ds.numberofstacks;
  if (par.particle[0] == 'd') {
    norm /= (maxj-minj)*(maxk-mink);
  } else if (options.normalize) {
    norm *= sum/h->Integral();
  }
  h->Scale(norm);

  TH1* hl = (TH1*)h->Clone("tmp1");
  TH1* hh = (TH1*)h->Clone("tmp2");
  if (error == TP_ERRORAREA  ||
      error == TP_ERRORLIMIT ||
      error == TP_AXISRANGE) {
    for (i=mini; i<maxi+1; ++i) {
      double err_sum = 0;
      for (j=minj; j<maxj+1; ++j) {
	for (k=mink; k<maxk+1; ++k) {
	  err_sum += h3e->GetBinContent(*x,*y,*z)*norm;
	}	
      }
      h->SetBinError(i, 3*err_sum);
      hl->SetBinContent(i, h->GetBinContent(i) - 3*err_sum);
      hh->SetBinContent(i, h->GetBinContent(i) + 3*err_sum);
    }
  }

  hl->SetLineColor(par.linecolor+10);
  hl->SetLineStyle(par.linestyle+style);
  hl->SetLineWidth(par.linewidth);

  hh->SetLineColor(par.linecolor+10);
  hh->SetLineStyle(par.linestyle+style);
  hh->SetLineWidth(par.linewidth);

  h->SetLineColor(par.linecolor);
  h->SetFillColor(par.linecolor+20);
  h->SetFillStyle(1001*error);
  h->SetLineStyle(par.linestyle+style);
  h->SetLineWidth(par.linewidth);
  h->SetMarkerSize(0);
  h->SetLabelFont(132, "xy");
  h->SetTitle(title.c_str());

  // Determine Y axis range
  if (error == TP_AXISRANGE) {
    TH1* bl = 0;
    TH1* bh = 0;
    if (options.errors && options.average) {
      bl = hl;
      bh = hh;
    } else {
      bl = h;
      bh = h;
    }
    if (! axisrange[gPad->GetNumber()].isinitialized) {
      axisrange[gPad->GetNumber()].max = bh->GetMaximum();
      axisrange[gPad->GetNumber()].min = bl->GetMinimum();
      axisrange[gPad->GetNumber()].isinitialized = true;
    } else {
      if (axisrange[gPad->GetNumber()].max < bh->GetMaximum()) {
	axisrange[gPad->GetNumber()].max = bh->GetMaximum();
      }
      if (axisrange[gPad->GetNumber()].min > bl->GetMinimum()) {
	axisrange[gPad->GetNumber()].min = bl->GetMinimum();
      }
    }
  }
  if (options.logy) {
    gPad->SetLogy(1);
  }

  if (axisrange[gPad->GetNumber()].isinitialized) {
    if (!options.logy) {
      h->GetYaxis()->SetRangeUser(
				  axisrange[gPad->GetNumber()].max
				  - 1.1*(axisrange[gPad->GetNumber()].max-axisrange[gPad->GetNumber()].min),
				  axisrange[gPad->GetNumber()].min
				  + 1.1*(axisrange[gPad->GetNumber()].max-axisrange[gPad->GetNumber()].min));
    } else {
      h->GetYaxis()->SetRangeUser(pow(10,
				      log10(axisrange[gPad->GetNumber()].max)
				      - 1.1*(log10(axisrange[gPad->GetNumber()].max)-log10(axisrange[gPad->GetNumber()].min))),
				  pow(10,
				      log10(axisrange[gPad->GetNumber()].min)
				      + 1.1*(log10(axisrange[gPad->GetNumber()].max)-log10(axisrange[gPad->GetNumber()].min))));
      if (axisrange[gPad->GetNumber()].min <= 0) {
	h->GetYaxis()->SetRangeUser(
				    (axisrange[gPad->GetNumber()].min
				     + 1.1*(axisrange[gPad->GetNumber()].max-axisrange[gPad->GetNumber()].min))*1.e-5,
				    axisrange[gPad->GetNumber()].min
				    + 1.1*(axisrange[gPad->GetNumber()].max-axisrange[gPad->GetNumber()].min));
      }
    }
  }
  
  if (error == TP_ERRORAREA && options.average) {
    //		h->DrawClone(string("e3"+string(startnew && options.newpage ? "":" same")).c_str());
    h->DrawCopy(string("e3"+string(options.firstplot ? "":" same")).c_str());
  } else if (error == TP_ERRORLIMIT && options.average) {
    hl->DrawCopy("hist same l");
    hh->DrawCopy("hist same l");
  } else if (error == TP_GRAPH) {
    //		h->DrawClone(string("hist l"+string(startnew && options.newpage ? "":" same")).c_str());
    h->DrawCopy(string("hist l"+string(options.firstplot ? "":" same")).c_str());
  }
	
  h->Delete();
  hl->Delete();
  hh->Delete();
  h3->Delete();
  h3e->Delete();
}


/* DrawOffsetGraph
 *
 */
void DrawOffsetGraph(Dataset& ds, Offset *offset, string id, string title, Parameters par, TPlotType error) {

  double xmin = ds.angle.GetZaxis()->GetXmin();
  double xmax = ds.angle.GetZaxis()->GetXmax();
  double ymin = 0, ymax = 0;
  bool first = true;

  TGraphErrors* e = new TGraphErrors;
  e->SetName(id.c_str());
  int j = 0;
  for (unsigned int i=0; i<offset->GetSize(); ++i) {
    if (offset->has[i]) {
      e->SetPoint(j,      (double)i/(double)ds.time.GetNbinsZ() *ds.time.GetZaxis()->GetXmax()
		  + (1-(double)i/(double)ds.time.GetNbinsZ())*ds.time.GetZaxis()->GetXmin(), offset->val[i]);
      e->SetPointError(j, 0, sqrt(offset->err[i])/(double)offset->has[i]);
      j ++;
    }
  }

  int pad = PADSX*PADSY-1+gPad->GetNumber();
  // Determine Y axis range
  if (error == TP_AXISRANGE) {
    for (int i = 0; i < e->GetN(); i++) {
      if (first) {
	ymin = e->GetY()[i] - e->GetEY()[i];
	ymax = e->GetY()[i] + e->GetEY()[i];
	first = false;
      } else {
	if (ymin > e->GetY()[i] - e->GetEY()[i]) ymin = e->GetY()[i] - e->GetEY()[i];
	if (ymax < e->GetY()[i] + e->GetEY()[i]) ymax = e->GetY()[i] + e->GetEY()[i];
      }
    }
    if (! axisrange[pad].isinitialized) {
      axisrange[pad].max = ymax;
      axisrange[pad].min = ymin;
      axisrange[pad].isinitialized = true;
    } else {
      if (axisrange[pad].max < ymax) {
	axisrange[pad].max = ymax;
      }
      if (axisrange[pad].min > ymin) {
	axisrange[pad].min = ymin;
      }
    }
  }
  e->SetLineColor(par.linecolor);
  e->SetFillColor(par.linecolor+20);
  e->SetFillStyle(1001);
  e->SetLineStyle(par.linestyle);
  e->SetLineWidth(par.linewidth);
  e->SetMarkerSize(0);
  e->GetXaxis()->SetRangeUser(xmin, xmax);
  if (axisrange[pad].isinitialized) {
    e->GetYaxis()->SetRangeUser(axisrange[pad].min, axisrange[pad].max);
  }
  e->GetXaxis()->SetLabelFont(132);
  e->GetYaxis()->SetLabelFont(132);
  double fontsize = 0.09;
  if (gPad->GetNumber() == 1) {
    fontsize = 0.09;
  } else {
    e->GetXaxis()->SetTitleSize(fontsize);
    e->GetXaxis()->SetTitle("#font[132]{log(#font[12]{E}/GeV)}");
  }
  e->GetXaxis()->SetLabelSize(fontsize);
  e->GetYaxis()->SetLabelSize(fontsize);
  e->SetTitle(title.c_str());

  if (error == TP_ERRORAREA) {
    cout << pad << "\t" << axisrange[pad].min << "\t" << axisrange[pad].max << endl;
    e->DrawClone(string("l3"+string(options.firstplot ? "a":"a same")).c_str());
  }
  if (error == TP_ERRORLIMIT) {
    //TArrayD hl = 
    //hl->DrawCopy("hist same l");
    //hh->DrawCopy("hist same l");
  }
  if (error == TP_GRAPH) {
    // cout << pad << "\t" << axisrange[pad].min << "\t" << axisrange[pad].max << endl;
    e->DrawClone(string("lx"+string(options.firstplot ? "a":" same")).c_str());
  }
  delete e;
  //delete h;
}


/* Plot
 * **********
 * Add a dataset or dataset's error to the plot
 */
void Plot(Dataset ds, Parameters par, TPlotType pass = TP_GRAPH) {
  
  if (debug)
    cout << "Plot            : Plotting data" << endl;

  TH1 *h;
  double t_delta, l_delta, a_delta, z_delta, e_delta;
  int t_min = 0;
  int t_max = ds.time .GetNbinsX();
  int l_min = 0;
  int l_max = ds.time .GetNbinsY();
  int a_min = 0;
  int a_max = ds.angle.GetNbinsX();
  int z_min = 0;
  int z_max = ds.angle.GetNbinsY();
  int e_min = 0;
  int e_max = ds.angle.GetNbinsZ();
  t_delta = (double)t_max/(double)options.tcuts;
  l_delta = (double)l_max/(double)options.lcuts;
  a_delta = (double)a_max/(double)options.acuts;
  z_delta = (double)z_max/(double)options.zcuts;
  e_delta = (double)e_max/(double)options.ecuts;

  //TH3D *errTl, *errTh, *errAl, *errAh;

  //errTl->Multiply(&ds.time, &ds.time, ds.numberofstacks, 1.);
  //errTh->Multiply(&ds.time, &ds.time, ds.numberofstacks, 1.);

  double sumtotal = ds.angle.Integral();

  stringstream s;

  // Use style specification as histogram id; usage of letter characters is
  // discouraged according to this website:
  // http://root.cern.ch/root/html516/TH3.html#TH3:Project3D

  string id = "";
  s << par.linecolor << par.linestyle << par.linewidth << pass;
  s >> id;

  for (int i=0; i<options.ecuts; ++i) {

    ds.time .GetZaxis()->SetRange(e_min+int(i*e_delta),e_min+int((i+1)*e_delta));
    ds.angle.GetZaxis()->SetRange(e_min+int(i*e_delta),e_min+int((i+1)*e_delta));

    for (int j=0; j<options.lcuts; ++j) {
      
      stringstream s;
      string id = "";
      s << par.linecolor << par.linestyle << par.linewidth << pass << i << j;
      s >> id;
      
      ds.time.GetYaxis()->SetRange(l_min+int(j*l_delta),l_min+int((j+1)*l_delta));
      
      // Draw time lag histogram
      c1.cd(1);
      h  = ds.time.Project3D(string("x0"+id).c_str());
      
      if (debug) 
	cout << "Plot            : Adding time lag plot" << endl;
      
      DrawGraph(ds, h, "Time lag w.r.t. light speed", sumtotal, par, i+j, pass, true);

      ds.time.GetYaxis()->SetRange(l_min, l_max);
    }

    for (int k=0; k<options.tcuts; ++k) {
      
      stringstream s;
      string id = "";
      s << par.linecolor << par.linestyle << par.linewidth << pass << i << k;
      s >> id;

      ds.time.GetXaxis()->SetRange(t_min+int(k*t_delta),t_min+int((k+1)*t_delta));

      // Draw lateral density histogram
      c1.cd(2);
      h  = ds.time.Project3D(string("y0"+id).c_str());
      
      if (debug)
	cout << "Plot            : Adding lateral density plot" << endl;
      
      DrawGraph(ds, h, "Lateral density", sumtotal, par, i+k, pass, true);

      ds.time.GetYaxis()->SetRange(t_min,t_max);
    }

    for (int j = 0; j < options.zcuts; j ++) {
      stringstream s;
      string id = "";
      s << par.linecolor << par.linestyle << par.linewidth << pass << i << j;
      s >> id;

      ds.angle.GetYaxis()->SetRange(z_min+int(j*z_delta),z_min+int((j+1)*z_delta));

      // Draw radial momentum angle histogram
      c1.cd(3);
      h = ds.angle.Project3D(string("x1"+id).c_str());
      
      if (debug)
	cout << "Plot            : Adding lateral momentum angle plot" << endl;
      
      DrawGraph(ds, h, "Momentum space angle to shower axis", sumtotal, par, i+j, pass, false);

      ds.angle.GetYaxis()->SetRange(z_min,z_max);
    }

    for (int k = 0; k < options.acuts; k ++) {
      stringstream s;
      string id = "";
      s << par.linecolor << par.linestyle << par.linewidth << pass << i << k;
      s >> id;

      ds.angle.GetXaxis()->SetRange(a_min+int(k*a_delta),a_min+int((k+1)*a_delta));

      // Draw momentum angle to shower histogram
      c1.cd(4);
      h = ds.angle.Project3D(string("y1"+id).c_str());
      
      if (debug)
	cout << "Plot            : Adding radial momentum angle plot" << endl;
      
      DrawGraph(ds, h, "Momentum azimuth angle", sumtotal, par, i+k, pass, false);

      ds.angle.GetXaxis()->SetRange(a_min,a_max);
    }

    ds.time. GetZaxis()->SetRange(e_min,e_max);
    ds.angle.GetZaxis()->SetRange(e_min,e_max);
  }

  if (options.lcuts > 1 || options.tcuts > 1) {
    for (int j = 0; j < options.lcuts; j ++) {
      ds.time.GetYaxis()->SetRange(l_min+int(j*l_delta),l_min+int((j+1)*l_delta));

      for (int k = 0; k < options.tcuts; k ++) {
	stringstream s;
	string id = "";
	s << par.linecolor << par.linestyle << par.linewidth << pass << 0 << j << k;
	s >> id;

	ds.time.GetXaxis()->SetRange(t_min+int(k*t_delta),t_min+int((k+1)*t_delta));

	// Draw kinetic energy histogram
	c1.cd(5);
	h = ds.time.Project3D(string("z0"+id).c_str());
	
	if (debug)
	  cout << "Plot            : Adding kinetic energy plot" << endl;
	
	DrawGraph(ds, h, "Kinetic energy", sumtotal, par, j+k, pass, true);

	ds.time.GetXaxis()->SetRange(t_min,t_max);
      }
      ds.time.GetYaxis()->SetRange(l_min,l_max);
    }
  } else {
    for (int j = 0; j < options.zcuts; j ++) {
      ds.angle.GetYaxis()->SetRange(z_min+int(j*z_delta),z_min+int((j+1)*z_delta));

      for (int k = 0; k < options.acuts; k ++) {
	stringstream s;
	string id = "";
	s << par.linecolor << par.linestyle << par.linewidth << pass << 0 << j << k;
	s >> id;

	ds.angle.GetXaxis()->SetRange(a_min+int(k*a_delta),a_min+int((k+1)*a_delta));

	// Draw kinetic energy histogram
	c1.cd(5);
	h = ds.angle.Project3D(string("z0"+id).c_str());
	
	if (debug)
	  cout << "Plot            : Adding kinetic energy plot" << endl;
	
	DrawGraph(ds, h, "Kinetic energy", sumtotal, par, j+k, pass, true);

	ds.angle.GetXaxis()->SetRange(a_min,a_max);
      }
      ds.angle.GetYaxis()->SetRange(z_min,z_max);
    }
  }

  // Draw offset histogram
  if (ds.hasoffset) {
    c1.cd(6)->cd(1);
    DrawOffsetGraph(ds, &ds.offsetX, string("ox"+id), "Offset", par, pass);
    c1.cd(6)->cd(2);
    DrawOffsetGraph(ds, &ds.offsetY, string("oy"+id), "", par, pass);
  }

  if (pass == TP_GRAPH || options.errors && pass == TP_ERRORAREA)
    options.firstplot = false;
}


/* Textout
 * **********
 * Add a dataset or dataset's error to the 
 */
void Textout(Dataset* ds, Parameters* par, int numsets) {
  
  if (debug) 
    cout << "Plot            : Plotting data" << endl;

  TH1 *h, *he;
  double t_delta, l_delta, a_delta, z_delta, e_delta;
  int set = 0;
  int t_min = 0;
  int t_max = ds[set].time .GetNbinsX();
  int l_min = 0;
  int l_max = ds[set].time .GetNbinsY();
  int a_min = 0;
  int a_max = ds[set].angle.GetNbinsX();
  int z_min = 0;
  int z_max = ds[set].angle.GetNbinsY();
  int e_min = 0;
  int e_max = ds[set].angle.GetNbinsZ();
  t_delta = (double)t_max/(double)options.tcuts;
  l_delta = (double)l_max/(double)options.lcuts;
  a_delta = (double)a_max/(double)options.acuts;
  z_delta = (double)z_max/(double)options.zcuts;
  e_delta = (double)e_max/(double)options.ecuts;

  // -------------------------------
  cout << endl << endl << "#Time lag histogram" << endl;
  cout << "#age\tpar\tecut\tlcut\tbin\tval\terr" << endl;

  for (int set = 0; set < numsets; set ++) {
    for (int i = 0; i < options.ecuts; i ++) {
      ds[set].time    .GetZaxis()->SetRange(e_min+int(i*e_delta),e_min+int((i+1)*e_delta));
      ds[set].time_err.GetZaxis()->SetRange(e_min+int(i*e_delta),e_min+int((i+1)*e_delta));

      for (int j = 0; j < options.lcuts; j ++) {
	ds[set].time    .GetYaxis()->SetRange(l_min+int(j*l_delta),l_min+int((j+1)*l_delta));
	ds[set].time_err.GetYaxis()->SetRange(l_min+int(j*l_delta),l_min+int((j+1)*l_delta));

	h  = ds[set].time    .Project3D("x0");
	he = ds[set].time_err.Project3D("x1");
					
	for (int bin = 1; bin < h->GetNbinsX()+1; bin ++) {
	  cout << par[set].age << "\t"
	       << ((par[set].particle[0] == 'p') + 1) << "\t"
	       << i << "\t"
	       << j << "\t"
	       << bin << "\t"
	       << h ->GetBinContent(bin) << "\t"
	       << he->GetBinContent(bin)*3. << endl;
	}
	cout << endl;

	h ->Delete();
	he->Delete();
      }
      ds[set].time    .GetYaxis()->SetRange(l_min,l_max);
      ds[set].time_err.GetYaxis()->SetRange(l_min,l_max);
    }
    ds[set].time    .GetZaxis()->SetRange(e_min,e_max);
    ds[set].time_err.GetZaxis()->SetRange(e_min,e_max);
  }
	
  // -------------------------------
  cout << endl << "#Lateral histogram" << endl;
  cout << "#age\tpar\tecut\ttcut\tbin\tval\terr" << endl;
	
  for (int set = 0; set < numsets; set ++) {
    for (int i = 0; i < options.ecuts; i ++) {
      ds[set].time    .GetZaxis()->SetRange(e_min+int(i*e_delta),e_min+int((i+1)*e_delta));
      ds[set].time_err.GetZaxis()->SetRange(e_min+int(i*e_delta),e_min+int((i+1)*e_delta));

      for (int j = 0; j < options.tcuts; j ++) {
	ds[set].time    .GetXaxis()->SetRange(t_min+int(j*t_delta),t_min+int((j+1)*t_delta));
	ds[set].time_err.GetXaxis()->SetRange(t_min+int(j*t_delta),t_min+int((j+1)*t_delta));

	h  = ds[set].time    .Project3D("y0");
	he = ds[set].time_err.Project3D("y1");
				
	for (int bin = 1; bin < h->GetNbinsX()+1; bin ++) {
	  cout << par[set].age << "\t"
	       << ((par[set].particle[0] == 'p') + 1) << "\t"
	       << i << "\t"
	       << j << "\t"
	       << bin << "\t"
	       << h ->GetBinContent(bin) << "\t"
	       << he->GetBinContent(bin)*3. << endl;
	}
	cout << endl;

	h ->Delete();
	he->Delete();
      }
      ds[set].time    .GetXaxis()->SetRange(t_min,t_max);
      ds[set].time_err.GetXaxis()->SetRange(t_min,t_max);
    }
    ds[set].time    .GetZaxis()->SetRange(e_min,e_max);
    ds[set].time_err.GetZaxis()->SetRange(e_min,e_max);
  }
	
  // -------------------------------
  cout << endl << "#Vertical momentum angle" << endl;
  cout << "#age\tpar\tecut\tzcut\tbin\tval\terr" << endl;
	
  for (int set = 0; set < numsets; set ++) {
    for (int i = 0; i < options.ecuts; i ++) {
      ds[set].angle    .GetZaxis()->SetRange(e_min+int(i*e_delta),e_min+int((i+1)*e_delta));
      ds[set].angle_err.GetZaxis()->SetRange(e_min+int(i*e_delta),e_min+int((i+1)*e_delta));

      for (int j = 0; j < options.zcuts; j ++) {
	ds[set].angle    .GetYaxis()->SetRange(z_min+int(j*z_delta),z_min+int((j+1)*z_delta));
	ds[set].angle_err.GetYaxis()->SetRange(z_min+int(j*z_delta),z_min+int((j+1)*z_delta));

	h  = ds[set].angle    .Project3D("x0");
	he = ds[set].angle_err.Project3D("x1");
				
	for (int bin = 1; bin < h->GetNbinsX()+1; bin ++) {
	  cout << par[set].age << "\t"
	       << ((par[set].particle[0] == 'p') + 1) << "\t"
	       << i << "\t"
	       << j << "\t"
	       << bin << "\t"
	       << h ->GetBinContent(bin) << "\t"
	       << he->GetBinContent(bin)*3. << endl;
	}
	cout << endl;

	h ->Delete();
	he->Delete();
      }
      ds[set].angle    .GetYaxis()->SetRange(z_min,z_max);
      ds[set].angle_err.GetYaxis()->SetRange(z_min,z_max);
    }
    ds[set].angle    .GetZaxis()->SetRange(e_min,e_max);
    ds[set].angle_err.GetZaxis()->SetRange(e_min,e_max);
  }
	
  // -------------------------------
  cout << endl << "#Horizontal momentum angle" << endl;
  cout << "#age\tpar\tecut\tacut\tbin\tval\terr" << endl;
	
  for (int set = 0; set < numsets; set ++) {
    for (int i = 0; i < options.ecuts; i ++) {
      ds[set].angle    .GetZaxis()->SetRange(e_min+int(i*e_delta),e_min+int((i+1)*e_delta));
      ds[set].angle_err.GetZaxis()->SetRange(e_min+int(i*e_delta),e_min+int((i+1)*e_delta));

      for (int j = 0; j < options.acuts; j ++) {
	ds[set].angle    .GetXaxis()->SetRange(a_min+int(j*a_delta),a_min+int((j+1)*a_delta));
	ds[set].angle_err.GetXaxis()->SetRange(a_min+int(j*a_delta),a_min+int((j+1)*a_delta));

	h  = ds[set].angle    .Project3D("y0");
	he = ds[set].angle_err.Project3D("y1");
				
	for (int bin = 1; bin < h->GetNbinsX()+1; bin ++) {
	  cout << par[set].age << "\t"
	       << ((par[set].particle[0] == 'p') + 1) << "\t"
	       << i << "\t"
	       << j << "\t"
	       << bin << "\t"
	       << h ->GetBinContent(bin) << "\t"
	       << he->GetBinContent(bin)*3. << endl;
	}
	cout << endl;

	h ->Delete();
	he->Delete();
      }
      ds[set].angle    .GetXaxis()->SetRange(a_min,a_max);
      ds[set].angle_err.GetXaxis()->SetRange(a_min,a_max);
    }
    ds[set].angle    .GetZaxis()->SetRange(e_min,e_max);
    ds[set].angle_err.GetZaxis()->SetRange(e_min,e_max);
  }
	
  // -------------------------------
  cout << endl << "#Kinetic energy" << endl;
	
  if (options.lcuts > 1 || options.tcuts > 1) {
    cout << "#age\tpar\tlcut\ttcut\tbin\tval\terr" << endl;
    for (int set = 0; set < numsets; set ++) {
      for (int j = 0; j < options.lcuts; j ++) {
	ds[set].time    .GetYaxis()->SetRange(l_min+int(j*l_delta),l_min+int((j+1)*l_delta));
	ds[set].time_err.GetYaxis()->SetRange(l_min+int(j*l_delta),l_min+int((j+1)*l_delta));

	for (int k = 0; k < options.tcuts; k ++) {
	  ds[set].angle    .GetXaxis()->SetRange(t_min+int(k*t_delta),t_min+int((k+1)*t_delta));
	  ds[set].angle_err.GetXaxis()->SetRange(t_min+int(k*t_delta),t_min+int((k+1)*t_delta));

	  h  = ds[set].time    .Project3D("z0");
	  he = ds[set].time_err.Project3D("z1");
				
	  for (int bin = 1; bin < h->GetNbinsX()+1; bin ++) {
	    cout << par[set].age << "\t"
		 << ((par[set].particle[0] == 'p') + 1) << "\t"
		 << j << "\t"
		 << k << "\t"
		 << bin << "\t"
		 << h ->GetBinContent(bin) << "\t"
		 << he->GetBinContent(bin)*3. << endl;
	  }
	  cout << endl;

	  h ->Delete();
	  he->Delete();
	}
	ds[set].time    .GetXaxis()->SetRange(t_min,t_max);
	ds[set].time_err.GetXaxis()->SetRange(t_min,t_max);
      }
      ds[set].time    .GetYaxis()->SetRange(l_min,l_max);
      ds[set].time_err.GetYaxis()->SetRange(l_min,l_max);
    }
  } else {
    cout << "#age\tpar\tzcut\tacut\tbin\tval\terr" << endl;
    for (int set = 0; set < numsets; set ++) {
      for (int j = 0; j < options.zcuts; j ++) {
	ds[set].angle    .GetYaxis()->SetRange(z_min+int(j*z_delta),z_min+int((j+1)*z_delta));
	ds[set].angle_err.GetYaxis()->SetRange(z_min+int(j*z_delta),z_min+int((j+1)*z_delta));

	for (int k = 0; k < options.acuts; k ++) {
	  ds[set].angle    .GetXaxis()->SetRange(a_min+int(k*a_delta),a_min+int((k+1)*a_delta));
	  ds[set].angle_err.GetXaxis()->SetRange(a_min+int(k*a_delta),a_min+int((k+1)*a_delta));

	  h  = ds[set].angle    .Project3D("z0");
	  he = ds[set].angle_err.Project3D("z1");
				
	  for (int bin = 1; bin < h->GetNbinsX()+1; bin ++) {
	    cout << par[set].age << "\t"
		 << ((par[set].particle[0] == 'p') + 1) << "\t"
		 << j << "\t"
		 << k << "\t"
		 << bin << "\t"
		 << h ->GetBinContent(bin) << "\t"
		 << he->GetBinContent(bin)*3. << endl;
	  }
	  cout << endl;

	  h ->Delete();
	  he->Delete();
	}
	ds[set].angle    .GetXaxis()->SetRange(a_min,a_max);
	ds[set].angle_err.GetXaxis()->SetRange(a_min,a_max);
      }
      ds[set].angle    .GetYaxis()->SetRange(z_min,z_max);
      ds[set].angle_err.GetYaxis()->SetRange(z_min,z_max);
    }
  }

  // -------------------------------
}

/*
 * Excess
 * ******
 *
 * calculated the charge excess
 */ 
Dataset Excess(Dataset& e, Dataset& p) {
  Dataset s = e;
  Dataset q = e;

  s.Add(p, true);
  q.Add(p, true, -1.);
	
  q.time .Divide(&s.time);
  q.angle.Divide(&s.angle);

  q.isallocated = true;
  q.hasoffset   = false;
  q.hasenergy   = false;
  q.haserror    = false;
  q.numberofstacks = 1;
  q.norm = 1.;

  return q;
}

/* Export
 * ******
 * Export the current plot canvas to file
 */
void Export() {
  
  string ps;
  if (options.postscript == 2) {
    ps = options.outfile;
    options.outfile.replace(options.outfile.length() -3, 3, ".pdf", 4);
  }
  
  if (debug) 
    cout << "Export          : Exporting to file " << options.outfile << endl;

  if (options.postscript) {
    c1.Update();
    options.psstream->Close();
    if (options.postscript == 2) {
      system(string("ps2pdfwr '"+ps+"' '"+options.outfile+"'").c_str());
      unlink(ps.c_str());
    }
  } else {
    c1.Print(options.outfile.c_str());
  }

  if (debug) cout << "Export          : Complete" << endl;
  else       cout << "Output written to " << options.outfile << endl;
}


/// -------------------------------------------------------------------
/// -------------------------------------------------------------------
/// -------------------------------------------------------------------

/* Main
 * ****
 * Main routine
 */
int main(int argc, const char** argv) {
  
  if (debug)
    cout << "rootplot        :" << endl;
  
  int err = 0;
  if ( ! (err = ProcessOptions(argc, argv)) ) {
    return err;
  }

  // set plotting style and prepare canvas
  InitPlot();
  
  if (debug)
    cout << "rootplot        : Averaging plots" << endl;
  
  Dataset avData[plots.size()];
  Parameters avParameters[plots.size()];
  int max = 0;
  
  // Read data from files
  int iPlot = 0;
  for (vector<Parameters>::iterator s = plots.begin();
       s != plots.end(); ++s) {
    
    Parameters pcp = *s;
    pcp.file = "";
    Dataset current = GetData(*s);
    if (pcp.particle[0] == 's') {
      pcp.particle[0] = 'p';
      Dataset extra = GetData(*s);
      pcp.particle[0] = 's';
      current.Add(extra, true);
    } else if (pcp.particle[0] == 'd') {
      pcp.particle[0] = 'p';
      Dataset extra = GetData(*s);
      pcp.particle[0] = 'd';
      current = Excess(current, extra);
    }

    int _av = 0;
    if (!options.average) _av = max;
    while (_av<max && avParameters[_av]!=pcp) 
      _av ++;
    
    if (debug)
      cout << "rootplot        : Plot " << iPlot << " to average " << _av << endl;

    if (_av == max) {
      avData[_av].Add(current);
      avParameters[_av] = pcp;
      max ++;
    } else {
      avData[_av].Add(current);
    }
    
    ++iPlot;
    
  } // end loop plots
  
  
  // Determine axis range for all plots
  int i = 0;
  while (i < max) {
    if (options.errors) 
      CalculateErrors(avData[i]);
    
    if (!options.textout) 
      Plot(avData[i], avParameters[i], TP_AXISRANGE);
    
    i++;
  }

  // Plot everything
  if (options.textout) {
    
    Textout(avData, avParameters, max);
    
  } else if (options.newpage) {
    
    i = 0;
    while (i < max) {
      int j = i;
      if (options.errors) {
	do {
	  Plot(avData[j], avParameters[j], TP_ERRORAREA);
	  j ++;
	} while (j < max && (j == 0 || avParameters[j].OnlyDiffersByParticle(avParameters[j-1])));
	j = i;
	do {
	  Plot(avData[j], avParameters[j], TP_ERRORLIMIT);
	  j ++;
	} while (j < max && (j == 0 || avParameters[j].OnlyDiffersByParticle(avParameters[j-1])));
	j = i;
      }
      do {
	Plot(avData[j], avParameters[j], TP_GRAPH);
	j ++;
      } while (j < max && (j == 0 || avParameters[j].OnlyDiffersByParticle(avParameters[j-1])));

      i = j;
      if (i < max-1) {
	c1.Update();
	if (options.postscript) options.psstream->NewPage();
	options.firstplot = true;
      }
    }
    Export();

  } else {

    // Make sure all errors are drawn BELOW actual plots
    // to keep it readible
    if (options.errors) {
      for (int i = 0; i < max; i ++) {
	Plot(avData[i], avParameters[i], TP_ERRORAREA);
      }
      for (int i = 0; i < max; i ++) {
	Plot(avData[i], avParameters[i], TP_ERRORLIMIT);
      }
    }
    for (int i = 0; i < max; i ++) {
      Plot(avData[i], avParameters[i]);
    }
    Export();
  }
	
  if (debug) 
    cout << "rootplot        : Exiting" << endl;

  return 0;
}

