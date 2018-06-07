#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"

#include "TFile.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TGraph2D.h"
#include "TMath.h"
//#include "lardataobj/RecoBase/SpacePoint.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TString.h"
#include "nusimdata/SimulationBase/MCParticle.h"

using namespace art;
using namespace std;

// make a poor-man's event display of mc -- for the ievcount'th event in the file
// Kirsty Duffy, Fermilab, Jan 16, 2017. Based on macro by
// Tom Junk, Fermilab, Oct 5, 2017.  Based on gallery demos by Marc Paterno.
//
// Usage:  setup dunetpc (or uboonecode) and gallery, run root.
//  .L centmc.C++
// centmc("inputfile.root","outputfile.root");
// the clip region is defined with mins and maxes.  A margin can be added around the outside.
// center the display on the first mc particle trajectory point


// some example default arguments
//void centmc(size_t ievcount=0, std::string const& filename="genie126nue_gen_g4_garfield_100.root", bool showneutrons=false,
//void centmc(size_t ievcount=0, std::string const& filename="g4garfield_numu100.root", bool showneutrons=false,
//       std::string tagstring="largeant",
//       double xclipmin=-250, double xclipmax=250, double yclipmin=-250, double yclipmax=250, double zclipmin=-250, double zclipmax=250,
//       double margin=0)
//void
//centmc(size_t ievcount=0, std::string const& filename="g4lar100nue.root", bool showneutrons=false,

void mctruth_evtdisplay(std::string const& filename, int ievcount=-999, double xcentre=-999, double ycentre=-999, double zcentre=-999, TString outdir = "./", std::string const& outfilename="mctruthevtdisplay_out.root", bool showneutrons=false,
       std::string tagstring="largeant",
       double xclipmin=-15, double xclipmax=15, double yclipmin=-15, double yclipmax=15, double zclipmin=-15, double zclipmax=15,
       double margin=0)
{

  InputTag mcptag{ tagstring };
  // Create a vector of length 1, containing the given filename.
  vector<string> filenames(1, filename);

  //TFile *fout = new TFile(outfilename.c_str(),"RECREATE");

  for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()) {

    auto const& mcparticles = *ev.getValidHandle<vector<simb::MCParticle>>(mcptag);
    if (!mcparticles.empty())
      {
	size_t icolor=1;

	if (ievcount != -999 && ievcount != ev.eventAuxiliary().event()) continue;


	TString savename = Form("evt_%d",ev.eventAuxiliary().event());
	TCanvas *c = new TCanvas(savename.Data(),"TGraph2D Event Display",0,0,800,800);

	// Make a dummy graph to define the axis scales and titles, and the plot title.

	TGraph2D *grt = new TGraph2D();
	grt->SetPoint(0,zclipmin-margin,xclipmin-margin,yclipmin-margin);
	grt->SetPoint(1,zclipmax+margin,xclipmax+margin,yclipmax+margin);
	grt->SetMarkerColor(0);
	grt->SetLineColor(0);
	grt->SetMarkerStyle(1);
	std::string titlestring=tagstring;
	titlestring += " MCParticle Display";
	grt->SetTitle(titlestring.c_str());
	grt->Draw("P");
	grt->GetXaxis()->SetTitle("Z");
	grt->GetXaxis()->SetTitleColor(4);
	grt->GetYaxis()->SetTitle("X");
	grt->GetYaxis()->SetTitleColor(4);
	grt->GetZaxis()->SetTitle("Y");
	grt->GetZaxis()->SetTitleColor(4);

	TGraph2D *gmu = new TGraph2D();
	gmu->SetLineColor(1);
	TGraph2D *gpi = new TGraph2D();
	gpi->SetLineColor(2);
	TGraph2D *gp = new TGraph2D();
	gp->SetLineColor(4);
	TGraph2D *ge = new TGraph2D();
	ge->SetLineColor(6);
	TLegend *l = new TLegend(0.75,0.75,0.88,0.88);
	l->AddEntry(gmu,"True #mu","l");
	l->AddEntry(gpi,"True #pi","l");
	l->AddEntry(gp,"True p","l");
	l->AddEntry(ge,"True e","l");
	l->Draw();

	TLorentzVector centerpos;
	bool first = false; // if true, centre display on first MCParticle trajectory point
	if (xcentre==-999 && ycentre==-999 && zcentre==-999){
	  first = true;
	}
	else{
	  centerpos.SetX(xcentre);
	  centerpos.SetY(ycentre);
	  centerpos.SetZ(zcentre);
	}

	for (size_t imcp=0;imcp<mcparticles.size(); ++imcp)
	  {
	    int ipdg = std::abs(mcparticles[imcp].PdgCode());
	    // plot only leptons, p, k, pi
	    int icolor=0;
	    if (ipdg == 11) icolor = 6;  // magenta for electrons
	    if (ipdg == 13) icolor = 1;  // black for muons
	    if (ipdg == 15) icolor = 5;  // taus
	    if (ipdg == 211) icolor = 2;  // red for pions
	    if (ipdg == 321) icolor = 7; // cyan for kaons
	    if (ipdg == 2212) icolor = 4; // protons are blue
	    if (ipdg == 2112 && showneutrons) icolor = 3; // neutrons are green

	    if (icolor>0)
	      {
		size_t numtpts = mcparticles[imcp].NumberTrajectoryPoints();

		// see if this particles has any points in our clip region
		bool isinclip = false;
		for (size_t ipt=0; ipt<numtpts; ++ipt)
		  {
		    TLorentzVector apos = mcparticles[imcp].Position(ipt);  // absolute position
		    if (first)
		      {
			first = false;
			centerpos = apos;
			cout << "found center pos: " << centerpos.X() << " " << centerpos.Y() << " " << centerpos.Z() << endl;
		      }
		    TLorentzVector pos = apos - centerpos;  // relative to the center position
		    double x = pos.X();
		    double y = pos.Y();
		    double z = pos.Z();
		    if (x<xclipmin || x>xclipmax) continue;
		    if (y<yclipmin || y>yclipmax) continue;
		    if (z<zclipmin || z>zclipmax) continue;
		    isinclip = true;
		  }

		if (isinclip)
		  {
		    bool prevpointinclip = false;
		    double xprev=0;
		    double yprev=0;
		    double zprev=0;

		    TGraph2D *gri = new TGraph2D();
		    int ipc=0;
		    for (size_t ipt=0; ipt<numtpts; ++ipt)
		      {
			TLorentzVector apos = mcparticles[imcp].Position(ipt);  // absolute position
			TLorentzVector pos = apos - centerpos;  // relative to the center position
			double x = pos.X();
			double y = pos.Y();
			double z = pos.Z();

			if (x>xclipmin && x<xclipmax &&
			    y>yclipmin && y<yclipmax &&
			    z>zclipmin && z<zclipmax)
			  {
			    gri->SetPoint(ipc,z,x,y);
			    ++ipc;
			    prevpointinclip = true;
			    xprev = x;
			    yprev = y;
			    zprev = z;
			  }
			else
			  {
			    if (prevpointinclip)  // current point not in clip but previous point is.  Still need to draw
			      // a line to the edge of the clip region in the direction of the current point outside.
			      {
				double xcos = x-xprev;
				double ycos = y-yprev;
				double zcos = z-zprev;
				double xedge = 0;   // new point
				double yedge = 0;
				double zedge = 0;
				double xec=0;   // candidate new point
				double yec=0;
				double zec=0;
				double cmm = sqrt(xcos*xcos+ycos*ycos+zcos*zcos);
				if (cmm>0)  // shouldn't happen -- two points on top of one another should not be on different sides of the clip
				  {
				    xcos /= cmm;
				    ycos /= cmm;
				    zcos /= cmm;

				    if (x>xclipmax)
				      {
					if (xcos != 0)
					  {
					    xec = xprev + xcos*(xclipmax-xprev)/xcos;  // should just be xclipmax
					    yec = yprev + ycos*(xclipmax-xprev)/xcos;
					    zec = zprev + zcos*(xclipmax-xprev)/xcos;
					    if (yec>yclipmin && yec<yclipmax && zec>zclipmin && zec<zclipmax)
					      {
						xedge = xec;
						yedge = yec;
						zedge = zec;
					      }
					  }
				      }
				    if (x<xclipmin)
				      {
					if (xcos != 0)
					  {
					    xec = xprev + xcos*(xclipmin-xprev)/xcos;  // should just be xclipmin
					    yec = yprev + ycos*(xclipmin-xprev)/xcos;
					    zec = zprev + zcos*(xclipmin-xprev)/xcos;
					    if (yec>yclipmin && yec<yclipmax && zec>zclipmin && zec<zclipmax)
					      {
						xedge = xec;
						yedge = yec;
						zedge = zec;
					      }
					  }
				      }

				    if (y>yclipmax)
				      {
					if (ycos != 0)
					  {
					    xec = xprev + xcos*(yclipmax-yprev)/ycos;  // should just be yclipmax
					    yec = yprev + ycos*(yclipmax-yprev)/ycos;
					    zec = zprev + zcos*(yclipmax-yprev)/ycos;
					    if (xec>xclipmin && xec<xclipmax && zec>zclipmin && zec<zclipmax)
					      {
						xedge = xec;
						yedge = yec;
						zedge = zec;
						  }
					  }
				      }
				    if (y<yclipmin)
				      {
					if (ycos != 0)
					  {
					    xec = xprev + xcos*(yclipmin-yprev)/ycos;  // should just be yclipmin
					    yec = yprev + ycos*(yclipmin-yprev)/ycos;
					    zec = zprev + zcos*(yclipmin-yprev)/ycos;
					    if (xec>xclipmin && xec<xclipmax && zec>zclipmin && zec<zclipmax)
					      {
						xedge = xec;
						yedge = yec;
						zedge = zec;
					      }
					  }
				      }

				    if (z>zclipmax)
				      {
					if (zcos != 0)
					  {
					    xec = xprev + xcos*(zclipmax-zprev)/zcos;  // should just be yclipmax
					    yec = yprev + ycos*(zclipmax-zprev)/zcos;
					    zec = zprev + zcos*(zclipmax-zprev)/zcos;
					    if (xec>xclipmin && xec<xclipmax && yec>yclipmin && yec<yclipmax)
					      {
						xedge = xec;
						yedge = yec;
						zedge = zec;
					      }
					  }
				      }
				    if (z<zclipmin)
				      {
					if (zcos != 0)
					  {
					    xec = xprev + xcos*(zclipmin-zprev)/zcos;  // should just be yclipmin
					    yec = yprev + ycos*(zclipmin-zprev)/zcos;
					    zec = zprev + zcos*(zclipmin-zprev)/zcos;
					    if (xec>xclipmin && xec<xclipmax && yec>yclipmin && yec<yclipmax)
					      {
						xedge = xec;
						yedge = yec;
						zedge = zec;
					      }
					  }
				      }

				    if (xedge != 0 || yedge != 0 || zedge != 0)
				      {

					gri->SetPoint(ipc,zedge,xedge,yedge);
					gri->SetMarkerColor(icolor);
					gri->SetLineColor(icolor);
					gri->SetMarkerStyle(1);
					gri->Draw("LINE,SAME");
					gri = new TGraph2D();
					ipc = 0;
				      }
				  }  // test to see if cmm == 0
			      }  // end test to see if prevpointinclip
			    prevpointinclip = false;
			  } // end block of current point not in clip

			    // an attempt to limit the sizes of the TGraphs.  Turns out it was just a ROOT drawing bug.

			    //if (ipc>400)
			    //  {
			    //    gri->SetMarkerColor(icolor);
			    //    gri->SetLineColor(icolor);
			    //    gri->SetMarkerStyle(1);
			    //    gri->Draw("LINE,SAME");
			    //    gri = new TGraph2D();
			    //    ipc = 0;
			    //    gri->SetPoint(ipc,z,x,y);
			    //    ++ipc;
			    // }
			    //cout << ipt << " " << x << " " << y << " " << z << endl;
		      }
		    gri->SetMarkerColor(icolor);
		    gri->SetLineColor(icolor);
		    gri->SetMarkerStyle(1);
		    if (gri->GetN()) gri->Draw("LINE,SAME");
		  }
	      }
	  }
	//c->Update();
	//fout->cd();
	//c->Write(savename.Data());

	TString printname = outdir+Form("/mctruth_run_%d_evt_%d.pdf",ev.eventAuxiliary().run(),ev.eventAuxiliary().event());
	c->Print(printname.Data());
      }

  }
}
