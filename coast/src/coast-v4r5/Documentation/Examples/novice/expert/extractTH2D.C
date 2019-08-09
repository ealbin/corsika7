#include <cmath>
using namespace std;

TH2D *extractX(TH3D *h3,const char *name,int ix,int ix_up=0)
{
    if (ix_up==0) ix_up=ix;

    int xbins=h3->GetXaxis()->GetNbins();
    int ybins=h3->GetYaxis()->GetNbins();
    int zbins=h3->GetZaxis()->GetNbins();

    double xmin=h3->GetXaxis()->GetXmin();
    double ymin=h3->GetYaxis()->GetXmin();
    double zmin=h3->GetZaxis()->GetXmin();

    double xmax=h3->GetXaxis()->GetXmax();
    double ymax=h3->GetYaxis()->GetXmax();
    double zmax=h3->GetZaxis()->GetXmax();


    printf("extraxtX from ix=%d up to ix_up=%d\n",ix,ix_up);
    printf("xbins=%4d \t ybins=%4d \t zbins=%4d\n",xbins,ybins,zbins);
    printf("xmin =%4.0f \t ymin =%4.0f \t zmin =%4.0f\n",xmin,ymin,zmin);
    printf("xmax =%4.0f \t ymax =%4.0f \t zmax =%4.0f\n",xmax,ymax,zmax);
    
   
    

    if (ix>xbins || ix<1) return NULL;

    TH2D *h2=new TH2D(name,name,ybins,ymin,ymax,zbins,zmin,zmax);
	for (int iy=1;iy<=ybins;iy++)
	    for (int iz=1;iz<=zbins;iz++)
	    {
		double sum=0;
		double sumerror2=0;
		for (int isum=ix;isum<=ix_up;isum++) {
		    sum+=h3->GetBinContent(isum,iy,iz);
		    sumerror2+=h3->GetBinError(isum,iy,iz)*h3->GetBinError(isum,iy,iz);
		
		}
		
		{
		h2->SetBinContent(iy,iz,sum);
		h2->SetBinError(iy,iz,sqrt(sumerror2));
		}



	    }



    return h2;
}




TH2D *extractY(const TH3D *h3,const char *name,int iy,int iy_up=0)
{  
    if (iy_up==0) iy_up=iy;

    int xbins=h3->GetXaxis()->GetNbins();
    int ybins=h3->GetYaxis()->GetNbins();
    int zbins=h3->GetZaxis()->GetNbins();

    double xmin=h3->GetXaxis()->GetXmin();
    double ymin=h3->GetYaxis()->GetXmin();
    double zmin=h3->GetZaxis()->GetXmin();

    double xmax=h3->GetXaxis()->GetXmax();
    double ymax=h3->GetYaxis()->GetXmax();
    double zmax=h3->GetZaxis()->GetXmax();

    if (iy>ybins+1 || iy<0) {printf("%d %d\n",iy,ybins);return NULL;}

    TH2D *h2=new TH2D(name,name,xbins,xmin,xmax,zbins,zmin,zmax);
    //TH2D *h2 = (TH2D*) h3->Clone ("name");
    //h2->Reset ("ICE");

    //for (int ix=1; ix<=ybins; ix++)
    for (int ix=0; ix<=xbins+1; ix++) {
	
	for (int iz=0; iz<=zbins+1; iz++)	{

	    double sum=0;
	    double sumerror2=0;

	    for (int isum=iy; isum<=iy_up; isum++) {
		
		double bin = h3->GetBinContent (ix,isum,iz);
		double err = h3->GetBinError (ix,isum,iz);
		double errN = sqrt (err*err + pow (h2->GetBinError(ix,iz),2));

		if (bin!=0 || err!=0) {
		    h2->Fill (h3->GetXaxis ()->GetBinCenter (ix), 
			      h3->GetZaxis ()->GetBinCenter (iz), bin);
		    h2->SetBinError (ix, iz, errN);
		}
		
	    }
	    
	    /*
	    if (sum!=0) {
		h2->SetBinContent(ix,iz,sum);
		h2->SetBinError(ix,iz,sqrt(sumerror2));
		}
	    */
	}
    }

    return h2;
}





#include <TGraph2D.h>
#include <TGraph2DErrors.h>

//TGraph2DErrors *TH2DToGraph2D(TH2D *h2,float wcut)
TGraph2D *TH2DToGraph2D(TH2D *h2,float wcut)
{


    int xbins=h2->GetXaxis()->GetNbins();
    int ybins=h2->GetYaxis()->GetNbins();
   

    double xmin=h2->GetXaxis()->GetXmin();
    double ymin=h2->GetYaxis()->GetXmin();
   

    double xmax=h2->GetXaxis()->GetXmax();
    double ymax=h2->GetYaxis()->GetXmax();
      

    double xwidth=h2->GetXaxis()->GetBinWidth(1);
    double ywidth=h2->GetYaxis()->GetBinWidth(1);
   


    int nbins=0;
	for (int iy=1;iy<=ybins;iy++)
	    for (int ix=1;ix<=xbins;ix++)
	    {
		double w=h2->GetBinContent(ix,iy);

		if (w>wcut) nbins++;
		
	    }

	int ibin=0;
//	TGraph2DErrors *gr=new TGraph2DErrors(nbins);
	TGraph2D *gr=new TGraph2D(nbins);

	for (int iy=1;iy<=ybins;iy++)
	    for (int ix=1;ix<=xbins;ix++)
	    {
		double w=h2->GetBinContent(ix,iy);

		if (w>wcut) 
		{
		double ew=1;sqrt(w);
		
		gr->SetPoint(ibin,xmin+(ix-0.5)*xwidth,ymin+(iy-0.5)*ywidth,w);
//		gr->SetPointError(ibin,0.5*xwidth,0.5*ywidth,ew);
//	gr->SetPointError(ibin,0,0,ew);
		ibin++;
		}

	    }




    return gr;
}


void TH2DCut(TH2D *h2,float wcut,int flag=1)
{


    int xbins=h2->GetXaxis()->GetNbins();
    int ybins=h2->GetYaxis()->GetNbins();
   
    int nbins=0;
	for (int iy=1;iy<=ybins;iy++)
	    for (int ix=1;ix<=xbins;ix++)
	    {
		double w=h2->GetBinContent(ix,iy);

		if (flag==1 && w<wcut)
		    { 
			h2->SetBinContent(ix,iy,0);
			h2->SetBinError(ix,iy,0);
		    }
		if (flag==-1 && w>wcut)
		    { 
			h2->SetBinContent(ix,iy,0);
			h2->SetBinError(ix,iy,0);
		    }
	    }

}


TH2D *TH2DGetChi(const char *name,TH2D *h2e,TH2D *h2)
{

    

    int xbins=h2->GetXaxis()->GetNbins();
    int ybins=h2->GetYaxis()->GetNbins();
   
    double xmin=h2->GetXaxis()->GetXmin();
    double ymin=h2->GetYaxis()->GetXmin();
   
    double xmax=h2->GetXaxis()->GetXmax();
    double ymax=h2->GetYaxis()->GetXmax();


    TH2D *hChi=new TH2D(name,name,xbins,xmin,xmax,ybins,ymin,ymax);
   
    int nbins=0;
	for (int iy=1;iy<=ybins;iy++)
	    for (int ix=1;ix<=xbins;ix++)
	    {
		double wh2=h2->GetBinContent(ix,iy);
		double wh2e=h2e->GetBinContent(ix,iy);
		double errorh2e=h2e->GetBinError(ix,iy);

		

		if (errorh2e<1) errorh2e=1;
		hChi->SetBinContent(ix,iy,pow((wh2-wh2e)/errorh2e,2));

		
	    }
	return hChi;
}
