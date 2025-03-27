#include "TCanvas.h"

// canvas
TCanvas* MyCanvas(const char* name, Int_t width = 650, Int_t height = 600) {
	auto* c = new TCanvas(name, "", width, height);
	c->SetFillColor(0);
	return c;
}

void SaveMyCanvas(TCanvas* canv, const char* name) {
	canv->Update();
	canv->RedrawAxis();
	gSystem->mkdir("Figures/", kTRUE);
	canv->SaveAs(Form("Figures/%s.pdf", name), "RECREATE");
	canv->SaveAs(Form("Figures/%s.png", name), "RECREATE");
	//canv->Close();
}

// graphs
