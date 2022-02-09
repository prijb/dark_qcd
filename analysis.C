#include <vector>



/*Code to compare the signal (tree) and QCD background (tree1) with an option to enable the trigger on either the signal or the background
  */



//Create functions for delta phi and delta eta
double del_phi(Double_t phi1, Double_t phi2) {
	Double_t pi = TMath::Pi();
	Double_t d_phi = phi1 - phi2;
	while (d_phi >= pi) d_phi -= 2 * pi;
	while (d_phi < -pi) d_phi += 2 * pi;
	return d_phi;
}
double del_eta(Double_t eta1, Double_t eta2) {
	Double_t d_eta = eta1 - eta2;
	return d_eta;
}

//Create a function to make a histogram
TH1D* make_hist(const char* hist_name, const char* hist_title, Int_t nbins, Double_t bin_start, Double_t bin_end, const char* xaxis, const char* yaxis, Color_t lcolor, Width_t lwidth) {
	TH1D* out_hist = new TH1D(hist_name, hist_title, nbins, bin_start, bin_end);
	out_hist->SetLineColor(lcolor);
	out_hist->SetLineWidth(lwidth);
	out_hist->GetXaxis()->SetTitle(xaxis);
	out_hist->GetYaxis()->SetTitle(yaxis);
	return out_hist;

}

//Create a function to draw stacked histograms
void hist_stack(const char* stack_name, std::vector<TH1D*> hist_vector) {

	if (hist_vector.size() == 2) {
		TCanvas* c1 = new TCanvas();
		const char* stack_title = hist_vector.at(0)->GetTitle();
		THStack* out_stack = new THStack(stack_name, stack_title);
		out_stack->Add(hist_vector.at(0));
		out_stack->Add(hist_vector.at(1),"S");
		out_stack->Draw("hist nostack");
		out_stack->GetXaxis()->SetTitle(hist_vector.at(0)->GetXaxis()->GetTitle());
		out_stack->GetYaxis()->SetTitle(hist_vector.at(0)->GetYaxis()->GetTitle());
		TLegend* out_legend = new TLegend(0.78, 0.695, 0.98, 0.775);
		out_legend->AddEntry(hist_vector.at(0)->GetName(), "Background", "l");
		out_legend->AddEntry(hist_vector.at(1)->GetName(), "Signal", "l");
		out_legend->Draw("Same");
		gPad->SetLogy();
	}
	else if (hist_vector.size() == 4) {
		TCanvas* c1 = new TCanvas();
		const char* stack_title = hist_vector.at(0)->GetTitle();
		THStack* out_stack = new THStack(stack_name, stack_title);
		out_stack->Add(hist_vector.at(0));
		out_stack->Add(hist_vector.at(1));
		out_stack->Add(hist_vector.at(2));
		out_stack->Add(hist_vector.at(3));
		out_stack->Draw("hist nostack");
		out_stack->GetXaxis()->SetTitle(hist_vector.at(0)->GetXaxis()->GetTitle());
		out_stack->GetYaxis()->SetTitle(hist_vector.at(0)->GetYaxis()->GetTitle());
		TLegend* out_legend = new TLegend(0.78, 0.695, 0.98, 0.775);
		out_legend->AddEntry(hist_vector.at(0)->GetName(), "Background leading", "l");
		out_legend->AddEntry(hist_vector.at(1)->GetName(), "Background subleading", "l");
		out_legend->AddEntry(hist_vector.at(2)->GetName(), "Signal leading", "l");
		out_legend->AddEntry(hist_vector.at(3)->GetName(), "Signal subleading", "l");
		out_legend->Draw("Same");
		gPad->SetLogy();
	}
}

//Create a class for a jet array which stores variables (pt_vec,...,index_vec) for n jets as n sized vectors 
class Jet_array {
	public:
		std::vector<double>* pt_vec = new std::vector<double>;
		std::vector<double>* e_vec = new std::vector<double>;
		std::vector<double>* echad_vec = new std::vector<double>;
		std::vector<double>* ehad_vec = new std::vector<double>;
		std::vector<double>* ecem_vec = new std::vector<double>;
		std::vector<double>* eem_vec = new std::vector<double>;
		std::vector<double>* phi_vec = new std::vector<double>;
		std::vector<double>* eta_vec = new std::vector<double>;
		std::vector<double>* index_vec = new std::vector<double>;
		Int_t max_ent;                                          //max_ent is the index of the jet with the largest transverse momentum

		~Jet_array() {
			delete pt_vec;
			delete e_vec;
			delete echad_vec;
			delete ehad_vec;
			delete ecem_vec;
			delete eem_vec;
			delete phi_vec;
			delete eta_vec;
			delete index_vec;
		}
};

//Create a class for a lightparticle array which stores variables (pt_vec,...,index_vec) for n generator particles as n sized vectors 
class Light_array {
	public:
		std::vector<double>* px_vec = new std::vector<double>;
		std::vector<double>* py_vec = new std::vector<double>;
		std::vector<double>* pt_vec = new std::vector<double>;
		std::vector<double>* eta_vec = new std::vector<double>;
		std::vector<double>* phi_vec = new std::vector<double>;
		std::vector<double>* pdgId_vec = new std::vector<double>;

		~Light_array(){
			delete px_vec;
                	delete py_vec;
                	delete pt_vec;
                	delete eta_vec;
                	delete phi_vec;
			delete pdgId_vec;
		}
};


//Create a function for jet clustering
std::vector<Light_array*>* cluster(Light_array* light_in, Jet_array* jet_in) {

	/*Define an output vector of two light arrays light_jet and light_rem.
	light_jet captures the particles clustered to reference jet (jet_in)
	light_rem captures the unclustered particles*/
	std::vector<Light_array*>* light_out = new std::vector<Light_array*>;
	Int_t nEnt = light_in->eta_vec->size();
	Light_array* light_jet = new Light_array;
	Light_array* light_rem = new Light_array;

	/*Perform the actual clustering: Takes jet_in and obtains deta and dphi with respect to the i'th particle in light_in.
	 If sqrt(deta**2+dphi**2) < 0.4, the i'th particle in light_in is considered to be clustered to jet_in and is placed in light_jet */
	for (Int_t i = 0; i < nEnt; i++) {
		Double_t dph = del_phi(jet_in->phi_vec->at(jet_in->max_ent), light_in->phi_vec->at(i));
		Double_t det = del_eta(jet_in->eta_vec->at(jet_in->max_ent), light_in->eta_vec->at(i));
		Double_t dr = TMath::Sqrt(dph * dph + det * det);
		if (dr < 0.4) {
			light_jet->px_vec->push_back(light_in->px_vec->at(i));
			light_jet->py_vec->push_back(light_in->py_vec->at(i));
			light_jet->pt_vec->push_back(light_in->pt_vec->at(i));
			light_jet->phi_vec->push_back(light_in->phi_vec->at(i));
			light_jet->eta_vec->push_back(light_in->eta_vec->at(i));
			light_jet->pdgId_vec->push_back(light_in->pdgId_vec->at(i));
		}
		else {
			light_rem->px_vec->push_back(light_in->px_vec->at(i));
			light_rem->py_vec->push_back(light_in->py_vec->at(i));
			light_rem->pt_vec->push_back(light_in->pt_vec->at(i));
			light_rem->phi_vec->push_back(light_in->phi_vec->at(i));
			light_rem->eta_vec->push_back(light_in->eta_vec->at(i));
			light_rem->pdgId_vec->push_back(light_in->pdgId_vec->at(i));
		}
	}
	light_out->push_back(light_jet);
	light_out->push_back(light_rem);
	return light_out;
}

//Create a function for removing maximum from a jet
Jet_array* jet_reduce(Jet_array* jet_in) {

	//Define the output jet_array
	Jet_array* jet_out = new Jet_array;
	Int_t nEnt_jet = jet_in->eta_vec->size();
	
	//Fill output Jet_array with elements that aren't the maximum
	for (Int_t i = 0; i < nEnt_jet; i++) {
		if (i != jet_in->max_ent) {
			jet_out->pt_vec->push_back(jet_in->pt_vec->at(i));
			jet_out->e_vec->push_back(jet_in->e_vec->at(i));
			jet_out->echad_vec->push_back(jet_in->echad_vec->at(i));
			jet_out->ehad_vec->push_back(jet_in->ehad_vec->at(i));
			jet_out->ecem_vec->push_back(jet_in->ecem_vec->at(i));
			jet_out->eem_vec->push_back(jet_in->eem_vec->at(i));
			jet_out->phi_vec->push_back(jet_in->phi_vec->at(i));
			jet_out->eta_vec->push_back(jet_in->eta_vec->at(i));
			jet_out->index_vec->push_back(jet_in->index_vec->at(i));
		}
	}
	return jet_out;
}

void analysis() {
	
	int sig_trigger = 0;
	int back_trigger = 0;

	//Part of code that considers the signal
	//Code to change file parameters
        std::string portal = "vector";
        std::string mass = "5";
        std::string ctau = "10";
        std::string xi = "1";
        std::string xl = "1";

        std::string fname_int = "../Ntuples/HiddenValley_"+portal+"_m_"+mass+"_ctau_"+ctau+"_xi_"+xi+"_xl_"+xl+"_privateMC_102X_DARKTIME_v5_generationNoFilter_forBP_2018_miniTrees.root";
        const char* fname = fname_int.c_str();

        cout << "Filename: " << fname << endl;

        //Part of code that considers the signal

        TFile* infile = new TFile(fname);
        TTree* tree = (TTree*)infile->Get("demo/timeTree");

	Int_t nEvent = tree->GetEntries();

	//Define entries for light_particle
	std::vector<double>* px_in = new std::vector<double>;
	std::vector<double>* py_in = new std::vector<double>;
	std::vector<double>* eta_in = new std::vector<double>;
	std::vector<double>* phi_in = new std::vector<double>;
	std::vector<double>* e_in = new std::vector<double>;
	std::vector<double>* pdgId_in = new std::vector<double>;
	std::vector<double>* vx_in = new std::vector<double>;
	std::vector<double>* vy_in = new std::vector<double>;
	std::vector<double>* vz_in = new std::vector<double>;
	//Define entries for muon
	std::vector<double>* px_muon = new std::vector<double>;
	std::vector<double>* py_muon = new std::vector<double>;
	std::vector<double>* eta_muon = new std::vector<double>;
	std::vector<double>* phi_muon = new std::vector<double>;
	std::vector<bool>* looseid_muon = new std::vector<bool>;
	//Define entries for PfJet
	std::vector<double>* pt_jet_in = new std::vector<double>;
	std::vector<double>* eta_jet_in = new std::vector<double>;
	std::vector<double>* phi_jet_in = new std::vector<double>;
	std::vector<double>* e_jet_in = new std::vector<double>;
	std::vector<double>* echad_jet_in = new std::vector<double>;
	std::vector<double>* ehad_jet_in = new std::vector<double>;
	std::vector<double>* ecem_jet_in = new std::vector<double>;
	std::vector<double>* eem_jet_in = new std::vector<double>;
	//Define entries for sv
	std::vector<double>* ipfjet_in = new std::vector<double>;
	std::vector<double>* svdeta_in = new std::vector<double>;
	std::vector<double>* svdphi_in = new std::vector<double>;
	std::vector<double>* svmass_in = new std::vector<double>;
	std::vector<double>* svntracks_in = new std::vector<double>;
	std::vector<double>* svchi2_in = new std::vector<double>;
	std::vector<double>* svndof_in = new std::vector<double>;
	std::vector<double>* svdxy_in = new std::vector<double>;
	std::vector<double>* svdxysig_in = new std::vector<double>;
	std::vector<double>* svd3d_in = new std::vector<double>;
	std::vector<double>* svd3dsig_in = new std::vector<double>;
	Bool_t hlt_jet_in, looseid;

	//Define everything relevant 
	Double_t px, py, vx, vy, vz, eta, phi, pt, e, rad, pdgId, pt_jet, eta_jet, phi_jet, e_jet , echad_jet, ehad_jet, ecem_jet, eem_jet, index_jet;

	//event and event_fin refer to the total and reduced event counts respectively
	Double_t event = 0;
	Double_t event_fin = 0;
	Double_t event_zero = 0;
	Double_t event_zero_jet = 0;
	Double_t event_sv = 0;

	//Jet number and jets with respective label
	Double_t n_jet = 0;
	Double_t n_jet_one = 0;
	Double_t n_jet_two = 0;
	Double_t n_jet_double = 0;

	//Set branch addresses 
	tree->SetBranchAddress("lightparticle_px", &px_in);
	tree->SetBranchAddress("lightparticle_py", &py_in);
	tree->SetBranchAddress("lightparticle_vx", &vx_in);
	tree->SetBranchAddress("lightparticle_vy", &vy_in);
	tree->SetBranchAddress("lightparticle_vz", &vz_in);
	tree->SetBranchAddress("lightparticle_eta", &eta_in);
	tree->SetBranchAddress("lightparticle_phi", &phi_in);
	tree->SetBranchAddress("lightparticle_e", &e_in);
	tree->SetBranchAddress("lightparticle_pdgId", &pdgId_in);
	tree->SetBranchAddress("muon_px", &px_muon);
	tree->SetBranchAddress("muon_py", &py_muon);
	tree->SetBranchAddress("muon_eta", &eta_muon);
	tree->SetBranchAddress("muon_phi", &phi_muon);
	tree->SetBranchAddress("muon_looseId", &looseid_muon);
	tree->SetBranchAddress("pfJetPt", &pt_jet_in);
	tree->SetBranchAddress("pfJetEta", &eta_jet_in);
	tree->SetBranchAddress("pfJetPhi", &phi_jet_in);
	tree->SetBranchAddress("pfJetE", &e_jet_in);
	tree->SetBranchAddress("pfJetChargedHadEnergy", &echad_jet_in);
	tree->SetBranchAddress("pfJetNeutralHadEnergy", &ehad_jet_in);
	tree->SetBranchAddress("pfJetChargedEmEnergy", &ecem_jet_in);
	tree->SetBranchAddress("pfJetNeutralEmEnergy", &eem_jet_in);
	tree->SetBranchAddress("sv_iPfJet", &ipfjet_in);
	tree->SetBranchAddress("sv_deta", &svdeta_in);
	tree->SetBranchAddress("sv_dphi", &svdphi_in);
	tree->SetBranchAddress("sv_mass", &svmass_in);
	tree->SetBranchAddress("sv_ntracks", &svntracks_in);
	tree->SetBranchAddress("sv_chi2", &svchi2_in);
	tree->SetBranchAddress("sv_ndof", &svndof_in);
	tree->SetBranchAddress("sv_dxy", &svdxy_in);
	tree->SetBranchAddress("sv_dxysig", &svdxysig_in);
	tree->SetBranchAddress("sv_d3d", &svd3d_in);
	tree->SetBranchAddress("sv_d3dsig", &svd3dsig_in);
	tree->SetBranchAddress("HLT_Mu9_IP6", &hlt_jet_in);

	//Set branch statuses
	tree->SetBranchStatus("*", 0);
	tree->SetBranchStatus("lightparticle_px", 1);
	tree->SetBranchStatus("lightparticle_py", 1);
	tree->SetBranchStatus("lightparticle_vx", 1);
	tree->SetBranchStatus("lightparticle_vy", 1);
	tree->SetBranchStatus("lightparticle_vz", 1);
	tree->SetBranchStatus("lightparticle_eta", 1);
	tree->SetBranchStatus("lightparticle_phi", 1);
	tree->SetBranchStatus("lightparticle_e", 1);
	tree->SetBranchStatus("lightparticle_pdgId", 1);
	tree->SetBranchStatus("muon_px", 1);
	tree->SetBranchStatus("muon_py", 1);
	tree->SetBranchStatus("muon_eta", 1);
	tree->SetBranchStatus("muon_phi", 1);
	tree->SetBranchStatus("muon_looseId", 1);
	tree->SetBranchStatus("pfJetPt", 1);
	tree->SetBranchStatus("pfJetEta", 1);
	tree->SetBranchStatus("pfJetPhi", 1);
	tree->SetBranchStatus("pfJetE", 1);
	tree->SetBranchStatus("pfJetChargedHadEnergy", 1);
	tree->SetBranchStatus("pfJetNeutralHadEnergy", 1);
	tree->SetBranchStatus("pfJetChargedEmEnergy", 1);
	tree->SetBranchStatus("pfJetNeutralEmEnergy", 1);
	tree->SetBranchStatus("sv_iPfJet", 1);
	tree->SetBranchStatus("sv_deta", 1);
	tree->SetBranchStatus("sv_dphi", 1);
	tree->SetBranchStatus("sv_mass", 1);
	tree->SetBranchStatus("sv_ntracks", 1);
	tree->SetBranchStatus("sv_chi2", 1);
	tree->SetBranchStatus("sv_ndof", 1);
	tree->SetBranchStatus("sv_dxy", 1);
	tree->SetBranchStatus("sv_dxysig", 1);
	tree->SetBranchStatus("sv_d3d", 1);
	tree->SetBranchStatus("sv_d3dsig", 1);
	tree->SetBranchStatus("HLT_Mu9_IP6", 1);

	//Make histograms for the total signal
	TH1D* hist_pdgId = make_hist("hist_pdgId", "Particle type", 100, -25, 25, "Particle type", "Number of particles", kBlack, 2);
	TH1D* hist_frac_num = make_hist("hist_frac_num", "Number of particles matched to pfjets", 100, 0., 20., "Number of particles", "Number of events", kBlack, 2);
	TH1D* hist_frac = make_hist("hist_frac", "Fraction of particles matched to pfjets", 100, 0., 1.1, "Fraction", "Number of particles", kBlack, 2);
	TH1D* hist_frac_muon_num = make_hist("hist_frac_muon_num", "Number of muons matched to pfjets", 100, 0., 20., "Number of muons", "Number of events", kBlack, 2);
	TH1D* hist_frac_muon = make_hist("hist_frac_muon", "Fraction of muons matched to pfjets", 100, 0., 1.1, "Fraction", "Number of events", kBlack, 2);
	TH1D* hist_eta_muon = make_hist("hist_eta_muon", "Eta distribution of each muon", 100, -2.8, 2.8, "Eta", "Number of muons", kBlack, 2);
	TH1D* hist_phi_muon = make_hist("hist_phi_muon", "Phi distribution of each muon", 100, -4.5, 4.5, "Phi", "Number of muons", kBlack, 2);
	TH1D* hist_pt_muon = make_hist("hist_pt_muon", "Transverse momentum distribution of muons", 100, 0., 500., "Transverse momentum (GeV)", "Number of muons", kBlack, 2);
	TH1D* hist_eta_jet = make_hist("hist_eta_jet", "Eta distribution of each jet", 100, -2.8, 2.8, "Eta", "Number of jets", kBlack, 2);
	TH1D* hist_deta_jet = make_hist("hist_deta_jet", "deta distribution of each jet", 100, 0., 5., "dEta", "Number of jet pairs", kBlack, 2);
	TH1D* hist_phi_jet = make_hist("hist_phi_jet", "Phi distribution of each jet", 100, -4.5, 4.5, "Phi", "Number of jets", kBlack, 2);
	TH1D* hist_dphi_jet = make_hist("hist_dphi_jet", "dphi distribution of each jet", 100, 0, 7., "dPhi", "Number of jet pairs", kBlack, 2);
	TH1D* hist_e = make_hist("hist_e", "Energy distribution of jets", 100, 0., 5000., "Energy (GeV)", "Number of jets", kBlack, 2);
	TH1D* hist_echad = make_hist("hist_echad", "Charged hadron energy distribution", 100, 0., 1.1, "Energy as a fraction of total jet energy", "Number of jets", kBlack, 2);
	TH1D* hist_ehad = make_hist("hist_ehad", "Neutral hadron energy distribution", 100, 0., 1.1, "Energy as a fraction of total jet energy", "Number of jets", kBlack, 2);
	TH1D* hist_ecem = make_hist("hist_ecem", "Charged EM energy distribution", 100, 0., 1.1, "Energy as a fraction of total jet energy", "Number of jets", kBlack, 2);
	TH1D* hist_eem = make_hist("hist_eem", "Neutral EM energy distribution", 100, 0., 1.1, "Energy as a fraction of total jet energy", "Number of jets", kBlack, 2);
	TH1D* hist_pt = make_hist("hist_pt", "Transverse momentum distribution of jets", 100, 0., 1000., "Transverse momentum (GeV)", "Number of jets", kBlack, 2);
	TH1D* hist_pt_ratio = make_hist("hist_pt_ratio", "Transverse momentum ratio of gen particles to jets", 100, 0., 5., "Transverse momentum ratio", "Number of jets", kBlack, 2);
	TH1D* hist_pt_ratio_muon = make_hist("hist_pt_ratio_muon", "Transverse momentum ratio of muons to jets", 100, 0., 5., "Transverse momentum ratio", "Number of jets", kBlack, 2);
	TH1D* hist_sv = make_hist("hist_sv", "Secondary vertex distribution", 100, 0., 20., "Number of secondary vertices matched", "Number of jets", kBlack, 2);
	TH1D* hist_svdeta= make_hist("hist_svdeta", "deta distribution", 100, -0.50, 0.50, "dEta", "Number of secondary vertices", kBlack, 2);
	TH1D* hist_svdphi = make_hist("hist_svdphi", "dphi distribution", 100, -0.50, 0.50, "dPhi", "Number of secondary vertices", kBlack, 2);
	TH1D* hist_svmass = make_hist("hist_svmass", "Mass distribution", 100, 0., 21.5, "Mass", "Number of secondary vertices", kBlack, 2);
	TH1D* hist_svntracks = make_hist("hist_svntracks", "Ntracks distribution", 100, 0., 12., "Ntracks", "Number of secondary vertices", kBlack, 2);
	TH1D* hist_svchi2 = make_hist("hist_svchi2", "Chi2 distribution", 100, 0., 21., "Chi2", "Number of secondary vertices", kBlack, 2);
	TH1D* hist_svndof = make_hist("hist_svndof", "ndof distribution", 100, 0., 18., "ndof", "Number of secondary vertices", kBlack, 2);
	TH1D* hist_svchi2div = make_hist("hist_svchi2div", "Chi2 and ndof comparison", 100, 0., 10., "Abs(chi2-ndof)", "Number of secondary vertices", kBlack, 2);
	TH1D* hist_svdxy = make_hist("hist_svdxy", "dxy distribution", 100, 0., 25., "dxy", "Number of secondary vertices", kBlack, 2);
	TH1D* hist_svdxysig = make_hist("hist_svdxysig", "dxysig distribution", 100, 0., 2000., "dxysig", "Number of secondary vertices", kBlack, 2);
	TH1D* hist_svd3d = make_hist("hist_svd3d", "d3d distribution", 100, 0., 136., "d3d", "Number of secondary vertices", kBlack, 2);
	TH1D* hist_svd3dsig = make_hist("hist_svd3dsig", "d3dsig distribution", 100, 0., 2900., "d3dsig", "Number of secondary vertices", kBlack, 2);
	TH1D* hist_num = make_hist("hist_num", "Clustered jet count distribution", 100, 0., 20., "Number of jets clustered", "Number of events", kBlack, 2);


	for (Int_t k = 0; k < nEvent; k++) {

		event += 1;
		tree->GetEntry(k);

		//Does everything only if it passes the trigger. Include || hlt_jet_in == 0 to include all events
		if (hlt_jet_in == 1) {
			//cout << "Event number:" << k + 1 << endl;
			event_fin += 1;

			//Create the input Jet_array
			Jet_array* jet_in = new Jet_array;
			Int_t nEnt_jet = pt_jet_in->size();

			//Fills jet_in with jets that fulfill cuts
			//Takes the eta and pt of each jet and places the appropriate jets in jet_in
			for (Int_t iEnt_jet_one = 0; iEnt_jet_one < nEnt_jet; iEnt_jet_one++) {
				if (TMath::Abs(eta_jet_in->at(iEnt_jet_one)) < 2.4 && pt_jet_in->at(iEnt_jet_one) > 20.) {
					jet_in->pt_vec->push_back(pt_jet_in->at(iEnt_jet_one));
					jet_in->eta_vec->push_back(eta_jet_in->at(iEnt_jet_one));
					jet_in->phi_vec->push_back(phi_jet_in->at(iEnt_jet_one));
					jet_in->e_vec->push_back(e_jet_in->at(iEnt_jet_one));
					jet_in->echad_vec->push_back(echad_jet_in->at(iEnt_jet_one));
					jet_in->ehad_vec->push_back(ehad_jet_in->at(iEnt_jet_one));
					jet_in->ecem_vec->push_back(ecem_jet_in->at(iEnt_jet_one));
					jet_in->eem_vec->push_back(eem_jet_in->at(iEnt_jet_one));
					jet_in->index_vec->push_back(iEnt_jet_one * 1.0);
				}
			}

			//Create the input lightparticle array
			//light_in collects all the chosen particles
			//light_out collects all the clustered from light_in 
			Light_array* light_in = new Light_array;
			Light_array* light_out = new Light_array;
			Int_t nEnt = px_in->size();

			//Fills light_in with detected particles with nEnt = number of particles in each event
			for (Int_t iEnt = 0; iEnt < nEnt; iEnt++) {

				//Get attributes of iEnt'th particle 
				px = px_in->at(iEnt);
				py = py_in->at(iEnt);
				vx = vx_in->at(iEnt);
				vy = vy_in->at(iEnt);
				vz = vz_in->at(iEnt);
				eta = eta_in->at(iEnt);
				phi = phi_in->at(iEnt);
				e = e_in->at(iEnt);
				pdgId = pdgId_in->at(iEnt);
				pt = TMath::Sqrt(px * px + py * py);
				rad = TMath::Sqrt(vx * vx + vy * vy);
				
				//Fills the attributes of the iEnt'th particle of the event into light_in if requirements are met
				if (pt > 2. && rad < 120. && TMath::Abs(vz) < 300. && TMath::Abs(eta) < 2.4 && TMath::Abs(pdgId) < 1000) {
					light_in->px_vec->push_back(px);
					light_in->py_vec->push_back(py);
					light_in->pt_vec->push_back(pt);
					light_in->eta_vec->push_back(eta);
					light_in->phi_vec->push_back(phi);
					light_in->pdgId_vec->push_back(pdgId);
				}
			}
			
			//Create the input muon array
			//muon_in collects all the chosen muons
			//muon_out collects all the clustered muons from muon_in
			Light_array* muon_in = new Light_array;
			Light_array* muon_out = new Light_array;
			Int_t nEnt_muon = px_muon->size();

			//Fills muon_in with detected particles with nEnt_muon = number of muons in each event
			for (Int_t iEnt_muon = 0; iEnt_muon < nEnt_muon; iEnt_muon++) {
				
				//Get attributes of iEnt_muon'th muon
				px = px_muon->at(iEnt_muon);
				py = py_muon->at(iEnt_muon);
				pt = TMath::Sqrt(px * px + py * py);
				eta = eta_muon->at(iEnt_muon);
				phi = phi_muon->at(iEnt_muon);
				looseid = looseid_muon->at(iEnt_muon);
				
				//Fills the attributes of the iEnt_muon'th muon of the event into muon_in if requirements are met
				if (pt > 5. && TMath::Abs(eta) < 2.4 && looseid == 1 ) {
					muon_in->px_vec->push_back(px);
					muon_in->py_vec->push_back(py);
					muon_in->pt_vec->push_back(pt);
					muon_in->eta_vec->push_back(eta);
					muon_in->phi_vec->push_back(phi);
					muon_in->pdgId_vec->push_back(13);
				}
			}

			//Initialise values used for clustering along with the sizes of the initial light_in and jet_in vectors
			Int_t nEnt_jet_one, nEnt_one;
			Int_t nEnt_muon_initial = muon_in->eta_vec->size();
			Int_t nEnt_initial = light_in->eta_vec->size();
			Int_t nEnt_jet_initial = jet_in->pt_vec->size();
			Double_t match = 0;
			Double_t match_muon = 0;
			Double_t match_muon_num = 0;
			Double_t jet_num = 0;
			Double_t px_clus_muon = 0;
			Double_t py_clus_muon = 0;
			Double_t pt_clus_muon;
			Int_t sv_size = ipfjet_in->size();

			//Records events with no particles
			if (nEnt_initial == 0) {
				event_zero += 1;
			}

			//Prevents clustering early on for any events with no jets
			if (nEnt_jet_initial == 0) {
				event_zero_jet += 1;
				continue;
			}

			//Cluster for as many elements in jet_in
			for (Int_t i = 0; i < nEnt_jet_initial; i++) {
				nEnt_jet_one = jet_in->pt_vec->size();    //Number of jets
				nEnt_one = light_in->eta_vec->size();     //Number of particles
				
				//Obtains the jet with the highest pt and stores its position in the array of jets
				if (nEnt_jet_one > 0) {
					Double_t pt_max = *std::max_element(jet_in->pt_vec->begin(), jet_in->pt_vec->end());
					Int_t pt_max_ent = std::max_element(jet_in->pt_vec->begin(), jet_in->pt_vec->end()) - jet_in->pt_vec->begin();
					jet_in->max_ent = pt_max_ent;
				}
				//Stops clustering if jets are exhausted
				else if (nEnt_jet_one == 0) {
					break;
				}
				

				Double_t sv_count = 0;
				//Sv clustering (even for jets without particles linked to them)
				for (Int_t j = 0; j < sv_size; j++) {
					Double_t sv_index = ipfjet_in->at(j);

					//sv_index is the index of the jet associated with the secondary vertex
					//If sv_index matches with the index of the leading pt jet, fill histogram with the vertex attributes
					if (sv_index == jet_in->index_vec->at(jet_in->max_ent) && svdxysig_in->at(j) > 50) {
						sv_count += 1;

						hist_svdeta->Fill(svdeta_in->at(j));
						hist_svdphi->Fill(svdphi_in->at(j));
						hist_svmass->Fill(svmass_in->at(j));
						hist_svntracks->Fill(svntracks_in->at(j));
						hist_svndof->Fill(svndof_in->at(j));
						hist_svchi2->Fill(svchi2_in->at(j));
						hist_svchi2div->Fill(TMath::Abs(svchi2_in->at(j) - svndof_in->at(j)));
						hist_svdxy->Fill(svdxy_in->at(j));
						hist_svdxysig->Fill(svdxysig_in->at(j));
						hist_svd3d->Fill(svd3d_in->at(j));
						hist_svd3dsig->Fill(svd3dsig_in->at(j));
					}
				}
				
				//Number of secondary vertices for the leading jet
				hist_sv->Fill(sv_count);

				//Get reduced jet for replacement later
				Jet_array* jet_new = jet_reduce(jet_in);
				
				//Particle clustering along with assignment of output
				//light_out is the first output of the cluster which gives the particles clustered to the leading jet
				//The second output of the cluster is the set of unclustered particles
				std::vector<Light_array*>* light_vector = cluster(light_in, jet_in);
				light_out = light_vector->at(0);

				//Muon clustering
				//Similar to light_out in behaviour
				std::vector<Light_array*>* muon_vector = cluster(muon_in, jet_in);
				muon_out = muon_vector->at(0);
				
				//Fills histogram with attributes of clustered muons
				//Also sums up px and py of all the clustered muons as a collective
				for (Int_t a = 0; a < muon_out->pt_vec->size(); a++) {
					px_clus_muon += muon_out->px_vec->at(a);
					py_clus_muon += muon_out->py_vec->at(a);
					hist_pt_muon->Fill(muon_out->pt_vec->at(a));
					hist_eta_muon->Fill(muon_out->eta_vec->at(a));
					hist_phi_muon->Fill(muon_out->phi_vec->at(a));
				}
				
				//Fills histogram with particle type of clustered particles
				for (Int_t b = 0; b < light_out->pt_vec->size(); b++) {
			
					hist_pdgId->Fill(light_out->pdgId_vec->at(b));
				}
				
				//Obtains total pt of clustered muons
				pt_clus_muon = TMath::Sqrt(px_clus_muon * px_clus_muon + py_clus_muon * py_clus_muon);

				//Compares total muon pt with jet pt as a fraction
				if (muon_out->pt_vec->size() > 0) {
					Double_t pt_ratio_muon = pt_clus_muon / jet_in->pt_vec->at(jet_in->max_ent);
					hist_pt_ratio_muon->Fill(pt_ratio_muon);
				}

				//Match gets the fraction of particles clustered with respect to the total number of particles from light_in
				if (nEnt_initial > 0) {
					match += (light_out->eta_vec->size() / (nEnt_initial * 1.0));
				}
				if (nEnt_muon_initial > 0) {
					match_muon += (muon_out->eta_vec->size() / (nEnt_muon_initial * 1.0));
					match_muon_num += (muon_out->eta_vec->size());
				}

				//Fills histogram with leading jet values
				n_jet += 1;
				hist_pt->Fill(jet_in->pt_vec->at(jet_in->max_ent));
				hist_e->Fill(jet_in->e_vec->at(jet_in->max_ent));
				hist_echad->Fill(jet_in->echad_vec->at(jet_in->max_ent)/ jet_in->e_vec->at(jet_in->max_ent));
				hist_ehad->Fill(jet_in->ehad_vec->at(jet_in->max_ent)/ jet_in->e_vec->at(jet_in->max_ent));
				hist_ecem->Fill(jet_in->ecem_vec->at(jet_in->max_ent)/ jet_in->e_vec->at(jet_in->max_ent));
				hist_eem->Fill(jet_in->eem_vec->at(jet_in->max_ent)/ jet_in->e_vec->at(jet_in->max_ent));
				hist_eta_jet->Fill(jet_in->eta_vec->at(jet_in->max_ent));
				hist_phi_jet->Fill(jet_in->phi_vec->at(jet_in->max_ent));


				//Clear out light_in and muon_in right before reassignment since reassigning pointers doesn't remove the initial data stored in light_in
				delete light_in;
				delete muon_in;

				//Reassign for next cluster with remainder particles and muons
				light_in = light_vector->at(1);
				muon_in = muon_vector->at(1);
				
				//Clear out light_vector and muon_vector since this will be called as a new pointer in the next cluster
				delete light_vector;
				delete muon_vector;

				//Clear out light_in and muon_in right before reassignment since reassigning pointers doesn't remove the initial data stored in light_in
				delete jet_in;
				jet_in = jet_new;

				//Records if the jet contains particles clustered around it
				if (light_out->eta_vec->size() > 0) {
					jet_num += 1;
				}

				delete light_out;
				delete muon_out;
			}
			//Above completes the clustering of the particles to all accepted jets for current event

			//Delete light_in, muon_in and jet_in since it will be called as 'new' in next event
			//jet_in = jet_new ensures that jet_new also gets deleted
			delete light_in;
			delete muon_in;
			delete jet_in;

			//Fill histogram with number of matched particles per jet and number of jets
			if (nEnt_initial > 0) {
				hist_frac->Fill(match);
				hist_num->Fill(jet_num);
			}
			if (nEnt_muon_initial > 0) {
				hist_frac_muon->Fill(match_muon);
				hist_frac_muon_num->Fill(match_muon_num);
			}

			
		}

	}
	//Clear vectors which will be filled out by branches
	delete px_in;
	delete py_in;
	delete eta_in;
	delete phi_in;
	delete e_in;
	delete pdgId_in;
	delete vx_in;
	delete vy_in;
	delete vz_in;
	delete px_muon;
	delete py_muon;
	delete eta_muon;
	delete phi_muon;
	delete looseid_muon;
	delete pt_jet_in;
	delete eta_jet_in;
	delete phi_jet_in;
	delete e_jet_in;
	delete echad_jet_in;
	delete ehad_jet_in;
	delete ecem_jet_in;
	delete eem_jet_in;
	delete ipfjet_in;
	delete svdeta_in;
	delete svdphi_in;
	delete svmass_in;
	delete svntracks_in;
	delete svchi2_in;
	delete svndof_in;
	delete svdxy_in;
	delete svdxysig_in;
	delete svd3d_in;
	delete svd3dsig_in;


	cout << "Total event count: " << event << endl;
	cout << "Accepted event count:" << event_fin << endl;
	cout << "Zero particle event count:" << event_zero << endl;
	cout << "Zero jet event count:" << event_zero_jet << endl;
	cout << "Total jet count:" << n_jet << endl;
	cout << "___________________________" << endl;
	cout << "" << endl;
	cout << "Percentage analysis: " << endl;
	cout << "Percentage of events accepted:                 " << (event_fin / event) * 100 << endl; cout << "Number of jets per accepted event:             " << (n_jet / event_fin) << endl;
	cout << "___________________________" << endl;
	cout << "" << endl;

	//Part of code that considers the background
	TFile* infile1 = new TFile("../Ntuples/QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_MINIAODSIM_v1_generationQCD.root");
        TTree* tree1 = (TTree*)infile1->Get("demo/timeTree");

	Int_t nEvent1 = tree1->GetEntries();
	Int_t batchsize = nEvent1;
	Int_t batchcount = nEvent1 / batchsize;
	//Int_t batchcount = 0;

	//Define entries for muon
	std::vector<double>* px_muon1 = new std::vector<double>;
	std::vector<double>* py_muon1 = new std::vector<double>;
	std::vector<double>* eta_muon1 = new std::vector<double>;
	std::vector<double>* phi_muon1 = new std::vector<double>;
	std::vector<bool>* looseid_muon1 = new std::vector<bool>;
	//Define entries for PfJet
	std::vector<double>* pt_jet_in1 = new std::vector<double>;
	std::vector<double>* eta_jet_in1 = new std::vector<double>;
	std::vector<double>* phi_jet_in1 = new std::vector<double>;
	std::vector<double>* e_jet_in1 = new std::vector<double>;
	std::vector<double>* echad_jet_in1 = new std::vector<double>;
	std::vector<double>* ehad_jet_in1 = new std::vector<double>;
	std::vector<double>* ecem_jet_in1 = new std::vector<double>;
	std::vector<double>* eem_jet_in1 = new std::vector<double>;
	//Define entries for sv
	std::vector<double>* ipfjet_in1 = new std::vector<double>;
	std::vector<double>* svdeta_in1 = new std::vector<double>;
	std::vector<double>* svdphi_in1 = new std::vector<double>;
	std::vector<double>* svmass_in1 = new std::vector<double>;
	std::vector<double>* svntracks_in1 = new std::vector<double>;
	std::vector<double>* svchi2_in1 = new std::vector<double>;
	std::vector<double>* svndof_in1 = new std::vector<double>;
	std::vector<double>* svdxy_in1 = new std::vector<double>;
	std::vector<double>* svdxysig_in1 = new std::vector<double>;
	std::vector<double>* svd3d_in1 = new std::vector<double>;
	std::vector<double>* svd3dsig_in1 = new std::vector<double>;
	Bool_t hlt_jet_in1, looseid1;

	//Define everything relevant 
	Double_t px1, py1, vx1, vy1, vz1, eta1, phi1, pt1, e1, rad1, pdgId1, pt_jet1, eta_jet1, phi_jet1, e_jet1, echad_jet1, ehad_jet1, ecem_jet1, eem_jet1, index_jet1;

	//event and event_fin refer to the total and reduced event counts respectively
	Double_t event1 = 0;
	Double_t event_fin1 = 0;
	Double_t event_zero_jet1 = 0;
	Double_t event_sv1 = 0;

	//Jet number and jets with respective label
	Double_t n_jet1 = 0;

	//Set branch addresses 
	tree1->SetBranchAddress("muon_px", &px_muon1);
	tree1->SetBranchAddress("muon_py", &py_muon1);
	tree1->SetBranchAddress("muon_eta", &eta_muon1);
	tree1->SetBranchAddress("muon_phi", &phi_muon1);
	tree1->SetBranchAddress("muon_looseId", &looseid_muon1);
	tree1->SetBranchAddress("pfJetPt", &pt_jet_in1);
	tree1->SetBranchAddress("pfJetEta", &eta_jet_in1);
	tree1->SetBranchAddress("pfJetPhi", &phi_jet_in1);
	tree1->SetBranchAddress("pfJetE", &e_jet_in1);
	tree1->SetBranchAddress("pfJetChargedHadEnergy", &echad_jet_in1);
	tree1->SetBranchAddress("pfJetNeutralHadEnergy", &ehad_jet_in1);
	tree1->SetBranchAddress("pfJetChargedEmEnergy", &ecem_jet_in1);

	//Set branch addresses 
	tree1->SetBranchAddress("muon_px", &px_muon1);
	tree1->SetBranchAddress("muon_py", &py_muon1);
	tree1->SetBranchAddress("muon_eta", &eta_muon1);
	tree1->SetBranchAddress("muon_phi", &phi_muon1);
	tree1->SetBranchAddress("muon_looseId", &looseid_muon1);
	tree1->SetBranchAddress("pfJetPt", &pt_jet_in1);
	tree1->SetBranchAddress("pfJetEta", &eta_jet_in1);
	tree1->SetBranchAddress("pfJetPhi", &phi_jet_in1);
	tree1->SetBranchAddress("pfJetE", &e_jet_in1);
	tree1->SetBranchAddress("pfJetChargedHadEnergy", &echad_jet_in1);
	tree1->SetBranchAddress("pfJetNeutralHadEnergy", &ehad_jet_in1);
	tree1->SetBranchAddress("pfJetChargedEmEnergy", &ecem_jet_in1);
	tree1->SetBranchAddress("pfJetNeutralEmEnergy", &eem_jet_in1);
	tree1->SetBranchAddress("sv_iPfJet", &ipfjet_in1);
	tree1->SetBranchAddress("sv_deta", &svdeta_in1);
	tree1->SetBranchAddress("sv_dphi", &svdphi_in1);
	tree1->SetBranchAddress("sv_mass", &svmass_in1);
	tree1->SetBranchAddress("sv_ntracks", &svntracks_in1);
	tree1->SetBranchAddress("sv_chi2", &svchi2_in1);
	tree1->SetBranchAddress("sv_ndof", &svndof_in1);
	tree1->SetBranchAddress("sv_dxy", &svdxy_in1);
	tree1->SetBranchAddress("sv_dxysig", &svdxysig_in1);
	tree1->SetBranchAddress("sv_d3d", &svd3d_in1);
	tree1->SetBranchAddress("sv_d3dsig", &svd3dsig_in1);
	tree1->SetBranchAddress("HLT_Mu9_IP6", &hlt_jet_in1);

	//Set branch statuses
	tree1->SetBranchStatus("*", 0);
	tree1->SetBranchStatus("muon_px", 1);
	tree1->SetBranchStatus("muon_py", 1);
	tree1->SetBranchStatus("muon_eta", 1);
	tree1->SetBranchStatus("muon_phi", 1);
	tree1->SetBranchStatus("muon_looseId", 1);
	tree1->SetBranchStatus("pfJetPt", 1);
	tree1->SetBranchStatus("pfJetEta", 1);
	tree1->SetBranchStatus("pfJetPhi", 1);
	tree1->SetBranchStatus("pfJetE", 1);
	tree1->SetBranchStatus("pfJetChargedHadEnergy", 1);
	tree1->SetBranchStatus("pfJetNeutralHadEnergy", 1);
	tree1->SetBranchStatus("pfJetChargedEmEnergy", 1);
	tree1->SetBranchStatus("pfJetNeutralEmEnergy", 1);
	tree1->SetBranchStatus("sv_iPfJet", 1);
	tree1->SetBranchStatus("sv_deta", 1);
	tree1->SetBranchStatus("sv_dphi", 1);
	tree1->SetBranchStatus("sv_mass", 1);
	tree1->SetBranchStatus("sv_ntracks", 1);
	tree1->SetBranchStatus("sv_chi2", 1);
	tree1->SetBranchStatus("sv_ndof", 1);
	tree1->SetBranchStatus("sv_dxy", 1);
	tree1->SetBranchStatus("sv_dxysig", 1);
	tree1->SetBranchStatus("sv_d3d", 1);
	tree1->SetBranchStatus("sv_d3dsig", 1);
	tree1->SetBranchStatus("HLT_Mu9_IP6", 1);

	//Make histograms for the total background
	TH1D* hist_frac_muon_num1 = make_hist("hist_frac_muon_num1", "Number of muons matched to pfjets", 100, 0., 20., "Number of muons", "Number of events", kGreen, 2);
	TH1D* hist_frac_muon1 = make_hist("hist_frac_muon1", "Fraction of muons matched to pfjets", 100, 0., 1.1, "Fraction", "Number of events", kGreen, 2);
	TH1D* hist_eta_muon1 = make_hist("hist_eta_muon1", "Eta distribution of each muon", 100, -2.8, 2.8, "Eta", "Number of muons", kGreen, 2);
	TH1D* hist_phi_muon1 = make_hist("hist_phi_muon1", "Phi distribution of each muon", 100, -4.5, 4.5, "Phi", "Number of muons", kGreen, 2);
	TH1D* hist_pt_muon1 = make_hist("hist_pt_muon1", "Transverse momentum distribution of muons", 100, 0., 500., "Transverse momentum (GeV)", "Number of muons", kGreen,2);
	TH1D* hist_eta_jet1 = make_hist("hist_eta_jet1", "Eta distribution of each jet", 100, -2.8, 2.8, "Eta", "Number of jets", kGreen, 2);
	TH1D* hist_deta_jet1 = make_hist("hist_deta_jet1", "deta distribution of each jet", 100, 0., 5., "dEta", "Number of jet pairs", kGreen, 2);
	TH1D* hist_phi_jet1 = make_hist("hist_phi_jet1", "Phi distribution of each jet", 100, -4.5, 4.5, "Phi", "Number of jets", kGreen, 2);
	TH1D* hist_dphi_jet1 = make_hist("hist_dphi_jet1", "dphi distribution of each jet", 100, 0, 7., "dPhi", "Number of jet pairs", kGreen, 2);
	TH1D* hist_e1 = make_hist("hist_e1", "Energy distribution of jets", 100, 0., 5000., "Energy (GeV)", "Number of jets", kGreen, 2);
	TH1D* hist_echad1 = make_hist("hist_echad1", "Charged hadron energy distribution", 100, 0., 1.1, "Energy as a fraction of total jet energy", "Number of jets", kGreen, 2);
	TH1D* hist_ehad1 = make_hist("hist_ehad1", "Neutral hadron energy distribution", 100, 0., 1.1, "Energy as a fraction of total jet energy", "Number of jets", kGreen, 2);
	TH1D* hist_ecem1 = make_hist("hist_ecem1", "Charged EM energy distribution", 100, 0., 1.1, "Energy as a fraction of total jet energy", "Number of jets", kGreen, 2);
	TH1D* hist_eem1 = make_hist("hist_eem1", "Neutral EM energy distribution", 100, 0., 1.1, "Energy as a fraction of total jet energy", "Number of jets", kGreen, 2);
	TH1D* hist_pt1 = make_hist("hist_pt1", "Transverse momentum distribution of jets", 100, 0., 1000., "Transverse momentum (GeV)", "Number of jets", kGreen, 2);
	TH1D* hist_pt_ratio_muon1 = make_hist("hist_pt_ratio_muon1", "Transverse momentum ratio of muons to jets", 100, 0., 5., "Transverse momentum ratio", "Number of jets", kGreen, 2);
	TH1D* hist_sv1 = make_hist("hist_sv1", "Secondary vertex distribution", 100, 0., 20., "Number of secondary vertices matched", "Number of jets", kGreen, 2);
	TH1D* hist_svdeta1 = make_hist("hist_svdeta1", "deta distribution", 100, -0.50, 0.50, "dEta", "Number of secondary vertices", kGreen, 2);
	TH1D* hist_svdphi1 = make_hist("hist_svdphi1", "dphi distribution", 100, -0.50, 0.50, "dPhi", "Number of secondary vertices", kGreen, 2);
	TH1D* hist_svmass1 = make_hist("hist_svmass1", "Mass distribution", 100, 0., 21.5, "Mass", "Number of secondary vertices", kGreen, 2);
	TH1D* hist_svntracks1 = make_hist("hist_svntracks1", "Ntracks distribution", 100, 0., 12., "Ntracks", "Number of secondary vertices", kGreen, 2);
	TH1D* hist_svchi21 = make_hist("hist_svchi21", "Chi2 distribution", 100, 0., 21., "Chi2", "Number of secondary vertices", kGreen, 2);
	TH1D* hist_svndof1 = make_hist("hist_svndof1", "ndof distribution", 100, 0., 18., "ndof", "Number of secondary vertices", kGreen, 2);
	TH1D* hist_svchi2div1 = make_hist("hist_svchi2div1", "Chi2 and ndof comparison", 100, 0., 10., "Abs(chi2-ndof)", "Number of secondary vertices", kGreen, 2);
	TH1D* hist_svdxy1 = make_hist("hist_svdxy1", "dxy distribution", 100, 0., 25., "dxy", "Number of secondary vertices", kGreen, 2);
	TH1D* hist_svdxysig1 = make_hist("hist_svdxysig1", "dxysig distribution", 100, 0., 2000., "dxysig", "Number of secondary vertices", kGreen, 2);
	TH1D* hist_svd3d1 = make_hist("hist_svd3d1", "d3d distribution", 100, 0., 136., "d3d", "Number of secondary vertices", kGreen, 2);
	TH1D* hist_svd3dsig1 = make_hist("hist_svd3dsig1", "d3dsig distribution", 100, 0., 2900., "d3dsig", "Number of secondary vertices", kGreen, 2);

	for (Int_t batch = 0; batch < batchcount; batch++) {
		cout << "Batch number:" << batch << endl;
		for (Int_t k = batch*batchsize; k < (batch+1)*batchsize; k++) {

			event1 += 1;
			tree1->GetEntry(k);

			//Does everything only if it passes the trigger. Include || hlt_jet_in == 0 to include all events
			if (hlt_jet_in1 == 1) {
				event_fin1 += 1;

				//cout << "Event number:" << k + 1 << endl;


				//Create the input Jet_array
				Jet_array* jet_in1 = new Jet_array;
				Int_t nEnt_jet1 = pt_jet_in1->size();

				//Fills jet_in1 with jets that fulfill cuts
				//Takes the eta and pt of each jet and places the appropriate jets in jet_in1
				for (Int_t iEnt_jet_one = 0; iEnt_jet_one < nEnt_jet1; iEnt_jet_one++) {
					if (TMath::Abs(eta_jet_in1->at(iEnt_jet_one)) < 2.4 && pt_jet_in1->at(iEnt_jet_one) > 20.) {
						jet_in1->pt_vec->push_back(pt_jet_in1->at(iEnt_jet_one));
						jet_in1->eta_vec->push_back(eta_jet_in1->at(iEnt_jet_one));
						jet_in1->phi_vec->push_back(phi_jet_in1->at(iEnt_jet_one));
						jet_in1->e_vec->push_back(e_jet_in1->at(iEnt_jet_one));
						jet_in1->echad_vec->push_back(echad_jet_in1->at(iEnt_jet_one));
						jet_in1->ehad_vec->push_back(ehad_jet_in1->at(iEnt_jet_one));
						jet_in1->ecem_vec->push_back(ecem_jet_in1->at(iEnt_jet_one));
						jet_in1->eem_vec->push_back(eem_jet_in1->at(iEnt_jet_one));
						jet_in1->index_vec->push_back(iEnt_jet_one * 1.0);
					}
				}

				//Create the input muon array
				//muon_in1 collects all the chosen muons
				//muon_out1 collects all the clustered muons from muon_in
				Light_array* muon_in1 = new Light_array;
				Light_array* muon_out1 = new Light_array;;
				Int_t nEnt_muon1 = px_muon1->size();

				//Fills muon_in with detected particles with nEnt_muon = number of muons in each event
				for (Int_t iEnt_muon1 = 0; iEnt_muon1 < nEnt_muon1; iEnt_muon1++) {
					
					//Get attributes of iEnt_muon'th muon
					px1 = px_muon1->at(iEnt_muon1);
					py1 = py_muon1->at(iEnt_muon1);
					pt1 = TMath::Sqrt(px1 * px1 + py1 * py1);
					eta1 = eta_muon1->at(iEnt_muon1);
					phi1 = phi_muon1->at(iEnt_muon1);
					looseid1 = looseid_muon1->at(iEnt_muon1);

					//Fills the attributes of the iEnt_muon'th muon of the event into muon_in if requirements are met
					if (pt1 > 5. && TMath::Abs(eta1) < 2.4 && looseid1 == 1) {
						muon_in1->px_vec->push_back(px1);
						muon_in1->py_vec->push_back(py1);
						muon_in1->pt_vec->push_back(pt1);
						muon_in1->eta_vec->push_back(eta1);
						muon_in1->phi_vec->push_back(phi1);
						muon_in1->pdgId_vec->push_back(13);
					}
				}



				//Initialise values used for clustering along with the sizes of the initial muon_in and jet_in vectors
				Int_t nEnt_jet_one1, nEnt_one1;
				Int_t nEnt_muon_initial1 = muon_in1->eta_vec->size();
				Int_t nEnt_jet_initial1 = jet_in1->pt_vec->size();
				Double_t match_muon1 = 0;
				Double_t match_muon_num1 = 0;
				Double_t px_clus_muon1 = 0;
				Double_t py_clus_muon1 = 0;
				Double_t pt_clus_muon1;
				Int_t sv_size1 = ipfjet_in1->size();

				//Prevents clustering early on for any events with no jets
				if (nEnt_jet_initial1 == 0) {
					event_zero_jet1 += 1;
					continue;
				}

				//Cluster for as many elements in jet_in
				for (Int_t i = 0; i < nEnt_jet_initial1; i++) {
					nEnt_jet_one1 = jet_in1->pt_vec->size();
					
					//Obtains the jet with the highest pt and stores its position in the array of jets
					if (nEnt_jet_one1 > 0) {
						Double_t pt_max1 = *std::max_element(jet_in1->pt_vec->begin(), jet_in1->pt_vec->end());
						Int_t pt_max_ent1 = std::max_element(jet_in1->pt_vec->begin(), jet_in1->pt_vec->end()) - jet_in1->pt_vec->begin();
						jet_in1->max_ent = pt_max_ent1;
					}
					//Stops clustering if jets are exhausted
					else if (nEnt_jet_one1 == 0) {
						break;
					}


					Double_t sv_count1 = 0;
					//Sv clustering (even for jets without muons linked to them)
					//sv_index is the index of the jet associated with the secondary vertex
					//If sv_index matches with the index of the leading pt jet, fill histogram with the vertex attributes
					for (Int_t j = 0; j < sv_size1; j++) {
						Double_t sv_index1 = ipfjet_in1->at(j);
						if (sv_index1 == jet_in1->index_vec->at(jet_in1->max_ent) && svdxysig_in1->at(j) > 50.) {
							sv_count1 += 1;

							hist_svdeta1->Fill(svdeta_in1->at(j));
							hist_svdphi1->Fill(svdphi_in1->at(j));
							hist_svmass1->Fill(svmass_in1->at(j));
							hist_svntracks1->Fill(svntracks_in1->at(j));
							hist_svchi21->Fill(svchi2_in1->at(j));
							hist_svndof1->Fill(svndof_in1->at(j));
							hist_svchi2div1->Fill(TMath::Abs(svchi2_in1->at(j) - svndof_in1->at(j)));
							hist_svdxy1->Fill(svdxy_in1->at(j));
							hist_svdxysig1->Fill(svdxysig_in1->at(j));
							hist_svd3d1->Fill(svd3d_in1->at(j));
							hist_svd3dsig1->Fill(svd3dsig_in1->at(j));
						}
					}
					hist_sv1->Fill(sv_count1);

					//Get reduced jet for replacement later
					Jet_array* jet_new1 = jet_reduce(jet_in1);

					//Muon clustering
					std::vector<Light_array*>* muon_vector1 = cluster(muon_in1, jet_in1);
					muon_out1 = muon_vector1->at(0);

					//Fills histogram with attributes of clustered muons
					//Also sums up px and py of all the clustered muons as a collective
					for (Int_t a = 0; a < muon_out1->pt_vec->size(); a++) {
						px_clus_muon1 += muon_out1->px_vec->at(a);
						py_clus_muon1 += muon_out1->py_vec->at(a);
						hist_pt_muon1->Fill(muon_out1->pt_vec->at(a));
						hist_eta_muon1->Fill(muon_out1->eta_vec->at(a));
						hist_phi_muon1->Fill(muon_out1->phi_vec->at(a));
					}
					if (nEnt_muon_initial1 > 0) {
						match_muon1 += (muon_out1->eta_vec->size() / (nEnt_muon_initial1 * 1.0));
						match_muon_num1 += (muon_out1->eta_vec->size());
					}

					pt_clus_muon1 = TMath::Sqrt(px_clus_muon1 * px_clus_muon1 + py_clus_muon1 * py_clus_muon1);

					//Compares total muon pt with jet pt as a fraction
					if (muon_out1->pt_vec->size() > 0) {
						Double_t pt_ratio_muon1 = pt_clus_muon1 / jet_in1->pt_vec->at(jet_in1->max_ent);
						hist_pt_ratio_muon1->Fill(pt_ratio_muon1);
					}
 
					//Fills histogram with leading jet values
					n_jet1 += 1;
					hist_pt1->Fill(jet_in1->pt_vec->at(jet_in1->max_ent));
					hist_e1->Fill(jet_in1->e_vec->at(jet_in1->max_ent));
					hist_echad1->Fill(jet_in1->echad_vec->at(jet_in1->max_ent) / jet_in1->e_vec->at(jet_in1->max_ent));
					hist_ehad1->Fill(jet_in1->ehad_vec->at(jet_in1->max_ent) / jet_in1->e_vec->at(jet_in1->max_ent));
					hist_ecem1->Fill(jet_in1->ecem_vec->at(jet_in1->max_ent) / jet_in1->e_vec->at(jet_in1->max_ent));
					hist_eem1->Fill(jet_in1->eem_vec->at(jet_in1->max_ent) / jet_in1->e_vec->at(jet_in1->max_ent));
					hist_eta_jet1->Fill(jet_in1->eta_vec->at(jet_in1->max_ent));
					hist_phi_jet1->Fill(jet_in1->phi_vec->at(jet_in1->max_ent));

					//Clear out muon_in right before reassignment since reassigning pointers doesn't remove the initial data stored in muon_in
					delete muon_in1;
					delete muon_out1;

					//Reassign for next cluster
					muon_in1 = muon_vector1->at(1);

					//Clear out jet_in right before reassignment to jet_new (jet_new can't be cleared since this pointer will be used in next cluster)
					delete muon_vector1;
					delete jet_in1;
					jet_in1 = jet_new1;
				}
				if (nEnt_muon_initial1 > 0) {
					hist_frac_muon1->Fill(match_muon1);
					hist_frac_muon_num1->Fill(match_muon_num1);
				}


				//Delete muon_in and jet_in since it will be called as 'new' in next event
				//jet_in = jet_new ensures that jet_new also gets deleted
				delete muon_in1;
				delete jet_in1;
			}
			else
				continue;
		}
	}
	cout << "Background" << endl;
	cout << "" << endl;
	cout << "Total event count: " << event1 << endl;
	cout << "Accepted event count:" << event_fin1 << endl;
	cout << "Zero jet event count:" << event_zero_jet1 << endl;
	cout << "Total jet count:" << n_jet1 << endl;
	cout << "___________________________" << endl;
	cout << "" << endl;
	cout << "Percentage analysis: " << endl;
	cout << "Percentage of events accepted:                 " << (event_fin1 / event1) * 100 << endl;
	cout << "Number of jets per accepted event:             " << (n_jet1 / event_fin1) << endl;
	
	delete px_muon1;
        delete py_muon1;
        delete eta_muon1;
        delete phi_muon1;
        delete looseid_muon1;
        delete pt_jet_in1;
        delete eta_jet_in1;
        delete phi_jet_in1;
        delete e_jet_in1;
        delete echad_jet_in1;
        delete ehad_jet_in1;
        delete ecem_jet_in1;
        delete eem_jet_in1;
        delete ipfjet_in1;
        delete svdeta_in1;
        delete svdphi_in1;
        delete svmass_in1;
        delete svntracks_in1;
        delete svchi2_in1;
        delete svndof_in1;
        delete svdxy_in1;
        delete svdxysig_in1;
        delete svd3d_in1;
        delete svd3dsig_in1;


	Double_t factor = 1.;

	hist_pdgId->Scale(factor / hist_pdgId->Integral());
	hist_frac->Scale(factor / hist_frac->Integral());
	hist_frac_muon->Scale(factor / hist_frac_muon->Integral());
	hist_frac_muon_num->Scale(factor / hist_frac_muon_num->Integral());
	hist_eta_muon->Scale(factor / hist_eta_muon->Integral());
	hist_phi_muon->Scale(factor / hist_phi_muon->Integral());
	hist_pt_muon->Scale(factor / hist_pt_muon->Integral());
	hist_eta_jet->Scale(factor / hist_eta_jet->Integral());
	hist_phi_jet->Scale(factor / hist_phi_jet->Integral());
	hist_e->Scale(factor / hist_e->Integral());
	hist_echad->Scale(factor / hist_echad->Integral());
	hist_ehad->Scale(factor / hist_ehad->Integral());
	hist_ecem->Scale(factor / hist_ecem->Integral());
	hist_eem->Scale(factor / hist_eem->Integral());
	hist_pt->Scale(factor / hist_pt->Integral());
	hist_pt_ratio_muon->Scale(factor / hist_pt_ratio_muon->Integral());
	hist_sv->Scale(factor / hist_sv->Integral());
	hist_svdeta->Scale(factor / hist_svdeta->Integral());
	hist_svdphi->Scale(factor / hist_svdphi->Integral());
	hist_svmass->Scale(factor / hist_svmass->Integral());
	hist_svntracks->Scale(factor / hist_svntracks->Integral());
	hist_svchi2->Scale(factor / hist_svchi2->Integral());
	hist_svndof->Scale(factor / hist_svndof->Integral());
	hist_svchi2div->Scale(factor / hist_svchi2div->Integral());
	hist_svdxy->Scale(factor / hist_svdxy->Integral());
	hist_svdxysig->Scale(factor / hist_svdxysig->Integral());
	hist_svd3d->Scale(factor / hist_svd3d->Integral());
	hist_svd3dsig->Scale(factor / hist_svd3dsig->Integral());
	hist_num->Scale(factor / hist_num->Integral());

	hist_frac_muon1->Scale(factor / hist_frac_muon1->Integral());
	hist_frac_muon_num1->Scale(factor / hist_frac_muon_num1->Integral());
	hist_eta_muon1->Scale(factor / hist_eta_muon1->Integral());
	hist_phi_muon1->Scale(factor / hist_phi_muon1->Integral());
	hist_pt_muon1->Scale(factor / hist_pt_muon1->Integral());
	hist_eta_jet1->Scale(factor / hist_eta_jet1->Integral());
	hist_phi_jet1->Scale(factor / hist_phi_jet1->Integral());
	hist_e1->Scale(factor / hist_e1->Integral());
	hist_echad1->Scale(factor / hist_echad1->Integral());
	hist_ehad1->Scale(factor / hist_ehad1->Integral());
	hist_ecem1->Scale(factor / hist_ecem1->Integral());
	hist_eem1->Scale(factor / hist_eem1->Integral());
	hist_pt1->Scale(factor / hist_pt1->Integral());
	hist_pt_ratio_muon1->Scale(factor / hist_pt_ratio_muon1->Integral());
	hist_sv1->Scale(factor / hist_sv1->Integral());
	hist_svdeta1->Scale(factor / hist_svdeta1->Integral());
	hist_svdphi1->Scale(factor / hist_svdphi1->Integral());
	hist_svmass1->Scale(factor / hist_svmass1->Integral());
	hist_svntracks1->Scale(factor / hist_svntracks1->Integral());
	hist_svchi21->Scale(factor / hist_svchi21->Integral());
	hist_svndof1->Scale(factor / hist_svndof1->Integral());
	hist_svchi2div1->Scale(factor / hist_svchi2div1->Integral());
	hist_svdxy1->Scale(factor / hist_svdxy1->Integral());
	hist_svdxysig1->Scale(factor / hist_svdxysig1->Integral());
	hist_svd3d1->Scale(factor / hist_svd3d1->Integral());
	hist_svd3dsig1->Scale(factor / hist_svd3dsig1->Integral());

	hist_frac_muon_num1->SetFillColor(3);
	hist_frac_muon1->SetFillColor(3);
	hist_sv1->SetFillColor(3);
	hist_svntracks1->SetFillColor(3);
	hist_svndof1->SetFillColor(3);

	hist_stack("stack_eta_muon", { hist_eta_muon1,hist_eta_muon });
	hist_stack("stack_phi_muon", { hist_phi_muon1,hist_phi_muon });
	hist_stack("stack_pt_muon", { hist_pt_muon1,hist_pt_muon });
	hist_stack("stack_eta_jet", { hist_eta_jet1,hist_eta_jet });
	hist_stack("stack_phi_jet", { hist_phi_jet1,hist_phi_jet });
	hist_stack("stack_frac_muon", { hist_frac_muon1,hist_frac_muon });
	hist_stack("stack_frac_muon_num", { hist_frac_muon_num1,hist_frac_muon_num });
	hist_stack("stack_e", { hist_e1,hist_e });
	hist_stack("stack_echad", { hist_echad1,hist_echad });
}
