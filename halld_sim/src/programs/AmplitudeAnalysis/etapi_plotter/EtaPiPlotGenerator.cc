#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "EtaPiPlotGenerator.h"
#include "IUAmpTools/Histogram1D.h"
#include "IUAmpTools/Kinematics.h"

EtaPiPlotGenerator::EtaPiPlotGenerator( const FitResults& results ) :
PlotGenerator( results )
{
	createHistograms();
}

EtaPiPlotGenerator::EtaPiPlotGenerator( ) :
PlotGenerator( )
{
	createHistograms();
}

void EtaPiPlotGenerator::createHistograms() {
  // calls to bookHistogram go here
  
  bookHistogram( kEtaPiMass, new Histogram1D( 86, 0.28, 2.0, "Metapi", "Invariant Mass of #eta #pi") );
  bookHistogram( kEtaCosTheta, new Histogram1D( 50, -1., 1., "cosTheta", "cos( #theta ) of Resonance Production") );

  bookHistogram( kPhiPi,  new Histogram1D( 50, -1*PI, PI, "PhiPi",  "#Phi_{#pi}" ) );
  bookHistogram( kPhiEta, new Histogram1D( 50, -1*PI, PI, "PhiPiEta", "#Phi_{#eta}" ) );
  bookHistogram( kPhi, new Histogram1D( 50, -1*PI, PI, "Phi", "#Phi" ) );
  bookHistogram( kphi, new Histogram1D( 50, -1*PI, PI, "phi", "#phi" ) );
  bookHistogram( kPsi, new Histogram1D( 50, -1*PI, PI, "psi", "#psi" ) );
  bookHistogram( kt, new Histogram1D( 100, 0, 2.00, "t", "-t" ) );
}

void
EtaPiPlotGenerator::projectEvent( Kinematics* kin, const string& reactionName ){
  //cout << cfgInfo()->amplitudeList( reactionName, "", "").at(0)->factors().at(0).at(5) << endl;
  double polAngle=stod(cfgInfo()->amplitudeList( reactionName, "", "").at(0)->factors().at(0).at(5));

  // Our reaction has the order 14 7 17 when using tree_to_amptools. so its proton pi0 eta 
  TLorentzVector beam   = kin->particle( 0 );
  TLorentzVector recoil = kin->particle( 1 );
  TLorentzVector p1 = kin->particle( 2 );
  TLorentzVector p2 = kin->particle( 3 );

  TLorentzVector resonance = p1 + p2; 
  TLorentzRotation resonanceBoost( -resonance.BoostVector() );

  TLorentzVector recoil_res = resonanceBoost * recoil;
  TLorentzVector p1_res = resonanceBoost * p1;
  TLorentzVector p2_res = resonanceBoost * p2;

  // normal to the production plane
  TVector3 y = (beam.Vect().Unit().Cross(-recoil.Vect().Unit())).Unit();
  
  //// choose helicity frame: z-axis opposite recoil proton in rho rest frame
  TVector3 z = -1. * recoil_res.Vect().Unit();
  TVector3 x = y.Cross(z).Unit();
  TVector3 angles(   (p2_res.Vect()).Dot(x),
                     (p2_res.Vect()).Dot(y),
                     (p2_res.Vect()).Dot(z) );
  

  // choose GJ frame (is this correct form?)
  // pg 78 of https://arxiv.org/pdf/1310.7498.pdf
  // The Gottfried-Jackson frame (GJ) is a frame where the resonance (X) is at rest. z
  // is in the direction of the beam and y is perpendicular to the production plane
  //TVector3 z = beam.Vect().Unit();
  //TVector3 x = y.Cross(z).Unit();
  //TVector3 angles(   (p2_res.Vect()).Dot(x),
  //                   (p2_res.Vect()).Dot(y),
  //                   (p2_res.Vect()).Dot(z) );

  GDouble cosTheta = angles.CosTheta();
  
  GDouble phi = angles.Phi();
  
  TVector3 eps(TMath::Cos(polAngle), TMath::Sin(polAngle), 0.0); // beam polarization vector
  GDouble Phi = atan2(y.Dot(eps), beam.Vect().Unit().Dot(eps.Cross(y)));

  GDouble psi = phi - Phi;
  if(psi < -1*PI) psi += 2*PI;
  if(psi > PI) psi -= 2*PI;

  // compute invariant t
  GDouble t = - 2* recoil.M() * (recoil.E()-recoil.M());

  // calls to fillHistogram go here
  
  fillHistogram( kEtaPiMass, ( resonance ).M() );
  
  fillHistogram( kEtaCosTheta, cosTheta );

  fillHistogram( kPhiPi,  p1.Phi() );
  fillHistogram( kPhiEta, p2.Phi() );
  fillHistogram( kPhi, Phi );
  fillHistogram( kphi, phi );

  fillHistogram( kPsi, psi );
  fillHistogram( kt, -t );      // fill with -t to make positive
}
