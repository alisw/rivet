// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief CEP of h+h- (h=pi,K,p) at sqrt(s)=200 GeV with forward proton tagging
  class STAR_2020_I1792394 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(STAR_2020_I1792394);


    /// @name Analysis methods
    //@{
    
    enum CENTRAL_PARTICLES_PID {_PION, _KAON, _PROTON, _nAllowedPids};
    enum PARTICLE_CHARGE { _PLUS, _MINUS, _nSigns };
    enum PARTICLE_DIRECTION { _E, _W, _nBeamDirections }; // E = negative p_z, W = positive p_Z
    
    const double minPt[_nAllowedPids] = { 0.2*GeV, 0.3*GeV, 0.4*GeV };
    const double maxMinPt[_nAllowedPids] = { 9e9*GeV, 0.7*GeV, 1.1*GeV };

    /// Book histograms and initialise projections before the run
    void init() {

      // all final-state particles
      const FinalState fs( Cuts::NOCUT );
      declare(fs, "FS_all");
      
      // all final-state particles within STAR acceptance for this
      // measurement (reconstructed in the TPC and TOF)
      Cut centralCuts = Cuts::abscharge > 0 
                     && Cuts::abseta < 0.7 
                     && Cuts::pT > 0.2*GeV
                     && (  Cuts::abspid == PID::PIPLUS
                        || Cuts::abspid == PID::KPLUS
                        || Cuts::abspid == PID::PROTON );
      const FinalState fs_central( centralCuts );
      declare(fs_central, "FS_central");
      
      // forward-scattered beam particles detectable in Roman Pots
      // Checking the ID is not needed
      Cut forwardCuts = Cuts::abscharge > 0 
                     && Cuts::abseta > 5.0; // inclusive cut to select forward particles
      const FinalState fs_forward(forwardCuts);
      declare(fs_forward, "FS_forward");


      // Book histograms with binning taken from HEPdata
      book(_h["m_pipi"], "d01-x01-y01");            _scaleFactor["m_pipi"] = 1.0;
      book(_h["m_kk"], "d02-x01-y01");              _scaleFactor["m_kk"] = 1.0;
      book(_h["m_ppbar"], "d03-x01-y01");           _scaleFactor["m_ppbar"] = 1.0e3;
      book(_h["y_pipi"], "d04-x01-y01");            _scaleFactor["y_pipi"] = 1.0;
      book(_h["y_kk"], "d05-x01-y01");              _scaleFactor["y_kk"] = 1.0;
      book(_h["y_ppbar"], "d06-x01-y01");           _scaleFactor["y_ppbar"] = 1.0e3;
      book(_h["deltaPhi_pipi"], "d07-x01-y01");     _scaleFactor["deltaPhi_pipi"] = 1.0;
      book(_h["deltaPhi_kk"], "d08-x01-y01");       _scaleFactor["deltaPhi_kk"] = 1.0e3;
      book(_h["deltaPhi_ppbar"], "d09-x01-y01");    _scaleFactor["deltaPhi_ppbar"] = 1.0e3;
      book(_h["tSum_pipi"], "d10-x01-y01");         _scaleFactor["tSum_pipi"] = 1.0;
      book(_h["tSum_kk"], "d11-x01-y01");           _scaleFactor["tSum_kk"] = 1.0;
      book(_h["tSum_ppbar"], "d12-x01-y01");        _scaleFactor["tSum_ppbar"] = 1.0e3;

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
        
      // Retrieve all final-state particles
      const FinalState & fs = apply<FinalState>(event, "FS_all");
      // Veto event if number of particles in the final state is different from 4
      if(fs.size() != 4)
          return;
      
      // Retrieve accepted centrally produced particles
      const FinalState & fs_central = apply<FinalState>(event, "FS_central");
      // Veto event if number of centrally produced particles is different from 2
      if(fs_central.size() != 2)
          return;
      
      // Retrieve forward-scattered particles
      const FinalState & fs_forward = apply<FinalState>(event, "FS_forward");
      // Veto event if number of forward particles is different from 2
      if(fs_forward.size() != 2)
          return;

      // Continue checking forward particles (intact beam particles)
      // Storing forward particles in an array with cell ID indicating the direction (p_z)
      bool forwardParticlesInFiducialRegion[_nBeamDirections] = {false};
      Particle forwardParticles[_nBeamDirections];
      Vector3 forwardParticle2Vec[_nBeamDirections];
      for(const Particle & p : fs_forward.particles()){
          const int dir = p.pz() > 0 ? _W : _E;
          forwardParticle2Vec[dir] = Vector3(p.px(), p.py(), 0.0);
          forwardParticles[dir] = p;
          forwardParticlesInFiducialRegion[dir] = p.px() > -0.2 
                                                && fabs(p.py()) > 0.2
                                                && fabs(p.py()) < 0.4
                                                && (pow( p.px() + 0.3, 2) + pow( p.py(), 2 )) < 0.25;
      }
      if( !forwardParticlesInFiducialRegion[_E] || !forwardParticlesInFiducialRegion[_W] )
          return;
      
      // Storing central particles in an array with cell ID indicating the charge
      Particle csParticles[_nSigns];
      int totalCharge = 0;
      for(const Particle & p : fs_central.particles()){
          csParticles[ p.charge()>0 ? _PLUS : _MINUS] = p;
          totalCharge += p.charge();
      }
      // Checking the charge conservation, just in case
      if( totalCharge != 0 )
          return;
      
      // Determine PID of the central pair
      const int pid = ( csParticles[_PLUS].pid()==PID::PIPLUS && csParticles[_MINUS].pid()==PID::PIMINUS ) ? _PION :
                     (( csParticles[_PLUS].pid()==PID::KPLUS  && csParticles[_MINUS].pid()==PID::KMINUS) ? _KAON : 
                     (( csParticles[_PLUS].pid()==PID::PROTON && csParticles[_MINUS].pid()==PID::ANTIPROTON ) ? _PROTON : _nAllowedPids) );
      // skip event if particles in a pair are of different ID (should not happen)
      if(pid == _nAllowedPids)
          return;
      
      // Checking if central particles pass selection (in principle important for KK and ppbar)
      bool centralParticlesWithinFiducialRegion = true;
      for(int i=0; i<_nSigns; ++i)
          if( csParticles[i].pT() < minPt[pid] || min(csParticles[i].pT(), csParticles[1-i].pT()) > maxMinPt[pid] ){
              centralParticlesWithinFiducialRegion = false;
              break;
          }
      if( !centralParticlesWithinFiducialRegion )
          return;
      
      
      //-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
      // At this point event satisfies the definition of the fiducial region for events accepted
      // in the CEP measurement at STAR at 200 GeV
      
      const FourMomentum centralState4Mom = csParticles[_PLUS].momentum() + csParticles[_MINUS].momentum();      
      const double invMass = centralState4Mom.mass()/GeV;
      const double rapidity = centralState4Mom.rapidity();
      const double deltaPhi = fabs( forwardParticle2Vec[_W].angle( forwardParticle2Vec[_E] ) )/degree;
      
      // We need beam particles to get momentum transfers
      /* // Fragment below did not work for Pythia, unfortunately (four momenta were [0,0,0,0]); using a workaround
      Particle beamParticles[_nBeamDirections];
      const ParticlePair & beams = Rivet::Beam().beams();
      beamParticles[_W] = beams.first;
      beamParticles[_E] = beams.second;
      */
      // workaround - at this point we know that process is exclusive (2 forward protons + 2 central state particles)
      // assume that both beams are of the same type and collision in symmetric (lab frame = c.m.s. frame)
      const double sqrt_s = (centralState4Mom + forwardParticles[_E].momentum() + forwardParticles[_W].momentum()).mass();
      const double beamParticleMass = forwardParticles[_W].momentum().mass(); 
      const double fabsPz = sqrt( sqrt_s*sqrt_s/4. - beamParticleMass*beamParticleMass );
      FourMomentum beamParticles4Mom[_nBeamDirections];
      beamParticles4Mom[_W] = FourMomentum( sqrt(beamParticleMass*beamParticleMass + fabsPz*fabsPz), 0., 0., fabsPz );
      beamParticles4Mom[_E] = FourMomentum( sqrt(beamParticleMass*beamParticleMass + fabsPz*fabsPz), 0., 0., -fabsPz );
      // end of workaround

      double t[_nBeamDirections];
      for(int dir=0; dir<_nBeamDirections; ++dir)
          t[dir] = (beamParticles4Mom[dir] - forwardParticles[dir].momentum()).mass2()/(GeV*GeV);
      const double tSum = fabs(t[_E] + t[_W]);
      
      if( pid==_PION ){
        _h["m_pipi"]->fill( invMass );
        _h["y_pipi"]->fill( rapidity );
        _h["deltaPhi_pipi"]->fill( deltaPhi );
        _h["tSum_pipi"]->fill( tSum );
      } else if( pid==_KAON ){
        _h["m_kk"]->fill( invMass );
        _h["y_kk"]->fill( rapidity );
        _h["deltaPhi_kk"]->fill( deltaPhi );
        _h["tSum_kk"]->fill( tSum );
      } else{
        _h["m_ppbar"]->fill( invMass );
        _h["y_ppbar"]->fill( rapidity );
        _h["deltaPhi_ppbar"]->fill( deltaPhi );
        _h["tSum_ppbar"]->fill( tSum );
      }


    }


    /// Normalise histograms etc., after the run
    void finalize() {

      const double scalingFactor = crossSection()/nanobarn/sumOfWeights();
      // scale to cross section
      for(auto &hist : _h)
        scale(hist.second, scalingFactor*_scaleFactor[hist.first]);
    }

    //@}


    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    map<string, double> _scaleFactor; // map with scale factors to ensure cross section units in agreement with HEPdata
    //@}


  };


  DECLARE_RIVET_PLUGIN(STAR_2020_I1792394);

}
