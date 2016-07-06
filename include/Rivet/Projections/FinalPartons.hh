// -*- C++ -*-
#ifndef RIVET_FinalPartons_HH
#define RIVET_FinalPartons_HH

#include "Rivet/Projections/FinalState.hh"

namespace Rivet {

    class FinalPartons : public FinalState {
        public:
            FinalPartons(Cut c=Cuts::open()) : FinalState(c) { }

            const Projection* clone() const {
                return new FinalPartons(*this);
            }

            void project(const Event& e);

        protected:
            bool accept(const Particle& p) const;
    };

}

#endif
