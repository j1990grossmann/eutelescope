#ifndef EUtelProcessorTrueHitAnalysis_H
#define EUTelProcessorTrueHitAnalysis_H 1

#include "EUTELESCOPE.h"

#include "marlin/Processor.h"

#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>

#include <IMPL/LCCollectionVec.h>

#include <array>
#include <cmath>
#include <map>
#include <string>
#include <vector>

namespace eutelescope {

  class EUTelProcessorTrueHitAnalysis : public marlin::Processor {

  public:
    virtual Processor *newProcessor() {
      return new EUTelProcessorTrueHitAnalysis;
    }

    virtual const std::string &name() const { return Processor::name(); }

    EUTelProcessorTrueHitAnalysis();

    virtual void init();

    virtual void processRunHeader(LCRunHeader *rhdr);

    virtual void processEvent(LCEvent *evt);

    virtual void check(LCEvent *evt);

    virtual void end();

    int findPairIndex(double const *a, std::vector<double const *> vec);

    void bookHistos();
    void fillHistos(LCEvent *event);

  protected:
    std::string _trueHitCollectionName;
    std::string _reconstructedHitCollectionName;

    int _iRun;
    int _iEvt;

  private:
    std::vector<int> _sensorIDVec;

    LCCollectionVec *_trueHitCollectionVec;
    LCCollectionVec *_reconstructedHitCollectionVec;

    std::array<std::map<int, AIDA::IHistogram1D *>, 8> _1DHistos;
    std::array<std::map<int, AIDA::IHistogram2D *>, 8> _2DHistos;
    std::vector<AIDA::IHistogram1D *> _xClustSize2Histos;
    std::vector<AIDA::IHistogram1D *> _yClustSize3Histos;

    void readCollections(LCEvent *event);
  };

  EUTelProcessorTrueHitAnalysis gEUTelProcessorTrueHitAnalysis;
}

#endif
