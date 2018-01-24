// Author Johannes Grossmann, HEPHY<mailto:Johannes.Grossmann@oeaw.ac.at>
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h"
#include "EUTelGeometryTelescopeGeoDescription.h"

#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelMAPSAHitMaker.h"
#include "EUTelRunHeaderImpl.h"

#include "EUTelBrickedClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelGenericSparseClusterImpl.h"
#include "EUTelGeometricClusterImpl.h"
#include "EUTelSimpleVirtualCluster.h"
#include "EUTelSparseClusterImpl.h"

#include "EUTelAlignmentConstant.h"
#include "EUTelExceptions.h"
#include "EUTelReferenceHit.h"

// marlin includes ".h"
#include "marlin/Global.h"
#include "marlin/Processor.h"

// aida includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ITree.h>
#include <marlin/AIDAProcessor.h>
#endif

// ROOT includes
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include <TVector3.h>
#endif


// lcio includes <.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackerPulseImpl.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTime.h>
#include <UTIL/LCTOOLS.h>

// system includes <>
#include <algorithm>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace std;
using namespace marlin;
using namespace eutelescope;

// definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string EUTelMAPSAHitMaker::_hitHistoLocalName       = "HitHistoLocal";
std::string EUTelMAPSAHitMaker::_hitHistoTelescopeName   = "HitHistoTelescope";
std::string EUTelMAPSAHitMaker::_densityPlotName         = "DensityPlot";
std::string EUTelMAPSAHitMaker::_clusterCenterHistoName  = "ClusterCenter";
std::string EUTelMAPSAHitMaker::_clusterCenterXHistoName = "ClusterCenterX";
std::string EUTelMAPSAHitMaker::_clusterCenterYHistoName = "ClusterCenterY";
#endif

EUTelMAPSAHitMaker::EUTelMAPSAHitMaker()
: Processor("EUTelMAPSAHitMaker"), _inputhitCollectionName(),
_hitCollectionName(), _referenceHitCollectionName("referenceHit"),
_referenceHitCollectionVec(), _wantLocalCoordinates(false),
_referenceHitLCIOFile("reference.slcio"), _iRun(0), _neg_shift(0), _iEvt(0), _total_evts_cms_ref_file(0),
_alreadyBookedSensorID(), _shift_pixel(), _aidaHistoMap(), _histogramSwitch(true) {
    // modify processor description
    _description = "EUTelMAPSAHitMaker is to Merge the LCIO TRACKERHIT Collections"
    "From the telescop with CMSPixel and MPA\nto the external "
    "frame of reference using the GEAR geometry description";
    
    //     registerInputCollection(LCIO::TRACKERPULSE, "PulseCollectionName",
    //                             "Input cluster collection name", _inputhitCollectionName,
    //                             std::string(""));
    registerInputCollection(LCIO::TRACKERHIT, "InputHitCollection",
                            "Input cluster collection name", _inputhitCollectionName,
                            std::string(""));
    
    registerOutputCollection(LCIO::TRACKERHIT, "HitCollectionName",
                             "Output hit collection name", _hitCollectionName,
                             std::string(""));
    registerProcessorParameter("MAPSA_Cluster_Root_File",
                               "Provide the root file with the clustered MAPSA Data",
                               _RootFilename,
                               std::string(""));
    registerProcessorParameter("MAPSA_Cluster_Shift_File",
                               "Provide the root file the Shiftvectors",
                               _ShiftFilename,
                               std::string(""));
    registerProcessorParameter("DataPath",
                               "Provide the path, where the root file resides",
                               _DataPath,
                               std::string(""));
    registerOptionalParameter(
        "EnableLocalCoordidates",
        "Hit coordinates are calculated in local reference frame of sensor",
        _wantLocalCoordinates, static_cast<bool>(false));
    
    registerOptionalParameter(
        "ReferenceCollection",
        "This is the name of the reference hit collection initialized in this "
        "processor. This collection provides the reference vector to correctly "
        "determine a plane corresponding to a global hit coordiante.",
        _referenceHitCollectionName, static_cast<string>("referenceHit"));
    
    registerOptionalParameter(
        "ReferenceHitFile",
        "This is the file where the reference hit collection is stored",
        _referenceHitLCIOFile, std::string("reference.slcio"));
}

void EUTelMAPSAHitMaker::init() {
    printParameters();
    if(_DataPath!="") {
        streamlog_out(MESSAGE0)<<"Data path = "<<_DataPath<<"\n";
    } else {
        streamlog_out(ERROR)<<"Data path not defined\n";
    }
    boost::filesystem::path dir (_DataPath);
    boost::filesystem::path file (_RootFilename);
    std::string file_shift=_RootFilename.substr(0,_RootFilename.length()-5)+"_jumps.root";
    boost::filesystem::path file1 (file_shift);
    boost::filesystem::path p = dir / file;
    boost::filesystem::path p1 = dir / file1;
    
    //   boost::filesystem::is_directory(data_dir);
    this->checkFile(p);
    this->checkFile(p1);
    _input_file = new TFile(p.string().c_str(),"READ");
    if(_input_file ->TObject::IsZombie()) {
        streamlog_out(ERROR)<<"Error opening file "<<p.string().c_str()<<std::endl;
        exit(1);
    }
    else
        streamlog_out(MESSAGE0)<<"Process file\t"<<p.string().c_str()<<std::endl;
    _reader = new TTreeReader("Clustertree", _input_file);
    _f_counter_cluster = new TTreeReaderValue<std::vector<CounterCluster>> (*_reader, "CounterCluster");
    // set to zero the run and event counters
    _iRun = 0;
    
    _input_file_corr = new TFile(p1.string().c_str(),"READ");
    if(_input_file_corr ->TObject::IsZombie()) {
        streamlog_out(ERROR)<<"Error opening file "<<p1.string().c_str()<<std::endl;
        exit(1);
    }
    else
        streamlog_out(MESSAGE0)<<"Open correlation_file\t"<<p1.string().c_str()<<std::endl;
    TTreeReader myReader("tree", _input_file_corr);
    TTreeReaderValue<std::vector<int>> shift_mpa(myReader, "shift_vec_mpa");
    TTreeReaderValue<std::vector<int>> shift_tel(myReader, "shift_vec_tel");
    myReader.Next();
    _shift_mpa   = static_cast<std::vector<int>>(*shift_mpa);
    _shift_pixel = static_cast<std::vector<int>>(*shift_tel);
    
    _goodflaghist = dynamic_cast<TH1S*>(_input_file_corr->Get("Sumflag"));
    _goodflaghist->SetDirectory(0);
    _input_file_corr->Close();
    
    
    geo::gGeometry().initializeTGeoDescription(EUTELESCOPE::GEOFILENAME,
                                               EUTELESCOPE::DUMPGEOROOT);
    
    _histogramSwitch = true;
    
    // only for global coord we need a refhit collection
    if (!_wantLocalCoordinates) {
        DumpReferenceHitDB();
    }
    _cms_ref_syncreader = std::unique_ptr<IO::LCReader>(IOIMPL::LCFactory::getInstance()->createLCReader());
    //     _cms_ref_syncreader->setReadCollectionNames("cms_ref_hit");
    StringVec lcioInputFiles_1 ; 
    Global::parameters->getStringVals("LCIOInputFiles" , lcioInputFiles_1 );
    
    for(const auto& file_lcio: lcioInputFiles_1){
        streamlog_out(MESSAGE0)<<"Processing the files "<<file_lcio<<"\n";
    }
    _cms_ref_syncreader->open(lcioInputFiles_1);
    _total_evts_cms_ref_file = _cms_ref_syncreader->getNumberOfEvents();
    _timer = new boost::timer::auto_cpu_timer();
}

void EUTelMAPSAHitMaker::DumpReferenceHitDB() {
    // create a reference hit collection file (DB)
    LCWriter *lcWriter = LCFactory::getInstance()->createLCWriter();
    try {
        lcWriter->open(_referenceHitLCIOFile, LCIO::WRITE_NEW);
    } catch (IOException &e) {
        streamlog_out(ERROR4) << e.what() << endl;
        exit(-1);
    }    
    streamlog_out(MESSAGE5) << "Writing to " << _referenceHitLCIOFile
    << std::endl;
    
    LCRunHeaderImpl *lcHeader = new LCRunHeaderImpl;
    lcHeader->setRunNumber(0);
    lcWriter->writeRunHeader(lcHeader);
    delete lcHeader;
    LCEventImpl *event = new LCEventImpl;
    event->setRunNumber(0);
    event->setEventNumber(0);
    LCTime *now = new LCTime;
    event->setTimeStamp(now->timeStamp());
    delete now;
    
    LCCollectionVec *referenceHitCollection =
    new LCCollectionVec(LCIO::LCGENERICOBJECT);
    double refVec[3];
    
    EVENT::IntVec sensorIDVec = geo::gGeometry().sensorIDsVec();
    
    for (EVENT::IntVec::iterator it = sensorIDVec.begin();
         it != sensorIDVec.end(); it++) {
        EUTelReferenceHit *refhit = new EUTelReferenceHit();
    
    int sensorID = *it;
    refhit->setSensorID(sensorID);
    refhit->setXOffset(geo::gGeometry().siPlaneXPosition(sensorID));
    refhit->setYOffset(geo::gGeometry().siPlaneYPosition(sensorID));
    refhit->setZOffset(geo::gGeometry().siPlaneZPosition(sensorID) +
    0.5 * geo::gGeometry().siPlaneZSize(sensorID));
    
    refVec[0] = 0.;
    refVec[1] = 0.;
    refVec[2] = 1.;
    
    double gRotation[3] = {0., 0., 0.};                         // not rotated
    gRotation[0] = geo::gGeometry().siPlaneZRotation(sensorID); // Euler alpha ;
    gRotation[1] = geo::gGeometry().siPlaneYRotation(sensorID); // Euler alpha ;
    gRotation[2] = geo::gGeometry().siPlaneXRotation(sensorID); // Euler alpha ;
    streamlog_out(DEBUG5) << "GEAR rotations: " << gRotation[0] << " "
    << gRotation[1] << " " << gRotation[2] << endl;
    gRotation[0] = gRotation[0] * 3.1415926 / 180.; //
    gRotation[1] = gRotation[1] * 3.1415926 / 180.; //
    gRotation[2] = gRotation[2] * 3.1415926 / 180.; //
    
    TVector3 _RotatedVector(refVec[0], refVec[1], refVec[2]);
    TVector3 _Xaxis(1.0, 0.0, 0.0);
    TVector3 _Yaxis(0.0, 1.0, 0.0);
    TVector3 _Zaxis(0.0, 0.0, 1.0);
    
    if (TMath::Abs(gRotation[2]) > 1e-6) {
        _RotatedVector.Rotate(gRotation[2], _Xaxis); // in ZY
    }
    if (TMath::Abs(gRotation[1]) > 1e-6) {
        _RotatedVector.Rotate(gRotation[1], _Yaxis); // in ZX
    }
    if (TMath::Abs(gRotation[0]) > 1e-6) {
        _RotatedVector.Rotate(gRotation[0], _Zaxis); // in XY
    }
    
    refhit->setAlpha(_RotatedVector[0]);
    refhit->setBeta(_RotatedVector[1]);
    refhit->setGamma(_RotatedVector[2]);
    referenceHitCollection->push_back(refhit);
         }
         event->addCollection(referenceHitCollection, _referenceHitCollectionName);
         
         lcWriter->writeEvent(event);
         delete event;
         lcWriter->close();
}

void EUTelMAPSAHitMaker::processRunHeader(LCRunHeader *rdr) {
    std::unique_ptr<EUTelRunHeaderImpl> header =
    std::make_unique<EUTelRunHeaderImpl>(rdr);
    header->addProcessor(type());
    
    // this is the right place also to check the geometry ID. This is a
    // unique number identifying each different geometry used at the
    // beam test. The same number should be saved in the run header and
    // in the xml file. If the numbers are different, instead of barely
    // quitting ask the user what to do.
    
    if (header->getGeoID() == 0)
        streamlog_out(WARNING0)
        << "The geometry ID in the run header is set to zero." << endl
        << "This may mean that the GeoID parameter was not set" << endl;
    
    // increment the run counter
    ++_iRun;
}

// void EUTelMAPSAHitMaker::readDataSource(int numEvents) {
void EUTelMAPSAHitMaker::processEvent(LCEvent *event) {
    _iEvt++;
    int runnumber = event->getRunNumber();
    int evtno_tel = event->getEventNumber();
    //  Merger flags and shifted evt nos
    bool goodflag = false;
    int mapsa_evt_no   = 0;
    int cms_ref_evt_no = 0;    
    // Initialize Old detector ID for histogramming and geometry constants
    int oldDetectorID = -100;    
    double xSize = 0., ySize = 0.;
    double resolutionX = 0., resolutionY = 0.;
    
    static const std::string mcpName_ref("cms_ref_hit");
    // This functions determines the correct shift in the datastreams.
    this->getShift(evtno_tel,goodflag,mapsa_evt_no,cms_ref_evt_no);
    if((_iEvt<100 && _iEvt%10==0 ) || (_iEvt%1000==0 && evtno_tel!=cms_ref_evt_no)){
        streamlog_out(MESSAGE0)<<"Tel event\tCMS_REF_event\tMapsaEvent\n"
        <<evtno_tel<<"\t"<<cms_ref_evt_no<<"\t"<<mapsa_evt_no<<"\n";
    }
    //     Prepare the output hit collection
    LCCollectionVec *hitCollection = nullptr;
    try {
        hitCollection = static_cast<LCCollectionVec *>(
            event->getCollection(_hitCollectionName));
    } catch (...) {
        hitCollection = new LCCollectionVec(LCIO::TRACKERHIT);
    }
    
    //  Input event from Marlin
    EUTelEventImpl *evt = static_cast<EUTelEventImpl *>(event);    
    if (evt->getEventType() == kEORE) {
        streamlog_out(DEBUG4) << "EORE found: nothing else to do." << endl;
        return;
    } else if (evt->getEventType() == kUNKNOWN) {
        streamlog_out(WARNING2) << "Event number " << evt->getEventNumber()
        << " in run " << evt->getRunNumber()
        << " is of unknown type. Continue considering it "
        "as a normal Data Event."
        << endl;
    }   
    //     Get the input hit collection
    //     Prepare the input hit collections
    const LCCollectionVec *input_hitCollection = nullptr;
    try {
        input_hitCollection = static_cast<LCCollectionVec *>(
            event->getCollection(_inputhitCollectionName));
    } catch (...) {
        streamlog_out(MESSAGE2) << "No input collection " << _inputhitCollectionName
        << " found on event " << event->getEventNumber()
        << " in run " << event->getRunNumber() << endl;
        return;
    }
    //     Decoder and encoder for SensorIDs
    UTIL::CellIDDecoder<TrackerHitImpl> hitCellDecoder(input_hitCollection);
    UTIL::CellIDEncoder<TrackerHitImpl> idHitEncoder(EUTELESCOPE::HITENCODING,hitCollection);
    //  Loop to copy the telescope hits in output hit collection
    for (int iHit = 0; iHit < input_hitCollection->getNumberOfElements();iHit++){
        auto pulseguard = std::make_unique<TrackerHitImpl>(*dynamic_cast<TrackerHitImpl *>(input_hitCollection->getElementAt(iHit)));
        int sensorID = hitCellDecoder(pulseguard.get())["sensorID"];
        auto globalhitprop = static_cast<int>(hitCellDecoder(pulseguard.get())["properties"]);
        idHitEncoder["sensorID"] = sensorID;
        //         Set the flag for local global coordinates
        idHitEncoder["properties"] = 0; // init
        if (!_wantLocalCoordinates)
            idHitEncoder["properties"] = kHitInGlobalCoord;
        // store values
        //         idHitEncoder.setCellID(pulseguard.get());
        hitCollection->emplace_back(pulseguard.get());
        auto * const pulse = pulseguard.release();
//         pulse->get
        if (sensorID != oldDetectorID) {
            oldDetectorID = sensorID;
            // check if the histos for this sensor ID have been booked already.
            if (_alreadyBookedSensorID.find(sensorID) ==
                _alreadyBookedSensorID.end()) {
                bookHistos(sensorID);
                }
                resolutionX = geo::gGeometry().siPlaneXResolution(sensorID); // mm
                resolutionY = geo::gGeometry().siPlaneYResolution(sensorID); // mm
                
                xSize = geo::gGeometry().siPlaneXSize(sensorID); // mm
                ySize = geo::gGeometry().siPlaneYSize(sensorID); // mm
                
        }
        // LOCAL coordinate system !!!!!!
        double telPos[3];
        const double* hitpos = pulse->getPosition();
        telPos[0] = hitpos[0];telPos[1] = hitpos[1];telPos[2] = hitpos[2];
        // We now plot the the hits in the EUTelescope local frame. This frame has the
        // coordinate centre at the sensor centre.
        #if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
        string tempHistoName;
        if (_histogramSwitch) {
            tempHistoName = _hitHistoLocalName + "_" + to_string(sensorID);
            if (AIDA::IHistogram2D *histo = dynamic_cast<AIDA::IHistogram2D *>(
                _aidaHistoMap[tempHistoName])) {
                histo->fill(hitpos[0], hitpos[1]);
                } else {
                    streamlog_out(ERROR1)
                    << "Not able to retrieve histogram pointer for " << tempHistoName
                    << ".\nDisabling histogramming from now on " << endl;
                    _histogramSwitch = false;
                }
        }
        #endif
        if (!_wantLocalCoordinates) {
            //
            // NOW !!
            // GLOBAL coordinate system !!!
            
            const double localPos[3] = {hitpos[0], hitpos[1], hitpos[2]};
            geo::gGeometry().local2Master(sensorID, localPos, telPos);
        }
        #if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
        if (_histogramSwitch) {
            tempHistoName = _hitHistoTelescopeName + "_" + to_string(sensorID);
            AIDA::IHistogram2D *histo2D =
            dynamic_cast<AIDA::IHistogram2D *>(_aidaHistoMap[tempHistoName]);
            if (histo2D) {
                histo2D->fill(hitpos[0], hitpos[1]);
            } else {
                streamlog_out(ERROR1)
                << "Not able to retrieve histogram pointer for " << tempHistoName
                << ".\nDisabling histogramming from now on " << endl;
                _histogramSwitch = false;
            }
        }
        #endif
        
    }
    //     End of the TelescopeHit Loop
    //     Now add the hits from the other subsystems
    if(goodflag==true){
        const LCCollection*    col_ref = nullptr;
        const LCEvent*    tmpevt = static_cast<LCEvent*>(_cms_ref_syncreader->readEvent(runnumber,cms_ref_evt_no));
        try {
            col_ref = static_cast<LCCollection*>(tmpevt->getCollection(mcpName_ref));
            UTIL::CellIDDecoder<TrackerHitImpl> hitCellDecoder_ref(col_ref);
            
            int nMcp = col_ref->getNumberOfElements();
            for( int j1 = 0 ; j1 < nMcp ; j1++ ) {
                //                 TrackerHitImpl* cms_refhit = dynamic_cast<TrackerHitImpl *>(col_ref->getElementAt(j1));
                auto cms_refhit = std::make_unique<TrackerHitImpl>(*dynamic_cast<TrackerHitImpl *>(col_ref->getElementAt(j1)));
                int sensorID = hitCellDecoder_ref(cms_refhit.get())["sensorID"];
                idHitEncoder["sensorID"] = sensorID;
                idHitEncoder["properties"] = 0; // init
                if (!_wantLocalCoordinates)
                    idHitEncoder["properties"] = kHitInGlobalCoord;
                // store values
                //         idHitEncoder.setCellID(pulse);        
                hitCollection->emplace_back(cms_refhit.get());
                auto* const pulse = cms_refhit.release();
                if (sensorID != oldDetectorID) {
                    oldDetectorID = sensorID;
                    // check if the histos for this sensor ID have been booked already.
                    if (_alreadyBookedSensorID.find(sensorID) ==
                        _alreadyBookedSensorID.end()) {
                        bookHistos(sensorID);
                        }
                        resolutionX = geo::gGeometry().siPlaneXResolution(sensorID); // mm
                        resolutionY = geo::gGeometry().siPlaneYResolution(sensorID); // mm
                        
                        xSize = geo::gGeometry().siPlaneXSize(sensorID); // mm
                        ySize = geo::gGeometry().siPlaneYSize(sensorID); // mm
                        
                }
                // LOCAL coordinate system !!!!!!
                double telPos[3];
                const double* hitpos = pulse->getPosition();
                telPos[0] = hitpos[0];telPos[1] = hitpos[1];telPos[2] = hitpos[2];
                // We now plot the the hits in the EUTelescope local frame. This frame has the
                // coordinate centre at the sensor centre.
                #if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
                string tempHistoName;
                if (_histogramSwitch) {
                    tempHistoName = _hitHistoLocalName + "_" + to_string(sensorID);
                    if (AIDA::IHistogram2D *histo = dynamic_cast<AIDA::IHistogram2D *>(
                        _aidaHistoMap[tempHistoName])) {
                        histo->fill(hitpos[0], hitpos[1]);
                        } else {
                            streamlog_out(ERROR1)
                            << "Not able to retrieve histogram pointer for " << tempHistoName
                            << ".\nDisabling histogramming from now on " << endl;
                            _histogramSwitch = false;
                        }
                }
                #endif
                if (!_wantLocalCoordinates) {
                    //
                    // NOW !!
                    // GLOBAL coordinate system !!!
                    
                    const double localPos[3] = {hitpos[0], hitpos[1], hitpos[2]};
                    geo::gGeometry().local2Master(sensorID, localPos, telPos);
                }
                #if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
                if (_histogramSwitch) {
                    tempHistoName = _hitHistoTelescopeName + "_" + to_string(sensorID);
                    AIDA::IHistogram2D *histo2D =
                    dynamic_cast<AIDA::IHistogram2D *>(_aidaHistoMap[tempHistoName]);
                    if (histo2D) {
                        histo2D->fill(hitpos[0], hitpos[1]);
                    } else {
                        streamlog_out(ERROR1)
                        << "Not able to retrieve histogram pointer for " << tempHistoName
                        << ".\nDisabling histogramming from now on " << endl;
                        _histogramSwitch = false;
                    }
                }
                #endif                
            }
        } catch(const std::exception& e){
            //             streamlog_out(MESSAGE0)<<e.what()<<std::endl;
        }
        //     }
        
        //     if(mapsa_evt_no>-1 && mapsa_evt_no<_reader->GetEntries(true)) {
        _reader->SetEntry(mapsa_evt_no);
        if(_reader->GetCurrentEntry()==1) {            
            streamlog_out(MESSAGE0)<<"Check mapsa event "<<_reader->GetCurrentEntry()<<" for valid MAPSA Data\n";
            if (!CheckValue(*_f_counter_cluster)) {
                streamlog_out(ERROR)<<"ill defined data\n";
            } else {
                streamlog_out(MESSAGE0)<<"Setup Successfully\n";
            }
        }
        //         if( _reader->GetCurrentEntry() <10){
        //             streamlog_out(DEBUG5)<<"TTreereader    Entry"<<setw(12)<<std::fixed<<std::right<<_reader->GetCurrentEntry()<<"\n";
        //             streamlog_out(DEBUG5)<<"Event Tel      Entry"<<setw(12)<<std::fixed<<std::right<<evt->getEventNumber()<<"\n";
        //             streamlog_out(DEBUG5)<<"Event CMSREF   Entry"<<setw(12)<<std::fixed<<std::right<<cms_ref_evt_no<<"\n";
        //         }
        for(const auto &cluster: (**_f_counter_cluster)) {
            static constexpr int sensorID = 101;
            if (sensorID != oldDetectorID) {
                oldDetectorID = sensorID;
                // check if the histos for this sensor ID have been booked already.
                if (_alreadyBookedSensorID.find(sensorID) ==
                    _alreadyBookedSensorID.end()) {
                    bookHistos(sensorID);
                    }
                    
                    resolutionX = geo::gGeometry().siPlaneXResolution(sensorID); // mm
                    resolutionY = geo::gGeometry().siPlaneYResolution(sensorID); // mm
                    
                    xSize = geo::gGeometry().siPlaneXSize(sensorID); // mm
                    ySize = geo::gGeometry().siPlaneYSize(sensorID); // mm
                    
            }
            
            //             double telPos[3];            
            double hitpos[3];            
            hitpos[0] = cluster.cog_x/1000. - xSize / 2.;
            hitpos[1] = -cluster.cog_y/1000. + ySize / 2.;
            hitpos[2] = 0.;
            // We now plot the the hits in the EUTelescope local frame. This frame has the
            // coordinate centre at the sensor centre.
            #if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
            string tempHistoName;
            if (_histogramSwitch) {
                tempHistoName = _hitHistoLocalName + "_" + to_string(sensorID);
                if (AIDA::IHistogram2D *histo = dynamic_cast<AIDA::IHistogram2D *>(
                    _aidaHistoMap[tempHistoName])) {
                    histo->fill(hitpos[0], hitpos[1]);
                    } else {
                        streamlog_out(ERROR1)
                        << "Not able to retrieve histogram pointer for " << tempHistoName
                        << ".\nDisabling histogramming from now on " << endl;
                        _histogramSwitch = false;
                    }
            }
            #endif
            
            if (!_wantLocalCoordinates) {
                //
                // NOW !!
                // GLOBAL coordinate system !!!
                
                const double localPos[3] = {hitpos[0], hitpos[1], hitpos[2]};
                //                 geo::gGeometry().local2Master(sensorID, localPos, telPos);
                geo::gGeometry().local2Master(sensorID, localPos, hitpos);
            }
            
            #if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
            if (_histogramSwitch) {
                tempHistoName = _hitHistoTelescopeName + "_" + to_string(sensorID);
                AIDA::IHistogram2D *histo2D =
                dynamic_cast<AIDA::IHistogram2D *>(_aidaHistoMap[tempHistoName]);
                if (histo2D) {
                    histo2D->fill(hitpos[0], hitpos[1]);
                } else {
                    streamlog_out(ERROR1)
                    << "Not able to retrieve histogram pointer for " << tempHistoName
                    << ".\nDisabling histogramming from now on " << endl;
                    _histogramSwitch = false;
                }
            }
            #endif
            // create the new hit for the MPA
            TrackerHitImpl *hit = new TrackerHitImpl;
            hit->setPosition(&hitpos[0]);
            float cov[TRKHITNCOVMATRIX] = {0., 0., 0., 0., 0., 0.};
            cov[0] = resolutionX * resolutionX; // cov(x,x)
            cov[2] = resolutionY * resolutionY; // cov(y,y)
            hit->setCovMatrix(cov);
            hit->setType(4);
            hit->setTime(0);
//             Abuse quality flag to store the mapsa event number
            hit->setQuality(mapsa_evt_no);
            idHitEncoder["sensorID"] = sensorID;
            idHitEncoder["properties"] = 0; // init
            if (!_wantLocalCoordinates)
                idHitEncoder["properties"] = kHitInGlobalCoord;
            idHitEncoder.setCellID(hit);
            // add the new hit to the hit collection
            hitCollection->emplace_back(hit);
        }
    }
    try {
        event->getCollection(_hitCollectionName);
    } catch (...) {
        event->addCollection(hitCollection, _hitCollectionName);
    }
    if(_iEvt<5){
        streamlog_out(MESSAGE0)<<"Hits input hit collection\n";
        //         LCTOOLS::printTrackerHits(input_hitCollection);
        streamlog_out(MESSAGE0)<<"Hits output hit collection\n";
        //         LCTOOLS::printTrackerHits(hitCollection);
    }
    if (isFirstEvent())
        _isFirstEvent = false;
}

void EUTelMAPSAHitMaker::end() {
    _goodflaghist->Delete();
    delete _f_counter_cluster;
    delete _reader;
    _input_file->Close();
    _input_file_corr->Close();
    delete _timer;
    _cms_ref_syncreader->close();
    streamlog_out(MESSAGE4) << "Successfully finished" << endl;
}

void EUTelMAPSAHitMaker::bookHistos(int sensorID) {
    
    #if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    
    string tempHistoName;
    string basePath = "plane_" + to_string(sensorID);
    AIDAProcessor::tree(this)->mkdir(basePath.c_str());
    basePath = basePath + "/";
    
    tempHistoName = _hitHistoLocalName + "_" + to_string(sensorID);
    
    // Note in the local frame the origin is at the centre of the sensor. So we
    // want this to look for hits in the -x/y direction, as well at the + axis.
    // We add and subtract a constant so we know for sure we can see all hits on
    // the histogram.
    double constant = 5;
    double xMin = -(geo::gGeometry().siPlaneXSize(sensorID) / 2) - constant;
    double xMax = (geo::gGeometry().siPlaneXSize(sensorID) / 2) + constant;
    
    double yMin = -(geo::gGeometry().siPlaneYSize(sensorID) / 2) - constant;
    double yMax = (geo::gGeometry().siPlaneYSize(sensorID) / 2) + constant;
    
    int xNBin = geo::gGeometry().siPlaneXNpixels(sensorID);
    int yNBin = geo::gGeometry().siPlaneYNpixels(sensorID);
    // Special Mapsa thing
    if(sensorID==101){
        xMin = -(geo::gGeometry().siPlaneXSize(sensorID) / 2) ;
        xMax = (geo::gGeometry().siPlaneXSize(sensorID) / 2)  ;
        
        yMin = -(geo::gGeometry().siPlaneYSize(sensorID) / 2) ;
        yMax = (geo::gGeometry().siPlaneYSize(sensorID) / 2)  ;
        //         Rounding Sqrt(12)
        xNBin = 4*geo::gGeometry().siPlaneXNpixels(sensorID);
        yNBin = 4*geo::gGeometry().siPlaneYNpixels(sensorID);
    }
    
    
    AIDA::IHistogram2D *hitHistoLocal =
    AIDAProcessor::histogramFactory(this)->createHistogram2D(
        (basePath + tempHistoName).c_str(), xNBin, xMin, xMax, yNBin, yMin,
                                                             yMax);
    if (hitHistoLocal) {
        hitHistoLocal->setTitle("Hit map in the detector local frame of reference");
        _aidaHistoMap.insert(make_pair(tempHistoName, hitHistoLocal));
    } else {
        streamlog_out(ERROR1) << "Problem booking the "
        << (basePath + tempHistoName) << ".\n"
        << "Very likely a problem with path name. Switching "
        "off histogramming and continue w/o"
        << endl;
        _histogramSwitch = false;
    }
    
    // 2 should be enough because it
    // means that the sensor is wrong
    // by all its size.
    double safetyFactor = 1.2;
    double xPosition = geo::gGeometry().siPlaneXPosition(sensorID);
    double yPosition = geo::gGeometry().siPlaneYPosition(sensorID);
    double xSize = geo::gGeometry().siPlaneXSize(sensorID);
    double ySize = geo::gGeometry().siPlaneYSize(sensorID);
    int xBin = geo::gGeometry().siPlaneXNpixels(sensorID);
    int yBin = geo::gGeometry().siPlaneYNpixels(sensorID);
    
    xMin = safetyFactor * (xPosition - (0.5 * xSize));
    xMax = safetyFactor * (xPosition + (0.5 * xSize));
    
    yMin = safetyFactor * (yPosition - (0.5 * ySize));
    yMax = safetyFactor * (yPosition + (0.5 * ySize));
    
    xNBin = static_cast<int>(safetyFactor * xBin);
    yNBin = static_cast<int>(safetyFactor * yBin);
    
    if(sensorID==101){
        double x_star = abs(geo::gGeometry().siPlaneRotation1(sensorID)*(0.5 * xSize)+geo::gGeometry().siPlaneRotation2(sensorID)*(0.5 * ySize)); // was -1 ;
        double y_star = abs(geo::gGeometry().siPlaneRotation3(sensorID)*(0.5 * xSize)+geo::gGeometry().siPlaneRotation4(sensorID)*(0.5 * ySize)); // was  0 ;
        xMin = (xPosition - x_star);
        xMax = (xPosition + x_star);        
        yMin = (yPosition - y_star);
        yMax = (yPosition + y_star);        
        xNBin = static_cast<int>(4*round(abs(geo::gGeometry().siPlaneXNpixels(sensorID)*geo::gGeometry().siPlaneRotation1(sensorID)+geo::gGeometry().siPlaneYNpixels(sensorID)*geo::gGeometry().siPlaneRotation2(sensorID))));
        yNBin = static_cast<int>(4*round(abs(geo::gGeometry().siPlaneXNpixels(sensorID)*geo::gGeometry().siPlaneRotation3(sensorID)+geo::gGeometry().siPlaneYNpixels(sensorID)*geo::gGeometry().siPlaneRotation4(sensorID))));
        streamlog_out(MESSAGE0)<<"2D Hist SensorID"<<std::setw(5)<<sensorID
        <<std::setw(10)<<std::right<<std::scientific<<std::setprecision(1)<<xMin
        <<std::setw(10)<<std::right<<std::scientific<<std::setprecision(1)<<xMax
        <<std::setw(10)<<std::right<<std::scientific<<std::setprecision(1)<<yMin
        <<std::setw(10)<<std::right<<std::scientific<<std::setprecision(1)<<yMax
        <<std::setw(10)<<std::right<<std::fixed<<xNBin
        <<std::setw(10)<<std::right<<std::fixed<<yNBin<<"\n";        
    }
    
    tempHistoName = _hitHistoTelescopeName + "_" + to_string(sensorID);
    AIDA::IHistogram2D *hitHistoTelescope =
    AIDAProcessor::histogramFactory(this)->createHistogram2D(
        (basePath + tempHistoName).c_str(), xNBin, xMin, xMax, yNBin, yMin,
                                                             yMax);
    
    if (hitHistoTelescope) {
        hitHistoTelescope->setTitle("Hit map in the telescope frame of reference");
        _aidaHistoMap.insert(make_pair(tempHistoName, hitHistoTelescope));
    } else {
        streamlog_out(ERROR1) << "Problem booking the "
        << (basePath + tempHistoName) << ".\n"
        << "Very likely a problem with path name. Switching "
        "off histogramming and continue w/o"
        << endl;
        _histogramSwitch = false;
    }    
    _alreadyBookedSensorID.insert(sensorID);
    
    #endif // AIDA
}
bool EUTelMAPSAHitMaker::CheckValue(const ROOT::Internal::TTreeReaderValueBase& value) {
    if (value.GetSetupStatus() < 0) {
        streamlog_out(ERROR) << "Error " << static_cast<int>(value.GetSetupStatus())
        << "setting up reader for " << value.GetBranchName() << '\n';
        return false;
    }
    return true;
}
void  EUTelMAPSAHitMaker::getShift(unsigned int evtno, bool& goodflag, int& mapsa_evt_no, int& cms_ref_evt_no){
    int tmp_bin=_goodflaghist->FindBin(evtno)-1;
    if(tmp_bin>=0 && _goodflaghist->GetBinContent(tmp_bin+1)>1.5){
        cms_ref_evt_no  = _shift_pixel.at(tmp_bin)+evtno;
        mapsa_evt_no    = _shift_pixel.at(tmp_bin)+evtno-_shift_mpa.at(tmp_bin);
        goodflag=true;
    } else {
        cms_ref_evt_no = -10;
        mapsa_evt_no   = -10;
        goodflag=false;
    }
    //     Cut the first 50 events (Due to the long time it took to start the telescope)
    if( mapsa_evt_no<50 || mapsa_evt_no>_reader->GetEntries(true) || cms_ref_evt_no >= _total_evts_cms_ref_file || cms_ref_evt_no<0 ){
        goodflag=false;
        cms_ref_evt_no = -20;
        mapsa_evt_no   = -20;
    }
    return;
}
void EUTelMAPSAHitMaker::checkFile(const boost::filesystem::path& p){
    try{
        // does p actually exist?
        if (boost::filesystem::exists(p)){
            // is p a regular file?
            if (boost::filesystem::is_regular_file(p)){
                streamlog_out(MESSAGE0) << p << " size is " << boost::filesystem::file_size(p) << '\n';
            } else if (boost::filesystem::is_directory(p)){      // is p a directory?
                streamlog_out(ERROR) << p << " is a directory containing:\n";
                copy(boost::filesystem::directory_iterator(p), boost::filesystem::directory_iterator(), // directory_iterator::value_type
                     ostream_iterator<boost::filesystem::directory_entry>(std::cout, "\n")); // is directory_entry, which is
                // converted to a path by the
                // path stream inserter
            } else {
                streamlog_out(MESSAGE0) << p << " exists, but is neither a regular file nor a directory\n";
            }
        } else {
            streamlog_out(MESSAGE0) << p << " does not exist\n";
        }
    } catch (const boost::filesystem::filesystem_error& ex){
        streamlog_out(ERROR) << ex.what() << '\n';
    }
}
