// Author: Wahid Redjeb - wahid.redjeb@cern.ch
// Date: 09/2021

#ifndef DataFormats_HGCalReco_SuperTrackster_h
#define DataFormats_HGCalReco_SuperTrackster_h

#include <array>
#include <vector>

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"

namespace ticl {
    class SuperTrackster {
        public:
            SuperTrackster()
            :   rawEnergy_(0.f),
                hasTrack(false){}

            SuperTrackster(const std::vector<uint32_t>& tracksters)
                :   tracksterIdxs_({tracksters}) {}
            
            SuperTrackster(const std::vector<uint32_t>& tracksters, const uint32_t& track)
                :   tracksterIdxs_({tracksters}),
                    trackIdx_(track),
                    hasTrack(true) {}

            inline const uint32_t trackIdx() const { return trackIdx_;}
            void setTrackPtr(const uint32_t trackIdx) { 
                trackIdx_ = trackIdx;
                hasTrack = true;
                }

            inline const std::vector<uint32_t> trackstersIdxs() const {return tracksterIdxs_;}
            void setTracksters(const std::vector<uint32_t>& tracksterIdxs) { tracksterIdxs_ = tracksterIdxs; }
            void addTrackster(const uint32_t tracksterIdx) { 
                tracksterIdxs_.push_back(tracksterIdx);
            }

            void setRawEnergy(float rawEnergy) { rawEnergy_ = rawEnergy;}
        private:
            float rawEnergy_;            
            std::vector<uint32_t> tracksterIdxs_; 
            uint32_t trackIdx_;
            bool hasTrack;

    };
}

#endif