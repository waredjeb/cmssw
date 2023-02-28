#include "FWCore/ParameterSet/interface/ValidatedPluginFactoryMacros.h"
#include "FWCore/ParameterSet/interface/ValidatedPluginMacros.h"
#include "LinkingAlgoFactory.h"
#include "LinkingAlgoByDirectionGeometric.h"
#include "LinkingAlgoByGNN.h"
#include "SmoothingAlgoByMLP.h"

EDM_REGISTER_VALIDATED_PLUGINFACTORY(LinkingAlgoFactory, "LinkingAlgoFactory");

DEFINE_EDM_VALIDATED_PLUGIN(LinkingAlgoFactory,
                            ticl::LinkingAlgoByDirectionGeometric,
                            "LinkingAlgoByDirectionGeometric");

DEFINE_EDM_VALIDATED_PLUGIN(LinkingAlgoFactory, ticl::LinkingAlgoByGNN, "LinkingAlgoByGNN");
DEFINE_EDM_VALIDATED_PLUGIN(LinkingAlgoFactory, ticl::SmoothingAlgoByMLP, "SmoothingAlgoByMLP");
