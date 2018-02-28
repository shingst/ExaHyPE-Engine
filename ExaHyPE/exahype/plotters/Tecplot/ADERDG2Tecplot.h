/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon 
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 * 
 * @authors: Sven Koeppel <koeppel@fias.uni-frankfurt.de>,
 *           Maurizio Tavelli <m.tavelli@unitn.it>
 **/
 
#ifndef _EXAHYPE_PLOTTERS_ADERDG_2_TECPLOT_H_
#define _EXAHYPE_PLOTTERS_ADERDG_2_TECPLOT_H_

#include "exahype/plotters/Plotter.h"

namespace exahype {
  namespace plotters {
    class ADERDG2Tecplot;
  }
}

/**
 * <h2>The Tecplot writer from the Trento code</h2>
 *
 * This plotter is a stub to connect the Fortran plotting code to ExaHyPE.
 *
 * Note, to enable this plotter, you need
 *
 *   export PROJECT_CFLAGS="-DTECPLOT"
 * 
 */
class exahype::plotters::ADERDG2Tecplot : public exahype::plotters::Plotter::Device {
 public:
  static tarch::logging::Log _log;
  
  int _orderPlusOne, _solverUnknowns, _writtenUnknowns;
  exahype::solvers::Solver::Type _solverType;

  static std::string getIdentifier();

  ADERDG2Tecplot(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing, exahype::solvers::Solver::Type type);

  virtual ~ADERDG2Tecplot();

  virtual void init(const std::string& filename, int orderPlusOne, int solverUnknowns, int writtenUnknowns, const std::string& select);

  virtual void plotPatch(
        const int cellDescriptionsIndex,
        const int element);

  virtual void startPlotting( double time );
  virtual void finishPlotting();
};

#endif/* _EXAHYPE_PLOTTERS_ADERDG_2_TECPLOT_H_ */
