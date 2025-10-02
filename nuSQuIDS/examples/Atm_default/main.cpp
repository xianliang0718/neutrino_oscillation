/******************************************************************************
 *    This program is free software: you can redistribute it and/or modify     *
 *   it under the terms of the GNU General Public License as published by      *
 *   the Free Software Foundation, either version 3 of the License, or         *
 *   (at your option) any later version.                                       *
 *                                                                             *
 *   This program is distributed in the hope that it will be useful,           *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *   GNU General Public License for more details.                              *
 *                                                                             *
 *   You should have received a copy of the GNU General Public License         *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.     *
 *                                                                             *
 *   Authors:                                                                  *
 *      Carlos Arguelles (University of Wisconsin Madison)                     *
 *         carguelles@icecube.wisc.edu                                         *
 *      Jordi Salvado (University of Wisconsin Madison)                        *
 *         jsalvado@icecube.wisc.edu                                           *
 *      Christopher Weaver (University of Wisconsin Madison)                   *
 *         chris.weaver@icecube.wisc.edu                                       *
 ******************************************************************************/

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <nuSQuIDS/nuSQuIDS.h>

/*
 * The problem of solving the propagation of the atmospheric neutrinos is and
 * energy and zenith dependent problem, for this we include the class nuSQUIDSAtm
 * that allows to solve a set of nuSUIDS energy dependent objects to take in to account the
 * zenith dependence.
 */

using namespace nusquids;

// If this is defined we are doing the sterile neutrino case

// Function that gives the initial flux, for this example we set it to 1.0
double flux_function(double enu, double cz)
{
  return 1.;
}

int main()
{
  // Units and constants class
  squids::Const units;
  // Number of neutrinos (3) standard 4 1-sterile
  // Number of neutrinos, 4 is the case with one sterile neutrino
  unsigned int numneu = 3;
  bool interactions = true;

  // Minimum and maximum values for the energy and cosine zenith, notice that the energy needs to have the
  // units, if they are omitted the input is in eV.
  double Emin = 1.e-3 * units.GeV;
  double Emax = 1.e3 * units.GeV;
  double czmin = -1.0;
  double czmax = 1.0;
  int E_num = 400;
  int cz_num = 200;
  // Declaration of the atmospheric object
  std::cout << "Begin: constructing nuSQuIDS-Atm object" << std::endl;
  nuSQUIDSAtm<> nus_atm(linspace(czmin, czmax, cz_num), logspace(Emin, Emax, E_num), numneu, both, interactions);
  std::cout << "End: constructing nuSQuIDS-Atm object" << std::endl;

  std::cout << "Begin: setting mixing angles." << std::endl;
  // set mixing angles, mass differences and cp phases
  nus_atm.Set_MixingAngle(0, 1, 0.563942);
  nus_atm.Set_MixingAngle(0, 2, 0.154085);
  nus_atm.Set_MixingAngle(1, 2, 0.785398);

  nus_atm.Set_SquareMassDifference(1, 7.65e-05);
  nus_atm.Set_SquareMassDifference(2, 0.00247);

  nus_atm.Set_CPPhase(0, 2, 0);
  if (numneu > 3)
  {
    nus_atm.Set_SquareMassDifference(3, -1.);
    nus_atm.Set_MixingAngle(1, 3, 0.160875);
  }
  std::cout << "End: setting mixing angles." << std::endl;

  // Setup integration precision
  nus_atm.Set_rel_error(1.0e-6);
  nus_atm.Set_abs_error(1.0e-6);
  nus_atm.Set_GSL_step(gsl_odeiv2_step_rk4);

  // Array that contains the values of the energies and cosine of the zenith, is the same length for every zenith
  auto e_range = nus_atm.GetERange();
  auto cz_range = nus_atm.GetCosthRange();

  std::cout << "Begin: setting initial state." << std::endl;

  // Construct the initial state, we set a flat spectra in zenith and log-energy
  marray<double, 4> inistate{nus_atm.GetNumCos(), nus_atm.GetNumE(), 2, numneu};
  std::fill(inistate.begin(), inistate.end(), 0);

  std::ifstream inputTXT("/home/xianliang/download/musairs/output1.txt");
  if (!inputTXT)
  {
    std::cout << "无法打开文件" << std::endl;
    return 1;
  }
  else
  {
    std::cout << "文件打开成功" << std::endl;
  }

  marray<double, 1> costh_marray = linspace(czmin, czmax, cz_num);
  marray<double, 1> enu_marray = logspace(log(Emin) - log(units.eV), log(Emax) - log(units.eV), E_num);

  std::string line;
  int cnt;
  int pdgIDValue;
  double Ek0Value;
  double ZenithValue;
  double WeightValue;
  double HeightValue;

  while (inputTXT >> pdgIDValue >> Ek0Value >> ZenithValue >> WeightValue >> HeightValue)
  {

    int rho = (pdgIDValue > 0) ? 0 : 1;     // Neutrino or antineutrino
    int flv = std::abs(pdgIDValue) / 2 - 6; // Flavor index
    // std::cout << "pdgID: " << pdgIDValue << " Ek0: " << Ek0Value << " Zenith: " << ZenithValue << " Weight: " << WeightValue << std::endl;

    auto cthit = std::lower_bound(costh_marray.begin(), costh_marray.end(), ZenithValue - 1);
    if (cthit >= costh_marray.end())
      std::cout << "SQUIDS::GetExpectation C ValueD : x value not in the array.";
    if (cthit != costh_marray.begin())
      cthit--;
    size_t cth_M = std::distance(costh_marray.begin(), cthit);

    double logE = log(Ek0Value * 1e9);
    auto logeit = std::lower_bound(enu_marray.begin(), enu_marray.end(), logE);
    if (logeit >= enu_marray.end())
      std::cout << "SQUIDS::GetExpectation E ValueD : x value not in the array.";
    if (logeit != enu_marray.begin())
      logeit--;
    size_t enu_M = std::distance(enu_marray.begin(), logeit);

    inistate[cth_M][enu_M][rho][flv] += WeightValue;
  }

  inputTXT.close();


  nus_atm.Set_initial_state(inistate, flavor);
  std::cout << "End: setting initial state." << std::endl;

  // Set to true the monitoring prgress bar and the vacuum oscillations
  nus_atm.Set_ProgressBar(true);
  nus_atm.Set_IncludeOscillations(true);

  // Here we do the evolution of all the states
  std::cout << "Begin: Evolution" << std::endl;
  nus_atm.EvolveState();
  std::cout << "End: Evolution" << std::endl;

  std::ofstream outFile("/home/xianliang/download/musairs/output2.txt");
  if (!outFile)
  {
    std::cerr << "无法打开输出文件" << std::endl;
    return 1;
  }

  double lEmin = log10(Emin);
  double lEmax = log10(Emax);

  for (double cz = czmin; cz < czmax; cz += (czmax - czmin) / (double)cz_num)
  {
    for (double lE = lEmin; lE < lEmax; lE += (lEmax - lEmin) / (double)E_num)
    {
      double E = pow(10.0, lE);
      for (int fl = 0; fl < numneu; fl++)
      {
        int pdgid = (fl + 6) * 2;
        outFile << pdgid << " " << E / units.GeV << " " << cz << " " << nus_atm.EvalFlavor(fl, cz, E, 0) << std::endl;
        outFile << -pdgid << " " << E / units.GeV << " " << cz << " " << nus_atm.EvalFlavor(fl, cz, E, 1) << std::endl;
      }
    }
  }
  return 0;
}