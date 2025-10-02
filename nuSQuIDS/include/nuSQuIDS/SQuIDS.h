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

/**
 \mainpage
 
 Introduction
 ------------
 
 Simple Quantum Integro-Differential Solver (SQuIDS) is a C++ code designed to 
 solve semi-analytically the evolution of a set of density matrices and scalar 
 functions. SQuIDS provides a base class from which users can derive new classes 
 to include new non-trivial terms from the right hand sides of density matrix 
 equations. The code was designed in the context of solving neutrino oscillation 
 problems, but can be applied in any problem that involves solving the quantum 
 evolution of a collection of particles with Hilbert space of dimension up to 6.
 
 Dependencies
 ------------
 
 SQuIDS is written using features of C++11 for efficiency, so it requires a
 reasonably new C++ compiler. GCC 4.8.1 or newer or Clang 3.3 or newer is
 recommended, but the code is written to be portable and should work with 
 other compilers as well. 
 
 The GNU Scientific Library (GSL),
 version 1.15 or newer is needed by SQuIDS to perform differential equation
 integration. The latest version can be found at 
 https://www.gnu.org/software/gsl/ .
 
 Finally, gnuplot (http://www.gnuplot.info) version 4.0 or newer is required 
 to generate plots from the example programs, but is not needed by the library 
 itself.
 
 Building
 --------
 
 First, run the included configure script: 
 
     ./configure
 
 It may be useful to run it with the `--help` flag in order to find out all
 available options. 
 
 After configuring the build, simply compile with the generated makefile:
 
     make
 
 You can verify that the library is working correctly by running the unit 
 tests with `make test`, and you can install the built library with 
 `make install`.
 
 Usage
 -----
 
 To learn about using the library, see the examples in the `examples` 
 subdirectory. 
 
 In general, implementing a system using this library requires defining a
 new subclass of the SQUIDS base class which encapsulates the problem, and
 using the SU_vector class to represent quantum mechanical states and 
 operators.
 
 */

#ifndef SQUIDS_H
#define SQUIDS_H

#if __cplusplus < 201103L
#error C++11 compiler required. Update your compiler and use the flag -std=c++11
#endif

#include "SUNalg.h"

#include <iosfwd>
#include <vector>
#include <memory>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

namespace squids{

///\brief SQuIDS main class
///
///density matrix kinetic equation solver
class SQuIDS {
 protected:
  
  ///\brief Structure that contains the node state
  struct SU_state
  {
    ///\brief Vector of SU(N) vectors that represents the quantum part of the state
    std::unique_ptr<SU_vector[]> rho;
    ///\brief Vector of scalars that represents the classic part of the state
    double* scalar; //not owned
  };

 private:
  bool CoherentRhoTerms,NonCoherentRhoTerms,OtherRhoTerms,GammaScalarTerms,OtherScalarTerms,AnyNumerics;
  bool is_init;
  bool adaptive_step;
 
  std::vector<double> x;
  double t;
  double t_ini;
  unsigned int nsteps;
  
  unsigned int size_rho;
  unsigned int size_state;
  
  std::unique_ptr<double[]> system;
  gsl_odeiv2_step_type const* step; //not owned
  gsl_odeiv2_system sys;
  
  double h;
  double h_min;
  double h_max;
  double abs_error;
  double rel_error;
  std::unique_ptr<SU_state[]> dstate;
  double* last_dstate_ptr;
  double* last_estate_ptr;
  
  //***************************************************************
  ///\brief Sets the evolution state and derivative system pointer for GSL use
  ///\param sp the backing storage for the state during evolution (estate)
  ///\param dp the backing storage for the derivative during evolution (dstate)
  void set_system_pointers(double* sp, double* dp);
  //interface function called by GSL
  friend int RHS(double ,const double*,double*,void*);
 
 protected:
  ///The number of nodes in the system
  unsigned int nx;
  ///The dimension of the hilbert space used for each density matrix
  unsigned int nsun;
  ///The number of density matrices per node
  unsigned int nrhos;
  ///The number of scalars per node
  unsigned int nscalars;
  ///contains constants and basis transformation
  Const params;
  ///the state of the system
  std::unique_ptr<SU_state[]> state;
  ///the state of the system during an evolution step
  std::unique_ptr<SU_state[]> estate;
  ///\brief Sets the current time of the system
  ///\param t_ Time to set.
  ///\warning Do not use this function unless you are setting the same time
  /// as the time to the system has already being evolved at.
  void Set_t(double t_) { t = t_; }
 public:
  //****************
  //Constructors
  //****************
  SQuIDS();
  //***************************************************************
  ///\brief Constructs a SQUIDS object
  ///
  ///\param nx Number of components in the array "x"
  ///\param dim Dimension of SU(n)
  ///\param nrho Number of density matrix in every "x" site
  ///\param nscalar Number of scalars in every "x" site
  ///\param ti initial value for the evolution parameter t
  SQuIDS(unsigned int nx, unsigned int dim, unsigned int nrho, unsigned int nscalar, double ti=0.0);

  ///\brief Move constructs a SQUIDS object from an existing object
  SQuIDS(SQuIDS&&);

  //***************************************************************
  virtual ~SQuIDS();
  
  //***************************************************************
  ///\brief Move assigns a SQUIDS object from an existing object
  SQuIDS& operator=(SQuIDS&&);

  //***************************************************************
  ///\brief Initializes a SQUIDS object
  ///
  ///\param nx Number of components in the array "x"
  ///\param dim Dimension of SU(n)
  ///\param nrho Number of density matrix in every "x" site
  ///\param nscalar Number of scalars in every "x" site
  ///\param ti initial value for the evolution parameter t
  void ini(unsigned int nx, unsigned int dim, unsigned int nrho, unsigned int nscalar, double ti=0.0);

  //***************************************************************
  ///\brief Get the number of nodes in the system
  unsigned int Get_nx() const{ return(nx); }
  
  ///\brief Get the number of density matrices per node
  unsigned int Get_nrhos() const{ return(nrhos); }
  
  ///\brief Get the number of scalars per node
  unsigned int Get_nscalars() const{ return(nscalars); }
  
  //***************************************************************
  ///\brief Set the range of values for the array "x"
  ///\param xini  x_min
  ///\param xend  x_max
  ///\param scale "log or lin" type of scale
  void Set_xrange(double xini, double xend, std::string scale);
  
  ///\brief Set the range of values for the array "x"
  ///\param xs  The x values to set
  ///\pre xs.size()==nx
  void Set_xrange(const std::vector<double>& xs);
  
  ///\brief Get the range of values for the array "x"
  const std::vector<double>& Get_xrange() const{ return(x); }

  //***************************************************************
  ///\brief Returns the closes position in the array x for the value given
  ///\param x value of x to look for
  unsigned int Get_i(double x) const;

  //***************************************************************
  ///\brief Returns de value in the position "i"
  ///\param i node position
  double Get_x(unsigned int i) const{return x[i];}

  //***************************************************************
  //virtual functions defined in the dervied class
  ///\brief H0 time independent evolution operator
  virtual SU_vector H0(double x, unsigned int irho) const{ return SU_vector(nsun);}
  ///\brief H1 time dependent evolution operator
  virtual SU_vector HI(unsigned int ix, unsigned int irho, double t) const{ return SU_vector(nsun);}
  ///\brief Attenuation and/or decoherence operator
  ///\param ix Index in the x-array
  ///\param t time
  virtual SU_vector GammaRho(unsigned int ix, unsigned int irho, double t) const{ return SU_vector(nsun);}
  ///\brief Function containing other possible operations, like non linear terms in rho
  ///or terms involving the scalar functions
  ///\param ix Index in the x-array
  ///\param t time
  virtual SU_vector InteractionsRho(unsigned int ix, unsigned int irho, double t) const{ return SU_vector(nsun);}
  ///\brief Attenuation for the scalar functions
  ///\param ix Index in the x-array
  ///\param t time
  virtual double GammaScalar(unsigned int ix, unsigned int irho, double t) const{return 0.0;}
  ///\brief Other possible interaction terms for the scalar functions.
  ///\param ix Index in the x-array
  ///\param t time
  virtual double InteractionsScalar(unsigned int ix, unsigned int irho, double t) const{return 0.0;}
  ///\brief Function to be evaluated before the derivative
  ///\param t time
  ///
  /// This function enables the user to perform operations or updates before the derivative.
  virtual void PreDerive(double t){}

  //***************************************************************
  ///\brief Computes the right hand side of the kinetic equation.
  ///\param t time
  void Derive(double t);

  //***************************************************************
  ///\brief Numerical evolution of the state using GSL
  ///\param dt evolution time interval.
  void Evolve(double dt);

  //***************************************************************
  //functions to set parameters in the object.
  // string -> parameter
  // other -> value
  //the possible parameters are:
  // ******* bool *************
  // CoherentInteractions
  // NonCoherentInteractions 
  // OtherInt
  // ScalarInteractions
  // AntiNeutrinos
  // ******* double *********
  // t -> time
  // Units -> time unit
  // The parameters in the struc "const" are also available from this functions
  // ******* int ************
  // nx
  // nsun
  // nrhos
  // nscalars
  // The parameters in the struc "const" are also available from this functions
  ///\brief Sets the GSL stepper function (numerical method)
  ///\param opt GSL step function
  void Set_GSL_step(gsl_odeiv2_step_type const* opt);

  ///\brief Turns on and off adaptive runge-kutta stepping
  ///\param opt If true: uses adaptive stepping, else: it does not.
  void Set_AdaptiveStep(bool opt);
  ///\brief Activate coherent interaction
  void Set_CoherentRhoTerms(bool opt);
  ///\brief Activate noncoherent interaction
  void Set_NonCoherentRhoTerms(bool opt);
  ///\brief Activate other SU_vector interactions
  void Set_OtherRhoTerms(bool opt);
  ///\brief Activate other scalar interactions
  void Set_GammaScalarTerms(bool opt);
  ///\brief Activate other scalar interactions
  void Set_OtherScalarTerms(bool opt);
  ///\brief If set to false will disable all numerics
  void Set_AnyNumerics(bool opt);
  ///\brief Set the minimum runge-kutta step
  void Set_h_min(double opt);
  ///\brief Set the maximum runge-kutta step
  void Set_h_max(double opt);
  ///\brief Set the initial runge-kutta step
  void Set_h(double opt);
  ///\brief Get the minimum runge-kutta step
  double Get_h_min() const;
  ///\brief Get the maximum runge-kutta step
  double Get_h_max() const;
  ///\brief Get the initial runge-kutta step
  double Get_h() const;
  ///\brief Set the numerical relative error
  void Set_rel_error(double opt);
  ///\brief Set the numerical absolute error
  void Set_abs_error(double opt);
  ///\brief Set the number of steps when not using adaptive stepping
  void Set_NumSteps(unsigned int opt);
   ///\brief Get the numerical relative error
  double Get_rel_error() const;
  ///\brief Get the numerical absolute error
  double Get_abs_error() const;
  ///\brief Get the number of steps when not using adaptive stepping
  double Get_NumSteps() const;

  //***************************************************************
  ///\brief Returns the expectation value for a given operator for a give state irho in a node ix.
  ///\param op operator
  ///\param irho index of rho
  ///\param ix index in the array "x"
  double GetExpectationValue(SU_vector op, unsigned int irho, unsigned int ix) const;

  //***************************************************************
  ///\brief Returns the expectation value for a given operator for a give state irho in a node ix when
  /// consindering averaging of oscillations. Oscillations are averaged if the phase is larger than scale
  /// the bool pointer is true for all scales that are averaged out.
  ///\param op operator
  ///\param irho index of rho
  ///\param ix index in the array "x"
  ///\param scale scale upon which oscillations will be averaged out
  ///\param avg bool array which is true for all scales that were averaged out
  double GetExpectationValue(SU_vector op, unsigned int nrh, unsigned int i, double scale, std::vector<bool>& avr) const;

  //***************************************************************
  ///\brief Returns the intermediate state using linear interpolation in "x"
  ///\param irho index of rho
  ///\param x value of x
  SU_vector GetIntermediateState(unsigned int irho, double x) const;

  //***************************************************************
  ///\brief Returns the expectation value for a given operator for the rho given by irho
  /// and using linear interpolation in "x"
  ///\param op operator 
  ///\param irho index of rho
  ///\param x value of x
  double GetExpectationValueD(const SU_vector& op, unsigned int irho, double x) const;

  //***************************************************************
  ///\brief Returns the expectation value for a given operator for a give state irho in a node ix when
  /// consindering averaging of oscillations. Oscillations are averaged if the phase is larger than scale
  /// the bool pointer is true for all scales that are averaged out. It uses linear interpolation in "x"
  ///\param op operator
  ///\param irho index of rho
  ///\param x value of x
  ///\param scale scale upon which oscillations will be averaged out
  ///\param avg bool array which is true for all scales that were averaged out
  double GetExpectationValueD(const SU_vector& op, unsigned int nrh, double x, double scale, std::vector<bool>& avr) const;

  ///This type encapsulates the temporary storage needed by GetExpectationValueD
  struct expectationValueDBuffer{
  private:
    SU_vector state, op;
    friend class SQuIDS;
  public:
    ///Construct temporary storage suitable for problems with the given dimension
    ///\param dim the problem dimension
    expectationValueDBuffer(unsigned int dim):
    state(dim),op(dim){}
  };

  ///\brief Returns the expectation value for a given operator for the rho given by irho
  /// and using linear interpolation in "x"
  ///\param op operator
  ///\param irho index of rho
  ///\param x value of x
  ///\param buf a buffer containing the necessary temporary storage. Must have
  ///           been initialized to the same dimension as the problem.
  double GetExpectationValueD(const SU_vector& op, unsigned int irho, double x, expectationValueDBuffer& buf) const;

 //***************************************************************
  ///\brief Returns the expectation value for a given operator for a give state irho in a node ix when
  /// consindering averaging of oscillations. Oscillations are averaged if the phase is larger than scale
  /// the bool pointer is true for all scales that are averaged out. It uses linear interpolation in "x"
  ///\param op operator
  ///\param irho index of rho
  ///\param x value of x
  ///\param scale scale upon which oscillations will be averaged out
  ///\param avg bool array which is true for all scales that were averaged out
  double GetExpectationValueD(const SU_vector& op, unsigned int nrh,  double x, expectationValueDBuffer& buf, double scale, std::vector<bool>& avr) const;

  ///\brief Returns the initial time of the system
  double Get_t_initial() const{ return(t_ini); }
  ///\brief Returns the current time of the system
  double Get_t() const{ return(t); }
  

  ///\brief Returns the parameter object for this system
  const Const& GetParams() const{ return(params); }
};

} //namespace squids
  
#endif
