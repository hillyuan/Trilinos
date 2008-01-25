// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER


/** \file
\brief  Unit test (CubatureDirect,CubatureTensor): correctness of
        integration of monomials for 3D reference cells.
\author Created by P. Bochev and D. Ridzal.
*/

#include "Intrepid_CubatureDirect.hpp"
#include "Intrepid_CubatureTensor.hpp"
#include "Intrepid_Utils.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

using namespace Intrepid;


/*
  Monomial evaluation.
    in 1D, for point p(x)    : x^xDeg
    in 2D, for point p(x,y)  : x^xDeg * y^yDeg
    in 3D, for point p(x,y,z): x^xDeg * y^yDeg * z^zDeg
*/
double computeMonomial(Point<double> p, int xDeg, int yDeg=0, int zDeg=0) {
  double val = 1.0;
  int polydeg[3];
  polydeg[0] = xDeg; polydeg[1] = yDeg; polydeg[2] = zDeg;
  for (int i=0; i<p.getDim(); i++) {
    val *= std::pow(p.getCoordinates()[i],polydeg[i]);
  }
  return val;
}


/*
  Computes integrals of monomials over a given reference cell.
*/
double computeIntegral(ECell cellType, int cubDegree, int xDeg, int yDeg, int zDeg) {

  Teuchos::RCP< Cubature<double> > myCub;  
  CubatureDirect<double> dCub;
  CubatureTensor<double> tCub;
  double val = 0.0;

  int ambientDim =  MultiCell<double>::getTopologicalDim(cellType);

  switch (cellType) {

    case CELL_TET:
        myCub = Teuchos::rcp(&dCub, false);
      break;

    case CELL_HEX:
        myCub = Teuchos::rcp(&tCub, false);
      break;

    case CELL_TRIPRISM:
        myCub = Teuchos::rcp(&tCub, false);
      break;

    default:
      TEST_FOR_EXCEPTION((cellType != CELL_TET) && (cellType != CELL_HEX) && (cellType != CELL_TRIPRISM),
                          std::invalid_argument,
                          ">>> ERROR (Unit Test -- Cubature -- 3D Monomial): Invalid cell type.");
  } // end switch

  int numCubPoints = myCub->getNumPoints(cellType, cubDegree);

  Teuchos::Array< Point<double> > cubPoints;
  Teuchos::Array<double> cubWeights;

  Point<double> tempPoint(ambientDim);
  cubPoints.assign(numCubPoints,tempPoint);
  cubWeights.assign(numCubPoints,0.0);

  myCub->getCubature(cubPoints, cubWeights, cellType, cubDegree);

  for (int i=0; i<numCubPoints; i++) {
    val += computeMonomial(cubPoints[i], xDeg, yDeg, zDeg)*cubWeights[i];
  }

  return val;
}


int main(int argc, char *argv[]) {

  // This little trick lets us print to std::cout only if
  // a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);
 
  *outStream \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|                 Unit Test (CubatureDirect,CubatureTensor)                   |\n" \
  << "|                                                                             |\n" \
  << "|     1) Computing integrals of monomials on reference cells in 3D            |\n" \
  << "|                - no BLAS, i.e. standard addition loops -                    |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n"\
  << "| TEST 1: integrals of monomials in 3D (non-BLAS version)                     |\n"\
  << "===============================================================================\n";

  // >>> ASSUMPTION: max polynomial degree integrated exactly is the same for
  // >>>             edges and triangles !!!
  // internal variables:
  int                                      errorFlag = 0;
  int                                      polyCt = 0;
  int                                      offset = 0;
  Teuchos::Array< Teuchos::Array<double> > testInt;
  Teuchos::Array< Teuchos::Array<double> > analyticInt;
  Teuchos::Array<double>                   tmparray(1);
  double                                   reltol = 1.0e+04 * INTREPID_TOL;
  int                                      numPoly = (INTREPID_MAX_CUBATURE_DEGREE_EDGE+1)*
                                                     (INTREPID_MAX_CUBATURE_DEGREE_EDGE+2)*
                                                     (INTREPID_MAX_CUBATURE_DEGREE_EDGE+3)/6;
  testInt.assign(numPoly, tmparray);
  analyticInt.assign(numPoly, tmparray);

  // get names of files with analytic values
  std::string basedir = "./data";
  std::stringstream namestream[3];
  std::string filename[3];
  namestream[0] << basedir << "/TET_integrals" << ".dat";
  namestream[0] >> filename[0];
  namestream[1] << basedir << "/HEX_integrals" << ".dat";
  namestream[1] >> filename[1];
  namestream[2] << basedir << "/TRIPRISM_integrals" << ".dat";
  namestream[2] >> filename[2];

  // reference cells tested
  ECell             testType[] = {CELL_TET, CELL_HEX, CELL_TRIPRISM};
  // format of data files with analytic values
  TypeOfExactData dataFormat[] = {INTREPID_UTILS_SCALAR, INTREPID_UTILS_FRACTION, INTREPID_UTILS_FRACTION};

  // compute and compare integrals
  try {
    for (int cellCt=0; cellCt < 3; cellCt++) {
      *outStream << "\nIntegrals of monomials on a reference " << MultiCell<double>::getCellName(testType[cellCt]) << ":\n";
      std::ifstream filecompare(&filename[cellCt][0]);
      // compute integrals
      for (int cubDeg=0; cubDeg <= INTREPID_MAX_CUBATURE_DEGREE_EDGE; cubDeg++) {
        polyCt = 0;
        testInt[cubDeg].resize((cubDeg+1)*(cubDeg+2)*(cubDeg+3)/6);
        for (int xDeg=0; xDeg <= cubDeg; xDeg++) {
          for (int yDeg=0; yDeg <= cubDeg-xDeg; yDeg++) {
            for (int zDeg=0; zDeg <= cubDeg-xDeg-yDeg; zDeg++) {
              testInt[cubDeg][polyCt] = computeIntegral(testType[cellCt], cubDeg, xDeg, yDeg, zDeg);
              polyCt++; 
            }
          }
        }
      }
      // get analytic values
      if (filecompare.is_open()) {
        getAnalytic(analyticInt, filecompare, dataFormat[cellCt]);
        // close file
        filecompare.close();
      }
      // perform comparison
      for (int cubDeg=0; cubDeg <= INTREPID_MAX_CUBATURE_DEGREE_EDGE; cubDeg++) {
        polyCt = 0;
        offset = 0;
        for (int xDeg=0; xDeg <= cubDeg; xDeg++) {
          for (int yDeg=0; yDeg <= cubDeg-xDeg; yDeg++) {
            for (int zDeg=0; zDeg <= cubDeg-xDeg-yDeg; zDeg++) {
              double abstol = ( analyticInt[polyCt+offset][0] == 0.0 ? reltol : std::fabs(reltol*analyticInt[polyCt+offset][0]) );
              double absdiff = std::fabs(analyticInt[polyCt+offset][0] - testInt[cubDeg][polyCt]);
              *outStream << "Cubature order " << std::setw(2) << std::left << cubDeg << " integrating "
                         << "x^" << std::setw(2) << std::left << xDeg << " * y^" << std::setw(2) << yDeg
                         << " * z^" << std::setw(2) << zDeg << ":" << "   "
                         << std::scientific << std::setprecision(16) << testInt[cubDeg][polyCt] << "   " << analyticInt[polyCt+offset][0] << "   "
                         << std::setprecision(4) << absdiff << "   " << "<?" << "   " << abstol << "\n";
              if (absdiff > abstol) {
                errorFlag++;
                *outStream << std::right << std::setw(118) << "^^^^---FAILURE!\n";
              }
              polyCt++;
            }
            offset = offset + INTREPID_MAX_CUBATURE_DEGREE_EDGE - cubDeg;
          }
          offset = offset + (INTREPID_MAX_CUBATURE_DEGREE_EDGE - cubDeg)*(INTREPID_MAX_CUBATURE_DEGREE_EDGE - cubDeg + 1)/2;
        }
        *outStream << "\n";
      }
      *outStream << "\n";
    }  // end for cellCt
  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1;
  };


  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);

  return errorFlag;
}
