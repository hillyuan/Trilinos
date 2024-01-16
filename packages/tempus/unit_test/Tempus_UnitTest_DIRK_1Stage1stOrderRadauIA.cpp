// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Tempus_UnitTest_RK_Utils.hpp"

namespace Tempus_Unit_Test {

using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::sublist;

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(DIRK_1Stage1stOrderRadauIA, Default_Construction)
{
  auto stepper = rcp(new Tempus::StepperDIRK_1Stage1stOrderRadauIA<double>());
  testDIRKAccessorsFullConstruction(stepper);

  // Test stepper properties.
  TEUCHOS_ASSERT(stepper->getOrder() == 1);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(DIRK_1Stage1stOrderRadauIA, StepperFactory_Construction)
{
  auto model = rcp(new Tempus_Test::SinCosModel<double>());
  testFactoryConstruction("RK Implicit 1 Stage 1st order Radau IA", model);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(DIRK_1Stage1stOrderRadauIA, AppAction)
{
  auto stepper = rcp(new Tempus::StepperDIRK_1Stage1stOrderRadauIA<double>());
  auto model   = rcp(new Tempus_Test::SinCosModel<double>());
  testRKAppAction(stepper, model, out, success);
}

}  // namespace Tempus_Unit_Test
