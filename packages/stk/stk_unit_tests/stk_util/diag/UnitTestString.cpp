// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include "gtest/gtest.h"
#include "stk_util/diag/String.hpp"      // for String, operator==, StringBase
#include "stk_util/diag/StringUtil.hpp"  // for case_strcmp
#include <string>                        // for string

TEST(StkString, case_insensitivity)
{
  sierra::String ABC("ABC");
  sierra::String abc("abc");
  bool strings_match_by_case_insensitivity = ABC == abc;
  EXPECT_TRUE(strings_match_by_case_insensitivity);
}

TEST(StkString, case_sensitive_comparison_with_sierra_string)
{
  sierra::String ABC("ABC");
  sierra::String abc("abc");
  bool strings_do_not_match = sierra::case_strcmp(ABC, abc);
  EXPECT_FALSE(strings_do_not_match);
}

TEST(StkString, case_sensitive_comparison_with_std_string)
{
  std::string ABC("ABC");
  std::string abc("abc");
  bool strings_do_not_match = sierra::case_strcmp(ABC, abc);
  EXPECT_FALSE(strings_do_not_match);
}

