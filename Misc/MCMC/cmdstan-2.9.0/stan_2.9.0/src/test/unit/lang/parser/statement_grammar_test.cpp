#include <gtest/gtest.h>
#include <exception>
#include <test/unit/lang/utility.hpp>

TEST(langParserStatement2Grammar, addConditionalCondition) {
  test_parsable("conditional_condition_good");
  test_throws("conditional_condition_bad_1",
              "conditions in if-else");
  test_throws("conditional_condition_bad_2",
              "conditions in if-else");
}

TEST(langParserStatementGrammar, validateIntExpr2) {
  test_parsable("validate_int_expr2_good");
  test_throws("validate_int_expr2_bad1",
              "expression denoting integer required");
  test_throws("validate_int_expr2_bad2",
              "expression denoting integer required");
  test_throws("validate_int_expr2_bad3",
              "expression denoting integer required");
  test_throws("validate_int_expr2_bad4",
              "expression denoting integer required");
}

TEST(langParserStatementGrammar, validateAllowSample) {
  test_throws("validate_allow_sample_bad1",
              "Sampling statements (~) and increment_log_prob() are");
  test_throws("validate_allow_sample_bad2",
              "Sampling statements (~) and increment_log_prob() are");
  test_throws("validate_allow_sample_bad3",
              "Sampling statements (~) and increment_log_prob() are");
}

TEST(langParserStatementGrammarDef, jacobianAdjustmentWarning) {
  test_parsable("validate_jacobian_warning_good");
  test_warning("validate_jacobian_warning1",
               "If so, you need to call increment_log_prob() with the log");
  test_warning("validate_jacobian_warning2",
               "If so, you need to call increment_log_prob() with the log");
  test_warning("validate_jacobian_warning3",
               "If so, you need to call increment_log_prob() with the log");
  test_warning("validate_jacobian_warning4",
               "If so, you need to call increment_log_prob() with the log");
  test_warning("validate_jacobian_warning5",
               "If so, you need to call increment_log_prob() with the log");
  test_warning("validate_jacobian_warning6",
               "If so, you need to call increment_log_prob() with the log");
}

TEST(langParserStatementGrammarDef, jacobianUserFacing) {
  test_warning("validate_jacobian_warning_user",
               "exp(y[1]) ~ normal(...)");
}

TEST(langParserStatementGrammarDef, comparisonsInBoundsTest) {
  test_parsable("validate_bounds_comparison");
  EXPECT_THROW(is_parsable("src/test/test-models/bad/lang/bad_bounds1.stan"),
               std::invalid_argument);
}

TEST(langParserStatementGrammar, validateAssignmentTypes) {
  test_throws("bad_var_assignment_type1",
              "mismatch in assignment");
  test_throws("bad_var_assignment_type2",
              "mismatch in assignment");
  test_throws("bad_var_assignment_vec_arr",
              "mismatch in assignment");
}

TEST(langParserStatementGrammar, assignRealToIntMessage) {
  test_throws("assign_real_to_int",
              "PARSER EXPECTED: <expression assignable to left-hand side>");
}

TEST(langParserStatementGrammar, useCdfWithSamplingNotation) {
  test_throws("cdf-sample",
              "CDF and CCDF functions may not be used with sampling notation.");
  test_throws("ccdf-sample",
              "CDF and CCDF functions may not be used with sampling notation.");
  test_throws("multiply_sample",
              "Only distribution names can be used with sampling (~) notation");
  test_throws("binomial_coefficient_sample",
              "Only distribution names can be used with sampling (~) notation");
}

