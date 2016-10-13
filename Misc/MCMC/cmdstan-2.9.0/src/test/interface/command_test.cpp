#include <stan/services/error_codes.hpp>
#include <cmdstan/command.hpp>
#include <gtest/gtest.h>
#include <string>
#include <test/utility.hpp>
#include <stdexcept>
#include <boost/math/policies/error_handling.hpp>

using cmdstan::test::convert_model_path;
using cmdstan::test::count_matches;
using cmdstan::test::multiple_command_separator;
using cmdstan::test::run_command;
using cmdstan::test::run_command_output;

TEST(StanUiCommand, countMatches) {
  EXPECT_EQ(-1, count_matches("", ""));
  EXPECT_EQ(-1, count_matches("", "abc"));

  EXPECT_EQ(0, count_matches("abc", ""));
  EXPECT_EQ(0, count_matches("abc", "ab"));
  EXPECT_EQ(0, count_matches("abc", "dab"));
  EXPECT_EQ(0, count_matches("abc", "abde"));

  EXPECT_EQ(0, count_matches("aa","a"));
  EXPECT_EQ(1, count_matches("aa","aa"));
  EXPECT_EQ(1, count_matches("aa","aaa"));
  EXPECT_EQ(2, count_matches("aa","aaaa"));
}

void test_sample_prints(const std::string& base_cmd) {
  std::string cmd(base_cmd);
  cmd += " num_samples=100 num_warmup=100";
  std::string cmd_output = run_command(cmd).output;
  // transformed data
  EXPECT_EQ(1, count_matches("x=", cmd_output)); 
  // transformed parameters
  EXPECT_TRUE(count_matches("z=", cmd_output) >= 200); 
  // model
  EXPECT_TRUE(count_matches("y=", cmd_output) >= 200);
  // generated quantities [only on saved iterations, should be num samples]
  EXPECT_TRUE(count_matches("w=", cmd_output) == 100);
}

void test_optimize_prints(const std::string& base_cmd) {
  std::string cmd(base_cmd);
  // cmd += " num_samples=100 num_warmup=100";
  std::string cmd_output = run_command(cmd).output;
  // transformed data
  EXPECT_EQ(1, count_matches("x=", cmd_output)); 
  // transformed parameters
  EXPECT_TRUE(count_matches("z=", cmd_output) >= 1); 
  // model
  EXPECT_TRUE(count_matches("y=", cmd_output) >= 1);
  // generated quantities [only on saved iterations, should be num samples]
  EXPECT_TRUE(count_matches("w=", cmd_output) == 1);
}

TEST(StanUiCommand, printReallyPrints) {
  std::vector<std::string> path_vector;
  path_vector.push_back("..");
  path_vector.push_back("src");
  path_vector.push_back("test");
  path_vector.push_back("test-models");
  path_vector.push_back("printer");

  std::string path = "cd test ";
  path += multiple_command_separator();
  path += " ";
  path += convert_model_path(path_vector);
  
  // SAMPLING
  // static HMC
  // + adapt
  test_sample_prints(path + " sample algorithm=hmc engine=static metric=unit_e adapt engaged=0");
  test_sample_prints(path + " sample algorithm=hmc engine=static metric=diag_e adapt engaged=0");
  test_sample_prints(path + " sample algorithm=hmc engine=static metric=dense_e adapt engaged=0");

  // - adapt
  test_sample_prints(path + " sample algorithm=hmc engine=static metric=unit_e adapt engaged=1");
  test_sample_prints(path + " sample algorithm=hmc engine=static metric=diag_e adapt engaged=1");
  test_sample_prints(path + " sample algorithm=hmc engine=static metric=dense_e adapt engaged=1");

  // NUTS
  // + adapt
  test_sample_prints(path + " sample algorithm=hmc engine=nuts metric=unit_e adapt engaged=0");
  test_sample_prints(path + " sample algorithm=hmc engine=nuts metric=diag_e adapt engaged=0");
  test_sample_prints(path + " sample algorithm=hmc engine=nuts metric=dense_e adapt engaged=0");

  // - adapt
  test_sample_prints(path + " sample algorithm=hmc engine=nuts metric=unit_e adapt engaged=1");
  test_sample_prints(path + " sample algorithm=hmc engine=nuts metric=diag_e adapt engaged=1");
  test_sample_prints(path + " sample algorithm=hmc engine=nuts metric=dense_e adapt engaged=1");

  // OPTIMIZATION
  test_optimize_prints(path + " optimize algorithm=newton");
  test_optimize_prints(path + " optimize algorithm=bfgs");
}

TEST(StanUiCommand, refresh_zero_ok) {
  std::vector<std::string> model_path;
  model_path.push_back("src");
  model_path.push_back("test");
  model_path.push_back("test-models");
  model_path.push_back("proper");
  
  std::string command = convert_model_path(model_path) + " sample num_samples=10 num_warmup=10 init=0 output refresh=0 file=test/output.csv";
  run_command_output out = run_command(command);
  EXPECT_EQ(int(stan::services::error_codes::OK), out.err_code);
  EXPECT_EQ(0, count_matches("Iteration:", out.output));
}

TEST(StanUiCommand, refresh_nonzero_ok) {
  std::vector<std::string> model_path;
  model_path.push_back("src");
  model_path.push_back("test");
  model_path.push_back("test-models");
  model_path.push_back("proper");
  
  std::string command = convert_model_path(model_path) + " sample num_samples=10 num_warmup=10 init=0 output refresh=1 file=test/output.csv";
  run_command_output out = run_command(command);
  EXPECT_EQ(int(stan::services::error_codes::OK), out.err_code);
  EXPECT_EQ(20, count_matches("Iteration:", out.output));
}

TEST(StanUiCommand, zero_init_value_fail) {
  std::string expected_message
    = "Rejecting initial value";

  std::vector<std::string> model_path;
  model_path.push_back("src");
  model_path.push_back("test");
  model_path.push_back("test-models");
  model_path.push_back("value_fail");

  std::string command = convert_model_path(model_path) + " sample init=0 output file=test/output.csv";
  run_command_output out = run_command(command);
  EXPECT_EQ(int(stan::services::error_codes::SOFTWARE), out.err_code);

  EXPECT_TRUE(out.header.length() > 0U);
  EXPECT_TRUE(out.body.length() > 0U);
  
  EXPECT_EQ(1, count_matches(expected_message, out.body))
    << "Failed running: " << out.command;
}

TEST(StanUiCommand, zero_init_domain_fail) {
  std::string expected_message
    = "Rejecting initial value";

  std::vector<std::string> model_path;
  model_path.push_back("src");
  model_path.push_back("test");
  model_path.push_back("test-models");
  model_path.push_back("domain_fail");
  
  std::string command = convert_model_path(model_path) + " sample init=0 output file=test/output.csv";
  
  run_command_output out = run_command(command);
  EXPECT_EQ(int(stan::services::error_codes::SOFTWARE), out.err_code);

  EXPECT_TRUE(out.header.length() > 0U);
  EXPECT_TRUE(out.body.length() > 0U);
  
  EXPECT_EQ(1, count_matches(expected_message, out.body))
    << "Failed running: " << out.command;
}

TEST(StanUiCommand, user_init_value_fail) {
  std::string expected_message
    = "Rejecting initial value";

  std::vector<std::string> model_path;
  model_path.push_back("src");
  model_path.push_back("test");
  model_path.push_back("test-models");
  model_path.push_back("value_fail");
  
  std::vector<std::string> init_path;
  init_path.push_back("src");
  init_path.push_back("test");
  init_path.push_back("test-models");
  init_path.push_back("value_fail.init.R");
  
  std::string command = convert_model_path(model_path)
    + " sample init=" + convert_model_path(init_path)
    + " output file=test/output.csv";

  run_command_output out = run_command(command);
  EXPECT_EQ(int(stan::services::error_codes::SOFTWARE), out.err_code);
  
  EXPECT_TRUE(out.header.length() > 0U);
  EXPECT_TRUE(out.body.length() > 0U);
  
  EXPECT_EQ(1, count_matches(expected_message, out.body))
    << "Failed running: " << out.command;
}

TEST(StanUiCommand, user_init_domain_fail) {
  std::string expected_message
    = "Rejecting initial value";

  std::vector<std::string> model_path;
  model_path.push_back("src");
  model_path.push_back("test");
  model_path.push_back("test-models");
  model_path.push_back("domain_fail");
  
  std::vector<std::string> init_path;
  init_path.push_back("src");
  init_path.push_back("test");
  init_path.push_back("test-models");
  init_path.push_back("domain_fail.init.R");
  
  std::string command = convert_model_path(model_path)
    + " sample init=" + convert_model_path(init_path)
    + " output file=test/output.csv";
  
  run_command_output out = run_command(command);
  EXPECT_EQ(int(stan::services::error_codes::SOFTWARE), out.err_code);
  
  EXPECT_TRUE(out.header.length() > 0U);
  EXPECT_TRUE(out.body.length() > 0U);
  
  EXPECT_EQ(1, count_matches(expected_message, out.body))
    << "Failed running: " << out.command;
}

TEST(StanUiCommand, CheckCommand_default) {
  std::vector<std::string> model_path;
  model_path.push_back("src");
  model_path.push_back("test");
  model_path.push_back("test-models");
  model_path.push_back("domain_fail"); // can use any model here
   
  std::string command = convert_model_path(model_path);
  run_command_output out = run_command(command);
  EXPECT_EQ(int(stan::services::error_codes::USAGE), out.err_code);
}

TEST(StanUiCommand, CheckCommand_help) {
  std::vector<std::string> model_path;
  model_path.push_back("src");
  model_path.push_back("test");
  model_path.push_back("test-models");
  model_path.push_back("domain_fail"); // can use any model here
  
  std::string command = convert_model_path(model_path) + " help";

  run_command_output out = run_command(command);
  EXPECT_EQ(int(stan::services::error_codes::OK), out.err_code);
}

TEST(StanUiCommand, CheckCommand_unrecognized_argument) {
  std::vector<std::string> model_path;
  model_path.push_back("src");
  model_path.push_back("test");
  model_path.push_back("test-models");
  model_path.push_back("domain_fail"); // can use any model here
  
  std::string command = convert_model_path(model_path) + " foo";

  run_command_output out = run_command(command);
  EXPECT_EQ(int(stan::services::error_codes::USAGE), out.err_code);
}

TEST(StanUiCommand, timing_info) {
  std::vector<std::string> model_path;
  model_path.push_back("src");
  model_path.push_back("test");
  model_path.push_back("test-models");
  model_path.push_back("proper");
  
  std::string command = convert_model_path(model_path) + " sample num_samples=10 num_warmup=10 init=0 output refresh=0 file=test/output.csv";
  run_command_output out = run_command(command);
  EXPECT_EQ(int(stan::services::error_codes::OK), out.err_code);
  
  std::fstream output_csv_stream("test/output.csv");
  std::stringstream output_sstream;
  output_sstream << output_csv_stream.rdbuf();
  output_csv_stream.close();
  std::string output = output_sstream.str();
  
  EXPECT_EQ(1, count_matches("#  Elapsed Time:", output));
  EXPECT_EQ(1, count_matches(" seconds (Warm-up)", output));
  EXPECT_EQ(1, count_matches(" seconds (Sampling)", output));
  EXPECT_EQ(1, count_matches(" seconds (Total)", output));
}

TEST(StanUiCommand, run_info) {
  std::vector<std::string> model_path;
  model_path.push_back("src");
  model_path.push_back("test");
  model_path.push_back("test-models");
  model_path.push_back("proper");
  
  std::string command = convert_model_path(model_path) + " sample num_samples=10 num_warmup=10 init=0 output refresh=0 file=test/output.csv";
  run_command_output out = run_command(command);
  EXPECT_EQ(int(stan::services::error_codes::OK), out.err_code);
  
  std::fstream output_csv_stream("test/output.csv");
  std::stringstream output_sstream;
  output_sstream << output_csv_stream.rdbuf();
  output_csv_stream.close();
  std::string output = output_sstream.str();
  
  EXPECT_EQ(1, count_matches("# method = sample", output));
  EXPECT_EQ(1, count_matches(" num_samples = 10", output));
  EXPECT_EQ(1, count_matches(" num_warmup = 10", output));
  EXPECT_EQ(1, count_matches(" init = 0", output));
}


// 
struct dummy_stepsize_adaptation {
  void set_mu(const double) {}
  void set_delta(const double) {}
  void set_gamma(const double) {}
  void set_kappa(const double) {}
  void set_t0(const double) {}
};

struct dummy_z {
  Eigen::VectorXd q;
};

template<class ExceptionType>
struct sampler {
  dummy_stepsize_adaptation _stepsize_adaptation;
  dummy_z _z;
  
  double get_nominal_stepsize() {
    return 0;
  }
  dummy_stepsize_adaptation get_stepsize_adaptation() {
    return _stepsize_adaptation;
  }
  
  void engage_adaptation() {}
  dummy_z z() {
    return _z;
  }
  
  void init_stepsize() {
    throw ExceptionType("throwing exception");
  }
};

template<typename T>
class StanUiCommandException : public ::testing::Test {

};
TYPED_TEST_CASE_P(StanUiCommandException);

TYPED_TEST_P(StanUiCommandException, init_adapt) {
  sampler<TypeParam> throwing_sampler;
  Eigen::VectorXd cont_params;
  
  EXPECT_FALSE(stan::services::sample::init_adapt(&throwing_sampler,
                                                  0, 0, 0, 0, cont_params,
                                                  &std::cout));
}

REGISTER_TYPED_TEST_CASE_P(StanUiCommandException,
                           init_adapt);

// exception types that can be thrown by Boost's math functions
typedef ::testing::Types<std::domain_error,
                         std::overflow_error,
                         std::underflow_error,
                         boost::math::rounding_error,
                         boost::math::evaluation_error> BoostExceptionTypes;

INSTANTIATE_TYPED_TEST_CASE_P(, StanUiCommandException, 
                              BoostExceptionTypes);

