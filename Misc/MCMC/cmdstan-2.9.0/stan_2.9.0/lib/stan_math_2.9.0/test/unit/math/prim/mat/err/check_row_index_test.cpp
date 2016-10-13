#include <stan/math/prim/mat/meta/get.hpp>
#include <stan/math/prim/arr/meta/get.hpp>
#include <stan/math/prim/mat/meta/length.hpp>
#include <stan/math/prim/mat/meta/is_vector.hpp>
#include <stan/math/prim/mat/meta/is_vector_like.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>
#include <stan/math/prim/mat/err/check_row_index.hpp>
#include <gtest/gtest.h>

TEST(ErrorHandlingMatrix, checkRowIndexMatrix) {
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> y;
  size_t i;
  
  i=2;
  y.resize(3,3);
  EXPECT_TRUE(stan::math::check_row_index("checkRowIndexMatrix",
                                                    "i", y, i));
  i=3;
  EXPECT_TRUE(stan::math::check_row_index("checkRowIndexMatrix",
                                                    "i", y, i));

  y.resize(2, 3);
  EXPECT_THROW(stan::math::check_row_index("checkRowIndexMatrix",
                                                     "i", y, i), 
               std::out_of_range);

  i=0;
  EXPECT_THROW(stan::math::check_row_index("checkRowIndexMatrix",
                                                     "i", y, i), 
               std::out_of_range);
}

TEST(ErrorHandlingMatrix, checkRowIndexMatrix_nan) {
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> y;
  size_t i;
  double nan = std::numeric_limits<double>::quiet_NaN();

  i=2;
  y.resize(3,3);
  y << nan, nan, nan,nan, nan, nan,nan, nan, nan;
  EXPECT_TRUE(stan::math::check_row_index("checkRowIndexMatrix",
                                                    "i", y, i));
  i=3;
  EXPECT_TRUE(stan::math::check_row_index("checkRowIndexMatrix",
                                                    "i", y, i));

  y.resize(2, 3);
  y << nan, nan, nan,nan, nan, nan;
  EXPECT_THROW(stan::math::check_row_index("checkRowIndexMatrix",
                                                     "i", y, i), 
               std::out_of_range);

  i=0;
  EXPECT_THROW(stan::math::check_row_index("checkRowIndexMatrix",
                                                     "i", y, i), 
               std::out_of_range);
}
