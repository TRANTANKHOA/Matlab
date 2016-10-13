#ifndef STAN_LANG_GRAMMARS_EXPRESSION07_GRAMMAR_DEF_HPP
#define STAN_LANG_GRAMMARS_EXPRESSION07_GRAMMAR_DEF_HPP

#include <boost/lexical_cast.hpp>
#include <boost/config/warning_disable.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/std_pair.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_function.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/qi_numeric.hpp>
#include <boost/spirit/include/classic_position_iterator.hpp>
#include <boost/spirit/include/support_multi_pass.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/variant/apply_visitor.hpp>
#include <boost/variant/recursive_variant.hpp>

#include <stan/lang/ast.hpp>
#include <stan/lang/grammars/whitespace_grammar.hpp>
#include <stan/lang/grammars/term_grammar.hpp>
#include <stan/lang/grammars/expression_grammar.hpp>
#include <stan/lang/grammars/expression07_grammar.hpp>

#include <cstddef>
#include <iomanip>
#include <iostream>
#include <istream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <stdexcept>

namespace stan {

  namespace lang {

    // see bare_type_grammar_def.hpp for original
    struct set_val2 {
      template <class> struct result;
      template <typename F, typename T1, typename T2>
      struct result<F(T1, T2)> { typedef void type; };
      template <typename T1, typename T2>
      void operator()(T1& lhs,
                      const T2& rhs) const {
        lhs = rhs;
      }
    };
    boost::phoenix::function<set_val2> set_val2_f;

    struct validate_expr_type3 {
      template <class> struct result;

      template <typename F, typename T1, typename T2, typename T3>
      struct result<F(T1, T2, T3)> { typedef void type; };

      void operator()(const expression& expr,
                      bool& pass,
                      std::ostream& error_msgs) const {
        pass = !expr.expression_type().is_ill_formed();
        if (!pass)
          error_msgs << "expression is ill formed" << std::endl;
      }
    };
    boost::phoenix::function<validate_expr_type3> validate_expr_type3_f;

    // FIXME: cut and paste from term grammar, having trouble w. includes
    struct set_fun_type3 {
      template <typename T1, typename T2>
      struct result { typedef fun type; };

      fun operator()(fun& fun,
                     std::ostream& error_msgs) const {
        std::vector<expr_type> arg_types;
        for (size_t i = 0; i < fun.args_.size(); ++i)
          arg_types.push_back(fun.args_[i].expression_type());
        fun.type_ = function_signatures::instance().get_result_type(fun.name_,
                                                                    arg_types,
                                                                    error_msgs);
        return fun;
      }
    };
    boost::phoenix::function<set_fun_type3> set_fun_type3_f;

    struct addition_expr3 {
      template <class> struct result;

      template <typename F, typename T1, typename T2, typename T3>
      struct result<F(T1, T2, T3)> { typedef void type; };

      void operator()(expression& expr1,
                      const expression& expr2,
                      std::ostream& error_msgs) const {
        if (expr1.expression_type().is_primitive()
            && expr2.expression_type().is_primitive()) {
          expr1 += expr2;
          return;
        }
        std::vector<expression> args;
        args.push_back(expr1);
        args.push_back(expr2);
        set_fun_type3 sft;
        fun f("add", args);
        sft(f, error_msgs);
        expr1 = expression(f);
      }
    };
    boost::phoenix::function<addition_expr3> addition3_f;


    struct subtraction_expr3 {
      template <class> struct result;
      template <typename F, typename T1, typename T2, typename T3>
      struct result<F(T1, T2, T3)> { typedef void type; };

      void operator()(expression& expr1,
                      const expression& expr2,
                      std::ostream& error_msgs) const {
        if (expr1.expression_type().is_primitive()
            && expr2.expression_type().is_primitive()) {
          expr1 -= expr2;
          return;
        }
        std::vector<expression> args;
        args.push_back(expr1);
        args.push_back(expr2);
        set_fun_type3 sft;
        fun f("subtract", args);
        sft(f, error_msgs);
        expr1 = expression(f);
      }
    };
    boost::phoenix::function<subtraction_expr3> subtraction3_f;



    template <typename Iterator>
    expression07_grammar<Iterator>::expression07_grammar(variable_map& var_map,
                                             std::stringstream& error_msgs,
                                             expression_grammar<Iterator>& eg)
      : expression07_grammar::base_type(expression07_r),
        var_map_(var_map),
        error_msgs_(error_msgs),
        term_g(var_map, error_msgs, eg) {
      using boost::spirit::qi::_1;
      using boost::spirit::qi::eps;
      using boost::spirit::qi::lit;
      using boost::spirit::qi::_pass;
      using boost::spirit::qi::_val;
      using boost::spirit::qi::labels::_r1;

      expression07_r.name("expression");
      expression07_r
        %=  term_g(_r1)
            [set_val2_f(_val, _1)]
        > *((lit('+')
             > term_g(_r1)
               [addition3_f(_val, _1, boost::phoenix::ref(error_msgs))])
            |
            (lit('-')
             > term_g(_r1)
             [subtraction3_f(_val, _1, boost::phoenix::ref(error_msgs))]))
        > eps[validate_expr_type3_f(_val, _pass,
                                    boost::phoenix::ref(error_msgs_))];
    }

  }
}
#endif
